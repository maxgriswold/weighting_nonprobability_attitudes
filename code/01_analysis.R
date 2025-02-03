library(survey)
library(betareg)
library(Hmisc) 

### Wrapper function to calculate relevant summary statistics
test.func <- function(all.dat, covars, outvars, G = 40, propens = TRUE, simul = FALSE, w = "WEIGHT2", samp.var = "S", alpha = .05) {
    c.dat <- sapply(all.dat, class)
    facs <- colnames(all.dat)[c.dat == "factor" | c.dat == "character"]
    facs <- intersect(covars, facs)
    out <- data.frame(DEFF = 1)
    tmp <- blend.sample(data = all.dat, samp.var = samp.var, G = G, covars = covars, facs = facs, w = w, propens = propens, simul = simul)
    out$DEFF <- mean(tmp[all.dat$S == 0, 1]^2, na.rm = TRUE)/mean(tmp[all.dat$S == 0, 1], na.rm = TRUE)^2
    out.test <- test.blend(all.dat, tmp, outvars, w.prob = "WEIGHT2")
    out$num.sig <- sum(out.test[, 6] <= alpha)
    out$ks.pval <- ks.test(out.test[, 6], "punif", 0, 1)$p.value
    out1 <- out
    out1$DEFF <- 1
    out1$num.sig <- sum(out.test[, 3] <= alpha)
    out1$ks.pval <- ks.test(out.test[, 3], "punif", 0, 1)$p.value
    rownames(out) <- " "
    return(list(output = out.test, info = out, info.unw = out1))
}

### Outer function used to calculate weights
blend.sample <- function(data, conv = NULL, samp.var = NULL, covars, facs = NULL, G = NULL, w = NULL, trim = NULL, simul = FALSE, propens = FALSE, verbose = FALSE) {

  if(length(facs) > 0 & !is.logical(facs)) {
    facs <- is.element(covars,  facs)
  }

  if (sum(G) == 0 | length(G) == 0) {
    G = 0
  } else if (G == TRUE) {
    G = min(table(strata))
  } 

  w.all = blend.sample.inner(data = data, conv = conv, samp.var = samp.var, covars = covars, facs = facs, w = w, trim = trim, simul = simul, propens = propens, verbose = verbose)

  if(G > 0) {

   if(length(conv) == 0) {

      rep.G <- assign.groups(strata = data[, samp.var],  G = G)

      out <- matrix(NA, NROW(data), G + 1)
      dimnames(out) <- list(names(w.all), c("All", paste0("Group", 1:G, sep = "")))
    
      out[, "All"] <- w.all

      for(g in 1:G) {
         tmp <- blend.sample.inner(data = data[rep.G != g, ], samp.var = samp.var, covars = covars, facs = facs, w = w, trim = trim, simul = simul, propens = propens, verbose = FALSE)
         out[rep.G != g, g + 1] <- tmp
      }

    } else {

      rep.G.prob <- assign.groups(strata = rep(1, NROW(data)),  G = G)
      rep.G.conv <- assign.groups(strata = rep(1, NROW(conv)),  G = G)

      out.prob <- matrix(NA, NROW(data), G+1)
      dimnames(out.prob) <- list(names(w.all$prob), c("All", paste0("Group", 1:G, sep = "")))

      out.conv <- matrix(NA, NROW(conv), G + 1)
      dimnames(out.prob) <- list(names(w.all$conv), c("All", paste0("Group", 1:G, sep = "")))

      out.prob[, "All"] <- w.all$prob
      out.conv[, "All"] <- w.all$conv

      for(g in 1:G) {
         tmp <- blend.sample.inner(data = data[rep.G.prob != g, ], conv=conv[rep.G.conv != g, ], samp.var = samp.var, covars = covars, facs = facs, w = w, trim = trim, simul = simul, propens = propens, verbose = FALSE)
         out.prob[rep.G.prob != g, g+1] <- tmp$prob
         out.conv[rep.G.conv != g, g+1] <- tmp$conv
      }
  
      out <- list(prob = out.prob, conv = out.conv)

    }

  } else {
    out <- w.all
  }

  return(out)

}


### Internal function for calculation of weights
blend.sample.inner <- function(data, conv = NULL, samp.var = NULL, covars, facs = NULL, w = NULL, trim = NULL, simul = FALSE, propens = FALSE, verbose = TRUE) {

  if(length(rownames(data)) == 0){rownames(data) = 1:NROW(data)}

  samps <- NULL
  if(length(conv) == 0){
    samps <- data[, samp.var] == 1
    names(samps) <- rownames(data)
    conv <- data[data[, samp.var] == 0, ]
    data <- data[data[, samp.var] == 1, ]
  }

  if(length(rownames(conv)) == 0){rownames(conv) = 1:NROW(conv)}

  if(length(covars) == 0) {
    data = as.data.frame(data)
    data$myIntercept <- 1
    if(length(conv) > 0){
      conv <- as.data.frame(conv)
      conv$myIntercept <- 1
    }
    covars <- "myIntercept"
  }

  form <- covars

  if(length(facs) > 0) {
    form[facs] <- paste0("factor(", covars[facs], ")")
  }

  form <- paste0("~ ", paste(form, collapse = " + "), sep = "")

  forma <- form
  if(length(w) == 1) {
    forma <- paste0(forma, " + ", w, sep = "")
  }

  form <- as.formula(form)
  forma <- as.formula(forma)

  prob1 <- model.matrix(forma, data = data)
  prob1 <- as.data.frame(prob1)

  conv1 <- model.matrix(form, data = conv)
  conv1 <- as.data.frame(conv1)

  if(length(w) == 1) {
    w <- prob1[, w]
    prob1 <- prob1[, -NCOL(prob1), drop = FALSE]
  } else if (length(w) == 0) {
    w <- rep(1, NROW(prob1))
  } else if (length(w) == NROW(prob1)){
    w  <-  w
  }

  w <- length(w) * w/sum(w)

  benchmarks <- colSums(w * prob1)

  names.tmp <- colnames(conv1)

  colnames(prob1) <- colnames(conv1) <- names(benchmarks) <- paste("X", 1:length(benchmarks), sep = "")
  covars1 <- names(benchmarks)

  all.dat <- rbind(prob1[, covars1, drop = FALSE], conv1[, covars1, drop = FALSE])

  if(!propens) {

    InitW <- NROW(prob1) * rep(1/NROW(conv1), NROW(conv1))

    form1 <- paste0("~ ", paste(colnames(conv1), collapse = " + "), " - 1", sep = "")
    form1 <- as.formula(form1)

    if(!simul){
      w.dat <- conv1
      all.w <- InitW
    } else {
      w.dat <- all.dat
      all.w <- c(w, rep(1, NROW(conv1)))
    }

    colnames(all.dat) <- names.tmp

    caldesign <- svydesign(ids = ~0, data = w.dat, weights = all.w)

    cali <- calibrate(design = caldesign, formula = form1, population = benchmarks, data = all.dat, calfun = "raking")

    myweights <- weights(cali)

    if(simul) {
      w  <-  myweights[1:NROW(prob1)]
      myweights <- myweights[(NROW(prob1) + 1):length(myweights)]
    }

  } else {

    P <- 258300000
    w <- P * w/sum(w)

    d <- 1/w

    prob2 <- data.frame(d = d, prob1)

    bad.col <- colspace(crossprod(as.matrix(prob2)))
    bad.col <- colnames(prob2)[bad.col]

    mod <- betareg(d~.-1, data = prob2[, !is.element(colnames(prob2), bad.col), drop = FALSE])

    d2 <- predict(mod, newdata = data.frame(conv1))

    d1 <- d

    d <- c(d1, d2)
    w.tmp <- 1/d
    w.tmp <- length(w.tmp) * w.tmp/sum(w.tmp)

    all.dat1 <- cbind(S = c(rep(0, NROW(prob1)), rep(1, NROW(conv1))), all.dat)

    mod1 <- glm(S~.-1, data = all.dat1, family = quasibinomial())
    coefs <- summary(mod1)$coefficients
    gamma <- mod1$fitted.values

    if(!simul) {
      q <- d * gamma/(1-gamma)
      w  <-  1/d[1:NROW(prob1)]
      myweights <- 1/q[(NROW(prob1) + 1):length(q)]
    } else {
      p <- d/(1-gamma)
      w  <-  1/p[1:NROW(prob1)]
      myweights <- 1/p[(NROW(prob1) + 1):length(p)]
    }

  }

  if(!simul){
    myweights <- sum(w) * myweights/sum(myweights)
    kappa <- sum(w) * sum(myweights^2)/(sum(w^2) * sum(myweights) + sum(w) * sum(myweights^2))
    w <- w * kappa
    myweights <- (1-kappa) * myweights
  } else {
    w <- (length(w) + length(myweights)) * w/(sum(w) + sum(myweights))
    myweights <- (length(w) + length(myweights)) * myweights/(sum(w) + sum(myweights))
  }

  if(length(trim) > 0) {
    myweights <- c(w, myweights)
    if(length(trim) == 1){trim <- c(trim, 1-trim)}
    cali <- trimWeights(cali,  upper  =  quantile(myweights, trim[2]),  lower  =  quantile(myweights, trim[1]))
    myweights <- weights(cali)
    w  <-  myweights[1:NROW(prob1)]
    myweights <- myweights[(NROW(prob1) + 1):length(myweights)]
  }

  names(myweights) <- rownames(conv1)
  names(w) <- rownames(prob1)

  if(length(samps) != 0) {
    blended <- rep(NA, length(samps))
    names(blended) <- names(samps)
    blended[names(w)] <- w
    blended[names(myweights)] <- myweights
  } else {
    blended <- list(conv = myweights, prob = w)
  }

  cali.tots <- colSums(myweights * conv1)
  tots.unw <- colSums(conv1)

  all.w <- c(w, myweights)
  balance.tab <- cbind(Benchmarks = benchmarks/NROW(prob1), Blended = colSums(all.w * all.dat)/sum(all.w), Conv = cali.tots/sum(myweights), Unw.Conv = tots.unw/NROW(conv1))

  names.tmp <- gsub("factor(", "", names.tmp, fixed = TRUE)
  names.tmp <- gsub(")", "", names.tmp, fixed = TRUE)
  names.tmp[names.tmp == "(Intercept"] <- "(Intercept)"

  rownames(balance.tab) <- names.tmp

  deffs <- c(Benchmarks = mean(c(w)^2)/mean(c(w))^2, Conv = mean(c(myweights)^2)/mean(c(myweights))^2, Blended = mean(c(w, myweights)^2)/mean(c(w, myweights))^2)
  balance.tab <- rbind(balance.tab, DEFF = NA)
  balance.tab["DEFF", names(deffs)] <- deffs
  if(!simul) {
    balance.tab["DEFF", "Unw.Conv"] <- kappa
  }

  if(verbose){print(balance.tab)}

  return(blended)

}

### Given weights, calculates results for outcome variables
test.blend <- function(data, weights, outvars, s.var = "S", w.prob = NULL) {

  drop <- rowMeans(!is.na(weights)) == 0
  data <- data[!drop, ]

  for(i in 1:length(outvars)) {
    if(substr(outvars[i], 1, 3) == "q22") {
      data[which(data[, outvars[i]] == 5), ] <- NA
    }
  }

  if(length(w.prob) > 0) {
    ws <- data[, w.prob]
    ws[is.na(ws)] <- sum(ws[!is.na(ws)]) * rep(1, sum(is.na(ws)))/sum(is.na(ws))
  } else {
    ws <- rep(1, NROW(data))
  }

  if(NCOL(weights) > 1) {
    weights <- weights[!drop, ]
    weights[is.na(weights)] <- 0
    caldesign1 <- svrepdesign(repweights = data.frame(weights[, colnames(weights) != "All"]), weights = weights[, "All"], data = data, type = "JK1")
  } else {
    weights <- weights[!drop]
    weights[is.na(weights)] <- 0
    caldesign1 <- svydesign(ids = ~0, data = data, weights = weights[, "All"])
  }

  caldesign2 <- svydesign(ids = ~0, data = data, weights = ws)

  out <- matrix(NA, length(outvars), 6)
  dimnames(out) <- list(outvars, c("prob.unw", "conv.unw", "p.val.unw", "prob.w", "conv.w", "p.val.w"))

  for(i in 1:length(outvars)) {

    n.vars <- length(table(data[, outvars[i]]))

    myregs0 <- lm(formula(paste(outvars[i], " ~ factor(", s.var, ") - 1", sep = "")), data = data, weights = ws)
    out[i, c(1, 2)] <- summary(myregs0)$coefficients[c(2, 1), "Estimate"]

    if(n.vars>2){
      myregs <- svyglm(formula(paste(outvars[i], " ~ factor(", s.var, ") - 1", sep = "")), design = caldesign1)
      myregs1 <- svyglm(formula(paste(outvars[i], " ~ ", s.var, sep = "")), design = caldesign1)
      myregs2 <- svyglm(formula(paste(outvars[i], " ~ ", s.var, sep = "")), design = caldesign2)
    } else {
      myregs <- svyglm(formula(paste(outvars[i], " ~ factor(", s.var, ") - 1", sep = "")), design = caldesign1, family = quasibinomial())
      myregs1 <- svyglm(formula(paste(outvars[i], " ~ ", s.var, sep = "")), design = caldesign1, family = quasibinomial())
      myregs2 <- svyglm(formula(paste(outvars[i], " ~ ", s.var, sep = "")), design = caldesign2, family = quasibinomial())
    }

    out[i, c(4, 5)] <- summary(myregs)$coefficients[c(2, 1), "Estimate"]
    out[i, 6] <- summary(myregs1)$coefficients[2, 4]
    out[i, 3] <- summary(myregs2)$coefficients[2, 4]

  }

  return(out)

}


### Assigns jackknife groups
assign.groups <- function(strata = NULL, n = length(strata), G = min(table(strata))) {

  # Assign replication groups for jackknife procedures
  states <- names(table(strata))

  Gs <- base::sample(1:G, G)
  Gs1 <- Gs

  rep.G <- rep(NA, n)

  for (i in 1:length(states)) {
    here <- which(strata == states[i])

    n <- length(here)
    samp <- base::sample(1:n, n)

    J <- floor(n/G)

    rep.G1 <- rep(NA, n)

    for (g in 1:G) {
      here1 <- samp[1:J]
      samp <- samp[-(1:J)]
      rep.G1[here1] <- g
    }

    if (length(is.na(rep.G1)) > 0) {
      here2 <- which(is.na(rep.G1))
      if (length(here2) > length(Gs1)) {
        Gs1 <- c(Gs1, Gs)
      }
      rep.G1[here2] <- Gs1[1:length(here2)]
      Gs1 <- Gs1[-(1:length(here2))]
    }

    rep.G[here] <- rep.G1
  }

  return(rep.G)
}


colspace <- function(A, r = nrow(A), k = ncol(A), verbose = TRUE) {

        ### Define RREF within colspace
        RREF <- function(X, ...) {

          GaussianElimination <- function(A, B, tol=sqrt(.Machine$double.eps),
                                        verbose=FALSE, fractions=FALSE){
            # A: coefficient matrix
            # B: right-hand side vector or matrix
            # tol: tolerance for checking for 0 pivot
            # verbose: if TRUE, print intermediate steps
            # fractions: try to express nonintegers as rational numbers
            # If B is absent returns the reduced row-echelon form of A.
            # If B is present, reduces A to RREF carrying B along.

            if ((!is.matrix(A)) || (!is.numeric(A)))
              stop("argument must be a numeric matrix")
            n <- nrow(A)
            m <- ncol(A)
            if (!missing(B)){
              B <- as.matrix(B)
              if (!(nrow(B) == nrow(A)) || !is.numeric(B))
                stop("argument must be numeric and must match the number of row of A")
              A <- cbind(A, B)
            }
            i <- j <- 1
            while (i <= n && j <= m){
              while (j <= m){
                currentColumn <- A[,j]
                currentColumn[1:n < i] <- 0
                # find maximum pivot in current column at or below current row
                which <- which.max(abs(currentColumn))
                pivot <- currentColumn[which]
                if (abs(pivot) <= tol) { # check for 0 pivot
                  j <- j + 1
                  next
                }
                if (which > i) A[c(i, which),] <- A[c(which, i),] # exchange rows
                A[i,] <- A[i,]/pivot # pivot
                row <- A[i,]
                A <- A - outer(A[,j], row) # sweep
                A[i,] <- row # restore current row
                if (verbose) if (fractions) print(fractions(A))
                else print(round(A, round(abs(log(tol,10)))))
                j <- j + 1
                break
              }
              i <- i + 1
            }
            # 0 rows to bottom
            zeros <- which(apply(A[, 1:m], 1, function(x) max(abs(x)) <= tol))
            if (length(zeros) > 0){
              zeroRows <- A[zeros, ]
              A <- A[-zeros, ]
              A <- rbind(A, zeroRows)
            }
            rownames(A) <- NULL
            if (fractions) fractions (A) else round(A, round(abs(log(tol, 10))))
          }

          GaussianElimination(X, ...)
          # returns the reduced row-echelon form of X
        }

        ## Returns the columns that create singularity in the matrix A
        if (k > r) {
          if (verbose){print("WARNING: More columns than rows in colspace.")}
        }
        rA <- RREF(A)
        pivot <- 1
        out <- NULL
        for (i in 2:k) {
          if (rA[(pivot + 1), i] == 1) {
            pivot <- pivot + 1
          } else {
            out <- c(out, i)
          }
        }
        return(out)
      }

### Calculate weighted variances of outcome variables
my.wtd.var <- function(data, outvars, S = "S", w) {

  s.cats <- unique(data[, S])
  s.cats <- s.cats[!is.na(s.cats)]

  out <- matrix(NA, length(outvars), length(s.cats))
  dimnames(out) <- list(outvars, s.cats)

  for(j in 1:length(outvars)) {

    for(k in 1:length(s.cats)) {

      here <- which(data[, S] == s.cats[k])

      out[j, k] <- wtd.var(data[here, outvars[j]], weights = w[here])

    }

  } 

  return(out)

}


outvars <- c("q17", "q18", "q22_1", "q22_2", "q22_3", "q22_4", "q22_5", "q22_6", "q23_1", "q23_2", "q23_3", "q23_4", "q23_5", "q23_6", "q24_1", "q24_2", "q24_3", "q24_4", "q24_5", "q25_1", "q25_2", "q25_3", "q25_4", "q25_5", "q25_6", "q25_7", "q35")

cov.all <- c('age' = "qage", 'gender' = "qgender", 'race' = "qrace1", 'income' = "qincome", 'political ideaology' = "qidea", 'vaccination' = "q53", 'education' = "qedu", 'like to try new products' = "q31_1", 'look for whats new' = "q31_2", 'time spent using technology' = "q32")

out.test <- test.func(all.dat, cov.all, outvars = outvars, G = 80)

ws.tmp <- all.dat$WEIGHT2
ws.tmp[is.na(ws.tmp)] <- 1
vars <- my.wtd.var(all.dat, outvars, w = ws.tmp)

cohens <- cbind(
  Before = abs(out.test$output[, "conv.unw"] - out.test$output[, "prob.unw"])/sqrt(vars[, 1]), 
  After = abs(out.test$output[, "conv.w"] - out.test$output[, "prob.w"])/sqrt(vars[, 1])
)
colMeans(cohens)

sens.tab <- matrix(NA, length(cov.all) + 1, length(out.test$info) + 1)
rownames(sens.tab) <- c("All", names(cov.all))
colnames(sens.tab) <- c(names(out.test$info), "Cohens")
sens.tab[1, ] <- c(as.matrix(out.test$info)[1, ], colMeans(cohens)[2])

for(i in 1:length(cov.all)) {
  out.test.tmp <- test.func(all.dat, cov.all[-i], outvars = outvars, G = 80)
  cohen.tmp <- mean(abs(out.test.tmp$output[, "conv.w"] - out.test.tmp$output[, "prob.w"])/sqrt(vars[, 1]))
  out.test.tmp <- c(as.matrix(out.test.tmp$info)[1, ], cohen.tmp)
  sens.tab[i + 1, ] <- out.test.tmp
}

