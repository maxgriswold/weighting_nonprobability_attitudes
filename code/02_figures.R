# Construct plot of blended score differences
# Max Griswold
# 12/7/23

library(data.table)
library(ggplot2)

setwd("C:/users/griswold/documents/GitHub/twitter-representative-pop/")

unit_scale <- T

df <- fread("blended_scores_v2.csv")
df <- df[, .(variable, prob.unw, conv.unw, p.val.unw, 
             prob.w, conv.w, p.val.w, 
             scale_l, scale_u, reverse)]

# Add on indicator for question:
df[1:8, question := "Views of Politicians"]
df[c(9:19, 27), question := "Policy Issues and Political Topics"]
df[20:26, question := "Critical Threats"]

df_unweighted <- df[, .(variable, prob.unw, conv.unw, p.val.unw, question, scale_l, scale_u, reverse)]
setnames(df_unweighted, c("conv.unw", "p.val.unw"), c("convenience_mean", "p_diff"))
df_unweighted[, weighted := "Unadjusted"]

df_weighted <- df[, .(variable, prob.unw, conv.w, p.val.w, question, scale_l, scale_u, reverse)]
setnames(df_weighted, c("conv.w", "p.val.w"), c("convenience_mean", "p_diff"))
df_weighted[, weighted := "Adjusted"]

df <- rbind(df_unweighted, df_weighted)

setnames(df, "prob.unw", "probability_mean")

# Reverse code the regulation and critical threat questions:
df[reverse == 1, `:=`(probability_mean = -1*probability_mean, convenience_mean = -1*convenience_mean)]

normalize <- function(mu, l, u){
  return(4*(mu - l)/(u - l) + 1)
}

calc_ci <- function(p, mu, sig = 0.05, upper = F){
  
  z <- qnorm(p/2)
  z_ci <- abs(qnorm(sig/2))
  
  se <- abs(mu/z)
  
  if (upper == F){
    ci <- mu - z_ci*se
  }else{
    ci <- mu + z_ci*se
  }
  
  return(ci)
  
}

if (unit_scale == T){
  df[, probability_mean := normalize(probability_mean, scale_l, scale_u)]
  df[, convenience_mean := normalize(convenience_mean, scale_l, scale_u)]
}

df[, mean_diff := probability_mean - convenience_mean]
df[, lower_diff := calc_ci(p_diff, mean_diff)]
df[, upper_diff := calc_ci(p_diff, mean_diff, upper = T)]

df[, sig := ifelse(lower_diff < 0 & upper_diff > 0, F, T)]

# Order the plot so the mean difference is cascading from largest to smallest
# difference

#order_vars <- rev(df[weighted == "Unweighted",][df[weighted == "Unweighted", order(mean_diff)]]$variable)

#df[, variable := factor(variable, levels = order_vars)]

plot_colors <- c("#969696", "#252525")

df[, variable := factor(variable, levels = unique(df$variable))]
df[, question := factor(question, levels = unique(df$question))]
df[, weighted := factor(weighted, levels = c("Unadjusted", "Adjusted"))]

if(unit_scale){
  scale_lims <- c(-2, 2)
}else{
  scale_lims <- c(-5, 5)
}

p <- ggplot(df, aes(x = mean_diff, y = variable)) +
      facet_grid(question ~ weighted, space = 'free', scale = 'free') +
      geom_vline(xintercept = 0, linetype = 2, alpha = 0.5, size = 0.4) +
      geom_linerange(aes(xmin = lower_diff, xmax = upper_diff), size = 0.8, alpha =  0.6) +
      geom_point(size = 2.2) +
      labs(title = "",
           x = "Mean Difference",
           y = "") +
      theme_bw() +
      scale_y_discrete(limits = rev) +
      lims(x = scale_lims) + 
      scale_color_manual(values = plot_colors) +
      coord_cartesian(clip = "off") +
      theme(plot.title = element_text(hjust = 0.5, family = 'sans', size = 16),
            strip.background = element_blank(),
            strip.text.x = element_text(family = 'sans', size = 12),
            legend.position = "none",
            legend.text = element_text(family = 'sans', size = 8),
            axis.ticks = element_line(linewidth = 0.5),
            axis.ticks.length = unit(5.6, "points"),
            axis.title = element_text(size = 11, family = 'sans'),
            axis.title.y = element_text(size = 11, family = 'sans', angle = 0),
            axis.text = element_text(size = 12, family = 'sans'),
            axis.text.x = element_text(size = 12, family = 'sans',
                                       margin = margin(t = 0, r = 0, b = 0, l = 0)),
            legend.title = element_text(family = 'sans', size = 10),
            strip.text.y = element_text(family = 'sans', size = 12), 
            plot.margin = margin(0.25, 0.25, 0.25, 1.25, unit = "in"),
            plot.caption = element_text(size = 6, hjust = 0))

subtitle <- ifelse(unit_scale, "_normalized", "")

ggsave(sprintf("adjusted_scores%s.pdf", subtitle), plot = p, device = pdf, width = 11.69, height = 8.27, units = 'in')
ggsave(sprintf("adjusted_scores%s.tiff", subtitle), plot = p, device = tiff, width = 8.27, height = 11.69, units = 'in', dpi = 600)
