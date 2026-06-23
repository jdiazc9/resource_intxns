library(ggplot2)
library(tune)
library(scales)

################################################################################
# AESTHETIC OPTIONS

mycolors <- setNames(
  c('#52c8f7',
    '#5283f7',
    '#56e016',
    '#2e780c'),
  c('1st order linear regression - Lasso',
    '2nd order linear regression - Lasso',
    '1st order linear regression - Ridge',
    '2nd order linear regression - Ridge')
)


################################################################################
# LIST RESULTS FILES

files <- list.files('../results/', full.names = T, pattern = 'linmod')

for (file in files) {

################################################################################
# LOAD DATA

load(file)
s <- eval_data$s

################################################################################
# PROCESS STATS

eval_stats <- do.call(rbind,
                      lapply(seq_along(s),
                             FUN = function(i) cbind(s = s[i],
                                                     eval_data[['results']][[i]]$stats)))

################################################################################
# PLOT: PREDICTED VS. OBSERVED, DIFFERENT METHODS

for (si_plot in c(1, 6, 14)) {
  
  plot_this <- do.call(rbind,
                       lapply(names(eval_data$results[[si_plot]]$y_pred),
                              FUN = function(m) {
                                data.frame(s = s[si_plot],
                                           method = m,
                                           y_true = eval_data$results[[si_plot]]$y_true[[1]],
                                           y_pred = eval_data$results[[si_plot]]$y_pred[[m]][[1]])
                              }))
  plot_this$method <- factor(plot_this$method,
                             levels = names(eval_data[['results']][[1]][['y_pred']]))
  
  ggplot(plot_this,
         aes(x = y_pred, y = y_true, color = method)) +
    geom_abline(slope = 1,
                intercept = 0,
                color = 'gray') +
    geom_point() +
    scale_x_continuous(name = 'Predicted fitness',
                       breaks = pretty_breaks(n = 2)) +
    scale_y_continuous(name = 'Measured fitness',
                       breaks = pretty_breaks(n = 2)) +
    scale_color_manual(name = 'Method',
                       values = as.character(mycolors)) +
    facet_wrap(~method,
               nrow = 1, labeller = as_labeller(
                 setNames(
                   c('1st order\nlinear regression\n(Lasso)',
                     '2nd order\nlinear regression\n(Lasso)',
                     '1st order\nlinear regression\n(Ridge)',
                     '2nd order\nlinear regression\n(Ridge)',
                     'FEE\nconcatenation'),
                   levels(plot_this$method)
                 )
               )) +
    coord_obs_pred() +
    theme_bw() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 16),
          legend.position = 'none') +
    ggtitle(paste(gsub('.RData', '', gsub('eval_linmod_', '', basename(file))),
                  '\n',
                  'training fraction = ', s[si_plot],
                  sep = ''))
  
  ggsave(filename = paste('../plots/linmods/', gsub('-', '', Sys.Date()), '-pred-vs-obs-s=', s[si_plot], '-', gsub('.RData', '', gsub('eval_linmod_', '', basename(file))), '.pdf', sep = ''),
         width = 220,
         height = 140,
         units = 'mm',
         limitsize = F)
  
}

################################################################################
# PLOT: R^2 VS. TRAINING FRACTION

eval_stats_plot <- do.call(data.frame,
                           aggregate(R2_xy ~ s + method,
                                     data = eval_stats,
                                     FUN = function(x) c(mean = mean(x),
                                                         median = median(x),
                                                         sd = sd(x),
                                                         q5 = pmax(0, as.numeric(quantile(x, 0.05))),
                                                         q95 = as.numeric(quantile(x, 0.95)))))
eval_stats_plot$s <- eval_stats_plot$s +
  0.05 * min(diff(s)) * setNames(seq(-1, 1, length.out = length(unique(eval_stats_plot$method))),
                                 names(eval_data[['results']][[1]][['y_pred']]))[eval_stats_plot$method]
eval_stats_plot$method <- factor(
  eval_stats_plot$method,
  levels = names(eval_data[['results']][[1]][['y_pred']])
)

ggplot(eval_stats_plot,
       aes(x = s, y = R2_xy.median, color = method,
           ymin = R2_xy.q5, ymax = R2_xy.q95)) +
  geom_line() +
  geom_pointrange(position = 'jitter') +
  scale_x_continuous(name = 'Fraction of observations\nin training set') +
  scale_y_continuous(name = expression(italic(R)^2),
                     limits = c(min(eval_stats_plot$R2_xy.q5) - 0.01, 1)) +
  scale_color_manual(name = 'Method',
                     values = as.character(mycolors),
                     labels = names(mycolors)) +
  theme_bw() +
  theme(aspect.ratio = 0.8,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  ggtitle(gsub('.RData', '', gsub('eval_linmod_', '', basename(file))))

ggplot(eval_stats_plot,
       aes(x = s, y = R2_xy.median, color = method, fill = method,
           ymin = R2_xy.q5, ymax = R2_xy.q95)) +
  geom_ribbon(alpha = 0.2,
              color = NA) +
  geom_line() +
  scale_x_continuous(name = 'Fraction of observations\nin training set',
                     expand = c(0, 0)) +
  scale_y_continuous(name = expression(italic(R)^2),
                     limits = c(min(eval_stats_plot$R2_xy.q5) - 0.01, 1)) +
  scale_color_manual(name = 'Method',
                     values = as.character(mycolors),
                     labels = names(mycolors)) +
  scale_fill_manual(name = 'Method',
                    values = as.character(mycolors),
                    labels = names(mycolors)) +
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  guides(color = 'none') +
  ggtitle(gsub('.RData', '', gsub('eval_linmod_', '', basename(file))))

ggsave(filename = paste('../plots/linmods/', gsub('-', '', Sys.Date()), '-R2-vs-fraction-', gsub('.RData', '', gsub('eval_linmod_', '', basename(file))), '.pdf', sep = ''),
       width = 180,
       height = 140,
       units = 'mm',
       limitsize = F)

}
