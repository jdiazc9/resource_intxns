rm(list = ls())
library(readxl)
library(glmnet)
source('ecoFunctions.R')
set.seed(0) # for reproducibility

# load data
files <- list.files('../data/', full.names = T)

################################################################################
# HELPER FUNCTIONS

makeRandomSplits <- function(N, m, reps = 100) { # here m is the number of in-sample observations
  
  lapply(1:reps,
         FUN = function(i) sample(N, N - m))
  
}

predict_lin2nd <- function(df, folds, alpha = 1) {
  
  # alpha = 1 -> Lasso regularization; alpha = 0 -> Ridge regularization
  
  y_pred <- vector('list', length(folds))
  
  X <- as.matrix(df[, -ncol(df)])
  y <- df[, ncol(df)]
  
  # pairwise interaction terms
  interaction_terms <- combn(ncol(X), 2, simplify = FALSE)
  
  XiXj <- do.call(cbind,
                  lapply(interaction_terms,
                         FUN = function(pair) X[, pair[1]] * X[, pair[2]]))
  
  # interaction names (for traceability)
  colnames(XiXj) <- sapply(interaction_terms,
                           FUN = function(pair) paste(colnames(df)[pair[1]], colnames(df)[pair[2]], sep = '_'))
  
  # full design matrix: main effects + interactions
  X_full <- cbind(X, XiXj)
  
  for (i in seq_along(folds)) {
    
    test_idx  <- folds[[i]]
    train_idx <- setdiff(1:256, folds[[i]])
    
    X_train <- X_full[train_idx, ]
    y_train <- y[train_idx]
    X_test  <- X_full[test_idx, ]
    
    # remove zero-variance columns of training set
    X_train <- X_train[, apply(X_train, 2, var) != 0]
    X_test <- X_test[, apply(X_train, 2, var) != 0]
    
    cv_fit <- cv.glmnet(X_train, y_train, alpha = alpha)
    y_pred[[i]] <- as.numeric(
      predict(cv_fit, newx = X_test, s = 'lambda.min')
    )
    
  }
  
  y_pred
  
}

################################################################################
# RUN PIPELINE

# run models and get stats
df <- do.call(rbind,
              lapply(files,
                     FUN = function(file) {
                       
                       df <- read.csv(file)
                       df <- df[, c(1:8, which(colnames(df) == 'mean'))]
                       colnames(df)[9] <- 'func'
                       
                       folds <- makeRandomSplits(2^8, 1 + 8 + 8*7/2)
                       folds[[101]] <- setdiff(1:2^8,
                                               which(rowSums(df[, 1:8]) <= 2)) # add one fold: all one- and two-resource environments
                       
                       y_true <- lapply(folds,
                                        FUN = function(fold) df[fold, ncol(df)])
                       
                       y_pred <- predict_lin2nd(df, folds)
                       
                       stats <- do.call(rbind,
                                        lapply(1:length(folds),
                                               FUN = function(i) {
                                                 
                                                 eval_mod <- lm(data.frame(y_true = y_true[[i]],
                                                                           y_pred = y_pred[[i]]),
                                                                formula = y_true ~ y_pred)
                                                 
                                                 data.frame(fold_id = i,
                                                            R2 = summary(eval_mod)$r.squared,
                                                            R2_xy = 1 - sum((y_true[[i]] - y_pred[[i]])^2) / sum((y_true[[i]] - mean(y_true[[i]]))^2),
                                                            RMSE = sqrt(mean((y_true[[i]] - y_pred[[i]])^2)))
                                                 
                                               }))
                       
                       return(
                         cbind(
                           isolate = gsub('_L8_resources_module_magda.csv', '', basename(file)),
                           stats
                         )
                       )
                       
                     }))

################################################################################
# PLOT

df$R2_xy <- sapply(df$R2_xy,
                   FUN = function(x) pmax(0, pmin(1, x)))

ggplot(df[df$fold_id != max(df$fold_id), ],
       aes(x = 0, y = R2_xy)) +
  geom_boxplot(outliers = FALSE,
               color = '#1d4070',
               fill = NA) +
  geom_jitter(width = 0.1,
              alpha = 0.2,
              shape = 16,
              color = '#1d4070') +
  geom_hline(data = df[df$fold_id == max(df$fold_id), ],
             aes(yintercept = R2_xy),
             linetype = 'dashed') +
  facet_wrap(~isolate, nrow = 1) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(name = expression(italic(R)^2)) +
  theme_bw() +
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = 'none')

ggsave(filename = paste('../plots/', gsub('-', '', Sys.Date()), '-lin2nd_lasso_onetwos.pdf', sep = ''),
       width = 220,
       height = 140,
       units = 'mm',
       limitsize = F)


