rm(list = ls())
library(readxl)
library(glmnet)
set.seed(0) # for reproducibility

# load data
files <- list.files('../data/', full.names = T)
file <- files[1]

################################################################################
# MAKE FOLDS FOR CROSS-VALIDATION (OR REPEATED RANDOM SPLITS)

makeCVfolds <- function(N, k = 5) {
  
  indices <- sample(N)  # shuffle 1:N
  folds <- vector('list', k)
  
  for (i in seq_len(k)) {
    folds[[i]] <- indices[seq(i, N, by = k)]
  }
  
  folds
  
}

makeRandomSplits <- function(N, f, reps = 100) { # here f is the fraction of in-sample observations
  
  lapply(1:reps,
         FUN = function(i) sample(N, round((1 - f)*N)))
  
}

################################################################################
# FIRST-ORDER MODEL; LASSO/RIDGE REGULARIZATION

predict_lin1st <- function(df, folds, alpha = 1) {
  
  # alpha = 1 -> Lasso regularization; alpha = 0 -> Ridge regularization
  
  y_pred <- vector('list', length(folds))  # will store out-of-sample predictions
  
  # separate in- and out-of-sample
  X <- as.matrix(df[, -ncol(df)])
  y <- df[, ncol(df)]
  
  for (i in seq_along(folds)) {
    
    test_idx  <- folds[[i]]
    train_idx <- setdiff(1:256, folds[[i]])
    
    X_train <- X[train_idx, ]
    y_train <- y[train_idx]
    X_test  <- X[test_idx, ]
    
    # remove zero-variance columns of training set
    X_train <- X_train[, apply(X_train, 2, var) != 0]
    X_test <- X_test[, apply(X_train, 2, var) != 0]
    
    # cross-validate lambda on the training set
    cv_fit <- cv.glmnet(X_train, y_train, alpha = alpha)
    
    # predict on test fold using lambda that minimizes CV error
    y_pred[[i]] <- as.numeric(
      predict(cv_fit, newx = X_test, s = 'lambda.min')
    )
    
  }
  
  y_pred
  
}

################################################################################
# SECOND-ORDER MODEL

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
    which_rm <- apply(X_train, 2, var) == 0
    X_train <- X_train[, !which_rm]
    X_test <- X_test[, !which_rm]
    
    cv_fit <- cv.glmnet(X_train, y_train, alpha = alpha)
    y_pred[[i]] <- as.numeric(
      predict(cv_fit, newx = X_test, s = 'lambda.min')
    )
    
  }
  
  y_pred
  
}

################################################################################
# WRAPPER FUNCTION: for a given set of folds, evaluate all methods

evalMethods <- function(df, s) {
  
  # the 's' parameter is interpreted as K for K-fold cross-validation if it is an integer,
  # and as a fraction f (fraction of in-sample observations) if it is in the range (0,1)
  
  # make folds
  if (s<1) {
    folds <- makeRandomSplits(nrow(df), f = s)
  } else {
    folds <- makeCVfolds(nrow(df), k = s)
  }
  
  # true fitness of out-of-samples
  y_true <- lapply(folds,
                   FUN = function(fold) df[fold, ncol(df)])
  
  # predictions
  ###
  # y_pred <- list(
  #   lin1st_lasso = predict_lin1st(df, folds, alpha = 1), # 1st order model, Lasso regularization
  #   lin2nd_lasso = predict_lin2nd(df, folds, alpha = 1), # 2nd order model, Lasso regularization
  #   lin1st_ridge = predict_lin1st(df, folds, alpha = 0), # 1st order model, Ridge regularization
  #   lin2nd_ridge = predict_lin2nd(df, folds, alpha = 0) # 2nd order model, Ridge regularization
  # )
  ###
  y_pred <- vector('list', length = 4)
  names(y_pred) <- c('lin1st_lasso',
                     'lin2nd_lasso',
                     'lin1st_ridge',
                     'lin2nd_ridge')
  for (m in names(y_pred)) y_pred[[m]] <- vector('list', length = length(folds))
  for (i in 1:length(folds)) {

    print(paste('fold', i))

    y_pred[[1]][[i]] <- predict_lin1st(df, folds[i], alpha = 1)[[1]]
    y_pred[[2]][[i]] <- predict_lin2nd(df, folds[i], alpha = 1)[[1]]
    y_pred[[3]][[i]] <- predict_lin1st(df, folds[i], alpha = 0)[[1]]
    y_pred[[4]][[i]] <- predict_lin2nd(df, folds[i], alpha = 0)[[1]]

  }
  ###
  
  # extract R2, RMSE
  stats <- do.call(rbind,
                   lapply(1:length(folds),
                          FUN = function(i) {

                            do.call(rbind,
                                    lapply(names(y_pred),
                                           FUN = function(m) {

                                             eval_mod <- lm(data.frame(y_true = y_true[[i]],
                                                                       y_pred = y_pred[[m]][[i]]),
                                                            formula = y_true ~ y_pred)

                                             data.frame(method = m,
                                                        fold_id = i,
                                                        R2 = summary(eval_mod)$r.squared,
                                                        R2_xy = 1 - sum((y_true[[i]] - y_pred[[m]][[i]])^2) / sum((y_true[[i]] - mean(y_true[[i]]))^2),
                                                        RMSE = sqrt(mean((y_true[[i]] - y_pred[[m]][[i]])^2)))

                                           }))

                          }))
  
  return(list(
    folds = folds,
    y_true = y_true,
    y_pred = y_pred,
    stats = stats
  ))
  
}

################################################################################
# RUN PIPELINE

# try different K
s <- seq(0.05, 0.7, by = 0.025)


for (file in files) {
  
  print(file)
  
  # load data
  df <- read.csv(file)
  df <- df[, c(1:8, which(colnames(df) == 'mean'))]
  colnames(df)[9] <- 'func'
  
  ###
  # eval_data <- list(s = s,
  #                   results = lapply(s,
  #                                    FUN = function(s_i) evalMethods(df, s_i)))
  ###
  eval_data <- list(s = s,
                    results = vector('list', length = length(s)))
  for (i in 1:length(s)) {
    print(s[i])
    eval_data$results[[i]] <- evalMethods(df, s[i])
  }
  ###
  
  save(eval_data,
       file = paste0('../results/eval_linmod_',
                     gsub('_L8_resources_module_magda.csv', '', basename(file)),
                     '.RData'))
  
}




