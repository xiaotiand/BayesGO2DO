
#' Cross-validated prediction function
#'
#' This is the prediction function with built-in cross validation.
#' @details
#' This function performs leave-one-out cross-validation and generates predicted values.
#'
#' @param res The model output from the Main() function
#' @param Y A list of response values in regressions
#' @param X A list of design matrices in regressions
#' @param burnin The number of burn-in iterations
#' @param ich The index of chain used for prediction (e.g., ich = 1)
#' @return A list of predicted survival times for all subjects
#'
#' @export

Predict = function(res, Y, X, burnin, nfold = NULL, ich = 1) {
  require(caret)
  require(mvtnorm)

  if(is.null(nfold)) nfold = 10
  GAMMAres = res$GAMMA_vec[[ich]][(burnin + 1):length(res$GAMMA_vec[[ich]])]
  Yaug = list()
  for(k in 1:K) {
    Yres = lapply(res$Y_vec[[ich]], function(x) x[[k]])
    Yres = Yres[(burnin + 1):length(res$Y_vec[[ich]])]
    Yaug[[k]] = Reduce("+", Yres) / length(Yres)
  }
  Ypred = lapply(1:K, function(k) rep(0, length(Yaug[[k]])))
  folds = lapply(Y, function(x) createFolds(x, k = nfold, list = FALSE, returnTrain = FALSE))

  for (k in 1:K) {
    for (m in 1:nfold) {
      Xtrain = X[[k]][folds[[k]] != m, ]
      Ytrain = Y[[k]][folds[[k]] != m]
      deltrain = delta[[k]][folds[[k]] != m]
      Xtest = X[[k]][folds[[k]] == m, ]
      Yhtest = Yaug[[k]][folds[[k]] == m]
      W = rep(0, length(GAMMAres))
      pred_list = list()
      for (i in 1:length(GAMMAres)) {
        gamma = which(GAMMAres[[i]][k, ] == 1)
        fit = glmnet(Xtrain[, gamma], as.matrix(data.frame(time = exp(Ytrain), status = deltrain)),
                     alpha = 1, family = "cox", intercept = FALSE)
        Bhat = as.numeric(coef(fit, s = fit$lambda[length(fit$lambda)]))
        if (length(gamma) == 1) {
          pred_list[[i]] = Xtest[, gamma] * Bhat
        } else {
          pred_list[[i]] = Xtest[, gamma] %*% Bhat
        }
        W[i] = dmvnorm(Yhtest, pred_list[[i]], sigma = diag(length(pred_list[[i]])))
      }
      pred_list = lapply(1:length(pred_list), function(i) pred_list[[i]])
      pred_list = Reduce("+", pred_list) / length(pred_list)
      Ypred[[k]][folds[[k]] == m] = pred_list
    }
  }
  return(Ypred)
}

