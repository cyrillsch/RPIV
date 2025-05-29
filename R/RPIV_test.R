
get_rf_predictions <- function(res_train, Ztrain, Ztest, regr_par){
  d <- NCOL(Ztrain)
  # reusing some code from the IVDML package to tune the random forest
  if(!exists("num.trees", regr_par)){
    regr_par$num.trees  <- 500
  }
  # number of mtry values to try, default 5
  if(!exists("num_mtry", regr_par)){
    regr_par$num_mtry <- 15
  }
  if(!exists("mtry", regr_par)){
    # default: try some selection of num_mtry integers in [min(d/3, sqrt(d)), 2*d/3]
    regr_par$mtry <- unique(round(seq(max(1, floor(min(sqrt(d),d/3))), ceiling(2*d/3), length.out = regr_par$num_mtry)))
  }
  # by default, we do not restrict max_depth, but regularize min_node_size
  if(!exists("max.depth", regr_par)){
    regr_par$max.depth <-  0
  }
  # number of min.node.size values to try, default: 15
  if(!exists("num_min.node.size", regr_par)){
    regr_par$num_min.node.size <- 15
  }
  # default: try some selection  of num_min.node.size integers from a logarithmic grid up to nX/4
  if(!exists("min.node.size", regr_par)){
    ntrain <- NROW(Ztrain)
    regr_par$min.node.size <- unique(round(5 * exp(seq(0, log(ntrain/20), length.out = regr_par$num_min.node.size))))
  }
  par_grid <- with(regr_par, expand.grid(num.trees = num.trees, mtry = mtry, max.depth = max.depth, min.node.size = min.node.size))
  Ztrain_data <- as.data.frame(Ztrain)
  colnames(Ztrain_data) <- paste("Z", 1:NCOL(Ztrain), sep = "")
  Ztrain_data$res <- res_train
  form <- "res ~ "
  for(j in 1:d){
    form <- paste(form, " + Z", j, sep = "")
  }

  min_OOB <- Inf
  fit_res <- NULL
  for(j in 1:NROW(par_grid)){
    fit_temp <- ranger::ranger(formula = formula(form), data = Ztrain_data, num.trees = par_grid[j, 1],
                       mtry = par_grid[j, 2], max.depth = par_grid[j, 3],
                       min.node.size = par_grid[j, 4], num.threads = 1)
    if(fit_temp$prediction.error < min_OOB){
      min_OOB <- fit_temp$prediction.error
      fit_res <- fit_temp
    }
  }
  pred_train <- fit_res$predictions
  Ztest_data <- as.data.frame(Ztest)
  colnames(Ztest_data) <- paste("Z", 1:NCOL(Ztrain), sep = "")
  pred_test <- predict(fit_res, data = Ztest_data)$predictions
  return(list(pred_train = pred_train, pred_test = pred_test))
}




get_w <- function(res_train, Ztrain, Ztest, upper_clip_quantile, regr_par){
  list_pred <- get_rf_predictions(res_train, Ztrain, Ztest, regr_par)
  pred_test <- list_pred$pred_test
  pred_train <- list_pred$pred_train
  if(upper_clip_quantile == 0){
    w <- sign(pred_test)
  } else {
    Kmax <- quantile(abs(pred_train), upper_clip_quantile)
    w <- sign(pred_test) * pmin(abs(pred_test), Kmax) / Kmax
  }
  return(w)
}


calc_sigmahatw2 <- function(variance_estimator, res, w, ZAw, clustering_test = NULL){
  if(variance_estimator == "homoskedastic"){
    return(mean((w+ZAw)^2) * mean(res^2))
  }
  if(variance_estimator == "heteroskedastic"){
    return(mean((w + ZAw)^2 * res^2) - mean(w * res)^2)
  }
  if(variance_estimator == "cluster"){
    Sg <- tapply((w+ ZAw) * res, clustering_test, sum)
    n0 <- length(w)
    return(sum(Sg^2) / n0 - n0/length(unique(clustering_test))*mean(w*res)^2)
  }
}



#' Residual Prediction Test for Linear Instrumental Variable Models
#'
#' Performs a hypothesis test for the well-specification of linear instrumental variable (IV) model.
#' More specifically, it tests the null-hypothesis
#' \eqn{H_0: \exists\beta\in \mathbb R^p \text{ s.t. } \mathbb E[Y-X^T\beta|Z] = 0.}
#' It uses sample splitting and a random forest to try to predict the two-stage
#' least-squares residuals from the instruments \eqn{Z}.
#'
#' @param Y A numeric vector. The outcome variable.
#' @param X A numeric matrix or vector. The endogenous explanatory variables.
#' @param C A numeric matrix, vector or `NULL`. The additional exogenous explanatory variables (optional).
#' @param Z A numeric matrix or vector. The instruments.
#' @param frac_A A numeric scalar between 0 and 1 or `NULL`. The fraction of the sample used for training (sample splitting). Default is `min(0.5, exp(1)/log(n))`, where `n` is the sample size.
#' @param gamma A non-negative scalar. If the variance estimator is less than gamma times the noise level (as estimated as by the mean of the squared residuals), gamma times the noise level is used as variance estimator.
#' @param variance_estimator Character string or vector. One or more of `"homoskedastic"`, `"heteroskedastic"`, `"cluster"`. Specifies the types of variance estimation used.
#' @param clustering A vector of cluster identifiers or `NULL`. Observations with the same value of `clustering` belong to the same cluster. Required if `variance_estimator` includes `"cluster"`.
#' @param upper_clip_quantile A scalar between 0 and 1. The estimated weight-function will be clipped at the corresponding quantile of the random forest predictions on the auxiliary sample. Use 0 to use the sign of the predictions. Default is 0.8.
#' @param regr_par A list of parameters passed to the random forest regression model. Supports `num.trees`, `num_mtry` (number of different mtry values to try out) or a vector `mtry`, a vector `max.depth`, `num_min.node.size` (number of different min.node.size values to try out) or a vector `min.node.size`.
#' @param fit_intercept Logical. Should an intercept be included in the model? Default is `TRUE`.
#'
#' @return If a single variance estimator is used, returns a list with:
#' \describe{
#'   \item{p_value}{p-value of the residual prediction test.}
#'   \item{test_statistic}{The value of the test statistic.}
#'   \item{var_fraction}{The estimated variance fraction, i.e., variance estimator divided by noise level estmate.}
#'   \item{T_null}{The value of the initial test statistic. If var_fraction >= gamma, it is equal to test_statistic, otherwise, it has larger absolute value.}
#'   \item{variance_estimator}{The variance estimator used.}
#' }
#' If multiple estimators are supplied, returns a named list of such results for each estimator.
#'
#' @details
#' The RPIV test splits the sample into an auxiliary and a main sample.
#' On the auxiliary sample, a random forest is used to predict the two-stage least squares residuals from the instruments.
#' The test statistic is the scalar product of the two-stage least-squares residuals with a clipped
#' and rescaled version of the learned function evaluated on the main sample divided by an estimator of its standard deviation.
#'
#' If `clustering` is supplied, sample splitting is done at cluster level (also for `variance_estimator` "homoskedastic"` or `"heteroskedastic"`).
#'
#'@examples
#' set.seed(1)
#' n <- 100
#' Z <- rnorm(n)
#' H <- rnorm(n)
#' C <- rnorm(n)
#' X <- Z + rnorm(n) + H
#' Y1 <- X - C - H + rnorm(n)
#' Y2 <- X - C - H + rnorm(n) + Z^2
#' RPIV_test(Y1, X, C, Z)
#' RPIV_test(Y2, X, C, Z)
#'
#'
#'
#' @references
#' Cyrill Scheidegger, Malte Londschien and Peter Bühlmann. A residual prediction test for the well-specification of linear instrumental variable models. Preprint, arXiv:????, 2025.
#'
#' @importFrom stats formula lm pnorm predict quantile
#' @import ranger
#' @export


RPIV_test <- function(Y, X, C = NULL, Z, frac_A = NULL, gamma = 0.05,
                      variance_estimator = "heteroskedastic", clustering = NULL,
                      upper_clip_quantile = 0.8, regr_par = list(), fit_intercept = TRUE){
  # Check input and convert to numeric or matrix
  Y <- try(as.numeric(Y), silent = TRUE)
  if (inherits(Y, "try-error")){
    stop("Y cannot be converted to a vector.")
  }
  N <- length(Y)
  matrix_ZXC <- function(var, var_name){
    if (!is.null(var)){
      mat <- try(as.matrix(var), silent = TRUE)
      if (inherits(mat, "try-error")) {
        stop(paste(var_name, " cannot be converted to a matrix. Make sure that it is a vector, matrix or data frame.", sep = ""))
      }
      if (nrow(mat) != N) {
        stop(paste("The number of rows in ", var_name, " must match the length of Y.", sep = ""))
      }
      return(mat)
    } else {
      return(NULL)
    }
  }
  Z <- matrix_ZXC(Z, "Z")
  X <- matrix_ZXC(X, "X")
  C <- matrix_ZXC(C, "C")
  if (!is.null(clustering)) {
    if (length(clustering) != N) stop("clustering must have the same length as Y.")
    if (!is.character(clustering) && !is.factor(clustering) && !is.numeric(clustering)) {
      stop("clustering must be a vector of identifiers (numeric, factor, or character).")
    }
    if (is.factor(clustering)){
      clustering <- as.numeric(clustering)
    }
    if (!("cluster" %in% variance_estimator)) {
      warning("Cluster structure is present but 'cluster' is not in variance_estimator.")
    }
  }
  if (!is.null(frac_A)) {
    if (!is.numeric(frac_A) || length(frac_A) != 1 || frac_A <= 0 || frac_A >= 1) {
      stop("frac_A must be NULL or a numeric scalar in (0, 1).")
    }
  }
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("gamma must be a non-negative numeric scalar.")
  }

  allowed_estimators <- c("homoskedastic", "heteroskedastic", "cluster")
  if (!all(variance_estimator %in% allowed_estimators)) {
    stop(sprintf("variance_estimator must be one or more of %s.", paste(allowed_estimators, collapse = ", ")))
  }

  if (!is.numeric(upper_clip_quantile) || length(upper_clip_quantile) != 1 || upper_clip_quantile < 0 || upper_clip_quantile > 1) {
    stop("upper_clip_quantile must be a numeric scalar in [0, 1].")
  }

  if (!is.list(regr_par)) {
    stop("regr_par must be a list.")
  }

  if (NCOL(X) > NCOL(Z)){
    stop("Z must have at least as many columns as X.")
  }

  # Check regr_par keys
  allowed_regr_par_keys <- c("num.trees", "num_mtry", "mtry", "max.depth", "num_min.node.size", "min.node.size")

  invalid_keys <- setdiff(names(regr_par), allowed_regr_par_keys)
  if (length(invalid_keys) > 0) {
    stop(sprintf(
      "regr_par contains invalid entries: %s. Allowed entries are: %s.",
      paste(invalid_keys, collapse = ", "),
      paste(allowed_regr_par_keys, collapse = ", ")
    ))
  }

  if (!is.logical(fit_intercept) || length(fit_intercept) != 1) {
    stop("fit_intercept must be a logical scalar (TRUE or FALSE).")
  }


  if(is.null(frac_A)){
    frac_A <- min(0.5, exp(1) / log(N))
  }
  # full matrix of regressors
  Xbar <- cbind(X, C)
  # full matrix of instruments
  Zbar <- cbind(Z, C)
  if(fit_intercept){
    Xbar <- cbind(rep(1, N), Xbar)
    Zbar <- cbind(rep(1, N), Zbar)
  }

  if(!is.null(clustering)){
    clusters <- unique(clustering)
    clusters_train <- sample(clusters, round(length(clusters) * frac_A))
    train <- which(clustering %in% clusters_train)
    clustering_test <- clustering[-train]
  } else {
    train <- sample(1:N, round(N * frac_A))
    clustering_test <- NULL
  }
  Ytrain <- Y[train]
  Ytest <- Y[-train]
  Xbartrain <- Xbar[train, ]
  Xbartest <- Xbar[-train, ]
  Zbartrain <- Zbar[train, ]
  Zbartest <- Zbar[-train, ]

  # TSLS on training set
  first_stage_train <- lm(Xbartrain ~ -1 + Zbartrain)
  fitted_train <- first_stage_train$fitted.values
  second_stage_train <- lm(Ytrain ~ -1 + fitted_train)
  beta_train <- second_stage_train$coefficients
  res_train <- Ytrain - Xbartrain %*% beta_train
  if(fit_intercept){
    w <- get_w(res_train, Zbartrain[, -1], Zbartest[, -1], upper_clip_quantile, regr_par)
  } else {
    w <- get_w(res_train, Zbartrain, Zbartest, upper_clip_quantile, regr_par)
  }
  # TSLS on test set
  first_stage_test <- lm(Xbartest ~ -1 + Zbartest)
  fitted_test <- first_stage_test$fitted.values
  second_stage_test <- lm(Ytest ~ -1 + fitted_test)
  beta_test <- second_stage_test$coefficients

  res_test <- Ytest - Xbartest %*% beta_test
  ZAw <- - fitted_test %*%  solve(t(fitted_test) %*% fitted_test, t(Xbartest) %*% w)

  n0 <- length(w)

  list_results <- list()
  for(v_e in variance_estimator){
    sigmahatw2 <- calc_sigmahatw2(v_e, res_test, w, ZAw, clustering_test)
    var_fraction <- sigmahatw2 / mean(res_test^2)
    T_null <- sum(w * res_test) / sqrt(n0 * sigmahatw2)
    if(var_fraction < gamma){
      sigmahatw2 <- gamma * mean(res_test^2)
    }
    test_statistic <- sum(w * res_test) / sqrt(n0 * sigmahatw2)
    pw <- 1 - pnorm(test_statistic)
    list_results[[v_e]] <- list(p_value = pw, test_statistic = test_statistic,
                                var_fraction = var_fraction, T_null = T_null,
                                variance_estimator = v_e)
  }
  if(length(variance_estimator) == 1){
    return(list_results[[1]])
  } else {
    return(list_results)
  }
}
