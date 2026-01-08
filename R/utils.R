tune_rf <- function(Utrain, Ztrain, Ztest, regr_par){
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
  Ztrain_data$U <- Utrain
  form <- "U ~ "
  for(j in 1:d){
    form <- paste(form, " + Z", j, sep = "")
  }
  OOBs <- rep(Inf, nrow(par_grid))
  ind_min <- NA
  fit_res <- NULL
  for(j in 1:NROW(par_grid)){
    fit_temp <- ranger::ranger(formula = formula(form), data = Ztrain_data, num.trees = par_grid[j, 1],
                               mtry = par_grid[j, 2], max.depth = par_grid[j, 3],
                               min.node.size = par_grid[j, 4], num.threads = 1)
    OOBs[j] <- fit_temp$prediction.error
    if(OOBs[j] == min(OOBs)){
      ind_min <- j
      fit_res <- fit_temp
    }
  }
  pred_train <- fit_res$predictions
  Ztest_data <- as.data.frame(Ztest)
  colnames(Ztest_data) <- paste("Z", 1:NCOL(Ztrain), sep = "")
  pred_test <- predict(fit_res, data = Ztest_data)$predictions
  par_opt <- list(num.trees = par_grid[ind_min, 1], mtry = par_grid[ind_min, 2], max.depth = par_grid[ind_min, 3], min.node.size = par_grid[ind_min, 4])
  return(list(pred_train = pred_train, pred_test = pred_test, par_opt = par_opt))
}

get_rf_predictions_from_tuned <- function(Utrain, Ztrain, Ztest, tuned_regr_par){
  d <- NCOL(Ztrain)
  Ztrain_data <- as.data.frame(Ztrain)
  colnames(Ztrain_data) <- paste("Z", 1:NCOL(Ztrain), sep = "")
  Ztrain_data$U <- Utrain
  form <- "U ~ "
  for(j in 1:d){
    form <- paste(form, " + Z", j, sep = "")
  }
  fit_res <- ranger::ranger(formula = formula(form), data = Ztrain_data, num.trees = tuned_regr_par$num.trees,
                            mtry = tuned_regr_par$mtry, max.depth = tuned_regr_par$max.depth,
                            min.node.size = tuned_regr_par$min.node.size, num.threads = 1)
  pred_train <- fit_res$predictions
  Ztest_data <- as.data.frame(Ztest)
  colnames(Ztest_data) <- paste("Z", 1:NCOL(Ztrain), sep = "")
  pred_test <- predict(fit_res, data = Ztest_data)$predictions
  return(list(pred_train = pred_train, pred_test = pred_test))
}
