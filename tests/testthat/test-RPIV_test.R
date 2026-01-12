test_that("RPIV_test returns correct structure for default arguments", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X - C  + rnorm(n)
  result <- RPIV_test(Y, X, C, Z)
  expect_type(result, "list")
  expect_named(result, c("p_value", "test_statistic", "var_fraction", "T_null", "variance_estimator"))
  expect_true(is.numeric(result$p_value))
  expect_true(is.numeric(result$test_statistic))
  expect_true(is.numeric(result$T_null))
  expect_true(is.numeric(result$var_fraction))
  expect_equal(result$variance_estimator, "heteroskedastic")
})

test_that("RPIV_test works with C = NULL", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)
  expect_no_error(RPIV_test(Y = Y, X = X, Z = Z))
})

test_that("RPIV_test handles multiple variance estimators", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)
  result <- RPIV_test(Y, X, Z = Z, variance_estimator = c("homoskedastic", "heteroskedastic"))
  expect_type(result, "list")
  expect_named(result, c("homoskedastic", "heteroskedastic"))
  expect_true(all(c("p_value", "test_statistic", "var_fraction", "T_null", "variance_estimator") %in%
                    names(result$homoskedastic)))
})

test_that("RPIV_test yields different test_statistic and T_null if gamma is large", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)
  result <- RPIV_test(Y, X, Z = Z, gamma = 1)
  expect_false(result$test_statistic == result$T_null)
})

test_that("RPIV_test yields same test_statistic and T_null if gamma is small", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)
  result <- RPIV_test(Y, X, Z = Z, gamma = 0.00001)
  expect_equal(result$test_statistic, result$T_null)
})

test_that("RPIV_test works with clustering", {
  set.seed(1)
  n <- 20
  cluster_ids <- rep(1:10, each = 2)
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X - C + rnorm(n)
  result <- RPIV_test(Y, X, C, Z, clustering = cluster_ids, variance_estimator = "cluster")
  expect_type(result, "list")
  expect_equal(result$variance_estimator, "cluster")
  cluster_ids <- as.factor(cluster_ids)
  result <- RPIV_test(Y, X, C, Z, clustering = cluster_ids, variance_estimator = "cluster")
  expect_type(result, "list")
  expect_equal(result$variance_estimator, "cluster")
  cluster_ids <- as.character(cluster_ids)
  result <- RPIV_test(Y, X, C, Z, clustering = cluster_ids, variance_estimator = "cluster")
  expect_type(result, "list")
  expect_equal(result$variance_estimator, "cluster")
})



test_that("RPIV_test works with various dimensional input and throws error in underidentified case", {
  set.seed(1)
  n <- 20
  Z1 <- rnorm(n)
  Z2 <- cbind(Z1, rnorm(n))
  Z3 <- cbind(Z2, rnorm(n))
  C <- cbind(rnorm(n), rnorm(n))
  X1 <- rnorm(n)
  X2 <- cbind(X1, rnorm(n))
  X3 <- cbind(X2, rnorm(n))
  Y <- rnorm(n)
  expect_error(RPIV_test(Y, X2, C, Z1))
  expect_error(RPIV_test(Y, X3, C, Z2))
  expect_no_error(RPIV_test(Y, X1, C, Z2))
  expect_no_error(RPIV_test(Y, X2, C, Z2))
  expect_no_error(RPIV_test(Y, X2, C, Z3))
  expect_no_error(RPIV_test(Y, X3, C, Z3))
})

test_that("RPIV_test throws errors on bad input", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- rnorm(n)
  Y <- rnorm(n)

  expect_error(RPIV_test(Y[1:19], X = X, Z = Z))
  expect_error(RPIV_test(Y, X, Z = matrix(Z, ncol = 1), frac_A = -0.1), "frac_A must be NULL or a numeric scalar in")
  expect_error(RPIV_test(Y, X, Z = Z, regr_par = list(num_trees = 5)), "regr_par contains invalid entries")
  expect_error(RPIV_test(Y, X, Z = Z, variance_estimator = "homoscedastic"), "variance_estimator must be one or more of")
})

test_that("RPIV_test works with fit_intercept = FALSE", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- rnorm(n)
  Y <- rnorm(n)

  result <- RPIV_test(Y, X, C, Z, fit_intercept = FALSE)

  expect_type(result, "list")
  expect_equal(result$variance_estimator, "heteroskedastic")

  expect_no_error(RPIV_test(Y, X, C = NULL, Z, fit_intercept = FALSE))
})


test_that("RPIV_test works with dataframes and vectors with attributes", {
  set.seed(1)
  n <- 20
  Z <- data.frame(cbind(rnorm(n), rnorm(n)))
  C <- data.frame(cbind(rnorm(n), rnorm(n)))
  X <- data.frame(cbind(rnorm(n), rnorm(n)))
  Y <- rnorm(n)
  attr(Y, "label") <- "testlabel"
  result <- RPIV_test(Y, X, C, Z)
  expect_type(result, "list")
  expect_equal(result$variance_estimator, "heteroskedastic")
})

test_that("RPIV_test raises error if Z is not convertible to vector", {
  set.seed(1)
  n <- 20
  Z <- list(a = c("a", "aa"), b = rnorm(n))
  C <- rnorm(n)
  X <- rnorm(n)
  Y <- rnorm(n)
  expect_error(RPIV_test(Y, X, C, Z))
})




test_that("weak_RPIV_test returns a function with correct output structure", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X - C + rnorm(n)

  test_fun <- weak_RPIV_test(Y, X, C, Z)

  expect_type(test_fun, "closure")

  result <- test_fun(beta = 1, type = "tune_and_fit")

  expect_type(result, "double")
  expect_named(result, "heteroskedastic")
  expect_true(is.numeric(result["heteroskedastic"]))
})

test_that("weak_RPIV_test works with C = NULL and/or fit_intercept = FALSE", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)

  test_fun1 <- weak_RPIV_test(Y = Y, X = X, Z = Z)
  test_fun2 <- weak_RPIV_test(Y = Y, X = X, C = C, Z = Z, fit_intercept = FALSE)
  test_fun3 <- weak_RPIV_test(Y = Y, X = X, Z = Z, fit_intercept = FALSE)


  expect_no_error(test_fun1(beta = 0, type = "tune_and_fit"))
  expect_no_error(test_fun2(beta = 0, type = "tune_and_fit"))
  expect_no_error(test_fun3(beta = 0, type = "tune_and_fit"))
})


test_that("weak_RPIV_test handles multiple variance estimators", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)

  test_fun <- weak_RPIV_test(
    Y, X, Z = Z,
    variance_estimator = c("homoskedastic", "heteroskedastic")
  )

  result <- test_fun(beta = 0, type = "tune_and_fit")

  expect_type(result, "double")
  expect_named(result, c("homoskedastic", "heteroskedastic"))
})


test_that("weak_RPIV_test works with clustering", {
  set.seed(1)
  n <- 20
  cluster_ids <- rep(1:10, each = 2)
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X - C + rnorm(n)

  test_fun <- weak_RPIV_test(
    Y, X, C, Z,
    clustering = cluster_ids,
    variance_estimator = "cluster"
  )
  result <- test_fun(beta = 1, type = "tune_and_fit")
  expect_named(result, "cluster")
})

test_that("weak_RPIV_test handles various dimensions and TSLS requirements", {
  set.seed(1)
  n <- 20
  Z1 <- rnorm(n)
  Z2 <- cbind(Z1, rnorm(n))
  Z3 <- cbind(Z2, rnorm(n))
  C <- cbind(rnorm(n), rnorm(n))
  X1 <- rnorm(n)
  X2 <- cbind(X1, rnorm(n))
  X3 <- cbind(X2, rnorm(n))
  Y <- rnorm(n)

  expect_error(weak_RPIV_test(Y, X2, C, Z1, fit_at_tsls = TRUE))
  expect_error(weak_RPIV_test(Y, X3, C, Z2, fit_at_tsls = TRUE))

  fun1 <- weak_RPIV_test(Y, X2, C, Z1, fit_at_tsls = FALSE)
  fun2 <- weak_RPIV_test(Y, X3, C, Z2, fit_at_tsls = FALSE)

  fun3 <- weak_RPIV_test(Y, X2, C, Z2)
  fun4 <- weak_RPIV_test(Y, X3, C, Z3)

  expect_error(fun1(1, "tune_and_fit"))
  expect_error(fun1(c(1,2), type = "fit"))
  expect_no_error(fun1(c(1,2), "tune_and_fit"))
  expect_no_error(fun2(c(1,2,3), "tune_and_fit"))

  expect_no_error(fun3(c(1,2), "fit"))
  expect_no_error(fun4(c(1,2,3), "recalculate"))
})


test_that("weak_RPIV_test throws errors on bad input", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- rnorm(n)
  Y <- rnorm(n)

  expect_error(weak_RPIV_test(Y[1:19], X = X, Z = Z))
  expect_error(
    weak_RPIV_test(Y, X, Z = Z, frac_A = -0.1),
    "frac_A must be NULL or a numeric scalar"
  )
  expect_error(
    weak_RPIV_test(Y, X, Z = Z, regr_par = list(num_trees = 5)),
    "regr_par contains invalid entries"
  )
  expect_error(
    weak_RPIV_test(Y, X, Z = Z, variance_estimator = "homoscedastic"),
    "variance_estimator must be one or more of"
  )
})

test_that("weak_RPIV_test works with dataframes and vectors with attributes", {
  set.seed(1)
  n <- 20
  Z <- data.frame(cbind(rnorm(n), rnorm(n)))
  C <- data.frame(cbind(rnorm(n), rnorm(n)))
  X <- data.frame(cbind(rnorm(n), rnorm(n)))
  Y <- rnorm(n)
  attr(Y, "label") <- "testlabel"

  test_fun <- weak_RPIV_test(Y, X, C, Z)
  result <- test_fun(beta = c(0, 0), type = "tune_and_fit")

  expect_type(result, "double")
})

test_that("weak_RPIV_test raises error if Z is not convertible to matrix", {
  set.seed(1)
  n <- 20
  Z <- list(a = c("a", "aa"), b = rnorm(n))
  C <- rnorm(n)
  X <- rnorm(n)
  Y <- rnorm(n)

  expect_error(weak_RPIV_test(Y, X, C, Z))
})

test_that("weak_RPIV_test rejects invalid type", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)

  test_fun <- weak_RPIV_test(Y, X, Z = Z)

  expect_error(
    test_fun(beta = 0, type = "invalid_type"),
    "type needs to be one of"
  )
})


test_that("weak_RPIV_test rejects invalid type", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X + rnorm(n)

  test_fun <- weak_RPIV_test(Y, X, Z = Z)

  expect_error(
    test_fun(beta = 0, type = "invalid_type"),
    "type needs to be one of"
  )
})


test_that("test_fun returned weak_RPIV_test is deterministic when type = 'recalculate'", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- rnorm(n)

  test_fun <- weak_RPIV_test(Y, X, Z = Z)

  T1 <- test_fun(beta = 0, type = "recalculate")
  T2 <- test_fun(beta = 0, type = "recalculate")

  expect_true(abs(T1 - T2) < 10^(-10))
})
