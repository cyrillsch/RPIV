test_that("RPIV_test returns correct structure for default arguments", {
  set.seed(1)
  n <- 20
  Z <- rnorm(n)
  C <- rnorm(n)
  X <- Z + rnorm(n)
  Y <- X - C  + rnorm(n)
  result <- RPIV_test(Y, X, C, Z)
  expect_type(result, "list")
  expect_named(result, c("p_value", "test_statistic", "keep_test_statistic", "var_fraction", "variance_estimator"))
  expect_true(is.numeric(result$p_value))
  expect_true(is.numeric(result$test_statistic))
  expect_true(is.logical(result$keep_test_statistic))
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
  expect_true(all(c("p_value", "test_statistic", "keep_test_statistic", "var_fraction", "variance_estimator") %in%
                    names(result$homoskedastic)))
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
