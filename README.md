
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RPIV

<!-- badges: start -->
<!-- badges: end -->

The RPIV package implements a residual prediction test for the well
specification of linear instrumental variable (IV) models, as presented
in Scheidegger, Londschien and Bühlmann (2025). For a response
$Y_i\in \mathbb R$, endogenous explanatory variables
$X_i\in \mathbb R^p$ and instruments $Z_i\in \mathbb R^d$ ($d\geq p$),
it tests the well-specification of the linear IV model (\`\`is the
linear IV model appropriate?’’). More formally, it tests the null
hypothesis
$$H_0: \exists \beta\in \mathbb R^p\text{ s.t. } \mathbb E[Y_i - X_i^T\beta|Z_i] = 0 \text{ a.s., } i =1,\ldots, n,$$
which is implied by the well-specification of the linear IV model (with
mean-independence assumption on the errors).

The model allows for additional exogenous explanatory variables
(\`\`exogenous controls’’) (denoted by $C$ in the R function) and an
intercept, which are added both to $X$ and $Z$.

For a detailed discussion of the method, we refer to Scheidegger,
Londschien and Bühlmann (2025). We now demonstrate, how the RPIV package
is used in practice.

## Installation

You can install the development version of IVDML from
[GitHub](https://github.com/) with

``` r
devtools::install_github("cyrillsch/RPIV")
#> Skipping install of 'RPIV' from a github remote, the SHA1 (cb085085) has not changed since last install.
#>   Use `force = TRUE` to force installation
```

## Example

This is a basic example presenting, how the well-specification of linear
IV models can be tested with the RPIV package. We simulate a dataset
with $n = 200$ observation and three responses.

``` r
set.seed(1)
n <- 200
C <- rnorm(n) # exogenous explanatory variable
Z <- cbind(rnorm(n), C + rnorm(n)) # instrumental variable
H <- rnorm(n) # hidden confounding
X <- Z[, 1] - Z[, 2] + rnorm(n) # endogenous explanatory variable
Y1 <- X - C + H + rnorm(n) # linear IV model
Y2 <- X - C + H + Z[, 1]^2 + rnorm(n) # invalid IV -> misspecified
Y3 <- 2 * sign(X - C) + H + rnorm(n) # nonlinear IV model -> misspecified
```

To apply the well-specification test to the three responses, we use the
function , which uses a heteroskedasticity robust variance estimator by
default.

``` r
library(RPIV)
result1 <- RPIV_test(Y = Y1, X = X, C = C, Z = Z)
result2 <- RPIV_test(Y = Y2, X = X, C = C, Z = Z)
result3 <- RPIV_test(Y = Y3, X = X, C = C, Z = Z)

result1$p_value
#> [1] 0.1575286
result2$p_value
#> [1] 0.0004228503
result3$p_value
#> [1] 0.005525054
```

We see that, indeed, well-specification is rejected at significance
level $\alpha = 0.05$ for the responses $Y_2$ and $Y_3$.

The RPIV package also supports cluster-robust inference. We simulate
data with 50 clusters of size 4, but the linear IV model is
well-specified otherwise.

``` r
set.seed(1)
n <- 200
clustering <- rep(1:50, length.out = n)
Z <- rep(rnorm(1:50), length.out = n) + 0.1 * rnorm(n)
H <- rep(rnorm(1:50), length.out = n) + 0.1 * rnorm(n)
X <- Z + rep(rnorm(1:50), length.out = n) + 0.1 * rnorm(n)
Y <- X + H + rep(rnorm(1:50), length.out = n) + 0.1 * rnorm(n)
```

We apply the test with three different variance estimators: assuming
homoskedasticity, robust to heteroskedasticity, robust to clustering.

``` r
result <- RPIV_test(Y = Y, X = X, C = NULL, Z = Z, variance_estimator = 
                      c("homoskedastic", "heteroskedastic", "cluster"), clustering = clustering)
result$homoskedastic$p_value
#> [1] 0.02844595
result$heteroskedastic$p_value
#> [1] 0.01728716
result$cluster$p_value
#> [1] 0.1347029
```

We see that only using the cluster-robust variance estimator does not
reject the null hypothesis at significance level $\alpha = 0.05$.

## More Examples

More examples can be found in Scheidegger, Londschien and Bühlmann
(2025) and the associated GithHub repository ????.

## References

Cyrill Scheidegger, Malte Londschien and Peter Bühlmann. A residual
prediction test for the well-specification of linear instrumental
variable models. Preprint, arXiv:????, 2025.
