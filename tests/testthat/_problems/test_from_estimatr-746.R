# Extracted from test_from_estimatr.R:746

# prequel ----------------------------------------------------------------------
library(estimatrZero)
rmcall <- function(x) {
  x$call <- NULL
  if (!is.null(x$terms)) attr(x$terms, ".Environment") <- NULL
  x
}
se_types    <- c("classical", "HC0", "HC1", "HC2", "HC3")
cr_se_types <- c("CR0", "stata", "CR2")

# test -------------------------------------------------------------------------
set.seed(43)
N <- 40
dat <- data.frame(
    Y  = rnorm(N), Y2 = rnorm(N),
    Z  = rbinom(N, 1, .5), X  = rnorm(N),
    B  = factor(rep(1:4, each = 10))
  )
ro  <- lm_robust(cbind(Y, Y2) ~ Z + X + factor(B), data = dat)
rfo <- lm_robust(cbind(Y, Y2) ~ Z + X, fixed_effects = ~B, data = dat)
