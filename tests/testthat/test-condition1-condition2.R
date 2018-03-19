context("Helper - Condition parsing for difference estimators")


test_that("Condition arguments behave as expected", {
  n <- 40
  dat <- data.frame(
    y = rnorm(n),
    bl = rep(1:5, each = 8),
    z = 1:4,
    ps = runif(n)
  )

  horvitz_thompson(y ~ z, data = dat, subset = z <= 2, condition_prs = ps)

  # Subsetting and just selecting two conditions
  expect_identical(
    horvitz_thompson(y ~ z, data = dat, subset = z <= 2, condition_prs = ps),
    horvitz_thompson(y ~ z, data = dat, condition1 = 1, condition2 = 2, condition_prs = ps)
  )

  expect_identical(
    difference_in_means(y ~ z, data = dat, subset = z <= 2),
    difference_in_means(y ~ z, data = dat, condition1 = 1, condition2 = 2)
  )

  expect_identical(
    difference_in_means(y ~ z, data = dat, subset = z <= 2, blocks = bl),
    difference_in_means(y ~ z, data = dat, condition1 = 1, condition2 = 2, blocks = bl)
  )

  # Subsetting and just selecting two conditions
  expect_identical(
    tidy(horvitz_thompson(
      y ~ z,
      data = dat,
      condition1 = 3,
      condition2 = 4,
      condition_prs = rep(0.5, nrow(dat))
    ))[c("estimate", "std.error")],
    tidy(horvitz_thompson(
      y ~ z,
      data = dat,
      condition1 = 4,
      condition2 = 3,
      condition_prs = rep(0.5, nrow(dat))
    ))[c("estimate", "std.error")] * c(-1, 1)
  )

  expect_identical(
    tidy(difference_in_means(
      y ~ z,
      data = dat,
      condition1 = 2,
      condition2 = 1
    ))[c("estimate", "std.error")],
    tidy(difference_in_means(
      y ~ z,
      data = dat,
      condition1 = 1,
      condition2 = 2
    ))[c("estimate", "std.error")] * c(-1, 1)
  )

  # Errors if not specifying both
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition1 = 4,
      condition_prs = ps
    ),
    "condition1"
  )
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition2 = 4,
      condition_prs = ps
    ),
    "condition1"
  )
  expect_error(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition_prs = ps
    ),
    "condition1"
  )

  expect_error(
    difference_in_means(
      y ~ z,
      data = dat,
      condition1 = 4
    ),
    "condition1"
  )
  expect_error(
    difference_in_means(
      y ~ z,
      data = dat,
      condition2 = 4
    ),
    "condition1"
  )
  expect_error(
    difference_in_means(
      y ~ z,
      data = dat
    ),
    "condition1"
  )


  # Specifying only one works with binary treatment
  dat$z <- c("Treated", "Control")
  expect_identical(
    difference_in_means(
      y ~ z,
      data = dat,
      condition1 = "Treated"
    ),
    difference_in_means(
      y ~ z,
      data = dat,
      condition1 = "Treated",
      condition2 = "Control"
    )
  )

  expect_identical(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition1 = "Treated"
    ),
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition1 = "Treated",
      condition2 = "Control"
    )
  )
  expect_identical(
    difference_in_means(
      y ~ z,
      data = dat,
      condition2 = "Treated"
    ),
    difference_in_means(
      y ~ z,
      data = dat,
      condition2 = "Treated",
      condition1 = "Control"
    )
  )

  expect_identical(
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition2 = "Treated"
    ),
    horvitz_thompson(
      y ~ z,
      data = dat,
      condition2 = "Treated",
      condition1 = "Control"
    )
  )

  # Works with factor
  dat$z <- factor(c("T", "C"))
  # Must pass string!
  difference_in_means(y ~ z, condition2 = "T", data = dat)
  # Errors if not found
  expect_error(
    difference_in_means(
      y ~ z,
      condition2 = 1,
      data = dat
    ),
    "`condition1` and `condition2` must be values found in the treatment"
  )

  dat$z <- 1
  expect_error(
    difference_in_means(y ~ z, data = dat),
    "Must have more than one value in treatment unless using Horvitz"
  )
})
