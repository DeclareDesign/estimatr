
context("na.omit_detailed")

df <- expand.grid(Y = c(1:5, NA), Z = c(LETTERS, NA))

stock <- na.omit(df)
detailed <- na.omit_detailed.data.frame(df)

stock_action <- attr(stock, "na.action")
detailed_action <- attr(detailed, "na.action")

test_that("Omits are the same", {
  expect_equal(
    as.vector(stock_action),
    as.vector(detailed_action)
  )
})

test_that("Row names are set correctly", {
  expect_identical(
    names(stock_action),
    names(detailed_action)
  )
})

