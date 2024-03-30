# added to .Rbuildignore to keep CRAN from complaining that broom
# isn't a SUGGESTed package
context("zzzbroom.R - .onAttach")

test_that(".onLoad message if old version of 'broom' is installed", {
  skip_if_not_installed("broom")
  skip_if(packageVersion("broom") > "0.5.0")
  library(broom)
  expect_message(
    .onLoad("estimatr", "estimatr"),
    "the `broom` package version 0.5.0 or earlier is loaded"
  )
})
