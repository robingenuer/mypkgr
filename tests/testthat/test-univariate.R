context("Univariate gaussian test")

test_that("correct result for univariate gaussian", {
  expect_equal(mvnpdf(x=matrix(1.96), mean=0, varcovM=diag(1), Log=FALSE)$y,
               dnorm(1.96))
  expect_equal(mvnpdf(x=matrix(c(1.96, -0.5), ncol = 2), mean=0, varcovM=diag(1), Log=FALSE)$y,
               dnorm(c(1.96, -0.5)))
})
