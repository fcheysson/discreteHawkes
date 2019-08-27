x <- hawkes(100, 1, 1, 2, 3.5)

test_that("intensity in 0", {
  expect_equal(intensity(x, 0), 3.5)
})

x <- hawkes(100, 1, 1, 2)

test_that("intensity in 0", {
  expect_equal(intensity(x, 0), 1)
})

test_that("compensator in T", {
  expect_lt(compensator(x, 100), 100+x$n/2)
})

y <- discrete(x, length = 100)