## ============================================================================
## This file contains tests for constraints. There are several important "edge"
## cases:
##
## 1.) s = NULL and s not equispaced.
## 2.) ps != (ms-1), for instance p = 0.
## 3.) different supports than c(0, 1).
##
test_that("Monotone (decreasing) constraints work.", {
  set.seed(1337)
  data = stats::rbeta(1000, 1, 9)
  expect_equal(discrepancy(polygram(data, s = 0,
                       monotone = "decreasing")), -3.669324, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = 1,
                       monotone = "decreasing")), -4.660745, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = 2,
                       monotone = "decreasing")), -4.715786, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = 3,
                       monotone = "decreasing")), -4.727753, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = quantile(data, (1:4)/5),
                       monotone = "decreasing")), -4.65472 , tolerance = 1e-6)
})

test_that("Monotone (increasing) constraints work.", {
  set.seed(1337)
  data = stats::rbeta(1000, 11, 1)*6 - 3
  expect_equal(discrepancy(polygram(data, s = 0, support = c(-3,3),
                                    monotone = "increasing")), -3.952639, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = 1, support = c(-3,3),
                                    monotone = "increasing")), -5.434724, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = 2, support = c(-3,3),
                                    monotone = "increasing")), -5.671504, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = 3, support = c(-3,3),
                                    monotone = "increasing")), -5.701314, tolerance = 1e-6)
  expect_equal(discrepancy(polygram(data, s = quantile(data, (1:4)/5), support = c(-3,3),
                                    monotone = "increasing")), -5.58013 , tolerance = 1e-4)
})

polygram(data, s = quantile(data, (1:4)/5), support = c(-3,3), monotone = "increasing", m = 0, p = -1) %>% plot


