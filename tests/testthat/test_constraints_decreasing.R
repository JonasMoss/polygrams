context("Monotone decreasing constraints work.")
set.seed(1337)
data = stats::rbeta(1000, 1, 9)
s = quantile(data, (1:4)/5)

expect_equal(discrepancy(polygram(data ~ decreasing, s = 0)),
             -3.669324,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ decreasing, s = 1)),
             -4.660745,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ decreasing, s = 2)),
             -4.715786,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ decreasing, s = 3)),
             -4.727753,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ decreasing, s = s)),
             -4.65472,
             tolerance = 1e-6)
