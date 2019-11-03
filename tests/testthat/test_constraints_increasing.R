context("Monotone increasing constraints work.")

set.seed(1337)
data = stats::rbeta(1000, 11, 1)*6 - 3
s = quantile(data, (1:4)/5)
support = c(-3,3)

expect_equal(discrepancy(polygram(data ~ increasing, s = 0, support = support)),
             -3.952639,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ increasing, s = 1, support = support)),
             -5.434724,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ increasing, s = 2, support = support)),
             -5.671504,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ increasing, s = 3, support = support)),
             -5.701314,
             tolerance = 1e-6)
expect_equal(discrepancy(polygram(data ~ increasing, s = s, support = support)),
             -5.58013,
             tolerance = 1e-6)

