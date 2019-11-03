context("utility_functions")

set.seed(313)
A = matrix(sample(2, 9, replace = TRUE), nrow = 2)
B = matrix(sample(2, 9, replace = TRUE), nrow = 3)
C = matrix(sample(2, 9, replace = TRUE), nrow = 4)

dsum = direct_sum(A, B, C)
expect_equal(dsum[1:2, 1:5], A)
expect_equal(dsum[3:5, 6:8], B)
expect_equal(dsum[6:9, 9:11], C)
expect_equal(dsum[4, 4], 0)
