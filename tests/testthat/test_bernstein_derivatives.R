context("bernstein derivatives")

lambda = c(0.4, 0.5, 0.1)
m = length(lambda) - 1

expect_equal(numDeriv::grad(function(x) dbernstein(x, lambda), x = 1),
             (-lambda[m] + lambda[m + 1])*(m + 1) * m, tolerance = 1e-8)

expect_equal(c(numDeriv::hessian(function(x) dbernstein(x, lambda), x = 1)),
             (lambda[m - 1] - 2 * lambda[m] + lambda[m + 1])*(m + 1)*m,
             tolerance = 1e-8)
expect_equal(c(numDeriv::hessian(function(x) dbernstein(x, lambda), x = 0)),
            (lambda[3] - 2*lambda[2] + lambda[1])*(m+1)*m,
            tolerance = 1e-3)
