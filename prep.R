                                        # P.1 a

hist(rnorm(n = 1000, mean = 0, sd = sqrt(2)), prob = TRUE)

mynorm = function(mu, sigma)
{
    function(x)
    {
        1/sqrt(2*pi*sigma^2)*exp(-(x - mu)^2/(2*sigma^2))
    }
}

f = mynorm(0, sqrt(2))

curve(expr = f, from = -4, to = 4, add = TRUE)

                                        # b
hist(
    qnorm(
        runif(1000, 0, 1)
    )
)

                                        # c

s = sqrt(10)*qnorm(runif(1000, 0, 1)) + 2

hist(s, prob = TRUE)

var(s)

                                        # d
?set.seed

# P.2
sample_point = function()
{
    c(runif(1, min = -1, max = 1), runif(1, min = -1, max = 1))
}

circle_check = function(p)
{
    sqrt(p[1]^2 + p[2]^2) <= 1
}

N = 10^(1:7)

pi_est = sapply(N, function(n) {
    4*sum(
          sapply(
              1:n,
              function(i) {
                  circle_check(sample_point())
              }
          ))/n
}
)
    
plot(log(N, base = 10), pi_est)
abline(h = pi)
