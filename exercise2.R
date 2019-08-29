# III.1 Metropolis-Hastings

# Model

# target
p = function(x) {
  sin(x)^2*exp(-abs(x))
}

# proposal
sigma = 0.1    # IF SIGMA IS TOO SMALL IT CANNOT JUMP BETWEEN THE TWO MODELS
rq = function(x_prime) {
  rnorm(n = 1, mean = x_prime, sd = sigma)
}

dq = function(x, x_prime) {
  dnorm(x = x, mean = x_prime, sd = sigma)
}


# Inference

M = 10000
x = numeric(M)

# Initialization
x[1] = runif(0, 1)

for (i in 1:M) {
  x_prime = rq(x[i])
  u       = runif(n = 1, min = 0, max = 1)
  alpha   = min(1, p(x_prime)         /  p(x[i])  * 
                   dq(x[i], x_prime)  / dq(x_prime, x[i])
               )
                   
                
  if (u <= alpha) { 
    x[i + 1] = x_prime
  } else
    x[i + 1] = x[i]
}

burnin = 1000
hist(x[burnin:M], breaks = 100, probability = TRUE, xlim = c(-10, 10))
curve(p, from = -10, to = 10, add = TRUE)


plot(p(seq(from = -10, to = 10, by = 0.01)))


# III.2 Gibbs sampling

mu    = c(7, 3)
sigma = matrix(
  c(0.3, 0.1, 
     0.1, 1),
 2)

muab = function(mu_a, mu_b, x_b) {
  mua + sigma[1, 2]/sigma(2, 2)*(x_b - mu_b)
}

sigmaab = function(sigma_a, sigma_b) {
  
}

p = function(c) { # samples condtional distribution given component c
  #mvrnorm(mu = mu, Sigma = sigma ) # from MASS package
  rnorm(n = 1, mu = muab(c), )
}

# Inference

x = list()
x[[1]] = c(0, 0)
