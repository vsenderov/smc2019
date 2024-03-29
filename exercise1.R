# I.1
# b)
N <- 10000
props <- rnorm(N, 2,2)

w <- (0<=props&props<=4)*1/4/dnorm(props,2,2)

library(plotrix)
weighted.hist(props,w, breaks = 30)

# c)
B <- 10000
N <- 10
means <- numeric(B)
medians <- numeric(B)
for(i in 1:B){
  props <- rnorm(N, 0,1)
  w <- (0<=props&props<=4)*1/4/dnorm(props, 0,1)
  # means[i] <- weighted.mean(props, w)
  means[i] <- mean(props*w)
  medians[i] <- median(props*w)
}
hist(means, breaks = 30)
hist(medians, breaks = 30)
mean(means, na.rm = TRUE)
mean(medians,na.rm = TRUE)

# e)
N <- 100
B <- 10000
means <- numeric(B)
for(i in 1:B){
  props <- runif(N, 0,4)
  w <- (0<=props&props<=4)/dunif(props, 0,4)
  means[i] <- mean(w)
}
hist(means)
mean(means)
var(means)

# f)
B <- 10000
N <- 1000
means <- numeric(B)
for(i in 1:B){
  props <- rnorm(N, 0, sd=1)
  w <- (0<=props&props<=4)/dnorm(props, 0,sd = 1)
  # means[i] <- weighted.mean(props, w)
  integral <-sum(props*w)
  z_hat <- sum(w)
  means[i] <- integral/z_hat
}
hist(means, breaks = 30)
mean(means, na.rm = TRUE)

# I.2
library(mvtnorm)
D_max <- 15
N <- 10000

prop_in_cube <- numeric(D_max-1)

for(D in 2:D_max){
  props <- rmvnorm(N, mean = rep(0,D),sigma = diag(rep(1, D)))
  pi_tilde <- apply(props, 1, function(x) all(abs(x)<=0.5))
  #q <- dmvnorm(props, mean = rep(0,D),sigma = diag(rep(1, D)))
  #w <- pi_tilde/q
  prop_in_cube[D-1] <-mean(pi_tilde)
}
  
plot(2:D_max, prop_in_cube, type = "l")

# I.3
# q)
D <- 1000
N <- 10
# without logs
props <- rmvnorm(N, mean = rep(0,D),sigma = diag(rep(4, D)))
pi_tilde <- apply(props, 1, function(x) dmvnorm(x, mean = rep(0,D), sigma = diag(rep(1, D))))
q <- apply(props, 1, function(x) dmvnorm(x, mean = rep(0,D), sigma = diag(rep(4, D))))
w <- pi_tilde/q
w_normalized <- w/sum(w)
w
w_normalized
sum(w_normalized)

# with logs
props <- rmvnorm(N, mean = rep(0,D),sigma = diag(rep(4, D)))
pi_tilde_log <- apply(props, 1, function(x) sum(dnorm(x, mean = 0, sd = 1, log = TRUE)))
q_log <- apply(props, 1, function(x) sum(dnorm(x, mean = 0, sd = 2, log = TRUE)))
w_log <- pi_tilde_log - q_log
w_log
exp(w_log)!=0

# normalize
w_bar <- exp(w_log-max(w_log))
w_bar 
w_normalized <- w_bar/sum(w_bar)
w_normalized
  

# I.4

y <- read.csv("seOMXlogreturns2012to2014.csv", header = FALSE)[,1]
theta <- c(0.98, 0.16, 0.7)
N <- 500
totalT <- 500
mean_states <- numeric(500)
mean_pred <- numeric(500)

# Initialize
x_old <- rnorm(N, mean = 0, sd = theta[2]*5)
w_log <- log(rep(1/N, N))
# Iterate through time
for(t in 1:totalT){
  # resample
  a <- sample(1:N, N, prob = exp(w_log - max(w_log)), replace = TRUE)
  # propagate
  x_new <- rnorm(N, mean = theta[1]*x_old[a], sd = theta[2])
  # Weight
  mean_pred[t] <- mean(x_new)
  w_log <- dnorm(y[t], mean = 0, sd = theta[3]*sqrt(exp(x_new)), log = TRUE)
  mean_states[t] <- weighted.mean(x_new, exp(w_log - max(w_log)))
  x_old <- x_new
}




plot(1:500, mean_states, type = "l", ylim = range(y))
lines(1:500, y, col = "blue")
lines(1:500, abs(y), col = "red")

plot(1:500, mean_pred, type = "l", ylim = range(y))
lines(1:500, mean_states, col = "blue")
mean_pred[1:499]-mean_states[2:500]



# 2.1 

y <- read.csv("seOMXlogreturns2012to2014.csv", header = FALSE)[,1]
theta <- c("phi" = 0.16,"sigma" =  0.7)
beta <-  seq(from = 0.1, to = 2, by = 0.1)


B <- 10 # number of particle filter runs
log_likelihood <- matrix(rep(0, length(beta)*B), nrow = length(beta), ncol = B)
for (j in 1:B)
{
  i = 1
  for (b in beta) {
    N <- 10
    totalT <- 500
    mean_states <- numeric(500)
    mean_pred <- numeric(500)
    
    # Initialize
    x_old <- rnorm(N, mean = 0, sd = theta[2]*5)
    w_log <- log(rep(1/N, N))
    # Iterate through time
    for(t in 1:totalT){
      # resample
      a <- sample(1:N, N, prob = exp(w_log - max(w_log)), replace = TRUE)
      
      # propagate
      x_new <- rnorm(N, mean = theta[1]*x_old[a], sd = theta[2])
      # Weight
      mean_pred[t] <- mean(x_new)
      w_log <- dnorm(y[t], mean = 0, sd = b*sqrt(exp(x_new)), log = TRUE)
      mean_states[t] <- weighted.mean(x_new, exp(w_log - max(w_log)))
      x_old <- x_new
      
      log_likelihood[i,j] = log_likelihood[i,j] + log(sum(exp((w_log)))) - log(N)
      
    }
    i = i + 1
  }  
}


boxplot(t(log_likelihood))



# 2.c

y <- read.csv("seOMXlogreturns2012to2014.csv", header = FALSE)[,1]
theta <- c("phi" = 0.16,"sigma" =  0.7)
beta <-  seq(from = 0.1, to = 2, by = 0.1)


B <- 10 # number of particle filter runs
log_likelihood <- matrix(rep(0, length(beta)*B), nrow = length(beta), ncol = B)
for (j in 1:B)
{
  i = 1
  for (b in beta) {
    N <- 200
    totalT <- 500
      mean_states <- numeric(500)
    mean_pred <- numeric(500)
    
    # Initialize
    x_old <- rnorm(N, mean = 0, sd = theta[2]*5)
    w_log <- log(rep(1/N, N))
    
    w_log_old <- w_log
    # Iterate through time
    for(t in 1:totalT){
      # resample
      #a <- sample(1:N, N, prob = exp(w_log - max(w_log)), replace = TRUE)
      a = 1:N
      # propagate
      x_new <- rnorm(N, mean = theta[1]*x_old[a], sd = theta[2])
      # Weight
      mean_pred[t] <- mean(x_new)
      w_log <- dnorm(y[t], mean = 0, sd = b*sqrt(exp(x_new)), log = TRUE) + w_log_old
      
      w_log_old <- w_log
      
      mean_states[t] <- weighted.mean(x_new, exp(w_log - max(w_log)))
      x_old <- x_new
      
      log_likelihood[i,j] = log_likelihood[i,j] + log(sum(exp((w_log)))) - log(N)
      
    }
    i = i + 1
  }  
}


boxplot( t(log_likelihood)  )




# II.2.b

# Model II Set-up
f = function(x_previous) {
  cos((x_previous))^2
}

C = 2

Q = 1

R = .01

# Guassian model with non-linear dynamics

proposal = function(x_previous, y) {
  K = Q * C * ( C * Q * C + R)^-1
  Sigma = sqrt( (1 - K * C) * Q ) 
  rnorm(n = N, mean = f(x_previous) + K*(y - C*f(x_previous)), Sigma)
}

weights = function(y,x_previous) {
  dnorm(y, mean = C * f(x_previous), sd = sqrt(C * Q * t(C) + R ))
}

# Fully Adapted Particle Filter - a variant of auxillary particle filter

# Initialization
T = 20
N = 10          # particles
x_true = numeric(T)
x = list()
y = list()
nu = list()
x_true[1] = 0 # We know the true starting state
x[[1]] = rep(x_true[1], N)

for (i in 2:T) {
  # compute new x and y
  x_true[i] = f(x_true[i-1]) + rnorm(n = 1, mean = 0, sd = sqrt(Q))
  
  y[[i]] = C * x_true[i] + rnorm(n = N, mean = 0, sd = sqrt(R))
  
  # resampling, vector of weights nu, weights = 1/n
  nu[[i-1]] = 1/N * weights(y[[i]], x[[i-1]])
  # ancestor indicies
  a = sample(x = 1:N, size = N, replace = TRUE, prob = nu[[i-1]])
  # propagation
  x[[i]] = proposal( x[[i-1]][a], y[[i]])
  
}


# Bootstrap particle filter
w = numeric(N)
x_bootstrap = list()
x_bootstrap[[1]] = rep(x_true[1], N) # we knnow the true starting state
w = rep(1/N, N) # iitial weights
x_weighted_mean = numeric(N)

for (i in 2:T) {
  a = sample(x = 1:N, size = N, replace = TRUE, prob = w) # ancestor indicesa
  x_bootstrap[[i]] = f(x_bootstrap[[i-1]][a]) + rnorm(n = N, mean = 0, sd = sqrt(Q)) # propagation step
  
  # y[[i]] = C * x_true[i] + rnorm(n = N, mean = 0, sd = sqrt(R))
  # y[[i]] has been computed previously
  
  w_unnormalized = dnorm(x = y[[i]], mean = C*x_bootstrap[[i]], sd = sqrt(R) )
  
  w = w_unnormalized / sum(w_unnormalized)
  
  x_weighted_mean[i] = sum(x_bootstrap[[i]]*w)
}

plot(x_true)
#lines(x_true, col = "red")


points(sapply(X = x, mean), col = "red")
lines(sapply(X = x, mean), col = "red")

points(x_weighted_mean, col = "green")
lines(x_weighted_mean, col = "green")


