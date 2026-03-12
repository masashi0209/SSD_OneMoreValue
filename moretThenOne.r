# 5th percentile of Standard Normal Distribution
K005 <- qnorm(0.05) 
# sample size
SS <- 5

# function definition
one_trial_delta <- function(mu, sigma, SS = 5, K005 = qnorm(0.05)) {
  # generate an original dataset
  preData <- rnorm(SS, mean = mu, sd = sigma)
  preM   <- mean(preData)
  preSD  <- sd(preData)
  preHC5 <- preM + preSD * K005
  # add two more
  addData  <- rnorm(2, mean = mu, sd = sigma)
  # if we wish to consider more than 2 values (n in general), replace this function to rnorm(n, mean = mu, sd = sigma)
  postData <- c(preData, addData)
  postM    <- mean(postData)
  postSD   <- sd(postData)
  postHC5  <- postM + postSD * K005
  # amount of change 
  deltaLogHC5 <- (postHC5 - preHC5) / preSD
  return(deltaLogHC5)
}

# set population parameters
mu    <- 3.141592      # pi
sigma <- 0.602214076   # Avogadro constant

# try a run
delta_example <- one_trial_delta(mu, sigma, SS, K005)
delta_example

# try many runs
n_iter <- 10000
HC5Change1 <- replicate(n_iter, one_trial_delta(mu, sigma, SS, K005))

# show result
hist(HC5Change1,
     breaks = 100, probability = TRUE,
     col = rgb(1.0, 0.725, 0.0, 0.6), border = "gray40",
     main = "Distribution of Δ log(HC5) (two values added simultaneously)",
     xlab = expression(Delta ~ log(HC5)))

# different parameters
mu2    <- 6.6743     # gravity constant
sigma2 <- 6.626      # Planck constant

# try with these different parameters
HC5Change2 <- replicate(n_iter, one_trial_delta(mu2, sigma2, SS, K005))

# show the second results
hist(HC5Change2,
     breaks = 100, probability = TRUE, add = TRUE,
     col = rgb(0.24, 0.60, 0.80, 0.5), border = NA)


legend("topright",
       legend = c("mu=3.141592, sigma=0.6022", "mu=6.6743, sigma=6.626"),
       fill   = c(rgb(1.0,0.725,0.0,0.6), rgb(0.24,0.60,0.80,0.5)),
       border = NA, bty = "n")

# these two are almot identical

# Once the distribution of dLogHC5/sn is obtained, the remaining task is simply to compute the probability that the change of interest occurs. To do so, one only needs to count the number of values of dLogHC5 that satisfy the specified condition.
# EXAMPLE
# now we have a data set of which standard deviation is sn = 2, and we want to know a probability that post-HC5 > 1.5*pre-HC5
# we have a distribution of deltaLogHC5/sn = HC5Change2 -> deltaLogHC5 = sn*HC5Change2 > log10(1.5) -> HC5Change2 > log10(1.5)/sn -> 

th <- log10(1.5)/2
n <- length(HC5Change2)
cnt <- sum(HC5Change2 > th)
p <- cnt / n
p
