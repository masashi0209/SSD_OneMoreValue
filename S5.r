# Function definitions
#Standard Deviation
S <- function(x, n) {
  sqrt((n - 1) / n + (x^2) / n)
}

# dLogHC5
DHC5 <- function(x, n, Kp) {
  x / sqrt(n * (n + 1)) + (S(x, n) - 1) * Kp
}


#parameters
# 5th percentile of SND
K005 <- qnorm(0.05, mean = 0, sd = 1)

# sample size, 5 as example
SS <- 5

#two examples with sn = 1, and sn = 0.5

x <- seq(-2, 4, length.out = 501)
y1 <- DHC5(x, SS, K005)
y2 <- 0.5 * y1
plot(
  x, y1, type = "l", lwd = 2, col = "steelblue",
  xlab = "x", ylab = "DHC5", main = "DHC5 and 0.5 * DHC5"
)
lines(x, y2, lwd = 2, col = "tomato")

# DHC5 has maximum at xm, and it is
xm <- sqrt(SS - 1) / sqrt(-1 + K005^2 + K005^2 * SS)

# the value of DHC5 at xm is 
maxDHC5 <- sqrt(SS - 1) / ( sqrt(SS * (1 + SS)) * sqrt(-1 + K005 ^2 * (1 + SS)) ) + K005 * ( -1 + sqrt( (K005^2 * (SS^2 - 1)) / ( SS * (-1 + K005^2 * (1 + SS)) ) ) )

# case 1: we want to know that post-HC5 is greater than pre-HC5
# need to find the intersections of DHC5 and 0
# intersections are derived analytically (do it yourself)
# intersection 1, lower
x1 <- -(K005 * sqrt(SS * (1 + SS)) + sqrt( K005^2 * (1 + SS) * (-1 + SS + K005^2 * (1 + SS)) ))/(-1 + K005^2 * (1 + SS) )
# intersection 2, higher
x2 <- (-K005 * sqrt(SS * (1 + SS)) + sqrt(K005^2 * (1 + SS) * (-1 + SS + K005^2 * (1 + SS))))/(-1 + K005^2 * (1 + SS))

# probability that post-HC5 is higher than pre-HC5
p <- pt(x2, df = SS - 1) - pt(x1, df = SS - 1)


# case 2: we want to know that the post-HC5 is greater than k*pre-HC5
# case 2-1: k<1
# we need to know the intersections of HC5 and log10(1/k) = -log10(k)
# intersections are derived analytically (do it yourself)
# since the intersections are complicated, we define a function
# case 2-1-1: sn = 1
xLower <- function(n, Kp, z) {
  num1 <- Kp * sqrt(n * (1 + n)) * log(10)
  num2 <- sqrt(n * (1 + n)) * log(z)
  # 中央の大きな sqrt の中身
  inner <- Kp^2 * (1 + n) * (
    (-1 + n + Kp^2 * (1 + n)) * log(10)^2 +
      2 * Kp * n * (1 + n) * log(10) * log(z) +
      n * (1 + n) * log(z)^2
  )
  num3 <- sqrt(inner)
  den  <- (-1 + Kp^2 * (1 + n)) * log(10)

  -( (num1 + num2 + num3) / den )
}

xHigher <- function(n, Kp, z) {
  num1 <- -Kp * sqrt(n * (1 + n)) * log(10)
  num2 <- -sqrt(n * (1 + n)) * log(z)
  inner <- Kp^2 * (1 + n) * (
    (-1 + n + Kp^2 * (1 + n)) * log(10)^2 +
      2 * Kp * n * (1 + n) * log(10) * log(z) +
      n * (1 + n) * log(z)^2
  )
  num3 <- sqrt(inner)
  den  <- (-1 + Kp^2 * (1 + n)) * log(10)

  (num1 + num2 + num3) / den
}

# intersection 1, lower, k=1/2 as example
k <- 1/2
x1 <- xLower(SS, K005, k)
x2 <- xHigher(SS, K005, k)

#integrating t-distribution from x1 to x2 gives the probability that post-HC5 > pre-HC5
p <- pt(x2, df = SS - 1) - pt(x1, df = SS - 1)

# case 2-1-2: sn not 1
sn <- 1.5
# convert value of k as 
kc <- k^(1/sn) 
# determine the integration interval and integrate t-distribution
x1 <- xLower(SS, K005, kc)
x2 <- xHigher(SS, K005, kc)
p <- pt(x2, df = SS - 1) - pt(x1, df = SS - 1)
# p is a probability that post-HC5 > k*pre-HC5, and if we want to know the prob. that post-HC5 < k*pre-HC5
pn <- 1 - p

# repeating the process for various sn, we have the left panel in Figure 4
k <- 1 / 2
sn_rep <- 1:18
#
psSn05 <- vector("list", length(sn_rep))
for (i in seq_along(sn_rep)) {
  sn <- 0.3 + 0.1 * (sn_rep[i] - 1)
  kc <- k^(1 / sn)  
  x1 <- xLower(SS, K005, kc)
  x2 <- xHigher(SS, K005, kc)
  p <- pt(x2, df = SS - 1) - pt(x1, df = SS - 1) 
  psSn05[[i]] <- c(s = sn, one_minus_p = 1 - p)
}
# I do not know what they are doing
psSn05 <- do.call(rbind, psSn05)
psSn05 <- as.data.frame(psSn05)
row.names(psSn05) <- NULL

plot(psSn05[, 1], psSn05[, 2],
     type = "l",  
     xlim = c(0.3, 2.0), 
     ylim = c(0, 0.366),  
     xlab = "sn", ylab = "p",
     frame.plot = TRUE  
)

# case 2-2: k > 1
k <- 1.5
# we need to know the min sn, which have a solution
minSn <- log10(k)/maxDHC5
# minSn is the starting point of the investigation

sn_rep <- 1:1000
#
psSn15 <- vector("list", length(sn_rep))
for (i in seq_along(sn_rep)) {
  sn <- minSn + 0.01 * (sn_rep[i] - 1) + 0.00001
  kc <- k^(1 / sn)  
  x1 <- xLower(SS, K005, kc)
  x2 <- xHigher(SS, K005, kc)
  p <- pt(x2, df = SS - 1) - pt(x1, df = SS - 1)
  psSn15[[i]] <- c(s = sn, p = p)
}
# I do not know what they are doing
psSn15 <- do.call(rbind, psSn15)
psSn15 <- as.data.frame(psSn15)
row.names(psSn15) <- NULL

plot(psSn15[, 1], psSn15[, 2],
     type = "l",  
     xlim = c(0.3, 2.0), 
     ylim = c(0, 0.6),  
     xlab = "sn", ylab = "p",
     frame.plot = TRUE  
)

