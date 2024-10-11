library(MASS)
library(triangle)

# Clear the workspace
rm(list=ls())

# Observed wickiup construction times:
wickiup_times <- c(0.5, 0.5, 1, 1, 2)

# Set the offset
wickiup_offset <- 0.25

# Adjust the data by subtracting the offset
adjusted_times <- wickiup_times - wickiup_offset

# Perform MLE to fit the Gamma distribution to the adjusted data
wickiup_gamma_fit <- fitdistr(adjusted_times, densfun = "gamma")

# Display the estimated parameters
print(wickiup_gamma_fit)

# Alpha and beta for gamma distribution
wickiup_alpha <- wickiup_gamma_fit$estimate["shape"]
# Scale parameter (not rate)
wickiup_beta  <- 1 / wickiup_gamma_fit$estimate["rate"]

# samples <- rgamma(10, shape = wickiup_alpha, scale = wickiup_beta)
# Make a histogram and density plot for the Wickiup construction time
# Make a histogram and density plot for the Wickiup construction time
#pdf('wickiup_density.pdf')
    W_vect <- seq(0, 4, by=0.01)
    hist(wickiup_times, freq=FALSE, breaks=seq(0.125, 4.125, by=.25),
         xlab='Labor Hours', ylab='Density', main="Wickiup construction time")
    lines(W_vect, dgamma(W_vect - wickiup_offset, shape=wickiup_alpha, scale=wickiup_beta), lwd=3)
#dev.off()

# H is the construction time for a tipi, which is drawn from a gamma
# distribution with a mean of 211.5 hours and a standard deviation of 0.25
# times the mean. The gamma distribution probability density function is
#
# f(x, alpha, beta) = beta^alpha * x^(alpha-1) * exp(-beta*x) / Gamma(alpha),
#
# where alpha is the shape parameter, beta is the rate parameter, and Gamma is
# the gamma function. The mean is mu_gamma = alpha/beta and the variance is
# var_gamma = alpha / beta^2. mu_gamma / var_gamma = beta. alpha and beta, in
# terms of the mean and standard deviation, are:
#
# alpha = mu_gamma^2 / var_gamma  = mu_gamma^2 / stdv_gamma^2
# beta = mu_gamma / var_gamma = mu_gamma / stdv_gamma^2
mean_H <- 211.5
stdv_H <- mean_H * .25
alpha_H <- (mean_H / stdv_H)^2
beta_H <- mean_H / (stdv_H^2)

# M is the total set-up and tear-down time of a tipi, which is drawn from a
# gamma distribution with a mean of .834 hours and a standard deviation of
# 0.25 times the mean.
mean_M <- .834
stdv_M <- mean_M * .25
alpha_M <- (mean_M / stdv_M)^2
beta_M <- mean_M / (stdv_M^2)

# W is the construction time of a wickiup, which is drawn from a symmetric
# triangular distribution (isosceles triangular distribution) ranging from
# 5 hours to 15 hours, which we represent by min_W and max_W. The mean (also
# the midpoint and mode) is 10 hours, which we represent by mid_W.
min_W <- 5
max_W <- 15
mid_W <- (min_W + max_W) / 2
# sample <- rtriangle(n = 1, a = a, b = b, c = c)

# rho is the annualized time discount factor, which is drawn from a symmetric
# triangular distribution that ranges from 8% to 20%.
min_rho <- 0.08
max_rho <- 0.20
mid_rho <- (min_rho + max_rho) / 2

# T is the time in days between moves. It is the x variable for most of our
# plots. We iterate over values of T then iterate over samples to identify the
# probability that each strategy (tipi or wickiup) is prefered per our time
# disccounting model. We sample 1000 times for each value of T.
T_vect <- seq(1, 365, by=1)
# Number of samples per T
S <- 1000
# Days per year
#Y <- 365.2425
Y <- 365
p_vect <- c()
for (T in T_vect) {
    num_tipi <- 0
    num_wick <- 0
    num_ties <- 0
    for (s in 1:S) {
        H <- rgamma(n = 1, shape = alpha_H, rate = beta_H)
        M <- rgamma(n = 1, shape = alpha_M, rate = beta_M)
        W <- rtriangle(n = 1, a = min_W, b = max_W, c = mid_W)
        rho <- rtriangle(n = 1, a = min_rho, b = max_rho, c = mid_rho)
        # Calculate r, the per period disccount factor
        r <- (1 + rho)^(T/Y) - 1

        # The time discounted labor cost of the tipi strategy
        c_tipi <- H + M/r
        # The time discounted labor cost of the wickiup strategy
        c_wick <- W/r

        if (c_tipi < c_wick) {
            num_tipi <- num_tipi + 1
        } else if(c_tipi == c_wick) {
            num_ties <- num_ties + 1
        } else {
            num_wick <- num_wick + 1
        }
    }
    p_vect <- c(p_vect, (num_tipi + num_ties/2)/S )
}

pdf('p_vs_T.pdf')
    plot(T_vect, p_vect, xlab='Days between moves', ylab='Probability Tipi is Preferred')
dev.off()

pdf('p_vs_yearly_moves.pdf')
    plot(Y/T_vect, p_vect, xlab='Moves per year', ylab='Probability Tipi is Preferred', xlim=c(0,10))
dev.off()