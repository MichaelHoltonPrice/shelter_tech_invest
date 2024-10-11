# This R script creates all of the results for the following scientific
# article:
#
# Welker and Price -- Homeward Bound: Exploring Decisions Involving Shelter and
#                     Transport
#
# To generate the results, run the script, after installing packages, if necessary.
#
# install.packages(c('MASS', 'triangle'))
# source('make_publication_results.R')

# Load necessary libraries
library(MASS)
library(triangle)

# Clear the workspace
rm(list=ls())

# Helper function to simulate shelter choice
simulate_shelter_choice <- function(T_vect, S, mean_H, stdv_H, mean_M, stdv_M,
                                    min_W, max_W, min_rho, max_rho, Y,
                                    fixed_rho = NULL) {
  # Simulates shelter choice probabilities and calculates indifference point
  #
  # Args:
  #   T_vect: Vector of time periods between moves (in days)
  #   S: Number of samples per T
  #   mean_H, stdv_H: Mean and standard deviation for tipi construction time
  #   mean_M, stdv_M: Mean and standard deviation for tipi setup/teardown time
  #   min_W, max_W: Min and max for wickiup construction time
  #   min_rho, max_rho: Min and max for annualized time discount factor
  #   Y: Number of days per year
  #   fixed_rho: Optional fixed value for rho (default: NULL)
  #
  # Returns:
  #   List containing p_vect, x_indif, and T_indif
  
  # Calculate gamma distribution parameters for H and M
  alpha_H <- (mean_H / stdv_H)^2
  beta_H <- mean_H / (stdv_H^2)
  alpha_M <- (mean_M / stdv_M)^2
  beta_M <- mean_M / (stdv_M^2)
  
  mid_W <- (min_W + max_W) / 2
  mid_rho <- (min_rho + max_rho) / 2
  
  p_vect <- c()
  mean_c_tipi <- c() 
  mean_c_wick <- c() 
  for (T in T_vect) {
    # Make S draws for each of H, M, W, and rho (if not fixed)
    H <- rgamma(n = S, shape = alpha_H, rate = beta_H)
    M <- rgamma(n = S, shape = alpha_M, rate = beta_M)
    W <- rtriangle(n = S, a = min_W, b = max_W, c = mid_W)
    
    if (is.null(fixed_rho)) {
      rho <- rtriangle(n = S, a = min_rho, b = max_rho, c = mid_rho)
    } else {
      rho <- rep(fixed_rho, S)
    }

    # Calculate r, the per period discount factor
    r <- (1 + rho)^(T/Y) - 1

    # Calculate time discounted labor costs
    c_tipi <- H + M/r
    c_wick <- W/r

    # Calculate probabilities
    num_tipi <- sum(c_tipi < c_wick)
    num_ties <- sum(c_tipi == c_wick)
    p_vect <- c(p_vect, (num_tipi + num_ties/2)/S)
    mean_c_tipi <- c(mean_c_tipi, mean(c_tipi))
    mean_c_wick <- c(mean_c_wick, mean(c_wick))
  }

  # Calculate indifference point
  moves_per_year <- Y / T_vect
  ind_indif <- which.min(abs(p_vect - .5))
  x_for_indif <- moves_per_year[ind_indif + c(-1, 0, 1)]
  y_for_indif <- p_vect[ind_indif + c(-1, 0, 1)]
  fit_for_indif <- lm(y_for_indif ~ x_for_indif)
  x_indif <- (0.5 - coef(fit_for_indif)[1]) / coef(fit_for_indif)[2]
  T_indif <- Y / x_indif

  # Calculate the mean values of the discounted cost


  return(list(p_vect = p_vect, x_indif = x_indif, T_indif = T_indif,
              mean_c_tipi = mean_c_tipi, mean_c_wick = mean_c_wick))
}

# Set parameter values.
#
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
# probability that each strategy (tipi or wickiup) is preferred per our time
# discounting model. We sample 1000 times for each value of T.
# Days per year
Y <- 365
T_vect <- seq(1, Y, by=1)
# Number of samples per T. Using 100,000 samples is enough to yield very smooth
# curves.
S <- 100000

# Run simulations for different rho values
results_variable <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_M,
                                            stdv_M, min_W, max_W, min_rho,
                                            max_rho, Y)
results_min <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_M,
                                       stdv_M, min_W, max_W, min_rho, max_rho,
                                       Y, fixed_rho = min_rho)
results_mid <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_M,
                                       stdv_M, min_W, max_W, min_rho, max_rho,
                                       Y, fixed_rho = (min_rho + max_rho) / 2)
results_max <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_M,
                                       stdv_M, min_W, max_W, min_rho, max_rho,
                                       Y, fixed_rho = max_rho)

# Make Figure 1
# We limit the plot x-limits to 10 moves per year or less, since by then the
# probability has already plateaued to 1.
moves_per_year <- Y / T_vect
ind_fig1 <- which(moves_per_year <= 10)
x_fig1 <- moves_per_year[ind_fig1]
y_fig1 <- results_variable$p_vect[ind_fig1]
pdf('p_vs_yearly_moves.pdf')
    plot(x_fig1, y_fig1,
    xlab='Moves per year',
    ylab='Probability Tipi is Preferred',
    type='l',
    lwd=3,
    xaxt='n')  # Suppress the default x-axis
    
    # Create a custom x-axis with x_indif included
    axis(1, at=c(0, 2, 4, 6, 8, 10, results_variable$x_indif), 
         labels=c("0", "2", "4", "6", "8", "10", 
                  sprintf("%.2f", results_variable$x_indif)))

    abline(v = results_variable$x_indif, col = "grey", lty = 2)
    abline(h = 0.5, col = "grey", lty = 2)
dev.off()

# Make Figure 2
pdf('p_vs_T.pdf')
    plot(T_vect, results_variable$p_vect,
    xlab='Days between moves',
    ylab='Probability Tipi is Preferred',
    type='l',
    lwd=3,
    xaxt='n')  # Suppress the default x-axis
    
    # Create a custom x-axis with T_indif included
    T_indif_text <- sprintf("%.1f", results_variable$T_indif)
    axis(1, at=c(0, 100, 200, 300), 
         labels=c("0", "100", "200", "300"))
    axis(1, at=c(results_variable$T_indif), 
         labels=c(T_indif_text),
         line=1)

    abline(v = results_variable$T_indif, col = "grey", lty = 2)
    abline(h = 0.5, col = "grey", lty = 2)
dev.off()

# Make Figure 3 (p_vect versus moves per year with fixed rho values)
pdf('p_vs_yearly_moves_fixed_rho.pdf')
plot(x_fig1, results_min$p_vect[ind_fig1],
     xlab='Moves per year', ylab='Probability Tipi is Preferred',
     type='l', lwd=2, col='black', xaxt='n', ylim=c(0,1))
lines(x_fig1, results_mid$p_vect[ind_fig1], col='gray50', lwd=2)
lines(x_fig1, results_max$p_vect[ind_fig1], col='gray70', lwd=2)

axis(1, at=c(0, 2, 4, 6, 8, 10))
abline(h = 0.5, col = "gray50", lty = 2)

legend("bottomright", 
       legend=c(paste("Min rho =", min_rho), 
                paste("Mid rho =", (min_rho + max_rho) / 2),
                paste("Max rho =", max_rho)),
       col=c("black", "gray50", "gray70"), 
       lty=c(1, 1, 1), lwd=2, cex=0.8)
dev.off()

# Make Figure 4 (mean values of c versus moves per year)
moves_per_year <- Y / T_vect
ind_fig4 <- which(moves_per_year <= 10)
x_fig4 <- moves_per_year[ind_fig4]

pdf('mean_c_vs_moves_per_year.pdf')
plot(x_fig4, results_variable$mean_c_tipi[ind_fig4],
     xlab='Moves per year',
     ylab='Mean Discounted Labor Cost',
     type='l', lwd=3, col='black',
     ylim=range(c(results_variable$mean_c_tipi[ind_fig4], 
                  results_variable$mean_c_wick[ind_fig4])))
lines(x_fig4, results_variable$mean_c_wick[ind_fig4], col='black', lwd=3,
      lty=3)

legend("topleft", 
       legend=c("Tipi", "Wickiup"),
       col="black", 
       lty=c(1, 3),
       lwd=3, cex=0.8)

dev.off()

# Print indifference points
cat("Indifference points (moves per year):\n")
cat("Variable rho:", results_variable$x_indif, "\n")
cat("Min rho:", results_min$x_indif, "\n")
cat("Mid rho:", results_mid$x_indif, "\n")
cat("Max rho:", results_max$x_indif, "\n")

# Print indifference points for T (days between moves)
cat("\nIndifference points (days between moves):\n")
cat("Variable rho:", results_variable$T_indif, "\n")
cat("Min rho:", results_min$T_indif, "\n")
cat("Mid rho:", results_mid$T_indif, "\n")
cat("Max rho:", results_max$T_indif, "\n")