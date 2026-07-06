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

# Set random number seed, which is an integer between 1 and 1,000,000 generated
# at random.org
set.seed(339652)

# Helper function to simulate shelter choice
simulate_shelter_choice <- function(T_vect, S, mean_H, stdv_H, mean_D, stdv_D,
                                    min_W, max_W, min_rho, max_rho, Y,
                                    fixed_rho = NULL) {
  # Simulates shelter choice probabilities and calculates indifference point
  #
  # Args:
  #   T_vect: Vector of time periods between moves (in days)
  #   S: Number of samples per T
  #   mean_H, stdv_H: Mean and standard deviation for tipi construction time
  #   mean_D, stdv_D: Mean and standard deviation for tipi setup/teardown time
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
  alpha_D <- (mean_D / stdv_D)^2
  beta_D <- mean_D / (stdv_D^2)
  
  mid_W <- (min_W + max_W) / 2
  mid_rho <- (min_rho + max_rho) / 2
  
  p_vect <- c()
  mean_c_tipi <- c() 
  mean_c_wick <- c() 
  for (T in T_vect) {
    # Make S draws for each of H, M, W, and rho (if not fixed)
    H <- rgamma(n = S, shape = alpha_H, rate = beta_H)
    D <- rgamma(n = S, shape = alpha_D, rate = beta_D)
    W <- rtriangle(n = S, a = min_W, b = max_W, c = mid_W)
    
    if (is.null(fixed_rho)) {
      rho <- rtriangle(n = S, a = min_rho, b = max_rho, c = mid_rho)
    } else {
      rho <- rep(fixed_rho, S)
    }

    # Calculate r, the per period discount factor
    r <- (1 + rho)^(T/Y) - 1

    # Calculate time discounted labor costs. Both shelters must provide shelter
    # from day 0. The tipi is built once (H at day 0) and re-erected at each
    # subsequent move (D at days T, 2T, ...), giving H + D/r. The wickiup is
    # rebuilt at the start of EVERY camp, including day 0 (W at days 0, T, 2T,
    # ...), giving a perpetuity that starts at day 0: W + W/r.
    c_tipi <- H + D/r
    c_wick <- W + W/r

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
# distribution with a mean of 226.5 hours and a standard deviation of 0.25
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
mean_H <- 226.5
stdv_H <- mean_H * .25
alpha_H <- (mean_H / stdv_H)^2
beta_H <- mean_H / (stdv_H^2)

# D is the total set-up and tear-down time of a tipi (the per-move cost), drawn from a
# gamma distribution with a mean of .834 hours and a standard deviation of
# 0.25 times the mean.
mean_D <- .834
stdv_D <- mean_D * .25
alpha_D <- (mean_D / stdv_D)^2
beta_D <- mean_D / (stdv_D^2)

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
# Number of samples per T. Using 1,000,000 samples is enough to yield stable
# numbers and curves.
S <- 1000000

# Run simulations for different rho values
results_variable <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_D,
                                            stdv_D, min_W, max_W, min_rho,
                                            max_rho, Y)
results_min <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_D,
                                       stdv_D, min_W, max_W, min_rho, max_rho,
                                       Y, fixed_rho = min_rho)
results_mid <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_D,
                                       stdv_D, min_W, max_W, min_rho, max_rho,
                                       Y, fixed_rho = (min_rho + max_rho) / 2)
results_max <- simulate_shelter_choice(T_vect, S, mean_H, stdv_H, mean_D,
                                       stdv_D, min_W, max_W, min_rho, max_rho,
                                       Y, fixed_rho = max_rho)

# Compute deterministic transition point. Setting c_tipi = c_wick with mean
# values (H + D/r = W + W/r) gives (1+rho)^(T/Y) = (H - D)/(H - W), hence:
T_determ <- Y * log((mean_H - mean_D)/(mean_H - mid_W)) / log(1 + mid_rho)
moves_per_year_determ <- Y / T_determ

# Make Figure 1
# We limit the plot x-limits to 10 moves per year or less, since by then the
# probability has already plateaued to 1.
moves_per_year <- Y / T_vect
ind_fig1 <- which(moves_per_year <= 10)
x_fig1 <- moves_per_year[ind_fig1]
y_fig1 <- results_variable$p_vect[ind_fig1]
pdf('p_vs_yearly_moves.pdf')
    plot(x_fig1, y_fig1,
    xlab='Moves per Year',
    ylab='Probability Tipi is Preferred',
    type='l',
    lwd=3,
    xaxt='n')  # Suppress the default x-axis
    
    # Create a custom x-axis with x_indif included
    axis(1, at=c(0, 2, 4, 6, 8, 10, results_variable$x_indif), 
         labels=c("0", "2", "4", "6", "8", "10", 
                  sprintf("%.2f", results_variable$x_indif)))

    abline(v = results_variable$x_indif, col = "grey")
    abline(h = 0.5, col = "grey")
    
    # Add vertical line at moves_per_year_determ
    abline(v = moves_per_year_determ, col = "grey", lty = 2)
    
    # Add legend
    legend("bottomright",
           legend = c("Simulated Indifference", "Deterministic Transition"),
           col = c("grey", "grey"),
           lty = c(1,2),
           lwd = 2,
           bty = "n")
dev.off()

# Make Figure 2
pdf('p_vs_T.pdf')
    plot(T_vect, results_variable$p_vect,
    xlab='Days between Moves',
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

    abline(v = results_variable$T_indif, col = "grey")
    abline(h = 0.5, col = "grey")
    
    # Add vertical line at T_determ
    abline(v = T_determ, col = "grey", lty = 2)
    
    # Add legend
    legend("topright",
           legend = c("Simulated Indifference", "Deterministic Transition"),
           col = c("grey", "grey"),
           lty = c(1,2),
           lwd = 2,
           bty = "n")
dev.off()

# Make Figure 3 (p_vect versus moves per year with fixed rho values)
pdf('p_vs_yearly_moves_fixed_rho.pdf')
plot(x_fig1, results_min$p_vect[ind_fig1],
     xlab='Moves per Year', ylab='Probability Tipi is Preferred',
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
     xlab='Moves per Year',
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
cat("Deterministic:", moves_per_year_determ, "\n")
cat("Min rho:", results_min$x_indif, "\n")
cat("Mid rho:", results_mid$x_indif, "\n")
cat("Max rho:", results_max$x_indif, "\n")

# Print indifference points for T (days between moves)
cat("\nIndifference points (days between moves):\n")
cat("Variable rho:", results_variable$T_indif, "\n")
cat("Deterministic:", T_determ, "\n")
cat("Min rho:", results_min$T_indif, "\n")
cat("Mid rho:", results_mid$T_indif, "\n")
cat("Max rho:", results_max$T_indif, "\n")

# Print mean costs for 4 moves per year
ind4 <- which.min(abs(moves_per_year - 4))
c_tipi4 <- results_variable$mean_c_tipi[ind4]
c_wick4 <- results_variable$mean_c_wick[ind4]
cat("c_tipi at 4 moves:", c_tipi4, "\n")
cat("c_wick at 4 moves:", c_wick4, "\n")
cat("Relative difference at 4 moves:", (c_wick4-c_tipi4)/c_tipi4, "\n")

# Add a plot of tipi diameter versus total number of hides for the supplement
# See:
# https://albertashistoricplaces.com/2024/03/13/tipis-bison-and-dogs-visualizing-an-archaeological-feature-in-southern-alberta/
tipi_diameters <- c(2.31, 4.60, 8.57)
tipi_outer_hides <- c(6, 16, 42)
tipi_inner_hides <- c(3, 8, 17)
tipi_total_hides <- tipi_outer_hides + tipi_inner_hides
# Do a linear fit
tipi_fit <- lm(tipi_total_hides ~ tipi_diameters)
tipi_fit_intercept <- as.numeric(tipi_fit$coefficients)[1]
tipi_fit_slope <- as.numeric(tipi_fit$coefficients)[2]

# Predict the diameter for a tipi with 15 total hides, of which 5 are for the
# inner lining
target_total_hides <- 15
predicted_diameter <- (target_total_hides - tipi_fit_intercept)/tipi_fit_slope
print(paste("Predicted tipi diameter for 15 total hides:", predicted_diameter))

# Plot the total number of hides against the tipi diameter, along with the fit
pdf('tipi_total_hides_versus_diameter.pdf')
plot(tipi_diameters, tipi_total_hides,
     xlab='Tipi Diameter (meters)',
     ylab='Total Number of Hides',
     type='p',
     pch=16,
     lwd=3,
     col='black')
abline(tipi_fit, col='grey', lwd=3)
# Add dashed lines for the predicted values
abline(h = target_total_hides, col='grey', lwd=2, lty=2)
abline(v = predicted_diameter, col='grey', lwd=2, lty=2)
dev.off()

# =====================================================================
# Sensitivity analysis: per-location vs. drawn-once variability in W
# =====================================================================
# The wickiup build cost W is the recurring per-move cost. Its value is
# uncertain, but that uncertainty can be of two kinds: (i) we may not know
# the typical build time (epistemic; a single value applies to a given
# decision-maker across all camps), or (ii) the build time may genuinely
# vary from camp to camp with local material availability (inherent
# variability; a fresh value at each camp). These have the same mean cost
# but different variance: the per-location component averages over the
# indefinite sequence of camps.
#
# Following Supplemental Data 2, we write the cost at camp k of sample s as
# W_{s,k} = mid_W + a_s + b_{s,k}, splitting the total variance var_W into a
# drawn-once part a_s (variance (1 - f) * var_W, shared across the sample's
# camps) and a per-location part b_{s,k} (variance f * var_W, independent
# across camps). The discounted sum of the per-location parts,
# sum_{k>=0} b_{s,k} z_s^k (the code's S), has mean 0 and variance
# f * var_W / (1 - z_s^2). We evaluate it by direct simulation: we draw the
# per-camp b_{s,k} and accumulate the discounted sum, truncating where the
# discount weight z_s^k makes the omitted tail negligible (< 1e-4 of its
# variance). The shared part enters through the exact geometric sum
# (mid_W + a_s) / (1 - z_s). f = 0 reproduces the main model (W drawn once);
# f = 1 makes W fully per-location. In the code the sample index s is implicit:
# each per-sample quantity (a, S, z, ...) is a vector over the Ns samples
# drawn at a given number of moves.
#
# The mean and variance of this cost are derived analytically in Supplemental
# Data 2 (Mathematics): the mean is independent of f (so the median indifference
# point does not move) and the variance falls with f (so the spread narrows).
# Those are the mechanism; the core result computed here is the full
# distribution of the cross-over (indifference) point, obtained by simulation.
# We sweep f and record its median and 95% central range. Because
# P(tipi preferred at m) = P(indifference point < m), these are read off the
# crossing points of the probability curve.

var_W  <- (max_W - min_W)^2 / 24   # variance of the symmetric triangular W
half_W <- (max_W - min_W) / 2      # half-width, for the centered shared part

# Probability tipi is preferred at a given number of moves per year, for a
# given per-location variance fraction f.
sens_P <- function(moves, f, Ns) {
  Td   <- Y / moves
  Hs   <- rgamma(Ns, shape = alpha_H, rate = beta_H)
  Ds   <- rgamma(Ns, shape = alpha_D, rate = beta_D)
  rhos <- rtriangle(Ns, a = min_rho, b = max_rho, c = mid_rho)
  z    <- (1 + rhos)^(-Td / Y)
  perp <- 1 / (1 - z)                          # sum_{k>=0} z^k
  a    <- sqrt(1 - f) * rtriangle(Ns, a = -half_W, b = half_W, c = 0)
  # Direct simulation of S = sum_{k>=0} b_k z^k (per-location fluctuations),
  # truncated where the omitted tail contributes < 1e-4 of Var(S). max(z) sets
  # the truncation so the tail is negligible for every sample.
  S <- numeric(Ns)
  if (f > 0) {
    Kmax <- ceiling(log(1e-4) / (2 * log(max(z))))
    zpow <- rep(1, Ns)                         # z^0
    for (k in 0:Kmax) {
      bk   <- sqrt(f) * rtriangle(Ns, a = -half_W, b = half_W, c = 0)
      S    <- S + bk * zpow
      zpow <- zpow * z
    }
  }
  c_wick_s <- (mid_W + a) * perp + S           # shared part (exact) + per-camp (S)
  c_tipi_s <- Hs + Ds * (perp - 1)             # perp - 1 = z/(1-z) = 1/r
  mean(c_tipi_s < c_wick_s)
}

f_vals     <- c(0, 0.25, 0.5, 0.75, 1)
moves_grid <- seq(0.7, 9, by = 0.1)   # spans the 2.5%-97.5% crossings for all f
Ns_sens    <- 50000

sens_summary <- data.frame(f = f_vals, median = NA, lo95 = NA, hi95 = NA, width = NA)
sens_curves  <- matrix(NA, nrow = length(moves_grid), ncol = length(f_vals))
for (j in seq_along(f_vals)) {
  Pj <- sapply(moves_grid, sens_P, f = f_vals[j], Ns = Ns_sens)
  sens_curves[, j] <- Pj
  qs <- approx(x = Pj, y = moves_grid, xout = c(0.025, 0.5, 0.975), ties = mean)$y
  sens_summary$lo95[j]   <- qs[1]
  sens_summary$median[j] <- qs[2]
  sens_summary$hi95[j]   <- qs[3]
  sens_summary$width[j]  <- qs[3] - qs[1]
}

cat("\nSensitivity to per-location vs drawn-once variability in W\n")
cat("(f = fraction of W's variance that is per-location; f=0 is the main model):\n")
print(round(sens_summary, 3))

# Figure: probability curves for f = 0 (drawn once) and f = 1 (per location)
pdf('sensitivity_W_pcurves.pdf')
plot(moves_grid, sens_curves[, 1], type = 'l', lwd = 3, col = 'black',
     xlab = 'Moves per Year', ylab = 'Probability Tipi is Preferred',
     xlim = c(0, 10), ylim = c(0, 1))
lines(moves_grid, sens_curves[, length(f_vals)], lwd = 3, lty = 2, col = 'grey40')
abline(h = 0.5, col = 'grey', lty = 3)
legend('bottomright', bty = 'n', lwd = 3, lty = c(1, 2), col = c('black', 'grey40'),
       legend = c('f = 0 (W drawn once, epistemic)',
                  'f = 1 (W drawn per location)'))
dev.off()

# Figure: median and 95% central range of the indifference point vs f
pdf('sensitivity_W_range.pdf')
plot(sens_summary$f, sens_summary$median, type = 'b', pch = 16, lwd = 2,
     ylim = range(c(sens_summary$lo95, sens_summary$hi95)),
     xlab = 'Fraction of W variance that is per-location (f)',
     ylab = 'Indifference point (moves per year)')
lines(sens_summary$f, sens_summary$lo95, type = 'b', pch = 1, lty = 2)
lines(sens_summary$f, sens_summary$hi95, type = 'b', pch = 1, lty = 2)
legend('right', bty = 'n', pch = c(16, 1), lty = c(1, 2),
       legend = c('median', '95% central range'))
dev.off()

# =====================================================================
# Comparative-statics figure: indifference point vs. cost profile
# =====================================================================
# The deterministic indifference point,
#   Mhat = log(1 + rho) / log((H - D)/(H - W)),
# plotted as a function of the wickiup build cost W (with H and D at their
# means), for three discount rates. Shows how the housing cost profile moves
# the mobility threshold, before the specific tipi/wickiup values are inserted.
# Uses the closed form only (no Monte Carlo).
Mhat_det <- function(H, D, W, rho) log(1 + rho) / log((H - D) / (H - W))
W_seq <- seq(5, 30, by = 0.1)
pdf('comparative_statics.pdf')
plot(W_seq, Mhat_det(mean_H, mean_D, W_seq, min_rho),
     type = 'l', lwd = 2, col = 'black', ylim = c(0, 8),
     xlab = 'Wickiup build cost W (hours)',
     ylab = 'Indifference point (moves per year)')
lines(W_seq, Mhat_det(mean_H, mean_D, W_seq, mid_rho), lwd = 2, col = 'grey40')
lines(W_seq, Mhat_det(mean_H, mean_D, W_seq, max_rho), lwd = 2, col = 'grey70')
legend('topright', bty = 'n', lwd = 2, col = c('black', 'grey40', 'grey70'),
       legend = c(paste('rho =', min_rho), paste('rho =', mid_rho),
                  paste('rho =', max_rho)))
dev.off()