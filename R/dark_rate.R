# Functions to Estimate Unobserved Illegal Activity

#' Initiate Parameters
#'
#' @param alpha set initial parameters for optimization: c(c, phi, eta, d, xi, alpha)
#'
#' @return vector of parameters
#' @examples
#' params0 <- set_inits();
#' @export
set_inits <- function(c = .5, phi = .5, eta = .5, d = 0.5, xi = .5, alpha = .5) {
  c(c = c, phi = phi, eta = eta, d = d, xi = xi, alpha = alpha)
}

# Objective Function
LF <- function(params, data) {
  pC <- pnorm(params['c'] + params['phi'] * data[, 'Z'] + params['eta'] * data[, 'R_prev'])
  pE <- pnorm(params['d'] + params['xi'] * data[, 'W'] + params['alpha'] * data[, 'R_prev'])
  log_likelihood <- mean(
    log((pC * pE)^(data[, 'R'] == 1) *
          ((1 - pC * pE)^(data[, 'R'] == 0)))
  )
  return(-log_likelihood)
}

LF_gradient <- function(params, data) {
  z_C <- params['c'] + params['phi'] * data[, 'Z'] + params['eta'] * data[, 'R_prev']
  z_E <- params['d'] + params['xi'] * data[, 'W'] + params['alpha'] * data[, 'R_prev']

  # CDFs (probabilities)
  log_pC <- pnorm(z_C, log.p = TRUE)
  log_pE <- pnorm(z_E, log.p = TRUE)
  pC <- exp(log_pC)
  pE <- exp(log_pE)

  # PDFs (densities) - these are the derivatives
  phi_C <- dnorm(z_C)  # ϕ(z_C) = d/dz Φ(z_C)
  phi_E <- dnorm(z_E)  # ϕ(z_E) = d/dz Φ(z_E)

  R <- data[, 'R']

  grad <- numeric(6)
  names(grad) <- c('c', 'phi', 'eta', 'd', 'xi', 'alpha')

  for(i in 1:length(R)) {
    if(R[i] == 1) {
      grad['c'] <- grad['c'] - phi_C[i] / pC[i]
      grad['phi'] <- grad['phi'] - phi_C[i] / pC[i] * data[i, 'Z']
      grad['eta'] <- grad['eta'] - phi_C[i] / pC[i] * data[i, 'R_prev']
      grad['d'] <- grad['d'] - phi_E[i] / pE[i]
      grad['xi'] <- grad['xi'] - phi_E[i] / pE[i] * data[i, 'W']
      grad['alpha'] <- grad['alpha'] - phi_E[i] / pE[i] * data[i, 'R_prev']
    } else {
      p_CE <- pC[i] * pE[i]
      factor <- p_CE / (1 - p_CE)
      grad['c'] <- grad['c'] + factor * phi_C[i] / pC[i]
      grad['phi'] <- grad['phi'] + factor * phi_C[i] / pC[i] * data[i, 'Z']
      grad['eta'] <- grad['eta'] + factor * phi_C[i] / pC[i] * data[i, 'R_prev']
      grad['d'] <- grad['d'] + factor * phi_E[i] / pE[i]
      grad['xi'] <- grad['xi'] + factor * phi_E[i] / pE[i] * data[i, 'W']
      grad['alpha'] <- grad['alpha'] + factor * phi_E[i] / pE[i] * data[i, 'R_prev']
    }
  }

  n <- nrow(data)
  grad <- grad / n

#  cat("Gradient: ", grad, "\n")
  return(grad)
}



#' Wrapper Function for optim.R
#'
#' Estimates parameter for objective function using optim.R
#' @param data Observed data
#' @param params Initial parameter vector for optimization
#' @param fct Objective function
#' @param grad Gradient function of objective function (optinal, enhances robustness and runtime of optimization algorithm)
#'
#' @return List of 6: par (estimates), value (optimized value of objective function), counts (number of calls to fn and gr), convergence (0 if converged), message, hessian
#' @examples
#' parm_sim0 <- set_parms()
#' parm_sim <- parm_sim0[3:length(parm_sim0)]
#' data0 <- sim_data(n, parm_sim0)
#' params0 <- set_inits()
#' result <- wrap_optim(data[,1:4], params0, LF)
#' @export
wrap_optim <- function(data, params, fct, grad=NULL){
  result <- optim(
    par = params,             # Initial parameters
    fn = fct,                    # Objective function
    gr = grad,                 # Faster and more robust if the gradient is provided
    data = data,               # Data passed to the function
    method = "BFGS",           # Gradient-based method (Alternatives: "Nelder-Mead", "CG", "L-BFGS-B")
    hessian = TRUE,
    control = list(
      maxit = 20000,           # Maximum iterations
      trace =1,                # Print progress
      ndeps = c(rep(1e-5, length(params))),  # Step size
      reltol = 1e-10           # Tolerance for convergence
    )
  )
  return(result)
}


get_results_parm <- function(fit, decs, n){
  hessian <- fit$hessian * n
  vcov_matrix <- solve(hessian)
  se <- sqrt(diag(vcov_matrix))
  ci <- data.frame(
    estimate = fit$par,
    lower = fit$par - 1.96*se,
    upper = fit$par + 1.96*se
  )
  results <- data.frame(
    estimates = round(fit$par, decs),
    SE = round(sqrt(diag(solve(hessian))), decs),
    ci.lower = round(ci$lower, decs),
    ci.upper = round(ci$upper, decs),
    #    p_value = 2 * (1 - pnorm(abs(results$Estimate/results$SE))),
    row.names = names(fit$par)
  )
}

get_C <- function(Z, R_prev, parm) {
  pC <- pnorm(parm['c'] + parm['phi'] * Z + parm['eta'] * R_prev)
}

get_E <- function(W, R_prev, parm) {
  pE <- pnorm(parm['d'] + parm['xi'] * W + parm['alpha'] * R_prev)
}

get_C_E <- function(data, parm){
  pC <- get_C(data[, 'Z'], data[, 'R_prev'], parm)
  pE <- get_E(data[, 'W'], data[, 'R_prev'], parm)
  result_C <- data.frame(estimates=mean(pC), SE=sd(pC), ci.lower=quantile(pC, 0.025), ci.upper=quantile(pC, 0.975), row.names='C_hat')
  result_E <- data.frame(estimates=mean(pE), SE=sd(pE), ci.lower=quantile(pE, 0.025), ci.upper=quantile(pE, 0.975), row.names='E_hat')
  return(rbind(result_E, result_C))
}

#' Combine results
#'
#' Combine results after parameter estimation
#' @param data Observed data
#' @param fit Return object of optim.R (list of 6)
#' @param decs Number of digits to round results
#'
#' @return Dataframe result table with estimates, standard errors and 95% confidence intervals
#' @examples
#' parm_sim0 <- set_parms()
#' parm_sim <- parm_sim0[3:length(parm_sim0)]
#' data0 <- sim_data(n, parm_sim0)
#' data <- data0[,1:4]
#' params0 <- set_inits()
#' result <- wrap_optim(data, params0, LF)
#' result_table <- get_results(data, result)
get_results <- function(data, fit, decs=4){
  results_parms <- get_results_parm(fit, decs, nrow(data))
  results_C_H <- get_C_E(data, fit$par)
  results_C_H <- round(results_C_H, decs)
  return(rbind(results_parms, results_C_H))
}

#' Calculate Dark Rate
#'
#' Calculate Dark Rate based on observed und estimated unobserved activity
#' @param observed binary vector (1 if observed event, 0 else)
#' @param estimate binary vector (1 if estimated event, 0 else)
#'
#' @return Decimal number DRC = 1 - (observed/estimate)
#' @examples
#' parm_sim0 <- set_parms()
#' parm_sim <- parm_sim0[3:length(parm_sim0)]
#' data0 <- sim_data(n, parm_sim0)
#' DRC <- get_dark_rate(sum(data0[,'R']), sum(data0[,'C']))
get_dark_rate <- function(observed, estimate) {
  return((estimate - observed)/estimate)
}
#get_dark_rate(6, 9)
