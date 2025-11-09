# Functions to Simulate Observed and Unobserved Illegal Activity

#' Set Parameters
#'
#' @param alpha Parameters for simulation: c(rho, tau, c, phi, eta, d, xi, alpha)
#'
#' @return vector of parameters
#' @examples
#' parm_sim0 <- set_parms();
#' @export
set_parms <- function(rho = 0, tau = 0, c = 0.4, phi = 0.3, eta = -0.2, d = 0.5, xi = 0.4, alpha = 0.2){
  c(rho = rho, tau = tau, c = c, phi = phi, eta = eta, d = d, xi = xi, alpha = alpha)
}

sim_R_prev <- function(v, parms, Z, W){
  as.numeric(v < parms['rho'] * Z + parms['tau'] * W)
}

sim_C <- function(u, parms, Z, R_prev){
  as.numeric(u < parms['c'] + parms['phi'] * Z + parms['eta'] * R_prev)
}
sim_E <- function(e, parms, W, R_prev){
  as.numeric(e < parms['d'] + parms['xi'] * W + parms['alpha'] * R_prev)
}

#' Simulate Crime
#'
#' Simulation of entities and time periods with observed and unobserved illegal cases
#' @param n Number of observations
#' @param parms Parameters for simulation: c(rho, tau, c_C, phi, eta, c_E, xi, alpha)
#'
#' @return simulated dataset as matrix
#' @examples
#' data <- sim_data();
#' @export
sim_data <- function(n, parms) {
  Z <- rbinom(n, 2, 0.5)
  W <- rbinom(n, 2, 0.5)
  v <- rnorm(n)
  u <- rnorm(n)
  e <- rnorm(n)
  R_prev <- sim_R_prev(v, parms, Z, W)
  C <- sim_C(u, parms, Z, R_prev)
  E <- sim_E(e, parms, W, R_prev)
  R <- E * C
  data <- cbind(Z, W, R_prev, R, C, E)
  return(data)
}
