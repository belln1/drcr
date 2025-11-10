rm(list = ls())
library(drcr)
# library(assertthat)
# source(file = "R/simulate_crime.R")
# source(file = "R/dark_rate.R")

set.seed(123) # Set Seed for Reproducibility

n <- 100000
decs <- 4

# Create Data
parm_sim0 <- set_parms()
parm_sim <- parm_sim0[3:length(parm_sim0)]
data0 <- sim_data(n, parm_sim0)
data <- data0[,1:4]

# Estimate Parameters
params0 <- set_inits()
result <- wrap_optim(data, params0, LF)
#result <- wrap_optim(data, params0, LF, grad=LF_gradient)

# Show Results
result_table <- get_results(data, result)
sim.values <- c(parm_sim, E_hat=mean(data0[,'E']), C_hat=mean(data0[,'C']))
result_table <- cbind(sim.values, result_table)
result_table

# Estimate Dark Rate: DRC = 1 - R/C
DRC <- get_dark_rate(mean(data[,'R']), result_table['C_hat',]$estimates)
DRC

################################################################################################



