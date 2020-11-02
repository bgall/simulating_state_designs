#################################################
# Define function: sim_xsection
#
# Description: simulates data for
# cross-sectional research designs with
# binary treatments and continuous or
# binary outcomes analyzed via linear 
# regression models. Useful for power
# calculations and analyses of the
# statistical properties of various
# research designs employed to study
# the effect of binary treatments.
#
# Arguments:
# n_states: number of total states
# n_treated: number of "treated" states
# y_control: mean outcome for control states
#            or mean probability outcome = 1
# y_treated: mean outcome for treated states
#            or mean probability outcome = 1
# y_sd_control: sd of outcome for control states
#               used for continuous outcomes
# y_sd_treated: sd of outcome for treated states
#               used for continuous outcomes
# sims: number of simulations
# cores: number of CPU cores to use for parallel processing
# response: either "binary" or "continuous"
#           OLS is always used
# keep_terms: terms from the model summary to keep
#             from each simulation
#################################################

# Load required packages
library(foreach)
library(doParallel)
library(broom)


# Define function
sim_xsection <- function(n_states = 50,
                         n_treated,
                         y_control,
                         y_treat,
                         y_sd_control,
                         y_sd_treated,
                         sims,
                         cores = 1,
                         response = "binary",
                         keep_terms) {
  ##############################################
  # TEST INPUTS, RETURN ERRORS AS NEEDED
  #############################################
  
  # Test the number of states treated is fewer or equal to total states
  if (any(n_treated > n_states,
          n_treated < 0))
    stop("Number of treated states must be > 0 and < # of total states.")
  
  # For binary DV, test inputs are in the range of the binomial
  # distribution parameters
  if (response == "binary") {
    if (any(y_control < 0,
            y_treat < 0,
            y_control > 1,
            y_treat > 1))
      stop("Response is binary. y_control and y_treat must be in [0,1]")
  }
  
  ##############################################
  # Generate quantities needed regardless
  # of model
  #############################################
  
  # Calculate number of control states
  n_control <- n_states - n_treated
  
  ##############################################
  # Set up parallel processing for simulation
  #############################################
  
  # Register a cluster
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  
  ##############################################
  # Draw outcomes for treated and controls:
  # assuming treatment is assigned as-if-randomly
  # across states.
  #############################################
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Binary DV
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (response == "binary") {
    sims_out <- foreach::foreach(i = 1:sims,
                                 .packages = c("dplyr", "broom")) %dopar% {
                                   # simulate data
                                   sim_data <-
                                     dplyr::bind_rows(data.frame(
                                       y = rbinom(
                                         n = n_treated,
                                         size = 1,
                                         prob = y_treat
                                       ),
                                       treated = 1
                                     ),
                                     data.frame(
                                       y = rbinom(
                                         n = n_control,
                                         size = 1,
                                         prob = y_control
                                       ),
                                       treated = 0
                                     ))
                                   
                                   # Fit linear probability model for each simulation
                                   # and keep only results for desired model parameters
                                   sim_data %>%
                                     lm(formula = y ~ treated, data = .) %>%
                                     broom::tidy() %>%
                                     dplyr::filter(term %in% keep_terms)
                                 }
    
    # Terminate cluster
    parallel::stopCluster(cl)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Continuous DV
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (response == "continuous") {
    sims_out <- foreach::foreach(i = 1:sims,
                                 .packages = c("dplyr", "broom")) %dopar% {
                                   # simulate data
                                   sim_data <-
                                     dplyr::bind_rows(data.frame(
                                       y = rnorm(n = n_treated, 
                                                 mean = y_treat, 
                                                 sd = y_sd_treated),
                                       treated = 1
                                     ),
                                     data.frame(
                                       y = rnorm(n = n_control, 
                                                 mean = y_control, 
                                                 sd = y_sd_control),
                                       treated = 0
                                     ))
                                   
                                   # Fit linear probability model for each simulation
                                   # and keep only results for desired model parameters
                                   sim_data %>%
                                     lm(formula = y ~ treated, data = .) %>%
                                     broom::tidy() %>%
                                     dplyr::filter(term %in% keep_terms)
                                 }
    
    # Terminate cluster
    parallel::stopCluster(cl)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Tidying regardless of DV type
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Cast simulation object into a single tidy
  # data frame of parameter estimates for all
  # simulated models
  sims_df <- sims_out %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(true_treat_minus_control = y_treat - y_control)
  
  # Return simulated terms
  sims_df
}

######################################
# Working examples
######################################

# Example unquoted parameter values
# for testing chunks within the function
#n_states = 50
#n_treated = 40
#y_control = 0.5
#y_treat = 0.2
#y_sd_control = 1
#y_sd_treated = 1
#sims = 10
#cores = 1
#response = "binary"
#keep_terms = "treated"

# Example simulation:
# binary DV, 2 cores, 10 simulation
# keep only the treatment variable results
sim_xsection(
  n_states = 50,
  n_treated = 40,
  y_control = 0.5,
  y_treat = 0.9,
  y_sd_control = 1,
  y_sd_treated = 1,
  sims = 10,
  cores = 2,
  response = "continuous",
  keep_terms = "treated"
)

# Example simulation:
# binary DV, 2 cores, 10 simulation
# two different assumed treatment effects
# returns a list.
lapply(
  X = c(0.7, 0.9),
  FUN = function(x) {
    sim_xsection(
      n_treat = 20,
      y_control = x,
      y_treat = 0.5,
      sim = 100,
      cores = 2,
      response = "binary",
      keep_terms = "treated"
    )
  }
)

# Now let's look at power and errors!

# Simulate treatment effects of 1% and 4%
# with 20 states treated and 30 control states
# with a linear probability model and 10,000
# simulations each.
smallfx <- lapply(
  X = c(0.51, 0.55),
  FUN = function(x) {
    sim_xsection(
      n_treat = 20,
      y_control = x,
      y_treat = 0.5,
      sim = 10000,
      cores = 2,
      response = "binary",
      keep_terms = "treated"
    )
  }
)

# Get departure from true value
# and create indicator for significant results
smallfx_df <- smallfx %>%
              dplyr::bind_rows() %>%
  dplyr::mutate(error_m = true_treat_minus_control - estimate,
                significant = ifelse(p.value < 0.05, 1, 0))

# Calculate power: share of simulations finding a 
# significant effect at p < 0.05
smallfx_df %>% 
  dplyr::group_by(true_treat_minus_control) %>%
  dplyr::summarise(power = mean(significant))

# Calculate average effect *conditional on a
# a statistically significant effect*
smallfx_df %>%
  dplyr::filter(significant == 1) %>%
  dplyr::group_by(true_treat_minus_control) %>%
  dplyr::summarise(mean_effect = mean(estimate))
