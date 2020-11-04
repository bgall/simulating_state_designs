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
# seed: random number generator seed for replication
#################################################

# Load required packages
library(foreach)
library(doParallel)
library(broom)
library(estimatr)

# Define function
sim_panel <- function(n_states = 50,
                      n_treated,
                      y_control,
                      y_treat,
                      y_sd_control,
                      y_sd_treated,
                      rand_intercept_mean = 0,
                      rand_intercept_sd,
                      sims,
                      years,
                      cores = 1,
                      response = "binary",
                      keep_terms = NA,
                      seed = 999) {
  ##############################################
  # Test inputs
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
  
  # Set seed
  set.seed(seed)
  
  # Calculate number of control states
  n_control <- n_states - n_treated
  
  # Generate a "skeleton" of state-years
  skeleton <- data.frame(
    state = as.factor(rep(1:n_states, times = years)),
    wave = rep(1:years, each = n_states)
  )
  
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
                                 .packages = c("dplyr", "broom", "estimatr")) %dopar% {
                                  
                                    # Add treatment status to skeleton
                                   # The number of (untreated) units will be determined
                                   # on an annual basis: the same number
                                   # of units are treated every year but which
                                   # units are treated is random
                                   template <- skeleton %>%
                                     dplyr::mutate(assign_rng = runif(n = n())) %>%
                                     dplyr::group_by(wave) %>%
                                     dplyr::mutate(rank_rng = rank(assign_rng)) %>%
                                     dplyr::ungroup() %>%
                                     dplyr::mutate(treated = dplyr::if_else(rank_rng <= n_treated, 1, 0))
                                   
                                   # Draw state intercepts, assuming constant across waves
                                   state_fx <-
                                     data.frame(
                                       state = as.factor(1:n_states),
                                       p_state = rnorm(n = n_states,
                                                       mean = rand_intercept_mean,
                                                       sd = rand_intercept_sd)
                                     )
                                   
                                   # Add state intercepts to template
                                   template <-
                                     dplyr::left_join(x = template,
                                                      y = state_fx,
                                                      by = "state")
                                   
                                   # Add effects based on treatment status
                                   # and calculate overall effects from constituents
                                   template <-
                                     template %>% dplyr::mutate(
                                       p_condition = dplyr::if_else(treated == 1, y_treat, y_control),
                                       p_sum = p_condition + p_state,
                                       p_actual = dplyr::case_when(p_sum > 1 ~ 1,
                                                                   p_sum < 0 ~ 0,
                                                                   TRUE ~ p_sum)
                                     )
                                   
                                   # Simulate outcomes using drawn
                                   # probabilities
                                   sim_data <-
                                     template %>% dplyr::mutate(y = rbinom(
                                       n = n(),
                                       size = 1,
                                       prob = p_actual
                                     ))
                                   
                                   # Fit linear probability model for each simulation
                                   # and keep only results for desired model parameters
                                 
                                   # Model WITH state fixed effects
                                    if(state_fe == TRUE) {
                                    results <- sim_data %>%
                                     estimatr::lm_robust(formula = y ~ treated, 
                                                         data = ., 
                                                         clusters = state,
                                                         se_type = "HC1",
                                                         fixed_effects = ~ state) %>%
                                                         broom::tidy()
                                  }
                                   
                                   # Model WITHOUT state fixed effects
                                   if(state_fe == TRUE) {
                                     results <- sim_data %>%
                                       estimatr::lm_robust(formula = y ~ treated, 
                                                           data = ., 
                                                           clusters = state,
                                                           se_type = "HC1") %>%
                                       broom::tidy()
                                   }
                                   
                                   if (!is.na(keep_terms))
                                     results %>% dplyr::filter(term %in% keep_terms)
                                 }
    
    # Terminate cluster
    parallel::stopCluster(cl)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Continuous DV
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  if (response == "continuous") {
    sims_out <- foreach::foreach(i = 1:sims,
                                 .packages = c("dplyr", "broom", "estimatr")) %dopar% {
                                  
                                    # Add treatment status to skeleton
                                   # The number of (untreated) units will be determined
                                   # on an annual basis: the same number
                                   # of units are treated every year but which
                                   # units are treated is random
                                   template <- skeleton %>%
                                     dplyr::mutate(assign_rng = runif(n = n())) %>%
                                     dplyr::group_by(wave) %>%
                                     dplyr::mutate(rank_rng = rank(assign_rng)) %>%
                                     dplyr::ungroup() %>%
                                     dplyr::mutate(treated = dplyr::if_else(rank_rng <= n_treated, 1, 0))
                                   
                                   # Draw state intercepts, assuming constant across waves
                                   state_fx <-
                                     data.frame(
                                       state = as.factor(1:n_states),
                                       p_state = rnorm(n = n_states,
                                                       mean = rand_intercept_mean,
                                                       sd = rand_intercept_sd)
                                     )
                                   
                                   # Add state intercepts to template
                                   template <-
                                     dplyr::left_join(x = template,
                                                      y = state_fx,
                                                      by = "state")
                                   
                                   # Add effects based on treatment status,
                                   # calculate overall effects from constituents,
                                   # and choose sd of outcome
                                   template <-
                                     template %>% dplyr::mutate(
                                       p_condition = dplyr::if_else(treated == 1, y_treat, y_control),
                                       p_sum = p_condition + p_state,
                                       p_actual = dplyr::case_when(p_sum > 1 ~ 1,
                                                                   p_sum < 0 ~ 0,
                                                                   TRUE ~ p_sum),
                                       y_sd_actual = dplyr::if_else(treated == 1, y_sd_treated, y_sd_control)
                                     )
                                   
                                   # Simulate outcomes using drawn
                                   # probabilities
                                   sim_data <-
                                     template %>% dplyr::mutate(y = rnorm(
                                       n = n(),
                                       mean = p_actual,
                                       sd = y_sd_actual
                                     ))
                                   
                                   # Fit linear probability model for each simulation
                                   # and keep only results for desired model parameters
                                   
                                   # Model WITH state fixed effects
                                   if(state_fe = TRUE){
                                   results <- sim_data %>%
                                     estimatr::lm_robust(formula = y ~ treated, 
                                                         data = ., 
                                                         clusters = state,
                                                         se_type = "HC1",
                                                         fixed_effects = ~ state) %>%
                                                         broom::tidy()
                                   }
                                   
                                   # Model WITHOUT state fixed effects
                                   if(state_fe = FALSE){
                                     results <- sim_data %>%
                                       estimatr::lm_robust(formula = y ~ treated, 
                                                           data = ., 
                                                           clusters = state,
                                                           se_type = "HC1") %>%
                                       broom::tidy()
                                   }
                                   
                                   if (!is.na(keep_terms))
                                     results %>% dplyr::filter(term %in% keep_terms)
                                 }
    
    # Terminate cluster
    parallel::stopCluster(cl)
  }
  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Tidying regardless of DV type
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Cast simulation object into a single tidy
# data frame of parameter estimates for all
# simulated models
sims_df <- sims_out %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(true_treat_minus_control = y_treat - y_control)

# Return selected model terms across
# all simulations
sims_df
}

######################################
# Working examples
######################################

# Example unquoted parameter values
# for testing chunks within the function
# n_states = 50
# n_treated = 10
# y_control = 0.20
# y_treat = 0.45
# y_sd_control = 1
# y_sd_treated = 1
# rand_intercept_mean = 0
# rand_intercept_sd = 0.2
# sims = 50
# years = 5
# cores = 2
# response = "binary"
# keep_terms = "treated"
# seed = 999

# Example simulation:
# binary DV, 2 cores, 99 simulations, 5 years
# keep only the treatment variable results
sim_panel(
  n_states = 50,
  n_treated = 10,
  y_control = 0.20,
  y_treat = 0.45,
  y_sd_control = 1,
  y_sd_treated = 1,
  rand_intercept_mean = 0,
  rand_intercept_sd = 0.2,
  sims = 99,
  years = 5,
  cores = 2,
  response = "binary",
  keep_terms = "treated",
  seed = 999
)

# Example simulation:
# binary DV, 2 cores, 10 simulation
# two different assumed treatment effects
# returns a list.
lapply(
  X = c(0.7, 0.9),
  FUN = function(x) {
    sim_panel(
      n_states = 50,
      n_treated = 10,
      y_control = 0.20,
      y_treat = x,
      y_sd_control = 1,
      y_sd_treated = 1,
      rand_intercept_mean = 0,
      rand_intercept_sd = 0.2,
      sims = 99,
      years = 5,
      cores = 2,
      response = "binary",
      keep_terms = "treated",
      seed = 999
    )
  }
)

# Now let's look at power and errors!

# Simulate treatment effects of 1% and 10%
# with 25 states treated and 25 control states
# in each year with 20 years of data and
# with a linear probability model and 1,000
# simulations each.
smallfx <- lapply(
  X = c(0.51, 0.60),
  FUN = function(x) {
    sim_panel(
      n_states = 50,
      n_treated = 25,
      y_control = 0.50,
      y_treat = x,
      y_sd_control = 1,
      y_sd_treated = 1,
      rand_intercept_mean = 0,
      rand_intercept_sd = 0.1,
      sims = 1000,
      years = 20,
      cores = 2,
      response = "binary",
      keep_terms = "treated",
      seed = 999
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
