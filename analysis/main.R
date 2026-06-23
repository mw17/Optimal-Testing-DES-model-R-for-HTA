###
rm(list=ls())

###

# ltc is one of: "T2DM", "HTN", "CKD"
#ltc <- "XXX"

# Test strategy
#test_strategy_name <- "YYY"

# PSA (used in filenames for outputs)
#psa <- ZZZ

# n_patients
#n_patients <- 1000

##################################

ltc <- "CKD"
test_strategy_name <- "B"
psa <- 1
n_patients <- 1000

##################################

# set seed
# Now using same seed for each strategy, but still using different seeds for each PSA and ltc
ltc_seed <- c("T2DM" = 1, "HTN" = 2, "CKD" = 3)
seed <- as.numeric(paste0(ltc_seed[ltc],
                          1,
                          sprintf('%05d', psa))) # assumes PSAs don't go above 99999
seed
set.seed(seed)

#######

#library(devtools)
#library(pkgload)

library(openxlsx2)
library(MASS)
library(condMVNorm)

###

# Source all files (except "unused" files) in R folder
files <- list.files("R")
files <- files[-grep("unused", files)]
for (file in files) { source(paste0("R/",file)) }


###

# Load input data
load("data/input_data_main.RData")
load("data/input_data_misc.RData")
load("data/input_data_utilities.RData")
load("data/input_data_costs.RData")
###

# Select relevant ltc input data from input_data_main
risk_factor_limits <- risk_factor_limits_ltc[[ltc]]
lme_parameters_list <- lme_parameters_list_ltc[[ltc]]
lme_ar_parameters_list <- lme_ar_parameters_list_ltc[[ltc]]
survival_coef_means <- survival_coef_means_ltc[[ltc]]
survival_coef_cov <- survival_coef_cov_ltc[[ltc]]

###
# Run this once per PSA sample
lme_predictors_sampled <- gen_lme_predictors(lme_parameters_list)
lme_ar_predictors_sampled <- gen_lme_ar_predictors(lme_ar_parameters_list)

# Run this once per PSA sample
survival_coefs_sampled <- gen_survival_coefs(survival_coef_means, survival_coef_cov)

#############
### Run the following to check that code is working and to get risk factor names ###
### Then tidy up but keep risk factor names ###

# Run this once per patient
patient_lme_predictors <- gen_patient_lme_predictors(lme_predictors_sampled)
patient_lme_ar_predictors <- gen_patient_lme_ar_predictors(lme_ar_predictors_sampled)

# Run this once per patient and time point
risk_factors_sampled <- gen_risk_factors_t0(ltc,
                                            patient_lme_predictors,
                                            age = 65, sex = 1)
log_rates_sampled <- gen_log_rates(ltc, survival_coefs_sampled, risk_factors_sampled)

# t > 0, e.g, t = 1
risk_factors_sampled <- gen_risk_factors(ltc,
                                         patient_lme_predictors,
                                         patient_lme_ar_predictors,
                                         age = 65, sex = 1, t = 1,
                                         y_prev = risk_factors_sampled, delta_t = 1)
log_rates_sampled <- gen_log_rates(ltc, survival_coefs_sampled, risk_factors_sampled)

risk_factor_names <- names(risk_factors_sampled)

rm(patient_lme_predictors, patient_lme_ar_predictors,
   risk_factors_sampled, log_rates_sampled)
#############

### Tests

# Get vector of test frequencies
xlsx_path_tests <- "data/testing_strategies.xlsx"
test_strategies <- wb_to_df(xlsx_path_tests, sheet = ltc, row_names = 1, col_names = 1)
test_frequencies <- test_strategies[test_strategy_name, ]

# Test names
test_names <- colnames(test_frequencies)
n_tests <- length(test_names)

# Create matrix of test times
colnames_test_times <- c("current_time", test_names, "test_group")
test_time_vector <- seq(from = 0.25, to = max_model_years, by = 0.25)
test_times <- matrix(nrow = length(test_time_vector), ncol = length(colnames_test_times),
                     dimnames = list(NULL, colnames_test_times))

test_times[, "current_time"] <- test_time_vector
for (test in test_names) {
  if (test_frequencies[test] == 0) { test_times[, test] <- 0
  } else if (test_frequencies[test] == 0.25) { test_times[, test] <- 1
  } else if (test_frequencies[test] == 0.5) { test_times[, test] <- rep(c(0, 1), 1000)[1:length(test_time_vector)]
  } else if (test_frequencies[test] == 1) { test_times[, test] <- rep(c(0, 0, 0, 1), 1000)[1:length(test_time_vector)]
  } else if (test_frequencies[test] == 2) { test_times[, test] <- rep(c(0, 0, 0, 0, 0, 0, 0, 1), 1000)[1:length(test_time_vector)]
  } else if (test_frequencies[test] == 3) { test_times[, test] <- rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1), 1000)[1:length(test_time_vector)]
  }
}

# Only keep rows in which tests occur
test_times <- test_times[which(rowSums(test_times[, grep("_test$", colnames(test_times))]) > 0), ]

# Need to group tests as can't have more than one event occurring at same time
# There can be overlap between groups
# The test group name uses the binary of the test indicators
test_times[, "test_group"] <- strtoi(apply(test_times[, test_names], 1, paste, collapse = ""), base = 2)

test_group_numbers <- unique(test_times[, "test_group"])
test_group_names <- paste0("test_group_", test_group_numbers)

###

# Create named vector of events
if (ltc == "T2DM") {
  event_names <- c(
    # Events
    "ihd", "t2dm_control", "stroke", "renal_failure_transplant", "b12_def", "amputation", "death",
    # Test groups
    test_group_names
  )
} else if (ltc == "HTN") {
  event_names <- c(
    # Events
    "ihd", "t2dm_control", "stroke", "renal_injury_specialist", "renal_failure_transplant", "death",
    # Test groups
    test_group_names
  )
} else if (ltc == "CKD") {
  event_names <- c(
    # Events
    "ihd", "t2dm_control", "stroke", "renal_injury_specialist", "renal_failure_transplant", "anaemia", "death",
    # Test groups
    test_group_names
  )
}

n_events <- length(event_names)

event_list <- 1:n_events
names(event_list) <- event_names

time_to_event_names <- paste0("time_to_", event_names)

array_colnames <- c("current_time", "time_since_previous_event", "current_event",
                    "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12",
                    time_to_event_names)

######################################################################
######################################################################

# Max number of events (inc. tests) (+ 1 for time 0)
n_times <- 4 * max_model_years + 1 # allows average 4 events (inc. tests) per year for max_model_years

time_to_event_array <- array(dim = c(n_times, length(array_colnames), n_patients),
                             dimnames = list(NULL, array_colnames, NULL))

# Create empty array to store risk factor values over time
risk_factor_array <- array(dim = c(n_times, length(risk_factor_names), n_patients),
                           dimnames = list(NULL, risk_factor_names, NULL))

### Initiate array values

time_to_event_array[1, "current_time", ] <- 0
time_to_event_array[1, "time_since_previous_event", ] <- 0    # setting as 0 so code doesn't break
time_to_event_array[1, "current_event", ] <- NA     # no event happens at time 0

# Test groups
for (test_group_number in test_group_numbers) {
  eval(parse(text = paste0(
    "time_to_event_array[1, \"time_to_test_group_", test_group_number,
    "\", ] <- min(test_times[which(test_times[, \"test_group\"] == ", test_group_number,
    "), \"current_time\"])"
  )))
}

# Sample age and sex for all patients
# Set seed again here, always using same set seed
seed <- 50
set.seed(seed)
age_n_patients <- rlnorm(n_patients, meanlog = baseline_age[[ltc]]$meanlog,
                         sdlog = baseline_age[[ltc]]$sdlog)
sex_n_patients <- rbinom(n_patients, size = 1, prob = baseline_sex_prob[[ltc]])

### Main loop starts here ###
for (i_patient in 1:n_patients) {
  
  age <- age_n_patients[i_patient]
  sex <- sex_n_patients[i_patient]
  
  # lme_predictors_sampled and lme_ar_predictors_sampled are same for all patients in a PSA sample
  
  # set seed again here, with a unique seed for each patient
  # this is so that these values are the same for each testing strategy
  seed <- as.numeric(paste0(sprintf('%05d', psa), # assumes PSAs don't go above 99999
                            sprintf('%04d', i_patient))) # assumes patients don't go above 9999
  set.seed(seed)
  
  # Run this once per patient
  patient_lme_predictors <- gen_patient_lme_predictors(lme_predictors_sampled)
  patient_lme_ar_predictors <- gen_patient_lme_ar_predictors(lme_ar_predictors_sampled)
  
  # Run this once per patient to initiate values
  risk_factors_sampled <- gen_risk_factors_t0(ltc,
                                              patient_lme_predictors,
                                              age = age, sex = sex)
  
  # Save risk factors to array
  risk_factor_array[1, , i_patient] <- unlist(risk_factors_sampled)
  
  # Calculate log rates
  log_rates_sampled <- gen_log_rates(ltc, survival_coefs_sampled, risk_factors_sampled)
  
  # Sample times to events
  time_to_event_array[1, "time_to_ihd", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$ihd)
  time_to_event_array[1, "time_to_t2dm_control", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$t2dm_control)
  time_to_event_array[1, "time_to_stroke", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$stroke)
  if (ltc %in% c("HTN", "CKD")) { time_to_event_array[1, "time_to_renal_injury_specialist", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$renal_injury_specialist) }
  time_to_event_array[1, "time_to_renal_failure_transplant", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$renal_failure_transplant)
  if (ltc == "CKD") { time_to_event_array[1, "time_to_anaemia", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$anaemia) }
  if (ltc == "T2DM") {
    time_to_event_array[1, "time_to_b12_def", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$b12_def)
    time_to_event_array[1, "time_to_amputation", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$amputation)
  }
  time_to_event_array[1, "time_to_death", i_patient] <- event_surv_func(event_log_rate = log_rates_sampled$death)
  
  for (i_time in 2 : dim(time_to_event_array)[1]) {
    
    # get current event index
    time_to_event_array[i_time, "current_event", i_patient] <-
      which.min(time_to_event_array[i_time - 1, time_to_event_names, i_patient])
    
    # update time_since_previous_event
    time_to_event_array[i_time, "time_since_previous_event", i_patient] <-
      time_to_event_array[i_time - 1,
                          time_to_event_names[time_to_event_array[i_time, "current_event", i_patient]],
                          i_patient]
    
    # update current_time
    time_to_event_array[i_time, "current_time", i_patient] <-
      time_to_event_array[i_time - 1, "current_time", i_patient] +
      time_to_event_array[i_time, "time_since_previous_event", i_patient]
    
    # stop this inner loop if current event is death
    if (time_to_event_array[i_time, "current_event", i_patient] ==
        which(event_names == "death")
    ) break
    
    # stop this inner loop if current time is greater than max_model_years
    # (and also set current event to NA since can't happen)
    if (time_to_event_array[i_time, "current_time", i_patient] > max_model_years
    ) {time_to_event_array[i_time, "current_event", i_patient] <- NA
      break}
    
    # Set actions to 0 so code doesn't break from NAs
    time_to_event_array[i_time, c("action_lph", "action_egfr", "action_ps",
                                  "action_hgb", "action_b12"), i_patient] <- 0
    
    # Save risk_factors at previous time point to y_prev, before any are overwritten
    y_prev <- risk_factors_sampled
    
    # Resample risk factors to get values for current time point
    # This needs to be done now as test results need to be known before any
    # interventions (e.g. changes in lme alpha) can be made
    risk_factors_sampled <- gen_risk_factors(ltc,
                                             patient_lme_predictors,
                                             patient_lme_ar_predictors,
                                             age = age, sex = sex,
                                             t = unname(time_to_event_array[i_time, "current_time", i_patient]),
                                             y_prev = y_prev,
                                             delta_t = time_to_event_array[i_time, "time_since_previous_event", i_patient])
    
    # Update underlying trajectory models in response to current event
    # Rather than rewriting model code to say e.g. 0.5 * alpha_i, am instead halving the alpha_i value as
    # gives the same result and is easier to code that way
    # Need to select the correct trajectory model (lme, lme_ar, fixed_ar) for each ltc/biomarker combination
    
    # Events that aren't tests
    # Trajectory modelling method is specific to each ltc
    # Also, included events are specific to each ltc
    
    if (!(time_to_event_array[i_time, "current_event", i_patient] %in% grep("test", event_names))) {
      
      if (ltc == "T2DM") {
        if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "ihd")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "t2dm_control")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "stroke")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "renal_failure_transplant")) {
          time_to_event_array[i_time, "action_egfr", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "b12_def")) {
          time_to_event_array[i_time, "action_b12", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "amputation")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        }
      } else if (ltc == "HTN") {
        if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "ihd")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "t2dm_control")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "stroke")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "renal_injury_specialist")) {
          time_to_event_array[i_time, "action_ps", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "renal_failure_transplant")) {
          time_to_event_array[i_time, "action_egfr", i_patient] <- 1
          
        }
      } else if (ltc == "CKD") {
        if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "ihd")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "t2dm_control")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "stroke")) {
          time_to_event_array[i_time, "action_lph", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "renal_injury_specialist")) {
          time_to_event_array[i_time, "action_ps", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "renal_failure_transplant")) {
          time_to_event_array[i_time, "action_egfr", i_patient] <- 1
          
        } else if (time_to_event_array[i_time, "current_event", i_patient] == which(event_names == "anaemia")) {
          time_to_event_array[i_time, "action_hgb", i_patient] <- 1
          
        }
      }
      
    } else {
      # Events that are tests
      # Only alter trajectories for risk factors included in current tests
      # Trajectory modelling method is specific to each ltc
      
      if (ltc == "T2DM") {
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "hba1c_test"] == 1) {
          if (risk_factors_sampled$hba1c > hba1c_threshold) {
            time_to_event_array[i_time, "action_lph", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "egfr_test"] == 1) {
          if (risk_factors_sampled$eGFR < egfr_threshold) {
            if ( sum(time_to_event_array[, "action_egfr", i_patient] == 1, na.rm = T) == 0 ) {
              time_to_event_array[i_time, "action_egfr", i_patient] <- 1
            }
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "lipid_test"] == 1) {
          # No trajectory for T2DM triglycerides
          if (risk_factors_sampled$ldl > ldl_threshold) {
            time_to_event_array[i_time, "action_lph", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "potassium_test"] == 1) {
          if ((risk_factors_sampled$potassium < potassium_threshold_lower) | (risk_factors_sampled$potassium > potassium_threshold_upper)) {
            time_to_event_array[i_time, "action_ps", i_patient] <- 1
          }
        }
        # can only be diagnosed with b12 once as treatment is permanent in model
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "vit_b12_test"] == 1) {
          if (risk_factors_sampled$b12 < b12_threshold) {
            if ( sum(time_to_event_array[, "action_b12", i_patient] == 1, na.rm = T) == 0 ) {
              time_to_event_array[i_time, "action_b12", i_patient] <- 1
            }
          }
        }
        
      } else if (ltc == "HTN") {
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "hba1c_test"] == 1) {
          if (risk_factors_sampled$hba1c > hba1c_threshold) {
            time_to_event_array[i_time, "action_lph", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "egfr_test"] == 1) {
          if (risk_factors_sampled$eGFR < egfr_threshold) {
            if ( sum(time_to_event_array[, "action_egfr", i_patient] == 1, na.rm = T) == 0 ) {
              time_to_event_array[i_time, "action_egfr", i_patient] <- 1
            }
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "lipid_test"] == 1) {
          if ((risk_factors_sampled$ldl > ldl_threshold) | (risk_factors_sampled$triglyceride > triglyceride_threshold)) {
            time_to_event_array[i_time, "action_lph", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "potassium_test"] == 1) {
          if ((risk_factors_sampled$potassium < potassium_threshold_lower) | (risk_factors_sampled$potassium > potassium_threshold_upper)) {
            time_to_event_array[i_time, "action_ps", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "sodium_test"] == 1) {
          if ((risk_factors_sampled$sodium < sodium_threshold_lower) | (risk_factors_sampled$sodium > sodium_threshold_upper)) {
            time_to_event_array[i_time, "action_ps", i_patient] <- 1
          }
        }
        
      } else if (ltc == "CKD") {
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "hba1c_test"] == 1) {
          if (risk_factors_sampled$hba1c > hba1c_threshold) {
            time_to_event_array[i_time, "action_lph", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "egfr_test"] == 1) {
          if (risk_factors_sampled$eGFR < egfr_threshold) {
            if ( sum(time_to_event_array[, "action_egfr", i_patient] == 1, na.rm = T) == 0 ) {
              time_to_event_array[i_time, "action_egfr", i_patient] <- 1
            }
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "lipid_test"] == 1) {
          if ((risk_factors_sampled$ldl > ldl_threshold) | (risk_factors_sampled$triglyceride > triglyceride_threshold)) {
            time_to_event_array[i_time, "action_lph", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "potassium_test"] == 1) {
          if ((risk_factors_sampled$potassium < potassium_threshold_lower) | (risk_factors_sampled$potassium > potassium_threshold_upper)) {
            time_to_event_array[i_time, "action_ps", i_patient] <- 1
          }
        }
        if (test_times[which(test_times[, "current_time"] == time_to_event_array[i_time, "current_time", i_patient]),
                       "fbc_test"] == 1) {
          if ((risk_factors_sampled$haemoglobin < haemoglobin_threshold_female & sex == 0) |
              (risk_factors_sampled$haemoglobin < haemoglobin_threshold_male & sex == 1)) {
            time_to_event_array[i_time, "action_hgb", i_patient] <- 1
          }
        }
        
      }
      
    }
    
    ##### Apply trajectory changes based on values for action_lph etc.
    if (time_to_event_array[i_time, "action_lph", i_patient] == 1) {
      lipid_and_hba1c_changes_output <- make_lipid_and_hba1c_changes(patient_lme_predictors, patient_lme_ar_predictors)
      patient_lme_predictors <- lipid_and_hba1c_changes_output$patient_lme_predictors
      patient_lme_ar_predictors <- lipid_and_hba1c_changes_output$patient_lme_ar_predictors
    }
    if (time_to_event_array[i_time, "action_egfr", i_patient] == 1) {
      eGFR_changes_output <- make_eGFR_changes(patient_lme_predictors, patient_lme_ar_predictors)
      patient_lme_predictors <- eGFR_changes_output$patient_lme_predictors
      patient_lme_ar_predictors <- eGFR_changes_output$patient_lme_ar_predictors
    }
    if (time_to_event_array[i_time, "action_ps", i_patient] == 1) {
      potassium_and_sodium_changes_output <- make_potassium_and_sodium_changes(patient_lme_predictors, patient_lme_ar_predictors)
      patient_lme_predictors <- potassium_and_sodium_changes_output$patient_lme_predictors
      patient_lme_ar_predictors <- potassium_and_sodium_changes_output$patient_lme_ar_predictors
    }
    if (time_to_event_array[i_time, "action_hgb", i_patient] == 1) {
      haemoglobin_changes_output <- make_haemoglobin_changes(patient_lme_predictors)
      patient_lme_predictors <- haemoglobin_changes_output$patient_lme_predictors
    }
    if (time_to_event_array[i_time, "action_b12", i_patient] == 1) {
      b12_changes_output <- make_b12_changes(patient_lme_predictors, patient_lme_ar_predictors)
      patient_lme_predictors <- b12_changes_output$patient_lme_predictors
      patient_lme_ar_predictors <- b12_changes_output$patient_lme_ar_predictors
    }
    
    # Some risk_factors_sampled values (eGFR and potassium) may have been updated by actions above since they
    # were originally sampled at this i_time
    
    # Save risk factors to array
    risk_factor_array[i_time, , i_patient] <- unlist(risk_factors_sampled)
    
    # Update times to non-test events
    log_rates_sampled <- gen_log_rates(ltc, survival_coefs_sampled, risk_factors_sampled)
    
    # Can only be diagnosed with ihd once
    if (which(event_names == "ihd") %in% time_to_event_array[, "current_event", i_patient]) {
      time_to_event_array[i_time, "time_to_ihd", i_patient] <- Inf
    } else {
      time_to_event_array[i_time, "time_to_ihd", i_patient] <-
        event_surv_func(event_log_rate = log_rates_sampled$ihd)
    }
    
    time_to_event_array[i_time, "time_to_t2dm_control", i_patient] <-
      event_surv_func(event_log_rate = log_rates_sampled$t2dm_control)
    
    time_to_event_array[i_time, "time_to_stroke", i_patient] <-
      event_surv_func(event_log_rate = log_rates_sampled$stroke)
    
    # Can only have AKI hospitalisation once since it has a lifetime QALY loss
    if (ltc %in% c("HTN", "CKD")) { 
      if (which(event_names == "renal_injury_specialist") %in% time_to_event_array[, "current_event", i_patient]) {
        time_to_event_array[i_time, "time_to_renal_injury_specialist", i_patient] <- Inf
      } else {
        time_to_event_array[i_time, "time_to_renal_injury_specialist", i_patient] <-
          event_surv_func(event_log_rate = log_rates_sampled$renal_injury_specialist)
      }
    }
    
    # Can only have end stage renal disease once, then has permanent impact on costs and QALYs
    # Also can't have ESRD event if have already had low eGFR test result.
    if (sum(time_to_event_array[, "action_egfr", i_patient] == 1, na.rm = T) == 0) {
      time_to_event_array[i_time, "time_to_renal_failure_transplant", i_patient] <-
        event_surv_func(event_log_rate = log_rates_sampled$renal_failure_transplant)
    } else {
      time_to_event_array[i_time, "time_to_renal_failure_transplant", i_patient] <- Inf
    }
    
    if (ltc == "CKD") {
      # Can't have anaemia event if less than 6 months since previous anaemia (event or test result) as
      # will still be on treatment
      # Using 0.50001 rather than 0.5, as otherwise can have anaemia event clashing with test times if testing every six months
      if (
        (time_to_event_array[i_time, "current_time", i_patient] - 
         max(time_to_event_array[
           time_to_event_array[, "action_hgb", i_patient] == 1,
           "current_time", i_patient], na.rm = T)) < 0.50001
      ) {
        time_to_event_array[i_time, "time_to_anaemia", i_patient] <-
          max(event_surv_func(event_log_rate = log_rates_sampled$anaemia),
              0.50001 - (
                time_to_event_array[i_time, "current_time", i_patient] - 
                  max(time_to_event_array[
                    time_to_event_array[, "action_hgb", i_patient] == 1,
                    "current_time", i_patient], na.rm = T)
              ))
      } else {
        time_to_event_array[i_time, "time_to_anaemia", i_patient] <-
          event_surv_func(event_log_rate = log_rates_sampled$anaemia)
      }
    }
    
    if (ltc == "T2DM") {
      
      # Can only have b12_def event once (from survival equation or from test result) as remain
      # on treatment permanently
      if (sum(time_to_event_array[, "action_b12", i_patient] == 1, na.rm = T) == 0) {
        time_to_event_array[i_time, "time_to_b12_def", i_patient] <-
          event_surv_func(event_log_rate = log_rates_sampled$b12_def)
      } else {
        time_to_event_array[i_time, "time_to_b12_def", i_patient] <- Inf
      }
      
      # Limit amputations to once, has permanent effect on QALYs
      if (which(event_names == "amputation") %in% time_to_event_array[, "current_event", i_patient]) {
        time_to_event_array[i_time, "time_to_amputation", i_patient] <- Inf
      } else {
        time_to_event_array[i_time, "time_to_amputation", i_patient] <-
          event_surv_func(event_log_rate = log_rates_sampled$amputation)
      }
      
      ### Limit AKI hospitalisation to once as the acute QALY loss represents lifetime QALY loss??? ###
      
      
    }
    
    time_to_event_array[i_time, "time_to_death", i_patient] <-
      event_surv_func(event_log_rate = log_rates_sampled$death)
    
    # Update times to tests
    for (test_group_number in test_group_numbers) {
      eval(parse(text = paste0(
        "time_to_event_array[i_time, \"time_to_test_group_",
        test_group_number,
        "\", i_patient] <- 
        min(test_times[
          which((test_times[, \"current_time\"] > time_to_event_array[i_time, \"current_time\", i_patient]) &
                  (test_times[, \"test_group\"] == ",
        test_group_number,
        ")),
          \"current_time\"]) -
        time_to_event_array[i_time, \"current_time\", i_patient]"
      )))
    }
    
  } # End loop over times
  
  # Remove values where time is greater than max_model_years
  if (time_to_event_array[i_time, "current_time", i_patient] > max_model_years) {
    time_to_event_array[i_time, c("current_event", "time_since_previous_event", "current_time"), i_patient] <- NA
  }
  
} # End loop over patients

###########

i_patient

###

seed <- as.numeric(paste0(ltc_seed[ltc],
                          2,
                          sprintf('%05d', psa))) # assumes PSAs don't go above 99999
set.seed(seed)
patient_qalys <- gen_patient_qalys(time_to_event_array)

seed <- as.numeric(paste0(ltc_seed[ltc],
                          3,
                          sprintf('%05d', psa))) # assumes PSAs don't go above 99999
set.seed(seed)
patient_costs <- gen_patient_costs(time_to_event_array)

patient_results <- cbind(patient_qalys, patient_costs)
fname <- paste0("results/patient_qalys_and_costs_20260203_", ltc, "_", test_strategy_name, "_psa", psa, ".txt")
#write.table(patient_results, file = fname, row.names = F, col.names = F)

# summary info
fname2 <- paste0("results/event_summary_20260203_", ltc, "_", test_strategy_name, "_psa", psa, ".txt")
#sink(fname2)
#source("analysis/trajectory_and_event_summary.R")
#:sink()
