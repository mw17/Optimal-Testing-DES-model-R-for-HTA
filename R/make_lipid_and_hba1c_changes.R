###

# Run this code if stroke/ihd/t2dm_control occurs or if lipids/hba1c go beyond thresholds

# HbA1c alpha_i multiplier starts at 0.95. Overall multiplier tends to about 0.9:
#for (n_actioned in 1:10) { print(1 - 0.1/(2 ^ n_actioned)) }
#temp <- 1
#for (n_actioned in 1:10) { temp <- temp * (1 - 0.1/(2 ^ n_actioned)) ; print(temp) }

make_lipid_and_hba1c_changes <- function(patient_lme_predictors, patient_lme_ar_predictors) {
  
  if (ltc == "T2DM") {
    
    n_actioned <- sum(time_to_event_array[, "action_lph", i_patient] == 1, na.rm = T)
    patient_lme_predictors$hba1c$alpha_i <- (1 - 0.1/(2 ^ n_actioned)) * patient_lme_predictors$hba1c$alpha_i
    rm(n_actioned)
    patient_lme_predictors$hba1c$beta_i <- 0
    
    patient_lme_predictors$cholesterol$alpha_i <-
      (survival_coefs_sampled$stroke["cholesterol"] * patient_lme_predictors$cholesterol$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["cholesterol"]
    patient_lme_predictors$cholesterol$beta_i <- 0
    
    patient_lme_predictors$hdl$alpha_i <-
      (survival_coefs_sampled$stroke["hdl"] * patient_lme_predictors$hdl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["hdl"]
    patient_lme_predictors$hdl$beta_i <- 0
    
    patient_lme_ar_predictors$ldl$alpha_i <-
      (survival_coefs_sampled$stroke["ldl"] * patient_lme_ar_predictors$ldl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["ldl"]
    patient_lme_ar_predictors$ldl$phi <- 0
    
    patient_lme_ar_predictors$hdl_ldl$alpha_i <-
      (survival_coefs_sampled$stroke["hdl_ldl"] * patient_lme_ar_predictors$hdl_ldl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["hdl_ldl"]
    patient_lme_ar_predictors$hdl_ldl$phi <- 0
    
  } else if (ltc == "HTN") {
    
    n_actioned <- sum(time_to_event_array[, "action_lph", i_patient] == 1, na.rm = T)
    patient_lme_predictors$hba1c$alpha_i <- (1 - 0.1/(2 ^ n_actioned)) * patient_lme_predictors$hba1c$alpha_i
    rm(n_actioned)
    patient_lme_predictors$hba1c$beta_i <- 0
    
    patient_lme_ar_predictors$cholesterol$alpha_i <-
      (survival_coefs_sampled$stroke["cholesterol"] * patient_lme_ar_predictors$cholesterol$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["cholesterol"]
    patient_lme_ar_predictors$cholesterol$phi <- 0
    
    patient_lme_predictors$hdl$alpha_i <-
      (survival_coefs_sampled$stroke["hdl"] * patient_lme_predictors$hdl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["hdl"]
    patient_lme_predictors$hdl$beta_i <- 0
    
    patient_lme_predictors$ldl$alpha_i <-
      (survival_coefs_sampled$stroke["ldl"] * patient_lme_predictors$ldl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["ldl"]
    patient_lme_predictors$ldl$beta_i <- 0
    
    patient_lme_ar_predictors$triglyceride$alpha_i <-
      (survival_coefs_sampled$stroke["triglyceride"] * patient_lme_ar_predictors$triglyceride$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["triglyceride"]
    patient_lme_ar_predictors$triglyceride$phi <- 0
    
    patient_lme_predictors$hdl_ldl$alpha_i <-
      (survival_coefs_sampled$stroke["hdl_ldl"] * patient_lme_predictors$hdl_ldl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["hdl_ldl"]
    patient_lme_predictors$hdl_ldl$beta_i <- 0
    
  } else if (ltc == "CKD") {
    
    n_actioned <- sum(time_to_event_array[, "action_lph", i_patient] == 1, na.rm = T)
    patient_lme_predictors$hba1c$alpha_i <- (1 - 0.1/(2 ^ n_actioned)) * patient_lme_predictors$hba1c$alpha_i
    rm(n_actioned)
    patient_lme_predictors$hba1c$beta_i <- 0
    
    patient_lme_predictors$cholesterol$alpha_i <-
      (survival_coefs_sampled$stroke["cholesterol"] * patient_lme_predictors$cholesterol$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["cholesterol"]
    patient_lme_predictors$cholesterol$beta_i <- 0
    
    patient_lme_ar_predictors$hdl$alpha_i <-
      (survival_coefs_sampled$stroke["hdl"] * patient_lme_ar_predictors$hdl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["hdl"]
    patient_lme_ar_predictors$hdl$phi <- 0
    
    patient_lme_predictors$ldl$alpha_i <-
      (survival_coefs_sampled$stroke["ldl"] * patient_lme_predictors$ldl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["ldl"]
    patient_lme_predictors$ldl$beta_i <- 0
    
    patient_lme_ar_predictors$triglyceride$alpha_i <-
      (survival_coefs_sampled$stroke["triglyceride"] * patient_lme_ar_predictors$triglyceride$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["triglyceride"]
    patient_lme_ar_predictors$triglyceride$phi <- 0
    
    patient_lme_predictors$hdl_ldl$alpha_i <-
      (survival_coefs_sampled$stroke["hdl_ldl"] * patient_lme_predictors$hdl_ldl$alpha_i +
         log(0.99)/5) / survival_coefs_sampled$stroke["hdl_ldl"]
    patient_lme_predictors$hdl_ldl$beta_i <- 0
    
  }
  
  return(list(patient_lme_predictors = patient_lme_predictors,
              patient_lme_ar_predictors = patient_lme_ar_predictors))
  
}


###############################################################################
###############################################################################

# How to alter risk factors to reduce stroke rate by 1%
# new_rate = 0.99 * old_rate
# exp(new_log_rate) = 0.99 * exp(old_log_rate)
# new_log_rate = log(0.99 * exp(old_log_rate))
# new_log_rate = log(0.99) + log(exp(old_log_rate))
# new_log_rate = log(0.99) + old_log_rate

# So, if we wanted this reduction in stroke rate to be entirely due to change in the risk factor ldl, then
# survival_coefs_sampled_event["ldl"] * risk_factors_sampled$ldl
# needs to increase by log(0.99)
# So, ldl_coef * ldl_rf_new = ldl_coef * ldl_rf_old + log(0.99)
# Hence, ldl_rf_new = (ldl_coef * ldl_rf_old + log(0.99)) / ldl_coef

# So, if we want this reduction in stroke rate to be spread across the five lipid measures then we do
# ldl_rf_new = (ldl_coef * ldl_rf_old + log(0.99)/5) / ldl_coef
# etc for each of the lipids
# And instead of changing the risk factors directly we will change the alphas

###############################################################################
###############################################################################
