###

# Run this code if anaemia occurs or if haemoglobin goes below thresholds (different for men and women)

# alpha_i multiplier starts at 1.05. Overall multiplier tends to about 1.1 :
#for (n_actioned in 1:10) { print(1 + 0.1/(2 ^ n_actioned))}
#temp <- 1
#for (n_actioned in 1:10) { temp <- temp * (1 + 0.1/(2 ^ n_actioned)) ; print(temp)}

make_haemoglobin_changes <- function(patient_lme_predictors) {
  
  if (ltc == "CKD") {
    
    patient_lme_predictors$haemoglobin$beta_i <- 0.5 * patient_lme_predictors$haemoglobin$beta_i
    
    n_actioned <- sum(time_to_event_array[, "action_hgb", i_patient] == 1, na.rm = T)
    patient_lme_predictors$haemoglobin$alpha_i <- (1 + 0.1/(2 ^ n_actioned)) * patient_lme_predictors$haemoglobin$alpha_i
    rm(n_actioned)
    
  }
  
  return(list(patient_lme_predictors = patient_lme_predictors))
  
}
