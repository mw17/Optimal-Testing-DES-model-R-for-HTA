###

# Run this code if b12_def occurs or if b12 goes below threshold (i.e. goes below 180)

# alpha_i multiplier starts at 1.05. Overall multiplier tends to about 1.1 :
#for (n_actioned in 1:10) { print(1 + 0.1/(2 ^ n_actioned))}
#temp <- 1
#for (n_actioned in 1:10) { temp <- temp * (1 + 0.1/(2 ^ n_actioned)) ; print(temp)}

make_b12_changes <- function(patient_lme_predictors, patient_lme_ar_predictors) {
  
  if (ltc == "T2DM") {
    
    patient_lme_predictors$b12$beta_i <- 0
    
    n_actioned <- sum(time_to_event_array[, "action_b12", i_patient] == 1, na.rm = T)
    patient_lme_predictors$b12$alpha_i <- (1 + 0.1/(2 ^ n_actioned)) * patient_lme_predictors$b12$alpha_i
    rm(n_actioned)
    
  }
  
  return(list(patient_lme_predictors = patient_lme_predictors,
              patient_lme_ar_predictors = patient_lme_ar_predictors))
  
}