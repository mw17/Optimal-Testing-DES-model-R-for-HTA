###

# Run this code if renal_failure/transplant occurs or if eGFR goes below 15

# Only allowed to happen once

# eGFR trajectory uses lme model for T2DM, and lme_ar model for HTN/CKD

make_eGFR_changes <- function(patient_lme_predictors, patient_lme_ar_predictors) {
  
  if (ltc == "T2DM") {
    
    patient_lme_predictors$eGFR$alpha_i <- 1.1 * patient_lme_predictors$eGFR$alpha_i
    
    patient_lme_predictors$eGFR$beta_i <- 0.5 * patient_lme_predictors$eGFR$beta_i
    
  } else {
    
    risk_factors_sampled$eGFR <- max(risk_factors_sampled$eGFR,
                                     rnorm(1, mean = new_eGFR$mean, sd = new_eGFR$sd))
    
  }
  
  return(list(patient_lme_predictors = patient_lme_predictors,
              patient_lme_ar_predictors = patient_lme_ar_predictors))
  
}
