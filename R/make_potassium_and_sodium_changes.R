###

# Run this code if renal_injury_specialist occurs or if potassium/sodium go outside thresholds
# Potassium ideally between 3.5 and 5.5
# Sodium ideally between 135 and 145

# Action is to reset intercept to initial value, and a 100% reduction in slope
# So for lme models just set beta_i to 0
# And for lme_ar models, set phi to 0 and reuse t=0 value for the y_prev

# potassium trajectory uses lme model for T2DM and HTN, and lme_ar model for CKD
# sodium trajectory uses lme model for HTN, and doesn't have trajectory for T2DM or CKD

# So for T2DM and HTN, nothing changes when this code is run >= 2 times as slopes are already 0

# NB Need to be careful about when this code is run within the wider model as this code overwrites
# risk_factors_sampled$potassium for CKD. This is fine currently, as any if statements using it are
# run before this code is run.

make_potassium_and_sodium_changes <- function(patient_lme_predictors, patient_lme_ar_predictors) {
  
  if (ltc == "T2DM") {
    patient_lme_predictors$potassium$beta_i <- 0
    
  } else if (ltc == "HTN") {
    patient_lme_predictors$potassium$beta_i <- 0
    
    patient_lme_predictors$sodium$beta_i <- 0
    
  } else if (ltc == "CKD") {
    
    lme_ar_predictors_sampled$potassium$phi <- 0
    risk_factors_sampled$potassium <- risk_factor_array[1, "potassium", i_patient]
    
  }
  
  return(list(patient_lme_predictors = patient_lme_predictors,
              patient_lme_ar_predictors = patient_lme_ar_predictors))
  
}
