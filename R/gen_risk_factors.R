###
gen_risk_factors <- function(ltc,
                             patient_lme_predictors,
                             patient_lme_ar_predictors,
                             age, sex, t, y_prev, delta_t) {
  
  risk_factors_sampled <- y_prev
  # For e.g. T2DM, the sodium, haemoglobin and triglyceride values aren't updated
  
  risk_factors_sampled$t <- t
  
  for (biomarker in names(patient_lme_predictors)) {
    eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                             " <- gen_lme_y(patient_lme_predictors = patient_lme_predictors$", biomarker, ", ",
                             "age = age, sex = sex, t = t)")))
  }
  
  for (biomarker in names(patient_lme_ar_predictors)) {
    eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                             " <- gen_lme_ar_y(patient_lme_ar_predictors = patient_lme_ar_predictors$", biomarker, ", ",
                             "age = age, sex = sex, ",
                             "y_prev = y_prev$", biomarker, ", ",
                             "delta_t = delta_t)")))
  }
  
  #for (biomarker in names(fixed_ar_predictors_sampled)) {
  #  eval(parse(text = paste0("risk_factors_sampled$", biomarker,
  #                           " <- gen_fixed_ar_y(fixed_ar_predictors_sampled = fixed_ar_predictors_sampled$", biomarker, ", ",
  #                           "age = age, sex = sex, ",
  #                           "y_prev = y_prev$", biomarker, ", ",
  #                           "delta_t = delta_t)")))
  #  
  #}
  
  # Force risk factor values to be within a biologically plausible range
  for (biomarker in names(risk_factors_sampled)[-(1:3)]) {
    
    # Don't apply risk_factor_limits min value to eGFR trajectory as trajectory
    # will instead be reset by corrective action of having kidney transplant, so
    # use min = 0 instead
    if (biomarker == "eGFR") {
      if (risk_factors_sampled$eGFR < 0) { risk_factors_sampled$eGFR <- 0 }
    } else {
      if ( eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                    " < risk_factor_limits$", biomarker, "$Min"))) ) {
        eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                 " <- risk_factor_limits$", biomarker, "$Min")))
      }
    }
    
    if ( eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                  " > risk_factor_limits$", biomarker, "$Max"))) ) {
      eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                               " <- risk_factor_limits$", biomarker, "$Max")))
    }
    
  }
  
  return(risk_factors_sampled)
}
