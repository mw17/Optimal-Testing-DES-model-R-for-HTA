###
gen_log_rates <- function(ltc, survival_coefs_sampled, risk_factors_sampled) {
  
  # This give the log rates for a single PSA sample and single patient
  # Patient has risk factors risk_factors_sampled

  # ihd
  survival_coefs_sampled_event <- survival_coefs_sampled$ihd
  ihd_log_rate <-
    survival_coefs_sampled_event["rate"] +
    survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
    survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
    survival_coefs_sampled_event["hba1c"] * risk_factors_sampled$hba1c +
    survival_coefs_sampled_event["cholesterol"] * risk_factors_sampled$cholesterol +
    survival_coefs_sampled_event["hdl"] * risk_factors_sampled$hdl +
    survival_coefs_sampled_event["ldl"] * risk_factors_sampled$ldl +
    survival_coefs_sampled_event["hdl_ldl"] * risk_factors_sampled$hdl_ldl +
    survival_coefs_sampled_event["triglyceride"] * risk_factors_sampled$triglyceride
  
  # t2dm_control
  ### What do we do about "genderother"? Ignore for now? ###
  survival_coefs_sampled_event <- survival_coefs_sampled$t2dm_control
  t2dm_control_log_rate <-
    survival_coefs_sampled_event["rate"] +
    survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
    survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
    survival_coefs_sampled_event["hba1c"] * risk_factors_sampled$hba1c
  
  # stroke
  survival_coefs_sampled_event <- survival_coefs_sampled$stroke
  stroke_log_rate <-
    survival_coefs_sampled_event["rate"] +
    survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
    survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
    survival_coefs_sampled_event["hba1c"] * risk_factors_sampled$hba1c +
    survival_coefs_sampled_event["cholesterol"] * risk_factors_sampled$cholesterol +
    survival_coefs_sampled_event["hdl"] * risk_factors_sampled$hdl +
    survival_coefs_sampled_event["ldl"] * risk_factors_sampled$ldl +
    survival_coefs_sampled_event["hdl_ldl"] * risk_factors_sampled$hdl_ldl +
    survival_coefs_sampled_event["triglyceride"] * risk_factors_sampled$triglyceride
  
  # renal_injury_specialist
  ### What do we do about "genderother"? Ignore for now? ###
  survival_coefs_sampled_event <- survival_coefs_sampled$renal_injury_specialist
  if (ltc == "HTN") {
    renal_injury_specialist_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["sodium"] * risk_factors_sampled$sodium +
      survival_coefs_sampled_event["potassium"] * risk_factors_sampled$potassium
  } else if (ltc == "CKD") {
    renal_injury_specialist_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["sodium"] * sodium_mean[ltc] + # using mean here as no trajectory
      survival_coefs_sampled_event["potassium"] * risk_factors_sampled$potassium
  }
  
  # renal_failure_transplant
  ### What do we do about "genderother"? Ignore for now? ###
  survival_coefs_sampled_event <- survival_coefs_sampled$renal_failure_transplant
  renal_failure_transplant_log_rate <-
    survival_coefs_sampled_event["rate"] +
    survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
    survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
    survival_coefs_sampled_event["eGFR"] * risk_factors_sampled$eGFR
  
  if (ltc == "CKD") {
    # anaemia
    survival_coefs_sampled_event <- survival_coefs_sampled$anaemia
    anaemia_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["haemoglobin"] * risk_factors_sampled$haemoglobin
  }
  
  if (ltc == "T2DM") {
    # b12_def
    survival_coefs_sampled_event <- survival_coefs_sampled$b12_def
    b12_def_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["b12"] * risk_factors_sampled$b12
    
    # amputation
    survival_coefs_sampled_event <- survival_coefs_sampled$amputation
    amputation_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["hba1c"] * risk_factors_sampled$hba1c +
      survival_coefs_sampled_event["cholesterol"] * risk_factors_sampled$cholesterol +
      survival_coefs_sampled_event["hdl"] * risk_factors_sampled$hdl +
      survival_coefs_sampled_event["ldl"] * risk_factors_sampled$ldl +
      survival_coefs_sampled_event["hdl_ldl"] * risk_factors_sampled$hdl_ldl +
      survival_coefs_sampled_event["triglyceride"] * risk_factors_sampled$triglyceride
  }
  
  # death
  survival_coefs_sampled_event <- survival_coefs_sampled$death
  if (ltc == "HTN") {
    death_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["eGFR"] * risk_factors_sampled$eGFR +
      survival_coefs_sampled_event["hba1c"] * risk_factors_sampled$hba1c +
      survival_coefs_sampled_event["cholesterol"] * risk_factors_sampled$cholesterol +
      survival_coefs_sampled_event["hdl"] * risk_factors_sampled$hdl +
      survival_coefs_sampled_event["ldl"] * risk_factors_sampled$ldl +
      survival_coefs_sampled_event["hdl_ldl"] * risk_factors_sampled$hdl_ldl +
      survival_coefs_sampled_event["triglyceride"] * risk_factors_sampled$triglyceride +
      survival_coefs_sampled_event["sodium"] * risk_factors_sampled$sodium +
      survival_coefs_sampled_event["potassium"] * risk_factors_sampled$potassium
  } else {
    death_log_rate <-
      survival_coefs_sampled_event["rate"] +
      survival_coefs_sampled_event["sex"] * risk_factors_sampled$sex +
      survival_coefs_sampled_event["age"] * risk_factors_sampled$age +
      survival_coefs_sampled_event["eGFR"] * risk_factors_sampled$eGFR +
      survival_coefs_sampled_event["hba1c"] * risk_factors_sampled$hba1c +
      survival_coefs_sampled_event["cholesterol"] * risk_factors_sampled$cholesterol +
      survival_coefs_sampled_event["hdl"] * risk_factors_sampled$hdl +
      survival_coefs_sampled_event["ldl"] * risk_factors_sampled$ldl +
      survival_coefs_sampled_event["hdl_ldl"] * risk_factors_sampled$hdl_ldl +
      survival_coefs_sampled_event["triglyceride"] * risk_factors_sampled$triglyceride +
      survival_coefs_sampled_event["sodium"] * sodium_mean[ltc] +  # using mean here as no trajectory
      survival_coefs_sampled_event["potassium"] * risk_factors_sampled$potassium
  }
  # return log rates
  if (ltc == "T2DM") {
    log_rates <- list(ihd = ihd_log_rate,
                      t2dm_control = t2dm_control_log_rate,
                      stroke = stroke_log_rate,
                      renal_failure_transplant = renal_failure_transplant_log_rate,
                      amputation = amputation_log_rate,
                      b12_def = b12_def_log_rate,
                      death = death_log_rate)
  } else if (ltc == "HTN") {
    log_rates <- list(ihd = ihd_log_rate,
                      t2dm_control = t2dm_control_log_rate,
                      stroke = stroke_log_rate,
                      renal_injury_specialist = renal_injury_specialist_log_rate,
                      renal_failure_transplant = renal_failure_transplant_log_rate,
                      death = death_log_rate)
  } else if (ltc =="CKD") {
    log_rates <- list(ihd = ihd_log_rate,
                      t2dm_control = t2dm_control_log_rate,
                      stroke = stroke_log_rate,
                      renal_injury_specialist = renal_injury_specialist_log_rate,
                      renal_failure_transplant = renal_failure_transplant_log_rate,
                      anaemia = anaemia_log_rate,
                      death = death_log_rate)
  }
  return(log_rates)
}
