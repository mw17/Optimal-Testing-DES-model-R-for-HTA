###
# Run this for a single patient once per biomarker and time point
gen_lme_y <- function(patient_lme_predictors,
                         age, sex, t) {
  
  # Convert t to days
  t_days <- t * 365.25
  
  # omegas are freshly sampled for each t
  omega_t <- rnorm(n = 1, mean = 0, sd = patient_lme_predictors$sigma_omega)     #r#
  
  lme_y <-
    patient_lme_predictors$alpha_i +
    t_days * patient_lme_predictors$beta_i +
    patient_lme_predictors$gamma_1 * age +
    patient_lme_predictors$gamma_2 * sex +
    omega_t
  
  return(unname(lme_y))
  
}
