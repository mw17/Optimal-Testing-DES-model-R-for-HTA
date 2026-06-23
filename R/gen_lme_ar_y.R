###
# Run this for a single patient once per biomarker and time point
gen_lme_ar_y <- function(patient_lme_ar_predictors,
                            age, sex, y_prev, delta_t) {
  
  # Convert delta_t to days
  delta_t_days <- delta_t * 365.25
  
  # omegas are freshly sampled for each t
  omega_t <- rnorm(n = 1, mean = 0, sd = patient_lme_ar_predictors$sigma_omega)     #r#
  
  lme_ar_y <-
    patient_lme_ar_predictors$alpha_i +
    patient_lme_ar_predictors$theta * y_prev +
    patient_lme_ar_predictors$phi * log(delta_t_days) +
    patient_lme_ar_predictors$gamma_1 * age +
    patient_lme_ar_predictors$gamma_2 * sex +
    omega_t
  
  return(unname(lme_ar_y))
  
}
