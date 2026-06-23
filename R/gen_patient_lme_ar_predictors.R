###
# Run this once per patient
gen_patient_lme_ar_predictors <- function(lme_ar_predictors_sampled) {
  
  # Sample epsilons
  # Combine alpha and epsilon_alpha_i to get alpha_i
  # This makes it easier to apply lifestyle changes later on
  
  lme_ar_epsilons_sampled <- list()
  patient_lme_ar_predictors <- list()
  
  for (biomarker in names(lme_ar_predictors_sampled)) {
    
    eval(parse(text = paste0("lme_ar_epsilons_sampled$", biomarker,
                             " <- rnorm(n = 1, mean = 0, sd = lme_ar_predictors_sampled[[\"",
                             biomarker, "\"]]$sigma_a)")))
    
    patient_lme_ar_predictors[[biomarker]]$alpha_i <-
      lme_ar_predictors_sampled[[biomarker]]$alpha + lme_ar_epsilons_sampled[[biomarker]]
    
    patient_lme_ar_predictors[[biomarker]]$theta <- lme_ar_predictors_sampled[[biomarker]]$theta
    patient_lme_ar_predictors[[biomarker]]$phi <- lme_ar_predictors_sampled[[biomarker]]$phi
    patient_lme_ar_predictors[[biomarker]]$gamma_1 <- lme_ar_predictors_sampled[[biomarker]]$gamma_1
    patient_lme_ar_predictors[[biomarker]]$gamma_2 <- lme_ar_predictors_sampled[[biomarker]]$gamma_2
    patient_lme_ar_predictors[[biomarker]]$sigma_omega <- lme_ar_predictors_sampled[[biomarker]]$sigma_omega
    
  }
  
  return(patient_lme_ar_predictors)
  
}
