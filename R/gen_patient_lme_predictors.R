###
# Run this once per patient
gen_patient_lme_predictors <- function(lme_predictors_sampled) {
  
  # Sample epsilons
  # Combine alpha and epsilon_alpha_i to get alpha_i
  # Combine beta and epsilon_beta_i to get beta_i
  # This makes it easier to apply lifestyle changes later on
  
  lme_epsilons_sampled <- list()
  patient_lme_predictors <- list()
  
  for (biomarker in names(lme_predictors_sampled)) {
    
    eval(parse(text = paste0("lme_epsilons_sampled$", biomarker,
                             " <- mvrnorm(n = 1, mu = c(0, 0), Sigma = lme_predictors_sampled[[\"",
                             biomarker, "\"]]$epsilons_cov)")))
    
    patient_lme_predictors[[biomarker]]$alpha_i <-
      lme_predictors_sampled[[biomarker]]$alpha + lme_epsilons_sampled[[biomarker]]["epsilon_alpha"]
    
    patient_lme_predictors[[biomarker]]$beta_i <-
      lme_predictors_sampled[[biomarker]]$beta + lme_epsilons_sampled[[biomarker]]["epsilon_beta"]
    
    patient_lme_predictors[[biomarker]]$gamma_1 <- lme_predictors_sampled[[biomarker]]$gamma_1
    patient_lme_predictors[[biomarker]]$gamma_2 <- lme_predictors_sampled[[biomarker]]$gamma_2
    patient_lme_predictors[[biomarker]]$sigma_omega <- lme_predictors_sampled[[biomarker]]$sigma_omega
    
  }

  return(patient_lme_predictors)

}
