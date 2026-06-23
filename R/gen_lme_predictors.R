###
# This function is run once per PSA sample, and output is reused across all patients in that sample
gen_lme_predictors <- function(lme_parameters_list) {
  
  lme_predictors_sampled <- list()
  
  if (ltc == "T2DM") {
    #biomarker_names <- c("eGFR", "hba1c", "potassium", "cholesterol", "hdl", "hdl_ldl", "triglyceride", "ldl", "b12")
    biomarker_names <-   c("eGFR", "hba1c", "potassium", "cholesterol", "hdl", "b12")
  } else if (ltc == "HTN") {
    #biomarker_names <- c("eGFR", "hba1c", "potassium", "sodium", "cholesterol", "hdl", "ldl", "triglyceride", "hdl_ldl")
    biomarker_names <- c("hba1c", "potassium", "sodium", "hdl", "ldl", "hdl_ldl")
  } else if (ltc == "CKD") {
    #biomarker_names <- c("eGFR", "hba1c", "potassium", "cholesterol", "hdl", "ldl", "hdl_ldl", "triglyceride", "haemoglobin")
    biomarker_names <- c("hba1c", "cholesterol", "ldl", "hdl_ldl", "haemoglobin")
  }
  
  for (biomarker in biomarker_names) {
    
    lme_parameters_object <- lme_parameters_list[[biomarker]]
    
    fixed_pars <- mvrnorm(n = 1, mu = lme_parameters_object$fixed_coef_means,      #r#
                          Sigma = lme_parameters_object$fixed_coef_cov)
    
    random_pars <- mvrnorm(n = 1, mu = lme_parameters_object$random_coef_means,     #r#
                           Sigma = lme_parameters_object$random_coef_cov)
    names(random_pars) <- paste0(names(random_pars), "_linear")
    
    sigma_a <- exp(random_pars["sigma_a_linear"])
    sigma_b <- exp(random_pars["sigma_b_linear"])
    rho <- 2 * (exp(random_pars["rho_linear"]) / (1 + exp(random_pars["rho_linear"]))) - 1
    sigma_omega <- exp(random_pars["sigma_omega_linear"])
    
    epsilons_cov <- matrix(data = c(sigma_a^2, rho*sigma_a*sigma_b, rho*sigma_a*sigma_b, sigma_b^2),
                           nrow = 2, 
                           dimnames = list(c("epsilon_alpha", "epsilon_beta"),
                                           c("epsilon_alpha", "epsilon_beta")))
    
    lme_predictors_sampled_biomarker <- list(
      alpha = fixed_pars["alpha"],
      beta = fixed_pars["beta"],
      gamma_1 = fixed_pars["gamma_1"],
      gamma_2 = fixed_pars["gamma_2"],
      sigma_a = sigma_a,
      sigma_b = sigma_b,
      rho = rho,
      sigma_omega = sigma_omega,
      #sigma_a_linear = random_pars[, "sigma_a_linear"],
      #sigma_b_linear = random_pars[, "sigma_b_linear"],
      #rho_linear = random_pars[, "rho_linear"],
      #sigma_omega_linear = random_pars[, "sigma_omega_linear"],
      epsilons_cov = epsilons_cov
    )
    
    eval(parse(text = paste0("lme_predictors_sampled$", biomarker, " <- lme_predictors_sampled_biomarker")))
    
  }
  return(lme_predictors_sampled)
}
