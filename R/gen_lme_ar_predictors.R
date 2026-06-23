###
# This function is run once per PSA sample, and output is reused across all patients in that sample
gen_lme_ar_predictors <- function(lme_ar_parameters_list) {
  
  lme_ar_predictors_sampled <- list()
  
  if (ltc == "T2DM") {
    #biomarker_names <- c("eGFR", "hba1c", "potassium", "cholesterol", "hdl", "hdl_ldl", "triglyceride", "ldl", "b12")
    biomarker_names <-   c("hdl_ldl", "ldl")
  } else if (ltc == "HTN") {
    #biomarker_names <- c("eGFR", "hba1c", "potassium", "sodium", "cholesterol", "hdl", "ldl", "triglyceride", "hdl_ldl")
    biomarker_names <- c("eGFR", "cholesterol", "triglyceride")
  } else if (ltc == "CKD") {
    #biomarker_names <- c("eGFR", "hba1c", "potassium", "cholesterol", "hdl", "ldl", "hdl_ldl", "triglyceride", "haemoglobin")
    biomarker_names <- c("eGFR", "potassium", "hdl", "triglyceride")
  }
  
  for (biomarker in biomarker_names) {
    
    lme_ar_parameters_object <- lme_ar_parameters_list[[biomarker]]
    
    fixed_pars <- mvrnorm(n = 1, mu = lme_ar_parameters_object$fixed_coef_means,     #r#
                          Sigma = lme_ar_parameters_object$fixed_coef_cov)
    
    random_pars <- mvrnorm(n = 1, mu = lme_ar_parameters_object$random_coef_means,     #r#
                           Sigma = lme_ar_parameters_object$random_coef_cov)
    names(random_pars) <- paste0(names(random_pars), "_linear")
    
    sigma_a <- exp(random_pars["sigma_a_linear"])
    sigma_omega <- exp(random_pars["sigma_omega_linear"])
    
    lme_ar_predictors_sampled_biomarker <- list(
      alpha = fixed_pars["alpha"],
      theta = fixed_pars["theta"],
      phi = fixed_pars["phi"],
      gamma_1 = fixed_pars["gamma_1"],
      gamma_2 = fixed_pars["gamma_2"],
      sigma_a = sigma_a,
      sigma_omega = sigma_omega
    )
    
    eval(parse(text = paste0("lme_ar_predictors_sampled$", biomarker, " <- lme_ar_predictors_sampled_biomarker")))
    
  }
  return(lme_ar_predictors_sampled)
}
