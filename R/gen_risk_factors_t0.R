###
gen_risk_factors_t0 <- function(ltc,
                                patient_lme_predictors,
                                age, sex) {
  
  # This function assumes that get_risk_factor_limits has already been run for this ltc
  
  t <- 0
  risk_factors_sampled <- list(age = age, sex = sex, t = t)
  
  for (biomarker in names(patient_lme_predictors)) {
    
    eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                             " <- gen_lme_y(patient_lme_predictors = patient_lme_predictors$", biomarker, ", ",
                             "age = age, sex = sex, t = t)")))
    
    # Force risk factor values to be within a biologically plausible range
    if ( eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                  " < risk_factor_limits$", biomarker, "$Min"))) ) {
      eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                               " <- risk_factor_limits$", biomarker, "$Min")))
    }
    if ( eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                  " > risk_factor_limits$", biomarker, "$Max"))) ) {
      eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                               " <- risk_factor_limits$", biomarker, "$Max")))
    }
    
  }
  
  # Autoregressive trajectory starting values
  rf_means <- rf_means_ltc[[ltc]]
  rf_cov <- rf_cov_ltc[[ltc]]
  rf_names <- names(rf_means)
  
  # names of risk factors that are on log scale
  log_rf_names <- rf_names[grep("_log", rf_names)]
  log_rf_names_unlogged <- gsub("_log", "", log_rf_names)
  
  # names of all risk factors, with _log removed if it exists
  rf_names <- gsub("_log", "", rf_names)
  
  names(rf_means) <- rf_names
  colnames(rf_cov) <- rownames(rf_cov) <- rf_names
  
  non_ar_biomarker_names <- c("age", names(patient_lme_predictors))
  given_ind <- which(rf_names %in% non_ar_biomarker_names)
  
  ar_biomarker_names <- setdiff(rf_names, non_ar_biomarker_names)
  dep_ind <- which(rf_names %in% ar_biomarker_names)
  
  x_given <- unlist(risk_factors_sampled)
  x_given <- x_given[non_ar_biomarker_names]
  x_given <- x_given[names(rf_means[given_ind])]    # ensuring order is correct
  
  for (x_name in names(x_given)) {
    if (x_name %in% log_rf_names_unlogged) {
      x_given[x_name] <- log(x_given[x_name])
    }
  }
  
  ar_y_t0_sampled <- rcmvnorm(n = 1, mean = rf_means, sigma = rf_cov,
                              dependent.ind = dep_ind, given.ind = given_ind,
                              X.given = x_given)     #r#
  ar_y_t0_sampled <- as.list(ar_y_t0_sampled)
  names(ar_y_t0_sampled) <- names(rf_means)[dep_ind]
  
  for (biomarker in names(ar_y_t0_sampled)) {
    
    if (biomarker %in% log_rf_names_unlogged) {
      ar_y_t0_sampled[biomarker] <- exp(unlist(ar_y_t0_sampled[biomarker]))
    }
    eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                             " <- ar_y_t0_sampled[[\"", biomarker, "\"]]")))
  }
  
  # Force risk factor values to be within a biologically plausible range
  for (biomarker in ar_biomarker_names) {
    if ( eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                  " < risk_factor_limits$", biomarker, "$Min"))) ) {
      eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                               " <- risk_factor_limits$", biomarker, "$Min")))
    }
    if ( eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                                  " > risk_factor_limits$", biomarker, "$Max"))) ) {
      eval(parse(text = paste0("risk_factors_sampled$", biomarker,
                               " <- risk_factor_limits$", biomarker, "$Max")))
    }
  }
  
  return(risk_factors_sampled)
}
