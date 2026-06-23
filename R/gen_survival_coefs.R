###
gen_survival_coefs <- function(survival_coef_means, survival_coef_cov) {
  
  survival_coefs_sampled <- list()
  
  for (event in names(survival_coef_means)) {
    
    survival_coef_means[[event]] <- survival_coef_means[[event]][c(colnames(survival_coef_cov[[event]]))]
    
    eval(parse(text = paste0("survival_coefs_sampled$", event,
                             " <- mvrnorm(n = 1, mu = survival_coef_means[[\"", event, "\"]],
                             Sigma = survival_coef_cov[[\"", event, "\"]])")))     #r#
  }
  return(survival_coefs_sampled)
}
