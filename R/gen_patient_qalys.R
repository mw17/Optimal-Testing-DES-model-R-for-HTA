###

gen_patient_qalys <- function(time_to_event_array) {
  
  # Are we taking age into account enough?
  
  # Use the time to event array, created by competing events code
  
  # Keep columns: current_time, time_since_previous_event, current_event,
  # action_lph, action_egfr, action_ps, action_hgb, action_b12
  # Add extra columns: current_event_qaly_loss, annual_qalys, qalys_since_previous_event,
  # current_event_qaly_loss_d, qalys_since_previous_event_d
  
  qalys_array_colnames <- c("current_time", "time_since_previous_event", "current_event",
                            "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12",
                            "current_event_qaly_loss", "annual_qalys", "qalys_since_previous_event",
                            "current_event_qaly_loss_d", "qalys_since_previous_event_d")
  
  qalys_array <- array(dim = c(n_times, length(qalys_array_colnames), n_patients),
                       dimnames = list(NULL, qalys_array_colnames, NULL))
  
  qalys_array[, c("current_time", "time_since_previous_event", "current_event",
                  "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12"), ] <-
    time_to_event_array[, c("current_time", "time_since_previous_event", "current_event",
                            "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12"), ]
  
  ###
  
  # Acute event QALYs (negative means there's a QALY loss)
  
  acute_event_qalys_pars <- acute_event_qalys_pars_ltc[[ltc]]
  
  acute_event_qalys <- matrix(0, nrow = n_patients, ncol = n_events,
                              dimnames = list(NULL, event_names))
  
  if (ltc == "T2DM") {
    acute_event_qalys[, "ihd"] <- rnorm(n_patients, mean = acute_event_qalys_pars$ihd$mean,
                                        sd = acute_event_qalys_pars$ihd$sd)
    acute_event_qalys[, "t2dm_control"] <- rnorm(n_patients, mean = acute_event_qalys_pars$t2dm_control$mean,
                                                 sd = acute_event_qalys_pars$t2dm_control$sd)
    acute_event_qalys[, "stroke"] <- rnorm(n_patients, mean = acute_event_qalys_pars$stroke$mean,
                                           sd = acute_event_qalys_pars$stroke$sd)
    acute_event_qalys[, "b12_def"] <- rnorm(n_patients, mean = acute_event_qalys_pars$b12_def$mean,
                                            sd = acute_event_qalys_pars$b12_def$sd)
    acute_event_qalys[, "amputation"] <- rnorm(n_patients, mean = acute_event_qalys_pars$amputation$mean,
                                               sd = acute_event_qalys_pars$amputation$sd)
    
  } else if (ltc == "HTN") {
    acute_event_qalys[, "ihd"] <- rnorm(n_patients, mean = acute_event_qalys_pars$ihd$mean,
                                        sd = acute_event_qalys_pars$ihd$sd)
    acute_event_qalys[, "t2dm_control"] <- rnorm(n_patients, mean = acute_event_qalys_pars$t2dm_control$mean,
                                                 sd = acute_event_qalys_pars$t2dm_control$sd)
    acute_event_qalys[, "renal_injury_specialist"] <- acute_event_qalys_pars$renal_injury_specialist
    
  } else if (ltc == "CKD") {
    acute_event_qalys[, "ihd"] <- rnorm(n_patients, mean = acute_event_qalys_pars$ihd$mean,
                                        sd = acute_event_qalys_pars$ihd$sd)
    acute_event_qalys[, "t2dm_control"] <- rnorm(n_patients, mean = acute_event_qalys_pars$t2dm_control$mean,
                                                 sd = acute_event_qalys_pars$t2dm_control$sd)
    acute_event_qalys[, "stroke"] <- rnorm(n_patients, mean = acute_event_qalys_pars$stroke$mean,
                                           sd = acute_event_qalys_pars$stroke$sd)
    acute_event_qalys[, "renal_injury_specialist"] <- acute_event_qalys_pars$renal_injury_specialist
    acute_event_qalys[, "anaemia"] <- rnorm(n_patients, mean = acute_event_qalys_pars$anaemia$mean,
                                            sd = acute_event_qalys_pars$anaemia$sd)
  }
  
  # Annual post event QALYs
  # Also including vector for post test eGFR
  
  # SE for T2DM no events was calculated as 0.079*2/3.92 = 0.04
  
  annual_post_event_qalys_pars <- annual_post_event_qalys_pars_ltc[[ltc]]
  
  if (ltc == "T2DM") {
    
    annual_post_event_qalys <- matrix(rep(rnorm(n_patients,
                                                mean = annual_post_event_qalys_pars$no_event$mean,
                                                sd = annual_post_event_qalys_pars$no_event$sd
    ), times = (n_events + 1)),
    nrow = n_patients, ncol = (n_events + 1),
    dimnames = list(NULL, c(event_names, "no_event")))
    
    annual_post_event_qalys[, "ihd"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$ihd$mean,
                                              sd = annual_post_event_qalys_pars$ihd$sd)
    annual_post_event_qalys[, "stroke"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$stroke$mean,
                                                 sd = annual_post_event_qalys_pars$stroke$sd)
    annual_post_event_qalys[, "renal_failure_transplant"] <-
      prop_transplant * rnorm(n_patients, mean = annual_post_event_qalys_pars$transplant$mean,
                              sd = annual_post_event_qalys_pars$transplant$sd) +
      (1 - prop_transplant) * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$dialysis$shape1,
                                    shape2 = annual_post_event_qalys_pars$dialysis$shape2)
    annual_post_event_qalys[, "amputation"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$amputation$mean,
                                                     sd = annual_post_event_qalys_pars$amputation$sd)
    annual_post_event_qalys[, "death"] <- annual_post_event_qalys_pars$death
    
    annual_post_eGFR_test_qalys <-
      prop_transplant_test * rnorm(n_patients, mean = annual_post_event_qalys_pars$transplant$mean,
                                   sd = annual_post_event_qalys_pars$transplant$sd) +
      (1 - prop_transplant_test) * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$dialysis$shape1,
                                         shape2 = annual_post_event_qalys_pars$dialysis$shape2)
    
  } else if (ltc == "HTN") {
    
    annual_post_event_qalys <- matrix(rep(pmin(rnorm(n_patients,
                                                     mean = annual_post_event_qalys_pars$no_event$mean,
                                                     sd = annual_post_event_qalys_pars$no_event$sd
    ), 1), times = (n_events + 1)),
    nrow = n_patients, ncol = (n_events + 1),
    dimnames = list(NULL, c(event_names, "no_event")))
    
    annual_post_event_qalys[, "ihd"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$ihd$mean,
                                              sd = annual_post_event_qalys_pars$ihd$sd)
    annual_post_event_qalys[, "stroke"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$stroke$mean,
                                                 sd = annual_post_event_qalys_pars$stroke$sd)
    annual_post_event_qalys[, "renal_failure_transplant"] <-
      prop_transplant * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$transplant$shape1,
                              shape2 = annual_post_event_qalys_pars$transplant$shape2) +
      (1 - prop_transplant) * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$dialysis$shape1,
                                    shape2 = annual_post_event_qalys_pars$dialysis$shape2)
    annual_post_event_qalys[, "death"] <- annual_post_event_qalys_pars$death
    
    annual_post_eGFR_test_qalys <-
      prop_transplant_test * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$transplant$shape1,
                                   shape2 = annual_post_event_qalys_pars$transplant$shape2) +
      (1 - prop_transplant_test) * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$dialysis$shape1,
                                         shape2 = annual_post_event_qalys_pars$dialysis$shape2)
    
  } else if (ltc == "CKD") {
    
    annual_post_event_qalys <- matrix(rep(rnorm(n_patients, mean = annual_post_event_qalys_pars$no_event$mean,
                                                sd = annual_post_event_qalys_pars$no_event$sd), times = (n_events + 1)),
                                      nrow = n_patients, ncol = (n_events + 1),
                                      dimnames = list(NULL, c(event_names, "no_event")))
    
    annual_post_event_qalys[, "ihd"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$ihd$mean,
                                              sd = annual_post_event_qalys_pars$ihd$sd)
    annual_post_event_qalys[, "stroke"] <- rnorm(n_patients, mean = annual_post_event_qalys_pars$stroke$mean,
                                                 sd = annual_post_event_qalys_pars$stroke$sd)
    annual_post_event_qalys[, "renal_failure_transplant"] <-
      prop_transplant * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$transplant$shape1,
                              shape2 = annual_post_event_qalys_pars$transplant$shape2) +
      (1 - prop_transplant) * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$dialysis$shape1,
                                    shape2 = annual_post_event_qalys_pars$dialysis$shape2)
    annual_post_event_qalys[, "death"] <- annual_post_event_qalys_pars$death
    
    annual_post_eGFR_test_qalys <-
      prop_transplant_test * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$transplant$shape1,
                                   shape2 = annual_post_event_qalys_pars$transplant$shape2) +
      (1 - prop_transplant_test) * rbeta(n_patients, shape1 = annual_post_event_qalys_pars$dialysis$shape1,
                                         shape2 = annual_post_event_qalys_pars$dialysis$shape2)
    
  }
  
  #######
  
  #qalys_array[1:5, , i_patient]
  
  for (i_patient in 1:n_patients) {
    
    # When i_time = 1
    qalys_array[1, "annual_qalys", i_patient] <- annual_post_event_qalys[i_patient, "no_event"]
    
    for (i_time in 2:n_times) {
      
      # Acute event qalys - from survival events
      qalys_array[i_time, "current_event_qaly_loss", i_patient] <-
        acute_event_qalys[i_patient, event_list[qalys_array[i_time, "current_event", i_patient]]]
      
      # Acute test result qalys - from test result events
      #if (grepl("test_group", names(event_list[qalys_array[i_time, "current_event", i_patient]]))) {
      #  
      #  if (!is.na(qalys_array[i_time, "action_egfr", i_patient])) {
      #    qalys_array[i_time, "current_event_qaly_loss", i_patient] <-
      #      acute_test_result_qalys[i_patient, "eGFR"]
      #  }
      #  # No acute event qalys for other actions
      #  
      #}
      
      # Get qalys since previous event by multiplying time since previous event by
      # previous annual qalys
      qalys_array[i_time, "qalys_since_previous_event", i_patient] <-
        qalys_array[i_time, "time_since_previous_event", i_patient] * 
        qalys_array[i_time - 1, "annual_qalys", i_patient]
      
      # Set the annual post event qalys to be the minimum of the existing ones and the new ones
      if (!is.na(event_list[qalys_array[i_time, "current_event", i_patient]])) {
        
        # survival events
        qalys_array[i_time, "annual_qalys", i_patient] <-
          min(qalys_array[i_time - 1, "annual_qalys", i_patient],
              annual_post_event_qalys[i_patient, event_list[qalys_array[i_time, "current_event", i_patient]]])
        
        # test result events
        # only test results events with an annual post event qaly change is low eGFR
        if (grepl("test_group", names(event_list[qalys_array[i_time, "current_event", i_patient]]))) {
          if (qalys_array[i_time, "action_egfr", i_patient] == 1) {
            qalys_array[i_time, "annual_qalys", i_patient] <-
              min(qalys_array[i_time - 1, "annual_qalys", i_patient],
                  annual_post_eGFR_test_qalys[i_patient])
          }
        }
        
      }
      
      # Discounted QALYs (see NICE TSD 15)
      
      qalys_array[i_time, "current_event_qaly_loss_d", i_patient] <-
        qalys_array[i_time, "current_event_qaly_loss", i_patient] *
        (1 + discount_rate) ^ -qalys_array[i_time, "current_time", i_patient]
      
      idr <- log(1 + discount_rate)
      
      qalys_array[i_time, "qalys_since_previous_event_d", i_patient] <- 
        qalys_array[i_time - 1, "annual_qalys", i_patient] *
        (exp(qalys_array[i_time, "current_time", i_patient] *(-idr)) -
           exp(qalys_array[i_time - 1, "current_time", i_patient] *(-idr))) /
        (-idr)
      
    } # end loop over times
    
  } # end loop over patients
  
  ##########
  
  # Total QALYs per patient
  patient_qalys <- matrix(nrow = n_patients, ncol = 2)
  colnames(patient_qalys) <- c("qalys", "qalys_d")
  for (i_patient in 1:n_patients) {
    patient_qalys[i_patient,] <-
      c(sum(qalys_array[, c("current_event_qaly_loss", "qalys_since_previous_event"), i_patient], na.rm = T),
        sum(qalys_array[, c("current_event_qaly_loss_d", "qalys_since_previous_event_d"), i_patient], na.rm = T))
  }
  
  return(patient_qalys)
  
}
