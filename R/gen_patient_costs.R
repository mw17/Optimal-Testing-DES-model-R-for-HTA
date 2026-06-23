###
gen_patient_costs <- function(time_to_event_array) {
  
  test_costs <- test_costs[test_names]
  
  # Test group costs
  
  test_group_costs <- numeric(0)
  
  for (test_group_number in test_group_numbers) {
    
    # Create test_indicators with test names in same order as test_costs
    test_indicators <- test_times[which(test_times[, "test_group"] == test_group_number),][1, test_names]
    
    # Sum test costs to get overall test group costs, and add nurse cost
    test_group_costs[paste0("test_group_", test_group_number)] <- sum(test_indicators * test_costs) + nurse_cost
    
  }
  
  # Acute event costs
  acute_event_costs <- matrix(0, nrow = n_patients, ncol = n_events,
                              dimnames = list(NULL, event_names))
  
  acute_event_costs[, "ihd"] <- runif(n_patients, min = acute_event_costs_pars$ihd$min,
                                      max = acute_event_costs_pars$ihd$max)
  acute_event_costs[, "t2dm_control"] <- acute_event_costs_pars$t2dm_control
  acute_event_costs[, "stroke"] <- rnorm(n_patients, mean = acute_event_costs_pars$stroke$mean,
                                         sd = acute_event_costs_pars$stroke$sd) 
  if (ltc != "T2DM") { acute_event_costs[, "renal_injury_specialist"] <-
    acute_event_costs_pars$renal_injury_specialist } # Kerr 2014
  acute_event_costs[, "renal_failure_transplant"] <-
    prop_transplant * acute_event_costs_pars$transplant  # transplant
  if (ltc == "CKD") { acute_event_costs[, "anaemia"] <- acute_event_costs_pars$anaemia }
  if (ltc == "T2DM") { acute_event_costs[, "b12_def"] <- acute_event_costs_pars$b12_def }
  if (ltc == "T2DM") { acute_event_costs[, "amputation"] <- acute_event_costs_pars$amputation }
  
  for (test_group_number in test_group_numbers) {
    acute_event_costs[, paste0("test_group_", test_group_number)] <-
      test_group_costs[paste0("test_group_", test_group_number)]
  }
  
  # Chronic event costs
  
  chronic_event_costs <- matrix(0, nrow = n_patients, ncol = n_events,
                                dimnames = list(NULL, event_names))
  
  chronic_event_costs[, "stroke"] <- rnorm(n_patients, mean = chronic_event_costs_pars$stroke$mean,
                                           sd = chronic_event_costs_pars$stroke$sd)
  chronic_event_costs[, "renal_failure_transplant"] <-
    prop_transplant * chronic_event_costs_pars$transplant +
    (1-prop_transplant) * chronic_event_costs_pars$dialysis
  
  if (ltc == "T2DM") { chronic_event_costs[, "b12_def"] <- chronic_event_costs_pars$b12_def }
  
  # Chronic action_lph cost (NB code assumes this cost is fixed)
  chronic_action_lph_cost <- chronic_event_costs_pars$action_lph
  
  # Costs for test result events
  acute_test_result_costs <- list()
  chronic_test_result_costs <- list()
  acute_test_result_costs$renal_failure_transplant <-
    rep(prop_transplant_test * acute_event_costs_pars$transplant, n_patients)  # transplant
  chronic_test_result_costs$renal_failure_transplant <-
    rep(prop_transplant_test * chronic_event_costs_pars$transplant +
    (1-prop_transplant_test) * chronic_event_costs_pars$dialysis, n_patients)
  if (ltc == "CKD") { acute_test_result_costs$anaemia <- acute_event_costs[, "anaemia"] }
  if (ltc == "T2DM") {
    acute_test_result_costs$b12_def <- acute_event_costs[, "b12_def"]
    chronic_test_result_costs$b12_def <- chronic_event_costs[, "b12_def"]
  }
  
  ################################################################
  
  # Now ready to create and fill in the costs array
  
  # Need to use the time to event array, created by competing events code
  # Keep columns: current_time, time_since_previous_event, current_event
  # Add extra columns: acute_cost, chronic_cost, acute_cost_d, chronic_cost_d
  
  costs_array_colnames <- c("current_time", "time_since_previous_event", "current_event",
                            "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12",
                            "acute_cost", "annual_chronic_cost", "chronic_cost", "acute_cost_d", "chronic_cost_d")
  
  costs_array <- array(dim = c(n_times, length(costs_array_colnames), n_patients),
                       dimnames = list(NULL, costs_array_colnames, NULL))
  
  costs_array[, c("current_time", "time_since_previous_event", "current_event",
                  "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12"), ] <-
    time_to_event_array[, c("current_time", "time_since_previous_event", "current_event",
                            "action_lph", "action_egfr", "action_ps", "action_hgb", "action_b12"), ]
  
  # Set acute_cost = 0, annual_chronic_cost = 0, and chronic_cost = 0 for all non-NA current_event rows,
  # as otherwise may have issues later with adding to NAs
  for (i_patient in 1:n_patients) {
    costs_array[which(!is.na(costs_array[, "current_event", i_patient])), "acute_cost", i_patient] <- 0
    costs_array[which(!is.na(costs_array[, "current_event", i_patient])), "annual_chronic_cost", i_patient] <- 0
    costs_array[which(!is.na(costs_array[, "current_event", i_patient])), "chronic_cost", i_patient] <- 0
  }

  for (i_patient in 1:n_patients) {
    
    # First time events occur (used for calculating chronic costs)
    # Create using survival events, then overwrite using actions to capture test result events where appropriate
    first_event_time <- list()
    for (event_number in event_list) {
      first_event_time[names(event_list[event_number])] <-
        costs_array[min(which(costs_array[, "current_event", i_patient] == event_number)),
                    "current_time", i_patient]
    }
    if (ltc == "T2DM") {
      first_event_time$b12_def <-
        costs_array[min(which(costs_array[, "action_b12", i_patient] == 1)), "current_time", i_patient]
    }
    
    # Can only have one of renal_failure_transplant and low eGFR test result events, but costs are different
    # so need to handle test result separately
    first_low_eGFR_test_time <-
      costs_array[min(which(
        (costs_array[, "action_egfr", i_patient] == 1) &
          (grepl("test_group", names(event_list[costs_array[i_time, "current_event", i_patient]]))))
      ), "current_time", i_patient]
    
    
    # action_lph is handled separately as it's triggered by multiple survival and test_result events
    first_action_lph_time <-
      costs_array[min(which(costs_array[, "action_lph", i_patient] == 1)), "current_time", i_patient]
    
    for (i_time in 2:n_times) {
      
      if(is.na(costs_array[i_time, "current_time", i_patient])) break
      
      ### Acute event costs
      
      # Survival events (including test costs, but not test result event costs)
      costs_array[i_time, "acute_cost", i_patient] <-
        acute_event_costs[i_patient, event_list[costs_array[i_time, "current_event", i_patient]]]
      
      # Test result events (have to add these as possible to have more than one at once)
      # Action column used as an indicator here
      if (time_to_event_array[i_time, "current_event", i_patient] %in% grep("test", event_names)) {
        costs_array[i_time, "acute_cost", i_patient] <-
          costs_array[i_time, "acute_cost", i_patient] + # add as test costs already saved here
          costs_array[i_time, "action_egfr", i_patient] * acute_test_result_costs$renal_failure_transplant[i_patient]
        if (ltc == "CKD") {
          costs_array[i_time, "acute_cost", i_patient] <-
            costs_array[i_time, "acute_cost", i_patient] +
            costs_array[i_time, "action_hgb", i_patient] * acute_test_result_costs$anaemia[i_patient]
        }
        else if (ltc == "T2DM") {
          # action_b12 only occurs once for survival event and test result combined, so not accidentally costing it twice here
          costs_array[i_time, "acute_cost", i_patient] <-
            costs_array[i_time, "acute_cost", i_patient] +
            costs_array[i_time, "action_b12", i_patient] * acute_test_result_costs$b12_def[i_patient]
        }
      }
      
      ### Chronic event costs
      
      # First get annual chronic event costs for all events until now including the current event and store this at i_time
      # Then calculate the chronic costs accrued since the previous event using the (i_time - 1) value
      
      for (event_number in event_list) {
        if (!is.na(first_event_time[[event_names[event_number]]])) {
          if (costs_array[i_time, "current_time", i_patient] >= first_event_time[[event_names[event_number]]]) {
            
            # Annual chronic event cost (for events prior to current event)
            costs_array[i_time, "annual_chronic_cost", i_patient] <-
              costs_array[i_time, "annual_chronic_cost", i_patient] +
              chronic_event_costs[i_patient, event_names[event_number]]
            
          }
        }
      }
      
      # Test result costs for low eGFR
      if (!is.na(first_low_eGFR_test_time)) {
        if (costs_array[i_time, "current_time", i_patient] >= first_low_eGFR_test_time) {
          costs_array[i_time, "annual_chronic_cost", i_patient] <-
            costs_array[i_time, "annual_chronic_cost", i_patient] + chronic_test_result_costs$renal_failure_transplant[i_patient]
        }
      }
      # Add action_lph chronic costs
      if (!is.na(first_action_lph_time)) {
        if (costs_array[i_time, "current_time", i_patient] >= first_action_lph_time) {
          costs_array[i_time, "annual_chronic_cost", i_patient] <-
            costs_array[i_time, "annual_chronic_cost", i_patient] + chronic_action_lph_cost
        }
      }
      
      # Time since previous event, multiplied by previous annual chronic event cost
      costs_array[i_time, "chronic_cost", i_patient] <-
        costs_array[i_time, "time_since_previous_event", i_patient] *
        costs_array[i_time - 1, "annual_chronic_cost", i_patient]

      # Discounted costs (see NICE TSD 15)
      
      costs_array[i_time, "acute_cost_d", i_patient] <-
        costs_array[i_time, "acute_cost", i_patient] *
        (1 + discount_rate) ^ -costs_array[i_time, "current_time", i_patient]
      
      idr <- log(1 + discount_rate)
      
      costs_array[i_time, "chronic_cost_d", i_patient] <- 
        costs_array[i_time - 1, "annual_chronic_cost", i_patient] *
        (exp(costs_array[i_time, "current_time", i_patient] *(-idr)) -
           exp(costs_array[i_time - 1, "current_time", i_patient] *(-idr))) /
        (-idr)
      
    } # end loop over times
  } # end loop over patients

  ##########
  
  # Total costs per patient
  patient_costs <- matrix(nrow = n_patients, ncol = 2)
  colnames(patient_costs) <- c("costs", "costs_d")
  for (i_patient in 1:n_patients) {
    patient_costs[i_patient,] <-
      c(sum(costs_array[, c("acute_cost", "chronic_cost"), i_patient], na.rm = T),
        sum(costs_array[, c("acute_cost_d", "chronic_cost_d"), i_patient], na.rm = T))
  }
  
  return(patient_costs)
  
}
