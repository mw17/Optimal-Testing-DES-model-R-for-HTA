###
# Event survival function
event_surv_func <- function(event_log_rate) {
  rexp(1, rate = exp(event_log_rate)) / 365.25     #r#
}

#############################################

# Dummy event
#dummy_event_log_rate <- -10

# Mean time to dummy_event (in years, using output from event_surv_func)
#dummy_event_sampled <- numeric(10000)
#for (i in 1:10000) { dummy_event_sampled[i] <- event_surv_func(event_log_rate = dummy_event_log_rate[1])}
#mean(dummy_event_sampled)

# Mean time to dummy_event (in years)
#(1 / exp(dummy_event_log_rate[1])) / 365.25

# Proportion who have dummy_event in 1 year
#1 - exp(-365.25 * exp(dummy_event_log_rate[1]))
