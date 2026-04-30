library(tidyverse)
# Tr 1 = Treatment 1 Dry; # Tr 2 = Treatment 2 Mid; # Tr 3 = Treatment 3 Wet

## 1. Low COI first

### Initial
LowInit <- 1
LowAfter_Tr1_gt <- 1
LowAfter_Tr2_gt <- 1
LowAfter_Tr3_gt <- 2

###  Runs
#### Exponential distribution of results -> rescaled to COI level's range -> split by binary classes
range_Low_Tr1 <- c(40, 0)
range_Low_Tr2 <- c(45, 5)
range_Low_Tr3 <- c(60, 15)

treatments <- c("Tr1", "Tr2", "Tr3")

running_sim_bin <- function(treatment, samp) {
  range_run <- get(paste0("range_Low_", treatment))
  
  run_Low <- rexp(samp, 1/2)
  run_Low_ord <- run_Low[order(-run_Low)]
  run_Low_scaled <- (range_run[1] - range_run[2]) * (run_Low_ord - min(run_Low_ord)) / (max(run_Low_ord) - min(run_Low_ord)) + range_run[2]
  
  print(paste("Assigning", treatment_run))
  return(run_Low_scaled)
}

# Setting the number of samples
nsamps <- c(5,10,15,20,100)

sim_results_view = data.frame(); sim_results_plot = data.frame()

for (samp in nsamps) {
  for (treatment in treatments){
    print(paste("Running", treatment, "with sample number", samp))
    treatment_run <- paste0("run_Low_", treatment,"_scaled")
    
    assign(treatment_run, running_sim_bin(treatment, samp))
    
    ### Binary (1 = R, 2 = S) significance
    LowAfter_sim <- ifelse (get(treatment_run) <= 10, 1, 2)
    
    ### Comparing results
    LowAfter_res <- ifelse (LowAfter_sim == get(paste0("LowAfter_", treatment,"_gt")), "True", "False") # Where True = the simulated result matches the ground-truth
    
    freq_table <- tidyr::spread(as.data.frame(table(LowAfter_res)), LowAfter_res, Freq)
    plot_table <- as.data.frame(table(LowAfter_res))
    if (!("False" %in% unique(LowAfter_res))) {
      freq_table$"False" <- 0
      new_row <- data.frame(LowAfter_res = "False", Freq = 0); plot_table <- rbind(plot_table, new_row)
    }
    print(cbind(treatment, freq_table))
    
    sim_results_view <- rbind(sim_results_view, cbind(treatment, freq_table, samp))
    sim_results_plot <- rbind(sim_results_plot, cbind(treatment, plot_table, samp))
  }
}

print(sim_results)

sim_results_plot <- sim_results_plot %>%
  group_by(samp) %>%
  mutate(Proportional = Freq / sum(Freq))

ggplot(sim_results_plot, aes(x=treatment, y=Proportional, fill=LowAfter_res)) + 
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~samp) +
  theme_minimal()

plot(density(run_Low_Tr1_scaled), xlab="Values", ylab="Density", col="red")
lines(density(run_Low_Tr2_scaled), xlab="Values", ylab="Density", col="blue")
lines(density(run_Low_Tr3_scaled), xlab="Values", ylab="Density", col="green")

# Quantitative

