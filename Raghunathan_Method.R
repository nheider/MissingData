
# Problem: How do you impute variables that are never observed together but their 
# conditional expectation is of interest? The dataset also contains a
# one extra var that is a surrogate of one of the vars of interest     
#
# Model: 
#  +---+---+---+
#  | X |   | Z |
#  |   | Y | Z |
#  +---+---+---+
#
# 
# Raghunathans Approach (Missing Data Analysis in Practice (2015): 156-9):

# set-up ----
library(tidyverse)
library(rstanarm)
library(MASS)

set.seed(123)

# functions ----
# generate full dataset and measurement error pattern of scenario b (p. 156) 

generate_data <- function(){
  
  sigma <- matrix(c(1, 1, 1.7, 1, 2, 1.5, 1.7, 1.5, 4), 3, 3)
  mu <- c(0,1,2)
  full_dataset <- mvrnorm(n = 1000, mu, sigma)
  
  colnames(full_dataset) <- c("X", "Y", "W")
  
  main_study <- full_dataset[1:900,]
  sub_study <- full_dataset[901:1000,]
  
  #delete all X Values of the main study 
  main_study[, 1] <- NA 
  
  # Scenario 2
  sub_study2 <- sub_study
  sub_study2[,2] <- NA
  scenario_2 <- rbind(main_study, sub_study2)
  
  output <- list(full_dataset, scenario_2)
  
  output 
} 

raghunathan <- function(scenario_b){
  main_study <- scenario_b[1:900,]
  sub_study <- scenario_b[901:1000,]

  first_blr <- stan_glm(X ~ W, data = sub_study, refresh = 0)

  a_0 <- median(as.data.frame(first_blr)$"(Intercept)") 
  a_1 <- median(as.data.frame(first_blr)$W) 
  tau_squared <- median(as.data.frame(first_blr)$sigma)^2 
    
  second_blr <- stan_glm(Y ~ W, data = main_study, refresh = 0)
 
  b_0 <- median(as.data.frame(second_blr)$"(Intercept)") 
  b_1 <- median(as.data.frame(second_blr)$W, n_draws)
  sigma_squared <- median(as.data.frame(second_blr)$sigma)^2 

# Based on p. 157 equation 8.1
  
  B_1 <- b_1 / a_1
  Sigma <- sigma_squared - B_1 * tau_squared
  B_0 <- b_0 - B_1 * a_0

# Generate values from (X|Y, W)

  psi_squared <- (B_1^2 / Sigma + 1/tau_squared)^(-1)

  Y_main_study <- main_study$Y
  W_main_study <- main_study$W

  muX <- numeric()
  
  for(i in 1:900){
    muX[i] <- psi_squared * (B_1 * (Y_main_study[i] - B_0) / Sigma + 
                     (a_0 + a_1*W_main_study[i]) / tau_squared)
  } 

  x_draws <- rnorm(900, muX, psi_squared)
  main_study$X <- x_draws
  
  main_study
} 

compare <- function(data_true, data_comp){
  
  reg_true <- lm(Y ~ X, data_true)
  reg_imp <- lm(Y ~ X, data_comp)
  
  output <- c(summary(reg_true)$r.squared, 
              summary(reg_imp)$r.squared, 
              reg_true$coefficient[2], 
              reg_imp$coefficient[2], 
              confint(reg_true)[2,1], 
              confint(reg_true)[2,2], 
              confint(reg_imp)[2,1], 
              confint(reg_imp)[2,2])
  output 
}

calculate_statistics <- function(comparison){
  B1_bias <- comparison[4] -  comparison[3]
  len_confint_true <- comparison[6] - comparison[5]
  len_confint_imp <- comparison[8] - comparison[7]
  
  confint_hit <- comparison[7] < comparison[3] & comparison[3] < comparison[8]
  
  output <- c(comparison[1], 
              comparison[2], 
              B1_bias, 
              len_confint_true, 
              len_confint_imp,
              confint_hit)
}

stats <- matrix(ncol = 6, nrow = 0)
colnames(stats)<- c("R_2_true", 
                    "R_2_imp", 
                    "B1_bias", 
                    "len_confint_true", 
                    "len_confint_imp", 
                    "confint_hit")

# simulation study ----
progress <- numeric()
j <- 1

for(i in 1:250){
  data_true <- as.data.frame(generate_data()[[1]])[1:900,]
  
  data_imp <- as.data.frame(generate_data()[[2]])
  data_comp <- raghunathan(data_imp)
  
  stats <- rbind(stats, calculate_statistics(compare(data_true, data_comp))) 
                 
  # progress indicator 
  progress[j] <- round(j/250, digits = 2)*100
  if(j >= 2){
    if(progress[j] > progress[j-1]){
      print(paste(progress[j], "% Completed", sep = ""))
    } 
    
  } else { 
    print(paste(progress[j], "% Completed", sep = ""))
  }               
  j <- j + 1
}

mean(as.data.frame(stats)$B1_bias) # very good 

mean(as.data.frame(stats)$R_2_true)
mean(as.data.frame(stats)$R_2_imp) # very good 

mean(as.data.frame(stats)$len_confint_true)
mean(as.data.frame(stats)$len_confint_imp) # very close to the true value 

mean(as.data.frame(stats)$confint_hit) # very bad 

# -> The confidence interval is larger than the true one but only slightly
# -> Overall it seems much to small for large ammount of imputed data

# -> The coverage is really bad (probably because of the small convidence interval)

# -> A possible fix would be to add some noise 

