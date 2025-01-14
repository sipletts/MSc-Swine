
#############################################
#                                           #                                                    
#            SIV Spillover Model            #
#             Pig + Human + Env             #
#                                           #
#            Status: working :D             #
#                                           #
#           Updated: 17.05.24               #
#                                           #                                                 
#############################################
packageVersion("deSolve")

#Define SEIR model
#Summary: Deterministic SEIR model with 2 latent states and 11 infectious states & environmental compartment & human compartment

#Update 08.05: sedimentation rate corrected for pig height
#Update 08.05: inactivation rate added to V instead of DR model 
#Update 14.05: corrected ci for quanta (0.8 ipv 1.25)
#Update 14.05: corrected for incorrect P_inf direct (was missing E-06 whoops)
#Update 15.05: added "High" contact rate based on the number of face contacts per day as a function of duration in barn 
#Update 15.05: corrected "Low" contact rate mistake (was multiplying by #hmans)
#Update 15.05: cleaned up code and added descriptions, and indicated which parameters are adjustable
#Update 17.05: added a simplified plot for pig dynamics

library(deSolve)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(knitr)

#MODEL
spillover_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Extracting compartments from the state vector
    Sp <- state[1]
    Epk <- state[2]  
    Epl <- state[3]  
    Ipk <- state[4:9] 
    Ipl <- state[10:14]  
    Rp <- state[15]
    
    Sh <- state[16]
    Ih <- state[17]
    Rh <- state[18]
    
    V <- state [19]
    
    # Calculate the total population
    Np <- Sp + Epk + Epl + sum(Ipk) + sum(Ipl) + Rp
    Nh <- Sh + Ih + Rh
    Nv <- V
    
    # Calculate the derivative of each compartment
    #PIGGIES
    dSp <- - beta * Sp * (sum(Ipk) + sum(Ipl)) / Np
    dEpk <- 0.475 * beta * Sp * (sum(Ipk) + sum(Ipl)) / Np - sigma1 * Epk
    dEpl <- (1-0.475) * beta * Sp * (sum(Ipk) + sum(Ipl)) / Np - sigma2 * Epl
    
    dIpk1 <- sigma1 * Epk - gamma * Ipk1
    dIpk2 <- gamma * Ipk1 - gamma * Ipk2
    dIpk3 <- gamma * Ipk2 - gamma * Ipk3
    dIpk4 <- gamma * Ipk3 - gamma * Ipk4
    dIpk5 <- gamma * Ipk4 - gamma * Ipk5
    dIpk6 <- gamma * Ipk5 - gamma * Ipk6
    
    dIpl1 <- sigma2 * Epl - gamma * Ipl1
    dIpl2 <- gamma * Ipl1 - gamma * Ipl2
    dIpl3 <- gamma * Ipl2 - gamma * Ipl3
    dIpl4 <- gamma * Ipl3 - gamma * Ipl4
    dIpl5 <- gamma * Ipl4 - gamma * Ipl5
    
    dRp <- gamma * Ipk6 + gamma * Ipl5
    
    #HUMANS    
    ######### ADJUSTABLE PARAMETER 1: CONTACT RATE HIGH VS LOW (turn off/on) ##########   
    
    #Close-contact DR Model pt1 (High contact rate):
    mucous_contacts <- 0.123   #per minute contact with mouth and/or nose
    contacts_by_human <- mucous_contacts * 25 #!!!match with duration!!! 
    
    #Close-contact DR Model pt1 (Low contact rate):
    #contacts_by_human <- 1
    
    probability_infected_contact <- (sum(Ipk) + sum(Ipl)) / Np
    contacts_infected_pig <- contacts_by_human * probability_infected_contact
    
    ################################################################################
    
    #Close-contact DR Model pt2 (P_inf): 
    #P_inf: separate R code: QMRA_2024 calculation
    
    #Close-contact beta: Beta = P_inf * contact rate
    beta_pk1 <- 1.5936E-06 * contacts_infected_pig  
    beta_pk2 <- 4.6581E-06 * contacts_infected_pig
    beta_pk3 <- 4.5791E-06 * contacts_infected_pig
    beta_pk4 <- 3.4521E-06 * contacts_infected_pig
    beta_pk5 <- 2.3251E-06 * contacts_infected_pig
    beta_pk6 <- 1.1982E-06 * contacts_infected_pig
    beta_pl1 <- 1.3760E-06 * contacts_infected_pig
    beta_pl2 <- 2.5228E-06 * contacts_infected_pig
    beta_pl3 <- 2.5228E-06 * contacts_infected_pig
    beta_pl4 <- 2.5228E-06 * contacts_infected_pig
    beta_pl5 <- 1.4360E-06 * contacts_infected_pig
    
    ######### ADJUSTABLE PARAMETER 2: DURATION HIGH VS LOW (turn 25mins/420mins) ##########  
    
    #Air DR Model pt1 (contact rate ie. inhales)
    #Volume of air inhaled per breath 
    inhale_volume  <- 0.035 #m3/min
    bpm <- 16 #breaths per min (ave 12-20)
    duration <- 420
    contact <- duration * bpm
    
    #Air DR Model pt2 (P_inf)
    volume_per_breath <- inhale_volume / bpm
    dose_air <- V * volume_per_breath #convert quanta back to TCID50 (4/5q = 1 TCID50; dose = TCID50 per breath!
    
    #Beta-Poisson DR model 
    a <- 0.581
    N_50 <- 9.45E+5
    P_human <- 1 - ((1 + dose_air * ((2^(1/a) - 1) / N_50))^-a)
    
    #Optional: Exponential DR Model -> similar results to beta-poisson
    #a <- 9.0E-8   #susceptibility/rate at which the response increases with the dose
    #P_human <- 1 - exp(-a * dose_air)
    
    #Air beta: P_inf * contact rate
    beta_v <- contact * P_human
    
    dSh <- -(beta_pk1 * Ipk1 + beta_pk2 * Ipk2 + beta_pk3 * Ipk3 +
               beta_pk4 * Ipk4 + beta_pk5 * Ipk5 + beta_pk6 * Ipk6 +
               beta_pl1 * Ipl1 + beta_pl2 * Ipl2 + beta_pl3 * Ipl3 +
               beta_pl4 * Ipl4 + beta_pl5 * Ipl5) * Sh - beta_v * Sh * V 
    dIh <- (beta_pk1 * Ipk1 + beta_pk2 * Ipk2 + beta_pk3 * Ipk3 +
              beta_pk4 * Ipk4 + beta_pk5 * Ipk5 + beta_pk6 * Ipk6 +
              beta_pl1 * Ipl1 + beta_pl2 * Ipl2 + beta_pl3 * Ipl3 +
              beta_pl4 * Ipl4 + beta_pl5 * Ipl5) * Sh - delta * Ih + beta_v * Sh * V
    dRh <- delta * Ih 
    
    #ENVIRONMENT 
    ######### ADJUSTABLE PARAMETER 3: ERq HIGH VS LOW (turn off/on) #######################
    #High:
    #dV <- q1 * (Ipk1 + Ipl1) + q2 * (Ipk2 + Ipl2) + q3 * (Ipk3 + Ipl3) + q3 * (Ipk4 + Ipl4) + q2 * (Ipk5 + Ipl5) + q1 * Ipk6 - s * V - v * V
    
    #Low: 
    dV <- (qa * Ipk1 + q1 *Ipl1 + qb * Ipk2 + q2 * Ipl2 + qc * Ipk3 + q3 * Ipl3 + qd * Ipk4 + q4 * Ipl4 + qe * Ipk5 + q5 * Ipl5 + qf * Ipk6) * 1.25 - s * V - v * V - mu * V
    
    ########################################################################################
    
    return(list(c(dSp, dEpk, dEpl, dIpk1, dIpk2, dIpk3, dIpk4, dIpk5, dIpk6,
                  dIpl1, dIpl2, dIpl3, dIpl4, dIpl5, dRp, dSh, dIh, dRh, dV)))
    
  })
}

# Initial state: 1 infectious (Ipk) pig, 2999 susceptible pigs, 10 employees
initial_state <- c(Sp = 2019, Epk = 0, Epl = 0, Ipk = c(1, 0, 0, 0, 0, 0), 
                   Ipl = c(0, 0, 0, 0, 0), Rp = 0, Sh = 10, Ih = 0, Rh = 0, V = 0)

#Parameters

########## ADJUSTABLE PARAMETERS: B (0.25 or 1.88) and EQr values ###################
parameters <- list(
  #pig
  beta = 1.88, 
  sigma1 = (1/2.48), 
  sigma2 = (1/1.16), 
  gamma = 1, 
  Np = 2020,
  #Human
  Nh = 10,
  delta = 1/8,
  #Env
  s = 21.12,
  v= 121.2,
  mu = 0.16,
  #High EQr
  #q1 = 97, q2 = 157.9, q3 = 194.65
  #Low EQr
  qa = 0.122, q1 = 0.085, qb = 17.97, q2 = 0.553, qc = 15.8, q3 = 0.553, qd = 2.516, q4 = 0.553, qe = 0.401, q5 = 0.094, qf = 0.064
) 

#Time vector (30 days: farrowing, 40 days: weaning, 100 days: finishing, Total: 170 days)
times <- seq(0, 100, by = 1)  

# Solve the differential equations
model <- ode(y = initial_state, times = times, func = spillover_model, parms = parameters)



#PLOTTING PER COMPARTMENT
model_df <- as.data.frame(model)
model_long <- model_df %>% pivot_longer(cols = -time, names_to = "compartment", values_to = "value")

#Humans 
human_compartments <- c("Sh", "Ih", "Rh")
ggplot(data = filter(model_long, compartment %in% human_compartments), aes(x = time, y = value, color = compartment)) +
  geom_line() +
  labs(title = "Human Dynamics", x = "time (days)", y = "# humans") +
  theme_minimal() 


#Pigs 
pig_compartments <- c("Sp", "Epk", "Epl", "Ipl1", "Ipl2", "Ipl3", "Ipl4", "Ipl5", "Ipk1", "Ipk2", "Ipk3", "Ipk4", "Ipk5", "Ipk6", "Rp")
ggplot(data = filter(model_long, compartment %in% pig_compartments), aes(x = time, y = value, color = compartment)) +
  geom_line() +
  labs(title = "Pig Dynamics - high beta", x = "time (days)", y = "# pigs") +
  theme_minimal()

#Pig : Total I 
model_df$Total_I <- rowSums(model_df[ , grepl("Ipk|Ipl", names(model_df))])
plot_data <- model_df %>%
  select(time, Total_I, Sp, Rp) %>%
  pivot_longer(cols = -time, names_to = "compartment", values_to = "value")
ggplot(plot_data, aes(x = time, y = value, color = compartment)) +
  geom_line() +
  labs(title = "Simplified Pig Dynamics",
       x = "Time (days)", y = "# pigs",
       color = "Compartment") +
  theme_minimal()

#Env
environment_compartments <- c("V")
ggplot(data = filter(model_long, compartment %in% environment_compartments), aes(x = time, y = value, color = compartment)) +
  geom_line() +
  labs(title = "Environment Dynamics - high beta, low q", x = "time (days)", y = "SIV in air (TCID50/m3)") +
  theme_minimal()

#V and total I pigs
model_df$Total_I <- rowSums(model_df[ , grepl("Ipk|Ipl", names(model_df))])
plot_data <- model_df %>%
  select(time, Total_I, V) %>%
  pivot_longer(cols = -time, names_to = "compartment", values_to = "value")
ggplot(plot_data, aes(x = time, y = value, color = compartment)) +
  geom_line() +
  labs(title = "Overlay of outbreak and SIV in air",
       x = "Time (days)", y = "Total I pigs / SIV in air (quanta)",
       color = "Compartment") +
  theme_minimal()


#Check timesteps 
print(model)


#Plotting for ppt: thicker lines 
#Pigs
library(ggplot2)
pig_compartments <- c("Sp", "Epk", "Epl", "Ipl1", "Ipl2", "Ipl3", "Ipl4", "Ipl5", "Ipk1", "Ipk2", "Ipk3", "Ipk4", "Ipk5", "Ipk6", "Rp")
ggplot(data = filter(model_long, compartment %in% pig_compartments), aes(x = time, y = value, color = compartment)) +
  geom_line(size = 1.5) +
  labs(title = "Pig Dynamics - high beta", x = "time (days)", y = "# pigs") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 20), axis.title = element_text(size = 15))

#Env
environment_compartments <- c("V")
ggplot(data = filter(model_long, compartment %in% environment_compartments), aes(x = time, y = value, color = compartment)) +
  geom_line(size = 1.5) +
  labs(title = "Environment Dynamics - low beta, low q", x = "time (days)", y = "SIV in air (TCID50/m3)") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 20), axis.title = element_text(size = 15))

#Humans
human_compartments <- c("Sh", "Ih", "Rh")
ggplot(data = filter(model_long, compartment %in% human_compartments), aes(x = time, y = value, color = compartment)) +
  geom_line(size = 1.5) +
  labs(title = "Human Dynamics", x = "time (days)", y = "# humans") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 20), axis.title = element_text(size = 15))

filtered_data <- filter(model_long, compartment %in% human_compartments)
print(kable(filtered_data, caption = "Values for Sh, Ih, and Rh compartments"))
