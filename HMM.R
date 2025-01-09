library(devtools)
#install_version("HMMCont", version = "1.0", repos = "https://cloud.r-project.org")
install_version("HMMpa", version = "1.0.1", repos = "https://cloud.r-project.org")

library(HiddenMarkov)
library(ggplot2)
library(dplyr)
library(bayesforecast)
library(forecast)
library(HMMpa)


hmmdata <- read.csv('z1_hr_mean.csv')
hmmdata <- hmmdata[order(hmmdata$seq),]

# Step 2: Data cleaning
hmmdata$day <- as.Date(hmmdata$day, format = "%d/%m/%Y") # Convert to Date format
hmmdata <- hmmdata %>% arrange(day, hr) # Ensure data is sorted by date and hour

mean_power <- hmmdata$mean_power

## exploratory data analysis ----------------------------------------------

par(mfrow=c(2,2))
# Step 3: Histogram of mean power
#png(filename = "pictures/histogram.png", width = 400, height = 350)
png("pictures/histogram.png", width=18, height=10, units="cm", res=300)
ggplot(hmmdata, aes(x = mean_power)) +
  geom_histogram(binwidth = 1000, position = "dodge", color = "white") +
  scale_color_brewer(palette = "Set1") +
  labs(title = "", x = "Mean Power (KWh)", y = "Frequency")
dev.off()

# step 4: autocorrelation plot
ts_data <- ts(hmmdata$mean_power, frequency = 24) # Assuming hourly data with daily seasonality
ggAcf(ts_data) +
  labs(title = "")

#png(filename = "pictures/acf.png", width = 400, height = 350)
png("pictures/acf.png", width=18, height=10, units="cm", res=300)
ggAcf(mean_power) +
  labs(title = "")
dev.off()

# Step 5: Time-series plot
#png(filename = "pictures/comptimeseries.png", width = 400, height = 400)
png("pictures/comptimeseries.png", width=18, height=10, units="cm", res=300)
ggplot(hmmdata, aes(x = as.POSIXct(paste(day, hr, sep = " "),
                                   format = "%Y-%m-%d %H"), y = mean_power)) +
  geom_line(color = "blue") +
  labs(title = "", x = "Time", y = "Mean Power (KWh)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_datetime(date_labels = "%b", date_breaks = "1 month")
  #scale_x_datetime(date_labels = "%d-%b", date_breaks = "1 day")
dev.off()

# Time-series plot zoomed in for the first 14 days
#png(filename = "pictures/zoomtimeseries.png", width = 400, height = 400)
png("pictures/zoomtimeseries.png", width=18, height=10, units="cm", res=300)
ggplot(
  hmmdata %>%
    filter(day >= as.Date("2017-01-01") & day <= as.Date("2017-01-14")), 
  aes(x = as.POSIXct(paste(day, hr, sep = " "), 
                     format = "%Y-%m-%d %H"), y = mean_power)) +
  geom_line(color = "blue") +
  labs(title = "", x = "Time", y = "Mean Power (KWh)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_datetime(date_labels = "%d-%b", date_breaks = "1 day")
dev.off()

## - explore time-based patterns
png("pictures/hourseries.png", width=18, height=10, units="cm", res=300)
ggplot(hmmdata, aes(x = as.factor(hr), y = mean_power)) +
  geom_boxplot() +
  labs(
    title = "Mean Power Consumption by Hour",
    x = "Hour of Day",
    y = "Mean Power Consumption (KWh)"
  ) +
  theme_minimal()
dev.off()


## analysis  ----------------------------------------------

# func_transitMat <- function() {
#   mat <- matrix(0.01, nrow = n, ncol = n)
#   diag(mat())
# }
# 
# func_delta <- function() {
#   if(n==1) return (1)
#   values <- rep(1/n, n)
#   rounded <- round(values, 3)
#   rounded[n] <- 1 - sum(rounded[-n])
#   
#   return (rounded)
# }
# 
# 
# func_pm <- function(n) {
#   set.seed(3745)
#   probs <- runif(n)
#   return(list(prob = probs))}


#### ----- alternative way -----------------------------

mean_power <- hmmdata$mean_power

# Define initial model parameters
num_states <- 2  # Change the number of states to experiment
init_prob <- rep(1 / num_states, num_states)  # Equal probability for initial states
transition_matrix <- matrix(1 / num_states, nrow = num_states, ncol = num_states)  # Equal probabilities for transitions
emission_means <- seq(min(mean_power), max(mean_power), length.out = num_states)
emission_sd <- rep(sd(mean_power), num_states)

# Define the HMM
hmm_model1 <- dthmm(
  x = mean_power,
  Pi = transition_matrix,
  delta = init_prob,
  pm = list(mean = emission_means, sd = emission_sd),
  distn = "norm"
)


# Fit the HMM using Baum-Welch algorithm
fitted_hmm <- BaumWelch(hmm_model1)


# Define a function to fit HMM and evaluate
fit_hmm1 <- function(data, num_states) {
  # Initial probabilities, transition matrix, and emission parameters
  init_prob <- rep(1 / num_states, num_states)  # Equal probability for initial states
  transition_matrix <- matrix(1 / num_states, nrow = num_states, ncol = num_states)  # Equal probabilities for transitions
  emission_means <- seq(min(mean_power), max(mean_power), length.out = num_states)
  emission_sd <- rep(sd(mean_power), num_states)
  
  # Define the HMM
  hmm_model <- dthmm(x = mean_power, Pi = transition_matrix,
    delta = init_prob, pm = list(mean = emission_means, sd = emission_sd),
    distn = "norm"
  )
  
  # Fit HMM using Baum-Welch algorithm
  fitted_model <- BaumWelch(hmm_model)
  
  return(fitted_model)
}

m=15
hmm_BW_n <- fit_hmm1(mean_power, num_states=m)

AIC_HMM(logL = hmm_BW_n$LL, m=m, k=2)

BIC_HMM(size = nrow(hmmdata), m=20, k=2, hmm_BW_n$LL)
Viterbi(hmm_BW_n) # State decoding
fit15res <- residuals(hmm_BW_n) # Residuals

#### plot residuals --------------------------------------------------------
fit15res_df1 <- data.frame(state = rep(paste0(15, " States"), each = 1416), 
                              resid = c(fit15res),
                              time = rep(hmmdata$day_in_year, 1)) %>% as_tibble() %>%
  mutate(state = factor(state, levels = c("15 States")))

(diag_st15_1 <- fit15res_df1 %>% ggplot(aes(x = time, y = resid)) +
    geom_hline(yintercept = 1.96, linetype="solid", color="gray70") +
    geom_hline(yintercept = -1.96, linetype="solid", color="gray70") +
    geom_hline(yintercept = 0, linetype="solid", color="gray70") +
    geom_hline(yintercept = 2.58, linetype="solid", color="gray70") +
    geom_hline(yintercept = -2.58, linetype="solid", color="gray70") +
    geom_point(size = 1.5, alpha = 0.6, color = "gray50") +
    coord_cartesian(ylim = c(-4, 4)) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = NULL))

(diag_st15_2 <- fit15res_df1 %>%
    ggplot(aes(x = resid)) +
    geom_histogram(aes(y = ..density..), fill = "gray50", color = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(fit15res_df1$resid), 
                              sd = sd(fit15res_df1$resid)), linewidth = 1) +
    theme_minimal() +
    labs(x = NULL, y = NULL))

(diag_st15_3 <-fit15res_df1 %>%
    ggplot(aes(sample = resid)) +
    stat_qq() +
    geom_abline() +
    coord_cartesian(ylim = c(-4, 4), xlim = c(-4, 4)) +
    theme_minimal() +
    labs(x = NULL, y = NULL, title = NULL))


#######---------consider hour column  --------------------------------------

hmmdata$hour_sin <- sin(2 * pi * hmmdata$hr / 24)
hmmdata$hour_cos <- cos(2 * pi * hmmdata$hr / 24)

observations <- cbind(mean_power, hmmdata$hour_sin, hmmdata$hour_cos)

hmm_model2 <- dthmm(
  x = observations,
  Pi = transition_matrix,
  delta = init_prob,
  pm = list(mean = list(c(mean(mean_power), 0, 0)), cov = diag(c(var(mean_power), 0.1, 0.1))),
  distn = "mvnorm",
  discrete = TRUE
)

fitted_hmm2 <- BaumWelch(hmm_model)
AIC_HMM(logL = fitted_hmm$LL, m=2, k=1)
BIC_HMM(size = nrow(hmmdata), m = 2, k = 1, fitted_hmm$LL)

