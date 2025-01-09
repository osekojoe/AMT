# Load required libraries
library(HiddenMarkov)

# Read the dataset
data <- read.csv("z1_hr_mean.csv")

# Inspect the dataset structure
str(data)

# Add a column for 'Day/Night'
data$time_period <- ifelse(data$hr >= 6 & data$hr <= 18, "Day", "Night")

# Split the dataset into day and night data
day_data <- data[data$time_period == "Day", ]
night_data <- data[data$time_period == "Night", ]


# Extract mean power consumption for modeling
day_consumption <- day_data$mean_power
night_consumption <- night_data$mean_power

# Define a function to fit HMM and evaluate
fit_hmm <- function(data, states) {
  # Initial probabilities, transition matrix, and emission parameters
  init_probs <- rep(1 / states, states)
  trans_probs <- matrix(1 / states, nrow = states, ncol = states)
  means <- seq(min(data), max(data), length.out = states)
  sds <- rep(sd(data), states)
  
  # Define HMM
  hmm_model <- dthmm(x = data, Pi = trans_probs, delta = init_probs,
    distn = "norm", pm = list(mean = means, sd = sds)
  )
  
  # Fit HMM using Baum-Welch algorithm
  fitted_model <- BaumWelch(hmm_model)
  
  return(fitted_model)
}

n=2
hmm_day_n <- fit_hmm(day_consumption, states=n)
AIC_HMM(logL = hmm_day_n$LL, m=n, k=2)
BIC_HMM(size = nrow(day_data), m=n, k=2, hmm_day_n$LL)

n=2
hmm_night_n <- fit_hmm(night_consumption, states=n)
AIC_HMM(logL = hmm_night_n$LL, m=n, k=2)
BIC_HMM(size = nrow(night_data), m=n, k=2, hmm_night_n$LL)


## ------------------------------------------------------------------------
# Visualize the decoded states
day_states <- Viterbi(hmm_day_n)
night_states <- Viterbi(hmm_night_n)

# Add decoded states to the dataset
day_data$state <- day_states
night_data$state <- night_states

# Plot states

ggplot(day_data, aes(x = hr, y = mean_power, color = as.factor(state))) +
  geom_line() +
  labs(title = "Daytime Power Consumption States", x = "Hour",
       y = "Mean Power (KWh)", color = "State")

ggplot(night_data, aes(x = hr, y = mean_power, color = as.factor(state))) +
  geom_line() +
  labs(title = "Nighttime Power Consumption States", x = "Hour", 
       y = "Mean Power (KWh)", color = "State")
