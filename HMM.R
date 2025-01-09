library(devtools)
#install_version("HMMCont", version = "1.0", repos = "https://cloud.r-project.org")
install_version("HMMpa", version = "1.0.1", repos = "https://cloud.r-project.org")

library(HiddenMarkov)
library(ggplot2)
library(dplyr)
library(bayesforecast)
library(forecast)
library(HMMpa)
library(gridExtra)


#----------  Set up theme
#===============================================================================
theme_minimal_adjstd <- function(...) {
  theme_minimal() +
    theme(
      plot.title = element_text(color = "gray0",  face = "bold", hjust = 0.5, size = 8),
      axis.line = element_line(linetype = "solid"),
      axis.text.x = element_text(color = "gray0"),
      axis.text.y = element_text(color = "gray0"),
      axis.title.x = element_text(color = "gray0", size = 6),
      axis.title.y = element_text(color = "gray0", size = 6),
      axis.ticks.y = element_line(),
      axis.ticks.x = element_line(),            panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold"),
      legend.direction = "horizontal",
      legend.position = "topright",
      legend.background = element_rect(fill = NA, color = NA),
      # legend.text = element_text(size = 14),
      legend.key.width = unit(2, "line"),
      strip.text = element_text(face = "bold"),
      strip.background = element_rect(fill = NA, color = NA)
    )
}

?element_text

#===============================================================================

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

gen_Tmatrix <- function(n) {
  # Create an n x n matrix filled with 0.01
  mat <- matrix(0.01, nrow = n, ncol = n)
  # Adjust the diagonal elements so that each row sums to 1
  diag(mat) <- 1 - 0.01 * (n - 1)
  
  return(mat)
}

# Test the function
print(gen_Tmatrix(3))

func_delta <- function(n) {
  if(n==1) return (1)
  values <- rep(1/n, n)
  rounded <- round(values, 3)
  rounded[n] <- 1 - sum(rounded[-n])

  return (rounded)
}
func_delta(6)

# func_pm <- function(n) {
#   set.seed(3745)
#   mean_power <- 
#   probs <- runif(n)
#   return(list(prob = probs))}

emission_means <- seq(min(mean_power), max(mean_power), length.out = num_states)
emission_sd <- rep(sd(mean_power), num_states)

hmmbs2 <- dthmm(hmmdata$mean_power, Pi=gen_Tmatrix(15), delta=func_delta(15),
                "norm", pm = list(mean = emission_means, sd = emission_sd))

summary(hmmbs2)
try <- BaumWelch(hmmbs2)

##############################################################################
#### ----- alternative way -----------------------------

mean_power <- hmmdata$mean_power

##### Define initial model parameters -------

gen_Tmatrix <- function(n) {
  mat <- matrix(0.01, nrow = n, ncol = n)
  diag(mat) <- 1 - 0.01 * (n - 1)
  return(mat)
}
func_delta <- function(n) {
  if(n==1) return (1)
  values <- rep(1/n, n)
  rounded <- round(values, 3)
  rounded[n] <- 1 - sum(rounded[-n])
  
  return (rounded)
}

gen_emisparams <- function(data, n) {
  set.seed(374)
  # Generate random means and standard deviations
  emission_means <- runif(n, min = min(data), max = max(data))
  emission_sd <- runif(n, min = 0.7 * sd(mean_power), max = sd(data))
  return(list(mean = emission_means, sd = emission_sd))
}

gen_emisparamsq <- function(data, n) {
  # Calculate quantiles and determine means as midpoints of quantile ranges
  quantiles <- quantile(data, probs = seq(0, 1, length.out = n + 1))
  emission_means <- (quantiles[-1] + quantiles[-length(quantiles)]) / 2
  
  # Calculate standard deviations within each quantile range
  emission_sd <- sapply(1:n, function(i) {
    range_data <- data[data >= quantiles[i] & data <= quantiles[i + 1]]
    if (length(range_data) > 1) {
      sd(range_data)  # Standard deviation within the range
    } else {
      0.1 * sd(data)  # Small default value if insufficient data
    }
  })
  
  # Return as a named list
  return(list(mean = emission_means, sd = emission_sd))
}


# Define the HMM
n = 15
hmm_model1 <- dthmm(
  x = mean_power,
  Pi = gen_Tmatrix(n),
  delta = func_delta(n),
  pm = gen_emisparamsq(hmmdata$mean_power, n), 
  distn = "norm"
)

fitted_hmm <- BaumWelch(hmm_model1, control = bwcontrol(maxiter = 1000))
AIC_HMM(logL = fitted_hmm$LL, m=n, k=2)
BIC_HMM(size = nrow(hmmdata), m=n, k=2, fitted_hmm$LL)

fit_hmm11 <- BaumWelch(hmm_model1, control = bwcontrol(maxiter = 1000))
fit_hmm15 <- BaumWelch(hmm_model1, control = bwcontrol(maxiter = 1000))

#summary(fitted_hmm)


Viterbi(hmm_BW_n) # State decoding
fit15res <- residuals(fit_hmm15) # Residuals
fit11res <- residuals(fit_hmm11)

#### plot residuals - state 15 -------------------------------------------------------
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
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL, title = NULL))

(diag_st15_2 <- fit15res_df1 %>%
    ggplot(aes(x = resid)) +
    geom_histogram(aes(y = ..density..), fill = "gray50", color = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(fit15res_df1$resid), 
                              sd = sd(fit15res_df1$resid)), linewidth = 1) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL))

(diag_st15_3 <-fit15res_df1 %>%
    ggplot(aes(sample = resid)) +
    stat_qq() +
    geom_abline() +
    coord_cartesian(ylim = c(-4, 4), xlim = c(-4, 4)) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL, title = NULL))


#######--------------------------------

#### plot residuals  state 11 --------------------------------------------------------
fit11res_df <- data.frame(state = rep(paste0(15, " States"), each = 1416), 
                           resid = c(fit11res),
                           time = rep(hmmdata$day_in_year, 1)) %>% as_tibble() %>%
  mutate(state = factor(state, levels = c("11 States")))

(diag_st11_1 <- fit11res_df %>% ggplot(aes(x = time, y = resid)) +
    geom_hline(yintercept = 1.96, linetype="solid", color="gray70") +
    geom_hline(yintercept = -1.96, linetype="solid", color="gray70") +
    geom_hline(yintercept = 0, linetype="solid", color="gray70") +
    geom_hline(yintercept = 2.58, linetype="solid", color="gray70") +
    geom_hline(yintercept = -2.58, linetype="solid", color="gray70") +
    geom_point(size = 1.5, alpha = 0.6, color = "gray50") +
    coord_cartesian(ylim = c(-4, 4)) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL, title = NULL))

(diag_st11_2 <- fit11res_df %>%
    ggplot(aes(x = resid)) +
    geom_histogram(aes(y = ..density..), fill = "gray50", color = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(fit15res_df1$resid), 
                              sd = sd(fit15res_df1$resid)), linewidth = 1) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL))

(diag_st11_3 <-fit11res_df %>%
    ggplot(aes(sample = resid)) +
    stat_qq() +
    geom_abline() +
    coord_cartesian(ylim = c(-4, 4), xlim = c(-4, 4)) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL, title = NULL))


png("pictures/hmmb15_diagns.png", units="cm", width = 20, height = 10, res = 300)
gridExtra::grid.arrange(diag_st15_1, diag_st15_2, diag_st15_3, ncol = 3)
dev.off()

png("pictures/hmmb11_diagns.png", units="cm", width = 20, height = 10, res = 300)
gridExtra::grid.arrange(diag_st11_1, diag_st11_2, diag_st11_3, ncol = 3)
dev.off()
