library(tidyverse) #  For data manupulation
library(magrittr)  # For the pipe opeerator.
library(HiddenMarkov)  # For working with HMMs.
library(depmixS4)  # For working with HMMs. 
library(xtable)  # For creating nicely formatted tables, typically for Latex or HTML, from R data structures.
library(forecast)  # For time series forecasting
library(gridExtra)  # To arrange multiple grid-based plots on a page.
library(lubridate)  # To do arithmetic with date-times.
library(HMMpa)  # Analysing Accelerometer Data Using HMMs.


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
# Read in data
dat <- read_csv("BroadShop.csv")

# Sort rows by systemcodenumber, and then by lastupdated
dat %<>% arrange(systemcodenumber, lastupdated)

# overview of data
glimpse(dat)

# freq table to show count of each unique value in a variable
table(dat$systemcodenumber)
table(dat$capacity)

# contingency (two-way freq) table for counts of each combination of variables
table(dat$systemcodenumber, dat$capacity)








#----------- Exploratory data analysis
#===============================================================================
broadstreet <- dat %>% filter(systemcodenumber == "Broad Street")
shopping <- dat %>% filter(systemcodenumber == "Shopping")


# AUTOCORRELATION FUNCTIONS
# The ggAcf function for plotting the autocorrelation function (ACF) in a ggplot2-style plot. 
(acf_bs <- ggAcf(broadstreet$percocc, ci = NULL, color = "gray20") + theme_minimal_adjstd() +
    labs(x = "Lag", y = "ACF", title = NULL))

(acf_ss <- ggAcf(shopping$percocc, ci = NULL, color = "gray20") + theme_minimal_adjstd() +
    labs(x = "Lag", y = "ACF", title = NULL))

# DISTRIBUTION PLOT (Histograms of observed data)

(hist_obs_bs <- broadstreet %>% ggplot(aes(x = percocc)) +
        geom_histogram(color = "white", fill = "gray60", show.legend = F, binwidth = 0.05) +
        scale_fill_brewer(palette = "Set1") + 
        theme_minimal_adjstd() +
        labs(x = "Relative occupancy", y = NULL, title = NULL))

(hist_obs_ss <- shopping %>% ggplot(aes(x = percocc)) +
        geom_histogram(color = "white", fill = "gray60", show.legend = F, binwidth = 0.05) +
        scale_fill_brewer(palette = "Set1") + 
        theme_minimal_adjstd() +
        labs(x = "Relative occupancy", y = NULL, title = NULL))

png("Figures_ak_bs/acf_bs.png", units="in", width = 9, height = 3.5, res = 300)
acf_bs
dev.off()

png("Figures_ak_ss/acf_ss.png", units="in", width = 9, height = 3.5, res = 300)
acf_ss
dev.off()

png("Figures_ak_bs/hist_obs_bs.png", units="in", width = 9, height = 3.5, res = 300)
hist_obs_bs
dev.off()

png("Figures_ak_ss/hist_obs_ss.png", units="in", width = 9, height = 3.5, res = 300)
hist_obs_ss
dev.off()
#########################
# TIME SERIES PLOTS
(ts_bs <- broadstreet %>% ggplot(aes(x = lastupdated, y = percocc)) +
    geom_line(show.legend = F, color = "gray20") +
    scale_color_brewer(palette = "Set1") +
    #facet_wrap(~systemcodenumber, ncol = 2, scales = "free") +
    theme_minimal_adjstd() +
    labs(y = "Relative occupancy", x = "Instance", y = NULL, title = NULL))


(ts_ss <- shopping %>% ggplot(aes(x = lastupdated, y = percocc)) +
        geom_line(show.legend = F, color = "gray20") +
        scale_color_brewer(palette = "Set1") +
        # facet_wrap(~systemcodenumber, ncol = 2, scales = "free") +
        theme_minimal_adjstd() +
        labs(y = "Relative occupancy", x = "Instance", y = NULL, title = NULL))


png("Figures_ak/ts_bs.png", units="in", width = 9, height = 3.5, res = 300)
ts_bs
dev.off()

##############################################################
# Time series zoom in plots
# Define start date
start_date <- as.Date("2016-10-04")
end_date <- start_date + days(15) - 1  # Subtracting 1 to exclude the 16th day


# Filter the data
zoomin_ts1 <- broadstreet %>%
    filter(as.Date(lastupdated) >= start_date & 
               as.Date(lastupdated) <= end_date)

zoomin_ts2 <- shopping %>%
    filter(as.Date(lastupdated) >= start_date & 
               as.Date(lastupdated) <= end_date)
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Determine timezone of your data
tz <- attr(broadstreet$lastupdated, "tzone") %||% "UTC"

# Generate custom breaks for each day using the determined timezone
morning_breaks <- seq(from = as.POSIXct(paste(start_date, "07:45:00"), tz = tz), 
                      to = as.POSIXct(paste(end_date, "07:45:00"), tz = tz), 
                      by = "2 day")

afternoon_breaks <- seq(from = as.POSIXct(paste(start_date, "16:45:00"), tz = tz), 
                        to = as.POSIXct(paste(end_date, "16:45:00"), tz = tz), 
                        by = "2 day")

# Combine the sequences
custom_breaks <- sort(c(morning_breaks, afternoon_breaks))
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



# Vectorized custom label function
custom_labels <- Vectorize(function(datetime) {
    if (is.na(datetime)) {
        return(NA)
    }
    
    day_label <- format(datetime, "%a")
    short_date <- format(datetime, "(%m/%d)")
    
    if (format(datetime, "%H:%M") == "07:45") {
        return(paste0("07:45 16:45\n", day_label, " ", short_date))
    } else {
        return("")
    }
}, vectorize.args = "datetime")

# Plot the data: Broadstreet
(zoomin_ts_bs <- zoomin_ts1 %>%
    ggplot(aes(x = lastupdated, y = percocc)) +
    geom_line(show.legend = FALSE, color = "gray10") +
    theme_minimal_adjstd() +
    theme(axis.text.x = element_text(size = 4)) +
    theme(axis.text.x = element_text(hjust = 0.5)) +
    labs(title = NULL) +
    scale_x_datetime(
        breaks = custom_breaks, 
        labels = custom_labels, 
        expand = c(0,0)
    ) +
    labs(x = "Instance", y = "Relative occupancy"))


# Plot the data: Shopping
(zoomin_ts_ss <- zoomin_ts2 %>%
    ggplot(aes(x = lastupdated, y = percocc)) +
    geom_line(show.legend = FALSE, color = "gray10") +
    theme_minimal_adjstd() +
    theme(axis.text.x = element_text(size = 4)) +
    theme(axis.text.x = element_text(hjust = 0.5)) +
    labs(title = NULL) +
    scale_x_datetime(
        breaks = custom_breaks, 
        labels = custom_labels, 
        expand = c(0,0)
    ) +
    labs(x = "Instance", y = "Relative occupancy"))

acf_bs
acf_ss
hist_obs_bs
hist_obs_ss
ts_bs
ts_ss
zoomin_ts_bs
zoomin_ts_ss
####################################################################################################
# The 8 graphs are arranged in two rows (each for each car park i.e., Broad sreet in 
# the top row and Shopping n the bottom row) and 4 columns (i.e., distributions/histograms 
# in the first column, auto-correlation functions in the second column, full time 
# series plots (for all-days’ times) in the third coolumn and zoomed-in time series 
# plots (for the first-15-days’ times) in the fourth/last column)
png("Figures_ak_bs/eda.png", units="in", width = 9, height = 5, res = 300)
gridExtra::grid.arrange(hist_obs_bs, acf_bs, ts_bs, zoomin_ts_bs, hist_obs_ss, acf_ss, ts_ss, zoomin_ts_ss, ncol = 4)
dev.off()

####################################################################################################
####################################################################################################
####################################################################################################
#                        END FOR NOW & End of EDA!
####################################################################################################
####################################################################################################
####################################################################################################
















########### ANALYSIS #########
#===============================================================================
#---------- Pi's (initial transition probability matrices): Pi with off-diagonal = 0.01
generate_Tmat <- function(n) {
    # Create an n x n matrix filled with 0.01
    mat <- matrix(0.01, nrow = n, ncol = n)
    
    # Adjust the diagonal elements so that each row sums to 1
    diag(mat) <- 1 - 0.01 * (n - 1)
    
    return(mat)
}

# Test the function
print(generate_Tmat(2))
print(generate_Tmat(5))

#===============================================================================
# initial state distribution, i.e., delta
generate_delta <- function(n) {
    # If n is 1, simply return 1
    if(n == 1) return(1)
    
    # Initial raw values
    values <- rep(1/n, n)
    
    # Rounding all but the last value
    rounded <- round(values, 3)
    rounded[n] <- 1 - sum(rounded[-n])
    
    return(rounded)
}

generate_delta(2)  # c(0.5, 0.5)
generate_delta(3)  # close to c(0.333, 0.333, 0.334) but sum is exactly 1
generate_delta(6)  # c(0.2, 0.2, 0.2, 0.2, 0.2)


#===============================================================================
# pm's: observation distribution parameters
generate_prob_list <- function(n) {
    # Set a seed for reproducibility
    set.seed(1720)
    
    # Generate a set of random probabilities from a uniform distribution
    probs <- runif(n)
    
    # Sort the probabilities to ensure they're increasing
    # sorted_probs <- sort(probs)
    
    return(list(prob = probs))
}

# Test the function
print(generate_prob_list(2))
print(generate_prob_list(15))


#===============================================================================
#----------- Broad Street 
# 2-15 states:
# the number of trials for the bin dist i.e., capacity of the car park
pn <- list(size=broadstreet$capacity)
?dthmm
# Then, we can set up the HMMs:
hmmbs2 <- dthmm(broadstreet$occupancy, generate_Tmat(2), generate_delta(2), "binom", generate_prob_list(2), pn, discrete=TRUE)
hmmbs3 <- dthmm(broadstreet$occupancy, generate_Tmat(3), generate_delta(3), "binom", generate_prob_list(3), pn, discrete=TRUE)
hmmbs4 <- dthmm(broadstreet$occupancy, generate_Tmat(4), generate_delta(4), "binom", generate_prob_list(4), pn, discrete=TRUE)
hmmbs5 <- dthmm(broadstreet$occupancy, generate_Tmat(5), generate_delta(5), "binom", generate_prob_list(5), pn, discrete=TRUE)
hmmbs6 <- dthmm(broadstreet$occupancy, generate_Tmat(6), generate_delta(6), "binom", generate_prob_list(6), pn, discrete=TRUE)
hmmbs7 <- dthmm(broadstreet$occupancy, generate_Tmat(7), generate_delta(7), "binom", generate_prob_list(7), pn, discrete=TRUE)
hmmbs8 <- dthmm(broadstreet$occupancy, generate_Tmat(8), generate_delta(8), "binom", generate_prob_list(8), pn, discrete=TRUE)
hmmbs9 <- dthmm(broadstreet$occupancy, generate_Tmat(9), generate_delta(9), "binom", generate_prob_list(9), pn, discrete=TRUE)
hmmbs10 <- dthmm(broadstreet$occupancy, generate_Tmat(10), generate_delta(10), "binom", generate_prob_list(10), pn, discrete=TRUE)
hmmbs11 <- dthmm(broadstreet$occupancy, generate_Tmat(11), generate_delta(11), "binom", generate_prob_list(11), pn, discrete=TRUE)
hmmbs12 <- dthmm(broadstreet$occupancy, generate_Tmat(12), generate_delta(12), "binom", generate_prob_list(12), pn, discrete=TRUE)
hmmbs13 <- dthmm(broadstreet$occupancy, generate_Tmat(13), generate_delta(13), "binom", generate_prob_list(13), pn, discrete=TRUE)
hmmbs14 <- dthmm(broadstreet$occupancy, generate_Tmat(14), generate_delta(14), "binom", generate_prob_list(14), pn, discrete=TRUE)
hmmbs15 <- dthmm(broadstreet$occupancy, generate_Tmat(15), generate_delta(15), "binom", generate_prob_list(15), pn, discrete=TRUE)


hmmbs16 <- dthmm(broadstreet$occupancy, generate_Tmat(16), generate_delta(16), "binom", generate_prob_list(16), pn, discrete=TRUE)
hmmbs17 <- dthmm(broadstreet$occupancy, generate_Tmat(17), generate_delta(17), "binom", generate_prob_list(17), pn, discrete=TRUE)
hmmbs18 <- dthmm(broadstreet$occupancy, generate_Tmat(18), generate_delta(18), "binom", generate_prob_list(18), pn, discrete=TRUE)
hmmbs19 <- dthmm(broadstreet$occupancy, generate_Tmat(19), generate_delta(19), "binom", generate_prob_list(19), pn, discrete=TRUE)
hmmbs20 <- dthmm(broadstreet$occupancy, generate_Tmat(20), generate_delta(20), "binom", generate_prob_list(20), pn, discrete=TRUE)
hmmbs21 <- dthmm(broadstreet$occupancy, generate_Tmat(21), generate_delta(21), "binom", generate_prob_list(21), pn, discrete=TRUE)
hmmbs22 <- dthmm(broadstreet$occupancy, generate_Tmat(22), generate_delta(22), "binom", generate_prob_list(22), pn, discrete=TRUE)
hmmbs23 <- dthmm(broadstreet$occupancy, generate_Tmat(23), generate_delta(23), "binom", generate_prob_list(23), pn, discrete=TRUE)
hmmbs24 <- dthmm(broadstreet$occupancy, generate_Tmat(24), generate_delta(24), "binom", generate_prob_list(24), pn, discrete=TRUE)
hmmbs25 <- dthmm(broadstreet$occupancy, generate_Tmat(25), generate_delta(25), "binom", generate_prob_list(25), pn, discrete=TRUE)
summary(hmmbs15_fit)
hist(residuals(hmmbs15_fit))
hmmbs15_u <- hmmbs15_fit$u


hmmbs2_fit <- BaumWelch(hmmbs2)
hmmbs3_fit <- BaumWelch(hmmbs3)
hmmbs4_fit <- BaumWelch(hmmbs4)
hmmbs5_fit <- BaumWelch(hmmbs5)
hmmbs6_fit <- BaumWelch(hmmbs6)
hmmbs7_fit <- BaumWelch(hmmbs7)
hmmbs8_fit <- BaumWelch(hmmbs8)
hmmbs9_fit <- BaumWelch(hmmbs9)
hmmbs10_fit <- BaumWelch(hmmbs10)
hmmbs11_fit <- BaumWelch(hmmbs11)
hmmbs12_fit <- BaumWelch(hmmbs12)
hmmbs13_fit <- BaumWelch(hmmbs13)
hmmbs14_fit <- BaumWelch(hmmbs14)
hmmbs15_fit <- BaumWelch(hmmbs15)

hmmbs16_fit <- BaumWelch(hmmbs16)
hmmbs17_fit <- BaumWelch(hmmbs17)
hmmbs18_fit <- BaumWelch(hmmbs18)
hmmbs19_fit <- BaumWelch(hmmbs19)
hmmbs20_fit <- BaumWelch(hmmbs20)
hmmbs21_fit <- BaumWelch(hmmbs21)
hmmbs22_fit <- BaumWelch(hmmbs22)
hmmbs23_fit <- BaumWelch(hmmbs23)
hmmbs24_fit <- BaumWelch(hmmbs24)
hmmbs25_fit <- BaumWelch(hmmbs25)

# Calculate AIC, BIC
# ?AIC_HMM
AIC_bs2 <- AIC_HMM(logL = hmmbs2_fit$LL, m = 2, k =1)
AIC_bs3 <- AIC_HMM(logL = hmmbs3_fit$LL, m = 3, k =1)
AIC_bs4 <- AIC_HMM(logL = hmmbs4_fit$LL, m = 4, k =1)
AIC_bs5 <- AIC_HMM(logL = hmmbs5_fit$LL, m = 5, k =1)
AIC_bs6 <- AIC_HMM(logL = hmmbs6_fit$LL, m = 6, k =1)
AIC_bs7 <- AIC_HMM(logL = hmmbs7_fit$LL, m = 7, k =1)
AIC_bs8 <- AIC_HMM(logL = hmmbs8_fit$LL, m = 8, k =1)
AIC_bs9 <- AIC_HMM(logL = hmmbs9_fit$LL, m = 9, k =1)
AIC_bs10 <- AIC_HMM(logL = hmmbs10_fit$LL, m = 10, k =1)
AIC_bs11 <- AIC_HMM(logL = hmmbs11_fit$LL, m = 11, k =1)
AIC_bs12 <- AIC_HMM(logL = hmmbs12_fit$LL, m = 12, k =1)
AIC_bs13 <- AIC_HMM(logL = hmmbs13_fit$LL, m = 13, k =1)
AIC_bs14 <- AIC_HMM(logL = hmmbs14_fit$LL, m = 14, k =1)
AIC_bs15 <- AIC_HMM(logL = hmmbs15_fit$LL, m = 15, k =1)
AIC_bs16 <- AIC_HMM(logL = hmmbs16_fit$LL, m = 16, k =1)
AIC_bs17 <- AIC_HMM(logL = hmmbs17_fit$LL, m = 17, k =1)
AIC_bs18 <- AIC_HMM(logL = hmmbs18_fit$LL, m = 18, k =1)

AIC_bs19 <- AIC_HMM(logL = hmmbs19_fit$LL, m = 19, k =1)
AIC_bs20 <- AIC_HMM(logL = hmmbs20_fit$LL, m = 20, k =1)
AIC_bs21 <- AIC_HMM(logL = hmmbs21_fit$LL, m = 21, k =1)
AIC_bs22 <- AIC_HMM(logL = hmmbs22_fit$LL, m = 22, k =1)
AIC_bs23 <- AIC_HMM(logL = hmmbs23_fit$LL, m = 23, k =1)
AIC_bs24 <- AIC_HMM(logL = hmmbs24_fit$LL, m = 24, k =1)
AIC_bs25 <- AIC_HMM(logL = hmmbs25_fit$LL, m = 25, k =1)


BIC_bs2 <- BIC_HMM(size = nrow(broadstreet), m = 2, k = 1, hmmbs2_fit$LL)
BIC_bs3 <- BIC_HMM(size = nrow(broadstreet), m = 3, k = 1, hmmbs3_fit$LL)
BIC_bs4 <- BIC_HMM(size = nrow(broadstreet), m = 4, k = 1, hmmbs4_fit$LL)
BIC_bs5 <- BIC_HMM(size = nrow(broadstreet), m = 5, k = 1, hmmbs5_fit$LL)
BIC_bs6 <- BIC_HMM(size = nrow(broadstreet), m = 6, k = 1, hmmbs6_fit$LL)
BIC_bs7 <- BIC_HMM(size = nrow(broadstreet), m = 7, k = 1, hmmbs7_fit$LL)
BIC_bs8 <- BIC_HMM(size = nrow(broadstreet), m = 8, k = 1, hmmbs8_fit$LL)
BIC_bs9 <- BIC_HMM(size = nrow(broadstreet), m = 9, k = 1, hmmbs9_fit$LL)
BIC_bs10 <- BIC_HMM(size = nrow(broadstreet), m = 10, k = 1, hmmbs10_fit$LL)
BIC_bs11 <- BIC_HMM(size = nrow(broadstreet), m = 11, k = 1, hmmbs11_fit$LL)
BIC_bs12 <- BIC_HMM(size = nrow(broadstreet), m = 12, k = 1, hmmbs12_fit$LL)
BIC_bs13 <- BIC_HMM(size = nrow(broadstreet), m = 13, k = 1, hmmbs13_fit$LL)
BIC_bs14 <- BIC_HMM(size = nrow(broadstreet), m = 14, k = 1, hmmbs14_fit$LL)
BIC_bs15 <- BIC_HMM(size = nrow(broadstreet), m = 15, k = 1, hmmbs15_fit$LL)
BIC_bs16 <- BIC_HMM(size = nrow(broadstreet), m = 16, k = 1, hmmbs16_fit$LL)
BIC_bs17 <- BIC_HMM(size = nrow(broadstreet), m = 17, k = 1, hmmbs17_fit$LL)
BIC_bs18 <- BIC_HMM(size = nrow(broadstreet), m = 18, k = 1, hmmbs18_fit$LL)



BIC_bs19 <- BIC_HMM(size = nrow(broadstreet), m = 19, k = 1, hmmbs19_fit$LL)
BIC_bs20 <- BIC_HMM(size = nrow(broadstreet), m = 20, k = 1, hmmbs20_fit$LL)
BIC_bs21 <- BIC_HMM(size = nrow(broadstreet), m = 21, k = 1, hmmbs21_fit$LL)
BIC_bs22 <- BIC_HMM(size = nrow(broadstreet), m = 22, k = 1, hmmbs22_fit$LL)
BIC_bs23 <- BIC_HMM(size = nrow(broadstreet), m = 23, k = 1, hmmbs23_fit$LL)
BIC_bs24 <- BIC_HMM(size = nrow(broadstreet), m = 24, k = 1, hmmbs24_fit$LL)
BIC_bs25 <- BIC_HMM(size = nrow(broadstreet), m = 25, k = 1, hmmbs25_fit$LL)

AIC_BIC_df_bs <- cbind(state = 2:18, 
                       k = c(5, 11, 19, 29, 41, 55, 71, 89, 109, 111, 222, 333, 444, 555, 666, 777, 888),
                       loglik = c(-hmmbs2_fit$LL, -hmmbs3_fit$LL, -hmmbs4_fit$LL, -hmmbs5_fit$LL,
                                  -hmmbs6_fit$LL, -hmmbs7_fit$LL, -hmmbs8_fit$LL, -hmmbs9_fit$LL, -hmmbs10_fit$LL,
                                  -hmmbs11_fit$LL, -hmmbs12_fit$LL, -hmmbs13_fit$LL, -hmmbs14_fit$LL, -hmmbs15_fit$LL,
                                  -hmmbs16_fit$LL, -hmmbs17_fit$LL, -hmmbs18_fit$LL),
                       AIC = c(AIC_bs2, AIC_bs3, AIC_bs4, AIC_bs5, AIC_bs6, AIC_bs7, AIC_bs8, AIC_bs9, AIC_bs10, AIC_bs11, AIC_bs12, AIC_bs13, AIC_bs14, AIC_bs15, AIC_bs16, AIC_bs17, AIC_bs18),
                       BIC = c(BIC_bs2, BIC_bs3, BIC_bs4, BIC_bs5, BIC_bs6, BIC_bs7, BIC_bs8, BIC_bs9, BIC_bs10, BIC_bs11, BIC_bs12, BIC_bs13, BIC_bs14, BIC_bs15, BIC_bs16, BIC_bs17, BIC_bs18))
xtable(AIC_BIC_df_bs)

# Making AIC_BIC_df_bs a dataframe
AIC_BIC_df_bs_df <- data.frame(AIC_BIC_df_bs)

# Extracting the state corresponding to the lowest AIC
(state_lowest_AIC <- AIC_BIC_df_bs_df$state[which.min(AIC_BIC_df_bs_df$AIC)])
(lowest_AIC <- min(AIC_BIC_df_bs_df$AIC))

# Extracting the state corresponding to the lowest BIC
(state_lowest_BIC <- AIC_BIC_df_bs_df$state[which.min(AIC_BIC_df_bs_df$BIC)])
(lowest_BIC <- min(AIC_BIC_df_bs_df$BIC))

cat("The state with the lowest AIC is:", state_lowest_AIC, "with AIC value:", lowest_AIC, "\n")
cat("The state with the lowest BIC is:", state_lowest_BIC, "with BIC value:", lowest_BIC)



###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
###@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@###
# Plot -State prdiction: Broad street
# # Filter the data
# # Define start date
start_date <- as.Date("2016-10-04")
end_date <- start_date + days(15) - 1  # Subtracting 1 to exclude the 16th day

zoomin_ts11 <- broadstreet %>%
    filter(as.Date(lastupdated) >= start_date & 
               as.Date(lastupdated) <= end_date)


zoomin_ts1 %>%
    ggplot(aes(x = lastupdated, y = percocc)) +
    geom_line(show.legend = FALSE, color = "gray10") +
    theme_minimal_adjstd() +
    labs(title = "Broadstreet (first 15 days)") +
    scale_x_datetime(
        breaks = custom_breaks, 
        labels = custom_labels, 
        expand = c(0,0)
    ) +
    labs(x = "Instance", y = "Relative occupancy") 




summary(hmmbs15_fit)
hmmbs15_u <- hmmbs15_fit$u
hmmbs15_fitu <- Viterbi(hmmbs15_fit)
hmmbs15_res <- residuals(hmmbs15_fit)

(state_pred_bs_15 <- broadstreet %>% mutate(viterbi = Viterbi(hmmbs15_fit)) %>%
        filter(as.Date(lastupdated) >= start_date & 
                   as.Date(lastupdated) <= end_date) %>%
        mutate(viterbi2 = case_when(
            viterbi == 1 ~ 0.5669913,
            viterbi == 2 ~ 0.8624401,
            viterbi == 3 ~ 0.9116224,
            viterbi == 4 ~ 0.0986054,
            viterbi == 5 ~ 0.9818833,
            viterbi == 6 ~ 0.7971906,
            viterbi == 7 ~ 0.6425767,
            viterbi == 8 ~ 0.2055812,
            viterbi == 9 ~ 0.2689094,
            viterbi == 10 ~ 0.8001092,
            viterbi == 11 ~ 0.1438808,
            viterbi == 12 ~ 0.7215477,
            viterbi == 13 ~ 0.3432669,
            viterbi == 14 ~ 0.4284895,
            TRUE         ~ 0.9520501  # This is the default condition, similar to the last 'else' in a nested ifelse
        )) %>%
        ggplot(aes(x = lastupdated)) +
        geom_hline(aes(yintercept = viterbi2), color = "gray80", linetype = 1) +
        geom_line(aes(y = percocc, group = 1), color = "gray70") +
        geom_point(aes(y = viterbi2), size = 5, shape = "-", color = "black") +
        #theme(axis.text.x = element_text(size = 5)) +
        # theme(axis.text.x = element_text(hjust = -0.5)) +
        # scale_color_manual(values = mycolorll) +
        scale_y_continuous(limits = c(0.05,1)) +
        theme_minimal_adjstd() +
        theme(
            axis.text.x = element_text(size = 6, hjust = 0.3),
            # axis.text.x = element_text(hjust = 0.3),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            #panel.grid.major.x = element_line(color = "gray95"),
            #panel.grid.minor.x = element_line(color = "gray95")
        ) +
        scale_x_datetime(
            breaks = custom_breaks, 
            labels = custom_labels, 
            expand = c(0,0)
        ) +
        labs(x = "Time", y = "Relative occupancy", color = "States"))


png("Figures_ak/state_pred_bs_15.png", units="in", width = 9, height = 3.5, res = 300)
state_pred_bs_15
dev.off()
############################################
############################################
############################################
############################################


summary(hmmbs15_fit)
Pi_matrix <- summary(hmmbs15_fit)$Pi
latex_code <- xtable(Pi_matrix)
print(latex_code, type = "latex")

# Adding the initial state distribution (delta) and the emission probabilities (pm$prob) as two additional rows at the bottom of your matrix before exporting to LaTeX, you can do this using the rbind() function.
# Extracting the matrices and vectors needed:
Pi_matrix <- summary(hmmbs15_fit)$Pi
delta <- summary(hmmbs15_fit)$delta
emission_probs <- summary(hmmbs15_fit)$pm$prob

# Binding these two vectors as rows to the Pi_matrix:
final_matrix <- rbind(Pi_matrix, delta = delta, emission_probs = emission_probs)

# Converting this enhanced matrix to LaTeX format:
latex_code <- xtable(final_matrix)
print(latex_code, type = "latex")

hmmbs15_u <- hmmbs15_fit$u
hmmbs15_fitu <- Viterbi(hmmbs15_fit)
hmmbs15_res <- residuals(hmmbs15_fit)

#---------- Diagnosis
#---------- Diagnosis for 15 states
hmmbs15_res_df1 <- data.frame(state = rep(paste0(15, " States"), each = 1312), 
                         resid = c(hmmbs15_res),
                         time = rep(broadstreet$lastupdated, 1)) %>% as_tibble() %>%
    mutate(state = factor(state, levels = c("15 States")))

hmmbs_res15_df <- hmmbs15_res_df1 %>% filter(state == "15 States")


(diag_st15_1 <- hmmbs_res15_df %>% ggplot(aes(x = time, y = resid)) +
    geom_hline(yintercept = 1.96, linetype="solid", color="gray70") +
    geom_hline(yintercept = -1.96, linetype="solid", color="gray70") +
    geom_hline(yintercept = 0, linetype="solid", color="gray70") +
    geom_hline(yintercept = 2.58, linetype="solid", color="gray70") +
    geom_hline(yintercept = -2.58, linetype="solid", color="gray70") +
    geom_point(size = 1.5, alpha = 0.6, color = "gray50") +
    coord_cartesian(ylim = c(-4, 4)) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL, title = NULL))


(diag_st15_2 <- hmmbs_res15_df %>%
    ggplot(aes(x = resid)) +
    geom_histogram(aes(y = ..density..), fill = "gray50", color = "white") +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(hmmbs_res15_df$resid), 
                              sd = sd(hmmbs_res15_df$resid)), linewidth = 1) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL))


(diag_st15_3 <-hmmbs_res15_df %>%
    ggplot(aes(sample = resid)) +
    stat_qq() +
    geom_abline() +
    coord_cartesian(ylim = c(-5, 5), xlim = c(-5, 5)) +
    theme_minimal_adjstd() +
    labs(x = NULL, y = NULL, title = NULL))


png("Figures_ak_bs/hmmbs15_diagns.png", units="in", width = 9, height = 5, res = 300)
gridExtra::grid.arrange(diag_st15_1, diag_st15_2, diag_st15_3, ncol = 1)
dev.off()

png("Figures_ak_bs/hmmbs15_state_pred_bs_ss_row.png", units="in", width = 9, height = 5, res = 300)
gridExtra::grid.arrange(state_pred_bs_15, state_pred_ss_15, ncol = 2)
dev.off()
png("Figures_ak_bs/hmmbs15_state_pred_bs_ss_col.png", units="in", width = 9, height = 5, res = 300)
gridExtra::grid.arrange(state_pred_bs_15, state_pred_ss_15, ncol = 1)
dev.off()

png("Figures_ak/diag_st15_3.png", units="in", width = 9, height = 5, res = 300)
diag_st15_3
dev.off()



#############################################################################
#############################################################################
#############################################################################
AIC_BIC_df2 <- cbind(state = 2:15, 
                     k = c(5,11,19,29,41,55,71,89,109,131,155,181,209,239),
                     loglik = c(-mod_sp2_fit$LL, -mod_sp3_fit$LL, -mod_sp4_fit$LL, -mod_sp5_fit$LL,
                                -mod_sp6_fit$LL, -mod_sp7_fit$LL, -mod_sp8_fit$LL, 
                                -mod_sp9_fit$LL, -mod_sp10_fit$LL, -mod_sp11_fit$LL,
                                -mod_sp12_fit$LL, -mod_sp13_fit$LL, -mod_sp14_fit$LL, -mod_sp15_fit$LL),
                     AIC = c(AIC_sp2, AIC_sp3, AIC_sp4, AIC_sp5, AIC_sp6, AIC_sp7, 
                             AIC_sp8, AIC_sp9, AIC_sp10,AIC_sp11, AIC_sp12, AIC_sp13,
                             AIC_sp14, AIC_sp15),
                     BIC = c(BIC_sp2, BIC_sp3, BIC_sp4, BIC_sp5, BIC_sp6, BIC_sp7, 
                             BIC_sp8, BIC_sp9, BIC_sp10, BIC_sp11, BIC_sp12, BIC_sp13,
                             BIC_sp14, BIC_sp15))
