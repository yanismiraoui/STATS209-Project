### SETTING UP THE ENVIR ONMENT ###
setwd("/Users/yanis/Documents/Stanford/STATS 209/Project/")
getwd()

library(estimatr)
library(tidyverse)
library(magrittr)
library(ggplot2)
library(DOS2)
library(optmatch)
library(RItools)
library(plyr)
library(rcbalance)
library(ggplot2)
library(grf)

source('utility.R') # for summarize_match function

### PREPROCESSING ###

# Load data csv file
data <- read.csv("./data/ABIDE_tab_compressed.csv", header = TRUE)
nrow(data)
colnames(data)
summary(data)
data$Z <- as.numeric(data$DX_GROUP) - 1
data$Z <- factor(data$Z)

# Plot the distribution of subjects diagnosed with ASD (using red) and controls (using blue)
ggplot(data, aes(x=Z, fill=Z)) + geom_bar() + scale_fill_manual(values=c("blue", "red")) + labs(x="Diagnosis", y="Number of subjects", title="Distribution of subjects diagnosed with ASD and controls", fill="Diagnosis", caption="Source: ABIDE")

# Number of students in each arm
table(data$Z)

# count number of na values in FIQ, VIQ and PIQ
sum(is.na(data$FIQ))
sum(is.na(data$VIQ))
sum(is.na(data$PIQ))

# get rid of rows with na values in FIQ, VIQ and PIQ
cat("Number of rows before cleaning: ", nrow(data), "\n")
data <- data[!is.na(data$FIQ),]
data <- data[!is.na(data$VIQ),]
data <- data[!is.na(data$PIQ),]
data <- data[data$FIQ != -9999,]
data <- data[data$VIQ != -9999,]
data <- data[data$PIQ != -9999,]
data <- data[data$AGE_AT_SCAN != -9999,]
data <- data[data$HANDEDNESS_CATEGORY != -9999,]
data <- data[data$CURRENT_MED_STATUS != -9999,]
cols<- c("Z", "FIQ", "VIQ", "PIQ", "SEX", "AGE_AT_SCAN", "HANDEDNESS_CATEGORY", "CURRENT_MED_STATUS")
# delete all the rows with na values in the columns of interest
data <- data[complete.cases(data[,cols]),]
cat("Number of rows after cleaning: ", nrow(data), "\n")

# Keep only the columns where we have MRI data
data <- data[!is.na(data$compressed_2_1),]
cat("Number of rows: ", nrow(data), "\n")

# get rid of rows with '`' in column CURRENT_MED_STATUS
data <- data[!grepl("`", data$CURRENT_MED_STATUS),]
cat("Number of rows: ", nrow(data), "\n")

# Covariates: SEX, AGE_AT_SCAN, HANDEDNESS, CURRENT_MED_STATUS, MRI_FEATURES
SIZE_COMPRESSED_MRI <- 3
# Create a vector of the columns with names f"compressed_{SIZE_COMPRESSED_MRI}_{i}" for i in 1:SIZE_COMPRESSED_MRI
cols_mri <- c()
for (i in 1:SIZE_COMPRESSED_MRI) {
  cols_mri <- c(cols_mri, paste0("compressed_", SIZE_COMPRESSED_MRI, "_", i))
}
cols_mri

# Covariates: SEX, AGE_AT_SCAN, HANDEDNESS, CURRENT_MED_STATUS, MRI_FEATURES
cols_covariates <- c("Z", "FIQ", "VIQ", "PIQ", "SEX", "AGE_AT_SCAN", "HANDEDNESS_CATEGORY", "CURRENT_MED_STATUS", cols_mri, "prop")
cols_covariates

### PROPENSITY SCORE MATCHING ###

# Analysis of the covariates
data_covs <- data[, cols_covariates]
data_covs$Z <- as.numeric(data_covs$Z)-1
data_covs$HANDEDNESS_CATEGORY <- as.numeric(factor(data_covs$HANDEDNESS_CATEGORY))-1
data_covs$CURRENT_MED_STATUS <- as.numeric(factor(data_covs$CURRENT_MED_STATUS))-1
xbal <- xBalance(Z ~ ., data=data_covs)
print(xbal)
par(mfrow=c(1,1))
ggplot(data_covs, aes(x=FIQ, fill=factor(Z))) + geom_histogram(alpha=0.5, position="identity", bins=20) + labs(title="FIQ", x="FIQ", y="Count") + scale_fill_discrete(name="Treatment")
ggplot(data_covs, aes(x=PIQ, fill=factor(Z))) + geom_histogram(alpha=0.5, position="identity", bins=20) + labs(title="PIQ", x="PIQ", y="Count") + scale_fill_discrete(name="Treatment")

# Compute propensity score
formula_covs <- Z ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + CURRENT_MED_STATUS + compressed_3_1 + compressed_3_2 + compressed_3_3
ps <- glm(formula_covs, data=data_covs, family=binomial())
data$prop <- ps$fitted.values

# Plot propensity score distributions 
par(mfrow=c(1,1))
ggplot(data_covs, aes(x=prop, fill=factor(Z))) + geom_density(alpha=0.5) + labs(title="Propensity Score", x="Propensity Score", y="Density") + scale_fill_discrete(name="Treatment")


### AIPW ###

# Set seed for reproducibility
set.seed(12345)

data_covs <- data[, cols_covariates]
n <- nrow(data_covs)
data_covs_treatment <- data_covs[data_covs$Z == 1,]
data_covs_control <- data_covs[data_covs$Z == 0,]
n_treatment_1 <- nrow(data_covs_treatment)%/%2
n_treatment_2 <- nrow(data_covs_treatment) - n_treatment_1
n_control_1 <- nrow(data_covs_control)%/%2
n_control_2 <- nrow(data_covs_control) - n_control_1
idx_treatment <- sample(1:nrow(data_covs_treatment), n_treatment_1)
idx_control <- sample(1:nrow(data_covs_control), n_control_1)
data_covs_1_treatment <- data_covs_treatment[idx_treatment,]
data_covs_2_treatment <- data_covs_treatment[-idx_treatment,]
data_covs_1_control <- data_covs_control[idx_control,]
data_covs_2_control <- data_covs_control[-idx_control,]

# Set the Y variable 
Y_name <- "FIQ"

X_1_treatment <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_1_treatment) %>% scale(center=T, scale=F)
ncol(X_1_treatment)
Y <- data_covs_1_treatment[, Y_name]
mu_1 <- regression_forest(X_1_treatment,Y)

X_1_control <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_1_control) %>% scale(center=T, scale=F)
ncol(X_1_control)
Y <- data_covs_1_control[, Y_name]
mu_1_control <- regression_forest(X_1_control,Y)

X_2_treatment <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_2_treatment) %>% scale(center=T, scale=F)
ncol(X_2_treatment)
Y <- data_covs_2_treatment[, Y_name]
mu_2 <- regression_forest(X_2_treatment,Y)

X_2_control <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_2_control) %>% scale(center=T, scale=F)
ncol(X_2_control)
Y <- data_covs_2_control[, Y_name]
mu_2_control <- regression_forest(X_2_control,Y)

mu_tilde_1_control <- predict(mu_2, X_1_control)+(1/n_treatment_1)*sum(data_covs_1_treatment[,Y_name]-predict(mu_2, X_1_treatment))
mu_tilde_1 <- predict(mu_2, X_1_treatment)+(1/n_control_1)*sum(data_covs_1_control[,Y_name]-predict(mu_2, X_1_control))

mu_tilde_2_control <- predict(mu_1, X_2_control)+(1/n_treatment_2)*sum(data_covs_2_treatment[,Y_name]-predict(mu_1, X_2_treatment))
mu_tilde_2 <- predict(mu_1, X_2_treatment)+(1/n_control_2)*sum(data_covs_2_control[,Y_name]-predict(mu_1, X_2_control))


mu_1_dr_1 <- (1/(n_treatment_1+n_control_1))*sum((data_covs_1_treatment[,Y_name]- mu_tilde_1)/data_covs_1_treatment$prop + mu_tilde_1)
mu_0_dr_1 <- (1/(n_treatment_1+n_control_1))*sum((data_covs_1_control[,Y_name]- mu_tilde_1_control)/data_covs_1_control$prop + mu_tilde_1_control)
tau_1 <- mu_1_dr_1 - mu_0_dr_1

mu_1_dr_2 <- (1/(n_treatment_2+n_control_2))*sum((data_covs_2_treatment[,Y_name]- mu_tilde_2)/data_covs_2_treatment$prop + mu_tilde_2)
mu_0_dr_2 <- (1/(n_treatment_2+n_control_2))*sum((data_covs_2_control[,Y_name]- mu_tilde_2_control)/data_covs_2_control$prop + mu_tilde_2_control)
tau_2 <- mu_1_dr_2 - mu_0_dr_2

tau <- ((n_control_1+n_treatment_1)/nrow(data_covs))*tau_1 + ((n_control_2+n_treatment_2)/nrow(data_covs))*tau_2
cat("tau: ", tau, "\n")




# AIPW estimation function
compute_AIPW <- function(data, Y_name, cols_covariates) {
    # Splitting data into treatment and control
    data_covs <- data[, cols_covariates]
    n <- nrow(data_covs)
    data_covs_treatment <- data_covs[data_covs$Z == 1,]
    data_covs_control <- data_covs[data_covs$Z == 0,]
    n_treatment_1 <- nrow(data_covs_treatment)%/%2
    n_treatment_2 <- nrow(data_covs_treatment) - n_treatment_1
    n_control_1 <- nrow(data_covs_control)%/%2
    n_control_2 <- nrow(data_covs_control) - n_control_1
    idx_treatment <- sample(1:nrow(data_covs_treatment), n_treatment_1)
    idx_control <- sample(1:nrow(data_covs_control), n_control_1)
    data_covs_1_treatment <- data_covs_treatment[idx_treatment,]
    data_covs_2_treatment <- data_covs_treatment[-idx_treatment,]
    data_covs_1_control <- data_covs_control[idx_control,]
    data_covs_2_control <- data_covs_control[-idx_control,]

    X_1_treatment <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)++compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_1_treatment) %>% scale(center=T, scale=F)
    ncol(X_1_treatment)
    Y <- data_covs_1_treatment[, Y_name]
    mu_1 <- regression_forest(X_1_treatment,Y)

    X_1_control <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)++compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_1_control) %>% scale(center=T, scale=F)
    ncol(X_1_control)
    Y <- data_covs_1_control[, Y_name]
    mu_1_control <- regression_forest(X_1_control,Y)

    X_2_treatment <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_2_treatment) %>% scale(center=T, scale=F)
    ncol(X_2_treatment)
    Y <- data_covs_2_treatment[, Y_name]
    mu_2 <- regression_forest(X_2_treatment,Y)

    X_2_control <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_2_control) %>% scale(center=T, scale=F)
    ncol(X_2_control)
    Y <- data_covs_2_control[, Y_name]
    mu_2_control <- regression_forest(X_2_control,Y)

    mu_tilde_1_control <- predict(mu_2, X_1_control)+(1/n_treatment_1)*sum(data_covs_1_treatment[,Y_name]-predict(mu_2, X_1_treatment))
    mu_tilde_1 <- predict(mu_2, X_1_treatment)+(1/n_control_1)*sum(data_covs_1_control[,Y_name]-predict(mu_2, X_1_control))

    mu_tilde_2_control <- predict(mu_1, X_2_control)+(1/n_treatment_2)*sum(data_covs_2_treatment[,Y_name]-predict(mu_1, X_2_treatment))
    mu_tilde_2 <- predict(mu_1, X_2_treatment)+(1/n_control_2)*sum(data_covs_2_control[,Y_name]-predict(mu_1, X_2_control))


    mu_1_dr_1 <- (1/(n_treatment_1+n_control_1))*sum((data_covs_1_treatment[,Y_name]- mu_tilde_1)/data_covs_1_treatment$prop + mu_tilde_1)
    mu_0_dr_1 <- (1/(n_treatment_1+n_control_1))*sum((data_covs_1_control[,Y_name]- mu_tilde_1_control)/data_covs_1_control$prop + mu_tilde_1_control)
    tau_1 <- mu_1_dr_1 - mu_0_dr_1

    mu_1_dr_2 <- (1/(n_treatment_2+n_control_2))*sum((data_covs_2_treatment[,Y_name]- mu_tilde_2)/data_covs_2_treatment$prop + mu_tilde_2)
    mu_0_dr_2 <- (1/(n_treatment_2+n_control_2))*sum((data_covs_2_control[,Y_name]- mu_tilde_2_control)/data_covs_2_control$prop + mu_tilde_2_control)
    tau_2 <- mu_1_dr_2 - mu_0_dr_2

    tau <- ((n_control_1+n_treatment_1)/nrow(data_covs))*tau_1 + ((n_control_2+n_treatment_2)/nrow(data_covs))*tau_2
  
    return(tau)
}

# Main Analysis
# Assuming 'data' is your dataframe and 'FIQ' is the outcome variable
Y_name <- "PIQ"
cols_covariates <- c("Z", "FIQ", "VIQ", "PIQ", "SEX", "AGE_AT_SCAN", "HANDEDNESS_CATEGORY", cols_mri, "prop")
tau <- compute_AIPW(data, Y_name, cols_covariates)
cat("AIPW estimate: ", tau, "\n")

# Bootstrapping for variance estimation
n_bootstraps <- 1000
bootstrap_estimates <- numeric(n_bootstraps)
for(b in 1:n_bootstraps) {
    # Resample data
    boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
    # Compute AIPW estimate with error handling
    result <- try(compute_AIPW(boot_data, Y_name, cols_covariates), silent = TRUE)

    if(class(result) == "try-error") {
        cat("Error in iteration", b, "\n")  # Optionally log error
        next  # Skip to the next iteration
    } else {
        bootstrap_estimates[b] <- result
    }
}

# Calculate variance of successful estimates
successful_estimates <- bootstrap_estimates[bootstrap_estimates != 0]
bootstrap_variance <- var(successful_estimates)
cat("Bootstrap variance: ", bootstrap_variance, "\n")

# Print the AIPW estimate
cat("AIPW estimate: ", tau, "\n")

# Print the 95% confidence interval
CI <- c(tau - 1.96*sqrt(bootstrap_variance), tau + 1.96*sqrt(bootstrap_variance))
cat("95% confidence interval: ", CI, "\n")


