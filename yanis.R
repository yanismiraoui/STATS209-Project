setwd("/Users/yanis/Documents/Stanford/STATS 209/Project/")
getwd()

library(estimatr)
library(tidyverse)
library(magrittr)
library(ggplot2)

# Load data csv file
data <- read.csv("./data/ABIDE_tab_compressed.csv", header = TRUE)
nrow(data)
colnames(data)
summary(data)
data$Z <- as.numeric(data$DX_GROUP) - 1
data$Z <- factor(data$Z)
# Covariates: SEX, AGE_AT_SCAN, HANDEDNESS, CURRENT_MED_STATUS, MRI_FEATURES
cols_covariates <- c("Z", "FIQ", "VIQ", "PIQ", "SEX", "AGE_AT_SCAN", "HANDEDNESS_CATEGORY", "CURRENT_MED_STATUS", cols_mri)
cols_covariates

# Plot the distribution of subjects diagnosed with ASD (using red) and controls (using blue)
ggplot(data, aes(x=Z, fill=Z)) + geom_bar() + scale_fill_manual(values=c("blue", "red")) + labs(x="Diagnosis", y="Number of subjects", title="Distribution of subjects diagnosed with ASD and controls", fill="Diagnosis", caption="Source: ABIDE")

# Number of students in each arm
table(data$Z)

# count number of na values in FIQ, VIQ and PIQ
sum(is.na(data$FIQ))
sum(is.na(data$VIQ))
sum(is.na(data$PIQ))

# get rid of rows with na values in FIQ, VIQ and PIQ
data <- data[!is.na(data$FIQ),]
data <- data[!is.na(data$VIQ),]
data <- data[!is.na(data$PIQ),]
data <- data[data$FIQ != -9999,]
data <- data[data$VIQ != -9999,]
data <- data[data$PIQ != -9999,]
data <- data[data$AGE_AT_SCAN != -9999,]
data <- data[data$HANDEDNESS_CATEGORY != -9999,]
data <- data[data$CURRENT_MED_STATUS != -9999,]
cat("Number of rows: ", nrow(data), "\n")

# Keep only the columns where we have MRI data
data <- data[!is.na(data$compressed_2_1),]
cat("Number of rows: ", nrow(data), "\n")

# get rid of rows with '`' in column CURRENT_MED_STATUS
data <- data[!grepl("`", data$CURRENT_MED_STATUS),]
cat("Number of rows: ", nrow(data), "\n")

table(data$Z)
ggplot(data, aes(x=Z, fill=Z)) + geom_bar() + scale_fill_manual(values=c("blue", "red")) + labs(x="Diagnosis", y="Number of subjects", title="Distribution of subjects diagnosed with ASD and controls", fill="Diagnosis", caption="Source: ABIDE")

# Plot distribution of scores for FIQ, VIQ and PIQ and for Z=0 and Z=1
ggplot(data, aes(x=FIQ, fill=Z)) + geom_density(alpha=0.5) + scale_fill_manual(values=c("blue", "red")) + labs(x="FIQ", y="Density", title="Distribution of FIQ scores", fill="Diagnosis", caption="Source: ABIDE")
ggplot(data, aes(x=VIQ, fill=Z)) + geom_density(alpha=0.5) + scale_fill_manual(values=c("blue", "red")) + labs(x="VIQ", y="Density", title="Distribution of FIQ scores", fill="Diagnosis", caption="Source: ABIDE")
ggplot(data, aes(x=PIQ, fill=Z)) + geom_density(alpha=0.5) + scale_fill_manual(values=c("blue", "red")) + labs(x="PIQ", y="Density", title="Distribution of FIQ scores", fill="Diagnosis", caption="Source: ABIDE")

# Covariates: SEX, AGE_AT_SCAN, HANDEDNESS, CURRENT_MED_STATUS, MRI_FEATURES
SIZE_COMPRESSED_MRI <- 3
# Create a vector of the columns with names f"compressed_{SIZE_COMPRESSED_MRI}_{i}" for i in 1:SIZE_COMPRESSED_MRI
cols_mri <- c()
for (i in 1:SIZE_COMPRESSED_MRI) {
  cols_mri <- c(cols_mri, paste0("compressed_", SIZE_COMPRESSED_MRI, "_", i))
}
cols_mri

# Using the covariates, we can create a linear model to predict the IQ of the subjects (FIQ, VIQ and PIQ)
# We use the ”basic controls” and Lin’s estimator to estimate the causal effect
data_covs <- data[, cols_covariates]
colnames(data_covs)
X <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs) %>% scale(center=T, scale=F)
Z <- data_covs$Z

Y_1 <- data_covs$FIQ
fit1 <- lm_robust(Y_1~Z+X+Z*X)
summary(fit1)

Y_2 <- data_covs$VIQ
fit2 <- lm_robust(Y_2~Z+X+Z*X)
summary(fit2)

Y_3 <- data_covs$PIQ
fit3 <- lm_robust(Y_3~Z+X+Z*X)
summary(fit3)





# Q6
library(grf)
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

X_1_treatment <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_1_treatment) %>% scale(center=T, scale=F)
ncol(X_1_treatment)
Y <- data_covs_1_treatment$FIQ
mu_1 <- regression_forest(X_1_treatment,Y)

X_1_control <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_1_control) %>% scale(center=T, scale=F)
ncol(X_1_control)
Y <- data_covs_1_control$FIQ
mu_1_control <- regression_forest(X_1_control,Y)

X_2_treatment <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_2_treatment) %>% scale(center=T, scale=F)
ncol(X_2_treatment)
Y <- data_covs_2_treatment$FIQ
mu_2 <- regression_forest(X_2_treatment,Y)

X_2_control <- model.matrix( ~ factor(SEX)+AGE_AT_SCAN+factor(HANDEDNESS_CATEGORY)+factor(CURRENT_MED_STATUS)+compressed_3_1+compressed_3_2+compressed_3_3 , data_covs_2_control) %>% scale(center=T, scale=F)
ncol(X_2_control)
Y <- data_covs_2_control$FIQ
mu_2_control <- regression_forest(X_2_control,Y)

mu_tilde_1_control <- predict(mu_2, X_1_control)+(1/n_treatment_1)*sum(data_covs_1_treatment$FIQ-predict(mu_2, X_1_treatment))
mu_tilde_1 <- predict(mu_2, X_1_treatment)+(1/n_control_1)*sum(data_covs_1_control$FIQ-predict(mu_2, X_1_control))
tau_1 <- (1/(n_control_1+n_treatment_1))*((sum(data_covs_1_treatment$FIQ) + sum(mu_tilde_1_control)) - (sum(data_covs_1_control$FIQ) + sum(mu_tilde_1)))

mu_tilde_2_control <- predict(mu_1, X_2_control)+(1/n_treatment_2)*sum(data_covs_2_treatment$FIQ-predict(mu_1, X_2_treatment))
mu_tilde_2 <- predict(mu_1, X_2_treatment)+(1/n_control_2)*sum(data_covs_2_control$FIQ-predict(mu_1, X_2_control))
tau_2 <- (1/(n_control_2+n_treatment_2))*((sum(data_covs_2_treatment$FIQ) + sum(mu_tilde_2_control)) - (sum(data_covs_2_control$FIQ) + sum(mu_tilde_2)))

tau <- ((n_control_1+n_treatment_1)/nrow(data_covs))*tau_1 + ((n_control_2+n_treatment_2)/nrow(data_covs))*tau_2
cat("tau: ", tau, "\n")
var_1 <- (1/n_treatment_1)*(1/(n_treatment_1-1))*sum((data_covs_1_treatment$FIQ-predict(mu_2, X_1_treatment))^2) + (1/n_control_1)*(1/(n_control_1-1))*sum((data_covs_1_control$FIQ-predict(mu_2, X_1_control))^2) + 1/(n_treatment_1+n_control_1)*var(mu_tilde_2[-nrow(mu_tilde_2),] - mu_tilde_1)
var_2 <- (1/n_treatment_2)*(1/(n_treatment_2-1))*sum((data_covs_2_treatment$FIQ-predict(mu_1, X_2_treatment))^2) + (1/n_control_2)*(1/(n_control_2-1))*sum((data_covs_2_control$FIQ-predict(mu_1, X_2_control))^2) + 1/(n_treatment_2+n_control_2)*var(mu_tilde_2[-nrow(mu_tilde_2),] - mu_tilde_1)
var <- ((n_treatment_1+n_control_1)/n)^2*var_1 + ((n_treatment_2+n_control_2)/n)^2*var_2

# Confidence interval for tau
CI <- c(tau - 1.96*sqrt(var), tau + 1.96*sqrt(var))
CI
