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

# Covariates: SEX, AGE_AT_SCAN, HANDEDNESS, CURRENT_MED_STATUS, MRI_FEATURES
cols_covariates <- c("Z", "FIQ", "VIQ", "PIQ", "SEX", "AGE_AT_SCAN", "HANDEDNESS_CATEGORY", "CURRENT_MED_STATUS", cols_mri)
cols_covariates

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

# Create random forest for predicting mu_0, mu_1, mu_2 and mu_3
cols <- c("female", "english", "hsgroup", "numcourses_nov1")
data_ssp <- data[data$ssp == 1| data$control == 1,]
data_sfp <- data[data$sfp == 1| data$control == 1,]
data_sfsp <- data[data$sfsp == 1| data$control == 1,]

# Split each dataset into 2 halves
n <- nrow(data_ssp)
data_ssp_treatment <- data_ssp[data_ssp$ssp == 1,]
data_ssp_control <- data_ssp[data_ssp$control == 1,]
n_treatment_1 <- nrow(data_ssp_treatment)%/%2
n_treatment_2 <- nrow(data_ssp_treatment) - n_treatment_1
n_control_1 <- nrow(data_ssp_control)%/%2
n_control_2 <- nrow(data_ssp_control) - n_control_1
idx_treatment <- sample(1:nrow(data_ssp_treatment), n_treatment_1)
idx_control <- sample(1:nrow(data_ssp_control), n_control_1)
data_ssp_1_treatment <- data_ssp_treatment[idx_treatment,]
data_ssp_2_treatment <- data_ssp_treatment[-idx_treatment,]
data_ssp_1_control <- data_ssp_control[idx_control,]
data_ssp_2_control <- data_ssp_control[-idx_control,]

X_ssp_1_treatment <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_ssp_1_treatment) %>% scale(center=T, scale=F)
Y <- data_ssp_1_treatment$grade_20059_fall
mu_1_ssp <- regression_forest(X_ssp_1_treatment,Y)

X_ssp_1_control <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_ssp_1_control) %>% scale(center=T, scale=F)
Y <- data_ssp_1_control$grade_20059_fall
mu_1_control <- regression_forest(X_ssp_1_control,Y)

X_ssp_2_treatment <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_ssp_2_treatment) %>% scale(center=T, scale=F)
Y <- data_ssp_2_treatment$grade_20059_fall
mu_2_ssp <- regression_forest(X_ssp_2_treatment,Y)

X_ssp_2_control <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_ssp_2_control) %>% scale(center=T, scale=F)
Y <- data_ssp_2_control$grade_20059_fall
mu_2_control <- regression_forest(X_ssp_2_control,Y)

mu_tilde_1_control <- predict(mu_2_ssp, X_ssp_1_control)+(1/n_treatment_1)*sum(data_ssp_1_treatment$grade_20059_fall-predict(mu_2_ssp, X_ssp_1_treatment))
mu_tilde_1_ssp <- predict(mu_2_ssp, X_ssp_1_treatment)+(1/n_control_1)*sum(data_ssp_1_control$grade_20059_fall-predict(mu_2_ssp, X_ssp_1_control))
tau_1 <- (1/(n_control_1+n_treatment_1))*((sum(data_ssp_1_treatment$grade_20059_fall) + sum(mu_tilde_1_control)) - (sum(data_ssp_1_control$grade_20059_fall) + sum(mu_tilde_1_ssp)))

mu_tilde_2_control <- predict(mu_1_ssp, X_ssp_2_control)+(1/n_treatment_2)*sum(data_ssp_2_treatment$grade_20059_fall-predict(mu_1_ssp, X_ssp_2_treatment))
mu_tilde_2_ssp <- predict(mu_1_ssp, X_ssp_2_treatment)+(1/n_control_2)*sum(data_ssp_2_control$grade_20059_fall-predict(mu_1_ssp, X_ssp_2_control))
tau_2 <- (1/(n_control_2+n_treatment_2))*((sum(data_ssp_2_treatment$grade_20059_fall) + sum(mu_tilde_2_control)) - (sum(data_ssp_2_control$grade_20059_fall) + sum(mu_tilde_2_ssp)))
tau_ssp <- ((n_control_1+n_treatment_1)/nrow(data_ssp))*tau_1 + ((n_control_2+n_treatment_2)/nrow(data_ssp))*tau_2
cat("tau_ssp: ", tau_ssp, "\n")
var_1 <- (1/n_treatment_1)*(1/(n_treatment_1-1))*sum((data_ssp_1_treatment$grade_20059_fall-predict(mu_2_ssp, X_ssp_1_treatment))^2) + (1/n_control_1)*(1/(n_control_1-1))*sum((data_ssp_1_control$grade_20059_fall-predict(mu_2_ssp, X_ssp_1_control))^2) + 1/(n_treatment_1+n_control_1)*var(mu_tilde_2_ssp[-nrow(mu_tilde_2_ssp),] - mu_tilde_1_ssp)
var_2 <- (1/n_treatment_2)*(1/(n_treatment_2-1))*sum((data_ssp_2_treatment$grade_20059_fall-predict(mu_1_ssp, X_ssp_2_treatment))^2) + (1/n_control_2)*(1/(n_control_2-1))*sum((data_ssp_2_control$grade_20059_fall-predict(mu_1_ssp, X_ssp_2_control))^2) + 1/(n_treatment_2+n_control_2)*var(mu_tilde_2_ssp[-nrow(mu_tilde_2_ssp),] - mu_tilde_1_ssp)
var <- ((n_treatment_1+n_control_1)/n)^2*var_1 + ((n_treatment_2+n_control_2)/n)^2*var_2

# Confidence interval for tau_ssp
CI <- c(tau_ssp - 1.96*sqrt(var), tau_ssp + 1.96*sqrt(var))
CI

# Doing the same thing as above but for sfp
data_sfp_treatment <- data_sfp[data_sfp$sfp == 1,]
data_sfp_control <- data_sfp[data_sfp$control == 1,]
n_treatment_1 <- nrow(data_sfp_treatment)%/%2
n_treatment_2 <- nrow(data_sfp_treatment) - n_treatment_1
n_control_1 <- nrow(data_sfp_control)%/%2
n_control_2 <- nrow(data_sfp_control) - n_control_1
idx_treatment <- sample(1:nrow(data_sfp_treatment), n_treatment_1)
idx_control <- sample(1:nrow(data_sfp_control), n_control_1)
data_sfp_1_treatment <- data_sfp_treatment[idx_treatment,]
data_sfp_2_treatment <- data_sfp_treatment[-idx_treatment,]
data_sfp_1_control <- data_sfp_control[idx_control,]
data_sfp_2_control <- data_sfp_control[-idx_control,]

X_sfp_1_treatment <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfp_1_treatment) %>% scale(center=T, scale=F)
Y <- data_sfp_1_treatment$grade_20059_fall
mu_1_sfp <- regression_forest(X_sfp_1_treatment,Y)

X_sfp_1_control <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfp_1_control) %>% scale(center=T, scale=F)
Y <- data_sfp_1_control$grade_20059_fall
mu_1_control <- regression_forest(X_sfp_1_control,Y, num.trees = 200)

X_sfp_2_treatment <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfp_2_treatment) %>% scale(center=T, scale=F)
Y <- data_sfp_2_treatment$grade_20059_fall
mu_2_sfp <- regression_forest(X_sfp_2_treatment,Y)

X_sfp_2_control <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfp_2_control) %>% scale(center=T, scale=F)
Y <- data_sfp_2_control$grade_20059_fall
mu_2_control <- regression_forest(X_sfp_2_control,Y)

mu_tilde_1_control <- predict(mu_2_sfp, X_sfp_1_control)+(1/n_treatment_1)*sum(data_sfp_1_treatment$grade_20059_fall-predict(mu_2_sfp, X_sfp_1_treatment))
mu_tilde_1_sfp <- predict(mu_2_sfp, X_sfp_1_treatment)+(1/n_control_1)*sum(data_sfp_1_control$grade_20059_fall-predict(mu_2_sfp, X_sfp_1_control))
tau_1 <- (1/(n_control_1+n_treatment_1))*((sum(data_sfp_1_treatment$grade_20059_fall) + sum(mu_tilde_1_control)) - (sum(data_sfp_1_control$grade_20059_fall) + sum(mu_tilde_1_sfp)))
mu_tilde_2_control <- predict(mu_1_sfp, X_sfp_2_control)+(1/n_treatment_2)*sum(data_sfp_2_treatment$grade_20059_fall-predict(mu_1_sfp, X_sfp_2_treatment))
mu_tilde_2_sfp <- predict(mu_1_sfp, X_sfp_2_treatment)+(1/n_control_2)*sum(data_sfp_2_control$grade_20059_fall-predict(mu_1_sfp, X_sfp_2_control))
tau_2 <- (1/(n_control_2+n_treatment_2))*((sum(data_sfp_2_treatment$grade_20059_fall) + sum(mu_tilde_2_control)) - (sum(data_sfp_2_control$grade_20059_fall) + sum(mu_tilde_2_sfp)))
tau_sfp <- ((n_control_1+n_treatment_1)/nrow(data_sfp))*tau_1 + ((n_control_2+n_treatment_2)/nrow(data_sfp))*tau_2
cat("tau_sfp: ", tau_sfp, "\n")

var_1 <- (1/n_treatment_1)*(1/(n_treatment_1-1))*sum((data_sfp_1_treatment$grade_20059_fall-predict(mu_2_sfp, X_sfp_1_treatment))^2) + (1/n_control_1)*(1/(n_control_1-1))*sum((data_sfp_1_control$grade_20059_fall-predict(mu_2_sfp, X_sfp_1_control))^2) + 1/(n_treatment_1+n_control_1)*var(mu_tilde_2_sfp[-nrow(mu_tilde_2_sfp),] - mu_tilde_1_sfp)
var_2 <- (1/n_treatment_2)*(1/(n_treatment_2-1))*sum((data_sfp_2_treatment$grade_20059_fall-predict(mu_1_sfp, X_sfp_2_treatment))^2) + (1/n_control_2)*(1/(n_control_2-1))*sum((data_sfp_2_control$grade_20059_fall-predict(mu_1_sfp, X_sfp_2_control))^2) + 1/(n_treatment_2+n_control_2)*var(mu_tilde_2_sfp[-nrow(mu_tilde_2_sfp),] - mu_tilde_1_sfp)
var <- ((n_treatment_1+n_control_1)/n)^2*var_1 + ((n_treatment_2+n_control_2)/n)^2*var_2

# Confidence interval for tau_sfp
CI <- c(tau_sfp - 1.96*sqrt(var), tau_sfp + 1.96*sqrt(var))
CI

# Doing the same thing as above but for sfsp
data_sfsp_treatment <- data_sfsp[data_sfsp$sfsp == 1,]
data_sfsp_control <- data_sfsp[data_sfsp$control == 1,]
n_treatment_1 <- nrow(data_sfsp_treatment)%/%2
n_treatment_2 <- nrow(data_sfsp_treatment) - n_treatment_1
n_control_1 <- nrow(data_sfsp_control)%/%2
n_control_2 <- nrow(data_sfsp_control) - n_control_1
idx_treatment <- sample(1:nrow(data_sfsp_treatment), n_treatment_1)
idx_control <- sample(1:nrow(data_sfsp_control), n_control_1)
data_sfsp_1_treatment <- data_sfsp_treatment[idx_treatment,]
data_sfsp_2_treatment <- data_sfsp_treatment[-idx_treatment,]
data_sfsp_1_control <- data_sfsp_control[idx_control,]
data_sfsp_2_control <- data_sfsp_control[-idx_control,]

X_sfsp_1_treatment <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfsp_1_treatment) %>% scale(center=T, scale=F)
Y <- data_sfsp_1_treatment$grade_20059_fall
mu_1_sfsp <- regression_forest(X_sfsp_1_treatment,Y)

X_sfsp_1_control <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfsp_1_control) %>% scale(center=T, scale=F)
Y <- data_sfsp_1_control$grade_20059_fall
mu_1_control <- regression_forest(X_sfsp_1_control,Y)

X_sfsp_2_treatment <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfsp_2_treatment) %>% scale(center=T, scale=F)
Y <- data_sfsp_2_treatment$grade_20059_fall
mu_2_sfsp <- regression_forest(X_sfsp_2_treatment,Y)

X_sfsp_2_control <- model.matrix( ~ factor(female)+factor(english)+factor(hsgroup)+numcourses_nov1, data_sfsp_2_control) %>% scale(center=T, scale=F)
Y <- data_sfsp_2_control$grade_20059_fall
mu_2_control <- regression_forest(X_sfsp_2_control,Y)

mu_tilde_1_control <- predict(mu_2_sfsp, X_sfsp_1_control)+(1/n_treatment_1)*sum(data_sfsp_1_treatment$grade_20059_fall-predict(mu_2_sfsp, X_sfsp_1_treatment))
mu_tilde_1_sfsp <- predict(mu_2_sfsp, X_sfsp_1_treatment)+(1/n_control_1)*sum(data_sfsp_1_control$grade_20059_fall-predict(mu_2_sfsp, X_sfsp_1_control))
tau_1 <- (1/(n_control_1+n_treatment_1))*((sum(data_sfsp_1_treatment$grade_20059_fall) + sum(mu_tilde_1_control)) - (sum(data_sfsp_1_control$grade_20059_fall) + sum(mu_tilde_1_sfsp)))
mu_tilde_2_control <- predict(mu_1_sfsp, X_sfsp_2_control)+(1/n_treatment_2)*sum(data_sfsp_2_treatment$grade_20059_fall-predict(mu_1_sfsp, X_sfsp_2_treatment))
mu_tilde_2_sfsp <- predict(mu_1_sfsp, X_sfsp_2_treatment)+(1/n_control_2)*sum(data_sfsp_2_control$grade_20059_fall-predict(mu_1_sfsp, X_sfsp_2_control))
tau_2 <- (1/(n_control_2+n_treatment_2))*((sum(data_sfsp_2_treatment$grade_20059_fall) + sum(mu_tilde_2_control)) - (sum(data_sfsp_2_control$grade_20059_fall) + sum(mu_tilde_2_sfsp)))
tau_sfsp <- ((n_control_1+n_treatment_1)/nrow(data_sfsp))*tau_1 + ((n_control_2+n_treatment_2)/nrow(data_sfsp))*tau_2
cat("tau_sfp: ", tau_sfsp, "\n")

var_1 <- (1/n_treatment_1)*(1/(n_treatment_1-1))*sum((data_sfsp_1_treatment$grade_20059_fall-predict(mu_2_sfsp, X_sfsp_1_treatment))^2) + (1/n_control_1)*(1/(n_control_1-1))*sum((data_sfsp_1_control$grade_20059_fall-predict(mu_2_sfsp, X_sfsp_1_control))^2) + 1/(n_treatment_1+n_control_1)*var(mu_tilde_2_sfsp[-nrow(mu_tilde_2_sfsp),] - mu_tilde_1_sfsp)
var_2 <- (1/n_treatment_2)*(1/(n_treatment_2-1))*sum((data_sfsp_2_treatment$grade_20059_fall-predict(mu_1_sfsp, X_sfsp_2_treatment))^2) + (1/n_control_2)*(1/(n_control_2-1))*sum((data_sfsp_2_control$grade_20059_fall-predict(mu_1_sfsp, X_sfsp_2_control))^2) + 1/(n_treatment_2+n_control_2)*var(mu_tilde_2_sfsp[-nrow(mu_tilde_2_sfsp),] - mu_tilde_1_sfsp)
var <- ((n_treatment_1+n_control_1)/n)^2*var_1 + ((n_treatment_2+n_control_2)/n)^2*var_2

# Confidence interval for tau_sfsp
CI <- c(tau_sfsp - 1.96*sqrt(var), tau_sfsp + 1.96*sqrt(var))
CI



