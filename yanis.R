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

### ANALYSIS OF DATA ###

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

### BASIC CONTROLS AND LIN'S ESTIMATOR ###

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


### BASIC CONTROL AND ML REGRESSION ADJUSTMENT ###

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


### PROPENSITY SCORE MATCHING ###

# Analysis of the covariates
data_covs <- data[, cols_covariates]
data_covs$Z <- as.numeric(data_covs$Z)-1
xbal <- xBalance(Z ~ ., data=data_covs)
print(xbal)
par(mfrow=c(1,1))
ggplot(data_covs, aes(x=FIQ, fill=factor(Z))) + geom_histogram(alpha=0.5, position="identity", bins=20) + labs(title="FIQ", x="FIQ", y="Count") + scale_fill_discrete(name="Treatment")
ggplot(data_covs, aes(x=PIQ, fill=factor(Z))) + geom_histogram(alpha=0.5, position="identity", bins=20) + labs(title="PIQ", x="PIQ", y="Count") + scale_fill_discrete(name="Treatment")

# Compute propensity score
formula_covs <- Z ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + CURRENT_MED_STATUS + compressed_3_1 + compressed_3_2 + compressed_3_3
ps <- glm(formula_covs, data=data_covs, family=binomial())
data_covs$prop <- ps$fitted.values

# Plot propensity score distributions 
par(mfrow=c(1,1))
hist(data_covs$prop[data_covs$Z==1], main="Propensity Score", xlab="Propensity Score", freq=FALSE, col="black")
hist(data_covs$prop[data_covs$Z==0], add=TRUE, freq=FALSE, col="red", alpha=0.5)
legend("topright", c("Treated","Control"), fill=c("black","red"), cex=2)
# same plot with ggplot2 density
ggplot(data_covs, aes(x=prop, fill=factor(Z))) + geom_density(alpha=0.5) + labs(title="Propensity Score", x="Propensity Score", y="Density") + scale_fill_discrete(name="Treatment")

# Matching on the propensity score
data_covs <- na.omit(data_covs)
mat.1 <- match_on(formula_covs, data=data_covs)
pairmatch.1 <- pairmatch(mat.1, data=data_covs)
summarize.1 <- summarize.match(data_covs, pairmatch.1)
print(pairmatch.1, grouped = TRUE)
formula_plot <- Z ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + CURRENT_MED_STATUS + compressed_3_1 + compressed_3_2 + compressed_3_3 -1
plot(xBalance(formula_plot,strata=list(unstrat=NULL, ms.2=~pairmatch.1), data=data_covs),ggplot = TRUE)

# Average absolute difference in propensity scores within matched pairs
mean(abs(summarize.1$prop.0 - summarize.1$prop.1))
# Maximum absolute difference in propensity score
max(abs(summarize.1$prop.0 - summarize.1$prop.1))

# Matching on the propensity score with a caliper
mat.2 <- addcaliper(mat.1, z=data_covs$Z, p=data_covs$prop, caliper=0.1)
pairmatch.2 <- pairmatch(mat.2, data=data_covs)
summarize.2 <- summarize.match(data_covs, pairmatch.2)
plot(xBalance(formula_plot, strata=list(unstrat=NULL, ms.2=~pairmatch.2),data=data_covs),ggplot = TRUE)

## Comparing (1) and (2)
par(mfrow=c(1,2))
plot(xBalance(formula_plot, strata=list(unstrat=NULL, ms.2=~pairmatch.1),data=data_covs),ggplot = TRUE)
plot(xBalance(formula_plot, strata=list(unstrat=NULL, ms.2=~pairmatch.2),data=data_covs),ggplot = TRUE)
# Average absolute difference in propensity scores within caliper-matched pairs
mean(abs(summarize.2$prop.0 - summarize.2$prop.1))
# Maximum absolute difference in propensity score
max(abs(summarize.2$prop.0 - summarize.2$prop.1))


# Matching on the propensity score with a caliper and a ratio
mat.3 <- addalmostexact(mat.2, z=data_covs$Z, f=data_covs$SEX, mult=10)
pairmatch.3 <- pairmatch(mat.3, data=data_covs)
summarize.3 <- summarize.match(data_covs, pairmatch.3)
plot(xBalance(formula_plot, strata=list(unstrat=NULL, ms.2=~pairmatch.3),data=data_covs),ggplot = TRUE)

# Comparing (2) and (3)
mean(abs(summarize.3$prop.0 - summarize.3$prop.1))
max(abs(summarize.3$prop.0 - summarize.3$prop.1))


# P-value from a FRT of Fisher’s sharp null f
B  <- 100000
T_obs <- mean((summarize.2$FIQ.1-summarize.2$FIQ.0))
T_perm <- rep(0, B)
p_value <- 0
n <- length(summarize.2$FIQ.1)
for (i in 1:B) {
    # Permute Z
    Z <- rbinom(n, 1, p=0.5)
    # Compute T_perm
    T_perm[i] <- mean((2*Z-1)*(summarize.2$FIQ.1-summarize.2$FIQ.0))
    # Update p-value
    if (T_perm[i] >= T_obs) {
        p_value <- p_value + 1
    }
}
p_value <- p_value / B
cat("p-value = ", p_value, "\n")

B  <- 100000
T_obs <- mean((summarize.2$PIQ.1-summarize.2$PIQ.0))
T_perm <- rep(0, B)
p_value <- 0
n <- length(summarize.2$PIQ.1)
for (i in 1:B) {
    # Permute Z
    Z <- rbinom(n, 1, p=0.5)
    # Compute T_perm
    T_perm[i] <- mean((2*Z-1)*(summarize.2$PIQ.1-summarize.2$PIQ.0))
    # Update p-value
    if (T_perm[i] >= T_obs) {
        p_value <- p_value + 1
    }
}
p_value <- p_value / B
cat("p-value = ", p_value, "\n")

B  <- 100000
T_obs <- mean((summarize.2$VIQ.1-summarize.2$VIQ.0))
T_perm <- rep(0, B)
p_value <- 0
n <- length(summarize.2$VIQ.1)
for (i in 1:B) {
    # Permute Z
    Z <- rbinom(n, 1, p=0.5)
    # Compute T_perm
    T_perm[i] <- mean((2*Z-1)*(summarize.2$VIQ.1-summarize.2$VIQ.0))
    # Update p-value
    if (T_perm[i] >= T_obs) {
        p_value <- p_value + 1
    }
}
p_value <- p_value / B
cat("p-value = ", p_value, "\n")




# Bias-corrected estimate of the average treatment effect on the treated
colnames(data_covs)
cols_0 <- c('SEX.0', 'AGE_AT_SCAN.0', 'HANDEDNESS_CATEGORY.0', 'CURRENT_MED_STATUS.0', 'compressed_3_1.0', 'compressed_3_2.0', 'compressed_3_3.0', 'prop.0', 'FIQ.0', 'FIQ.1')
cols_1 <- c('SEX.1', 'AGE_AT_SCAN.1', 'HANDEDNESS_CATEGORY.1', 'CURRENT_MED_STATUS.1', 'compressed_3_1.1', 'compressed_3_2.1', 'compressed_3_3.1', 'prop.1', 'FIQ.0', 'FIQ.1')
control <- summarize.2[, cols_0]
colnames(control) <- c('SEX', 'AGE_AT_SCAN', 'HANDEDNESS_CATEGORY', 'CURRENT_MED_STATUS', 'compressed_3_1', 'compressed_3_2', 'compressed_3_3', 'prop', 'FIQ.0', 'FIQ.1')
treated <- summarize.2[, cols_1]
colnames(treated) <- c('SEX', 'AGE_AT_SCAN', 'HANDEDNESS_CATEGORY', 'CURRENT_MED_STATUS', 'compressed_3_1', 'compressed_3_2', 'compressed_3_3', 'prop', 'FIQ.0', 'FIQ.1')

mu_hat_0 <- lm(FIQ.0 ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + CURRENT_MED_STATUS + compressed_3_1 + compressed_3_2 + compressed_3_3, data=control)
mu_hat_1 <- lm(FIQ.1 ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + CURRENT_MED_STATUS + compressed_3_1 + compressed_3_2 + compressed_3_3, data=treated)

bias <- predict(mu_hat_0, treated) - predict(mu_hat_0, control)
tau <-  mean(summarize.2$FIQ.1 - summarize.2$FIQ.0)
tau_corr  <- tau - mean(bias)
cat("Bias-corrected estimate:", tau_corr, "\n")

# Variance estimate
V <- 1/nrow(summarize.2)**2 * (sum((treated$FIQ.1 - predict(mu_hat_1, treated))**2) + sum((control$FIQ.0 - predict(mu_hat_0, control))**2))
cat("Variance:", V)

