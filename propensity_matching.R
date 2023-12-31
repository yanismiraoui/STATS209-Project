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
cols_covariates <- c("Z", "FIQ", "VIQ", "PIQ", "SEX", "AGE_AT_SCAN", "HANDEDNESS_CATEGORY", cols_mri)
cols_covariates

### PROPENSITY SCORE MATCHING ###

# Analysis of the covariates
data_covs <- data[, cols_covariates]
data_covs$Z <- as.numeric(data_covs$Z)-1
data_covs$HANDEDNESS_CATEGORY <- as.numeric(factor(data_covs$HANDEDNESS_CATEGORY))-1
xbal <- xBalance(Z ~ ., data=data_covs)
print(xbal)
par(mfrow=c(1,1))
ggplot(data_covs, aes(x=FIQ, fill=factor(Z))) + geom_histogram(alpha=0.5, position="identity", bins=20) + labs(title="FIQ", x="FIQ", y="Count") + scale_fill_discrete(name="Treatment")
ggplot(data_covs, aes(x=PIQ, fill=factor(Z))) + geom_histogram(alpha=0.5, position="identity", bins=20) + labs(title="PIQ", x="PIQ", y="Count") + scale_fill_discrete(name="Treatment")

# Compute propensity score
formula_covs <- Z ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + compressed_3_1 + compressed_3_2 + compressed_3_3
ps <- glm(formula_covs, data=data_covs, family=binomial())
data_covs$prop <- ps$fitted.values

# Plot propensity score distributions 
par(mfrow=c(1,1))

ggplot() + 
  geom_histogram(aes(x = data_covs[data_covs$Z == 1, "prop"], y = after_stat(density), 
                     fill = "Treatment"), 
                 binwidth = 0.01, 
                 color = "black", 
                 alpha = 0.5) +
  geom_histogram(aes(x = data_covs[data_covs$Z == 0, "prop"], y = after_stat(density), 
                     fill = "Control"), 
                 binwidth = 0.01,  
                 color = "black", 
                 alpha = 0.5) +
  labs(x = "Propensit", y = "Frequency") + 
  scale_fill_manual(values = c("Treatment" = "firebrick", "Control" = "lightblue")) +
  guides(fill = guide_legend(title = "Group"))


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
formula_plot <- Z ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + compressed_3_1 + compressed_3_2 + compressed_3_3 -1
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
plot(xBalance(formula_plot, strata=list(unadj=NULL, matches=~pairmatch.1),data=data_covs),ggplot = TRUE)
plot(xBalance(formula_plot, strata=list(unadj=NULL, matches=~pairmatch.2),data=data_covs),ggplot = TRUE)
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
cols_0 <- c('SEX.0', 'AGE_AT_SCAN.0', 'HANDEDNESS_CATEGORY.0', 'compressed_3_1.0', 'compressed_3_2.0', 'compressed_3_3.0', 'prop.0', 'FIQ.0', 'FIQ.1')
cols_1 <- c('SEX.1', 'AGE_AT_SCAN.1', 'HANDEDNESS_CATEGORY.1', 'compressed_3_1.1', 'compressed_3_2.1', 'compressed_3_3.1', 'prop.1', 'FIQ.0', 'FIQ.1')
control <- summarize.2[, cols_0]
colnames(control) <- c('SEX', 'AGE_AT_SCAN', 'HANDEDNESS_CATEGORY', 'compressed_3_1', 'compressed_3_2', 'compressed_3_3', 'prop', 'FIQ.0', 'FIQ.1')
treated <- summarize.2[, cols_1]
colnames(treated) <- c('SEX', 'AGE_AT_SCAN', 'HANDEDNESS_CATEGORY', 'compressed_3_1', 'compressed_3_2', 'compressed_3_3', 'prop', 'FIQ.0', 'FIQ.1')

mu_hat_0 <- lm(FIQ.0 ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + compressed_3_1 + compressed_3_2 + compressed_3_3, data=control)
mu_hat_1 <- lm(FIQ.1 ~ SEX + AGE_AT_SCAN + HANDEDNESS_CATEGORY + compressed_3_1 + compressed_3_2 + compressed_3_3, data=treated)

bias <- predict(mu_hat_0, treated) - predict(mu_hat_0, control)
tau <-  mean(summarize.2$FIQ.1 - summarize.2$FIQ.0)
tau_corr  <- tau - mean(bias)
cat("Bias-corrected estimate:", tau_corr, "\n")

# Variance estimate
V <- 1/nrow(summarize.2)**2 * (sum((treated$FIQ.1 - predict(mu_hat_1, treated))**2) + sum((control$FIQ.0 - predict(mu_hat_0, control))**2))
cat("Variance:", V)

# Confidence interval
CI <- c(tau_corr - 1.96*sqrt(V), tau_corr + 1.96*sqrt(V))
cat("Confidence interval:", CI, "\n")







