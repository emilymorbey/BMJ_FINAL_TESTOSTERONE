


### THIS CODE PROVIDES THE MENDELIAN RANDOMISATION ANALYSIS OF TESTOSTERONE AND SHBG LEVELS TO CARDIOVASCULAR DISEASE RISK

### 1. Mendelian Randomisation of Testosterone to Cardiovascular Disease Risk in Men 
### 2. Mendelian Randomisation of Tesotsterone to Cardiovascular Disease Risk in Women
### 3. Mendelian Randomisation of SHBG to Cardiovascular Disease Risk in Men
### 4. Mendelian Randomisation of SHBG to Cardiovascular Disease Risk in Women
### 5. MR-PRESSO analysis of Testosterone to Cardiovascular Disease Risk
### 6. MULTIVARIABLE MR FOR DISCOVERY OF POTENTIAL MEDIATORS OF TESTOSTERONE CAD RELATIONSHIP
### 7. Phenotyping and Survival Analysis of Testosterone in UK Biobank 








#####################################    1. Mendelian Randomisation of Testosterone to Cardiovascular Disease risk in Men    ########################################



library(tidyverse)
library(readxl)
library(MendelianRandomization)
library(data.table)

setwd("C:/Users/emorb/OneDrive - University of Cambridge/PhD/MR/Testosterone_CAD_MR/Testosterone CAD MR R files/TestosteroneCAD/")


#### with free t weights

t <- read.table("R scripts/MR HARMONISATION/RERUNNING MEDIATORS/male_testosterone_cluster.txt", header = TRUE)
cad_hits <- read.table("R scripts/MR HARMONISATION/RERUNNING MEDIATORS/male_testosterone_cluster_and_cad.txt", header = TRUE)

merged <- merge(t, cad_hits, by=c("CHR", "BP"))

## HARMONISATION

merged <- merged %>% 
  select(SNP.x, 
         T_effect=ALLELE1, 
         T_other=ALLELE0, 
         T_beta=BETA, 
         T_se=SE, 
         C_effect=reference_allele, 
         C_other=other_allele, 
         C_beta=male_beta, 
         C_se=male_se, 
         C_P=male_p_value)



merged <- merged %>% 
  mutate(T_abs_beta=abs(T_beta))

merged <- merged %>% 
  mutate(T_inc_allele=ifelse(T_beta<0, T_other, T_effect))

merged <- merged %>% 
  mutate(C_beta_harmonised=ifelse(C_effect!=T_inc_allele, C_beta*-1, C_beta))

merged$C_se <- as.numeric(merged$C_se)
IVW_weights <- merged$C_se^-2
IVW <- lm(C_beta_harmonised ~ T_abs_beta -1, weights = IVW_weights, data=merged)
summary(IVW)


plot(merged$T_abs_beta, merged$C_beta_harmonised)
abline(IVW, col = "red")


MRObject = mr_input(bx = merged$T_abs_beta, bxse = merged$T_se, 
                    by = merged$C_beta_harmonised, byse = merged$C_se, snps = merged$SNP.x)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)


png("leave_one_out_plot.png", width = 3000, height = 6000, res = 300)
mr_loo(MRObject)
dev.off()

###### MALE TESTOSTERONE ####################################################


# setting up the basic plot ########################################


plot(merged$T_abs_beta, merged$C_beta_harmonised, pch = 16, cex = 0.7,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male testosterone effect on CAD (including outlier)")

# Add error bars
segments(
  x0 = merged$T_abs_beta,
  y0 = merged$C_beta_harmonised - merged$C_se,
  x1 = merged$T_abs_beta,
  y1 = merged$C_beta_harmonised + merged$C_se,
  col = "black"
)

segments(
  x0 = merged$T_abs_beta - merged$T_se, 
  y0 = merged$C_beta_harmonised,
  x1 = merged$T_abs_beta + merged$T_se,
  y1 = merged$C_beta_harmonised,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

merged$C_se <- as.numeric(merged$C_se)
IVW_weights <- merged$C_se^-2 
inverse_weighted_LR <- lm(merged$C_beta_harmonised ~ merged$T_abs_beta- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = 0.003, b = 0.011, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("bottomleft", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))







library(devtools)


# devtools::install_github("rondolab/MR-PRESSO", force = TRUE)
library(MRPRESSO)

# run the M_TESTSOTERONE_CAD script before running this 

merged <- as.data.frame(merged)
mr_presso(BetaOutcome = "C_beta_harmonised", BetaExposure = "T_abs_beta", SdOutcome = "C_se", SdExposure = "T_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = merged, NbDistribution = 3500,  SignifThreshold = 0.05)



### removing outliers and performing analysis again 





#### with free t weights

t <- read.table("R scripts/MR HARMONISATION/RERUNNING MEDIATORS/male_testosterone_cluster.txt", header = TRUE)
cad_hits <- read.table("R scripts/MR HARMONISATION/RERUNNING MEDIATORS/male_testosterone_cluster_and_cad.txt", header = TRUE)

cad_hits %>%  filter(rsid_ukb=="rs56196860")


merged <- merge(t, cad_hits, by=c("CHR", "BP"))



## HARMONISATION

merged <- merged %>% 
  select(SNP.x, 
         T_effect=ALLELE1, 
         T_other=ALLELE0, 
         T_beta=BETA, 
         T_se=SE, 
         C_effect=reference_allele, 
         C_other=other_allele, 
         C_beta=male_beta, 
         C_se=male_se, 
         C_P=male_p_value)

merged <- merged[!merged$SNP.x == "rs56196860", ]

write_tsv(merged, "merged.tsv")

merged <- merged %>% 
  mutate(T_abs_beta=abs(T_beta))

merged <- merged %>% 
  mutate(T_inc_allele=ifelse(T_beta<0, T_other, T_effect))

merged <- merged %>% 
  mutate(C_beta_harmonised=ifelse(C_effect!=T_inc_allele, C_beta*-1, C_beta))

merged$C_se <- as.numeric(merged$C_se)
IVW_weights <- merged$C_se^-2
IVW <- lm(C_beta_harmonised ~ T_abs_beta -1, weights = IVW_weights, data=merged)
summary(IVW)


plot(merged$T_abs_beta, merged$C_beta_harmonised)
abline(IVW, col = "red")


MRObject = mr_input(bx = merged$T_abs_beta, bxse = merged$T_se, 
                    by = merged$C_beta_harmonised, byse = merged$C_se, snps = merged$SNP.x)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)


###### MALE TESTOSTERONE ####################################################


# setting up the basic plot ########################################


plot(merged$T_abs_beta, merged$C_beta_harmonised, pch = 16, cex = 0.7,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male testosterone effect on CAD (excluding outlier)")

# Add error bars
segments(
  x0 = merged$T_abs_beta,
  y0 = merged$C_beta_harmonised - merged$C_se,
  x1 = merged$T_abs_beta,
  y1 = merged$C_beta_harmonised + merged$C_se,
  col = "black"
)

segments(
  x0 = merged$T_abs_beta - merged$T_se, 
  y0 = merged$C_beta_harmonised,
  x1 = merged$T_abs_beta + merged$T_se,
  y1 = merged$C_beta_harmonised,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

merged$C_se <- as.numeric(merged$C_se)
IVW_weights <- merged$C_se^-2 
inverse_weighted_LR <- lm(merged$C_beta_harmonised ~ merged$T_abs_beta- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = -0.002, b = 0.218, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("bottomright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))






#### RUNNING IN WOMEN



t <- read.table("R scripts/MR HARMONISATION/RERUNNING MEDIATORS/female_testosterone_cluster.txt", header = TRUE)
cad_hits <- read.table("R scripts/MR HARMONISATION/RERUNNING MEDIATORS/female_testosterone_cluster_cad_info.txt", header = TRUE)




merged <- merge(t, cad_hits, by=c("CHR", "BP"))



## HARMONISATION

merged <- merged %>% 
  select(SNP.x, 
         T_effect=ALLELE1, 
         T_other=ALLELE0, 
         T_beta=BETA, 
         T_se=SE, 
         C_effect=reference_allele, 
         C_other=other_allele, 
         C_beta=female_beta, 
         C_se=female_se, 
         C_P=male_p_value)

merged <- merged[!merged$SNP == "rs56196860", ]



merged <- merged %>% 
  mutate(T_abs_beta=abs(T_beta))

merged <- merged %>% 
  mutate(T_inc_allele=ifelse(T_beta<0, T_other, T_effect))

merged <- merged %>% 
  mutate(C_beta_harmonised=ifelse(C_effect!=T_inc_allele, C_beta*-1, C_beta))

merged$C_se <- as.numeric(merged$C_se)
IVW_weights <- merged$C_se^-2
IVW <- lm(C_beta_harmonised ~ T_abs_beta -1, weights = IVW_weights, data=merged)
summary(IVW)


plot(merged$T_abs_beta, merged$C_beta_harmonised)
abline(IVW, col = "red")


MRObject = mr_input(bx = merged$T_abs_beta, bxse = merged$T_se, 
                    by = merged$C_beta_harmonised, byse = merged$C_se, snps = merged$SNP.x)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


mr_plot(MRObject, interactive=FALSE, labels=TRUE)
mr_forest(MRObject, ordered=TRUE)
mr_loo(MRObject)
mr_funnel(MRObject)

merged <- as.data.frame(merged)
mr_presso(BetaOutcome = "C_beta_harmonised", BetaExposure = "T_abs_beta", SdOutcome = "C_se", SdExposure = "T_se", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = merged, NbDistribution = 3500,  SignifThreshold = 0.05)



plot(merged$T_abs_beta, merged$C_beta_harmonised, pch = 16, cex = 0.7,
     xlab = "SNP effect on Testosterone",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female testosterone effect on CAD")

# Add error bars
segments(
  x0 = merged$T_abs_beta,
  y0 = merged$C_beta_harmonised - merged$C_se,
  x1 = merged$T_abs_beta,
  y1 = merged$C_beta_harmonised + merged$C_se,
  col = "black"
)

segments(
  x0 = merged$T_abs_beta - merged$T_se, 
  y0 = merged$C_beta_harmonised,
  x1 = merged$T_abs_beta + merged$T_se,
  y1 = merged$C_beta_harmonised,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

merged$C_se <- as.numeric(merged$C_se)
IVW_weights <- merged$C_se^-2 
inverse_weighted_LR <- lm(merged$C_beta_harmonised ~ merged$T_abs_beta- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = -0.002, b = 0.058, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("bottomright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))


#####################################    3. Mendelian Randomisation of SHBG to Cardiovascular Disease risk in Men            ########################################


library(tidyverse)
library(readxl)
library(MendelianRandomization)
library(openxlsx)


################################################################################

# HARMONISATION AND MR

##################################################################################

# looking at the allele matching and frequencies etc.
M_SHBG_proxies_output <- read_excel("not found inputs/SNPs_M_SHBG_AND_CAD.xlsx", sheet = "T&P E&O")
M_SHBG_proxies_output <- M_SHBG_proxies_output[-1,]


allele_matching <- select(M_SHBG_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "male_beta", "male_se" )

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_SHBG = "SNP",
    ALLELE1_SHBG = "ALLELE1",
    ALLELE0_SHBG = "ALLELE0",
    A1FREQ_SHBG = "A1FREQ",
    BETA_SHBG = "BETA",
    SE_SHBG = "SE",
    Effect_allele_CAD = "reference_allele",
    Other_allele_CAD = "other_allele",
    eaf_CAD = "eaf",
    male_beta_CAD = "male_beta",
    male_se_CAD = "male_se"
  )

# identify trait increasing allele for SHBG

allele_matching$SHBG_inc_allele <- if_else(allele_matching$BETA_SHBG<0, allele_matching$ALLELE0_SHBG, 
                                             allele_matching$ALLELE1_SHBG)

allele_matching$BETA_SHBG <- as.numeric(allele_matching$BETA_SHBG)
allele_matching$ABS_BETA_SHBG <- abs(allele_matching$BETA_SHBG)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$male_beta_CAD <- as.numeric(allele_matching$male_beta_CAD)

allele_matching$HARM_MALE_BETA_CAD <- if_else(allele_matching$SHBG_inc_allele!=allele_matching$Effect_allele_CAD,
                                       allele_matching$male_beta_CAD*-1, allele_matching$male_beta_CAD)




plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD)
M_SHBG_proxies_output$male_se <- as.numeric(M_SHBG_proxies_output$male_se)
IVW_weights <- M_SHBG_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model






M_SHBG_proxies_output$male_beta <- as.numeric(M_SHBG_proxies_output$male_beta)

allele_matching$ABS_BETA_SHBG <- as.numeric(allele_matching$ABS_BETA_SHBG)
allele_matching$male_se_CAD <- as.numeric(allele_matching$male_se_CAD)
allele_matching$SE_SHBG <- as.numeric(allele_matching$SE_SHBG)

MRObject = mr_input(bx = allele_matching$ABS_BETA_SHBG, bxse = allele_matching$SE_SHBG, 
                    by = allele_matching$HARM_MALE_BETA_CAD, byse = allele_matching$male_se_CAD)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

mr_plot(MRObject, orientate=TRUE, line="ivw")
mr_plot(mr_allmethods(MRObject, method="ivw"))


mr_allmethods(MRObject)
mr_plot(mr_allmethods(MRObject))

plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male SHBG")  # Replace with your desired y-axis label

M_SHBG_proxies_output$male_se <- as.numeric(M_SHBG_proxies_output$male_se)
IVW_weights <- M_SHBG_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")


# looking for SHBG SNP
SHBG_SNPS <- allele_matching %>%
  filter(SNP_SHBG=="rs1799941")

snp_x <- 0.12
snp_y <- 0.0074
snp_label <- "rs1799941"



plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_MALE_BETA_CAD,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Male SHBG")  # Replace with your desired y-axis label

M_SHBG_proxies_output$male_se <- as.numeric(M_SHBG_proxies_output$male_se)
IVW_weights <- M_SHBG_proxies_output$male_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_MALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")


text(snp_x, snp_y, snp_label, col = "blue", pos = 1, cex = 0.7)







#####################################    4. Mendelian Randomisation of SHBG to Cardiovascular Disease risk in Women            ########################################



################################################################################

# HARMONISATION AND MR

##################################################################################

library(tidyverse)
library(readxl)
library(MendelianRandomization)

# looking at the allele matching and frequencies etc.
F_SHBG_proxies_output <- read_excel("not found inputs/SNPs_F_SHBG_AND_CAD.xlsx", sheet = "T&P E&O")
F_SHBG_proxies_output <- F_SHBG_proxies_output[-1,]


allele_matching <- select(F_SHBG_proxies_output, "SNP", "ALLELE1", "ALLELE0", "A1FREQ", "BETA", "SE", "reference_allele", "other_allele", "eaf", "female_beta", "female_se" )


### removing the outlying SNP

allele_matching <- allele_matching[!allele_matching$SNP == "rs56196860", ]

# renaming the columns for ease of use 
allele_matching <- allele_matching %>%
  rename(
    SNP_SHBG = "SNP",
    ALLELE1_SHBG = "ALLELE1",
    ALLELE0_SHBG = "ALLELE0",
    A1FREQ_SHBG = "A1FREQ",
    BETA_SHBG = "BETA",
    SE_SHBG = "SE",
    effect_allele_CAD = "reference_allele",
    other_allele_CAD = "other_allele",
    eaf_CAD = "eaf",
    female_beta_CAD = "female_beta",
    female_se_CAD = "female_se"
  )

# identify trait increasing allele for SHBG

allele_matching$SHBG_inc_allele <- if_else(allele_matching$BETA_SHBG<0, allele_matching$ALLELE0_SHBG, 
                                           allele_matching$ALLELE1_SHBG)

allele_matching$BETA_SHBG <- as.numeric(allele_matching$BETA_SHBG)
allele_matching$ABS_BETA_SHBG <- abs(allele_matching$BETA_SHBG)

# harmonising so the effect alleles for CAD and SHBG are the same
# changing the betas here 

allele_matching$female_beta_CAD <- as.numeric(allele_matching$female_beta_CAD)

allele_matching$HARM_FEMALE_BETA_CAD <- if_else(allele_matching$SHBG_inc_allele!=allele_matching$effect_allele_CAD,
                                              allele_matching$female_beta_CAD*-1, allele_matching$female_beta_CAD)




plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_FEMALE_BETA_CAD)
F_SHBG_proxies_output$female_se <- as.numeric(F_SHBG_proxies_output$female_se)
IVW_weights <- F_SHBG_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red")
summary_model <- summary(inverse_weighted_LR)
summary_model



F_SHBG_proxies_output$female_beta <- as.numeric(F_SHBG_proxies_output$female_beta)

allele_matching$ABS_BETA_SHBG <- as.numeric(allele_matching$ABS_BETA_SHBG)
allele_matching$female_se_CAD <- as.numeric(allele_matching$female_se_CAD)
allele_matching$SE_SHBG <- as.numeric(allele_matching$SE_SHBG)

MRObject = mr_input(bx = allele_matching$ABS_BETA_SHBG, bxse = allele_matching$SE_SHBG, 
                    by = allele_matching$HARM_FEMALE_BETA_CAD, byse = allele_matching$female_se_CAD, snps = allele_matching$SNP_SHBG)

mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


mr_allmethods(MRObject)
mr_plot(MRObject, error = TRUE, line = "allmethods", interactive = FALSE, labels=FALSE)

mr_plot(mr_allmethods(MRObject))



allele_matching %>%
  filter(SNP_SHBG=="rs1799941")


plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_FEMALE_BETA_CAD,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female SHBG")  # Replace with your desired y-axis label

F_SHBG_proxies_output$female_se <- as.numeric(F_SHBG_proxies_output$female_se)
IVW_weights <- F_SHBG_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG - 1, weights = IVW_weights)
summary(inverse_weighted_LR)

abline(inverse_weighted_LR, col = "red")





plot(allele_matching$ABS_BETA_SHBG, allele_matching$HARM_FEMALE_BETA_CAD, pch = 16, cex = 0.7,
     xlab = "SNP effect on SHBG",  # Replace with your desired x-axis label
     ylab = "SNP effect on CAD",
     main = "Female SHBG")

# Add error bars
segments(
  x0 = allele_matching$ABS_BETA_SHBG,
  y0 = allele_matching$HARM_FEMALE_BETA_CAD - allele_matching$female_se_CAD,
  x1 = allele_matching$ABS_BETA_SHBG,
  y1 = allele_matching$HARM_FEMALE_BETA_CAD + allele_matching$female_se_CAD,
  col = "black"
)

segments(
  x0 = allele_matching$ABS_BETA_SHBG - allele_matching$SE_SHBG, 
  y0 = allele_matching$HARM_FEMALE_BETA_CAD,
  x1 = allele_matching$ABS_BETA_SHBG + allele_matching$SE_SHBG, 
  y1 = allele_matching$HARM_FEMALE_BETA_CAD,
  col = "black"
)

# adding the lines of the different models ############################
# IVW

F_SHBG_proxies_output$female_se <- as.numeric(F_SHBG_proxies_output$female_se)
IVW_weights <- F_SHBG_proxies_output$female_se^-2 
inverse_weighted_LR <- lm(allele_matching$HARM_FEMALE_BETA_CAD ~ allele_matching$ABS_BETA_SHBG- 1 ,weights=IVW_weights)
summary(inverse_weighted_LR)
abline(inverse_weighted_LR, col="red", lwd=1.6)


# EGGER
abline(a = -0.004, b = -0.128, col = "blue", lty = 1, lwd=1.6)

# MEDIAN 

legend("topright", legend = c("IWV method", "MR-Egger method"),
       col = c("red", "blue"), lty = c(1, 1), lwd = c(1.6, 1.6))



snp_x <- 0.12
snp_y <- -0.01
snp_label <- "rs1799941"

text(snp_x, snp_y, snp_label, col = "purple", pos = 1, cex = 0.7)



#####################################    5. MR-PRESSO analysis of Testosterone to Cardiovascular Disease Risk                ########################################


library(devtools)
devtools::install_github("rondolab/MR-PRESSO", force = TRUE)
library(MRPRESSO)

# run the M_TESTSOTERONE_CAD script before running this 

allele_matching <- as.data.frame(allele_matching)
mr_presso(BetaOutcome = "HARM_MALE_BETA_CAD", BetaExposure = "ABS_BETA_T", SdOutcome = "male_se_CAD", SdExposure = "SE_T", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = allele_matching, NbDistribution = 3000,  SignifThreshold = 0.05)


################################### 6. MULTIVARIABLE MR FOR DISCOVERY OF POTENTIAL MEDIATORS OF TESTOSTERONE CAD RELATIONSHIP #############################################


#### MULTIVARIABLE MR ####


library(data.table)
library(tidyverse)
library(dplyr)
library(MendelianRandomization)
library(ggplot2)


## reading in reference file to get the chrpos info for the testosterone snps 

ref <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Katherine/sandbox/apps/funcanno_v1.0/ukb_mfi_AllChr_v3.txt")
ref <- ref %>% 
  rename(CHR=V1,
         POS=V2,
         ALLELE0=V3,
         ALLELE1=V4,
         ID=V5)





### reading in the testosterone SNPs to try and match them and pull out their chrpos info

testosterone_cluster_snps <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/M_Testosterone_signals.tsv")
testosterone_cluster_snps <- testosterone_cluster_snps %>% 
  rename(ID=Signal)



### merge them with the reference set 

merged <- left_join(testosterone_cluster_snps, ref, by="ID")



### they all matched except the ones on the x chromosome 

### reading in the whole file and checking the allele info 

testosterone_cluster <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/M_Testosterone_cluster.tsv")

testosterone_cluster <- testosterone_cluster %>% 
  rename(ID=Signal)

merged_2 <- left_join(testosterone_cluster, merged, by="ID")
merged_2 <- merged_2 %>% 
  select(ID, Trait_raising, Other_allele, Weight, SE_weight, CHR, POS,ALLELE0, ALLELE1)

merged_2 <- merged_2 %>% 
  mutate(CHRPOS=paste0(CHR,":",POS))


write_tsv(merged_2, "/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_cluster_with_chrpos.tsv")

### added in the chrpos info for the x chromosome snps 

testosterone_cluster <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_cluster_with_chrpos.tsv")

testosterone_cluster <- testosterone_cluster %>% 
  select(ID, Trait_raising, Other_allele, Weight, SE_weight, CHR, POS)


### need to pull out the weights from the freeT gwas instead 

freet <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/GWAS_outputs/FreeT_Male_AllChrs.txt")

freet <- freet %>% 
  rename(POS=BP)

freet$CHR <- as.character(freet$CHR)
testosterone_cluster$CHR <- as.character(testosterone_cluster$CHR)

testosterone_cluster_free_t_weights <- merge(testosterone_cluster, freet,  by=c("CHR", "POS"))

testosterone_cluster_free_t_weights <- testosterone_cluster_free_t_weights %>% 
  select(CHR, POS, ID, Trait_raising, Other_allele, Weight, SE_weight, SNP, ALLELE1, ALLELE0, BETA, SE)

testosterone_cluster_free_t_weights <- testosterone_cluster_free_t_weights %>% 
  mutate(CHRPOS=paste0(CHR,":",POS))

### now investigating the HDL and LDL cholesterol files 
### this is to check the snps are in the right format 

hdl <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/HDL_INV_EUR_HRC_1KGP3_others_MALE.meta.singlevar.results.gz")
ldl <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/LDL_INV_EUR_HRC_1KGP3_others_MALE.meta.singlevar.results.gz")


### good, we have chrpos info so we can match on that
### grepping out the files in linux based on the chrpos info


## reading back in the grepped info to match with the testosterone chrpos info and run the MR 

hdl_grepped <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_snps_in_HDL_info.tsv")
hdl_grepped <- hdl_grepped %>% 
  mutate(ALLELE1=ALT,
         ALLELE0=REF)

hdl_grepped <- hdl_grepped %>% 
  filter(rsID!="rs56196860")

hdl_grepped <- hdl_grepped %>% 
  rename(CHR=CHROM,
         POS=POS_b37)

testosterone_cluster$POS <- as.integer(testosterone_cluster$POS)
hdl_grepped$CHR <- as.character(hdl_grepped$CHR)
testosterone_cluster$CHR <- as.character(testosterone_cluster$CHR)


merged <- merge(testosterone_cluster_free_t_weights, hdl_grepped, by=c("CHR", "POS"))

duplicates <- merged[duplicated(merged$ID), ]
merged <- merged[!(merged$ID %in% merged$ID[duplicated(merged$ID)]), ]


merged <- merged %>% 
  rename(effect_allele_t = ALLELE1.x, 
         other_allele_t = ALLELE0.y,
         beta_t = BETA, 
         se_t = SE.x,
         effect_allele_hdl = ALT,
         other_allele_hdl = REF,
         beta_hdl = EFFECT_SIZE,
         se_hdl = SE.y)


merged_selected <- merged %>% 
  select(ID, effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_hdl, other_allele_hdl, beta_hdl,se_hdl)


merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(hdl_beta_harmonised=ifelse(effect_allele_hdl!=t_inc_allele, beta_hdl*-1, beta_hdl))

merged_selected$se_hdl <- as.numeric(merged_selected$se_hdl)
IVW_weights <- merged_selected$se_hdl^-2
IVW <- lm(hdl_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)


plot(merged_selected$t_abs_beta, merged_selected$hdl_beta_harmonised)
abline(IVW, col = "red")


plot(merged_selected$t_abs_beta, merged_selected$hdl_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "HDL-c",
     main = "MR of testosterone to HDL-c")
abline(IVW, col = "red")






MRObject = mr_input(bx = merged_selected$t_abs_beta, bxse = merged_selected$se_t, 
                    by = merged_selected$hdl_beta_harmonised, byse = merged_selected$se_hdl, snps = merged_selected$ID)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)



hdl_info <- merged_selected %>% 
  select(ID, t_abs_beta,t_inc_allele, hdl_beta_harmonised, se_hdl, se_t)




#### REPEATING FOR LDL

## reading back in the grepped info to match with the testosterone chrpos info and run the MR 

ldl_grepped <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_snps_in_LDL_info.tsv")
ldl_grepped <- ldl_grepped %>% 
  mutate(ALLELE1=ALT,
         ALLELE0=REF)

ldl_grepped <- ldl_grepped %>% 
  filter(rsID!="rs56196860")

ldl_grepped <- ldl_grepped %>% 
  rename(CHR=CHROM,
         POS=POS_b37)

ldl_grepped$CHR <- as.character(ldl_grepped$CHR)

merged <- merge(testosterone_cluster_free_t_weights, ldl_grepped, by=c("CHR", "POS"))

duplicates <- merged[duplicated(merged$ID), ]
merged <- merged[!(merged$ID %in% merged$ID[duplicated(merged$ID)]), ]

merged <- merged %>% 
  rename(effect_allele_t = ALLELE1.x, 
         other_allele_t = ALLELE0.y,
         beta_t = BETA, 
         se_t = SE.x,
         effect_allele_ldl = ALT,
         other_allele_ldl = REF,
         beta_ldl = EFFECT_SIZE,
         se_ldl = SE.y)

merged_selected <- merged %>% 
  select(ID, effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_ldl, other_allele_ldl, beta_ldl,se_ldl)

merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(ldl_beta_harmonised=ifelse(effect_allele_ldl!=t_inc_allele, beta_ldl*-1, beta_ldl))

merged_selected$se_ldl <- as.numeric(merged_selected$se_ldl)
IVW_weights <- merged_selected$se_ldl^-2
IVW <- lm(ldl_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)

plot(merged_selected$t_abs_beta, merged_selected$ldl_beta_harmonised)
abline(IVW, col = "red")


plot(merged_selected$t_abs_beta, merged_selected$ldl_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "LDL-c",
     main = "MR of testosterone to LDL-c")
abline(IVW, col = "red")

MRObject = mr_input(bx = merged_selected$t_abs_beta, bxse = merged_selected$se_t, 
                    by = merged_selected$ldl_beta_harmonised, byse = merged_selected$se_ldl, snps = merged_selected$ID)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

ldl_info <- merged_selected %>% 
  select(ID, t_abs_beta,t_inc_allele, ldl_beta_harmonised, se_ldl, se_t)





### repeating for triglycerides - bear in mind these are log transformed 



## reading back in the grepped info to match with the testosterone chrpos info and run the MR 

tg_grepped <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_in_tg.tsv")
tg_grepped <- tg_grepped %>% 
  mutate(ALLELE1=ALT,
         ALLELE0=REF)

tg_grepped <- tg_grepped %>% 
  filter(rsID!="rs56196860")

tg_grepped <- tg_grepped %>% 
  rename(CHR=CHROM,
         POS=POS_b37)

tg_grepped$CHR <- as.character(tg_grepped$CHR)

merged <- merge(testosterone_cluster_free_t_weights, tg_grepped, by=c("CHR", "POS"))

duplicates <- merged[duplicated(merged$ID), ]
merged <- merged[!(merged$ID %in% merged$ID[duplicated(merged$ID)]), ]


merged <- merged %>% 
  rename(effect_allele_t = ALLELE1.x, 
         other_allele_t = ALLELE0.y,
         beta_t = BETA, 
         se_t = SE.x,
         effect_allele_tg = ALT,
         other_allele_tg = REF,
         beta_tg = EFFECT_SIZE,
         se_tg = SE.y)


merged_selected <- merged %>% 
  select(ID, effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_tg, other_allele_tg, beta_tg, se_tg)


merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(tg_beta_harmonised=ifelse(effect_allele_tg!=t_inc_allele, beta_tg*-1, beta_tg))

merged_selected$se_tg <- as.numeric(merged_selected$se_tg)
IVW_weights <- merged_selected$se_tg^-2
IVW <- lm(tg_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)


plot(merged_selected$t_abs_beta, merged_selected$tg_beta_harmonised)
abline(IVW, col = "red")


plot(merged_selected$t_abs_beta, merged_selected$tg_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "Triglycerides",
     main = "MR of testosterone to Triglycerides")
abline(IVW, col = "red")






MRObject = mr_input(bx = merged_selected$t_abs_beta, bxse = merged_selected$se_t, 
                    by = merged_selected$tg_beta_harmonised, byse = merged_selected$se_tg, snps = merged_selected$ID)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)




#### it is actually blood pressure that is the likely mediator so we need to look at that
## also looked up on chromosome and position info 

dbp_grepped <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/dbp_found.tsv")
dbp_grepped <- dbp_grepped %>% 
  rename(CHR=chromosome,
         POS=base_pair_location)


dbp_grepped$CHR <- as.character(dbp_grepped$CHR)
dbp_grepped <- dbp_grepped %>% 
  filter(rs_id!="rs56196860")


merged <- merge(testosterone_cluster_free_t_weights, dbp_grepped, by=c("CHR", "POS"))

merged <- merged[grepl("[AGTC]", merged$effect_allele), ]

merged <- merged[!merged$ID == "rs56196860", ]
duplicates <- merged[duplicated(merged$ID), ]
merged <- merged[!(merged$ID %in% merged$ID[duplicated(merged$ID)]), ]

merged <- merged %>% 
  rename(effect_allele_t = ALLELE1, 
         other_allele_t = ALLELE0,
         beta_t = BETA, 
         se_t = SE,
         effect_allele_dbp = effect_allele,
         other_allele_dbp = other_allele,
         beta_dbp = beta,
         se_dbp = standard_error)


merged_selected <- merged %>% 
  select(ID, CHR, POS,  effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_dbp, other_allele_dbp, beta_dbp,se_dbp)


merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(dbp_beta_harmonised=ifelse(effect_allele_dbp!=t_inc_allele, beta_dbp*-1, beta_dbp))

merged_selected$se_dbp <- as.numeric(merged_selected$se_dbp)
IVW_weights <- merged_selected$se_dbp^-2
IVW <- lm(dbp_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)


plot(merged_selected$t_abs_beta, merged_selected$dbp_beta_harmonised)
abline(IVW, col = "red")


plot(merged_selected$t_abs_beta, merged_selected$dbp_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "DBP",
     main = "MR of testosterone to DBP")
abline(IVW, col = "red")


MRObject = mr_input(bx = merged_selected$t_abs_beta, bxse = merged_selected$se_t, 
                    by = merged_selected$dbp_beta_harmonised, byse = merged_selected$se_dbp, snps = merged_selected$ID)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)

dbp_info <- merged_selected %>% 
  select(ID, CHR,POS, t_abs_beta,t_inc_allele, dbp_beta_harmonised, se_dbp, se_t)

#### it is actually blood pressure that is the likely mediator so we need to look at that
## also looked up on chromosome and position info 

sbp_grepped <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/sbp_found.tsv")
sbp_grepped <- sbp_grepped %>% 
  rename(CHR=chromosome,
         POS=base_pair_location)

sbp_grepped <- sbp_grepped %>% 
  filter(rs_id!="rs56196860")

sbp_grepped$CHR <- as.character(sbp_grepped$CHR)

merged <- merge(testosterone_cluster_free_t_weights, sbp_grepped, by=c("CHR", "POS"))
merged <- merged[grepl("[AGTC]", merged$effect_allele), ]

merged <- merged[!merged$ID == "rs56196860", ]
duplicates <- merged[duplicated(merged$ID), ]
merged <- merged[!(merged$ID %in% merged$ID[duplicated(merged$ID)]), ]

merged <- merged %>% 
  rename(effect_allele_t = ALLELE1, 
         other_allele_t = ALLELE0,
         beta_t = BETA, 
         se_t = SE,
         effect_allele_sbp = effect_allele,
         other_allele_sbp = other_allele,
         beta_sbp = beta,
         se_sbp = standard_error)


merged_selected <- merged %>% 
  select(ID, effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_sbp, other_allele_sbp, beta_sbp,se_sbp)

merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(sbp_beta_harmonised=ifelse(effect_allele_sbp!=t_inc_allele, beta_sbp*-1, beta_sbp))

merged_selected$se_sbp <- as.numeric(merged_selected$se_sbp)
IVW_weights <- merged_selected$se_sbp^-2
IVW <- lm(sbp_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)

plot(merged_selected$t_abs_beta, merged_selected$sbp_beta_harmonised)
abline(IVW, col = "red")

plot(merged_selected$t_abs_beta, merged_selected$sbp_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "SBP",
     main = "MR of testosterone to SBP")
abline(IVW, col = "red")

MRObject = mr_input(bx = merged_selected$t_abs_beta, bxse = merged_selected$se_t, 
                    by = merged_selected$sbp_beta_harmonised, byse = merged_selected$se_sbp, snps = merged_selected$ID)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)



#### multivariable MR with DBP


### need to pull out the CAD info 

testosterone <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_cluster_with_chrpos_free_t_weights.tsv")
cad_grepped <- fread("/rfs/project/rfs-mpB3sSsgAn4/Studies/People/Emily/Testosterone_CAD_MR/Multivariable_MR/testosterone_snps_in_cad_with_proxies.tsv")


merged <- merge(testosterone, cad_grepped, by=c("CHR", "BP"))
merged <- merged[grepl("[AGTC]", merged$reference_allele), ]
merged <- merged %>% 
  filter(SNP.x!="rs56196860")


merged <- merged %>% 
  rename(effect_allele_t = ALLELE1, 
         other_allele_t = ALLELE0,
         beta_t = BETA, 
         se_t = SE,
         effect_allele_cad = reference_allele,
         other_allele_cad = other_allele,
         beta_cad = male_beta,
         se_cad = male_se)


merged_selected <- merged %>% 
  select(SNP.x, CHR, BP, effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_cad, other_allele_cad, beta_cad,se_cad)

merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(cad_beta_harmonised=ifelse(effect_allele_cad!=t_inc_allele, beta_cad*-1, beta_cad))

merged_selected$se_cad <- as.numeric(merged_selected$se_cad)
IVW_weights <- merged_selected$se_cad^-2
IVW <- lm(cad_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)


plot(merged_selected$t_abs_beta, merged_selected$cad_beta_harmonised)
abline(IVW, col = "red")


plot(merged_selected$t_abs_beta, merged_selected$cad_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "CAD",
     main = "MR of testosterone to DBP")
abline(IVW, col = "red")



MRObject = mr_input(bx = merged_selected$t_abs_beta, bxse = merged_selected$se_t, 
                    by = merged_selected$cad_beta_harmonised, byse = merged_selected$se_cad, snps = merged_selected$ID)
mr_ivw(MRObject)
mr_egger(MRObject)
mr_median(MRObject)


####  calculating the sample size that would be needed in a randomised controlled trial to generate this odds ratio 

### OR = 1.10


library(powerMediation)


# Run sample size calculation for continuous predictor
SSizeLogisticCon(0.073, 1.168, 0.05, 0.9)


### doing it again but removing the outlying SNP


merged <- merge(testosterone, cad_grepped, by=c("CHR", "BP"))

merged <- merged[grepl("[AGTC]", merged$reference_allele), ]


merged <- merged %>% 
  rename(effect_allele_t = ALLELE1, 
         other_allele_t = ALLELE0,
         beta_t = BETA, 
         se_t = SE,
         effect_allele_cad = reference_allele,
         other_allele_cad = other_allele,
         beta_cad = male_beta,
         se_cad = male_se)


merged <- merged[!merged$SNP.x == "rs56196860", ]


merged_selected <- merged %>% 
  select(SNP.x, CHR, BP, effect_allele_t, other_allele_t, beta_t, se_t, effect_allele_cad, other_allele_cad, beta_cad,se_cad)


merged_selected <- merged_selected %>% 
  mutate(t_abs_beta=abs(beta_t))

merged_selected <- merged_selected %>% 
  mutate(t_inc_allele=ifelse(beta_t<0, other_allele_t, effect_allele_t))

merged_selected <- merged_selected %>% 
  mutate(cad_beta_harmonised=ifelse(effect_allele_cad!=t_inc_allele, beta_cad*-1, beta_cad))

merged_selected$se_cad <- as.numeric(merged_selected$se_cad)
IVW_weights <- merged_selected$se_cad^-2
IVW <- lm(cad_beta_harmonised ~ t_abs_beta -1, weights = IVW_weights, data=merged_selected)
summary(IVW)


plot(merged_selected$t_abs_beta, merged_selected$cad_beta_harmonised)
abline(IVW, col = "red")


plot(merged_selected$t_abs_beta, merged_selected$cad_beta_harmonised,
     xlab = "Testosterone", 
     ylab = "CAD",
     main = "MR of testosterone to DBP")
abline(IVW, col = "red")


cad_info <- merged_selected %>% 
  select(ID=SNP.x, CHR, POS=BP, t_inc_allele, t_abs_beta,se_t,cad_beta_harmonised, se_cad)

#### 

### merging all this info
cad_info$CHR <- as.character(cad_info$CHR)

t_dbp_cad <- merge(cad_info, dbp_info, by=c("CHR", "POS"))
# Define the weights as 1 / variance
weights <- 1 / (t_dbp_cad$se_cad^2)

# Run Weighted Least Squares regression
mv_ivw <- lm(cad_beta_harmonised ~ t_abs_beta.x + dbp_beta_harmonised, weights = weights, data = t_dbp_cad)

# Display results
summary(mv_ivw)




### trying for sbp instead 

t_sbp_cad <- merge(sbp_info, cad_info, by="ID")

# Define the weights as 1 / variance
weights <- 1 / (t_sbp_cad$se_cad^2)

# Run Weighted Least Squares regression
mv_ivw <- lm(cad_beta_harmonised ~ t_abs_beta.x + sbp_beta_harmonised, weights = weights, data = t_sbp_cad)

# Display results
summary(mv_ivw)

t_sbp_dbp_cad <- merge(cad_info, sbp_info, by="ID", all.x = TRUE)
t_sbp_dbp_cad <- t_sbp_dbp_cad %>% 
  select(ID, t_inc_allele=t_inc_allele.x, t_abs_beta=t_abs_beta.x, se_t=se_t.x, cad_beta_harmonised, se_cad, sbp_beta_harmonised, se_sbp)

t_sbp_dbp_cad <- merge(t_sbp_dbp_cad, dbp_info, by="ID", all.x=TRUE)

# Define the weights as 1 / variance
weights <- 1 / (t_sbp_dbp_cad$se_cad^2)

# Run Weighted Least Squares regression
mv_ivw <- lm(cad_beta_harmonised ~ t_abs_beta.x + sbp_beta_harmonised + dbp_beta_harmonised, weights = weights, data = t_sbp_dbp_cad)

# Display results
summary(mv_ivw)








#### HDL and LDLD

cad_ldl <- merge(cad_info, ldl_info, by="ID")

# Define the weights as 1 / variance
weights <- 1 / (cad_ldl$se_cad^2)

# Run Weighted Least Squares regression
mv_ivw <- lm(cad_beta_harmonised ~ t_abs_beta.x + ldl_beta_harmonised, weights = weights, data = cad_ldl)

# Display results
summary(mv_ivw)

cad_hdl <- merge(cad_info, hdl_info, by="ID")

# Define the weights as 1 / variance
weights <- 1 / (cad_hdl$se_cad^2)

# Run Weighted Least Squares regression
mv_ivw <- lm(cad_beta_harmonised ~ t_abs_beta.x + hdl_beta_harmonised, weights = weights, data = cad_hdl)

# Display results
summary(mv_ivw)









#####################################    7. Phenotyping and Survival Analysis of Testosterone in UK Biobank                  ########################################
#####################################    This section has to be carried out in the RAP 
#####################################    as EIDs of participants can only remain in the RAP and are required 
#####################################    for phenotypic analysis 

#####################################    7. Phenotyping and Survival Analysis of Testosterone in UK Biobank                  ########################################



library(tidyverse)
library(dplyr)
library(survival)
library(lubridate)
library(stringr)
library(cli)
library(data.table)
library(ggplot2)

#############            STEP 1: DATA CLEANING              #####################  


# including the survivor variables like age of recruitment etc. 

surv <- read.csv("data_participant_surv_2.csv")

View(surv)

colnames(surv) <- c("IID", "T", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "DATEASSESSMENT", "LTF")





################################################################################
######################### DATING AND PHENOTYPING CAD PHENOTYPES  ###############
################################################################################

## read in data on ICD_10 codes 
## icd10_1 is the participants and the dates of their diagnosis of specific 
## types of conditions 
## this dates them by category of condition rather than the very specific 
## condition types that are listed in the ICD_10 codes themselves 
## icd10_2 is a row per diagnosis of each individual 
## so each individual may have as many rows as they have conditions
## these are described by their ICD_10 code


icd10_1 <- read.csv("male_icd10_1.csv")
icd10_2 <- read.csv("male_icd10_2_hesin_diag.csv")

colnames(icd10_2) <- c("dnx_hesin_diag_id", "eid", "diag_icd10")



# assign the patterns of characters that we want to search for in the 
# ICD_10 codes 
# all of the cardiovascular conditions begin with an I and the 
# subsequent letters describe the more specific subtypes 


patterns <- c("I21", "I22", "I23","I24.1", "I25.2")


# now we are filtering the long list of conditions which has a row 
# for each individual and just keeping the ones with the cardiovascular 
# ICD_10 patterns 


relevant_icd10 <- icd10_2 %>%
  filter(grepl(paste(patterns, collapse = "|"), diag_icd10))



# then merge with the file which has all the dates to make it possible 
# to date the cardiovascular conditions 
# the different columns from the ICD_10_1 file represent 
# the different types of conditions 
# within them is the date which the person was diagnosed with the condition 


icd10_dates <- merge(icd10_1, relevant_icd10, by = "eid")
icd10_dates$X131298.0.0 <- as.Date(icd10_dates$X131298.0.0)
icd10_dates$X131300.0.0 <- as.Date(icd10_dates$X131300.0.0)
icd10_dates$X131302.0.0 <- as.Date(icd10_dates$X131302.0.0)
icd10_dates$X131304.0.0 <- as.Date(icd10_dates$X131304.0.0)
icd10_dates$X131306.0.0 <- as.Date(icd10_dates$X131306.0.0)


# now it is possible to find the minimum value of the dates of the 
# different diagnoses of types of CAD conditions 
# so we can use the pmin function to select the earliest date
# this date is required as we want the earliest diagnosis of CAD 
# to be used as the date of the event in the survival analysis 


icd10_dates$earliest_cad_date <- pmin(
  icd10_dates$X131298.0.0,
  icd10_dates$X131300.0.0,
  icd10_dates$X131302.0.0,
  icd10_dates$X131304.0.0,
  icd10_dates$X131306.0.0,
  na.rm = TRUE
)

icd10_dates <- icd10_dates %>% 
  rename(IID=eid)

icd10_dates_unique <- icd10_dates %>%
  group_by(IID) %>%
  slice_min(order_by = earliest_cad_date, with_ties = FALSE) %>%
  ungroup()


icd10_dates <- icd10_dates_unique

# checking if there are any NAs to see if anything has gone wrong in 
# finding the minimum value 
# now for every CAD condition recorded there should be the date associated with 
# the earliest diagnosis 

sum(is.na(icd10_dates$earliest_cad_date))

# changing the name of the ID column so it is possible to merge with the next 
# table

names(icd10_dates)[names(icd10_dates) == "eid"] <- "IID"







################################################################################
############# phenotyping CAD cases ############################################
################################################################################


## phenotypes 1 contains the list of participants and their icd10 and icd9
## diagnoses 
## there are very few ICD9 diagnoses because they are an older form
## i believe ICD9 codes were just used in scotland
## phenotypes 2 has testosterone levels 
## also has all of the recorded operations for these individuals 
## and any self reported illness 

phenotypes <- read.csv("new_cad_1.csv")
phenotypes2 <- read.csv("new_cad_2.csv")


# here we are going to collapse the icd10 data so that there is 
# not one row per condition per individual but so they are all 
# in one row for that individual

phenotypes_collapsed <- phenotypes %>%
  group_by(eid) %>%
  summarize(
    diag_icd10 = paste(diag_icd10, collapse = ", "),
    diag_icd9 = paste(diag_icd9, collapse = "")
  ) %>%
  ungroup()


# then we are going to pull out any CAD related conditions based on their 
# icd10 code and create a new column called CAD_ICD10 which places a 1
# if any CAD conditions were present in their ICD10 list 
# and a 0 if there were not any 

phenotypes_collapsed <- phenotypes_collapsed %>%
  mutate(CAD_ICD10 = if_else(grepl("I21|I22|I23|I24.1|I25.2", phenotypes_collapsed$diag_icd10), 1, 0))


# now doing the same for ICD9 codes 
# the ones listed in this code are the way CAD is recorded in ICD9


phenotypes_collapsed <- phenotypes_collapsed %>%
  mutate("CAD_ICD9" = if_else(grepl("^(413|414|434|436)", phenotypes_collapsed$diag_icd9), 1, 0))


## now we have coded whether or not individuals have or do not have CAD 
## as defined by a long list of ICD10 and ICD9 codes 
## now we are going to merge this with the file which has data on operations 
## and self reported CAD 


phenotypes_all <- merge(phenotypes_collapsed, phenotypes2, by = "eid")

## and rename the columns so we can understand them 

names(phenotypes_all)[names(phenotypes_all) == "X30850.0.0"] <- "T"
names(phenotypes_all)[names(phenotypes_all) == "X41272.0.0"] <- "OPS"
names(phenotypes_all)[names(phenotypes_all) == "X20002.0.0"] <- "Self-report"
names(phenotypes_all)[names(phenotypes_all) == "eid"] <- "IID"



## adding the ICD_10 dates file onto the phenotypes_all file 

phenotypes_all <- merge(phenotypes_all, icd10_dates, by = "IID", all.x =TRUE)


## removing all the columns for the different types of CAD 
## as we do not need these anymore


phenotypes_all <- subset(phenotypes_all, select = -c(X131306.0.0, X131304.0.0, 
                                                     X131302.0.0, X131300.0.0, 
                                                     X131298.0.0))


table(phenotypes_all$CAD_ICD9)
####### LOCATING AND DATING OPERATIONS #########################################


## now we are going to pull out all of the operations which are associated 
## with CAD - these fall under this long list of codes 
## if people have had these operations, they get a 1, if not they get a 0


phenotypes_all <- phenotypes_all %>%
  mutate("CAD_OP" = if_else(grepl("K40.1|K40.2|K40.3|K40.4|K40.8|K40.9|K41.1|
  K41.2|K41.3|K41.4|K41.8|K41.9|K45.1|K45.2|K45.3|K45.4|K45.5|K49.1|K49.2|K49.8|K49.9|K50.2|K75.1|K75.2|K75.3|K75.4|K75.8|K75.9", phenotypes_all$OPS), 1, 0))


## now we are going to find out what the dates of these operations were
## this is using a similar method as we used for the ICD10 data


## this file has data on the operations of each individual and when these operations
## happened

op_dates <- read.csv("operations_hesin_oper.csv")
colnames(op_dates) <- c("dnx_hesin_oper_id", "eid", "oper4", "opdate")


patterns <- c("K40.1|K40.2|K40.3|K40.4|K40.8|K40.9|K41.1|
  K41.2|K41.3|K41.4|K41.8|K41.9|K45.1|K45.2|K45.3|K45.4|K45.5|K49.1|K49.2|K49.8|K49.9|K50.2|K75.1|K75.2|K75.3|K75.4|K75.8|K75.9")


## now filtering the operations data to keep only the operations associated 
## with CAD


relevant_ops <- op_dates %>%
  filter(grepl(paste(patterns, collapse = "|"), oper4))


## renaming the ID file so we can merge with the phenotypes file

names(relevant_ops)[names(relevant_ops) == "eid"] <- "IID"

# Convert opdate to Date format if it's not already
relevant_ops$opdate <- as.Date(relevant_ops$opdate)

# Keep only the earliest operation per individual
relevant_ops <- relevant_ops %>%
  group_by(IID) %>%
  slice_min(order_by = opdate, with_ties = FALSE) %>%
  ungroup()


## merging the operations dates to the phenotypes_all file

phenotypes_all <- merge(phenotypes_all, relevant_ops, by = "IID", all.x = TRUE)

## removing redundant columns 

phenotypes_all <- subset(phenotypes_all, select = -c(dnx_hesin_oper_id))

duplicated_iids <- phenotypes_all %>% 
  filter(duplicated(IID) | duplicated(IID, fromLast = TRUE))

# View the duplicated rows
print(duplicated_iids)



########### GATHERING ALL CASES OF CAD INTO A SINGLE BINARY PHENOTYPE ##########


## setting all the CAD binary outcomes as numerics 

phenotypes_all$CAD_ICD10 <- as.numeric(phenotypes_all$CAD_ICD10)
phenotypes_all$CAD_ICD9 <- as.numeric(phenotypes_all$CAD_ICD9)
phenotypes_all$CAD_OP <- as.numeric(phenotypes_all$CAD_OP)


## telling R that if there is a 1 in any of these 3 columns to put a
## 1 in our new CADBIN column 
## this CADBIN column has a 1 if CAD has been identified by either ICD10, ICD9
## or operations codes 


phenotypes_all$CADBIN <- as.numeric(rowSums(phenotypes_all[, c("CAD_ICD10", "CAD_OP", "CAD_ICD9")]) > 0)




## now we are merging the earliest CAD date column from the ICD_10 data 
## and the operation date for those who had the operation, and selecting the 
## first instance

phenotypes_all$earliest_cad_date_all <- pmin(phenotypes_all$earliest_cad_date,
                                             phenotypes_all$opdate, na.rm = TRUE)



## then checking if there is anyone that does not have cad and has a date suggesting 
## they have cad 

sum(any(phenotypes_all$CADBIN == 0 & !is.na(phenotypes_all$earliest_cad_date_all)))

## now checking if there is anyone that does have cad but does not have a date 

sum(any(phenotypes_all$CADBIN == 1 & is.na(phenotypes_all$earliest_cad_date_all)))

phenotypes_condensed <- phenotypes_all %>% select("IID", "T", "CADBIN", "earliest_cad_date_all")


# Check for duplicated IIDs
duplicated_iids <- phenotypes_condensed %>% 
  filter(duplicated(IID) | duplicated(IID, fromLast = TRUE))

# View the duplicated rows
print(duplicated_iids)

## now merging our CAD survivorship info with the surv file which we 
## loaded in first and contains all the relevant covariates 


surv <- merge(phenotypes_condensed, surv, by = "IID", all.x = TRUE)
names(surv)[names(surv) == "T.x"] <- "T"
surv <- surv %>% select(-"T.y")

colnames(surv)


########### CALCULATING SURVIVOR VARIABLES - DATES ETC. #########################

# if they have a date in their lost to follow up column, place a 1
# otherwise, leave as 0 - now we have a binary column which says whether 
# someone was lost to follow up or not

surv$LTFBIN <- ifelse(surv$LTF == "" ,0,1)

# adding censoring date as the current date 
surv$censdate <- Sys.Date()
surv$censdate[surv$LTFBIN == 1] <- surv$LTF[surv$LTFBIN == 1]

table(surv$LTFBIN)

# Extract month and day of the earliest recorded CAD instance 
surv$cadmonth <- month(surv$earliest_cad_date_all)
surv$cadday <- day(surv$earliest_cad_date_all)
surv$cadyear <- year(surv$earliest_cad_date_all)
surv$censyear <- year(surv$censdate)

surv_with_dates <- surv

# changing months to numbers for month of birth
surv$monthbirthnum <- as.integer(factor(surv$MONTHBIRTH, levels = month.name))

# creating a censoring variable = anyone who does not come up as a CAD case
surv$censored <- ifelse(surv$CADBIN == 0, 1, 0)

table(surv$censored)



# creating a year of recruitment variable - this will be inaccurate 
surv$year_of_recruitment <- year(surv$DATEASSESSMENT)
surv$year_of_CAD <- year(surv$earliest_cad_date_all)

# creating a time to event variable
surv$timetoCADyears <- surv$cadyear-surv$year_of_recruitment


# creating a time to censoring variable 
surv$timetoCENSORyears <- ifelse(surv$censored == 1, surv$censyear - surv$year_of_recruitment, NA)

# creating a general time to event variable 
surv$timetoEVENT <- ifelse(is.na(surv$timetoCENSOR), surv$timetoCAD, surv$timetoCENSOR)

# create a variable for testosterone deficiency 
surv$testosterone_deficiency <- ifelse(surv$T < 12, 1, 0)

# need complete cases for testosterone deficiency 
# surv <- surv[complete.cases(surv$T.x), ]

any(duplicated(surv$IID))
any(duplicated(surv$IID))
surv <- surv[!duplicated(surv$IID), ]
surv <- surv[!duplicated(surv$IID), ]

surv$timetoCAD2 <- difftime(surv$earliest_cad_date_all, surv$timetoEVENT, units = "weeks")


# select relevant columns 

surv_key_variables <- surv %>% select("IID", "T", "CADBIN",
                                      "timetoEVENT")


# writing this surv file out so it can be used in the RAP to run the actual
# models

write.csv(surv_key_variables, "CAD_SURV.csv", row.names = TRUE)



# reading in the surv file that we have just written out 
# and removing some redundant columns 

surv <- read.csv("CAD_SURV.csv")



################################################################################
################# ADDING IN COVARIATES #########################################
################################################################################



###### ALL COVARIATES READ IN  ################################################

medications <- read.csv("medications_participant.csv")
sociodemographics <- read.csv("sociodemographics_participant.csv")
other_illness <- read.csv("other_diseases_participant.csv")




###### DIABETES ################################################################

diabetes <- read.csv("diabetes3_participant.csv")


colnames(diabetes) <- c("IID", "SELFREPORT", "MEDICATION", "DOCTOR", "HBA1C", "ICD10", "ICD9")


#### type 1

diabetes <- diabetes %>%
  mutate(TYPE1DIAB = if_else(grepl("E10|O240", diabetes$ICD10) | 
                               grepl("1222", diabetes$SELFREPORT)|
                               grepl("25001|25011|25021|25031|25041|25051|25061|25071|25081|25091|25003|25013|25023|25033|25043|25053|25063|25073|25083|25093", diabetes$ICD9), 1, 0))

table(diabetes$TYPE1DIAB)


#### type 2

diabetes <- diabetes %>%
  mutate(TYPE2DIAB = if_else(
    grepl("E11|O241", ICD10) | 
      grepl("1223|1220", SELFREPORT) |
      grepl("25000|25010|25020|25030|25040|25050|25060|25070|25080|25090|25002|25012|25022|25032|25042|25052|25062|25072|25082|25092", ICD9) |
      grepl("1140868902|1140874646|1140874674|1140874718|1140874744|1140883066|1140884600|1141152590|1141157284|1141168660|1141171646|1141173882|1141189090", MEDICATION) |
      replace_na(HBA1C > 48, FALSE), 
    1, 0
  ))

table(diabetes$TYPE2DIAB)



###### OTHER ILLNESS ###########################################################


colnames(other_illness) <- c("IID", "ICD10", "ICD9", "SELFREPORT", "MEDICATION")


###### arthritis

other_illness <- other_illness %>%
  mutate(ARTHRITIS = if_else(grepl("M05|M06", other_illness$ICD10) | 
                               grepl("1464", other_illness$SELFREPORT)|
                               grepl("714", other_illness$ICD9), 1, 0))


##### afib 

other_illness <- other_illness %>%
  mutate(AFIB = if_else(grepl("I48", other_illness$ICD10) | 
                          grepl("1471|1483", other_illness$SELFREPORT)|
                          grepl("4273|4720", other_illness$ICD9), 1, 0))


##### chronic kidney disease

other_illness <- other_illness %>%
  mutate(KIDNEY_DISEASE = if_else(grepl("N183|N184|N185", other_illness$ICD10) | 
                                    grepl("1192|1519|1609", other_illness$SELFREPORT)|
                                    grepl("5853|5855|5810|5820|5900|V420|V451", other_illness$ICD9), 1, 0))




##### migraine

other_illness <- other_illness %>%
  mutate(MIGRAINE = if_else(grepl("G43|G440|N943", other_illness$ICD10) | 
                              grepl("1265", other_illness$SELFREPORT)|
                              grepl("346", other_illness$ICD9), 1, 0))


##### SLE 

other_illness <- other_illness %>%
  mutate(SLE = if_else(grepl("M32", other_illness$ICD10) | 
                         grepl("1381", other_illness$SELFREPORT)|
                         grepl("7100", other_illness$ICD9), 1, 0))



##### MENTAL ILLNESS

other_illness <- other_illness %>%
  mutate(MENTAL_ILLNESS = if_else(grepl("F03|F068|F09|F20|F22|F23|F259|F28|F29|F31|F39|F53|F333", other_illness$ICD10) | 
                                    grepl("1289|1291", other_illness$SELFREPORT)|
                                    grepl("295|298|296", other_illness$ICD9), 1, 0))



##### ED

other_illness <- other_illness %>%
  mutate(ED = if_else(grepl("N484", other_illness$ICD10) | 
                        grepl("1518", other_illness$SELFREPORT)|
                        grepl("60784", other_illness$ICD9)|
                        grepl("1141168936|1141168948|1141168944|1141168946|1140869100|1140883010", other_illness$MEDICATION), 1, 0))




##################### MEDICATIONS #############################################

medications <- read.csv("medications_participant.csv")

colnames(medications) <- c("IID", "MEDICATION", "BPMEDS")


medications <- medications %>%
  mutate(HYPERTENSION = if_else(grepl("1140860192|1140860292|1140860696|
                                      1140860728|1140860750|1140860806|1140860882|
                                      1140860904|1140861088|1140861190|1140861276|
                                      1140866072|1140866078|1140866090|1140866102|
                                      1140866108|1140866122|1140866138|1140866156|
                                      1140866162|1140866724|1140866738|1140868618|
                                      1140872568|1140874706|1140874744|1140875808|
                                      1140879758|1140879760|1140879762|1140879802|
                                      1140879806|1140879810|1140879818|1140879822|
                                      1140879826|1140879830|1140879834|1140879842|
                                      1140879866|1140884298|1140888552|1140888556|
                                      1140888560|1140888646|1140909706|1140910442|
                                      1140910614|1140916356|1140923272|1140923336|
                                      1140923404|1140923712|1140926778|1140928226|
                                      1141145660|1141146126|1141152998|1141153026|
                                      1141164276|1141165470|1141166006|1141169516|
                                      1141171336|1141180592|1141180772|1141180778|
                                      1141184722|1141193282|1141194794|1141194810", medications$MEDICATION) | 
                                  grepl("2", medications$BPMEDS), 1, 0))



##### corticosteroids 


medications <- medications %>%
  mutate(CORTICOSTEROIDS = if_else(grepl("1140874790|1140874816|
1140874896.00|1140874930|1140874976|1141145782|1141173346", medications$MEDICATION) , 1, 0))


table(medications$CORTICOSTEROIDS)



##### antipsychotics 


medications <- medications %>%
  mutate(ANTIPSYCHOTICS = if_else(grepl("1140867420|1140867444|1140927956|
                                        1140928916|1141152848|1141153490|
                                        1141169714|1141195974", medications$MEDICATION) , 1, 0))


table(medications$ANTIPSYCHOTICS)







# SOCIODEMOGRAPHIC VARIABLES 

sociodemographics <- read.csv("sociodemographics_participant.csv")


### smoking


# smoking - this one is more tricky 

sociodemographics$exsmoker <- ifelse(sociodemographics$Ever.smoked...Instance.0 == "1" & sociodemographics$Current.tobacco.smoking...Instance.0=="0" ,"2",NA)
sociodemographics$nonsmoker <- ifelse(sociodemographics$Ever.smoked...Instance.0 == "0" ,"1",NA)
sociodemographics$NUM_CIGS_DAILY <- as.numeric(sociodemographics$Number.of.cigarettes.currently.smoked.daily..current.cigarette.smokers....Instance.0)
sociodemographics$SmokingCategory <- ifelse(sociodemographics$NUM_CIGS_DAILY < 10 & sociodemographics$NUM_CIGS_DAILY > 0, "3",
                                            ifelse(sociodemographics$NUM_CIGS_DAILY >= 10 & sociodemographics$NUM_CIGS_DAILY < 20, "4", 
                                                   ifelse(sociodemographics$NUM_CIGS_DAILY >= 20, "5", NA)))



# Remove NA values and replace with empty strings
sociodemographics$exsmoker[is.na(sociodemographics$exsmoker)] <- ""
sociodemographics$nonsmoker[is.na(sociodemographics$nonsmoker)] <- ""
sociodemographics$SmokingCategory[is.na(sociodemographics$SmokingCategory)] <- ""

# Combine columns into UKBBSMOKING with no white space
sociodemographics$UKBBSMOKING <- paste0(sociodemographics$exsmoker, sociodemographics$nonsmoker, sociodemographics$SmokingCategory)


sociodemographics$UKBBSMOKING <- as.factor(sociodemographics$UKBBSMOKING)

sociodemographics %>%
  mutate(UKBBSMOKING = factor(UKBBSMOKING,
                              levels = c("","1","2","3","4","5")))



sociodemographics <- sociodemographics %>% select(-"nonsmoker", -"SmokingCategory", -"NUM_CIGS_DAILY")





colnames(sociodemographics) <- c("IID", "AGERECRUIT", "MONTHBIRTH", "YEARBIRTH", "DEPRIVATION", "LOST_TO_FOLLOW_UP", "ETHNICITY", 
                                 "BMI", "EVERSMOKED", "SMOKINGSTATUS", "CURRENTSMOKING", "CIGSDAILY", "SBP1", "DATEASSESSMENT", 
                                 "SBP2", "CHOLESTEROL", "HDL", "ILLFATH", "ILLMOTH", "ILLSIBS", "EXSMOKER", "UKBBSMOKING")





##### BLOOD PRESSURE 

invalid_columns <- which(names(sociodemographics) == "" | is.na(names(sociodemographics)))
print(invalid_columns)

names(sociodemographics)[invalid_columns] <- paste0("InvalidName", seq_along(invalid_columns))

sociodemographics <- sociodemographics %>%
  mutate(SBP = (SBP1 + SBP2) / 2)


sociodemographics <- sociodemographics %>%
  rowwise() %>%
  mutate(SBP_SD = sd(c(SBP1, SBP2))) %>%
  ungroup()







##### CHOLESTEROL HDL RATIO


sociodemographics <- sociodemographics %>%
  mutate(CHOLESTEROLTOHDL = CHOLESTEROL / HDL)






##### ethnicity 


sociodemographics <- sociodemographics %>%
  mutate(ETHNICITY_CATEGORY = case_when(
    ETHNICITY %in% c(1, 1001, 1002, 1003, -1, -3) ~ 1,  # White or not stated
    ETHNICITY == 3001 ~ 2,  # Indian
    ETHNICITY == 3002 ~ 3,  # Pakistani
    ETHNICITY == 3003 ~ 4,  # Bangladeshi
    ETHNICITY == 3004 ~ 5,  # Other Asian
    ETHNICITY == 4001 ~ 6,  # Black Caribbean
    ETHNICITY == 4002 ~ 7,  # Black African
    ETHNICITY == 5 ~ 8,  # Chinese
    ETHNICITY == 6 ~ 9,  # Other ethnic group
    TRUE ~ NA_real_  # Default case, if none of the conditions match
  ))



sociodemographics <- sociodemographics %>%
  mutate(ETHNICITY_CATEGORY = factor(ETHNICITY_CATEGORY, levels = 1:9, labels = c(
    "White or not stated",
    "Indian",
    "Pakistani",
    "Bangladeshi",
    "Other Asian",
    "Black Caribbean",
    "Black African",
    "Chinese",
    "Other ethnic group"
  )))




##### ill family member

sociodemographics <- sociodemographics %>%
  mutate(FAMHISTORY = if_else(
    grepl("\\b1\\b", ILLFATH) | grepl("\\b1\\b", ILLMOTH) | grepl("\\b1\\b", ILLSIBS),
    1, 0
  ))



#### selecting the important variables before merging 

medications <- medications %>% select("IID", "HYPERTENSION", "CORTICOSTEROIDS", "ANTIPSYCHOTICS")
diabetes <- diabetes %>% select("IID", "TYPE1DIAB", "TYPE2DIAB")
other_illness <- other_illness %>% select(-c("MEDICATION", "SELFREPORT"))


###### MERGING ALL OF THESE TABLES 

COVARIATES <- merge(diabetes, medications, by = "IID")
COVARIATES <- merge(COVARIATES, other_illness, by = "IID")
COVARIATES <- merge(COVARIATES, sociodemographics, by = "IID")


###### merging with the survival data

COMPLETE_DATA <- merge(surv, COVARIATES, by = "IID", all.x=TRUE)


COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("CIGSDAILY"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("EXSMOKER"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("ICD10", "ICD9"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("LOST_TO_FOLLOW_UP"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("SMOKINGSTATUS", "CURRENTSMOKING"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("X"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("ILLFATH", "ILLMOTH", "ILLSIBS"))
COMPLETE_DATA <- COMPLETE_DATA %>% select(-c("SBP1", "SBP2"))

any(duplicated(surv_with_dates$IID))
any(duplicated(COMPLETE_DATA$IID))
surv_with_dates <- surv_with_dates[!duplicated(surv_with_dates$IID), ]
COMPLETE_DATA <- COMPLETE_DATA[!duplicated(COMPLETE_DATA$IID), ]



REMOVE_CAD_PRIOR_TO_ASSESSMENT <- merge(surv_with_dates, COMPLETE_DATA, by = "IID", all.x=TRUE)


REMOVE_CAD_PRIOR_TO_ASSESSMENT$DATEASSESSMENT.x <- as.Date(REMOVE_CAD_PRIOR_TO_ASSESSMENT$DATEASSESSMENT.x)

table(REMOVE_CAD_PRIOR_TO_ASSESSMENT$earliest_cad_date_all <= REMOVE_CAD_PRIOR_TO_ASSESSMENT$DATEASSESSMENT.x)

table(REMOVE_CAD_PRIOR_TO_ASSESSMENT$CADBIN.x)



dates <- read.csv("attendance_participant.csv")

names(dates)[names(dates) == "Participant.ID"] <- "IID"

dat <- merge(REMOVE_CAD_PRIOR_TO_ASSESSMENT, dates, by = "IID")

names(dat)[names(dat) == "Date.of.attending.assessment.centre...Instance.0"] <- "ASSESSMENT_DATE"


cleaned_data <- dat %>%
  filter(is.na(earliest_cad_date_all) | is.na(ASSESSMENT_DATE) | earliest_cad_date_all >= ASSESSMENT_DATE)



cleaned_data <- cleaned_data %>% select(-c("MONTHBIRTH.x"))
cleaned_data <- cleaned_data %>% select(-c("YEARBIRTH.x"))
cleaned_data <- cleaned_data %>% select(-c("DATEASSESSMENT.x"))
cleaned_data <- cleaned_data %>% select(-c("AGERECRUIT.y"))
cleaned_data <- cleaned_data %>% select(-c("LTF"))
cleaned_data <- cleaned_data %>% select(-c("LTFBIN"))
cleaned_data <- cleaned_data %>% select(-c("censdate"))
cleaned_data <- cleaned_data %>% select(-c("cadmonth"))
cleaned_data <- cleaned_data %>% select(-c("cadday"))
cleaned_data <- cleaned_data %>% select(-c("cadyear"))
cleaned_data <- cleaned_data %>% select(-c("censyear"))
cleaned_data <- cleaned_data %>% select(-c("T.y"))
cleaned_data <- cleaned_data %>% select(-c("CADBIN.y"))
cleaned_data <- cleaned_data %>% select(-c("MONTHBIRTH.y"))
cleaned_data <- cleaned_data %>% select(-c("YEARBIRTH.y"))


table(cleaned_data$ASSESSMENT_DATE > cleaned_data$earliest_cad_date_all)


cleaned_data$got_CAD <- !is.na(cleaned_data$earliest_cad_date_all)
cleaned_data$date_comparison <- ifelse(cleaned_data$got_CAD, 
                                       ifelse(cleaned_data$earliest_cad_date_all < cleaned_data$ASSESSMENT_DATE, "Before", "After"),
                                       NA)


# Count the number of participants who never got CAD
num_na_CAD <- sum(is.na(cleaned_data$earliest_cad_date_all))
sum(cleaned_data$CADBIN.x==1)


table(cleaned_data$date_comparison, useNA = "ifany")



table(cleaned_data$CADBIN.x)


cleaned_data <- cleaned_data %>% select(-c("earliest_cad_date_all"))
cleaned_data <- cleaned_data %>% select(-c("date_comparison"))


complete_data <- cleaned_data[complete.cases(cleaned_data), ]

table(complete_data$CADBIN.x)

nrow(complete_data)

qrisksurv <- complete_data

names(qrisksurv)[names(qrisksurv) == "T.x"] <- "T"
names(qrisksurv)[names(qrisksurv) == "CADBIN.x"] <- "CADBIN"


# fixing the ethnicity labelling 



unique(complete_data$ETHNICITY_CATEGORY)

table(complete_data$ETHNICITY_CATEGORY)

complete_data$ETHNICITY_CATEGORY <- as.factor(complete_data$ETHNICITY_CATEGORY)


complete_data$ETHNICITY_CATEGORY <- relevel(complete_data$ETHNICITY_CATEGORY, ref = "White or not stated")



# creating quartiles of the testosterone distribution

quartiles <- quantile(qrisksurv$T, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
qrisksurv$testosterone_quartile <- cut(qrisksurv$T, breaks = quartiles,
                                       labels = c("lower", "lower middle", "upper middle", "upper"),
                                       include.lowest = TRUE)






# age categories 

qrisksurv$age_group <- cut(qrisksurv$AGERECRUIT.x,
                           breaks = c(40, 45, 50, 55, 60, 65, 70, Inf),
                           labels = c("40-45", "45-50", "50-55", "55-60", "60-65", "65-70", "70+"),
                           right = FALSE)


complete_data <- qrisksurv



### smoking variable


complete_data$UKBBSMOKING <- factor(complete_data$UKBBSMOKING, levels = c("1", "2", "3", "4", "5"))
complete_data$UKBBSMOKING <- relevel(complete_data$UKBBSMOKING, ref = "1")

# deficient, sufficient, and high groups 

# Define cutoff points and labels
cut_points <- c(-Inf, 12, 18, 25, 30)
labels <- c("deficient", "sufficient", "high", "very high")

# Create categorical variable for testosterone categories
complete_data$testosterone_category <- cut(complete_data$T, breaks = cut_points, labels = labels, include.lowest = TRUE)

complete_data$testosterone_category <- relevel(complete_data$testosterone_category, ref = "sufficient")

complete_data <- complete_data %>% 
  mutate(testosterone_binary=ifelse(T<12, "Deficient", "Sufficient"))


#### DESCRIPTIVE TABLE

# Install if not already installed

library(tableone)


complete_data <- complete_data %>% 
  rename(AGE_RECRUIT=AGERECRUIT.x)

# Define covariates and factors
vars <- c("AGE_RECRUIT", "CADBIN", "BMI", "TYPE2DIAB", "HYPERTENSION", "CORTICOSTEROIDS",
          "ANTIPSYCHOTICS", "ARTHRITIS", "AFIB", "KIDNEY_DISEASE", "MIGRAINE",
          "SLE", "MENTAL_ILLNESS", "ED", "DEPRIVATION", "CHOLESTEROL",
          "HDL", "SBP", "SBP_SD", "CHOLESTEROLTOHDL", "FAMHISTORY", "EVERSMOKED")

factorVars <- c("CADBIN", "TYPE2DIAB", "HYPERTENSION", "CORTICOSTEROIDS", "ANTIPSYCHOTICS", 
                "ARTHRITIS", "AFIB", "KIDNEY_DISEASE", "MIGRAINE", "SLE", 
                "MENTAL_ILLNESS", "ED", "FAMHISTORY", "EVERSMOKED")

# Create table
table1 <- CreateTableOne(vars = vars, strata = "testosterone_binary", 
                         data = complete_data, factorVars = factorVars)

# Convert to data frame
table_df <- as.data.frame(print(table1, printToggle = FALSE, noSpaces = TRUE))

# View or export
head(table_df)

table_df <- table_df %>% 
  select(Deficient, Sufficient, p)


table_df_out <- rownames_to_column(table_df, var = "Variable")

# Now write to .tsv file
write_tsv(table_df_out, "descriptive_phenotypic_info.tsv")







# SURVIVAL ANALYSIS 


complete_data_cases <- complete_data %>% 
  filter(got_CAD==1)
# making adjustments - just looking at deficient and sufficient 

complete_data$testosterone_binary <- ifelse(complete_data$T < 12, "deficient", "sufficient")

complete_data$testosterone_binary<- factor(complete_data$testosterone_binary, levels = c("deficient", "sufficient"))

# Reorder levels so that "sufficient" is the reference category
complete_data$testosterone_binary <- relevel(complete_data$testosterone_binary, ref = "sufficient")


# Create survival object
CADsurv <- Surv(time = complete_data$timetoEVENT, event = complete_data$CADBIN)

# Fit Kaplan-Meier survival curves
km_fit <- survfit(CADsurv ~ testosterone_binary, data = complete_data)



plot(
  km_fit,
  col = c("blue", "red"),           # Colors for the curves
  lty = 1:2,                        # Line types for the curves
  xlab = "Time to Event",
  ylab = "Survival Probability",
  main = "Kaplan-Meier Survival Curves"
)



# Function to calculate rates and CIs
calculate_rates <- function(complete_data, age_group) {
  group_data <- complete_data %>% filter(age_group == !!age_group)
  person_years <- sum(group_data$timetoEVENT)
  incident_cases <- sum(group_data$CADBIN)
  rate_per_1000 <- (incident_cases / person_years) * 1000
  se_rate <- sqrt(incident_cases) / person_years * 1000
  lower_ci <- rate_per_1000 - 1.96 * se_rate
  upper_ci <- rate_per_1000 + 1.96 * se_rate
  c(Rate_per_1000 = rate_per_1000, Lower_CI = lower_ci, Upper_CI = upper_ci, Person_Years = person_years, Incident_Cases = incident_cases)
}

# List of age groups
age_groups <- levels(complete_data$age_group)

# Calculate for each group
results <- data.frame()
for (age_group in age_groups) {
  rates <- calculate_rates(complete_data, age_group)
  results <- rbind(results, data.frame(Age_Group = age_group, t(rates)))
}

# View results
print(results)

write.csv(results, "cardiovascular_disease_rates.csv", row.names = FALSE)

table(complete_data$CADBIN)

cox <- coxph(CADsurv~testosterone_binary, data = complete_data)
summary(cox)

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$BMI)
summary(cox) 

cox <- coxph(CADsurv~complete_data$BMI)
summary(cox)

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$TYPE1DIAB)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$TYPE2DIAB)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$HYPERTENSION)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$CORTICOSTEROIDS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ANTIPSYCHOTICS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ARTHRITIS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$AFIB)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$KIDNEY_DISEASE)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$MIGRAINE)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$SLE)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$MENTAL_ILLNESS)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ED)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$AGERECRUIT)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$DEPRIVATION)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$ETHNICITY_CATEGORY)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$BMI)
summary(cox) 


cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$UKBBSMOKING)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$SBP)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$SBP_SD)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$CHOLESTEROLTOHDL)
summary(cox) 

cox <- coxph(CADsurv~complete_data$testosterone_binary+complete_data$FAMHISTORY)
summary(cox) 




# Define a function to extract Cox model results from the summary
extract_cox_summary <- function(model, model_name) {
  summary_model <- summary(model)
  coef <- summary_model$coefficients
  confint <- summary_model$conf.int
  results <- data.frame(
    model = model_name,
    hazard_ratio = coef[, "exp(coef)"],
    p.value = coef[, "Pr(>|z|)"],
    conf.low = confint[, "lower .95"],
    conf.high = confint[, "upper .95"]
  )
  return(results)
}

# Define a list of models with their names
models <- list(
  "testosterone_binary" = coxph(CADsurv ~ testosterone_binary, data = complete_data),
  "testosterone_binary + BMI" = coxph(CADsurv ~ testosterone_binary + BMI, data = complete_data),
  "testosterone_binary + TYPE1DIAB" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$TYPE1DIAB),
  "testosterone_binary + TYPE2DIAB" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$TYPE2DIAB),
  "testosterone_binary + HYPERTENSION" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$HYPERTENSION),
  "testosterone_binary + CORTICOSTEROIDS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$CORTICOSTEROIDS),
  "testosterone_binary + ANTIPSYCHOTICS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ANTIPSYCHOTICS),
  "testosterone_binary + ARTHRITIS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ARTHRITIS),
  "testosterone_binary + AFIB" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$AFIB),
  "testosterone_binary + KIDNEY_DISEASE" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$KIDNEY_DISEASE),
  "testosterone_binary + MIGRAINE" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$MIGRAINE),
  "testosterone_binary + SLE" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$SLE),
  "testosterone_binary + MENTAL_ILLNESS" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$MENTAL_ILLNESS),
  "testosterone_binary + ED" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ED),
  "testosterone_binary + AGERECRUIT" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$AGE_RECRUIT),
  "testosterone_binary + DEPRIVATION" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$DEPRIVATION),
  "testosterone_binary + ETHNICITY_CATEGORY" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$ETHNICITY_CATEGORY),
  "testosterone_binary + UKBBSMOKING" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$UKBBSMOKING),
  "testosterone_binary + SBP" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$SBP),
  "testosterone_binary + SBP_SD" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$SBP_SD),
  "testosterone_binary + CHOLESTEROLTOHDL" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$CHOLESTEROLTOHDL),
  "testosterone_binary + FAMHISTORY" = coxph(CADsurv ~ complete_data$testosterone_binary + complete_data$FAMHISTORY)
)

# Extract results for each model and combine them into a single table
combined_results <- bind_rows(lapply(names(models), function(model_name) {
  model <- models[[model_name]]
  extract_cox_summary(model, model_name)
}))

# Print the combined results
print(combined_results)

# Save the results to a CSV file (optional)
write.csv(combined_results, file = "separate_mediators.csv", row.names = TRUE)



##### ADJUSTING FOR EVERYTHING 

# Fit Cox proportional hazards model with multiple predictors
cox <- coxph(
  formula = CADsurv ~ testosterone_binary + BMI + TYPE1DIAB + TYPE2DIAB + HYPERTENSION +
    CORTICOSTEROIDS + ANTIPSYCHOTICS + ARTHRITIS + AFIB + KIDNEY_DISEASE +
    MIGRAINE + SLE + MENTAL_ILLNESS + ED + AGE_RECRUIT + DEPRIVATION +
    ETHNICITY_CATEGORY + UKBBSMOKING + SBP + SBP_SD + CHOLESTEROLTOHDL +
    FAMHISTORY,
  data = complete_data
)

# Summarize the model
summary(cox)



cox_summary <- summary(cox)


results <- data.frame(
  exp_coef = cox_summary$coefficients[, "exp(coef)"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  lower_ci = cox_summary$conf.int[, 1],  # Lower 95% CI
  upper_ci = cox_summary$conf.int[, 2]   # Upper 95% CI
)

write.csv(results, file = "cox_model_results.csv", row.names = TRUE)





