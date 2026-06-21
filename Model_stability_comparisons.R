############################################
###  Compare model repetitions: setup    ###
############################################


#clean workspace
rm(list = ls())

#load packages
library(Hmsc); library(coda)

#list WDs
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Data'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_4th_models'
wd_plots <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_4th_models/Plots'
wd_model_comparison <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_4th_models/Model_comparison'

#setwd
setwd(wd_models)

#check files in model folder
list.files(wd_models)



########################################
###  Load m2 repetitions             ###
########################################



#load m2 rep1 objects
m2_rep1 <- readRDS('m2_spatial_latent.rds')
mpost_m2_rep1 <- readRDS('mpost_m2_spatial_latent.rds')
postBeta_m2_rep1 <- readRDS('postBeta_m2_spatial_latent.rds')
Omega_m2_rep1 <- readRDS('Omega_m2_spatial_latent.rds')


#load m2 rep2 objects
m2_rep2 <- readRDS('m2_rep2_spatial_latent.rds')
mpost_m2_rep2 <- readRDS('mpost_m2_rep2_spatial_latent.rds')
postBeta_m2_rep2 <- readRDS('postBeta_m2_rep2_spatial_latent.rds')
Omega_m2_rep2 <- readRDS('Omega_m2_rep2_spatial_latent.rds')


#load m2 rep3 objects
m2_rep3 <- readRDS('m2_rep3_spatial_latent.rds')
mpost_m2_rep3 <- readRDS('mpost_m2_rep3_spatial_latent.rds')
postBeta_m2_rep3 <- readRDS('postBeta_m2_rep3_spatial_latent.rds')
Omega_m2_rep3 <- readRDS('Omega_m2_rep3_spatial_latent.rds')



########################################
###  Extract m2 comparison matrices  ###
########################################



#extract beta mean matrices
Beta_mean_m2_rep1 <- postBeta_m2_rep1$mean
Beta_mean_m2_rep2 <- postBeta_m2_rep2$mean
Beta_mean_m2_rep3 <- postBeta_m2_rep3$mean


#extract Omega mean matrices
Omega_mean_m2_rep1 <- Omega_m2_rep1[[1]]$mean
Omega_mean_m2_rep2 <- Omega_m2_rep2[[1]]$mean
Omega_mean_m2_rep3 <- Omega_m2_rep3[[1]]$mean


#extract Omega support matrices
Omega_support_m2_rep1 <- Omega_m2_rep1[[1]]$support
Omega_support_m2_rep2 <- Omega_m2_rep2[[1]]$support
Omega_support_m2_rep3 <- Omega_m2_rep3[[1]]$support



########################################
###  Compare m2 beta means            ###
########################################



#compare beta posterior means across repetitions
cor_beta_m2_rep1_rep2 <- cor(as.vector(Beta_mean_m2_rep1),
                             as.vector(Beta_mean_m2_rep2))

cor_beta_m2_rep1_rep3 <- cor(as.vector(Beta_mean_m2_rep1),
                             as.vector(Beta_mean_m2_rep3))

cor_beta_m2_rep2_rep3 <- cor(as.vector(Beta_mean_m2_rep2),
                             as.vector(Beta_mean_m2_rep3))


#combine results
beta_m2_correlations <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  correlation = c(cor_beta_m2_rep1_rep2,
                  cor_beta_m2_rep1_rep3,
                  cor_beta_m2_rep2_rep3)
)

#inspect results
beta_m2_correlations



########################################
###  Beta mean differences            ###
########################################



mad_beta_m2_rep1_rep2 <- mean(abs(Beta_mean_m2_rep1 -
                                    Beta_mean_m2_rep2))

mad_beta_m2_rep1_rep3 <- mean(abs(Beta_mean_m2_rep1 -
                                    Beta_mean_m2_rep3))

mad_beta_m2_rep2_rep3 <- mean(abs(Beta_mean_m2_rep2 -
                                    Beta_mean_m2_rep3))


beta_m2_differences <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  mean_absolute_difference =
    c(mad_beta_m2_rep1_rep2,
      mad_beta_m2_rep1_rep3,
      mad_beta_m2_rep2_rep3)
)

beta_m2_differences


# Beta posterior means are highly reproducible across independent
# MCMC runs (r = 0.98–0.99).
# 
# Mean absolute differences are small
# (0.06–0.11 beta units),
# suggesting environmental response estimates are robust to
# MCMC stochasticity.



################################################
###  Compare Omega means: species pairs only ###
################################################


#extract lower triangle values only
#
#this removes:
#
#1. the diagonal
#2. the duplicated upper triangle
#
#so we compare only real species-pair associations

Omega_mean_m2_rep1_pairs <- Omega_mean_m2_rep1[lower.tri(Omega_mean_m2_rep1)]
Omega_mean_m2_rep2_pairs <- Omega_mean_m2_rep2[lower.tri(Omega_mean_m2_rep2)]
Omega_mean_m2_rep3_pairs <- Omega_mean_m2_rep3[lower.tri(Omega_mean_m2_rep3)]


#compare Omega posterior means across repetitions
cor_omega_pairs_m2_rep1_rep2 <- cor(Omega_mean_m2_rep1_pairs,
                                    Omega_mean_m2_rep2_pairs)

cor_omega_pairs_m2_rep1_rep3 <- cor(Omega_mean_m2_rep1_pairs,
                                    Omega_mean_m2_rep3_pairs)

cor_omega_pairs_m2_rep2_rep3 <- cor(Omega_mean_m2_rep2_pairs,
                                    Omega_mean_m2_rep3_pairs)


#combine results
omega_pairs_m2_correlations <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  correlation = c(cor_omega_pairs_m2_rep1_rep2,
                  cor_omega_pairs_m2_rep1_rep3,
                  cor_omega_pairs_m2_rep2_rep3)
)

#inspect results
omega_pairs_m2_correlations



###########################################
###  Omega mean absolute differences    ###
###########################################


mad_omega_m2_rep1_rep2 <- mean(abs(Omega_mean_m2_rep1_pairs -
                                     Omega_mean_m2_rep2_pairs))

mad_omega_m2_rep1_rep3 <- mean(abs(Omega_mean_m2_rep1_pairs -
                                     Omega_mean_m2_rep3_pairs))

mad_omega_m2_rep2_rep3 <- mean(abs(Omega_mean_m2_rep2_pairs -
                                     Omega_mean_m2_rep3_pairs))


omega_m2_differences <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  mean_absolute_difference =
    c(mad_omega_m2_rep1_rep2,
      mad_omega_m2_rep1_rep3,
      mad_omega_m2_rep2_rep3)
)

#inspect results
omega_m2_differences



#################################################
###  Omega sign agreement: rep1 vs rep2       ###
#################################################



#extract signs of residual associations
sign_rep1 <- sign(Omega_mean_m2_rep1_pairs)
sign_rep2 <- sign(Omega_mean_m2_rep2_pairs)

#calculate agreement
same_sign_rep1_rep2 <- mean(sign_rep1 == sign_rep2)

same_sign_rep1_rep2



#################################################
###  Distribution of Omega magnitudes         ###
#################################################



#proportion of species pairs with moderate residual associations

mean(abs(Omega_mean_m2_rep1_pairs) > 0.1)
mean(abs(Omega_mean_m2_rep1_pairs) > 0.2)
mean(abs(Omega_mean_m2_rep1_pairs) > 0.3)

mean(abs(Omega_mean_m2_rep2_pairs) > 0.1)
mean(abs(Omega_mean_m2_rep2_pairs) > 0.2)
mean(abs(Omega_mean_m2_rep2_pairs) > 0.3)

mean(abs(Omega_mean_m2_rep3_pairs) > 0.1)
mean(abs(Omega_mean_m2_rep3_pairs) > 0.2)
mean(abs(Omega_mean_m2_rep3_pairs) > 0.3)



#################################
###  Omega sign agreement     ###
#################################


sign_rep1 <- sign(Omega_mean_m2_rep1_pairs)
sign_rep2 <- sign(Omega_mean_m2_rep2_pairs)
sign_rep3 <- sign(Omega_mean_m2_rep3_pairs)

same_sign_rep1_rep2 <- mean(sign_rep1 == sign_rep2)

same_sign_rep1_rep3 <- mean(sign_rep1 == sign_rep3)

same_sign_rep2_rep3 <- mean(sign_rep2 == sign_rep3)

omega_sign_agreement <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  proportion_same_sign = c(same_sign_rep1_rep2,
                           same_sign_rep1_rep3,
                           same_sign_rep2_rep3)
)

omega_sign_agreement



#################################################
###  Extract Omega support species pairs      ###
#################################################



Omega_support_m2_rep1_pairs <- Omega_support_m2_rep1[
  lower.tri(Omega_support_m2_rep1)]

Omega_support_m2_rep2_pairs <- Omega_support_m2_rep2[
  lower.tri(Omega_support_m2_rep2)]

Omega_support_m2_rep3_pairs <- Omega_support_m2_rep3[
  lower.tri(Omega_support_m2_rep3)]



summary(Omega_support_m2_rep1_pairs)

summary(Omega_support_m2_rep2_pairs)

summary(Omega_support_m2_rep3_pairs)



#################################################
###  Omega support correlations               ###
#################################################



cor_support_m2_rep1_rep2 <- cor(
  Omega_support_m2_rep1_pairs,
  Omega_support_m2_rep2_pairs)

cor_support_m2_rep1_rep3 <- cor(
  Omega_support_m2_rep1_pairs,
  Omega_support_m2_rep3_pairs)

cor_support_m2_rep2_rep3 <- cor(
  Omega_support_m2_rep2_pairs,
  Omega_support_m2_rep3_pairs)


omega_support_correlations <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  correlation = c(cor_support_m2_rep1_rep2,
                  cor_support_m2_rep1_rep3,
                  cor_support_m2_rep2_rep3))

omega_support_correlations



#############################################
###  Omega support absolute differences   ###
#############################################

mad_support_m2_rep1_rep2 <- mean(
  abs(Omega_support_m2_rep1_pairs -
        Omega_support_m2_rep2_pairs))

mad_support_m2_rep1_rep3 <- mean(
  abs(Omega_support_m2_rep1_pairs -
        Omega_support_m2_rep3_pairs))

mad_support_m2_rep2_rep3 <- mean(
  abs(Omega_support_m2_rep2_pairs -
        Omega_support_m2_rep3_pairs))

omega_support_differences <- data.frame(
  comparison = c('rep1_vs_rep2',
                 'rep1_vs_rep3',
                 'rep2_vs_rep3'),
  mean_absolute_difference =
    c(mad_support_m2_rep1_rep2,
      mad_support_m2_rep1_rep3,
      mad_support_m2_rep2_rep3))

omega_support_differences


# Under the current MCMC settings, environmental responses are highly reproducible across repetitions. Residual association estimates are substantially less stable. Posterior support values show moderate reproducibility, with average differences below 0.1 support units. These results provide a baseline against which longer MCMC runs can be evaluated.


########################################
###  Summarise m2 stability results   ###
########################################


#combine all stability metrics into one table
m2_stability_summary <- data.frame(
  model = rep('m2_spatial_latent', 3),
  comparison = beta_m2_correlations$comparison,
  
  beta_mean_correlation = beta_m2_correlations$correlation,
  beta_mean_absolute_difference =
    beta_m2_differences$mean_absolute_difference,
  
  omega_mean_correlation =
    omega_pairs_m2_correlations$correlation,
  omega_mean_absolute_difference =
    omega_m2_differences$mean_absolute_difference,
  
  omega_sign_agreement =
    omega_sign_agreement$proportion_same_sign,
  
  omega_support_correlation =
    omega_support_correlations$correlation,
  omega_support_absolute_difference =
    omega_support_differences$mean_absolute_difference
)

#inspect summary table
m2_stability_summary

#setwd
setwd(wd_model_comparison)

#save stability summary
write.csv(m2_stability_summary,
          'm2_stability_summary.csv',
          row.names = FALSE)


