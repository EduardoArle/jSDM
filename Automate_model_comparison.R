############################################
###  Compare model repetitions: setup    ###
############################################


#clean workspace
rm(list = ls())

#load packages
library(Hmsc)

#list WDs
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Sensitivity analysis jSDMs reptiles/Models'
wd_model_comparison <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Sensitivity analysis jSDMs reptiles/Stability_summaries'

#setwd
setwd(wd_models)

#check files in model folder
list.files(wd_models)



############################################
###  Choose model family to compare      ###
############################################



#model family to compare
#
#for now:
#
#m2_spatial_latent:
#full latent-variable model with all retained taxa
#
#later we will change this to m3_clean_spatial_latent
model_label <- 'clean_spatial_latent_3000samples_3000transient_10thin_4chains'

#file patterns for this model family
postBeta_pattern <- paste0('^postBeta_',
                           model_label,
                           '_rep[0-9]+\\.rds$')

Omega_pattern <- paste0('^Omega_',
                        model_label,
                        '_rep[0-9]+\\.rds$')



#################################
###  Find model output files  ###
#################################



#find postBeta files for selected model family
postBeta_files <- list.files(wd_models,
                             pattern = postBeta_pattern,
                             full.names = FALSE)

#find Omega files for selected model family
Omega_files <- list.files(wd_models,
                          pattern = Omega_pattern,
                          full.names = FALSE)

#inspect files
postBeta_files
Omega_files



###################################
###  Order files by repetition  ###
###################################



#function to identify repetition number from file name
#
#files without "_rep" are treated as rep1
get_rep_number <- function(file_name){
  
  if(grepl('_rep[0-9]+', file_name)){
    
    rep_number <- sub('.*_rep([0-9]+).*',
                      '\\1',
                      file_name)
    
    rep_number <- as.numeric(rep_number)
    
  } else {
    
    rep_number <- 1
  }
  
  return(rep_number)
}


#extract repetition numbers
postBeta_reps <- sapply(postBeta_files,
                        get_rep_number)

Omega_reps <- sapply(Omega_files,
                     get_rep_number)


#order files by repetition number
postBeta_files <- postBeta_files[order(postBeta_reps)]
Omega_files <- Omega_files[order(Omega_reps)]


#inspect ordered files
postBeta_files
Omega_files



############################################
###  Load model outputs                   ###
############################################



#load postBeta objects
postBeta_list <- lapply(postBeta_files,
                        readRDS)

#load Omega objects
Omega_list <- lapply(Omega_files,
                     readRDS)


#check number of repetitions found
length(postBeta_list)

length(Omega_list)



############################################
###  Extract comparison matrices          ###
############################################



#extract beta mean matrices from all repetitions
Beta_mean_list <- lapply(postBeta_list,
                         function(x) x$mean)


#extract Omega mean matrices from all repetitions
Omega_mean_list <- lapply(Omega_list,
                          function(x) x[[1]]$mean)


#extract Omega support matrices from all repetitions
Omega_support_list <- lapply(Omega_list,
                             function(x) x[[1]]$support)


#check dimensions
lapply(Beta_mean_list, dim)

lapply(Omega_mean_list, dim)

lapply(Omega_support_list, dim)



############################################
###  Create repetition comparisons        ###
############################################



#number of repetitions
n_reps <- length(Beta_mean_list)

#repetition labels
rep_labels <- paste0('rep', 1:n_reps)

#all pairwise comparisons among repetitions
rep_comparisons <- combn(rep_labels,
                         2,
                         simplify = FALSE)

#inspect comparisons
rep_comparisons



##################################
###  Test comparison indexing  ###
##################################



#take first comparison as a test
test_comp <- rep_comparisons[[1]]

#find list positions
i <- match(test_comp[1], rep_labels)
j <- match(test_comp[2], rep_labels)

#inspect
test_comp
i
j



################################
###  Compare beta means      ###
################################



#empty object to store results
beta_results <- data.frame()


#loop over repetition comparisons
for(k in seq_along(rep_comparisons)){
  
  #get comparison labels
  comp <- rep_comparisons[[k]]
  
  #get list positions
  i <- match(comp[1], rep_labels)
  j <- match(comp[2], rep_labels)
  
  #calculate correlation
  beta_cor <- cor(as.vector(Beta_mean_list[[i]]),
                  as.vector(Beta_mean_list[[j]]))
  
  #calculate mean absolute difference
  beta_mad <- mean(abs(Beta_mean_list[[i]] -
                         Beta_mean_list[[j]]))
  
  #store results
  beta_results <- rbind(beta_results,
                        data.frame(comparison = paste(comp[1],
                                                      comp[2],
                                                      sep = '_vs_'),
                                   beta_mean_correlation = beta_cor,
                                   beta_mean_absolute_difference = beta_mad))
}


#inspect results
beta_results



############################################
###  Compare Omega means                  ###
############################################



#empty object to store results
omega_mean_results <- data.frame()


#loop over repetition comparisons
for(k in seq_along(rep_comparisons)){
  
  #get comparison labels
  comp <- rep_comparisons[[k]]
  
  #get list positions
  i <- match(comp[1], rep_labels)
  j <- match(comp[2], rep_labels)
  
  #extract species-pair values only
  omega_i <- Omega_mean_list[[i]][lower.tri(Omega_mean_list[[i]])]
  omega_j <- Omega_mean_list[[j]][lower.tri(Omega_mean_list[[j]])]
  
  #calculate correlation
  omega_cor <- cor(omega_i,
                   omega_j)
  
  #calculate mean absolute difference
  omega_mad <- mean(abs(omega_i -
                          omega_j))
  
  #store results
  omega_mean_results <- rbind(omega_mean_results,
                              data.frame(comparison = paste(comp[1],
                                                            comp[2],
                                                            sep = '_vs_'),
                                         omega_mean_correlation = omega_cor,
                                         omega_mean_absolute_difference = omega_mad))
}


#inspect results
omega_mean_results



############################################
###  Compare Omega signs                  ###
############################################



#empty object to store results
omega_sign_results <- data.frame()


#loop over repetition comparisons
for(k in seq_along(rep_comparisons)){
  
  #get comparison labels
  comp <- rep_comparisons[[k]]
  
  #get list positions
  i <- match(comp[1], rep_labels)
  j <- match(comp[2], rep_labels)
  
  #extract species-pair values only
  omega_i <- Omega_mean_list[[i]][lower.tri(Omega_mean_list[[i]])]
  omega_j <- Omega_mean_list[[j]][lower.tri(Omega_mean_list[[j]])]
  
  #calculate sign agreement
  sign_agreement <- mean(sign(omega_i) == sign(omega_j))
  
  #store results
  omega_sign_results <- rbind(omega_sign_results,
                              data.frame(comparison = paste(comp[1],
                                                            comp[2],
                                                            sep = '_vs_'),
                                         omega_sign_agreement = sign_agreement))
}


#inspect results
omega_sign_results



#####################################
###  Compare Omega support        ###
#####################################



#empty object to store results
omega_support_results <- data.frame()


#loop over repetition comparisons
for(k in seq_along(rep_comparisons)){
  
  #get comparison labels
  comp <- rep_comparisons[[k]]
  
  #get list positions
  i <- match(comp[1], rep_labels)
  j <- match(comp[2], rep_labels)
  
  #extract species-pair support values only
  support_i <- Omega_support_list[[i]][lower.tri(Omega_support_list[[i]])]
  support_j <- Omega_support_list[[j]][lower.tri(Omega_support_list[[j]])]
  
  #calculate correlation
  support_cor <- cor(support_i,
                     support_j)
  
  #calculate mean absolute difference
  support_mad <- mean(abs(support_i -
                            support_j))
  
  #store results
  omega_support_results <- rbind(omega_support_results,
                                 data.frame(comparison = paste(comp[1],
                                                               comp[2],
                                                               sep = '_vs_'),
                                            omega_support_correlation = support_cor,
                                            omega_support_absolute_difference = support_mad))
}


#inspect results
omega_support_results



############################################
###  Summarise stability results          ###
############################################



#combine all stability metrics into one table
stability_summary <- merge(beta_results,
                           omega_mean_results,
                           by = 'comparison')

stability_summary <- merge(stability_summary,
                           omega_sign_results,
                           by = 'comparison')

stability_summary <- merge(stability_summary,
                           omega_support_results,
                           by = 'comparison')


#add model label
stability_summary$model <- model_label

#reorder columns
stability_summary <- stability_summary[, c('model',
                                           'comparison',
                                           'beta_mean_correlation',
                                           'beta_mean_absolute_difference',
                                           'omega_mean_correlation',
                                           'omega_mean_absolute_difference',
                                           'omega_sign_agreement',
                                           'omega_support_correlation',
                                           'omega_support_absolute_difference')]

#inspect summary table
stability_summary


#setwd
setwd(wd_model_comparison)

#save stability summary
write.csv(stability_summary,
          paste0(model_label, '_stability_summary.csv'),
          row.names = FALSE)

