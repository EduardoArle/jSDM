########################################
###  Compare ESS diagnostics: setup  ###
########################################


#clean workspace
rm(list = ls())

#load packages
library(Hmsc); library(coda)

#list WDs
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Sensitivity analysis jSDMs reptiles/Models'
wd_ess_diagnostics <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/HMSC/Model_stability_tests/ESS_diagnostics'



#####################
###  Functions    ###
#####################



#function to extract model settings from mpost file name
extract_model_settings <- function(file_name){
  
  #extract samples
  samples <- sub('.*_([0-9]+)samples_.*',
                 '\\1',
                 file_name)
  
  samples <- as.numeric(samples)
  
  
  #extract transient
  transient <- sub('.*samples_([0-9]+)transient_.*',
                   '\\1',
                   file_name)
  
  transient <- as.numeric(transient)
  
  
  #extract thin
  thin <- sub('.*transient_([0-9]+)thin_.*',
              '\\1',
              file_name)
  
  thin <- as.numeric(thin)
  
  
  #extract number of chains
  nChains <- sub('.*thin_([0-9]+)chains_.*',
                 '\\1',
                 file_name)
  
  nChains <- as.numeric(nChains)
  
  
  #extract repetition number
  rep <- sub('.*_rep([0-9]+)\\.rds$',
             '\\1',
             file_name)
  
  rep <- as.numeric(rep)
  
  
  #create model label without mpost_ and without repetition
  model <- sub('^mpost_',
               '',
               file_name)
  
  model <- sub('_rep[0-9]+\\.rds$',
               '',
               model)
  
  
  #return settings
  return(data.frame(file_name = file_name,
                    model = model,
                    samples = samples,
                    transient = transient,
                    thin = thin,
                    nChains = nChains,
                    rep = rep))
}

#function to calculate ESS summary
calculate_ESS_summary <- function(ess_values){
  
  #summarise ESS values
  ess_summary <- data.frame(ess_min = min(ess_values),
                            ess_q1 = quantile(ess_values, 0.25),
                            ess_median = median(ess_values),
                            ess_mean = mean(ess_values),
                            ess_q3 = quantile(ess_values, 0.75),
                            ess_max = max(ess_values))
  
  #return summary
  return(ess_summary)
}

#function to calculate relative ESS
calculate_relative_ESS <- function(ess_summary,
                                   max_ESS){
  
  #calculate relative ESS metrics
  ess_summary$relative_ess_min <- ess_summary$ess_min / max_ESS
  
  ess_summary$relative_ess_q1 <- ess_summary$ess_q1 / max_ESS
  
  ess_summary$relative_ess_median <- ess_summary$ess_median / max_ESS
  
  ess_summary$relative_ess_mean <- ess_summary$ess_mean / max_ESS
  
  ess_summary$relative_ess_q3 <- ess_summary$ess_q3 / max_ESS
  
  ess_summary$relative_ess_max <- ess_summary$ess_max / max_ESS
  
  #return table
  return(ess_summary)
}

#function to prepare ESS output for slides
prepare_ESS_slide_output <- function(ess_summary,
                                     max_ESS){
  
  #create simplified ESS table
  ess_output <- data.frame(ess_median = ess_summary$ess_median,
               relative_ess_median = ess_summary$ess_median / max_ESS)
  
  #return table
  return(ess_output)
}


#################################
###  Find mpost files         ###
#################################



#setwd
setwd(wd_models)

#check files in model folder
list.files(wd_models)

#find mpost files
mpost_files <- list.files(wd_models,
                pattern = '^mpost_clean_spatial_latent_.*_rep[0-9]+\\.rds$',
                full.names = FALSE)

#inspect mpost files
mpost_files

#check number of mpost files
length(mpost_files)



############################################
###  Create metadata table                ###
############################################



#extract settings from all mpost files
mpost_metadata <- do.call(rbind,
                          lapply(mpost_files,
                                 extract_model_settings))

#order metadata table
mpost_metadata <- mpost_metadata[order(mpost_metadata$samples,
                                       mpost_metadata$nChains,
                                       mpost_metadata$rep), ]

#reset row names
rownames(mpost_metadata) <- NULL

#inspect metadata table
mpost_metadata

#check dimensions
dim(mpost_metadata)



############################################
###  Calculate maximum possible ESS       ###
############################################



#calculate maximum possible ESS
#
#in HMSC, samples = retained samples per chain
mpost_metadata$max_ESS <- mpost_metadata$samples *
  mpost_metadata$nChains

#inspect maximum possible ESS
mpost_metadata[, c('model',
                   'rep',
                   'samples',
                   'nChains',
                   'max_ESS')]



#############################################
###  Calculate ESS for all model outputs  ###
#############################################


#empty object to store ESS results
ess_results <- data.frame()


#loop over mpost files
for(i in 1:nrow(mpost_metadata)){
  
  #load mpost object
  mpost_i <- readRDS(mpost_metadata$file_name[i])
  
  #calculate maximum possible ESS
  max_ESS_i <- mpost_metadata$max_ESS[i]
  
  
  ##################
  ###  Beta ESS  ###
  ##################
  
  
  #calculate Beta ESS
  beta_ess <- effectiveSize(mpost_i$Beta)
  
  #summarise Beta ESS
  beta_summary <- calculate_ESS_summary(beta_ess)
  
  #prepare Beta output
  beta_output <- prepare_ESS_slide_output(beta_summary,
                                          max_ESS = max_ESS_i)
  

  ###################
  ###  Omega ESS  ###
  ###################
  
  
  #calculate Omega ESS
  omega_ess <- effectiveSize(mpost_i$Omega[[1]])
  
  #summarise Omega ESS
  omega_summary <- calculate_ESS_summary(omega_ess)
  
  #prepare Omega output
  omega_output <- prepare_ESS_slide_output(omega_summary,
                                           max_ESS = max_ESS_i)
  
  
  ########################
  ###  Store results   ###
  ########################
  
  
  #store results
  ess_results <- rbind(ess_results,
                       data.frame(model = mpost_metadata$model[i],
                                  rep = mpost_metadata$rep[i],
                                  samples = mpost_metadata$samples[i],
                                  transient = mpost_metadata$transient[i],
                                  thin = mpost_metadata$thin[i],
                                  nChains = mpost_metadata$nChains[i],
                                  max_ESS = max_ESS_i,
                                  beta_median_ESS = beta_output$ess_median,
                                  beta_relative_median_ESS = beta_output$relative_ess_median,
                                  omega_median_ESS = omega_output$ess_median,
                                  omega_relative_median_ESS = omega_output$relative_ess_median))
}

#inspect results
ess_results

#setwd
setwd(wd_ess_diagnostics)

#save ESS results
write.csv(ess_results,
          'ESS_diagnostics_all_models.csv',
          row.names = FALSE)


