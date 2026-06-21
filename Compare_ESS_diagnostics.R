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
