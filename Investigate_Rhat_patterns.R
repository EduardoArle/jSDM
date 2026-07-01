#####################################
###  Detailed Rhat diagnostics    ###
#####################################



#load packages
library(Hmsc); library(coda)

#list WDs
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Data'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Sensitivity analysis jSDMs reptiles/Models'
wd_output <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/HMSC/Model_stability_tests/Rhat_detailed_diagnostics'
wd_beta_rhat_species_plots <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/HMSC/Model_stability_tests/Beta_Rhat_species_plots'
wd_omega_rhat_species_plots <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/HMSC/Model_stability_tests/Omega_Rhat_species_plots'
  


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

#function to extract species name from Beta parameter names
extract_beta_species <- function(parameter_names){
  
  #extract species name
  species <- sub('.*\\), ([A-Za-z0-9]+) \\(S[0-9]+\\).*',
                 '\\1',
                 parameter_names)
  
  #return species
  return(species)
}

#function to extract predictor name from Beta parameter names
extract_beta_predictor <- function(parameter_names){
  
  #extract predictor name
  predictor <- sub('B\\[([^ ]+).*',
                   '\\1',
                   parameter_names)
  
  #return predictor
  return(predictor)
}

#function to plot Beta Rhat across settings for one species
plot_beta_rhat_species <- function(beta_rhat_all, species_name){
  
  #define setting order
  setting_order <- unique(beta_rhat_all$setting[order(beta_rhat_all$samples,
                                                      beta_rhat_all$nChains)])
  
  #define predictor order
  predictor_order <- c('(Intercept)', 'T_min', 'Precipitation', 'Elev_MEAN',
                       'Barrenness', 'Ndvi_MEAN', 'DISTURB')
  
  #subset species
  beta_sp <- beta_rhat_all[beta_rhat_all$species == species_name, ]
  
  #set layout
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  #define repetition colours
  rep_cols <- c('lightblue', 'pink', 'lightgreen')
  
  #loop over predictors
  for(i in 1:length(predictor_order))
  {
    
    #subset predictor
    beta_pred <- beta_sp[beta_sp$predictor == predictor_order[i], ]
    
    #set y axis limit
    y_max <- max(beta_pred$rhat, na.rm = TRUE)
    
    #create empty plot
    plot(1:length(setting_order), rep(NA, length(setting_order)),
         ylim = c(1, y_max), xaxt = 'n', xlab = '',
         ylab = 'Beta Rhat', main = predictor_order[i])
    
    #add x axis
    axis(1, at = 1:length(setting_order), labels = F)

    #add angled labels
    for(j in 1:length(setting_order))
    {
      
      text(x = j, #y = par('usr')[3] - 0.5,
           y = par('usr')[3] - 0.12 * diff(par('usr')[3:4]),
           labels = setting_order[j], srt = 35,
           adj = 1, xpd = NA, cex = 0.8)
      }
    
    #add reference lines
    abline(h = 1.1, lty = 2)
    
    #create matrix to store Rhat values
    rhat_mat <- matrix(NA, nrow = length(setting_order),
                       ncol = length(unique(beta_pred$rep)))
    
    #loop over repetitions
    for(j in unique(beta_pred$rep))
    {
      
      #subset repetition
      beta_rep <- beta_pred[beta_pred$rep == j, ]
      
      #order by setting
      beta_rep <- beta_rep[match(setting_order, beta_rep$setting), ]
      
      #add line
      lines(1:length(setting_order), beta_rep$rhat, type = 'l',
            col = rep_cols[j], lwd = 1.5)
      
      #store Rhat values
      rhat_mat[, j] <- beta_rep$rhat
    }
    
    #add mean line
    lines(1:length(setting_order), rowMeans(rhat_mat, na.rm = TRUE),
          col = 'grey60', lwd = 2)

  }
  
  #create empty panel for legend
  plot.new()
  
  #add legend
  legend('center', legend = c('Rep 1', 'Rep 2', 'Rep 3', 'Mean', 'Rhat = 1.1'),
         col = c(rep_cols, 'grey60', 'black'), lty = c(1, 1, 1, 1, 2),
         lwd = c(1.5, 1.5, 1.5, 1.75, 1.5), bty = 'n', cex = 0.8,
         y.intersp = 0.8)
  
  #add species title
  mtext(species_name, outer = TRUE, line = 1, cex = 1.4)

}

#function to plot Omega Rhat across settings for one species
plot_omega_rhat_species <- function(omega_rhat_all, species_name){
  
  #define setting order
  setting_order <- unique(omega_rhat_all$setting[order(omega_rhat_all$samples,
                                                       omega_rhat_all$nChains)])
  
  #subset species pairs
  omega_sp <- omega_rhat_all[omega_rhat_all$species1 == species_name |
                               omega_rhat_all$species2 == species_name, ]
  
  #define pair order
  pair_order <- unique(omega_sp$pair)
  
  #set layout
  par(mfrow = c(4, 5), mar = c(4, 4, 3, 1), oma = c(0, 0, 3, 0))
  
  #define repetition colours
  rep_cols <- c('lightblue', 'pink', 'lightgreen')
  
  #loop over pairs
  for(i in 1:length(pair_order))
  {
    
    #subset pair
    omega_pair <- omega_sp[omega_sp$pair == pair_order[i], ]
    
    #remove duplicated symmetric entries
    omega_pair <- omega_pair[!duplicated(omega_pair[,c('pair','setting',
                                                       'rep')]), ]
    
    #set y axis limit
    y_max <- max(omega_pair$rhat, na.rm = TRUE)
    
    #create empty plot
    plot(1:length(setting_order), rep(NA, length(setting_order)),
         ylim = c(1, y_max), xaxt = 'n', xlab = '',
         ylab = 'Omega Rhat', main = pair_order[i])
    
    #add x axis
    axis(1, at = 1:length(setting_order), labels = F)
    
    #add angled labels
    for(j in 1:length(setting_order))
    {
      
      text(x = j, y = par('usr')[3] - 0.12 * diff(par('usr')[3:4]),
           labels = setting_order[j], srt = 35,
           adj = 1, xpd = NA, cex = 0.8)
    }
    
    #add reference line
    abline(h = 1.1, lty = 2)
    
    #create matrix to store Rhat values
    rhat_mat <- matrix(NA, nrow = length(setting_order),
                       ncol = length(unique(omega_pair$rep)))
    
    #loop over repetitions
    for(j in unique(omega_pair$rep))
    {
      
      #subset repetition
      omega_rep <- omega_pair[omega_pair$rep == j, ]
      
      #order by setting
      omega_rep <- omega_rep[match(setting_order, omega_rep$setting), ]
      
      #add line
      lines(1:length(setting_order), omega_rep$rhat, type = 'l',
            col = rep_cols[j], lwd = 1.5)
      
      #store Rhat values
      rhat_mat[, j] <- omega_rep$rhat
    }
    
    #add mean line
    lines(1:length(setting_order), rowMeans(rhat_mat, na.rm = TRUE),
          col = 'grey60', lwd = 2)
  }
  
  #create empty panel for legend
  plot.new()
  
  #add legend
  legend('center', legend = c('Rep 1', 'Rep 2', 'Rep 3', 'Mean', 'Rhat = 1.1'),
         col = c(rep_cols, 'grey60', 'black'), lty = c(1, 1, 1, 1, 2),
         lwd = c(1.5, 1.5, 1.5, 1.75, 1.5), bty = 'n', cex = 1,
         y.intersp = 0.7)
  
  #add species title
  mtext(species_name, outer = TRUE, line = 1, cex = 1.4)
}



#################################
###  Find mpost files         ###
#################################



#setwd
setwd(wd_models)

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



########################################
###  Extract all Beta Rhat values    ###
########################################



#empty object to store Beta Rhat values
beta_rhat_all <- data.frame()

#loop over mpost files
for(i in 1:nrow(mpost_metadata)){
  
  cat('\n')
  cat('=================================\n')
  cat('Processing model', i, 'of', nrow(mpost_metadata), '\n')
  cat('=================================\n')
  
  
  #load mpost object
  mpost_i <- readRDS(mpost_metadata$file_name[i])
  
  
  #calculate Beta Rhat
  beta_rhat_i <- gelman.diag(mpost_i$Beta,
                             multivariate = FALSE)
  
  
  #extract Rhat values
  beta_rhat_values <- beta_rhat_i$psrf[, 'Point est.']
  
  
  #create table for current model
  beta_rhat_table_i <- data.frame(
    model = mpost_metadata$model[i],
    samples = mpost_metadata$samples[i],
    transient = mpost_metadata$transient[i],
    thin = mpost_metadata$thin[i],
    nChains = mpost_metadata$nChains[i],
    rep = mpost_metadata$rep[i],
    parameter = names(beta_rhat_values),
    rhat = as.numeric(beta_rhat_values)
  )
  
  
  #append to main table
  beta_rhat_all <- rbind(beta_rhat_all,
                         beta_rhat_table_i)
}


#extract species names
beta_rhat_all$species <- extract_beta_species(beta_rhat_all$parameter)

#extract predictor names
beta_rhat_all$predictor <- extract_beta_predictor(beta_rhat_all$parameter)

#inspect table
head(beta_rhat_all)

#create short model setting label
beta_rhat_all$setting <- paste0(beta_rhat_all$samples,
                                '/',
                                beta_rhat_all$transient,
                                '/',
                                beta_rhat_all$nChains)

#inspect labels
unique(beta_rhat_all$setting)



##################################################
###  Beta Rhat table by repetition             ###
##################################################



#create table with one column per repetition
beta_rhat_reps <- reshape(beta_rhat_all[, c('species',
                                            'predictor',
                                            'setting',
                                            'rep',
                                            'rhat')],
                          idvar = c('species',
                                    'predictor',
                                    'setting'),
                          timevar = 'rep',
                          direction = 'wide')

#rename repetition columns
colnames(beta_rhat_reps) <- sub('rhat.',
                                'rep',
                                colnames(beta_rhat_reps),
                                fixed = TRUE)

#define desired setting order
beta_rhat_reps$setting <- factor(beta_rhat_reps$setting,
                                 levels = c('250/250/2',
                                            '500/500/2',
                                            '1000/1000/2',
                                            '1000/1000/4',
                                            '3000/3000/4'))

#sort table
beta_rhat_reps <- beta_rhat_reps[
  order(beta_rhat_reps$species,
        beta_rhat_reps$predictor,
        beta_rhat_reps$setting),
]

#reset row names
rownames(beta_rhat_reps) <- NULL

#inspect table
head(beta_rhat_reps)



###################################################
###  Create wide Beta Rhat notebook table       ###
###################################################



#create setting-repetition label
beta_rhat_all$setting_rep <- paste0(beta_rhat_all$setting,
                                    '_rep',
                                    beta_rhat_all$rep)

#create wide table with one column per setting and repetition
beta_rhat_wide <- reshape(beta_rhat_all[, c('species',
                                            'predictor',
                                            'setting_rep',
                                            'rhat')],
                          idvar = c('species',
                                    'predictor'),
                          timevar = 'setting_rep',
                          direction = 'wide')

#clean column names
colnames(beta_rhat_wide) <- sub('rhat.',
                                '',
                                colnames(beta_rhat_wide),
                                fixed = TRUE)

#inspect table
head(beta_rhat_wide)

#check dimensions
dim(beta_rhat_wide)



##############################
###  Save Beta Rhat table  ###
##############################



#set output directory
setwd(wd_output)

#save full long table
write.csv(beta_rhat_all,
          'beta_rhat_all_long.csv',
          row.names = FALSE)

#save repetition table
write.csv(beta_rhat_reps,
          'beta_rhat_by_repetition.csv',
          row.names = FALSE)

#save wide notebook table
write.csv(beta_rhat_wide,
          'beta_rhat_wide_notebook.csv',
          row.names = FALSE)



#####################################
###  Save Beta Rhat species plots ###
#####################################



#set working directory
setwd(wd_beta_rhat_species_plots)

#define species names
species_names <- unique(beta_rhat_all$species)

#loop over species
for(i in 1:length(species_names))
{
  
  #define species
  species_i <- species_names[i]
  
  #save plot
  png(paste0('beta_rhat_', species_i, '.png'),
      width = 2200, height = 1800, res = 200)
  
  #plot species
  plot_beta_rhat_species(beta_rhat_all, species_i)
  
  #close device
  dev.off()
  
}



#################################
###  Calculate Omega Rhat     ###
#################################



#empty object to store results
omega_rhat_all <- data.frame()

#loop over models
for(i in 1:nrow(mpost_metadata))
{
  
  #read posterior
  setwd(wd_models)
  mpost_i <- readRDS(mpost_metadata$file_name[i])
  
  #calculate Omega Rhat
  omega_rhat_i <- gelman.diag(mpost_i$Omega[[1]],
                              multivariate = FALSE)
  
  #store results
  omega_rhat_all <- rbind(omega_rhat_all,
                          data.frame(model = mpost_metadata$model[i],
                                     samples = mpost_metadata$samples[i],
                                     transient = mpost_metadata$transient[i],
                                     thin = mpost_metadata$thin[i],
                                     nChains = mpost_metadata$nChains[i],
                                     rep = mpost_metadata$rep[i],
                                     parameter = rownames(omega_rhat_i$psrf),
                                     rhat = omega_rhat_i$psrf[, 'Point est.']))
  
  print(i)
}



##################################
###  Parse Omega parameters    ###
##################################



#extract first species
omega_rhat_all$species1 <- sub('Omega1\\[(.*) \\(S[0-9]+\\), .*',
                               '\\1',
                               omega_rhat_all$parameter)

#extract second species
omega_rhat_all$species2 <- sub('Omega1\\[.*, (.*) \\(S[0-9]+\\)\\]',
                               '\\1',
                               omega_rhat_all$parameter)

#inspect results
head(omega_rhat_all[, c('parameter', 'species1', 'species2', 'rhat')])



################################
###  Create Omega pair names ###
################################



#create unordered pair name
omega_rhat_all$pair <- 
  ifelse(omega_rhat_all$species1 <= omega_rhat_all$species2,
      paste(omega_rhat_all$species1, omega_rhat_all$species2, sep = ' - '),
      paste(omega_rhat_all$species2, omega_rhat_all$species1, sep = ' - '))

#inspect results
head(omega_rhat_all[, c('species1', 'species2', 'pair', 'rhat')])

#create setting label
omega_rhat_all$setting <- paste0(omega_rhat_all$samples, '/',
                                 omega_rhat_all$transient, '/',
                                 omega_rhat_all$nChains)



######################################
###  Save Omega Rhat species plots ###
######################################



#set working directory
setwd(wd_omega_rhat_species_plots)

#define species names
species_names <- sort(unique(c(omega_rhat_all$species1,
                               omega_rhat_all$species2)))

#loop over species
for(i in 1:length(species_names))
{
  
  #define species
  species_i <- species_names[i]
  
  #save plot
  png(paste0('omega_rhat_', species_i, '.png'),
      width = 2600,
      height = 2000,
      res = 200)
  
  #plot species
  plot_omega_rhat_species(omega_rhat_all,
                          species_i)
  
  #close device
  dev.off()
  
  print(species_i)
}


