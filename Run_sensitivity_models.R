############################################
###  Run sensitivity models: setup       ###
############################################


#load packages
library(Hmsc)

#list WDs
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Data'
wd_sensitivity_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Sensitivity analysis jSDMs reptiles/Models'

#setwd
setwd(wd_data)

#load main reptile table
reptiles <- read.csv('Phd_survey_samples_2024_withEgrets.csv')

#check data
dim(reptiles)
names(reptiles)



#####################################
###  Sensitivity settings         ###
#####################################



#model configuration label
#
#this label will be used in file names
model_label <- 'clean_spatial_latent_500samples_500transient_10thin_2chains'

#MCMC settings
samples <- 1000

transient <- 1000

thin <- 10

nChains <- 2

#number of independent repetitions
n_repetitions <- 10



############################################
###  Species matrix (clean version)      ###
############################################



#create community matrix
#
#reptile species plus cattle egret counts
Y <- reptiles[, c(15:38, 100)]

#rename egret column
colnames(Y)[colnames(Y) == 'BubibisSmpl'] <- 'CattleEgret'

#remove rare taxa
occ_sites <- colSums(Y > 0)

Y <- Y[, occ_sites >= 10]

#remove problematic taxa
#
#LizUnIdent:
#unresolved taxonomic identity
#
#Tesgra:
#turtle species with distinct ecology
Y <- Y[, !colnames(Y) %in% c('LizUnIdent',
                             'Tesgra')]

#convert to matrix
Y <- as.matrix(Y)

#check dimensions
dim(Y)

#check taxa
colnames(Y)



############################################
###  Environmental matrix                ###
############################################



#create environmental matrix
XData <- reptiles[, c('T_min',
                      'Precipitation',
                      'Elev_MEAN',
                      'Barrenness',
                      'Ndvi_MEAN',
                      'DISTURB')]

#standardise predictors
XScaled <- scale(XData)

#convert back to data frame
XScaled <- as.data.frame(XScaled)

#check dimensions
dim(XScaled)

#check means
colMeans(XScaled)

#check standard deviations
apply(XScaled, 2, sd)



############################################
###  Study design and coordinates        ###
############################################



#create study design
studyDesign <- data.frame(sample = as.factor(reptiles$Sample))

#extract coordinates
coords <- as.matrix(reptiles[, c('Lon_sample',
                                 'Lat_sample')])

#match coordinates to samples
rownames(coords) <- reptiles$Sample

#create spatial random level
rL.site <- HmscRandomLevel(sData = coords)

#check dimensions
dim(coords)

#inspect first rows
head(studyDesign)

head(coords)



#########################
###  Model formula    ###
#########################



#define model formula
XFormula <- ~ T_min +
  Precipitation +
  Elev_MEAN +
  Barrenness +
  Ndvi_MEAN +
  DISTURB

#inspect formula
XFormula



#################################
###  Create base HMSC model   ###
#################################



#create base model
m_base <- Hmsc(Y = Y,
               XData = XScaled,
               XFormula = XFormula,
               distr = 'lognormal poisson',
               studyDesign = studyDesign,
               ranLevels = list(sample = rL.site))

#inspect model
m_base



##################################
###  Run model repetitions     ###
##################################



#loop over repetitions
for(i in 4:n_repetitions){
  
  cat('\n')
  cat('=================================\n')
  cat('Running repetition', i, '\n')
  cat('=================================\n')
  
  #create repetition label
  rep_label <- paste0(model_label,
                      '_rep',
                      i)
  
  #show repetition label
  print(rep_label)
  
  #create output file names
  model_file <- paste0(rep_label,
                       '.rds')
  
  postBeta_file <- paste0('postBeta_',
                          rep_label,
                          '.rds')
  
  mpost_file <- paste0('mpost_',
                       rep_label,
                       '.rds')
  
  Omega_file <- paste0('Omega_',
                       rep_label,
                       '.rds')
  
  #show file names
  print(model_file)
  print(postBeta_file)
  print(mpost_file)
  print(Omega_file)
  
  
  
  ##########################################
  ### Fit model                          ###
  ##########################################
  
  
  
  m_fit <- sampleMcmc(m_base,
                      samples = samples,
                      thin = thin,
                      transient = transient,
                      nChains = nChains,
                      verbose = 1000)
  
  
  ##########################################
  ### Save fitted model                  ###
  ##########################################
  
  
  setwd(wd_sensitivity_models)
  
  saveRDS(m_fit,
          model_file)
  
  
  
  ##########################################
  ### Extract and save model outputs     ###
  ##########################################
  
  
  
  #extract posterior beta estimates
  postBeta_fit <- getPostEstimate(m_fit,
                                  parName = 'Beta')
  
  saveRDS(postBeta_fit,
          postBeta_file)
  
  
  #convert posterior samples to coda objects
  mpost_fit <- convertToCodaObject(m_fit)
  
  saveRDS(mpost_fit,
          mpost_file)
  
  
  #compute residual associations
  Omega_fit <- computeAssociations(m_fit)
  
  saveRDS(Omega_fit,
          Omega_file)
  
}


model_label

