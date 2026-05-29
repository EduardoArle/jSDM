############################
###  Load reptile data   ###
############################


#clean workspace
rm(list = ls())

#load packages
library(Hmsc)

#list WDs
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Data'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_3rd_models'
wd_plots <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_3rd_models/Plots'

#load main reptile table
setwd(wd_data)
reptiles <- read.csv('Phd_survey_samples_2024_withEgrets.csv')

#check table dimensions
dim(reptiles)

#check column names
names(reptiles)

#look at the first rows
head(reptiles)



############################
###  Inspect columns     ###
############################



#show column names with their position
data.frame(column_number = 1:ncol(reptiles), column_name = names(reptiles))



############################
### Species matrix (Y)   ###
############################



#create community matrix

#reptile species plus cattle egret counts, BubibisSmpl is included as an additional species
#rather than as an environmental predictor
Y <- reptiles[, c(15:38, 100)]

#rename egret column for clarity in community plots
colnames(Y)[colnames(Y) == 'BubibisSmpl'] <- 'CattleEgret'

#check dimensions
dim(Y)

#look at first rows
head(Y)

#calculate total abundance of each species
colSums(Y)

#calculate number of occupied sites per species
colSums(Y > 0)


############################
### Filter rare species  ###
############################



#calculate the number of occupied sites per taxon
#this counts how many sites have abundance greater than zero
#taxa with very few occupied sites provide little information for estimating
#environmental responses and residual associations
occ_sites <- colSums(Y > 0)

#inspect prevalence before filtering
occ_sites

#keep only taxa occurring in at least 10 sites
#this filters rare reptile taxa while retaining CattleEgret
#because it occurs in enough sites
Y <- Y[, occ_sites >= 10]

#check which taxa remain after filtering
names(Y)

#check dimensions of the filtered community matrix
dim(Y)

#check prevalence again after filtering
colSums(Y > 0)


#remove problematic taxa suggested during exploratory discussion
#
#LizUnIdent:
#unresolved taxonomic identity
#
#Tesgra:
#turtle species with potentially distinct ecology
#
#for now, we create an alternative community matrix
#without modifying the original Y object
#
#importantly, CattleEgret remains in Y_m3
#because it is now treated as part of the community response matrix
Y_m3 <- Y[, !colnames(Y) %in% c('LizUnIdent',
                                'Tesgra')]

#check remaining taxa
colnames(Y_m3)

#check dimensions
dim(Y_m3)



############################################
### Raw species–environment associations ###
############################################



#select environmental predictors included in models
env_vars <- c('T_min',
              'Precipitation',
              'Elev_MEAN',
              'Barrenness',
              'Ndvi_MEAN',
              'DISTURB')

#extract environmental matrix
X_env <- reptiles[, env_vars]

#extract species abundance matrix
Y_species <- Y

#compute Spearman correlations between species abundances
#and environmental predictors
species_env_cor <- matrix(NA,
                          nrow = length(env_vars),
                          ncol = ncol(Y_species))

#assign row and column names
rownames(species_env_cor) <- env_vars
colnames(species_env_cor) <- colnames(Y_species)

#compute correlations
for(i in 1:length(env_vars)){
  
  for(j in 1:ncol(Y_species)){
    
    species_env_cor[i, j] <- cor(X_env[, i],
                                 Y_species[, j],
                                 method = 'spearman',
                                 use = 'complete.obs')
  }
}



###################################
### Prepare correlation heatmap ###
###################################



#transpose for plotting
species_env_cor_plot <- t(species_env_cor)

#define colour scale
species_env_zlim <- max(abs(species_env_cor_plot),
                        na.rm = TRUE)

species_env_cols <- colorRampPalette(c('darkorchid4',
                                       'white',
                                       'darkgreen'))(100)

species_env_legend_vals <- seq(-species_env_zlim,
                               species_env_zlim,
                               length.out = 100)




#########################################
### Raw species-species correlations  ###
#########################################



#compute raw species-species correlations
#
#Spearman is used because abundance data are count-based,
#zero-heavy, and not normally distributed
species_cor <- cor(Y,
                   method = 'spearman',
                   use = 'pairwise.complete.obs')

#prepare matrix for visualisation
#
#remove diagonal and upper triangle
species_cor_plot <- species_cor

species_cor_plot[upper.tri(species_cor_plot,
                           diag = FALSE)] <- NA

#define colour scale
species_cor_zlim <- max(abs(species_cor_plot),
                        na.rm = TRUE)

species_cor_cols <- colorRampPalette(c('darkorchid4',
                                       'white',
                                       'darkgreen'))(100)

species_cor_legend_vals <- seq(-species_cor_zlim,
                               species_cor_zlim,
                               length.out = 100)



###############################
### Plot species prevalence ###
###############################



#calculate species prevalence after filtering
species_prev <- colSums(Y > 0)

#plot prevalence of retained species
#
#this visualises how many sites each retained species occupies
#after applying the filtering threshold
#
#the dashed line indicates the minimum prevalence threshold used
par(mar = c(8, 5, 2, 1))

barplot(species_prev,
        las = 2,
        ylab = 'Occupied sites',
        main = 'Retained community taxa')

abline(h = 10,
       lty = 2)



############################
### Environmental matrix ###
############################



#create environmental matrix (X)
#this object contains the environmental predictors that will be used
#to explain variation in reptile community composition across sites
XData <- reptiles[, c('T_min',
                      'Precipitation',
                      'Elev_MEAN',
                      'Barrenness',
                      'Ndvi_MEAN',
                      'DISTURB')]

#rationale for the selected predictors:
#
#T_min:
#represents thermal limitation and broad climatic conditions
#
#Precipitation:
#represents water availability and climatic moisture
#
#Elev_MEAN:
#captures topographic variation and associated environmental gradients
#
#Barrenness:
#describes habitat openness and exposed substrate
#
#Ndvi_MEAN:
#proxy for vegetation productivity and greenness
#
#DISTURB:
#represents anthropogenic disturbance intensity

#importantly, we are NOT including:
#
#PreyLizAbun
#HerpAbun
#SpeciesRich
#
#because these variables are derived from the reptile community itself
#including them would create circularity in the model

#inspect first rows of environmental matrix
#this allows us to verify that variables were imported correctly
head(XData)

#summarise environmental variables
#this helps us inspect:
#
#variable ranges
#possible missing values
#differences in scale among predictors
#potentially problematic values
summary(XData)



############################
### Standardise XData    ###
############################



#environmental predictors are currently measured on very different scales
#
#for example:
#
#Precipitation may vary by hundreds of units
#Ndvi_MEAN varies roughly between 0 and 1
#Elev_MEAN may vary by hundreds of metres
#
#if predictors remain on different scales:
#
#coefficients become difficult to compare
#some variables may dominate numerically
#model fitting may become unstable
#MCMC mixing may become poorer
#
#therefore, it is standard practice in jSDMs and hierarchical models
#to centre and scale predictors before modelling

#standardise environmental predictors
#
#each variable will now have:
#
#mean = 0
#standard deviation = 1
#
#conceptually:
#
#X_scaled = (X - mean(X)) / sd(X)
#
#this means the model now interprets predictors in units of
#standard deviations rather than original measurement units
XScaled <- scale(XData)

#convert scaled matrix back into a data frame
#
#scale() returns a matrix object
#converting back to data.frame makes later manipulation easier
XScaled <- as.data.frame(XScaled)

#inspect first rows of scaled predictors
#
#values are now dimensionless and centred around zero
head(XScaled)

#check means of standardised predictors
#
#values should be extremely close to zero
#small numerical deviations are normal
colMeans(XScaled)

#check standard deviations of standardised predictors
#
#all values should now be 1
apply(XScaled, 2, sd)

#important conceptual point:
#
#after standardisation:
#
#a coefficient no longer means:
#"effect per one millimetre of rain"
#
#instead, it means:
#"effect per one standard deviation increase in precipitation"
#
#this makes effect sizes across predictors much more comparable



######################
### Study design   ###
######################



#each row in Y and XScaled corresponds to one sampling site
#
#Hmsc needs a study design object that identifies the sampling units
#this becomes especially important later when we add random effects
#such as spatial structure, locality, region, or repeated sampling design

#create site identifiers
#
#we use the Sample column because it identifies the sampling units
studyDesign <- data.frame(sample = as.factor(reptiles$Sample))

#inspect first rows
head(studyDesign)

#check number of sampling units
#
#this should match the number of rows in Y and XScaled
nrow(studyDesign)
nrow(Y)
nrow(XScaled)

#important conceptual point:
#
#for now, this object is very simple
#it only says:
#
#"each row is a sampling site"
#
#later, if we add random effects, HMSC will use this object
#to connect observations to those random effects



############################
### Model setup          ###
############################



#convert Y to a matrix
#
#Hmsc expects the response data to be a matrix
#rows are sites
#columns are species
Y <- as.matrix(Y)

#define model formula
#
#this says that species abundances are modelled as a function of:
#
#minimum temperature
#precipitation
#elevation
#egret field counts
#barrenness
#NDVI
#disturbance
#
#at this stage, we are fitting only environmental responses
#we are not yet adding latent variables or spatial random effects
XFormula <- ~ T_min +
  Precipitation +
  Elev_MEAN +
  Barrenness +
  Ndvi_MEAN +
  DISTURB


#############################
### Fix spatial matching  ###
#############################



#check whether sample identifiers are unique
length(unique(reptiles$Sample))
nrow(reptiles)

#recreate study design
#
#each row is one sampling site
studyDesign <- data.frame(sample = as.factor(reptiles$Sample))

#create coordinate matrix again
coords <- as.matrix(reptiles[, c('Lon_sample', 'Lat_sample')])

#give coordinates row names that match the sample names
#
#this is essential because HMSC needs to know which coordinates belong
#to which sampling unit
rownames(coords) <- levels(studyDesign$sample)

#check matching
head(rownames(coords))
head(levels(studyDesign$sample))

#create spatial random level again
rL.site <- HmscRandomLevel(sData = coords)

#rebuild latent-variable model
m2 <- Hmsc(Y = Y,
           XData = XScaled,
           XFormula = XFormula,
           distr = 'poisson',
           studyDesign = studyDesign,
           ranLevels = list(sample = rL.site))

#inspect model
m2



############################
### First HMSC model     ###
############################



#before fitting the model, we first define its structure
#
#this includes:
#
#Y = the species abundance matrix
#XData = the environmental predictors
#XFormula = the model formula
#distr = the response distribution


#create HMSC model object
#
#distr = 'poisson' because the species data are counts
#
#conceptually, this model asks:
#
#"how does each reptile species respond to the measured environmental
#predictors (including egret abundance/activity?) == NO
#
#this is still a relatively simple model
#that is good, because we want to understand the baseline first
m1 <- Hmsc(Y = Y,
           XData = XScaled,
           XFormula = XFormula,
           distr = 'poisson',
           studyDesign = studyDesign)

#inspect model object
m1



#############################
### Fit first HMSC model  ###
#############################



#fit the first model using MCMC sampling
#
#this is the step where the model is actually estimated
#
#because this is our first test run, we use a small number of samples
#the goal is NOT to get final results yet
#the goal is simply to check that:
#
#the data are correctly formatted
#the model runs without errors
#the workflow is working

m1 <- sampleMcmc(m1,
                 samples = 100,
                 thin = 10,
                 transient = 100,
                 nChains = 2,
                 verbose = 100)

#what these arguments mean:
#
#samples = 100
#number of posterior samples saved per chain
#this is very low and only suitable for testing
#
#thin = 10
#keeps one sample every 10 MCMC iterations
#
#transient = 100
#burn-in period
#early iterations are discarded because the chain is still adapting
#
#nChains = 2
#runs two independent chains
#later we use multiple chains to assess convergence
#
#verbose = 100
#print progress every 100 iterations

#important conceptual point:
#
#this is NOT the final model
#
#this is only a smoke test
#if it runs, we know that the model structure, response matrix,
#predictor matrix, and study design are compatible



############################
### Basic MCMC checks    ###
############################



#check MCMC convergence diagnostics
#
#HMSC provides posterior support and convergence summaries
#but first we inspect the effective sample sizes and potential scale reduction
#factor for model parameters
#
#effective sample size tells us how much independent information we have
#from the MCMC samples
#
#potential scale reduction factor checks whether different chains are
#mixing towards the same posterior distribution

mpost <- convertToCodaObject(m1)

#effective sample size for beta parameters
#
#beta parameters are the environmental responses of each species
effectiveSize(mpost$Beta)

#Gelman-Rubin diagnostic for beta parameters
#
#values close to 1 are good
#with this tiny test model, we do NOT expect perfect convergence
gelman.diag(mpost$Beta, multivariate = FALSE)

#important conceptual point:
#
#because we only used:
#
#100 samples
#2 chains
#short burn-in
#
#these diagnostics are only a technical check
#not a serious convergence assessment
#
#for now, we only want to know:
#
#"did the model run and are the outputs accessible?"



#############################
### Fit longer test model ###
#############################



#the first model was only a smoke test
#
#now we increase the MCMC settings to get a more meaningful preliminary fit
#
#this is still NOT the final model
#but it should behave much better than the tiny first run

m1 <- sampleMcmc(m1,
                 samples = 1000,
                 thin = 10,
                 transient = 1000,
                 nChains = 4,
                 verbose = 1000)

#setwd
setwd(wd_models)

#save the environment-only model
#
#this avoids having to rerun the MCMC every time we reopen R
#this is especially important because HMSC models can take a long time to fit
saveRDS(m1, 'm1_environment_only.rds')

#extract posterior beta estimates for m1
#
#beta coefficients describe species responses to environmental predictors
postBeta_m1 <- getPostEstimate(m1, parName = 'Beta')

#save beta estimates
setwd(wd_models)
saveRDS(postBeta_m1, 'postBeta_m1_environment_only.rds')

#samples = 1000
#saves 1000 posterior samples per chain
#
#thin = 10
#keeps one sample every 10 iterations
#
#transient = 1000
#discards the first 1000 iterations as burn-in
#
#nChains = 4
#runs four chains, which gives better convergence diagnostics
#
#verbose = 1000
#print progress less frequently

setwd(wd_models)
mpost <- convertToCodaObject(m1)

#save posterior objects used for diagnostics
#
#these are useful because they allow us to inspect convergence again
#without recreating the coda objects every time
saveRDS(mpost, 'mpost_m1_environment_only.rds')


#inspect available beta parameter names
#
#this helps us choose parameters for visualising chain behaviour
colnames(as.matrix(mpost$Beta[[1]]))

effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta, multivariate = FALSE)



#########################
### Visualise ESS m1  ###
#########################



#extract effective sample sizes for beta parameters from m1
ess_beta_m1 <- effectiveSize(mpost$Beta)

#visualise ESS distribution for m1
par(mar = c(5, 5, 2, 1))

hist(ess_beta_m1,
     breaks = 20,
     main = 'Effective sample size of beta parameters: m1',
     xlab = 'ESS')



############################
###  Species responses   ###
############################



#extract posterior mean beta coefficients
#
#beta coefficients describe species responses to environmental predictors
#
#positive values:
#species tends to increase with predictor
#
#negative values:
#species tends to decrease with predictor
#
#values near zero:
#weak or uncertain relationship

postBeta <- getPostEstimate(m1, parName = 'Beta')

#inspect dimensions
dim(postBeta$mean)

#inspect posterior mean coefficients
postBeta$mean




#important conceptual point:
#
#at this stage, species responses are still being explained only
#by the measured environmental predictors
#
#we have NOT yet added:
#
#latent variables
#residual species associations
#spatial random effects
#species interaction structure
#
#therefore, species are currently modelled as conditionally independent
#given the environment
#
#conceptually, this is still relatively close to a stacked SDM framework
#
#this is useful because it allows us to first understand:
#
#how species respond to measured environmental gradients
#which predictors appear important
#whether the model captures sensible ecological patterns
#
#before introducing latent structure and covariance among species
#
#this step is extremely important because once latent variables are added,
#interpretation becomes much more difficult
#
#latent variables can absorb:
#
#missing environmental predictors
#spatial structure
#sampling artefacts
#shared habitat preferences
#biotic interactions
#or combinations of all of these
#
#therefore, understanding the environment-only model first
#provides an essential ecological baseline



################################
### Visualise beta matrix m1 ###
################################



#plot posterior mean beta coefficients
#
#rows = environmental predictors
#columns = community taxa
#
#colour represents the direction and magnitude of species responses
#
#positive values:
#species tends to increase with predictor
#
#negative values:
#species tends to decrease with predictor
#
#values near zero:
#weak or uncertain response
#
#important:
#this plot shows posterior mean beta values only
#it does not show uncertainty around those estimates

#manually define predictor names
pred_names <- c('Intercept',
                'T_min',
                'Precipitation',
                'Elev_MEAN',
                'Barrenness',
                'Ndvi_MEAN',
                'DISTURB')

#add predictor names to beta matrix
rownames(postBeta$mean) <- pred_names

#define colour scale
#
#blue = negative response
#white = weak or near-zero response
#red = positive response
#
#the colour scale is symmetrical around zero
#this makes positive and negative responses visually comparable
zlim <- max(abs(postBeta$mean))

cols <- colorRampPalette(c('blue', 'white', 'red'))(100)



#############################
### Visualise in R window ###
#############################



#show beta heatmap in R plotting window
#
#layout creates two plotting areas:
#
#1 = heatmap
#2 = colour legend
layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main heatmap
par(mar = c(5, 2, 3, 4))

image(t(postBeta$mean),
      axes = FALSE,
      col = cols,
      zlim = c(-zlim, zlim),
      main = 'Community-wide environmental responses')

#add species names on x axis
axis(1,
     at = seq(0, 1, length.out = ncol(postBeta$mean)),
     labels = colnames(postBeta$mean),
     las = 2,
     cex.axis = 0.7)

#add predictor names on y axis
axis(2,
     at = seq(0, 1, length.out = length(pred_names)),
     labels = rev(pred_names),
     las = 2,
     cex.axis = 0.8)

#colour legend
par(mar = c(5, 2, 3, 4))

legend_vals <- seq(-zlim, zlim, length.out = 100)

image(x = 1,
      y = legend_vals,
      z = matrix(legend_vals, nrow = 1),
      col = cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Posterior mean beta',
      side = 4,
      line = 2.8,
      cex = 0.9)

#important conceptual point:
#
#this plot is one of the core outputs of many jSDMs
#
#it summarises how each species responds to each environmental predictor
#
#species with similar response profiles may:
#
#share habitat preferences
#respond similarly to climate
#show similar predator sensitivity
#occupy similar ecological niches
#
#however, at this stage:
#
#similar responses do NOT imply interactions
#
#two species may appear similar simply because:
#
#they respond similarly to climate
#they prefer similar habitats
#they co-occur under the same environmental conditions
#
#this distinction becomes extremely important later when interpreting
#residual covariance and latent-variable structure
#
#one especially interesting predictor here is:
#
#EgrtPredPress
#
#because it allows us to inspect whether reptile species differ
#in their apparent responses to cattle egret abundance/activity



#############################
### Spatial random level  ###
#############################



#create coordinate matrix
#
#rows = sampling sites
#columns = spatial coordinates
#
#these coordinates will allow HMSC to model spatial latent structure
coords <- as.matrix(reptiles[, c('Lon_sample',
                                 'Lat_sample')])

#inspect first rows
head(coords)

#create spatial random level
#
#this allows the model to estimate latent spatial structure among sites
#
#conceptually:
#
#sites that are geographically close may resemble each other
#for reasons not fully captured by the measured predictors
#
#the latent variables attempt to capture this residual spatial structure
rL.site <- HmscRandomLevel(sData = coords)

#important conceptual point:
#
#we are NOT explicitly telling the model:
#
#"these species interact"
#
#instead, we are allowing the model to estimate residual structure
#among sites after accounting for measured environmental predictors
#
#this residual structure may reflect:
#
#missing predictors
#spatial autocorrelation
#species interactions
#shared habitat preferences
#sampling artefacts
#or combinations of these
#
#therefore, latent variables should be interpreted cautiously
#
#they are extremely useful statistically
#but often ambiguous biologically



#############################
### Latent-variable model ###
#############################



#now we build a new model that includes:
#
#environmental predictors
#PLUS
#latent spatial structure
#
#this is where the model becomes a "true" latent-variable jSDM

#create latent-variable HMSC model
m2 <- Hmsc(Y = Y,
           XData = XScaled,
           XFormula = XFormula,
           distr = 'poisson',
           studyDesign = studyDesign,
           ranLevels = list(sample = rL.site))

#important conceptual point:
#
#the environmental component still estimates species responses
#to measured predictors such as:
#
#temperature
#precipitation
#elevation
#barrenness
#NDVI
#disturbance
#
#however, the model can now ALSO estimate additional latent structure
#among sampling sites
#
#conceptually, the model now asks:
#
#"after accounting for measured environmental predictors,
#is there still structured community variation remaining?"
#
#if yes, latent variables attempt to capture that remaining structure
#
#this is extremely important because ecological datasets almost always
#contain unmeasured structure
#
#examples:
#
#missing habitat predictors
#microclimate
#historical effects
#dispersal limitation
#observer bias
#biotic interactions
#spatial autocorrelation
#
#one of the central challenges in jSDMs is that latent variables can
#absorb ALL of these simultaneously
#
#therefore:
#
#latent structure is statistically real
#but biologically ambiguous

#inspect model object
m2



#############################
### Fit latent model      ###
#############################



#fit the latent-variable model
#
#as before, this is still a preliminary model
#not a final publishable run
#
#because latent-variable models are more complex than environment-only models
#they may take longer and may need more MCMC iterations to converge well

m2 <- sampleMcmc(m2,
                 samples = 1000,
                 thin = 10,
                 transient = 1000,
                 nChains = 4,
                 verbose = 1000)

#setwd
setwd(wd_models)

#save the latent-variable model
#
#this model is slower because it includes spatial latent structure
#saving it allows us to inspect diagnostics, coefficients, and associations later
#without repeating the expensive sampling step
saveRDS(m2, 'm2_spatial_latent.rds')

#extract posterior beta estimates for m2
#
#this allows us to compare environmental responses before and after
#adding latent spatial structure
postBeta_m2 <- getPostEstimate(m2, parName = 'Beta')

#save beta estimates
setwd(wd_models)
saveRDS(postBeta_m2, 'postBeta_m2_spatial_latent.rds')

#convert posterior samples to coda objects
#
#this lets us inspect MCMC diagnostics
mpost2 <- convertToCodaObject(m2)

#save posterior objects used for diagnostics
#
#these are useful because they allow us to inspect convergence again
#without recreating the coda objects every time

setwd(wd_models)
saveRDS(mpost2, 'mpost_m2_spatial_latent.rds')

#inspect available beta parameter names
#
#this helps us choose parameters for visualising chain behaviour
colnames(as.matrix(mpost2$Beta[[1]]))


#check effective sample size for beta parameters
#
#these are the environmental responses of each species
effectiveSize(mpost2$Beta)

#check Gelman diagnostics for beta parameters
#
#values close to 1 are good
#with latent-variable models, convergence may be harder than in m1
gelman.diag(mpost2$Beta, multivariate = FALSE)

#important conceptual point:
#
#we are first checking convergence for beta parameters only
#
#this tells us whether the estimated environmental responses are behaving
#reasonably across chains
#
#later, we can also inspect latent variables and residual associations
#
#but one step at a time:
#
#first ask:
#
#"did adding latent structure make the environmental responses unstable?"



#########################
### Visualise ESS m2  ###
#########################



#extract effective sample sizes for beta parameters from m2
ess_beta_m2 <- effectiveSize(mpost2$Beta)

#visualise ESS distribution for m2
par(mfrow = c(1,1))
par(mar = c(5, 5, 2, 1))

hist(ess_beta_m2,
     breaks = 20,
     main = 'Effective sample size of beta parameters: m2',
     xlab = 'ESS')



############################
### Interpretation notes ###
############################



#model m1:
#
#environment-only model
#
#this model estimates species responses to measured predictors only
#there are no latent variables and no residual association structure
#
#conceptually, this is close to a stacked SDM:
#
#species abundances ~ measured environment
#
#this model is useful as a baseline because it tells us what can be explained
#by the variables we explicitly included

#model m2:
#
#spatial latent-variable jSDM
#
#this model includes the same environmental predictors as m1
#but also includes spatial latent structure through the random level
#
#conceptually:
#
#species abundances ~ measured environment + residual spatial structure
#
#the latent structure captures remaining spatial/community variation
#after accounting for the measured predictors

#important warning:
#
#latent variables do NOT automatically represent species interactions
#
#they can absorb:
#
#missing environmental predictors
#spatial autocorrelation
#sampling artefacts
#shared habitat preferences
#dispersal limitation
#biotic interactions
#or combinations of these
#
#therefore, latent structure is statistically useful
#but biologically ambiguous

#convergence interpretation:
#
#Gelman values close to 1 indicate that chains are mixing similarly
#values clearly above 1 suggest poor convergence
#
#in this preliminary m2 model, some parameters converge reasonably well
#especially for common species
#
#rarer species tend to show poorer convergence
#because the model has less information to estimate their responses
#
#latent-variable models are also harder to fit because environmental effects
#and latent spatial structure can compete to explain the same variation

#practical decision:
#
#for now, m2 is useful for learning and exploration
#but it should not be treated as a final publishable model
#
#before any serious ecological interpretation, we would need:
#
#longer MCMC chains
#stronger convergence checks
#possibly additional model simplification
#possibly additional filtering of rare species
#model comparison
#posterior uncertainty around coefficients and associations



################################
### Visualise beta matrix m2 ###
################################



#plot posterior mean beta coefficients from m2
#
#m2 includes the same environmental predictors as m1
#but also includes spatial latent structure
#
#this plot allows comparison with m1:
#
#which environmental responses remain similar?
#which weaken after adding latent structure?
#which change direction?
#
#important:
#this plot shows posterior mean beta values only
#it does not show uncertainty around those estimates

#add predictor names to beta matrix
rownames(postBeta_m2$mean) <- pred_names

#define colour scale
#
#we use the same colour scale as the m1 heatmap
#so that the two plots are visually comparable
#
#blue = negative response
#white = weak or near-zero response
#red = positive response
zlim_m2 <- zlim

cols_m2 <- cols



#############################
### Visualise in R window ###
#############################



#show beta heatmap in R plotting window
layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main heatmap
par(mar = c(10, 10, 3, 1))

image(t(postBeta_m2$mean),
      axes = FALSE,
      col = cols_m2,
      zlim = c(-zlim_m2, zlim_m2),
      main = 'Community-wide environmental responses')

#add species names on x axis
axis(1,
     at = seq(0, 1, length.out = ncol(postBeta_m2$mean)),
     labels = colnames(postBeta_m2$mean),
     las = 2,
     cex.axis = 0.7)

#add predictor names on y axis
axis(2,
     at = seq(0, 1, length.out = length(pred_names)),
     labels = rev(pred_names),
     las = 2,
     cex.axis = 0.8)

#colour legend
par(mar = c(10, 2, 3, 4))

image(x = 1,
      y = legend_vals,
      z = matrix(legend_vals, nrow = 1),
      col = cols_m2,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Posterior mean beta',
      side = 4,
      line = 2.8,
      cex = 0.9)



#############################
### Residual associations ###
#############################



#compute residual species associations from m2
#
#Omega describes residual association structure among species
#after accounting for:
#
#measured environmental predictors
#PLUS
#latent spatial/community structure
#
#important:
#
#Omega does NOT automatically mean species interactions
#
#positive residual association:
#two species tend to co-occur more than expected after accounting for predictors
#
#negative residual association:
#two species tend to co-occur less than expected after accounting for predictors
#
#but this may reflect:
#
#missing predictors
#shared habitat preferences
#spatial structure
#sampling effects
#biotic interactions
#or combinations of these

Omega <- computeAssociations(m2)

#inspect object structure
str(Omega)

#inspect available random levels
names(Omega)

#extract posterior mean residual association matrix
Omega_mean <- Omega[[1]]$mean

#extract posterior support for association signs
Omega_support <- Omega[[1]]$support

#inspect dimensions
dim(Omega_mean)

#inspect first rows
round(Omega_mean[1:5, 1:5], 2)
round(Omega_support[1:5, 1:5], 2)



#########################
### Visualise Omega   ###
#########################



#plot posterior mean residual associations
#
#rows and columns = reptile species
#
#positive values:
#species pairs co-occur more than expected after accounting for predictors
#
#negative values:
#species pairs co-occur less than expected after accounting for predictors
#
#values near zero:
#weak residual association
#
#important:
#Omega does NOT automatically represent direct species interactions

#prepare Omega matrix for visualisation
#
#the diagonal is always 1 because each species is perfectly associated
#with itself
#
#Omega is symmetric, so the upper triangle duplicates information
#
#therefore, we remove the diagonal and upper triangle, but keep all
#lower-triangle associations visible
Omega_plot <- Omega_mean

Omega_plot[upper.tri(Omega_plot, diag = FALSE)] <- NA

#define colour scale
omega_zlim <- max(abs(Omega_plot), na.rm = TRUE)

omega_cols <- colorRampPalette(c('blue', 'white', 'red'))(100)

#plot Omega matrix with colour legend
layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main Omega heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_plot[nrow(Omega_plot):1, ]),
      axes = FALSE,
      col = omega_cols,
      zlim = c(-omega_zlim, omega_zlim),
      main = 'Residual community associations')

#add species names on x axis
axis(1,
     at = seq(0, 1, length.out = ncol(Omega_plot)),
     labels = colnames(Omega_plot),
     las = 2,
     cex.axis = 0.7)

#add species names on y axis
axis(2,
     at = seq(0, 1, length.out = nrow(Omega_plot)),
     labels = rev(rownames(Omega_plot)),
     las = 2,
     cex.axis = 0.7)


#colour legend
par(mar = c(6, 1, 2, 4))

omega_legend_vals <- seq(-omega_zlim, omega_zlim, length.out = 100)

image(x = 1,
      y = omega_legend_vals,
      z = matrix(omega_legend_vals, nrow = 1),
      col = omega_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Residual association',
      side = 4,
      line = 2.8,
      cex = 0.9)

#reset layout
layout(1)



#############################
### Clean latent model    ###
#############################



#create alternative latent-variable model
#
#this model excludes:
#
#LizUnIdent
#Tesgra
#
#the goal is to inspect whether:
#
#residual associations
#environmental responses
#latent structure
#and convergence behaviour
#
#change after removing these taxa

#make sure Y_m3 is a matrix
Y_m3 <- as.matrix(Y_m3)

m3 <- Hmsc(Y = Y_m3,
           XData = XScaled,
           XFormula = XFormula,
           distr = 'poisson',
           studyDesign = studyDesign,
           ranLevels = list(sample = rL.site))

#inspect model object
m3



#############################
### Fit clean latent model ###
#############################



#fit the clean latent-variable model
#
#m3 uses the same structure as m2
#but excludes:
#
#LizUnIdent
#Tesgra
#
#we keep the same MCMC settings as m2
#so the models remain comparable

m3 <- sampleMcmc(m3,
                 samples = 1000,
                 thin = 10,
                 transient = 1000,
                 nChains = 4,
                 verbose = 1000)


#setwd
setwd(wd_models)


#save the clean latent-variable model
#
#this model excludes LizUnIdent and Tesgra
saveRDS(m3, 'm3_clean_spatial_latent.rds')


#extract posterior beta estimates for m3
#
#this allows comparison after removing problematic taxa
postBeta_m3 <- getPostEstimate(m3, parName = 'Beta')

#save beta estimates
setwd(wd_models)
saveRDS(postBeta_m3, 'postBeta_m3_clean_spatial_latent.rds')


#convert posterior samples to coda objects
mpost3 <- convertToCodaObject(m3)

#save posterior objects used for diagnostics
#
#these are useful because they allow us to inspect convergence again
#without recreating the coda objects every time

setwd(wd_models)
saveRDS(mpost3, 'mpost_m3_clean_spatial_latent.rds')

#check effective sample size for beta parameters
effectiveSize(mpost3$Beta)

#check Gelman diagnostics for beta parameters
gelman.diag(mpost3$Beta, multivariate = FALSE)



#########################
### Visualise ESS m3  ###
#########################



#extract effective sample sizes for beta parameters from m3
ess_beta_m3 <- effectiveSize(mpost3$Beta)

#visualise ESS distribution for m3
par(mar = c(5, 5, 2, 1))

hist(ess_beta_m3,
     breaks = 20,
     main = 'Effective sample size of beta parameters: m3',
     xlab = 'ESS')



##############################
### Fit longer clean model ###
##############################



#fit longer version of m3
#
#m3_long uses the same structure as m3
#but longer chains, stronger thinning, and longer burn-in
#
#goal:
#test whether remaining convergence issues improve with longer sampling

m3_long <- sampleMcmc(m3,
                      samples = 3000,
                      thin = 20,
                      transient = 3000,
                      nChains = 4,
                      verbose = 1000)

#setwd
setwd(wd_models)

#save longer clean latent-variable model
saveRDS(m3_long, 'm3_long_clean_spatial_latent.rds')

#convert posterior samples to coda objects
mpost3_long <- convertToCodaObject(m3_long)

#save posterior diagnostic object
saveRDS(mpost3_long, 'mpost_m3_long_clean_spatial_latent.rds')

#extract effective sample sizes for beta parameters
ess_beta_m3_long <- effectiveSize(mpost3_long$Beta)

#check Gelman diagnostics for beta parameters
gelman.diag(mpost3_long$Beta, multivariate = FALSE)



################################
### Residual associations m3 ###
################################



#compute residual species associations from m3
#
#this allows comparison of residual covariance structure
#before and after removing potentially problematic taxa
Omega_m3 <- computeAssociations(m3)

#inspect object structure
str(Omega_m3)

#inspect available random levels
names(Omega_m3)

#extract posterior mean residual association matrix
Omega_mean_m3 <- Omega_m3[[1]]$mean

#extract posterior support for association signs
Omega_support_m3 <- Omega_m3[[1]]$support

#inspect dimensions
dim(Omega_mean_m3)

#inspect first rows
round(Omega_mean_m3[1:5, 1:5], 2)
round(Omega_support_m3[1:5, 1:5], 2)



#############################
### Prepare Omega plot m3 ###
#############################



#prepare Omega matrix for visualisation
#
#remove diagonal and upper triangle
Omega_plot_m3 <- Omega_mean_m3

Omega_plot_m3[upper.tri(Omega_plot_m3,
                        diag = FALSE)] <- NA

#define colour scale
omega_zlim_m3 <- max(abs(Omega_plot_m3),
                     na.rm = TRUE)

omega_legend_vals_m3 <- seq(-omega_zlim_m3,
                            omega_zlim_m3,
                            length.out = 100)



#######################################
### Residual associations m3_long   ###
#######################################



#compute residual species associations from m3_long
#
#this tests whether residual covariance structure
#is stable after longer posterior sampling
Omega_m3_long <- computeAssociations(m3_long)

#extract posterior mean residual association matrix
Omega_mean_m3_long <- Omega_m3_long[[1]]$mean

#extract posterior support for association signs
Omega_support_m3_long <- Omega_m3_long[[1]]$support



##################################
### Prepare Omega plot m3_long ###
##################################



#prepare Omega matrix for visualisation
#
#remove diagonal and upper triangle
Omega_plot_m3_long <- Omega_mean_m3_long

Omega_plot_m3_long[upper.tri(Omega_plot_m3_long,
                             diag = FALSE)] <- NA

#define colour scale
omega_zlim_m3_long <- max(abs(Omega_plot_m3_long),
                          na.rm = TRUE)

omega_legend_vals_m3_long <- seq(-omega_zlim_m3_long,
                                 omega_zlim_m3_long,
                                 length.out = 100)



############################
### Summary for meeting  ###
############################



#this section creates simple summaries of what was done today
#the aim is to have clear outputs to discuss with Yoni
#
#today's workflow:
#
#1. loaded reptile community data
#2. created species abundance matrix
#3. removed species occurring in fewer than 10 sites
#4. created environmental predictor matrix
#5. standardised predictors
#6. fitted environment-only HMSC model
#7. fitted spatial latent-variable HMSC model
#8. inspected preliminary convergence diagnostics

#summarise retained species prevalence
species_summary <- data.frame(
  species = colnames(Y),
  occupied_sites = colSums(Y > 0),
  total_abundance = colSums(Y)
)

#save species summary
write.csv(species_summary,
          'species_summary_retained_species.csv',
          row.names = FALSE)

#summarise predictors used in the model
predictor_summary <- data.frame(
  predictor = colnames(XData),
  mean_raw = apply(XData, 2, mean, na.rm = TRUE),
  sd_raw = apply(XData, 2, sd, na.rm = TRUE),
  min_raw = apply(XData, 2, min, na.rm = TRUE),
  max_raw = apply(XData, 2, max, na.rm = TRUE)
)

#save predictor summary
write.csv(predictor_summary,
          'predictor_summary_first_model.csv',
          row.names = FALSE)

#extract beta estimates from both models
postBeta_m1 <- getPostEstimate(m1, parName = 'Beta')
postBeta_m2 <- getPostEstimate(m2, parName = 'Beta')

#manually define predictor names
#these match the rows of the beta matrix
pred_names <- c('Intercept',
                'T_min',
                'Precipitation',
                'Elev_MEAN',
                'BubibisSmpl',
                'Barrenness',
                'Ndvi_MEAN',
                'DISTURB')

#add names to beta matrices
rownames(postBeta_m1$mean) <- pred_names
rownames(postBeta_m2$mean) <- pred_names

#extract egret predation responses from both models
egret_beta_summary <- data.frame(
  species = colnames(Y),
  m1_environment_only = postBeta_m1$mean['BubibisSmpl', ],
  m2_spatial_latent = postBeta_m2$mean['BubibisSmpl', ],
  difference_m2_minus_m1 =
    postBeta_m2$mean['BubibisSmpl', ] -
    postBeta_m1$mean['BubibisSmpl', ]
)

#save egret response summary
write.csv(egret_beta_summary, 'egret_beta_summary_m1_vs_m2.csv', row.names = FALSE)



##################################
### Add values to heatmaps     ###
##################################



add_matrix_values <- function(mat,
                              digits = 2,
                              cex = 0.45,
                              threshold = 0.45){
  
  x_at <- seq(0, 1, length.out = ncol(mat))
  y_at <- seq(0, 1, length.out = nrow(mat))
  
  for(i in seq_len(nrow(mat))){
    
    for(j in seq_len(ncol(mat))){
      
      if(is.finite(mat[i, j])){
        
        text_col <- ifelse(abs(mat[i, j]) >= threshold,
                           'white',
                           'black')
        
        text(x = x_at[j],
             y = y_at[nrow(mat) - i + 1],
             labels = round(mat[i, j], digits),
             cex = cex,
             col = text_col)
      }
    }
  }
}



############################
### Save final figures   ###
############################



#save prevalence plot for discussion slides
setwd(wd_plots)
png('species_prevalence_retained_taxa.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 2, 1))

barplot(species_prev,
        las = 2,
        ylab = 'Occupied sites',
        main = 'Retained community taxa')

abline(h = 10,
       lty = 2)

dev.off()



#############################################
### Save raw species correlation plot      ###
#############################################



#save raw species-species correlation plot
png('raw_species_species_correlations.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main correlation heatmap
par(mar = c(8, 8, 3, 1))

image(t(species_cor_plot[nrow(species_cor_plot):1, ]),
      axes = FALSE,
      col = species_cor_cols,
      zlim = c(-species_cor_zlim,
               species_cor_zlim),
      main = 'Raw species co-occurrence structure')

add_matrix_values(species_cor_plot,
                  digits = 2,
                  cex = 0.4,
                  threshold = 0.4)

axis(1,
     at = seq(0, 1,
              length.out = ncol(species_cor_plot)),
     labels = colnames(species_cor_plot),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1,
              length.out = nrow(species_cor_plot)),
     labels = rev(rownames(species_cor_plot)),
     las = 2,
     cex.axis = 0.7)

#colour legend
par(mar = c(6, 1, 2, 4))

image(x = 1,
      y = species_cor_legend_vals,
      z = matrix(species_cor_legend_vals,
                 nrow = 1),
      col = species_cor_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Spearman correlation',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



########################################################
### Save raw species-environment correlations        ###
########################################################


setwd(wd_plots)
#save raw species-environment association heatmap
png('raw_species_environment_correlations.png',
    width = 2200,
    height = 1600,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main heatmap
par(mar = c(10, 10, 3, 1))

image(t(species_env_cor_plot[nrow(species_env_cor_plot):1, ]),
      axes = FALSE,
      col = species_env_cols,
      zlim = c(-species_env_zlim,
               species_env_zlim),
      main = 'Raw species-environment associations')

add_matrix_values(species_env_cor_plot,
                  digits = 2,
                  cex = 0.45,
                  threshold = 0.4)

#species names
axis(1,
     at = seq(0, 1,
              length.out = ncol(species_env_cor_plot)),
     labels = colnames(species_env_cor_plot),
     las = 2,
     cex.axis = 0.7)

#environmental predictors
axis(2,
     at = seq(0, 1,
              length.out = nrow(species_env_cor_plot)),
     labels = rev(rownames(species_env_cor_plot)),
     las = 2,
     cex.axis = 0.8)

#colour legend
par(mar = c(10, 2, 3, 4))

image(x = 1,
      y = species_env_legend_vals,
      z = matrix(species_env_legend_vals,
                 nrow = 1),
      col = species_env_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Spearman correlation',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



#########################################
### Save beta matrix plots for slides ###
#########################################



#set plot output folder
setwd(wd_plots)

#ensure predictor names are attached to both beta matrices
rownames(postBeta_m1$mean) <- pred_names
rownames(postBeta_m2$mean) <- pred_names

#define shared colour scale for m1 and m2
#
#using the same colour scale makes the two heatmaps directly comparable
beta_zlim <- max(abs(c(postBeta_m1$mean,
                       postBeta_m2$mean)))

beta_cols <- colorRampPalette(c('blue', 'white', 'red'))(100)

beta_legend_vals <- seq(-beta_zlim, beta_zlim, length.out = 100)



#save beta matrix m1 plot for discussion slides
png('beta_matrix_heatmap_m1.png',
    width = 2200,
    height = 1600,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

par(mar = c(10, 10, 3, 1))

beta_plot_m1 <- t(postBeta_m1$mean)

image(t(beta_plot_m1[nrow(beta_plot_m1):1, ]),
      axes = FALSE,
      col = beta_cols,
      zlim = c(-beta_zlim, beta_zlim),
      main = 'Community-wide environmental responses')

add_matrix_values(beta_plot_m1,
                  digits = 2,
                  cex = 0.45,
                  threshold = 1.5)

axis(1,
     at = seq(0, 1, length.out = ncol(beta_plot_m1)),
     labels = colnames(beta_plot_m1),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = nrow(beta_plot_m1)),
     labels = rev(rownames(beta_plot_m1)),
     las = 2,
     cex.axis = 0.8)

par(mar = c(10, 2, 3, 4))

image(x = 1,
      y = beta_legend_vals,
      z = matrix(beta_legend_vals, nrow = 1),
      col = beta_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Posterior mean beta',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



#save beta matrix m2 plot for discussion slides
png('beta_matrix_heatmap_m2.png',
    width = 2200,
    height = 1600,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

par(mar = c(10, 10, 3, 1))

beta_plot_m2 <- t(postBeta_m2$mean)

image(t(beta_plot_m2[nrow(beta_plot_m2):1, ]),
      axes = FALSE,
      col = beta_cols,
      zlim = c(-beta_zlim, beta_zlim),
      main = 'Community-wide environmental responses')

add_matrix_values(beta_plot_m2,
                  digits = 2,
                  cex = 0.45,
                  threshold = 1.5)

axis(1,
     at = seq(0, 1, length.out = ncol(beta_plot_m2)),
     labels = colnames(beta_plot_m2),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = nrow(beta_plot_m2)),
     labels = rev(rownames(beta_plot_m2)),
     las = 2,
     cex.axis = 0.8)

par(mar = c(10, 2, 3, 4))

image(x = 1,
      y = beta_legend_vals,
      z = matrix(beta_legend_vals, nrow = 1),
      col = beta_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Posterior mean beta',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)


#save beta matrix m3 plot for discussion slides
png('beta_matrix_heatmap_m3.png',
    width = 2200,
    height = 1600,
    res = 200)

#ensure predictor names are attached to m3 beta matrix
rownames(postBeta_m3$mean) <- pred_names

#use shared colour scale across m1, m2 and m3
beta_zlim_m1_m2_m3 <- max(abs(c(postBeta_m1$mean,
                                postBeta_m2$mean,
                                postBeta_m3$mean)))

beta_legend_vals_m1_m2_m3 <- seq(-beta_zlim_m1_m2_m3,
                                  beta_zlim_m1_m2_m3,
                                  length.out = 100)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

par(mar = c(10, 10, 3, 1))

beta_plot_m3 <- t(postBeta_m3$mean)

image(t(beta_plot_m3[nrow(beta_plot_m3):1, ]),
      axes = FALSE,
      col = beta_cols,
      zlim = c(-beta_zlim_m1_m2_m3,
               beta_zlim_m1_m2_m3),
      main = 'Community-wide environmental responses')

add_matrix_values(beta_plot_m3,
                  digits = 2,
                  cex = 0.45,
                  threshold = 1.5)

axis(1,
     at = seq(0, 1,
              length.out = ncol(beta_plot_m3)),
     labels = colnames(beta_plot_m3),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1,
              length.out = nrow(beta_plot_m3)),
     labels = rev(rownames(beta_plot_m3)),
     las = 2,
     cex.axis = 0.8)

par(mar = c(10, 2, 3, 4))

image(x = 1,
      y = beta_legend_vals_m1_m2_m3,
      z = matrix(beta_legend_vals_m1_m2_m3, nrow = 1),
      col = beta_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Posterior mean beta',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



############################################
### Save BubibisSmpl plots for slides  ###
############################################



#set shared y-axis limits for m1 and m2
#
#this makes the two barplots visually comparable
#and avoids clipping extreme values
egret_ylim <- range(c(egret_beta,
                      egret_beta_m2))

#add small visual buffer
egret_ylim <- egret_ylim + c(-0.05, 0.05)



#save EgrtPredPress m1 plot for discussion slides
png('egret_beta_species_responses_m1.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 2, 1))

barplot(egret_beta,
        las = 2,
        ylim = egret_ylim,
        ylab = 'Posterior mean beta',
        main = 'Species responses to BubibisSmpl')

abline(h = 0)

dev.off()



#save EgrtPredPress m2 plot for discussion slides
png('egret_beta_species_responses_m2.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 2, 1))

barplot(egret_beta_m2,
        las = 2,
        ylim = egret_ylim,
        ylab = 'Posterior mean beta',
        main = 'Species responses to BubibisSmpl')

abline(h = 0)

dev.off()



#################################################
### Save m1 vs m2 comparison plot for slides  ###
#################################################



#calculate coefficient shifts after adding latent structure
egret_delta <- egret_beta_m2 - egret_beta

#save comparison plot
png('egret_beta_m2_minus_m1.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 3, 1))

barplot(egret_delta,
        las = 2,
        ylab = 'Change in posterior mean beta',
        main = 'Effect of latent structure on BubibisSmpl responses')

abline(h = 0,
       lty = 2)

dev.off()



############################################################
### Save BubibisSmpl plots for m2 and m3 comparison    ###
############################################################



#species retained in m3
species_m3 <- colnames(postBeta_m3$mean)

#extract BubibisSmpl responses from m2,
#restricted to species retained in m3
egret_beta_m2_restricted <- postBeta_m2$mean['BubibisSmpl',
                                             species_m3]

#extract BubibisSmpls responses from m3
egret_beta_m3 <- postBeta_m3$mean['BubibisSmpl', ]

#calculate effect of taxon removal
egret_delta_m3_minus_m2 <- egret_beta_m3 - egret_beta_m2_restricted

#shared y-axis limits for m2 restricted and m3
egret_ylim_m2_m3 <- range(c(egret_beta_m2_restricted,
                            egret_beta_m3))

egret_ylim_m2_m3 <- egret_ylim_m2_m3 + c(-0.05, 0.05)



#save EgrtPredPress m2 restricted plot
png('egret_beta_species_responses_m2_restricted.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 2, 1))

barplot(egret_beta_m2_restricted,
        las = 2,
        ylim = egret_ylim_m2_m3,
        ylab = 'Posterior mean beta',
        main = 'Species responses to BubibisSmpl')

abline(h = 0)

dev.off()



#save BubibisSmpl m3 plot
png('egret_beta_species_responses_m3.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 2, 1))

barplot(egret_beta_m3,
        las = 2,
        ylim = egret_ylim_m2_m3,
        ylab = 'Posterior mean beta',
        main = 'Species responses to BubibisSmpl')

abline(h = 0)

dev.off()



#save m3 minus m2 comparison plot
png('egret_beta_m3_minus_m2.png',
    width = 1800,
    height = 1200,
    res = 200)

par(mar = c(8, 5, 3, 1))

barplot(egret_delta_m3_minus_m2,
        las = 2,
        ylab = 'Change in posterior mean beta',
        main = 'Effect of taxon removal on BubibisSmpl responses')

abline(h = 0,
       lty = 2)

dev.off()



######################################
### Save Omega m2 plot for slides  ###
######################################



#save residual association plot for discussion slides
png('omega_residual_species_associations.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main Omega heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_plot[nrow(Omega_plot):1, ]),
      axes = FALSE,
      col = omega_cols,
      zlim = c(-omega_zlim, omega_zlim),
      main = 'Residual species associations')

add_matrix_values(Omega_plot,
                  digits = 2,
                  cex = 0.4,
                  threshold = 0.5)

axis(1,
     at = seq(0, 1, length.out = ncol(Omega_plot)),
     labels = colnames(Omega_plot),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = nrow(Omega_plot)),
     labels = rev(rownames(Omega_plot)),
     las = 2,
     cex.axis = 0.7)

#colour legend
par(mar = c(6, 1, 2, 4))

image(x = 1,
      y = omega_legend_vals,
      z = matrix(omega_legend_vals, nrow = 1),
      col = omega_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Residual association',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



######################################
### Save Omega m3 plot for slides  ###
######################################



#save residual association plot for m3
png('omega_residual_species_associations_m3.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main Omega heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_plot_m3[nrow(Omega_plot_m3):1, ]),
      axes = FALSE,
      col = omega_cols,
      zlim = c(-omega_zlim_m3, omega_zlim_m3),
      main = 'Residual species associations')

add_matrix_values(Omega_plot_m3,
                  digits = 2,
                  cex = 0.4,
                  threshold = 0.5)

axis(1,
     at = seq(0, 1, length.out = ncol(Omega_plot_m3)),
     labels = colnames(Omega_plot_m3),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = nrow(Omega_plot_m3)),
     labels = rev(rownames(Omega_plot_m3)),
     las = 2,
     cex.axis = 0.7)

#colour legend
par(mar = c(6, 1, 2, 4))

image(x = 1,
      y = omega_legend_vals_m3,
      z = matrix(omega_legend_vals_m3, nrow = 1),
      col = omega_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Residual association',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



###########################################
### Save Omega m3_long plot for slides  ###
###########################################



#save residual association plot for m3_long
png('omega_residual_species_associations_m3_long.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main Omega heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_plot_m3_long[nrow(Omega_plot_m3_long):1, ]),
      axes = FALSE,
      col = omega_cols,
      zlim = c(-omega_zlim_m3_long,
               omega_zlim_m3_long),
      main = 'Residual species associations')

add_matrix_values(Omega_plot_m3_long,
                  digits = 2,
                  cex = 0.4,
                  threshold = 0.5)

axis(1,
     at = seq(0, 1, length.out = ncol(Omega_plot_m3_long)),
     labels = colnames(Omega_plot_m3_long),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = nrow(Omega_plot_m3_long)),
     labels = rev(rownames(Omega_plot_m3_long)),
     las = 2,
     cex.axis = 0.7)

#colour legend
par(mar = c(6, 1, 2, 4))

image(x = 1,
      y = omega_legend_vals_m3_long,
      z = matrix(omega_legend_vals_m3_long, nrow = 1),
      col = omega_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Residual association',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



#############################################
### Difference in residual associations   ###
#############################################



#species retained in m3
species_m3 <- colnames(Omega_mean_m3)

#restrict m2 residual association matrix to species retained in m3
Omega_mean_m2_restricted <- Omega_mean[species_m3,
                                       species_m3]

#calculate change in residual associations after taxon removal
#
#positive values:
#residual association is stronger in m3 than in m2
#
#negative values:
#residual association is weaker in m3 than in m2
Omega_diff_m3_minus_m2 <- Omega_mean_m3 - Omega_mean_m2_restricted

#prepare difference matrix for visualisation
#
#remove diagonal and upper triangle
Omega_diff_plot <- Omega_diff_m3_minus_m2

Omega_diff_plot[upper.tri(Omega_diff_plot,
                          diag = FALSE)] <- NA

#define colour scale for difference plot
omega_diff_zlim <- max(abs(Omega_diff_plot),
                       na.rm = TRUE)

omega_diff_cols <- colorRampPalette(c('orange',
                                      'white',
                                      'darkgreen'))(100)

omega_diff_legend_vals <- seq(-omega_diff_zlim,
                              omega_diff_zlim,
                              length.out = 100)



#########################################################
### Save Omega difference plot for slides             ###
#########################################################



#save difference in residual associations between m3 and m2
png('omega_difference_m3_minus_m2.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main difference heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_diff_plot[nrow(Omega_diff_plot):1, ]),
      axes = FALSE,
      col = omega_diff_cols,
      zlim = c(-omega_diff_zlim,
               omega_diff_zlim),
      main = 'Change in residual species associations')

add_matrix_values(Omega_diff_plot,
                  digits = 2,
                  cex = 0.4,
                  threshold = 0.3)

axis(1,
     at = seq(0, 1,
              length.out = ncol(Omega_diff_plot)),
     labels = colnames(Omega_diff_plot),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1,
              length.out = nrow(Omega_diff_plot)),
     labels = rev(rownames(Omega_diff_plot)),
     las = 2,
     cex.axis = 0.7)

#colour legend
par(mar = c(6, 1, 2, 4))

image(x = 1,
      y = omega_diff_legend_vals,
      z = matrix(omega_diff_legend_vals,
                 nrow = 1),
      col = omega_diff_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Change in residual association',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



###########################################################
### Difference in residual associations: m3_long vs m3 ###
###########################################################



#calculate change in residual associations after longer sampling
#
#positive values:
#residual association is stronger in m3_long than in m3
#
#negative values:
#residual association is weaker in m3_long than in m3
Omega_diff_m3_long_minus_m3 <- Omega_mean_m3_long - Omega_mean_m3

#prepare difference matrix for visualisation
#
#remove diagonal and upper triangle
Omega_diff_plot_m3_long <- Omega_diff_m3_long_minus_m3

Omega_diff_plot_m3_long[upper.tri(Omega_diff_plot_m3_long,
                                  diag = FALSE)] <- NA

#define colour scale for difference plot
omega_diff_zlim_m3_long <- max(abs(Omega_diff_plot_m3_long),
                               na.rm = TRUE)

omega_diff_legend_vals_m3_long <- seq(-omega_diff_zlim_m3_long,
                                      omega_diff_zlim_m3_long,
                                      length.out = 100)


########################################################
### Save Omega difference plot for m3_long minus m3  ###
########################################################


 
#save difference in residual associations between m3_long and m3
png('omega_difference_m3_long_minus_m3.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main difference heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_diff_plot_m3_long[nrow(Omega_diff_plot_m3_long):1, ]),
      axes = FALSE,
      col = omega_diff_cols,
      zlim = c(-omega_diff_zlim_m3_long,
               omega_diff_zlim_m3_long),
      main = 'Change in residual species associations')

add_matrix_values(Omega_diff_plot_m3_long,
                  digits = 2,
                  cex = 0.4,
                  threshold = 0.3)

axis(1,
     at = seq(0, 1,
              length.out = ncol(Omega_diff_plot_m3_long)),
     labels = colnames(Omega_diff_plot_m3_long),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1,
              length.out = nrow(Omega_diff_plot_m3_long)),
     labels = rev(rownames(Omega_diff_plot_m3_long)),
     las = 2,
     cex.axis = 0.7)

#colour legend
par(mar = c(6, 1, 2, 4))

image(x = 1,
      y = omega_diff_legend_vals_m3_long,
      z = matrix(omega_diff_legend_vals_m3_long,
                 nrow = 1),
      col = omega_diff_cols,
      axes = FALSE,
      xlab = '',
      ylab = '')

axis(4,
     las = 2,
     cex.axis = 0.8)

mtext('Change in residual association',
      side = 4,
      line = 2.8,
      cex = 0.9)

dev.off()

layout(1)



###################################
### Save ESS plots for slides   ###
###################################



#save ESS comparison plot for m1 and m2
setwd(wd_plots)

png('ess_beta_m1_m2_histograms.png',
    width = 2200,
    height = 1200,
    res = 200)

par(mfrow = c(1, 2),
    mar = c(5, 5, 3, 1),
    oma = c(0, 0, 0, 0))

hist(ess_beta_m1,
     breaks = seq(0, 650, by = 25),
     xlim = c(0, 650),
     ylim = c(0, 45),
     yaxs = 'i',
     main = 'm1',
     xlab = 'ESS',
     ylab = 'Frequency')

hist(ess_beta_m2,
     breaks = seq(0, 650, by = 25),
     xlim = c(0, 650),
     ylim = c(0, 45),
     yaxs = 'i',
     main = 'm2',
     xlab = 'ESS',
     ylab = 'Frequency')

dev.off()

#reset plotting layout
par(mfrow = c(1, 1))



#####################################
### Save ESS plot for m2 and m3   ###
#####################################



#save ESS comparison plot for m2 and m3
setwd(wd_plots)

png('ess_beta_m2_m3_histograms.png',
    width = 2200,
    height = 1200,
    res = 200)

par(mfrow = c(1, 2),
    mar = c(5, 5, 3, 1),
    oma = c(0, 0, 0, 0))

hist(ess_beta_m2,
     breaks = seq(0, 300, by = 15),
     xlim = c(0, 300),
     ylim = c(0, 45),
     yaxs = 'i',
     main = 'm2',
     xlab = 'ESS',
     ylab = 'Frequency')

hist(ess_beta_m3,
     breaks = seq(0, 300, by = 15),
     xlim = c(0, 300),
     ylim = c(0, 45),
     yaxs = 'i',
     main = 'm3',
     xlab = 'ESS',
     ylab = 'Frequency')

dev.off()

#reset plotting layout
par(mfrow = c(1, 1))



##########################################
### Save ESS plot for m3 and m3_long   ###
##########################################



#save ESS comparison plot for m3 and m3_long
setwd(wd_plots)

png('ess_beta_m3_m3long_histograms.png',
    width = 2200,
    height = 1200,
    res = 200)

#shared x-axis range
ess_xlim_long <- c(0, 1500)

#shared histogram breaks
ess_breaks_long <- seq(0, 1500, by = 50)

#calculate shared y-axis limit
hist_m3 <- hist(ess_beta_m3,
                breaks = ess_breaks_long,
                plot = FALSE)

hist_m3_long <- hist(ess_beta_m3_long,
                     breaks = ess_breaks_long,
                     plot = FALSE)

ess_ylim_long <- c(0,
                   max(c(hist_m3$counts,
                         hist_m3_long$counts)))

par(mfrow = c(1, 2),
    mar = c(5, 5, 3, 1),
    oma = c(0, 0, 0, 0))

hist(ess_beta_m3,
     breaks = ess_breaks_long,
     xlim = ess_xlim_long,
     ylim = ess_ylim_long,
     yaxs = 'i',
     main = 'm3',
     xlab = 'ESS',
     ylab = 'Frequency')

hist(ess_beta_m3_long,
     breaks = ess_breaks_long,
     xlim = ess_xlim_long,
     ylim = ess_ylim_long,
     yaxs = 'i',
     main = 'm3_long',
     xlab = 'ESS',
     ylab = 'Frequency')

dev.off()

#reset plotting layout
par(mfrow = c(1, 1))



###################################
### Save trace plots for slides ###
###################################



#select parameter to inspect
#
#this is the BubibisSmpl response for Eumsch
trace_param <- 'B[BubibisSmpl (C5), Eumsch (S4)]'

#check that the parameter exists in all models
trace_param %in% colnames(as.matrix(mpost$Beta[[1]]))
trace_param %in% colnames(as.matrix(mpost2$Beta[[1]]))
trace_param %in% colnames(as.matrix(mpost3$Beta[[1]]))

#extract chains from m1
trace_m1 <- as.matrix(mpost$Beta)[, trace_param]

#extract chains from m2
trace_m2 <- as.matrix(mpost2$Beta)[, trace_param]

#extract chains from m3
trace_m3 <- as.matrix(mpost3$Beta)[, trace_param]


#shared y-axis limits for m1, m2 and m3
trace_ylim <- range(c(trace_m1,
                      trace_m2,
                      trace_m3))

#save trace plot comparison for discussion slides
setwd(wd_plots)

png('trace_beta_egret_eumsch_m1_m2_m3.png',
    width = 2200,
    height = 1500,
    res = 200)

par(mfrow = c(3, 1),
    mar = c(3, 5, 2, 1))

#trace plot m1
plot(trace_m1,
     type = 'l',
     ylim = trace_ylim,
     xlab = '',
     ylab = 'Posterior beta',
     main = 'm1',
     xaxt = 'n')

#trace plot m2
plot(trace_m2,
     type = 'l',
     ylim = trace_ylim,
     xlab = '',
     ylab = 'Posterior beta',
     main = 'm2',
     xaxt = 'n')

#trace plot m3
plot(trace_m3,
     type = 'l',
     ylim = trace_ylim,
     xlab = 'MCMC iteration',
     ylab = 'Posterior beta',
     main = 'm3')

dev.off()

#reset layout
par(mfrow = c(1, 1))



############################################
### Save trace plots for m3 and m3_long  ###
############################################



#check that the parameter exists in m3_long
trace_param %in% colnames(as.matrix(mpost3_long$Beta[[1]]))

#extract chain from m3_long
trace_m3_long <- as.matrix(mpost3_long$Beta)[, trace_param]

#shared y-axis limits for m3 and m3_long
trace_ylim_m3_long <- range(c(trace_m3,
                              trace_m3_long))

#shared x-axis limits
trace_xlim_m3_long <- c(1,
                        length(trace_m3_long))

#save trace plot comparison for m3 and m3_long
setwd(wd_plots)

png('trace_beta_bubibis_eumsch_m3_m3long.png',
    width = 2200,
    height = 1200,
    res = 200)

par(mfrow = c(2, 1),
    mar = c(3, 5, 2, 1))

#trace plot m3
plot(seq_along(trace_m3),
     trace_m3,
     type = 'l',
     xlim = trace_xlim_m3_long,
     ylim = trace_ylim_m3_long,
     xlab = '',
     ylab = 'Posterior beta',
     main = 'm3',
     xaxt = 'n')

#trace plot m3_long
plot(seq_along(trace_m3_long),
     trace_m3_long,
     type = 'l',
     xlim = trace_xlim_m3_long,
     ylim = trace_ylim_m3_long,
     xlab = 'MCMC iteration',
     ylab = 'Posterior beta',
     main = 'm3_long')

dev.off()

#reset layout
par(mfrow = c(1, 1))



###############################################
### Prepare Gelman diagnostics for plotting ###
###############################################



#calculate Gelman diagnostics for beta parameters
gelman_m1 <- gelman.diag(mpost$Beta,
                         multivariate = FALSE)

gelman_m2 <- gelman.diag(mpost2$Beta,
                         multivariate = FALSE)

gelman_m3 <- gelman.diag(mpost3$Beta,
                         multivariate = FALSE)

gelman_m3_long <- gelman.diag(mpost3_long$Beta,
                              multivariate = FALSE)

#extract point estimates
gelman_beta_m1 <- gelman_m1$psrf[, 'Point est.']
gelman_beta_m2 <- gelman_m2$psrf[, 'Point est.']
gelman_beta_m3 <- gelman_m3$psrf[, 'Point est.']
gelman_beta_m3_long <- gelman_m3_long$psrf[, 'Point est.']



###########################################
### Save Gelman plots for slides        ###
###########################################



#shared x-axis limits
gelman_xlim <- c(0.95, 5)

#shared histogram breaks
gelman_breaks <- seq(0.95, 5, by = 0.05)

#shared y-axis limit
gelman_hist_m1 <- hist(gelman_beta_m1,
                       breaks = gelman_breaks,
                       plot = FALSE)

gelman_hist_m2 <- hist(gelman_beta_m2,
                       breaks = gelman_breaks,
                       plot = FALSE)

gelman_hist_m3 <- hist(gelman_beta_m3,
                       breaks = gelman_breaks,
                       plot = FALSE)

gelman_hist_m3_long <- hist(gelman_beta_m3_long,
                            breaks = gelman_breaks,
                            plot = FALSE)

gelman_ylim <- c(0,
                 max(c(gelman_hist_m1$counts,
                       gelman_hist_m2$counts,
                       gelman_hist_m3$counts,
                       gelman_hist_m3_long$counts)))



#save plot
setwd(wd_plots)

png('gelman_beta_all_models.png',
    width = 2400,
    height = 1800,
    res = 200)

par(mfrow = c(2, 2),
    mar = c(5, 5, 3, 1))

#m1
hist(gelman_beta_m1,
     breaks = gelman_breaks,
     xlim = gelman_xlim,
     ylim = gelman_ylim,
     yaxs = 'i',
     main = 'm1',
     xlab = 'Gelman-Rubin PSRF',
     ylab = 'Frequency')

abline(v = 1.1,
       lty = 2)

#m2
hist(gelman_beta_m2,
     breaks = gelman_breaks,
     xlim = gelman_xlim,
     ylim = gelman_ylim,
     yaxs = 'i',
     main = 'm2',
     xlab = 'Gelman-Rubin PSRF',
     ylab = 'Frequency')

abline(v = 1.1,
       lty = 2)

#m3
hist(gelman_beta_m3,
     breaks = gelman_breaks,
     xlim = gelman_xlim,
     ylim = gelman_ylim,
     yaxs = 'i',
     main = 'm3',
     xlab = 'Gelman-Rubin PSRF',
     ylab = 'Frequency')

abline(v = 1.1,
       lty = 2)

#m3_long
hist(gelman_beta_m3_long,
     breaks = gelman_breaks,
     xlim = gelman_xlim,
     ylim = gelman_ylim,
     yaxs = 'i',
     main = 'm3_long',
     xlab = 'Gelman-Rubin PSRF',
     ylab = 'Frequency')

abline(v = 1.1,
       lty = 2)

dev.off()

#reset layout
par(mfrow = c(1, 1))



#############################
### Notes for discussion  ###
#############################



#points to discuss with Avi:
#
#1. Is a Poisson model appropriate for these count data,
#or should we later consider overdispersion or another distribution?
#
#2. Is the 10 occupied-site threshold reasonable for the first model?
#
#3. Should LizUnIdent be kept, given that it is not a resolved species?
#
#4. Which predictors are ecologically justified in the first model?
#
#5. Is EgrtPredPress best treated as a predictor of reptile abundance,
#or should it be interpreted more cautiously as an index correlated with
#habitat and landscape structure?
#
#6. Does adding spatial latent structure change the estimated responses
#to EgrtPredPress?
#
#7. What would be the clearest biological question for this dataset:
#
#community response to egret abundance/activity?
#species-specific vulnerability to egret predation?
#residual associations after accounting for predator pressure?
#comparison between observed egret diet and reptile community composition?







