############################
###  Load reptile data   ###
############################


#clean workspace
rm(list = ls())

#load packages
library(Hmsc)

#list WDs
wd_data <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Data'
wd_models <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_models'
wd_plots <- '/Users/carloseduardoaribeiro/Documents/Post-doc_2/Egret_reptile_models/Plots'

#load main reptile table
setwd(wd_data)
reptiles <- read.csv('NorthernIsraelReptiles.csv')

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



#create species matrix
Y <- reptiles[, 24:47]

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



#calculate the number of occupied sites per species
#this counts how many sites have abundance greater than zero
#species with very few occupied sites provide little information for estimating
#environmental responses and species associations
occ_sites <- colSums(Y > 0)

#inspect prevalence before filtering
#this lets us see which species are common, intermediate, or very rare
occ_sites

#keep only species occurring in at least 10 sites
#this is a conservative first threshold to make the first model more stable
#rare species can be added back later, but for now we prioritise learning the workflow
Y <- Y[, occ_sites >= 10]

#check which species remain after filtering
names(Y)

#check dimensions of the filtered community matrix
#rows should still be sites, columns are now the retained species
dim(Y)

#check prevalence again after filtering
colSums(Y > 0)



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
        main = 'Retained reptile taxa')

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
                      'EgrtPredPress',
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
#EgrtPredPress:
#our biologically most interesting variable
#represents estimated cattle egret predation pressure
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
### Standardise XData     ###
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
#egret predation pressure
#barrenness
#NDVI
#disturbance
#
#at this stage, we are fitting only environmental responses
#we are not yet adding latent variables or spatial random effects
XFormula <- ~ T_min +
  Precipitation +
  Elev_MEAN +
  EgrtPredPress +
  Barrenness +
  Ndvi_MEAN +
  DISTURB

#create HMSC model object
#
#distr = 'poisson' because the species data are counts
#
#conceptually, this model asks:
#
#"how does each reptile species respond to the measured environmental
#predictors, including egret predation pressure?"
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



############################
### Fit first HMSC model  ###
############################



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

mpost <- convertToCodaObject(m1)

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
#columns = reptile species
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
                'EgrtPredPress',
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
par(mar = c(10, 10, 3, 1))

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
par(mar = c(10, 2, 3, 4))

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
#in their apparent responses to cattle egret predation pressure



##############################
### Egret predation effect ###
##############################


#extract posterior mean responses to egret predation pressure
#
#row 5 corresponds to EgrtPredPress because:
#
#1 = intercept
#2 = T_min
#3 = Precipitation
#4 = Elev_MEAN
#5 = EgrtPredPress
#
#positive beta:
#species tends to occur more where egret predation pressure is higher
#
#negative beta:
#species tends to occur less where egret predation pressure is higher

egret_beta <- postBeta$mean[5, ]

#inspect values
egret_beta

#plot species responses to egret predation pressure
par(mar = c(8, 5, 2, 1))

barplot(egret_beta,
        las = 2,
        ylab = 'Posterior mean beta',
        main = 'Species responses to egret predation pressure')

#abline at zero
#
#species above zero:
#positive association with egret predation pressure
#
#species below zero:
#negative association with egret predation pressure
abline(h = 0)


#important conceptual point:
#
#these are CONDITIONAL environmental responses
#
#the model is estimating species responses to egret predation pressure
#while simultaneously accounting for the other predictors in the model
#
#therefore, this is not simply a raw correlation
#
#however, we still must be very careful with interpretation
#
#a negative response does NOT automatically prove direct predation effects
#
#the pattern could also emerge because:
#
#egrets prefer certain habitats
#egret pressure correlates with disturbance
#egret pressure correlates with vegetation openness
#some reptile species avoid habitats frequented by egrets
#
#therefore:
#
#these coefficients represent statistical associations conditional on
#the included predictors
#
#not definitive causal ecological mechanisms



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
#egret predation pressure
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

#convert posterior samples to coda objects
#
#this lets us inspect MCMC diagnostics
mpost2 <- convertToCodaObject(m2)

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
par(mar = c(5, 5, 2, 1))

hist(ess_beta_m2,
     breaks = 20,
     main = 'Effective sample size of beta parameters: m2',
     xlab = 'ESS')



###############################
### Visualise trace plots   ###
###############################



#select parameter to inspect
trace_param <- 'B[EgrtPredPress (C5), Eumsch (S4)]'

#extract chains from m1
trace_m1 <- as.matrix(mpost$Beta)[, trace_param]

#extract chains from m2
trace_m2 <- as.matrix(mpost2$Beta)[, trace_param]



#############################
### Visualise in R window ###
#############################



#plot trace behaviour for m1 and m2
par(mfrow = c(2, 1),
    mar = c(4, 5, 2, 1))

#trace plot m1
plot(trace_m1,
     type = 'l',
     xlab = 'MCMC iteration',
     ylab = 'Posterior beta',
     main = 'm1')

#trace plot m2
plot(trace_m2,
     type = 'l',
     xlab = 'MCMC iteration',
     ylab = 'Posterior beta',
     main = 'm2')

#reset layout
par(mfrow = c(1, 1))



############################
### Save fitted models   ###
############################



#setwd
setwd(wd_models)

#save the environment-only model
#
#this avoids having to rerun the MCMC every time we reopen R
#this is especially important because HMSC models can take a long time to fit
saveRDS(m1, 'm1_environment_only.rds')

#save the latent-variable model
#
#this model is slower because it includes spatial latent structure
#saving it allows us to inspect diagnostics, coefficients, and associations later
#without repeating the expensive sampling step
saveRDS(m2, 'm2_spatial_latent.rds')



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
#possibly fewer predictors
#possibly additional filtering of rare species
#model comparison
#posterior uncertainty around coefficients and associations



#########################
### Save diagnostics  ###
#########################



#save posterior objects used for diagnostics
#
#these are useful because they allow us to inspect convergence again
#without recreating the coda objects every time
saveRDS(mpost, 'mpost_m1_environment_only.rds')
saveRDS(mpost2, 'mpost_m2_spatial_latent.rds')



############################
### Extract trace chains ###
############################



#inspect available beta parameter names
#
#this helps us choose parameters for visualising chain behaviour
colnames(as.matrix(mpost$Beta[[1]]))

colnames(as.matrix(mpost2$Beta[[1]]))



############################
### Save beta estimates  ###
############################



#extract posterior beta estimates for m1
#
#beta coefficients describe species responses to environmental predictors
postBeta_m1 <- getPostEstimate(m1, parName = 'Beta')

#extract posterior beta estimates for m2
#
#this allows us to compare environmental responses before and after
#adding latent spatial structure
postBeta_m2 <- getPostEstimate(m2, parName = 'Beta')

#save beta estimates
saveRDS(postBeta_m1, 'postBeta_m1_environment_only.rds')
saveRDS(postBeta_m2, 'postBeta_m2_spatial_latent.rds')



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



#################################
### Egret predation effect m2 ###
#################################



#extract posterior mean responses to egret predation pressure
#from the latent-variable model (m2)
#
#this allows direct comparison with m1
#
#positive beta:
#species tends to occur more where egret predation pressure is higher
#
#negative beta:
#species tends to occur less where egret predation pressure is higher
#
#importantly:
#
#these responses are now estimated AFTER accounting for:
#
#environmental predictors
#PLUS
#latent spatial/community structure

egret_beta_m2 <- postBeta_m2$mean[5, ]

#inspect values
egret_beta_m2

#plot species responses to egret predation pressure
par(mar = c(8, 5, 2, 1))

barplot(egret_beta_m2,
        las = 2,
        ylim = c(-0.6, 0.2),
        ylab = 'Posterior mean beta',
        main = 'Species responses to EgrtPredPress')

#abline at zero
abline(h = 0)


#important conceptual point:
#
#these responses are estimated conditionally on:
#
#measured environmental predictors
#PLUS
#latent spatial/community structure
#
#therefore, differences between m1 and m2 may indicate that:
#
#some apparent environmental responses in m1
#were partly associated with residual spatial structure
#
#this is one of the key reasons latent-variable jSDMs
#can substantially change coefficient interpretation



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

Omega_plot[upper.tri(Omega_plot, diag = TRUE)] <- NA

#define colour scale
omega_zlim <- max(abs(Omega_plot), na.rm = TRUE)

omega_cols <- colorRampPalette(c('blue', 'white', 'red'))(100)

#plot Omega matrix with colour legend
layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main Omega heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_plot),
      axes = FALSE,
      col = omega_cols,
      zlim = c(-omega_zlim, omega_zlim),
      main = 'Residual species associations')

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
                'EgrtPredPress',
                'Barrenness',
                'Ndvi_MEAN',
                'DISTURB')

#add names to beta matrices
rownames(postBeta_m1$mean) <- pred_names
rownames(postBeta_m2$mean) <- pred_names

#extract egret predation responses from both models
egret_beta_summary <- data.frame(
  species = colnames(Y),
  m1_environment_only = postBeta_m1$mean['EgrtPredPress', ],
  m2_spatial_latent = postBeta_m2$mean['EgrtPredPress', ],
  difference_m2_minus_m1 =
    postBeta_m2$mean['EgrtPredPress', ] -
    postBeta_m1$mean['EgrtPredPress', ]
)

#save egret response summary
write.csv(egret_beta_summary, 'egret_beta_summary_m1_vs_m2.csv', row.names = FALSE)



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
        main = 'Retained reptile taxa')

abline(h = 10,
       lty = 2)

dev.off()



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
png('beta_matrix_heatmap.png',
    width = 2200,
    height = 1600,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

par(mar = c(10, 10, 3, 1))

image(t(postBeta_m1$mean),
      axes = FALSE,
      col = beta_cols,
      zlim = c(-beta_zlim, beta_zlim),
      main = 'Community-wide environmental responses')

axis(1,
     at = seq(0, 1, length.out = ncol(postBeta_m1$mean)),
     labels = colnames(postBeta_m1$mean),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = length(pred_names)),
     labels = rev(pred_names),
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

image(t(postBeta_m2$mean),
      axes = FALSE,
      col = beta_cols,
      zlim = c(-beta_zlim, beta_zlim),
      main = 'Community-wide environmental responses')

axis(1,
     at = seq(0, 1, length.out = ncol(postBeta_m2$mean)),
     labels = colnames(postBeta_m2$mean),
     las = 2,
     cex.axis = 0.7)

axis(2,
     at = seq(0, 1, length.out = length(pred_names)),
     labels = rev(pred_names),
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



############################################
### Save EgrtPredPress plots for slides  ###
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
        main = 'Species responses to EgrtPredPress')

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
        main = 'Species responses to EgrtPredPress')

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
        main = 'Effect of latent structure on EgrtPredPress responses')

abline(h = 0,
       lty = 2)

dev.off()



###################################
### Save Omega plot for slides  ###
###################################



#save residual association plot for discussion slides
png('omega_residual_species_associations.png',
    width = 2200,
    height = 1400,
    res = 200)

layout(matrix(c(1, 2), nrow = 1),
       widths = c(5, 1))

#main Omega heatmap
par(mar = c(8, 8, 3, 1))

image(t(Omega_plot),
      axes = FALSE,
      col = omega_cols,
      zlim = c(-omega_zlim, omega_zlim),
      main = 'Residual species associations')

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



###################################
###  Save ESS plot for slides   ###
###################################



#shared x-axis limits for m1 and m2
ess_xlim <- range(c(ess_beta_m1,
                    ess_beta_m2))

#save ESS comparison plot for discussion slides
setwd(wd_plots)

png('ess_beta_m1_m2_histograms.png',
    width = 2200,
    height = 1200,
    res = 200)

par(mfrow = c(1, 2),
    mar = c(5, 5, 3, 1))

hist(ess_beta_m1,
     breaks = 20,
     xlim = ess_xlim,
     main = 'm1',
     xlab = 'ESS')

hist(ess_beta_m2,
     breaks = 20,
     xlim = ess_xlim,
     main = 'm2',
     xlab = 'ESS')

dev.off()

#reset plotting layout
par(mfrow = c(1, 1))



###################################
### Save trace plots for slides ###
###################################



#save trace plot comparison for discussion slides
png('trace_beta_egret_eumsch_m1_m2.png',
    width = 2200,
    height = 1200,
    res = 200)

#shared y-axis limits for m1 and m2
trace_ylim <- range(c(trace_m1,
                      trace_m2))

par(mfrow = c(2, 1),
    mar = c(5, 5, 2, 1))

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
     xlab = 'MCMC iteration',
     ylab = 'Posterior beta',
     main = 'm2')

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
#community response to egret predation pressure?
#species-specific vulnerability to egret predation?
#residual associations after accounting for predator pressure?
#comparison between observed egret diet and reptile community composition?







