#load packages
library(Hmsc); library(corrplot)


## Explore available HMSC example data ##

#list datasets available in the Hmsc package
data(package = "Hmsc")


## Load example dataset ##

#load simulated teaching dataset included in Hmsc
data(TD)

#list objects currently loaded in the workspace
ls()

#inspect the structure of the TD object
str(TD, max.level = 2)

#TD is a list containing simulated community data, covariates,
#traits, phylogeny, random levels, and a fitted Hmsc model


## Extract fitted Hmsc model ##

#extract the fitted model stored inside TD
m = TD$m

#this model has already been fitted, so we can inspect outputs directly


## Convert posterior samples for MCMC diagnostics ##

#convert posterior samples to coda format
mpost = convertToCodaObject(m)

#coda format is used for diagnostics such as effective sample size
#and Gelman-Rubin convergence statistics


## Check MCMC convergence for Beta parameters ##

#effective sample size of beta parameters
effectiveSize(mpost$Beta)

#larger values indicate better MCMC mixing

## Interpretation for this model ##

#most values here are below 50, meaning posterior sampling is weak
#this happens because the example model included in TD$m was fitted
#with very short MCMC chains for speed in the tutorial

#therefore Beta estimates and plots are still useful for learning
#the workflow, but not reliable enough for scientific conclusions

#in real analyses, we would increase the number of posterior samples
#to obtain effective sample sizes above ~100–200


#Gelman-Rubin diagnostic for beta parameters
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

#psrf values close to 1 indicate good agreement between chains

#rule of thumb:
#~1.00 excellent convergence
#<1.05 very good
#<1.10 acceptable
#>1.20 poor convergence


## Interpretation for this model ##

#several Beta parameters show psrf values above 1.2,
#indicating that chains did not converge well

#this is expected because the example model TD$m was fitted
#with very short MCMC sampling for demonstration purposes

#therefore posterior estimates are still useful for learning
#the workflow, but not reliable for formal inference

#in real analyses we would increase the number of iterations
#until psrf values approach 1


## Evaluate explanatory power ##

#compute fitted predictions for the same data used to fit the model
preds = computePredictedValues(m)

#evaluate model fit
evaluateModelFit(hM = m, predY = preds)

#RMSE measures prediction error (smaller values indicate better fit)

#AUC measures how well the model separates presences from absences
#values close to 1 indicate excellent discrimination ability

#TjurR2 measures the difference between predicted probabilities
#at presences versus absences (larger values indicate better fit)


## Interpretation for this model ##

#AUC values close to 1 indicate that the model almost perfectly
#distinguishes presence from absence for all species

#TjurR2 values between ~0.68 and ~0.82 indicate strong explanatory power

#RMSE values are small, indicating low prediction error

#these metrics represent explanatory power because predictions
#were evaluated on the same data used to fit the model

#predictive power would require cross-validation


## Extract posterior estimates of environmental responses ##

#retrieve posterior summaries of Beta parameters
postBeta = getPostEstimate(m, parName = "Beta")

#postBeta$mean shows average effect size of each predictor on each species

#postBeta$support shows posterior probability that the effect is positive

#postBeta$supportNeg shows posterior probability that the effect is negative


## Interpretation for this model ##

#predictor x1 shows strong positive support (support = 1.00)
#for all species, indicating that x1 increases occurrence probability

#predictor x2o shows strong negative support for species sp_003 and sp_004
#and weaker evidence of a negative effect for the other species

#intercept terms are negative for all species and describe baseline
#occurrence probability when predictors equal zero


## Plot environmental responses ##

#plot posterior support for each Beta coefficient
plotBeta(m,
         post = postBeta,
         param = "Support",
         supportLevel = 0.95)

#rows correspond to species
#columns correspond to environmental predictors

#red indicates strong posterior support for a positive effect
#blue indicates strong posterior support for a negative effect
#white indicates weak or uncertain support

#x2o represents the categorical predictor x2 after conversion to a dummy variable
#it indicates the effect of category "o" relative to the reference category "c"


## Interpretation for this model ##

#predictor x1 shows strong positive effects for all species
#indicating that occurrence probability increases with x1

#predictor x2o shows strong negative effects for species sp_003
#and sp_004, indicating lower occurrence probability in category "o"
#relative to the reference category "c"; effects on the other species
#are weaker and not strongly supported

#intercept terms are mostly negative and represent baseline
#occurrence probability when predictors equal zero and x2 = "c"


## Compute variance partitioning ##

#calculate how explained variation is attributed to model components
VP = computeVariancePartitioning(m)

#VP stores the proportion of explained variation assigned to each component

VP

#VP$vals shows proportion of explained variation for each species
#attributed to each predictor and random effect


## Interpretation for this model ##

#predictor x1 explains most of the explained variation for species
#sp_003, sp_004, and sp_001, indicating strong environmental control by x1

#predictor x2 explains only a small proportion of the explained variation
#for all species

#plot-level random effect explains a large proportion of the explained
#variation for species sp_002, suggesting strong spatially structured
#residual variation not captured by the environmental predictors

#sample-level random effect contributes very little explained variation
#for all species

#overall the model explains approximately 73% of the variation across species


## Plot variance partitioning results ##

#visualise how explained variation is divided among predictors
plotVariancePartitioning(m, VP = VP)

#each bar corresponds to one species

#colours represent contributions of environmental predictors
#and random effects to the explained variation


## Compute residual species associations ##

#calculate residual association matrix between species
OmegaCor = computeAssociations(m)

#OmegaCor stores correlation structure at each random-effect level

OmegaCor[[1]]


## Interpretation of residual species associations ##

#residual correlations range between -1 and 1

#values close to zero indicate little or no association after accounting
#for environmental predictors and random effects

#all residual correlations in this model are small (close to zero),
#indicating weak unexplained co-occurrence structure between species

#posterior support values are far from the 0.95 threshold,
#so none of the residual associations are strongly supported

#this suggests that most species co-occurrence patterns are already
#explained by environmental predictors rather than biotic interactions


## Plot residual species associations with strong posterior support ##

#select posterior support threshold
supportLevel = 0.95

#retain only correlations with strong support
toPlot =
  ((OmegaCor[[1]]$support > supportLevel) +
     (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

#plot residual correlation matrix
corrplot(toPlot,
         method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         title = paste("random effect level:", m$rLNames[1]),
         mar = c(0,0,1,0))

#blue indicates negative residual association
#red indicates positive residual association
#white indicates no strongly supported association


## Interpretation of residual species association plot ##

#the plot shows residual correlations between species after accounting
#for environmental predictors (x1, x2) and random effects (sample, plot)

#only correlations with strong posterior support (≥ 0.95 or ≤ 0.05)
#are displayed in colour

#in this model, no species pairs show strongly supported residual
#associations, so the plot appears mostly white

#this indicates that most co-occurrence patterns among species are
#already explained by environmental predictors and spatial structure

#there is no evidence for additional structured residual associations
#that could suggest strong biotic interactions or missing predictors

#overall, community structure in this dataset is primarily driven by
#environmental variation rather than unexplained species-to-species
#dependencies
