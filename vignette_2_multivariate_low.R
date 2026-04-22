#load packages
library(Hmsc); library(corrplot)

#set seed
set.seed(1)


## Generate simulated community data ##

n = 100                                 #number of sampling units (sites)
x1 = rnorm(n)                           #first environmental covariate
x2 = rnorm(n)                           #second environmental covariate
XData = data.frame(x1 = x1, x2 = x2)    #predictor data used in the model

alpha = c(0,0,0,0,0)                    #true intercept for each species
beta1 = c(1,1,-1,-1,0)                  #true effect of x1 on each species
beta2 = c(1,-1,1,-1,0)                  #true effect of x2 on each species
sigma = c(1,1,1,1,1)                    #residual variation for each species

L = matrix(NA, nrow = n, ncol = 5)      #linear predictor for each species
Y = matrix(NA, nrow = n, ncol = 5)      #observed response matrix (sites x species)

for (j in 1:5){
  L[,j] = alpha[j] + beta1[j]*x1 + beta2[j]*x2   #deterministic response of species j
  Y[,j] = L[,j] + rnorm(n, sd = sigma[j])         #observed response = predictor + residual variation
}

colnames(Y) = paste0("sp", 1:5)         #species names

#species 1 and 2 respond positively to x1
#species 3 and 4 respond negatively to x1
#species responses to x2 differ among species
#species 5 does not respond to either predictor


## Construct multivariate HMSC model ##

#fit HMSC model with additive effects of x1 and x2 on all species
m = Hmsc(Y = Y, XData = XData, XFormula = ~x1+x2)

#the same environmental predictors are used for all species
#HMSC will estimate species-specific intercepts and slopes


## Fit multivariate HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               verbose = verbose)

#model now contains posterior samples of species-specific regression coefficients


## Define MCMC sampling settings ##

#number of independent MCMC chains
nChains = 2

#TRUE = fast test run (for learning syntax), FALSE = full run (reliable estimates)
test.run = FALSE

#fast run (not reliable parameter estimates)
if (test.run){
  
  thin = 1                         #keep every sample
  samples = 10                     #number of posterior samples per chain
  transient = 5                    #burn-in period (discard early unstable samples)
  verbose = 0                      #suppress progress output
  
  #full run (reproduces results shown in the vignette)
} else {
  
  thin = 10                        #keep every 10th sample to reduce autocorrelation
  samples = 1000                   #number of posterior samples per chain
  transient = 500*thin             #burn-in period
  verbose = 0                      #suppress progress output
}


## Fit multivariate HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               nParallel = nChains,
               verbose = verbose)

#model now contains posterior samples of species-specific regression coefficients


## Convert posterior samples into coda format ##

#convert fitted model output for MCMC diagnostics
mpost = convertToCodaObject(m)


## Check effective sample sizes of regression parameters ##

#effective sample size for each beta parameter
effectiveSize(mpost$Beta)

#larger values indicate better mixing of the MCMC chains


## Check Gelman–Rubin convergence diagnostic ##

#potential scale reduction factor (psrf) for each beta parameter
gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

#values close to 1 indicate convergence across chains


## Visualise convergence diagnostics ##

#plot distributions of effective sample sizes and psrf values
par(mfrow = c(1,2))

hist(effectiveSize(mpost$Beta),
     main = "ess(beta)")

hist(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf,
     main = "psrf(beta)")


## Evaluate explanatory power of multivariate model ##

#posterior predicted values for each species at each sampling unit
preds = computePredictedValues(m)

#evaluate model fit
evaluateModelFit(hM = m, predY = preds)

#R2 = proportion of variance explained for each species
#RMSE = prediction error for each species


## Evaluate predictive power with cross-validation ##

#create two-fold cross-validation partition
partition = createPartition(m, nfolds = 2)

#predict responses for held-out data
preds = computePredictedValues(m,
                               partition = partition,
                               nParallel = nChains)

#evaluate predictive performance
evaluateModelFit(hM = m, predY = preds)

#returns RMSE and R2 for each species based on cross-validation


## Extract posterior estimates of regression coefficients (Beta) ##

#retrieve posterior summaries of Beta parameters
postBeta = getPostEstimate(m, parName = "Beta")

#plot support for each coefficient
plotBeta(m,
         post = postBeta,
         param = "Support",
         supportLevel = 0.95)

#red = positive effect with strong posterior support
#blue = negative effect with strong posterior support
#white = no supported effect


## Estimate species-to-species residual associations via latent factors ##

#create study design describing sampling units (one row per observation)
studyDesign = data.frame(sample = as.factor(1:n))

#define random effect at sampling-unit level
rL = HmscRandomLevel(units = studyDesign$sample)

#refit model including latent random effect
m = Hmsc(Y = Y,
         XData = XData,
         XFormula = ~x1 + x2,
         studyDesign = studyDesign,
         ranLevels = list(sample = rL))

#run MCMC sampling to estimate parameters including residual associations
m = sampleMcmc(m,
               thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               nParallel = nChains,
               verbose = verbose)


## Conceptual note: why add a sampling-unit random effect? ##

#adding a random effect at the sampling-unit level introduces latent variables
#that capture residual co-variation among species after accounting for predictors

#this allows the model to estimate species-to-species residual associations

#in a multivariate model, these latent factors describe shared structure across species
#for example: species that tend to occur together (positive association)
#or species that tend to avoid each other (negative association)

#these associations are stored in the residual correlation matrix (Omega)

#biologically, residual associations may reflect:
#biotic interactions (competition, facilitation)
#missing environmental predictors
#dispersal limitation
#unmeasured habitat structure

#this step turns the model from a multivariate regression into a joint SDM


## Check MCMC convergence for both Beta and Omega parameters ##

#convert posterior samples to coda object
mpost = convertToCodaObject(m)

#store current plotting settings
oldpar = par(no.readonly = TRUE)

#set plotting layout
par(mfrow = c(2,2))

#effective sample size for environmental effects
hist(effectiveSize(mpost$Beta), main = "ess(beta)")

#Gelman-Rubin diagnostic for environmental effects
hist(gelman.diag(mpost$Beta, multivariate = FALSE)$psrf,
     main = "psrf(beta)")

#effective sample size for species association matrix (Omega)
hist(effectiveSize(mpost$Omega[[1]]), main = "ess(omega)")

#Gelman-Rubin diagnostic for Omega
hist(gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf,
     main = "psrf(omega)")

#restore plotting settings
par(oldpar)


## Re-check posterior estimates of regression coefficients (Beta) ##

#retrieve posterior summaries of Beta parameters
postBeta = getPostEstimate(m, parName = "Beta")

#plot support for each coefficient after adding latent random effect
plotBeta(m,
         post = postBeta,
         param = "Support",
         supportLevel = 0.95)

#Beta estimates should remain consistent with the model without random effects


## Extract species-to-species residual correlations ##

#compute residual association matrix
OmegaCor = computeAssociations(m)

#set posterior support threshold
supportLevel = 0.95

#retain only strongly supported correlations
toPlot =
  ((OmegaCor[[1]]$support > supportLevel) +
     (OmegaCor[[1]]$support < (1 - supportLevel)) > 0) *
  OmegaCor[[1]]$mean

#plot correlation matrix
corrplot(toPlot,
         method = "color",
         col = colorRampPalette(c("blue","white","red"))(200),
         title = paste("random effect level:", m$rLNames[1]),
         mar = c(0,0,1,0))

#only strongly supported residual correlations are shown
#diagonal values are always 1 (within-species correlation)
#off-diagonal values represent residual associations between species