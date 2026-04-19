#load packages
library(Hmsc)

#set seed
set.seed(1)

## Generate simulated data ##

n = 50                          #number of sampling units (sites)
x = rnorm(n)                    #environmental covariate
alpha = 0                       #true intercept
beta = 1                        #true effect of x on the response
sigma = 1                       #residual variation (unexplained noise)
L = alpha + beta*x              #linear predictor (deterministic part of the model)
y = L + rnorm(n, sd = sigma)    #observed response = predictor + residual variation

plot(x, y, las=1)                   
 


## Fit linear model and estimate parameters ##

#combine predictor and response into a data frame
df = data.frame(x, y)    

#fit linear regression: estimate intercept (alpha) and slope (beta)
m.lm = lm(y ~ x, data = df)     

#compare estimated values with true parameters used to simulate the data
summary(m.lm)

# note #

#lm estimates the relationship between x and y assuming residual variation is independent


## Construct HMSC model ##

#convert response vector to matrix format required by HMSC
Y = as.matrix(y)     

#environmental predictor data
XData = data.frame(x = x)           

#define model structure (response explained by predictor x)
m = Hmsc(Y = Y, XData = XData, XFormula = ~x)   

#default response distribution in HMSC is Gaussian
#so this model is equivalent to Bayesian linear regression

#same model written explicitly
#m = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "normal")

## Fit model with Bayesian MCMC sampling ##

#sampleMcmc estimates posterior distributions of model parameters instead of single point estimates

#key settings:
#nChains = number of independent sampling chains
#samples = number of posterior samples retained
#transient = burn-in period removed before analysis
#verbose = controls progress output


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
  verbose = 5                      #print progress every 5 iterations

#full run (reproduces results shown in the vignette)
} else {
  
  thin = 5                         #keep every 5th sample to reduce autocorrelation
  samples = 1000                   #number of posterior samples per chain
  transient = 500*thin             #burn-in period
  verbose = 500*thin               #print progress every 2500 iterations
}



## Fit HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               verbose = verbose)

#after this step, m is no longer just the model structure: it now contains fitted parameter estimates

#(Bayesian posterior samples, not single point estimates like in lm)


## Extract posterior estimates ##

#convert fitted model output into coda format for inspection
mpost = convertToCodaObject(m)      

#summarise posterior estimates of regression coefficients (intercept and slope)
summary(mpost$Beta)                

#Beta contains the environmental effects estimated by the model

#posterior mean of intercept ≈ 0 (true value used in simulation)
#posterior mean of slope ≈ 1 (true value used in simulation)
#95% credible interval of slope does not include 0 → predictor x has a clear positive effect


## Evaluate model fit ##

#posterior predicted values for each observation
preds = computePredictedValues(m)          

#compute RMSE and R² for model predictions
evaluateModelFit(hM = m, predY = preds)    

#RMSE = prediction error (expected ≈ sigma used in simulation)
#R2 = proportion of variance explained by predictor x


## Check MCMC convergence ##

plot(mpost$Beta)   #trace plots for regression coefficients across MCMC chains

#good trace plots: chains overlap, mix well, and show no trend across iterations


## Quantitative MCMC convergence diagnostics ##

#effective n of independent posterior samples (values close to 2000 = excellent)
effectiveSize(mpost$Beta)    

#PSRF ≈ 1 indicates both chains converged to same posterior distribution
gelman.diag(mpost$Beta, multivariate=FALSE)$psrf  

#values close to 1 indicate good convergence; values >1.1 would suggest problems

#effective sample sizes are high, indicating little autocorrelation among posterior samples
#PSRF values are ~1, indicating excellent agreement between the two MCMC chains


## Check assumptions of linear model (lm) ##

nres.lm = rstandard(m.lm)        #standardised residuals from lm model
preds.lm = fitted.values(m.lm)   #predicted values from lm model

par(mfrow=c(1,2))                #plot two panels side by side
hist(nres.lm, las = 1)           #check residual normality
plot(preds.lm, nres.lm, las = 1) #check homoscedasticity (constant variance)
abline(a=0,b=0)                  #reference line at zero residual


## Check assumptions of HMSC model ##

preds.mean = apply(preds, FUN=mean, MARGIN=1)  #posterior mean predss per obs
nres = scale(y - preds.mean)                   #standardised residuals

par(mfrow=c(1,2))
hist(nres, las = 1)                            #check residual normality
plot(preds.mean,nres, las = 1)                  #check homoscedasticity
abline(a=0,b=0)


## Generate simulated presence–absence data ##

#simulate binary response using latent-variable formulation of probit regression
#presence occurs when linear predictor + noise exceeds threshold (0)
y = 1*(L + rnorm(n, sd = 1) > 0)

plot(x, y, las=1)


## Construct probit HMSC model ##

#convert response vector to matrix format required by HMSC
Y = as.matrix(y)

#define probit regression model (presence–absence response)
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "probit")

#probit models probability of presence as a function of predictor x


## Fit probit HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               nParallel = nChains,
               verbose = verbose)

#model now contains posterior samples of regression parameters


## Extract posterior estimates ##

mpost = convertToCodaObject(m)

summary(mpost$Beta)

#Beta now represents effect of predictor x on probability of presence
#positive slope indicates higher x increases probability species is present

#model structure identical to Gaussian example, but response distribution is now binary (probit link)


## Quantitative MCMC convergence diagnostics ##

effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta, multivariate=FALSE)$psrf

#effective sample sizes slightly smaller than Gaussian case (binary data contain less information)
#PSRF values ≈ 1 indicate good agreement between chains and reliable convergence


## Evaluate model fit (probit model) ##

#posterior predicted probability of presence for each observation
preds = computePredictedValues(m)

#evaluate predictive performance of presence–absence model
evaluateModelFit(hM = m, predY = preds)

#RMSE = prediction error (less interpretable for binary responses than for Gaussian models)

#AUC = ability of model to discriminate presences from absences (values closer to 1 indicate better performance)

#TjurR2 = proportion of variance explained for binary-response models

#for binary responses, model fit is evaluated using AUC and TjurR2 instead of classical R²

#binary response data contain less information than continuous responses, so parameter uncertainty is larger


## Generate simulated count data ##

#simulate count response using Poisson regression formulation
#expected count increases with linear predictor L

y = rpois(n, lambda = exp(L))

plot(x, y, las = 1)

#Poisson model assumes response variable consists of non-negative integer counts (0,1,2,...)


## Construct Poisson HMSC model ##

#convert response vector to matrix format required by HMSC
Y = as.matrix(y)

#define Poisson regression model (count response)
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "poisson")

#Poisson regression models expected count as a function of predictor x


## Fit Poisson HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               nParallel = nChains,
               verbose = verbose)

#model now contains posterior samples of regression parameters


## Extract posterior estimates ##

mpost = convertToCodaObject(m)

summary(mpost$Beta)

#Beta now represents effect of predictor x on expected count
#positive slope indicates higher x increases expected number of individuals


## Quantitative MCMC convergence diagnostics ##

effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

#effective sample sizes slightly < than Gaussian case (count data contain less information)
#PSRF values ≈ 1 indicate good agreement between chains and reliable convergence


## Evaluate model fit (Poisson model) ##

#posterior predicted counts for each observation
preds = computePredictedValues(m)

#evaluate predictive performance of count-response model
evaluateModelFit(hM = m, predY = preds)

#RMSE = prediction error on count scale
#R2 = proportion of variance explained by predictor x

#Poisson regression assumes mean and variance of counts are approximately equal
#if variance >> mean, lognormal Poisson model is usually more appropriate


## Generate simulated lognormal Poisson count data ##

#standard Poisson regression assumes mean ≈ variance
#ecological count data often show extra variation (overdispersion)

#lognormal Poisson model allows additional unexplained variability
#by adding Gaussian noise to the linear predictor before exponentiation

y = rpois(n, lambda = exp(L + rnorm(n, sd = 2)))

plot(x, y, las = 1)

#counts now vary more strongly around the expected relationship with x
#this better reflects realistic ecological abundance data

#lognormal Poisson plays a role similar to negative binomial regression
#but is the overdispersed count model currently implemented in HMSC

#extra Gaussian noise represents unmeasured ecological drivers
#in jSDMs this becomes latent structure shared across species


## Construct lognormal Poisson HMSC model ##

#convert response vector to matrix format required by HMSC
Y = as.matrix(y)

#define lognormal Poisson regression model for overdispersed count data
m = Hmsc(Y = Y, XData = XData, XFormula = ~x, distr = "lognormal poisson")

#lognormal Poisson model allows extra variation around expected counts
#this is often more realistic for ecological abundance data than standard Poisson


## Fit lognormal Poisson HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               nParallel = nChains,
               verbose = verbose)

#model now contains posterior samples of regression parameters


## Extract posterior estimates ##

mpost = convertToCodaObject(m)

summary(mpost$Beta)

#Beta now represents effect of predictor x on expected abundance
#positive slope indicates higher x increases expected counts
#model also allows additional unexplained variation around this relationship


## Quantitative MCMC convergence diagnostics ##

effectiveSize(mpost$Beta)

gelman.diag(mpost$Beta, multivariate = FALSE)$psrf

#mixing is typically worse than in Gaussian and simpler Poisson models
#because the model is more flexible and harder to sample efficiently


## Evaluate model fit (lognormal Poisson model) ##

#use realised posterior predictions, not expected values
#this is required for occurrence and conditional-abundance fit measures

preds = computePredictedValues(m, expected = FALSE)

#evaluate predictive performance of overdispersed count-response model
evaluateModelFit(hM = m, predY = preds)

#RMSE and SR2 = overall fit for count data
#O.AUC, O.TjurR2, O.RMSE = fit for occurrence component (presence vs absence)
#C.SR2 and C.RMSE = fit for abundance conditional on presence

#lognormal Poisson model separates three questions:
#can the model predict counts overall
#can it predict occurrence
#can it predict abundance where the species is present


## Generate simulated hierarchical data ##

#number of sampling units
n = 100

#environmental predictor
x = rnorm(n)

#true fixed-effect parameters
alpha = 0
beta = 1
sigma = 1

#linear predictor (fixed-effect component)
L = alpha + beta*x

#number of plots
np = 10

#standard deviation among plots
sigma.plot = 1

#ID of each sampling unit
sample.id = 1:n

#assign each sampling unit to a plot
plot.id = sample(1:np, n, replace = TRUE)

#simulate plot-level random intercepts
ap = rnorm(np, sd = sigma.plot)

#assign plot effect to each sampling unit
a = ap[plot.id]

#observed response = fixed effect + plot effect + residual noise
y = L + a + rnorm(n, sd = sigma)

#convert plot IDs to factor
plot.id = as.factor(plot.id)

#visualise response by plot
plot(x, y, col = plot.id, las = 1)


## Prepare data for HMSC ##

#environmental predictor data
XData = data.frame(x = x)

#convert response vector to matrix format required by HMSC
Y = as.matrix(y)


## Construct hierarchical HMSC model ##

#define study design: each sampling unit belongs to one plot
studyDesign = data.frame(sample = as.factor(sample.id),
                         plot = as.factor(plot.id))

#define random level corresponding to plots
rL = HmscRandomLevel(units = studyDesign$plot)

#fit Gaussian HMSC model with plot-level random intercept
m = Hmsc(Y = Y,
         XData = XData,
         XFormula = ~x,
         studyDesign = studyDesign,
         ranLevels = list(plot = rL))


## Fit hierarchical HMSC model ##

m = sampleMcmc(m, thin = thin,
               samples = samples,
               transient = transient,
               nChains = nChains,
               nParallel = nChains,
               verbose = verbose)

#model now contains posterior samples for both fixed effects and plot-level random effects


## Evaluate explanatory power of hierarchical model ##

#predict fitted values using same data as model fitting
preds = computePredictedValues(m)

#evaluate model fit
MF = evaluateModelFit(hM = m, predY = preds)

#extract R2
MF$R2


## Cross-validation by sampling unit ##

#assign each sampling unit to one of two folds
partition = createPartition(m, nfolds = 2, column = "sample")

#compute cross-validated predictions
preds = computePredictedValues(m, partition = partition, nParallel = nChains)

#evaluate predictive performance
MF = evaluateModelFit(hM = m, predY = preds)

#extract predictive R2
MF$R2

#this evaluates prediction for new sampling units from plots already represented in training data


## Cross-validation by plot ##

#assign each plot to one of two folds
partition = createPartition(m, nfolds = 2, column = "plot")

#inspect plot-level partitioning
t(cbind(plot.id, partition)[1:15,])

#this evaluates prediction for sampling units belonging to entirely new plots

#compute cross-validated predictions
preds = computePredictedValues(m, partition = partition, nParallel = nChains)

#evaluate predictive performance
MF = evaluateModelFit(hM = m, predY = preds)

#extract predictive R2
MF$R2

#this evaluates prediction for sampling units in entirely new plots
#predictive performance is lower because the model cannot estimate the random intercept for unseen plots
#predictions therefore rely mainly on fixed effects


## Generate simulated spatially structured data ##

#standard deviation of spatial effect
sigma.spatial = 2

#range parameter controlling how quickly spatial correlation decays with distance
alpha.spatial = 0.5

#create unique ID for each sampling unit
sample.id = rep(NA, n)
for (i in 1:n){
  sample.id[i] = paste0("location_", as.character(i))
}
sample.id = as.factor(sample.id)

#simulate 2D coordinates for sampling units
xycoords = matrix(runif(2*n), ncol = 2)
rownames(xycoords) = sample.id
colnames(xycoords) = c("x-coordinate", "y-coordinate")

#simulate spatially autocorrelated random effect
a = MASS::mvrnorm(mu = rep(0, n),
                  Sigma = sigma.spatial^2*exp(-as.matrix(dist(xycoords))/alpha.spatial))

#observed response = fixed effect + spatial effect + residual noise
y = L + a + rnorm(n, sd = sigma)

#convert response vector to matrix format required by HMSC
Y = as.matrix(y)


## Visualise spatial structure ##

colfunc = colorRampPalette(c("cyan", "red"))
ncols = 100
cols = colfunc(ncols)

par(mfrow = c(1,2))
for (i in 1:2){
  if (i == 1) value = x
  if (i == 2) value = y
  value = value - min(value)
  value = 1 + (ncols - 1)*value/max(value)
  plot(xycoords[,1], xycoords[,2], col = cols[value], pch = 16,
       main = c("x", "y")[i], asp = 1)
}
