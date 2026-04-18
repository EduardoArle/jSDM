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

#RMSE = prediction error (expected ≈ sigma used in simulation)
#R2 = proportion of variance explained by predictor x


## Check MCMC convergence ##

plot(mpost$Beta)   #trace plots for regression coefficients across MCMC chains

#good trace plots: chains overlap, mix well, and show no trend across iterations

