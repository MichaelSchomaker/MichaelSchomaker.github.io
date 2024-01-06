library(CICI)
data(EFV)

# 1) STANDARD PARAMETRIC G-FORMULA
est <- gformula(X=EFV,
                Lnodes  = c("adherence.1","weight.1",
                            "adherence.2","weight.2",
                            "adherence.3","weight.3",
                            "adherence.4","weight.4"
                ),
                Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                abar=seq(0,10,1)
)

est        # print estimates
plot(est)  # plot estimates


# Bootstrap and CI (and using parallelization)
est_CI <- gformula(X=EFV,
                   Lnodes  = c("adherence.1","weight.1",
                               "adherence.2","weight.2",
                               "adherence.3","weight.3",
                               "adherence.4","weight.4"
                   ),
                   Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                   Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                   abar=seq(0,10,1), B=20, ncores=7 
)

est_CI
plot(est_CI, CI=T, time.points=c(1,5)) # check ?plot.gformula

# Other estimands than E(Y^(abar))
est <- gformula(X=EFV,
                Lnodes  = c("adherence.1","weight.1",
                            "adherence.2","weight.2",
                            "adherence.3","weight.3",
                            "adherence.4","weight.4"
                ),
                Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                abar=seq(0,2,1), ret=T
)

est # default
custom.measure(est, fun=prop,categ=1) # P(Y^a=1) - > identical results
custom.measure(est, fun=prop,categ=0) # P(Y^a=0)
custom.measure(est, fun=prop, categ=0, cond="sex==1") # conditional on sex=1
# does not make sense here, just for illustration
custom.measure(est, fun=quantile, probs=0.1) 


# Model specification
m <- make.model.formulas(X=EFV,
                         Lnodes  = c("adherence.1","weight.1",
                                     "adherence.2","weight.2",
                                     "adherence.3","weight.3",
                                     "adherence.4","weight.4"
                         ),
                         Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                         Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4")
)
m$model.names # all models potentially relevant for gformula(), given *full* past

# Update model formulas (automated): screening with LASSO 
glmnet.formulas <-  model.formulas.update(m$model.names, EFV)


# use these models in gformula
est <- gformula(X=EFV,
                Lnodes  = c("adherence.1","weight.1",
                            "adherence.2","weight.2",
                            "adherence.3","weight.3",
                            "adherence.4","weight.4"
                ),
                Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                Yform=glmnet.formulas$Ynames, Lform=glmnet.formulas$Lnames,
                abar=seq(0,2,1)
)
est

# Update model formulas (manual) and include penalized splines:
# ? fit.updated.formulas
# ? model.update



# Multiple Imputation
# suppose the following subsets were actually multiply imputed data (M=2)
EFV_1 <- EFV[1:2500,]
EFV_2 <- EFV[2501:5000,]

# first: conduct analysis on each imputed data set. Set ret=T.
m1 <- gformula(X=EFV_1,
               Lnodes  = c("adherence.1","weight.1",
                           "adherence.2","weight.2",
                           "adherence.3","weight.3",
                           "adherence.4","weight.4"
               ),
               Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
               Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
               abar=seq(0,5,1), verbose=F, ret=T
)

m2 <- gformula(X=EFV_2,
               Lnodes  = c("adherence.1","weight.1",
                           "adherence.2","weight.2",
                           "adherence.3","weight.3",
                           "adherence.4","weight.4"
               ),
               Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
               Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
               abar=seq(0,5,1), verbose=F, ret=T
)

# second combine results
m_imp <- mi.boot(list(m1,m2), mean) # uses MI rules & returns 'gformula' object
plot(m_imp)

# custom estimand: evaluate probability of suppression (Y=0), among females
m_imp2 <- mi.boot(list(m1,m2), prop, categ=0, cond="sex==1")
plot(m_imp2)


# Read more about other useful options in the details section of ?gformula:
# survival settings
# diagnostics
# custom interventions
# parallelization using the options 'ncores'
# track progress, even under parallelizations with 'prog'
# For reproducability use 'seed'


#### WEIGHETD SEQUENTIAL G-FORMULA
library(CICIplus) 

# Seq. g-formula
est2 <- sgf(X=EFV,
            Lnodes  = c("adherence.1","weight.1",
                        "adherence.2","weight.2",
                        "adherence.3","weight.3",
                        "adherence.4","weight.4"
            ),
            Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
            Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
            abar=seq(0,10,1)
)
est2
plot(est2)



# Weighted seq. g-formula to address positivity violations

# Step 1: calculate weights (add baseline variables to Lnodes)
w <- calc.weights(X=EFV, 
                         Lnodes  = c("sex", "metabolic", "log_age",
                                     "NRTI" ,"weight.0",
                                     "adherence.1","weight.1",
                                     "adherence.2","weight.2",
                                     "adherence.3","weight.3",
                                     "adherence.4","weight.4"),
                         Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
                         Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
                         d.method="binning", abar=seq(0,10,1),
                         c=c(0.01,0.001)
)
summary(w)
# d.method alternatives: parametric density or HAL (for small n)
# check ?calc.weights


# Step 2: sequential g-formula with outcome weights
est3 <- sgf(X=EFV,
              Lnodes  = c("adherence.1","weight.1",
                          "adherence.2","weight.2",
                          "adherence.3","weight.3",
                          "adherence.4","weight.4"
              ),
              Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
              Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
              Yweights = w$`0.01`,
              abar=seq(0,10,1)
)



# Realistic setup
# (requires user-written learners which are available upon request)

library(SuperLearner) 
est5 <- sgf(X=EFV,
            Lnodes  = c("adherence.1","weight.1",
                        "adherence.2","weight.2",
                        "adherence.3","weight.3",
                        "adherence.4","weight.4"
            ),
            Ynodes  = c("VL.0","VL.1","VL.2","VL.3","VL.4"),
            Anodes  = c("efv.0","efv.1","efv.2","efv.3","efv.4"),
            Yweights = w$`0.01`,
            SL.library = list(c("SL.glm","screen.glmnet_nVar_1_4_10_150"),
                              c("SL.mgcv_base", "screen.cramersv_3")),
            SL.export = c("SL.mgcv_base","screen.glmnet_nVar_base",
                          "screen.cramersv_base","predict.SL.mgcv"),
            abar=seq(0,10,1), B=200, ncores=7, seed=41188, prog="C:/temp",
            cvControl = list(V=2), # SL option
            calc.support=TRUE
)
est5                                    # estimates
est5$SL.weights                         # Super Learner summary
plot(est5, time.points = c(1,5), CI=T)  # plot with CIs

