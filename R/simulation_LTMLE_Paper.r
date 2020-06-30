###################################################################################################################################################
# Supplementary Material to Schomaker M, Luque Fernandez MA, Leroy V, Davies MA.                                                                  #
# Using Longitudinal Targeted Maximum Likelihood Estimation in Complex Settings with Dynamic Interventions.                                       #
# Statistics in Medicine. 2019;38:4888-911.                                                                                                       #
#                                                                                                                                                 #
# Simulation from Section 6                                                                                                                       #
# speed can be improved by i) reduction of learner ser 3 (see below) and ii) reduction of simulation runs (see below)                             #
#                                                                                                                                                 #
###################################################################################################################################################


# Set working directory
setwd('C:/Users/schomakm/Dropbox/Documents/Research_Statistics/Code_Reproduce/LTMLE')
#setwd("//home//mschomaker//Code_Reproduce//Bootstrap")

directory <- getwd()

library(simcausal)       # currently only on CRAN archive (June 2020)
library(SuperLearner)
library(ltmle)
library(ggplot2)

set.seed(241280)

###########################
# Data generating process #
###########################


# definition: left-truncated normal distribution
rnorm_trunc <- function(n, mean, sd, minval = 0, maxval = 10000,
                        min.low = 0, max.low = 50, min.high = 5000, max.high = 10000)
{                                                                                      
  out <- rnorm(n = n, mean = mean, sd = sd)
  minval <- minval[1]; min1 <- min.low[1]; max1 <- max.low[1]
  maxval <- maxval[1]; min2 <- min.high[1]; max2 <- max.high[1]
  leq.zero <- length(out[out <= minval])
  geq.max <- length(out[out >= maxval])
  out[out <= minval] <- runif(n = leq.zero, min = min1, max = max1)
  out[out >= maxval] <- runif(n = geq.max, min = min2, max = max2)
  out
}

#define the number of time points
t.end <- 12

#initialize the dag
D <- DAG.empty()

#everything at baseline (t = 0)
D.base <- D +
  node("V1",                                      #region -> 0 = west, 1 = southern 
       t = 0,
       distr = "rbern",
       prob = 4392/5826) +
  node("V2",                                      #sex -> 0 = female, 1 = male 
       t = 0,
       distr = "rbern",
       prob = ifelse(V1[0] == 1, 2222/4392, 758/1434)) +
  node("V3",                                      #age 
       t = 0,
       distr = "runif",
       min = 1,
       max = 5) +
  node("L1",                                      #cd4-count 
       t = 0,
       distr = "rnorm_trunc",
       mean = ifelse(V1[0] == 1, 650, 720),
       sd = ifelse(V1[0] == 1, 350, 400),
       minval = 0, maxval = 10000,
       min.low = 0, max.low = 50, min.high = 5000, max.high = 10000) +
  node("L1scaled",                                #auxilliary-variable 
       t = 0,
       distr = "rnorm",
       mean = (L1[0]-671.7468)/(10*352.2788)+1,
       sd = 0) +
  node("L2",                                      #cd4% 
       t = 0,
       distr = "rnorm_trunc",
       mean = .16 + .05 * (L1[0] - 650)/650,
       sd = .07,
       minval = .06, maxval = .8,
       min.low = .03, max.low = .09, min.high = .7, max.high = .8) +
  node("L2scaled",                                #auxilliary-variable 
       t = 0,
       distr = "rnorm",
       mean = (L2[0]-0.1648594)/(10*0.06980332)+1,
       sd = 0) +
  node("L3",                                      #waz 
       t = 0,
       distr = "rnorm_trunc",
       mean = ifelse(V1[0] == 1, - 1.65 + .1 * V3[0] + .05 * (L1[0] - 650)/650 + .05 * (L2[0] - 16)/16,
                     - 2.05 + .1 * V3[0] + .05 * (L1[0] - 650)/650 + .05 * (L2[0] - 16)/16),
       sd = 1,
       minval = -5, maxval = 5,
       min.low = -10, max.low = -3, min.high = 3, max.high = 10) +
  node("A",                                       # ART
       t = 0,
       distr = "rbern",
       prob = 0) +
  node("C",                                      # censoring (death?)
       t = 0,
       distr = "rbern",
       prob = 0,
       EFU = T) +
  node("Y",                                      # HAZ
       t = 0,
       distr = "rnorm_trunc",
       mean = -2.6 + .1 * I(V3[0] > 2) + .3 * I(V1[0] == 0) + (L3[0] + 1.45),
       sd = 1.1,
       minval = -5, maxval = 5,
       min.low = -10, max.low = -3, min.high = 3, max.high = 10)

#time-dependent variables at later time-points
D <- D.base +
  node("L1",                                      #cd4-count
       t = 1:4,
       distr = "rnorm_trunc",
       mean = 13*log(t * (1034-662)/8) + L1[t-1] + 2 * L2[t-1] + 2 * L3[t-1] + 2.5 * A[t-1],
       sd = 50,
       minval = 0, maxval = 10000,
       min.low = 0, max.low = 50, min.high = 5000, max.high = 10000) +
  node("L1",                                      #cd4-count
       t = 5:8,
       distr = "rnorm_trunc",
       mean = 4*log(t * (1034-662)/8) + L1[t-1] + 2 * L2[t-1] + 2 * L3[t-1] + 2.5 * A[t-1],
       sd = 50,
       minval = 0, maxval = 10000,
       min.low = 0, max.low = 50, min.high = 5000, max.high = 10000) +
  node("L1",                                      #cd4-count
       t = 9:t.end,
       distr = "rnorm_trunc",
       mean = L1[t-1] + 2 * L2[t-1] + 2 * L3[t-1] + 2.5 * A[t-1],
       sd = 50,
       minval = 0, maxval = 10000,
       min.low = 0, max.low = 50, min.high = 5000, max.high = 10000) +
  node("L2",                                      #cd4%
       t = 1:t.end,
       distr = "rnorm_trunc",
       mean = L2[t-1] + .0003 * (L1[t]-L1[t-1]) + .0005 * (L3[t-1]) + .0005 * A[t-1] * L1scaled[0],
       sd = .02,
       minval = .06, maxval = .8,
       min.low = .03, max.low = .09, min.high = .7, max.high = .8) +
  node("L3",                                      #waz
       t = 1:t.end,
       distr = "rnorm_trunc",
       mean = L3[t-1] + .0017 * (L1[t] - L1[t-1]) + .2 * (L2[t] - L2[t-1]) + .005 * A[t-1] * L2scaled[0],
       sd = .5,
       minval = -5, maxval = 5,
       min.low = -10, max.low = -3, min.high = 3, max.high = 10) +
  node("A",                                       #art
       t = 1:t.end,
       distr = "rbern",
       prob = ifelse(A[t-1] == 1, 1, plogis(-2.4 + .015 * (750 - L1[t]) + 5 * (.2 - L2[t]) - .8 * L3[t] + .8 * t))) +
  node("C",                                      # censoring
       t = 1:t.end,
       distr = "rbern",
       prob = plogis(-6 + .01 * (750 - L1[t]) + 1 * (.2 - L2[t]) - .65 * L3[t] - A[t]),
       EFU = T)

D <- D +
  node("Y",                                      #haz
       t = 1:t.end,
       distr = "rnorm_trunc",
       mean = Y[t-1] +
         .00005 * (L1[t] - L1[t-1]) - 0.000001 * ((L1[t] - L1[t-1])*sqrt(L1scaled[0]))^2 +
         .01 * (L2[t] - L2[t-1]) - .0001 * ((L2[t] - L2[t-1])*sqrt(L2scaled[0]))^2 +
         .07 * ((L3[t]-L3[t-1])*(L3[0]+1.5135)) - .001 * ((L3[t]-L3[t-1])*(L3[0]+1.5135))^2 +
         .005 * A[t] + .075 * A[t-1] + .05 * A[t]*A[t-1] ,
       sd = .01,
       minval = -5, maxval = 5,
       min.low = -10, max.low = -3, min.high = 3, max.high = 10)




#specify (dynamic) interventions
Dset <- set.DAG(D)
int.a <- c(node("A", t = 1:t.end, distr = "rbern",
              prob = ifelse(A[t-1] == 1, 1, ifelse(L1[t] < theta1 | L2[t] < theta2 | L3[t] < theta3, 1, 0))),
           node("C", t = 1:t.end, distr = "rbern", prob = 0))

D.dyn1 <- Dset + action("A_th1", nodes = int.a, theta1 = 10000, theta2 = .99, theta3 = 10) # treat all (A=1)
D.dyn2 <- Dset + action("A_th2", nodes = int.a, theta1 = 750, theta2 = .25, theta3 = -2) # treat 750/25/-2
D.dyn3 <- Dset + action("A_th3", nodes = int.a, theta1 = 350, theta2 = .15, theta3 = -2) # treat 350/15/-2
D.dyn4 <- Dset + action("A_th4", nodes = int.a, theta1 = -1, theta2 = -.1, theta3 = -11) # treat never (A=0)

# Intervention on the DAG given the specified interventions
dat1 <- simcausal::sim(DAG = D.dyn1, actions = "A_th1", n = 1000000, rndseed = 7693)
dat2 <- simcausal::sim(DAG = D.dyn2, actions = "A_th2", n = 1000000, rndseed = 7693)
dat3 <- simcausal::sim(DAG = D.dyn3, actions = "A_th3", n = 1000000, rndseed = 7693)
dat4 <- simcausal::sim(DAG = D.dyn4, actions = "A_th4", n = 1000000, rndseed = 7693)


true <- matrix(NA,ncol=4,nrow=13)
rownames(true) <- 0:12
colnames(true) <- c("immediate ART", "<750; <25%; < -2", "<350; <15%; < -2", "no ART")

true[1,1] <- mean(((dat1[["A_th1"]])$Y_0))
true[2,1] <- mean(((dat1[["A_th1"]])$Y_1))
true[3,1] <- mean(((dat1[["A_th1"]])$Y_2))
true[4,1] <- mean(((dat1[["A_th1"]])$Y_3))
true[5,1] <- mean(((dat1[["A_th1"]])$Y_4))
true[6,1] <- mean(((dat1[["A_th1"]])$Y_5))
true[7,1] <- mean(((dat1[["A_th1"]])$Y_6))
true[8,1] <- mean(((dat1[["A_th1"]])$Y_7))
true[9,1] <- mean(((dat1[["A_th1"]])$Y_8))
true[10,1] <- mean(((dat1[["A_th1"]])$Y_9))
true[11,1] <- mean(((dat1[["A_th1"]])$Y_10))
true[12,1] <- mean(((dat1[["A_th1"]])$Y_11))
true[13,1] <- mean(((dat1[["A_th1"]])$Y_12))

true[1,2] <- mean(((dat2[["A_th2"]])$Y_0))
true[2,2] <- mean(((dat2[["A_th2"]])$Y_1))
true[3,2] <- mean(((dat2[["A_th2"]])$Y_2))
true[4,2] <- mean(((dat2[["A_th2"]])$Y_3))
true[5,2] <- mean(((dat2[["A_th2"]])$Y_4))
true[6,2] <- mean(((dat2[["A_th2"]])$Y_5))
true[7,2] <- mean(((dat2[["A_th2"]])$Y_6))
true[8,2] <- mean(((dat2[["A_th2"]])$Y_7))
true[9,2] <- mean(((dat2[["A_th2"]])$Y_8))
true[10,2] <- mean(((dat2[["A_th2"]])$Y_9))
true[11,2] <- mean(((dat2[["A_th2"]])$Y_10))
true[12,2] <- mean(((dat2[["A_th2"]])$Y_11))
true[13,2] <- mean(((dat2[["A_th2"]])$Y_12))

true[1,3] <- mean(((dat3[["A_th3"]])$Y_0))
true[2,3] <- mean(((dat3[["A_th3"]])$Y_1))
true[3,3] <- mean(((dat3[["A_th3"]])$Y_2))
true[4,3] <- mean(((dat3[["A_th3"]])$Y_3))
true[5,3] <- mean(((dat3[["A_th3"]])$Y_4))
true[6,3] <- mean(((dat3[["A_th3"]])$Y_5))
true[7,3] <- mean(((dat3[["A_th3"]])$Y_6))
true[8,3] <- mean(((dat3[["A_th3"]])$Y_7))
true[9,3] <- mean(((dat3[["A_th3"]])$Y_8))
true[10,3] <- mean(((dat3[["A_th3"]])$Y_9))
true[11,3] <- mean(((dat3[["A_th3"]])$Y_10))
true[12,3] <- mean(((dat3[["A_th3"]])$Y_11))
true[13,3] <- mean(((dat3[["A_th3"]])$Y_12))

true[1,4] <- mean(((dat4[["A_th4"]])$Y_0))
true[2,4] <- mean(((dat4[["A_th4"]])$Y_1))
true[3,4] <- mean(((dat4[["A_th4"]])$Y_2))
true[4,4] <- mean(((dat4[["A_th4"]])$Y_3))
true[5,4] <- mean(((dat4[["A_th4"]])$Y_4))
true[6,4] <- mean(((dat4[["A_th4"]])$Y_5))
true[7,4] <- mean(((dat4[["A_th4"]])$Y_6))
true[8,4] <- mean(((dat4[["A_th4"]])$Y_7))
true[9,4] <- mean(((dat4[["A_th4"]])$Y_8))
true[10,4] <- mean(((dat4[["A_th4"]])$Y_9))
true[11,4] <- mean(((dat4[["A_th4"]])$Y_10))
true[12,4] <- mean(((dat4[["A_th4"]])$Y_11))
true[13,4] <- mean(((dat4[["A_th4"]])$Y_12))

write.csv(true,"true.csv")


##############
# Simulation #
##############

# Main Setup #

#################################################################
R <- 1000  # Number of simulation runs; note: very computer intensive
           # if learner set 3 is modified, speed can be improved  

N <- c(200,600,1000)

learner1 <- list(Q=c("SL.glm"),
                   g=c("SL.glm"))
learner2 <- list(Q=c("SL.glm", "SL.mean", "SL.glm.interaction"),
                   g=c("SL.glm", "SL.mean" , "SL.glm.interaction"))
learner3 <- list(Q=c("SL.glm","SL.mean", "SL.stepAIC", "SL.glm.interaction", "SL.gam", "SL.bayesglm"),
                  g=c("SL.glm","SL.mean", "SL.stepAIC", "SL.glm.interaction", "SL.gam", "SL.bayesglm")
                  )

mylibrary <- list(learner1 , learner2, learner3)       
mylibrary_names <- c("learner1", "learner2", "learner3")            
                  
gbounds <- c(0.01,0.025,0.04)

#################################################################
ptm <- proc.time()
allcombinations <- apply(expand.grid(N, gbounds, mylibrary_names),1,paste, collapse=" ")

results_always_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_always_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_750_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_750_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_350_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_350_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_never_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_never_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_CP_always_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_CP_always_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_CP_750_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_CP_750_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_CP_350_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_CP_350_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_CP_never_6  <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_CP_never_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_ss_always_6 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_ss_750_6 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_ss_350_6 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_ss_never_6 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_ss_always_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_ss_750_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_ss_350_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))
results_ss_never_12 <- matrix(NA,ncol=length(allcombinations), nrow=R, dimnames=list(NULL,allcombinations))

results_SL_g_always_12 <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_g_750_12    <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_g_350_12    <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_g_never_12  <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))

results_SL_g_always_6  <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_g_750_6     <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_g_350_6     <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_g_never_6   <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))

results_SL_Q_always_12 <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_Q_750_12    <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_Q_350_12    <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_Q_never_12  <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))

results_SL_Q_always_6  <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_Q_750_6     <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_Q_350_6     <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))
results_SL_Q_never_6   <- matrix(NA,ncol=length(mylibrary[length(mylibrary)][[1]][[1]]), nrow=R, dimnames=list(NULL,mylibrary[length(mylibrary)][[1]][[1]]))

# extract SL weights
ew  <- function(myfit){return(myfit[,2])}
ewg <- function(myfit){
if(is.character(myfit[[1]])){return(NULL)}else{
return(myfit[,2])}
}

###############################
myQform <- c(Y_1  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1  + L3_1  + A_1 ",
             Y_2  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1  + L3_1  + A_1  + Y_1  + L1_1  + L2_1  + L3_1  + A_1",
             Y_3  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2  + L3_2  + A_2  + Y_2  + L1_2  + L2_2  + L3_2  + A_2",
             Y_4  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3  + L3_3  + A_3  + Y_3  + L1_3  + L2_3  + L3_3  + A_3",
             Y_5  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4  + L3_4  + A_4  + Y_4  + L1_4  + L2_4  + L3_4  + A_4",
             Y_6  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5  + L3_5  + A_5  + Y_5  + L1_5  + L2_5  + L3_5  + A_5",
             Y_7  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6  + L3_6  + A_6  + Y_6  + L1_6  + L2_6  + L3_6  + A_6",
             Y_8  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7  + L3_7  + A_7  + Y_7  + L1_7  + L2_7  + L3_7  + A_7",
             Y_9  = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8  + L3_8  + A_8  + Y_8  + L1_8  + L2_8  + L3_8  + A_8",
             Y_10 = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9  + L3_9  + A_9  + Y_9  + L1_9  + L2_9  + L3_9  + A_9",
             Y_11 = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10 + L3_10 + A_10 + Y_10 + L1_10 + L2_10 + L3_10 + A_10",
             Y_12 = "Q.kplus1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11 + L3_11 + A_11 + Y_11 + L1_11 + L2_11 + L3_11 + A_11"
    )
    
mygform <- c("A_1  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1",  
             "C_1  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1",
             "A_2  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1   + Y_1   + L1_2  + L2_2  + L3_2", 
             "C_2  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1   + Y_1   + L1_2  + L2_2  + L3_2  + A_2",
             "A_3  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2   + L3_2   + A_2   + Y_2   + L1_3  + L2_3  + L3_3",  
             "C_3  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2   + L3_2   + A_2   + Y_2   + L1_3  + L2_3  + L3_3  + A_3", 
             "A_4  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3   + L3_3   + A_3   + Y_3   + L1_4  + L2_4  + L3_4",  
             "C_4  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3   + L3_3   + A_3   + Y_3   + L1_4  + L2_4  + L3_4  + A_4", 
             "A_5  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4   + L3_4   + A_4   + Y_4   + L1_5  + L2_5  + L3_5",  
             "C_5  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4   + L3_4   + A_4   + Y_4   + L1_5  + L2_5  + L3_5  + A_5", 
             "A_6  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5   + L3_5   + A_5   + Y_5   + L1_6  + L2_6  + L3_6",  
             "C_6  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5   + L3_5   + A_5   + Y_5   + L1_6  + L2_6  + L3_6  + A_6", 
             "A_7  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6   + L3_6   + A_6   + Y_6   + L1_7  + L2_7  + L3_7",  
             "C_7  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6   + L3_6   + A_6   + Y_6   + L1_7  + L2_7  + L3_7  + A_7", 
             "A_8  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7   + L3_7   + A_7   + Y_7   + L1_8  + L2_8  + L3_8",  
             "C_8  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7   + L3_7   + A_7   + Y_7   + L1_8  + L2_8  + L3_8  + A_8", 
             "A_9  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8   + L3_8   + A_8   + Y_8   + L1_9  + L2_9  + L3_9",  
             "C_9  ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8   + L3_8   + A_8   + Y_8   + L1_9  + L2_9  + L3_9  + A_9", 
             "A_10 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9   + L3_9   + A_9   + Y_9   + L1_10 + L2_10 + L3_10", 
             "C_10 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9   + L3_9   + A_9   + Y_9   + L1_10 + L2_10 + L3_10 + A_10", 
             "A_11 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10  + L3_10  + A_10  + Y_10  + L1_11 + L2_11 + L3_11", 
             "C_11 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10  + L3_10  + A_10  + Y_10  + L1_11 + L2_11 + L3_11 + A_11 ",
             "A_12 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11  + L3_11  + A_11  + Y_11  + L1_12 + L2_12 + L3_12", 
             "C_12 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11  + L3_11  + A_11  + Y_11  + L1_12 + L2_12 + L3_12 + A_12"
             )

###############################

####################
# start simulation #
####################

for(r in 1:R)suppressWarnings(try({
#
ptm.ib <- proc.time()
#
cat(paste("This is simulation run number",r, "\n"))
#
for(n in N){
cat(paste("use n=",n,"; \n"))
simdat <- simcausal::sim(DAG = Dset, n = n)
simdat <-  simdat[,-c(1,6,8,10,11)]
colnames(simdat)[1:7] <- c("b1","b2","b3","b4","b5","b6","b7")
simdat[, c(grep("C",colnames(simdat)))][simdat[, c(grep("C",colnames(simdat)))]==1] <- 2
simdat[, c(grep("C",colnames(simdat)))][simdat[, c(grep("C",colnames(simdat)))]==0] <- 1
simdat[, c(grep("C",colnames(simdat)))][simdat[, c(grep("C",colnames(simdat)))]==2] <- 0
for(j in c(grep("C",colnames(simdat)))){simdat[, j] <- BinaryToCensoring(is.uncensored=simdat[, j])}
#
cd4s_350.1  <- (simdat$L1_1  < 350  | simdat$L2_1  < 0.15 | simdat$L3_1  < -2)
cd4s_350.2  <- (simdat$L1_2  < 350  | simdat$L2_2  < 0.15 | simdat$L3_2  < -2)   | cd4s_350.1
cd4s_350.3  <- (simdat$L1_3  < 350  | simdat$L2_3  < 0.15 | simdat$L3_3  < -2)   | cd4s_350.2
cd4s_350.4  <- (simdat$L1_4  < 350  | simdat$L2_4  < 0.15 | simdat$L3_4  < -2)   | cd4s_350.3
cd4s_350.5  <- (simdat$L1_5  < 350  | simdat$L2_5  < 0.15 | simdat$L3_5  < -2)   | cd4s_350.4
cd4s_350.6  <- (simdat$L1_6  < 350  | simdat$L2_6  < 0.15 | simdat$L3_6  < -2)   | cd4s_350.5
cd4s_350.7  <- (simdat$L1_7  < 350  | simdat$L2_7  < 0.15 | simdat$L3_7  < -2)   | cd4s_350.6
cd4s_350.8  <- (simdat$L1_8  < 350  | simdat$L2_8  < 0.15 | simdat$L3_8  < -2)   | cd4s_350.7
cd4s_350.9  <- (simdat$L1_9  < 350  | simdat$L2_9  < 0.15 | simdat$L3_9  < -2)   | cd4s_350.8
cd4s_350.10 <- (simdat$L1_10 < 350  | simdat$L2_10 < 0.15 | simdat$L3_10 < -2)   | cd4s_350.9
cd4s_350.11 <- (simdat$L1_11 < 350  | simdat$L2_11 < 0.15 | simdat$L3_11 < -2)   | cd4s_350.10
cd4s_350.12 <- (simdat$L1_12 < 350  | simdat$L2_12 < 0.15 | simdat$L3_12 < -2)   | cd4s_350.11
abar350s <- as.matrix(cbind(cd4s_350.1,cd4s_350.2,cd4s_350.3,cd4s_350.4,cd4s_350.5,cd4s_350.6,cd4s_350.7,cd4s_350.8,cd4s_350.9,cd4s_350.10,cd4s_350.11,cd4s_350.12))
abar350s[is.na(abar350s)] <- 0
#
cd4s_750.1  <- (simdat$L1_1  < 750  | simdat$L2_1  < 0.25 | simdat$L3_1  < -2)
cd4s_750.2  <- (simdat$L1_2  < 750  | simdat$L2_2  < 0.25 | simdat$L3_2  < -2)   | cd4s_750.1
cd4s_750.3  <- (simdat$L1_3  < 750  | simdat$L2_3  < 0.25 | simdat$L3_3  < -2)   | cd4s_750.2
cd4s_750.4  <- (simdat$L1_4  < 750  | simdat$L2_4  < 0.25 | simdat$L3_4  < -2)   | cd4s_750.3
cd4s_750.5  <- (simdat$L1_5  < 750  | simdat$L2_5  < 0.25 | simdat$L3_5  < -2)   | cd4s_750.4
cd4s_750.6  <- (simdat$L1_6  < 750  | simdat$L2_6  < 0.25 | simdat$L3_6  < -2)   | cd4s_750.5
cd4s_750.7  <- (simdat$L1_7  < 750  | simdat$L2_7  < 0.25 | simdat$L3_7  < -2)   | cd4s_750.6
cd4s_750.8  <- (simdat$L1_8  < 750  | simdat$L2_8  < 0.25 | simdat$L3_8  < -2)   | cd4s_750.7
cd4s_750.9  <- (simdat$L1_9  < 750  | simdat$L2_9  < 0.25 | simdat$L3_9  < -2)   | cd4s_750.8
cd4s_750.10 <- (simdat$L1_10 < 750  | simdat$L2_10 < 0.25 | simdat$L3_10 < -2)   | cd4s_750.9
cd4s_750.11 <- (simdat$L1_11 < 750  | simdat$L2_11 < 0.25 | simdat$L3_11 < -2)   | cd4s_750.10
cd4s_750.12 <- (simdat$L1_12 < 750  | simdat$L2_12 < 0.25 | simdat$L3_12 < -2)   | cd4s_750.11
abar750s <- as.matrix(cbind(cd4s_750.1,cd4s_750.2,cd4s_750.3,cd4s_750.4,cd4s_750.5,cd4s_750.6,cd4s_750.7,cd4s_750.8,cd4s_750.9,cd4s_750.10,cd4s_750.11,cd4s_750.12))
abar750s[is.na(abar750s)] <- 0
#
for(g in gbounds){
  lc <- 0
  for(l in mylibrary){
#
lc <- lc+1
cat(paste("g=",g,"; l=",mylibrary_names[lc],"; "))
myindex <- c(seq(1:length(allcombinations)))[apply(rbind(grepl(paste(g),allcombinations),grepl(mylibrary_names[lc],allcombinations),grepl(paste(n),allcombinations)),2,all)]
#
write.csv(paste("Simulation number",r, ";", "use n=",n,"; ", "g=",g,"; l=",mylibrary_names[lc],"; "), file=paste(directory,"/runs.csv",sep=""))
#
cat("t=6; ")
tmle_6_always  <- suppressMessages(suppressWarnings(ltmle(simdat[,1:43],
                  Anodes=c(grep("A_",colnames(simdat[,1:43]))),
                  Cnodes=c(grep("C",colnames(simdat[,1:43]))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat[,1:43])),grep("L2_",colnames(simdat[,1:43])),grep("L3_",colnames(simdat[,1:43])))),
                  Ynodes=c(grep("Y_",colnames(simdat[,1:43]))), Yrange=c(-10,10),
                  Qform=myQform[1:6], gform=mygform[1:12], abar=rep(1,6), stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))
                  
tmle_6_750  <-  suppressMessages(suppressWarnings(ltmle(simdat[,1:43],
                  Anodes=c(grep("A_",colnames(simdat[,1:43]))),
                  Cnodes=c(grep("C",colnames(simdat[,1:43]))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat[,1:43])),grep("L2_",colnames(simdat[,1:43])),grep("L3_",colnames(simdat[,1:43])))),
                  Ynodes=c(grep("Y_",colnames(simdat[,1:43]))), Yrange=c(-10,10),
                  Qform=myQform[1:6], gform=mygform[1:12], abar=abar750s[,1:6], stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))
  
tmle_6_350  <-  suppressMessages(suppressWarnings(ltmle(simdat[,1:43],
                  Anodes=c(grep("A_",colnames(simdat[,1:43]))),
                  Cnodes=c(grep("C",colnames(simdat[,1:43]))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat[,1:43])),grep("L2_",colnames(simdat[,1:43])),grep("L3_",colnames(simdat[,1:43])))),
                  Ynodes=c(grep("Y_",colnames(simdat[,1:43]))), Yrange=c(-10,10),
                  Qform=myQform[1:6], gform=mygform[1:12], abar=abar350s[,1:6], stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))

tmle_6_never  <-  suppressMessages(suppressWarnings(ltmle(simdat[,1:43],
                  Anodes=c(grep("A_",colnames(simdat[,1:43]))),
                  Cnodes=c(grep("C",colnames(simdat[,1:43]))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat[,1:43])),grep("L2_",colnames(simdat[,1:43])),grep("L3_",colnames(simdat[,1:43])))),
                  Ynodes=c(grep("Y_",colnames(simdat[,1:43]))), Yrange=c(-10,10),
                  Qform=myQform[1:6], gform=mygform[1:12], abar=rep(0,6), stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))

#
cat("t=12 \n")
tmle_12_always  <- suppressMessages(suppressWarnings(ltmle(simdat,
                  Anodes=c(grep("A_",colnames(simdat))),
                  Cnodes=c(grep("C",colnames(simdat))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat)),grep("L2_",colnames(simdat)),grep("L3_",colnames(simdat)))),
                  Ynodes=c(grep("Y_",colnames(simdat))), Yrange=c(-10,10),
                  Qform=myQform, gform=mygform, abar=rep(1,12), stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))

tmle_12_750  <-  suppressMessages(suppressWarnings(ltmle(simdat,
                  Anodes=c(grep("A_",colnames(simdat))),
                  Cnodes=c(grep("C",colnames(simdat))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat)),grep("L2_",colnames(simdat)),grep("L3_",colnames(simdat)))),
                  Ynodes=c(grep("Y_",colnames(simdat))), Yrange=c(-10,10),
                  Qform=myQform, gform=mygform, abar=abar750s[,1:12], stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))

tmle_12_350  <-  suppressMessages(suppressWarnings(ltmle(simdat,
                  Anodes=c(grep("A_",colnames(simdat))),
                  Cnodes=c(grep("C",colnames(simdat))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat)),grep("L2_",colnames(simdat)),grep("L3_",colnames(simdat)))),
                  Ynodes=c(grep("Y_",colnames(simdat))), Yrange=c(-10,10),
                  Qform=myQform, gform=mygform, abar=abar350s[,1:12], stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))

tmle_12_never  <-  suppressMessages(suppressWarnings(ltmle(simdat,
                  Anodes=c(grep("A_",colnames(simdat))),
                  Cnodes=c(grep("C",colnames(simdat))),
                  Lnodes=sort(c(grep("L1_",colnames(simdat)),grep("L2_",colnames(simdat)),grep("L3_",colnames(simdat)))),
                  Ynodes=c(grep("Y_",colnames(simdat))), Yrange=c(-10,10),
                  Qform=myQform, gform=mygform, abar=rep(0,12), stratify=FALSE, variance.method="ic",
                  SL.library=l, estimate.time=F, gbounds=c(g,1))))

# 
results_always_6[r,myindex] <-  tmle_6_always$est[1]
results_CP_always_6[r,myindex] <-  as.numeric((summary(tmle_6_always)$treatment$CI[1] < true[7,1]) & (summary(tmle_6_always)$treatment$CI[2] > true[7,1]))
results_always_12[r,myindex] <-  tmle_12_always$est[1]
results_CP_always_12[r,myindex] <-  as.numeric((summary(tmle_12_always)$treatment$CI[1] < true[13,1]) & (summary(tmle_12_always)$treatment$CI[2] > true[13,1]))

results_750_6[r,myindex] <-  tmle_6_750$est[1]
results_CP_750_6[r,myindex] <-  as.numeric((summary(tmle_6_750)$treatment$CI[1] < true[7,2]) & (summary(tmle_6_750)$treatment$CI[2] > true[7,2]))
results_750_12[r,myindex] <-  tmle_12_750$est[1]
results_CP_750_12[r,myindex] <-  as.numeric((summary(tmle_12_750)$treatment$CI[1] < true[13,2]) & (summary(tmle_12_750)$treatment$CI[2] > true[13,2]))

results_350_6[r,myindex] <-  tmle_6_350$est[1]
results_CP_350_6[r,myindex] <-  as.numeric((summary(tmle_6_350)$treatment$CI[1] < true[7,3]) & (summary(tmle_6_350)$treatment$CI[2] > true[7,3]))
results_350_12[r,myindex] <-  tmle_12_350$est[1]
results_CP_350_12[r,myindex] <-  as.numeric((summary(tmle_12_350)$treatment$CI[1] < true[13,3]) & (summary(tmle_12_350)$treatment$CI[2] > true[13,3]))

results_never_6[r,myindex] <-  tmle_6_never$est[1]
results_CP_never_6[r,myindex] <-  as.numeric((summary(tmle_6_never)$treatment$CI[1] < true[7,4]) & (summary(tmle_6_never)$treatment$CI[2] > true[7,4]))
results_never_12[r,myindex] <-  tmle_12_never$est[1]
results_CP_never_12[r,myindex] <-  as.numeric((summary(tmle_12_never)$treatment$CI[1] < true[13,4]) & (summary(tmle_12_never)$treatment$CI[2] > true[13,4]))

results_ss_always_12[r,myindex] <- tmle_12_always$fit$Qstar$Y_12$n
results_ss_750_12[r,myindex] <- tmle_12_750$fit$Qstar$Y_12$n
results_ss_350_12[r,myindex] <- tmle_12_350$fit$Qstar$Y_12$n
results_ss_never_12[r,myindex] <- tmle_12_never$fit$Qstar$Y_12$n

results_ss_always_6[r,myindex] <- tmle_6_always$fit$Qstar$Y_6$n
results_ss_750_6[r,myindex] <- tmle_6_750$fit$Qstar$Y_6$n
results_ss_350_6[r,myindex] <- tmle_6_350$fit$Qstar$Y_6$n
results_ss_never_6[r,myindex] <- tmle_6_never$fit$Qstar$Y_6$n

if(lc==length(mylibrary) & g == gbounds[1] & n==N[length(N)])try({
results_SL_Q_always_12[r,]    <-  apply(t(matrix(unlist(lapply(tmle_12_always$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)    
results_SL_Q_750_12[r,]       <-  apply(t(matrix(unlist(lapply(tmle_12_750$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)       
results_SL_Q_350_12[r,]       <-  apply(t(matrix(unlist(lapply(tmle_12_350$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)       
results_SL_Q_never_12[r,]     <-  apply(t(matrix(unlist(lapply(tmle_12_never$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)     
                      
results_SL_Q_always_6[r,]     <-  apply(t(matrix(unlist(lapply(tmle_6_always$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
results_SL_Q_750_6[r,]        <-  apply(t(matrix(unlist(lapply(tmle_6_750$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
results_SL_Q_350_6[r,]        <-  apply(t(matrix(unlist(lapply(tmle_6_350$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
results_SL_Q_never_6[r,]      <-  apply(t(matrix(unlist(lapply(tmle_6_never$fit$Q, ew)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)

results_SL_g_always_12[r,]    <-  apply(t(matrix(unlist(lapply(tmle_12_always$fit$g, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)    
results_SL_g_750_12[r,]       <-  apply(t(matrix(unlist(lapply(tmle_12_750$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)       
results_SL_g_350_12[r,]       <-  apply(t(matrix(unlist(lapply(tmle_12_350$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)       
results_SL_g_never_12[r,]     <-  apply(t(matrix(unlist(lapply(tmle_12_never$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)     
                      
results_SL_g_always_6[r,]     <-  apply(t(matrix(unlist(lapply(tmle_6_always$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
results_SL_g_750_6[r,]        <-  apply(t(matrix(unlist(lapply(tmle_6_750$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
results_SL_g_350_6[r,]        <-  apply(t(matrix(unlist(lapply(tmle_6_350$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
results_SL_g_never_6[r,]      <-  apply(t(matrix(unlist(lapply(tmle_6_never$fit$Q, ewg)),nrow=length(l$Q),dimnames=list(l[[1]],NULL))),2,mean)
})

#
}}

#
tempresults <- list(results_always_6,results_always_12,results_750_6,results_750_12,results_350_6,results_350_12,results_never_6,results_never_12,
                    results_CP_always_6,results_CP_always_12,results_CP_750_6,results_CP_750_12,results_CP_350_6,results_CP_350_12,results_CP_never_6,results_CP_never_12,                                                                                                                                      
                    results_ss_always_12,results_ss_750_12,results_ss_350_12,results_ss_never_12,
                    results_SL_Q_always_12, results_SL_Q_750_12, results_SL_Q_350_12, results_SL_Q_never_12, results_SL_Q_always_6, results_SL_Q_750_6, results_SL_Q_350_6, results_SL_Q_never_6, 
                    results_SL_g_always_12, results_SL_g_750_12, results_SL_g_350_12, results_SL_g_never_12, results_SL_g_always_6, results_SL_g_750_6, results_SL_g_350_6, results_SL_g_never_6
                    )
save(tempresults, file=paste(directory,"/tempresults.Rdata",sep=""))
#
}
#
ptm.ib2 <- proc.time()
ibt <- round(((ptm.ib2-ptm.ib)/60)[1],digits=2)
if(r==1){cat(paste("The simulation will run for about another", R*ibt-ibt, "minutes \n"))}
timetogo <-  (paste("The simulation will run for about another", R*ibt-ibt, "minutes \n"))
write.csv(timetogo, file=paste(directory,"/timetogo.csv",sep=""))
#
}))

BIAS6     <- matrix(NA,nrow=4,ncol=length(allcombinations),dimnames=list(c("always","750","350","never"),allcombinations))
COVERAGE6 <- matrix(NA,nrow=4,ncol=length(allcombinations),dimnames=list(c("always","750","350","never"),allcombinations))
BIAS12     <- matrix(NA,nrow=4,ncol=length(allcombinations),dimnames=list(c("always","750","350","never"),allcombinations))
COVERAGE12 <- matrix(NA,nrow=4,ncol=length(allcombinations),dimnames=list(c("always","750","350","never"),allcombinations))

BIAS6[1,] <- apply(results_always_6,2,mean)-true[7,1]
BIAS6[2,] <- apply(results_750_6,2,mean)-true[7,2]
BIAS6[3,] <- apply(results_350_6,2,mean)-true[7,3]
BIAS6[4,] <- apply(results_never_6,2,mean)-true[7,4]
BIAS6 <- abs(BIAS6)
BIAS12[1,] <- apply(results_always_12,2,mean)-true[13,1]
BIAS12[2,] <- apply(results_750_12,2,mean)-true[13,2]
BIAS12[3,] <- apply(results_350_12,2,mean)-true[13,3]
BIAS12[4,] <- apply(results_never_12,2,mean)-true[13,4]
BIAS12 <- abs(BIAS12)

COVERAGE6[1,] <- apply(results_CP_always_6,2,mean)
COVERAGE6[2,] <- apply(results_CP_750_6,2,mean)
COVERAGE6[3,] <- apply(results_CP_350_6,2,mean)
COVERAGE6[4,] <- apply(results_CP_never_6,2,mean)
COVERAGE12[1,] <- apply(results_CP_always_12,2,mean)
COVERAGE12[2,] <- apply(results_CP_750_12,2,mean)
COVERAGE12[3,] <- apply(results_CP_350_12,2,mean)
COVERAGE12[4,] <- apply(results_CP_never_12,2,mean)

RESULTS12 <-  data.frame(cbind(rbind(expand.grid(N, gbounds, mylibrary_names),expand.grid(N, gbounds, mylibrary_names),expand.grid(N, gbounds, mylibrary_names),expand.grid(N, gbounds, mylibrary_names)),
                   c(rep("(1) always",length(allcombinations)),rep("(2) 750/25%/-2",length(allcombinations)),rep("(3) 350/15%/-2",length(allcombinations)),rep("(4) never",length(allcombinations))),
                   c(t(BIAS12)), c(t(COVERAGE12)), rep(12,length(allcombinations)*4)
                   ))
colnames(RESULTS12) <- c("n","g","Learner","Intervention","Bias","Coverage","t") 
RESULTS6 <-  data.frame(cbind(rbind(expand.grid(N, gbounds, mylibrary_names),expand.grid(N, gbounds, mylibrary_names),expand.grid(N, gbounds, mylibrary_names),expand.grid(N, gbounds, mylibrary_names)),
                   c(rep("(1) always",length(allcombinations)),rep("(2) 750/25%/-2",length(allcombinations)),rep("(3) 350/15%/-2",length(allcombinations)),rep("(4) never",length(allcombinations))),
                   c(t(BIAS6)), c(t(COVERAGE6)), rep(6,length(allcombinations)*4)
                   ))
colnames(RESULTS6) <- c("n","g","Learner","Intervention","Bias","Coverage","t") 



pdf(file=paste(directory,"/Figure4_Paper_t6_and_t12.pdf",sep=""), width=9)  # reported in the paper for t=12
contour1 <- ggplot(RESULTS6, aes(x=g, y=n))
contour2 <- contour1 + geom_raster(aes(fill = Bias), interpolate = TRUE) + scale_fill_gradient(low="green", high="red")  + facet_grid(Learner ~ Intervention)
contour3 <- contour2 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour3 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="|Bias|")) + ggtitle("Bias at t=6") )
contour4 <- ggplot(RESULTS12, aes(x=g, y=n))
contour5 <- contour4 + geom_raster(aes(fill = Bias), interpolate = TRUE) + scale_fill_gradient(low="green", high="red")  + facet_grid(Learner ~ Intervention)
contour6 <- contour5 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour6 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="|Bias|")) + ggtitle("Bias at t=12") )
contour7 <- ggplot(RESULTS6, aes(x=g, y=n))
contour8 <- contour7 + geom_raster(aes(fill = Coverage), interpolate = TRUE) + scale_fill_gradient2(midpoint=0.95, mid="green")  + facet_grid(Learner ~ Intervention)
contour9 <- contour8 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour9 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="Coverage")) + ggtitle("Coverage Probability for t=6") )
contour10 <- ggplot(RESULTS12, aes(x=g, y=n))
contour11 <- contour10 + geom_raster(aes(fill = Coverage), interpolate = TRUE) + scale_fill_gradient2(midpoint=0.95, mid="green")  + facet_grid(Learner ~ Intervention)
contour12 <- contour11 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour12 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="Coverage")) + ggtitle("Coverage Probability for t=12") )
dev.off()

pdf(file=paste(directory,"/Figure4_Paper_noshading_t6_and_t12.pdf",sep=""), width=9)    # as above, but without interpolating colours
contour1 <- ggplot(RESULTS6, aes(x=g, y=n))
contour2 <- contour1  + geom_tile(aes(fill=Bias)) + scale_fill_gradient(low="green", high="red")  + facet_grid(Learner ~ Intervention)
contour3 <- contour2 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour3 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="|Bias|")) + ggtitle("Bias at t=6") )
contour4 <- ggplot(RESULTS12, aes(x=g, y=n))
contour5 <- contour4 + geom_tile(aes(fill=Bias)) + scale_fill_gradient(low="green", high="red")  + facet_grid(Learner ~ Intervention)
contour6 <- contour5 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour6 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="|Bias|")) + ggtitle("Bias at t=12") )
contour7 <- ggplot(RESULTS6, aes(x=g, y=n))
contour8 <- contour7 + geom_tile(aes(fill=Coverage)) + scale_fill_gradient2(midpoint=0.95, mid="green")  + facet_grid(Learner ~ Intervention)
contour9 <- contour8 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour9 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="Coverage")) + ggtitle("Coverage Probability for t=6") )
contour10 <- ggplot(RESULTS12, aes(x=g, y=n))
contour11 <- contour10 + geom_tile(aes(fill=Coverage)) + scale_fill_gradient2(midpoint=0.95, mid="green")  + facet_grid(Learner ~ Intervention)
contour12 <- contour11 + theme_bw()  + scale_x_continuous("truncation level",breaks=c(0,0.02,0.04)) + scale_y_continuous("sample size")  + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
plot(contour12 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="Coverage")) + ggtitle("Coverage Probability for t=12") )
dev.off()

SLweights <- c(apply(results_SL_Q_always_12,2,mean), apply(results_SL_Q_750_12,2,mean), apply(results_SL_Q_350_12,2,mean),
      apply(results_SL_Q_never_12,2,mean), apply(results_SL_Q_always_6,2,mean), apply(results_SL_Q_750_6,2,mean),
      apply(results_SL_Q_350_6,2,mean), apply(results_SL_Q_never_6,2,mean), apply(results_SL_g_always_12,2,mean),
      apply(results_SL_g_750_12,2,mean), apply(results_SL_g_350_12,2,mean), apply(results_SL_g_never_12,2,mean),
      apply(results_SL_g_always_6,2,mean), apply(results_SL_g_750_6,2,mean), apply(results_SL_g_350_6,2,mean),
      apply(results_SL_g_never_6,2,mean))
lw <-  length(mylibrary[length(mylibrary)][[1]]$Q)
SLW <- data.frame(SLweights, names(SLweights), c(rep("Q-model",length(SLweights)/2),rep("g-model",length(SLweights)/2)), c(rep(12,length(SLweights)/4),rep(6,length(SLweights)/4),rep(12,length(SLweights)/4),rep(6,length(SLweights)/4)),
                  rep(c(rep("1) always",lw),rep("2) 750/25/-2",lw),rep("3) 350/15/-2",lw),rep("4) never",lw)),4))
colnames(SLW) <- c("weight","learner","model","time", "Intervention")
 
pdf(file=paste(directory,"/SL_summary.pdf",sep=""), width=9) # super learner summary, not reported in paper
slp1 <- ggplot(SLW, aes(x=model, y=weight))     
slp2 <- slp1  + geom_col(aes(fill = learner)) + facet_grid(Intervention ~ time) + theme_bw()
slp3 <- slp2 + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
slp4 <- slp3 + guides(fill = guide_legend(keywidth = 2, keyheight = 2, title="Learner")) + ggtitle(paste("Weights of super learners for n=",N[length(N)], "and g=", gbounds[1])) 
plot(slp4)
dev.off()

sams <- c(t(rbind(apply(results_ss_always_12,2,mean),apply(results_ss_750_12,2,mean),apply(results_ss_350_12,2,mean),apply(results_ss_never_12,2,mean)))[1:3,])
SS <- data.frame(sams,rep(N,4),c("1) always","1) always","1) always","2) 750/25-2","2) 750/25-2","2) 750/25-2","3) 350/15/-2","3) 350/15/-2","3) 350/15/-2","4) never","4) never","4) never"))
colnames(SS) <- c("sample_size","original","intervention")

pdf(file=paste(directory,"/sample_size.pdf",sep=""), width=7.5) # "usable" sample size, i.e. how many remain uncensored and follow rule of interest; not reported in paper
ss1 <- ggplot(SS, aes(x=intervention, y=sample_size)) 
ss2 <- ss1 + geom_col() + facet_grid(original ~ .) + theme_bw()
ss3 <- ss2 + theme(axis.title.x = element_text(size=16), axis.text.x = element_text(size=16),axis.title.y = element_text(size=16, angle = 90), axis.text.y = element_text(size=16), legend.text =  element_text(size=16), legend.title =  element_text(size=16, face = "bold", hjust = 0),legend.position =   "right")
ss4 <- ss3 + ggtitle("Usable sample size at t=12, stratified by intervention and baseline sample size") + scale_y_continuous("sample size at t=12") 
plot(ss4)
dev.off()

#
results <- list(results_always_6,results_always_12,results_750_6,results_750_12,results_350_6,results_350_12,results_never_6,results_never_12,
                    results_CP_always_6,results_CP_always_12,results_CP_750_6,results_CP_750_12,results_CP_350_6,results_CP_350_12,results_CP_never_6,results_CP_never_12,                                                                                                                                      
                    results_ss_always_12,results_ss_750_12,results_ss_350_12,results_ss_never_12,
                    results_SL_Q_always_12, results_SL_Q_750_12, results_SL_Q_350_12, results_SL_Q_never_12, results_SL_Q_always_6, results_SL_Q_750_6, results_SL_Q_350_6, results_SL_Q_never_6, 
                    results_SL_g_always_12, results_SL_g_750_12, results_SL_g_350_12, results_SL_g_never_12, results_SL_g_always_6, results_SL_g_750_6, results_SL_g_350_6, results_SL_g_never_6
                    )
save(results, file=paste(directory,"/results.Rdata",sep=""))
#

# How long did it take?
ptm2 <-  proc.time()
simulationsdauer <- ptm2-ptm
simulationsdauer <- (simulationsdauer/60)
simulationsdauer <- round(simulationsdauer[1]/60, digits=2)
cat(paste("The simulation time was", simulationsdauer, "hours \n"))



###############################
# Calculate support - Table 3 #
###############################

N=1000000
simdat <- simcausal::sim(DAG = Dset, n = N) 
simdat <-  simdat[,-c(1,6,8,10,11)]
colnames(simdat)[1:7] <- c("b1","b2","b3","b4","b5","b6","b7")


fillup<-function(myvec){
for(i in 2:length(myvec)){
if(is.na(myvec)[i]){myvec[i]<-myvec[i-1]} 
}
return(myvec)
}
eligible <- function(myvec){
myvec[is.na(myvec)] <- 0
i <- 1
elig <- length(myvec)+1
while(i <= length(myvec)){
if(myvec[i]==1){elig<-i
i <- length(myvec)+1
}else{i <- i+1}
}
return(elig)
}
actuallytreated <- simdat[,c(grep("A_",colnames(simdat)))]

d11  <- (simdat$L1_1  < 350  | simdat$L2_1  < 0.15 | simdat$L3_1  < -2)   
d12  <- (simdat$L1_2  < 350  | simdat$L2_2  < 0.15 | simdat$L3_2  < -2) 
d13  <- (simdat$L1_3  < 350  | simdat$L2_3  < 0.15 | simdat$L3_3  < -2) 
d14  <- (simdat$L1_4  < 350  | simdat$L2_4  < 0.15 | simdat$L3_4  < -2) 
d15  <- (simdat$L1_5  < 350  | simdat$L2_5  < 0.15 | simdat$L3_5  < -2) 
d16  <- (simdat$L1_6  < 350  | simdat$L2_6  < 0.15 | simdat$L3_6  < -2) 
d17  <- (simdat$L1_7  < 350  | simdat$L2_7  < 0.15 | simdat$L3_7  < -2) 
d18  <- (simdat$L1_8  < 350  | simdat$L2_8  < 0.15 | simdat$L3_8  < -2) 
d19  <- (simdat$L1_9  < 350  | simdat$L2_9  < 0.15 | simdat$L3_9  < -2) 
d110 <- (simdat$L1_10  < 350  | simdat$L2_10  < 0.15 | simdat$L3_10  < -2)
d111 <- (simdat$L1_11  < 350  | simdat$L2_11  < 0.15 | simdat$L3_11  < -2)
d112 <- (simdat$L1_12  < 350  | simdat$L2_12  < 0.15 | simdat$L3_12  < -2)
d1rule <- apply(cbind(d11,d12,d13,d14,d15,d16,d17,d18,d19,d110,d111,d112),2,as.numeric)
dsum <- t(apply(cbind(d11,d12,d13,d14,d15,d16,d17,d18,d19,d110,d111,d112),1,cumsum))
dsum[dsum!=1] <- 0
eg <- apply(dsum,1,eligible)
shouldbetreated_d1 <- matrix(NA,nrow=N,ncol=13)
for(i in 1:N){shouldbetreated_d1[i,eg[i]]<-1}
shouldbetreated_d1[is.na(shouldbetreated_d1)[,1],1]<-0
shouldbetreated_d1 <- t(apply(shouldbetreated_d1,1,fillup))
posp <- actuallytreated==shouldbetreated_d1 
posp1  <- as.numeric(posp[,1])
posp2  <- as.numeric(posp[,2])       
posp3  <- as.numeric(posp[,3])       
posp4  <- as.numeric(posp[,4])       
posp5  <- as.numeric(posp[,5])       
posp6  <- as.numeric(posp[,6])       
posp7  <- as.numeric(posp[,7])       
posp8  <- as.numeric(posp[,8])       
posp9  <- as.numeric(posp[,9])       
posp10 <- as.numeric(posp[,10])       
posp11 <- as.numeric(posp[,11])       
posp12 <- as.numeric(posp[,12])       
p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- p7 <- p8 <- p9 <- p10 <- p11 <- p12 <- rep(NA,length(posp1))
evdat <- cbind(simdat,posp1,posp2,posp3,posp4,posp5,posp6,posp7,posp8,posp9,posp10,posp11,posp12)

p1[is.na(evdat$posp1)==F]  <- predict(glm(posp1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1, family=binomial,data=evdat),type="response")
p2[is.na(evdat$posp2)==F & evdat$posp1==1]  <- predict(glm(posp2 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1   + Y_1   + L1_2  + L2_2  + L3_2, family=binomial,data=evdat[evdat$posp1==1,]),type="response")
p3[is.na(evdat$posp3)==F & evdat$posp2==1]  <- predict(glm(posp3 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2   + L3_2   + A_2   + Y_2   + L1_3  + L2_3  + L3_3, family=binomial,data=evdat[evdat$posp2==1,]),type="response")
p4[is.na(evdat$posp4)==F & evdat$posp3==1]  <- predict(glm(posp4 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3   + L3_3   + A_3   + Y_3   + L1_4  + L2_4  + L3_4, family=binomial,data=evdat[evdat$posp3==1,]),type="response")
p5[is.na(evdat$posp5)==F & evdat$posp4==1]  <- predict(glm(posp5 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4   + L3_4   + A_4   + Y_4   + L1_5  + L2_5  + L3_5, family=binomial,data=evdat[evdat$posp4==1,]),type="response")
p6[is.na(evdat$posp6)==F & evdat$posp5==1]  <- predict(glm(posp6 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5   + L3_5   + A_5   + Y_5   + L1_6  + L2_6  + L3_6, family=binomial,data=evdat[evdat$posp5==1,]),type="response")
p7[is.na(evdat$posp7)==F & evdat$posp6==1]  <- predict(glm(posp7 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6   + L3_6   + A_6   + Y_6   + L1_7  + L2_7  + L3_7, family=binomial,data=evdat[evdat$posp6==1,]),type="response")
p8[is.na(evdat$posp8)==F & evdat$posp7==1]  <- predict(glm(posp8 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7   + L3_7   + A_7   + Y_7   + L1_8  + L2_8  + L3_8, family=binomial,data=evdat[evdat$posp7==1,]),type="response")
p9[is.na(evdat$posp9)==F & evdat$posp8==1]  <- predict(glm(posp9 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8   + L3_8   + A_8   + Y_8   + L1_9  + L2_9  + L3_9, family=binomial,data=evdat[evdat$posp8==1,]),type="response")
p10[is.na(evdat$posp10)==F & evdat$posp9==1] <- predict(glm(posp10 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9   + L3_9   + A_9   + Y_9   + L1_10 + L2_10 + L3_10, family=binomial,data=evdat[evdat$posp9==1,]),type="response")
p11[is.na(evdat$posp11)==F & evdat$posp10==1] <- predict(glm(posp11 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10  + L3_10  + A_10  + Y_10  + L1_11 + L2_11 + L3_11, family=binomial,data=evdat[evdat$posp10==1,]),type="response")
p12[is.na(evdat$posp12)==F & evdat$posp11==1] <- predict(glm(posp12 ~  b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11  + L3_11  + A_11  + Y_11  + L1_12 + L2_12 + L3_12, family=binomial,data=evdat[evdat$posp11==1,]),type="response")

p350 <- t(apply(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),1,cumprod))

#
d21 <- (simdat$L1_1  < 750  | simdat$L2_1  < 0.25 | simdat$L3_1  < -2)   
d22 <- (simdat$L1_2  < 750  | simdat$L2_2  < 0.25 | simdat$L3_2  < -2) 
d23 <- (simdat$L1_3  < 750  | simdat$L2_3  < 0.25 | simdat$L3_3  < -2) 
d24 <- (simdat$L1_4  < 750  | simdat$L2_4  < 0.25 | simdat$L3_4  < -2) 
d25 <- (simdat$L1_5  < 750  | simdat$L2_5  < 0.25 | simdat$L3_5  < -2) 
d26 <- (simdat$L1_6  < 750  | simdat$L2_6  < 0.25 | simdat$L3_6  < -2) 
d27 <- (simdat$L1_7  < 750  | simdat$L2_7  < 0.25 | simdat$L3_7  < -2) 
d28 <- (simdat$L1_8  < 750  | simdat$L2_8  < 0.25 | simdat$L3_8  < -2) 
d29 <- (simdat$L1_9  < 750  | simdat$L2_9  < 0.25 | simdat$L3_9  < -2) 
d210 <- (simdat$L1_10  < 750  | simdat$L2_10  < 0.25 | simdat$L3_10  < -2)
d211 <- (simdat$L1_11  < 750  | simdat$L2_11  < 0.25 | simdat$L3_11  < -2)
d212 <- (simdat$L1_12  < 750  | simdat$L2_12  < 0.25 | simdat$L3_12  < -2)
dsum2 <- t(apply(cbind(d21,d22,d23,d24,d25,d26,d27,d28,d29,d210,d211,d212),1,cumsum))
dsum2[dsum2!=1] <- 0
eg <- apply(dsum2,1,eligible)
shouldbetreated_d2 <- matrix(NA,nrow=N,ncol=13)
for(i in 1:N){shouldbetreated_d2[i,eg[i]]<-1}
shouldbetreated_d2[is.na(shouldbetreated_d2)[,1],1]<-0
shouldbetreated_d2 <- t(apply(shouldbetreated_d2,1,fillup))
posp <- actuallytreated==shouldbetreated_d2 
posp1  <- as.numeric(posp[,1])
posp2  <- as.numeric(posp[,2])       
posp3  <- as.numeric(posp[,3])       
posp4  <- as.numeric(posp[,4])       
posp5  <- as.numeric(posp[,5])       
posp6  <- as.numeric(posp[,6])       
posp7  <- as.numeric(posp[,7])       
posp8  <- as.numeric(posp[,8])       
posp9  <- as.numeric(posp[,9])       
posp10 <- as.numeric(posp[,10])       
posp11 <- as.numeric(posp[,11])       
posp12 <- as.numeric(posp[,12])       
p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- p7 <- p8 <- p9 <- p10 <- p11 <- p12 <- rep(NA,length(posp1))
evdat <- cbind(simdat,posp1,posp2,posp3,posp4,posp5,posp6,posp7,posp8,posp9,posp10,posp11,posp12)

p1[is.na(evdat$posp1)==F]  <- predict(glm(posp1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1, family=binomial,data=evdat),type="response")
p2[is.na(evdat$posp2)==F & evdat$posp1==1]  <- predict(glm(posp2 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1   + Y_1   + L1_2  + L2_2  + L3_2, family=binomial,data=evdat[evdat$posp1==1,]),type="response")
p3[is.na(evdat$posp3)==F & evdat$posp2==1]  <- predict(glm(posp3 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2   + L3_2   + A_2   + Y_2   + L1_3  + L2_3  + L3_3, family=binomial,data=evdat[evdat$posp2==1,]),type="response")
p4[is.na(evdat$posp4)==F & evdat$posp3==1]  <- predict(glm(posp4 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3   + L3_3   + A_3   + Y_3   + L1_4  + L2_4  + L3_4, family=binomial,data=evdat[evdat$posp3==1,]),type="response")
p5[is.na(evdat$posp5)==F & evdat$posp4==1]  <- predict(glm(posp5 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4   + L3_4   + A_4   + Y_4   + L1_5  + L2_5  + L3_5, family=binomial,data=evdat[evdat$posp4==1,]),type="response")
p6[is.na(evdat$posp6)==F & evdat$posp5==1]  <- predict(glm(posp6 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5   + L3_5   + A_5   + Y_5   + L1_6  + L2_6  + L3_6, family=binomial,data=evdat[evdat$posp5==1,]),type="response")
p7[is.na(evdat$posp7)==F & evdat$posp6==1]  <- predict(glm(posp7 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6   + L3_6   + A_6   + Y_6   + L1_7  + L2_7  + L3_7, family=binomial,data=evdat[evdat$posp6==1,]),type="response")
p8[is.na(evdat$posp8)==F & evdat$posp7==1]  <- predict(glm(posp8 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7   + L3_7   + A_7   + Y_7   + L1_8  + L2_8  + L3_8, family=binomial,data=evdat[evdat$posp7==1,]),type="response")
p9[is.na(evdat$posp9)==F & evdat$posp8==1]  <- predict(glm(posp9 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8   + L3_8   + A_8   + Y_8   + L1_9  + L2_9  + L3_9, family=binomial,data=evdat[evdat$posp8==1,]),type="response")
p10[is.na(evdat$posp10)==F & evdat$posp9==1] <- predict(glm(posp10 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9   + L3_9   + A_9   + Y_9   + L1_10 + L2_10 + L3_10, family=binomial,data=evdat[evdat$posp9==1,]),type="response")
p11[is.na(evdat$posp11)==F & evdat$posp10==1] <- predict(glm(posp11 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10  + L3_10  + A_10  + Y_10  + L1_11 + L2_11 + L3_11, family=binomial,data=evdat[evdat$posp10==1,]),type="response")
p12[is.na(evdat$posp12)==F & evdat$posp11==1] <- predict(glm(posp12 ~  b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11  + L3_11  + A_11  + Y_11  + L1_12 + L2_12 + L3_12, family=binomial,data=evdat[evdat$posp11==1,]),type="response")

p750 <- t(apply(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),1,cumprod))


###############
shouldbetreated_d3 <- matrix(1,ncol=12,nrow=N)
posp <- actuallytreated==shouldbetreated_d3 
posp1  <- as.numeric(posp[,1])
posp2  <- as.numeric(posp[,2])       
posp3  <- as.numeric(posp[,3])       
posp4  <- as.numeric(posp[,4])       
posp5  <- as.numeric(posp[,5])       
posp6  <- as.numeric(posp[,6])       
posp7  <- as.numeric(posp[,7])       
posp8  <- as.numeric(posp[,8])       
posp9  <- as.numeric(posp[,9])       
posp10 <- as.numeric(posp[,10])       
posp11 <- as.numeric(posp[,11])       
posp12 <- as.numeric(posp[,12])       
p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- p7 <- p8 <- p9 <- p10 <- p11 <- p12 <- rep(NA,length(posp1))
evdat <- cbind(simdat,posp1,posp2,posp3,posp4,posp5,posp6,posp7,posp8,posp9,posp10,posp11,posp12)

p1[is.na(evdat$posp1)==F]  <- predict(glm(posp1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1, family=binomial,data=evdat),type="response")
p2[is.na(evdat$posp2)==F & evdat$posp1==1]  <- predict(glm(posp2 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1   + Y_1   + L1_2  + L2_2  + L3_2, family=binomial,data=evdat[evdat$posp1==1,]),type="response")
p3[is.na(evdat$posp3)==F & evdat$posp2==1]  <- predict(glm(posp3 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2   + L3_2   + A_2   + Y_2   + L1_3  + L2_3  + L3_3, family=binomial,data=evdat[evdat$posp2==1,]),type="response")
p4[is.na(evdat$posp4)==F & evdat$posp3==1]  <- predict(glm(posp4 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3   + L3_3   + A_3   + Y_3   + L1_4  + L2_4  + L3_4, family=binomial,data=evdat[evdat$posp3==1,]),type="response")
p5[is.na(evdat$posp5)==F & evdat$posp4==1]  <- predict(glm(posp5 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4   + L3_4   + A_4   + Y_4   + L1_5  + L2_5  + L3_5, family=binomial,data=evdat[evdat$posp4==1,]),type="response")
p6[is.na(evdat$posp6)==F & evdat$posp5==1]  <- predict(glm(posp6 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5   + L3_5   + A_5   + Y_5   + L1_6  + L2_6  + L3_6, family=binomial,data=evdat[evdat$posp5==1,]),type="response")
p7[is.na(evdat$posp7)==F & evdat$posp6==1]  <- predict(glm(posp7 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6   + L3_6   + A_6   + Y_6   + L1_7  + L2_7  + L3_7, family=binomial,data=evdat[evdat$posp6==1,]),type="response")
p8[is.na(evdat$posp8)==F & evdat$posp7==1]  <- predict(glm(posp8 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7   + L3_7   + A_7   + Y_7   + L1_8  + L2_8  + L3_8, family=binomial,data=evdat[evdat$posp7==1,]),type="response")
p9[is.na(evdat$posp9)==F & evdat$posp8==1]  <- predict(glm(posp9 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8   + L3_8   + A_8   + Y_8   + L1_9  + L2_9  + L3_9, family=binomial,data=evdat[evdat$posp8==1,]),type="response")
p10[is.na(evdat$posp10)==F & evdat$posp9==1] <- predict(glm(posp10 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9   + L3_9   + A_9   + Y_9   + L1_10 + L2_10 + L3_10, family=binomial,data=evdat[evdat$posp9==1,]),type="response")
p11[is.na(evdat$posp11)==F & evdat$posp10==1] <- predict(glm(posp11 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10  + L3_10  + A_10  + Y_10  + L1_11 + L2_11 + L3_11, family=binomial,data=evdat[evdat$posp10==1,]),type="response")
p12[is.na(evdat$posp12)==F & evdat$posp11==1] <- predict(glm(posp12 ~  b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11  + L3_11  + A_11  + Y_11  + L1_12 + L2_12 + L3_12, family=binomial,data=evdat[evdat$posp11==1,]),type="response")

palways <- t(apply(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),1,cumprod))

####
shouldbetreated_d4 <- matrix(0,ncol=12,nrow=N)
posp <- actuallytreated==shouldbetreated_d4 
posp1  <- as.numeric(posp[,1])
posp2  <- as.numeric(posp[,2])       
posp3  <- as.numeric(posp[,3])       
posp4  <- as.numeric(posp[,4])       
posp5  <- as.numeric(posp[,5])       
posp6  <- as.numeric(posp[,6])       
posp7  <- as.numeric(posp[,7])       
posp8  <- as.numeric(posp[,8])       
posp9  <- as.numeric(posp[,9])       
posp10 <- as.numeric(posp[,10])       
posp11 <- as.numeric(posp[,11])       
posp12 <- as.numeric(posp[,12])       
p1 <- p2 <- p3 <- p4 <- p5 <- p6 <- p7 <- p8 <- p9 <- p10 <- p11 <- p12 <- rep(NA,length(posp1))
evdat <- cbind(simdat,posp1,posp2,posp3,posp4,posp5,posp6,posp7,posp8,posp9,posp10,posp11,posp12)

p1[is.na(evdat$posp1)==F]  <- predict(glm(posp1 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1, family=binomial,data=evdat),type="response")
p2[is.na(evdat$posp2)==F & evdat$posp1==1]  <- predict(glm(posp2 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_1  + L2_1   + L3_1   + A_1   + Y_1   + L1_2  + L2_2  + L3_2, family=binomial,data=evdat[evdat$posp1==1,]),type="response")
p3[is.na(evdat$posp3)==F & evdat$posp2==1]  <- predict(glm(posp3 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_2  + L2_2   + L3_2   + A_2   + Y_2   + L1_3  + L2_3  + L3_3, family=binomial,data=evdat[evdat$posp2==1,]),type="response")
p4[is.na(evdat$posp4)==F & evdat$posp3==1]  <- predict(glm(posp4 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_3  + L2_3   + L3_3   + A_3   + Y_3   + L1_4  + L2_4  + L3_4, family=binomial,data=evdat[evdat$posp3==1,]),type="response")
p5[is.na(evdat$posp5)==F & evdat$posp4==1]  <- predict(glm(posp5 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_4  + L2_4   + L3_4   + A_4   + Y_4   + L1_5  + L2_5  + L3_5, family=binomial,data=evdat[evdat$posp4==1,]),type="response")
p6[is.na(evdat$posp6)==F & evdat$posp5==1]  <- predict(glm(posp6 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_5  + L2_5   + L3_5   + A_5   + Y_5   + L1_6  + L2_6  + L3_6, family=binomial,data=evdat[evdat$posp5==1,]),type="response")
p7[is.na(evdat$posp7)==F & evdat$posp6==1]  <- predict(glm(posp7 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_6  + L2_6   + L3_6   + A_6   + Y_6   + L1_7  + L2_7  + L3_7, family=binomial,data=evdat[evdat$posp6==1,]),type="response")
p8[is.na(evdat$posp8)==F & evdat$posp7==1]  <- predict(glm(posp8 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_7  + L2_7   + L3_7   + A_7   + Y_7   + L1_8  + L2_8  + L3_8, family=binomial,data=evdat[evdat$posp7==1,]),type="response")
p9[is.na(evdat$posp9)==F & evdat$posp8==1]  <- predict(glm(posp9 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_8  + L2_8   + L3_8   + A_8   + Y_8   + L1_9  + L2_9  + L3_9, family=binomial,data=evdat[evdat$posp8==1,]),type="response")
p10[is.na(evdat$posp10)==F & evdat$posp9==1] <- predict(glm(posp10 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_9  + L2_9   + L3_9   + A_9   + Y_9   + L1_10 + L2_10 + L3_10, family=binomial,data=evdat[evdat$posp9==1,]),type="response")
p11[is.na(evdat$posp11)==F & evdat$posp10==1] <- predict(glm(posp11 ~ b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_10 + L2_10  + L3_10  + A_10  + Y_10  + L1_11 + L2_11 + L3_11, family=binomial,data=evdat[evdat$posp10==1,]),type="response")
p12[is.na(evdat$posp12)==F & evdat$posp11==1] <- predict(glm(posp12 ~  b1 + b2 + b3 + b4 + b5 + b6 + b7 + L1_11 + L2_11  + L3_11  + A_11  + Y_11  + L1_12 + L2_12 + L3_12, family=binomial,data=evdat[evdat$posp11==1,]),type="response")

pnever <- t(apply(cbind(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),1,cumprod))


plot(density(na.omit(p750)),col="red",xlim=c(0,0.05),ylim=c(0,0.5))
lines(density(na.omit(p350)))
lines(density(na.omit(palways)),col="blue")
lines(density(na.omit(pnever)),col="green")

# TABLE 3 #
mean(na.omit(p350[,12]) < 0.025)
mean(na.omit(p750[,12]) < 0.025)
mean(na.omit(palways[,12]) < 0.025)
mean(na.omit(pnever[,12]) < 0.025)













