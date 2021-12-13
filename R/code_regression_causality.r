####################
# 0) Load Packages #
####################

library(simcausal)  # causal data simulation
library(doParallel) # for parallelization

#################################
# 1) Data-Generating Processes  #
#################################

# a) simple reference setting
M <- DAG.empty()
M <- M +
  node("L",
       distr = "rnorm",
       mean=1, sd=1) +
  node("A",
       distr = "rbern",
       prob = plogis(-0.5 + 2*L)) +
  node("Y",
       distr = "rnorm",
      mean= 2 + A + 3*L)
Mset <- set.DAG(M)   # setup 1
true1 <- 1 # true effect for A

# b) setup for incorrect model specification
M2 <- DAG.empty()
M2 <- M2 +
  node("L",
       distr = "rnorm",
       mean=1, sd=1) +
  node("A",
       distr = "rbern",
       prob = plogis(-0.5 + 2*L )) +
  node("Y",
       distr = "rnorm",
      mean= 2 + A + 0.5*L^2)
Mset2 <- set.DAG(M2) # setup 2
true2 <- 1 # true effect for A

# c) setup for collider bias
M3 <- DAG.empty()
M3 <- M3 +
  node("L",
       distr = "rnorm",
       mean=1, sd=1) +
  node("A",
       distr = "rbern",
       prob = plogis(-0.5 + 2*L)) +
  node("Y",
       distr = "rnorm",
      mean= 2 + A + 3*L) +
  node("L2",                 # L2 = collider
       distr = "rnorm",
       mean=Y*A, sd=1)
Mset3 <- set.DAG(M3) # setup 3
true3 <- 1 # true effect for A


# d) effect modification
M4 <- DAG.empty()
M4 <- M4 +
  node("L",
       distr = "rnorm",
       mean=1, sd=1) +
  node("A",
       distr = "rbern",
       prob = plogis(-0.5 + 2*L)) +
  node("Y",
       distr = "rnorm",
      mean= 2 + A + 3*L + A*L)
Mset4 <- set.DAG(M4) # setup 4
# true marginal ATE
a4.1 <- node("A", distr = "rbern", prob = 1)  # set A=1 and A=0 (intervene)
a4.0 <- node("A", distr = "rbern", prob = 0) 
Mset4 <- Mset4 + action("a4.0", nodes = a4.0) + action("a4.1", nodes = a4.1) # post-intervention DAG
int.dat4 <- simcausal::sim(DAG = Mset4, actions = c("a4.1", "a4.0"), n = 1000000, rndseed = 345)
Mset4 <- set.targetE(Mset4, outcome = "Y", param = "a4.1-a4.0")
true4 <- eval.target(Mset4, data = int.dat4)$res # true effect (ATE) for A

# as d); but with randomized treatment assignment
M4b <- DAG.empty()
M4b <- M4b +
  node("L",
       distr = "rnorm",
       mean=1, sd=1) +
  node("A",            # random treatment assignment
       distr = "rbern",
       prob = plogis(-0.5)) +
  node("Y",
       distr = "rnorm",
      mean= 2 + A + 3*L + A*L)
Mset4b <- set.DAG(M4b) # setup 4
# true marginal ATE
Mset4b <- Mset4b + action("a4.0", nodes = a4.0) + action("a4.1", nodes = a4.1) 
int.dat4b <- simcausal::sim(DAG = Mset4b, actions = c("a4.1", "a4.0"),
                            n = 1000000, rndseed = 345)
Mset4b <- set.targetE(Mset4b, outcome = "Y", param = "a4.1-a4.0")
true4b <- eval.target(Mset4b, data = int.dat4b)$res # true effect for A

# e) collapsibility
M5 <- DAG.empty()
M5 <- M5 +
  node("L",
       distr = "rnorm",
       mean=1, sd=1) +
  node("A",
       distr = "rbern",
       prob =  0.5) +   # randomized experiment
  node("Y",             # binary outcome
       distr = "rbern",
        prob = plogis(A+L))
Mset5 <- set.DAG(M5) # setup 5
# true marginal effects
a5.1 <- node("A", distr = "rbern", prob = 1); a5.0 <- node("A", distr = "rbern", prob = 0) # set A=1 and A=0 (intervene)
Mset5 <- Mset5 + action("a5.0", nodes = a5.0) + action("a5.1", nodes = a5.1)
int.dat5  <- simcausal::sim(DAG = Mset5, actions = c("a5.1", "a5.0"), n = 1000000, rndseed = 345)
Mset5     <- set.targetE(Mset5, outcome = "Y", param = "a5.1"); true5.1   <- eval.target(Mset5, data = int.dat5)$res # P(Y(had A=1)=1)
Mset5     <- set.targetE(Mset5, outcome = "Y", param = "a5.0"); true5.2   <- eval.target(Mset5, data = int.dat5)$res # P(Y(had A=0)=1)
true5_ATE <- true5.1-true5.2                                  # true ATE
true5_OR <- ((true5.1)/(1-true5.1))/((true5.2)/(1-true5.2))   # true MOR
true5_logOR <- log(((true5.1)/(1-true5.1))/((true5.2)/(1-true5.2))) # true log MOR

# f) mediation
M6 <- DAG.empty()
M6 <- M6 +
  node("A",
       distr = "rbern",
       prob = plogis(-0.5)) +
  node("M",                           # M = mediator
       distr = "rbern",
       prob = plogis(0.5 - 2*A)) +
  node("Y",
       distr = "rnorm",
      mean= 2 + M + A)
Mset6 <- set.DAG(M6) # setup 6
# true marginal effect estimates
a6.1 <- node("A", distr = "rbern", prob = 1)
a6.0 <- node("A", distr = "rbern", prob = 0)
Mset6 <- Mset6 + action("a6.0", nodes = a6.0) + action("a6.1", nodes = a6.1)
int.dat6 <- simcausal::sim(DAG = Mset6, actions = c("a6.1", "a6.0"), n = 1000000, rndseed = 345)
Mset6 <- set.targetE(Mset6, outcome = "Y", param = "a6.1-a6.0")
true6 <- eval.target(Mset6, data = int.dat6)$res  # true ATE

# g)
M7 <- DAG.empty()
M7 <- M7 +
  node("L1",
       distr = "rnorm",
       mean=1, sd=1) +
  node("L2",
       distr = "rnorm",
       mean=-1, sd=1) +
  node("A",
       distr = "rbern",
       prob = plogis(-0.5 + 2*L1 + L2)) +
  node("CA",             # Prob. that A is missing
       distr = "rbern",
       prob = plogis(1.5+0.5*L1 + 0.5*L2)) +
  node("CL2",            # Prob. that L2 is missing
       distr = "rbern",
       prob = plogis(1.5-0.75*L1 + 0.75*A)) +
  node("Y",
       distr = "rnorm",
      mean= 2 + A)  +
  node("CY",             # Prob. that Y is missing
       distr = "rbern",
       prob = plogis(1.5+0.25*L1 + 0.25*L2 + 0.5*A))
Mset7 <- set.DAG(M7) # setup 7
true7 <- 1           # true ATE

####################
# 2) simulation    #
####################

runs <- 10000  # number of simulation runs
N    <- 1000   # sample size

# Simulation loop
cl <- makeCluster(max(detectCores()-1,1)) # use all cores of computer, minus 1
registerDoParallel(cl)                    # start parallelization

sim <- foreach(r = 1:runs, .combine=rbind, .packages="simcausal") %dopar% {

# draw data
simdat   <- sim(DAG = Mset,   n = N, verbose=F)
simdat2  <- sim(DAG = Mset2,  n = N, verbose=F)
simdat3  <- sim(DAG = Mset3,  n = N, verbose=F)
simdat4  <- sim(DAG = Mset4,  n = N, verbose=F)
simdat4b <- sim(DAG = Mset4b, n = N, verbose=F)
simdat5  <- sim(DAG = Mset5,  n = N, verbose=F)
simdat6  <- sim(DAG = Mset6,  n = N, verbose=F)
simdat7  <- sim(DAG = Mset7,  n = N, verbose=F)
# models for the different setups
m1    <- lm(Y~A+L,data=simdat)
m2    <- lm(Y~A+L,data=simdat2)
m3    <- lm(Y~A+L+L2,data=simdat3)
m4    <- lm(Y~A+L,data=simdat4)
m4b   <- lm(Y~A+L,data=simdat4b)
m5    <- glm(Y~A+L,data=simdat5,family=binomial)
m5.2  <- glm(Y~A,data=simdat5,family=binomial)
m6    <- lm(Y~A+M,data=simdat6)
m6.2  <- lm(Y~A,data=simdat6)
m7    <- lm(Y~A+L1+L2,
                  data=simdat7[simdat7$CY==1 & simdat7$CA==1 & simdat7$CL2==1,])
# store treatment estimates
results <- c(coef(m1)[2], coef(m2)[2], coef(m3)[2], coef(m4)[2], coef(m4b)[2],
             coef(m5)[2], coef(m5.2)[2], coef(m6)[2], coef(m6.2)[2],coef(m7)[2])
}

stopCluster(cl)   # stop parallelization

#########################
# 3) Evaluating Results #
#########################

# Bias
truth <- c(true1,true2,true3,true4,true4b,true5_logOR,true5_logOR,true6,true6,true7)
BIAS <- apply(sim,2,mean)-truth
names(BIAS) <- c("ATE, setup 1:\n simple \n","ATE, setup 2:\n incorrect MS \n",
                 "ATE, setup 3:\n collider \n structure \n",
                 "ATE, setup 4:\n effect \n modification \n (no random.) \n",
                 "ATE, setup 4b:\n effect \n modification \n (random.)\n",
                 "MOR, setup 5:\n collapsibility, \n condititional \n",
                 "MOR, setup 5:\n collapsibility, \n crude \n",
                 "ATE, setup 6:\n mediation,\n conditional \n",
                 "ATE, setup 6:\n mediation,\n crude \n",
                 "ATE, setup 7:\n MNAR+CC \n \n")
BIAS

# Visualize Bias
library(ggplot2)
Bias_data <- as.data.frame(matrix(NA, nrow=length(BIAS),ncol=2))
Bias_data[,1] <- BIAS
Bias_data[,2] <- names(BIAS)
colnames(Bias_data) <- c("Bias","Setup")

pdf(file=paste(getwd(),"/BIAS.pdf",sep=""), width=12)
ggplot(as.data.frame(Bias_data), aes(x=Setup,y=Bias,size=I(2.5))) + geom_point() +
       theme_bw() +
       scale_x_discrete("Simulation Setup") +
       scale_y_continuous("Bias", breaks=seq(-1,1,0.25))  +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size=12),axis.title.y = element_text(size=12, angle = 90), axis.text.y = element_text(size=12))+
  geom_vline(aes(xintercept = 8.5), size = .25, linetype = "dashed") +
  geom_hline(aes(yintercept = 0)) +
  annotate("text", x = 1.5, y = 0.75, label = "Bias with respect to ATE") +
  annotate("text", x = 9.5, y = 0.75, label = "Bias with respect to log MOR")+
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = 0.7, ymax = 0.8,  alpha = .2)  +
  annotate("rect", xmin = 8.6, xmax = 10.4, ymin = 0.7, ymax = 0.8,  alpha = .2)
dev.off()