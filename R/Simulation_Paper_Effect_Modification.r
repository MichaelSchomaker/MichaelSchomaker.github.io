#########################################################################################################################################
# Luque-Fernandez MA, Redondo-Sanchez D, Schomaker M.                                                                                   #
# Effect Modification and Collapsibility in Evaluations of Public Health Interventions.                                                 #
# American Journal of Public Health. 2019;109(3):e12-e3.                                                                                #
#                                                                                                                                       #
# File for simulation reported in the letter                                                                                            #
# A full summary of letter, answer, answer to answer and code is available at. https://github.com/migariane/HETMOR-Causal-Inference     #
#########################################################################################################################################

# Copyright (c) 2018  <Miguel Angel Luque-Fernandez, Daniel Redondo-Sanchez and Michael Schomaker>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

##############################################################################################

library(simcausal)      # currently only on CRAN Aarchive (as of 4 June 2020)
set.seed(24121980)

###############################################
# Data generating process for cancer example  #
###############################################

# Using library simcausal to simulate the data 
M <- DAG.empty()
M <- M + 
  node("w1", # age (0/1); 1 -> Advanced age
       distr = "rbern", prob = 0.6) + 
  node("w2", # ses (1/2/3/4/5); Advanced age, higher probability of belonging to higher socioeconomic status
       distr = "rcat.b1",
       probs = { 
         plogis(-3.10 - 0.05 * w1)    # 1 = upper middle class, 4%                         
         plogis(-1.25 - 0.04 * w1)    # 2 = middle class, 22%                              
         plogis(-0.05 - 0.03 * w1)    # 3 = lower middle  49%                              
         plogis(-1.65 - 0.02 * w1)    # 4 = working class  16%                             
         }) +                         # 5= last class comes from "normalizing", about 9%   
  node("w3", # comorbidities (1/2/3/4/5)
       distr = "rcat.b1",
       probs = {
         plogis(-0.5 + 0.010 * w1 + 0.15 * w2)      # 1 = none
         plogis(-1.8 + 0.005 * w1 + 0.10 * w2)      # 2 = Hypertension
         plogis(-1.6 + 0.010 * w1 + 0.12 * w2)      # 3 = Diabetes II
         plogis(-1.5 + 0.010 * w1 + 0.15 * w2)      # 4 = other
         plogis(-2.5 + 0.020 * w1 + 0.20 * w2)      # 5 = alzheimer
         plogis(-2.3 + 0.015 * w1 + 0.15 * w2)      # 6 = previous heart attack
         }) +  
  node("w4", # cancer stage (1/2/3/4), again last stage comes from "normalizing"
       distr = "rcat.b1",
       probs = {
         plogis(-0.4 + 0.02 * w1 - 0.050 * w2)    
         plogis(-1.0 + 0.03 * w1 - 0.055 * w2)
         plogis(-1.2 + 0.01 * w1 - 0.040 * w2)
         }) + 
  node("a", # treatment: a=1 -> dual therapy; a=0 -> mono therapy
       distr = "rbern",
       prob  = plogis(-2.3 + 0.05 * w2 + 0.25 * w3 + w4)) 

# 4 setups for the outcome under different levels of effect modification: no effect modifcation, mild EM, normal EM, and strong EM      
M1 <- M + node("y", distr = "rbern", prob = 1 - plogis(4.25 + 2.9 * a - 0.65 * w1 - 0.25 * w2 - 0.2 * w3 - 0.75 * w4))
M2 <- M + node("y", distr = "rbern", prob = 1 - plogis(4.25 + 2.9 * a - 0.65 * w1 - 0.25 * w2 - 0.2 * w3 - 0.75 * w4 - 0.10 * a * w2 - 0.20 * a * w3))
M3 <- M + node("y", distr = "rbern", prob = 1 - plogis(4.25 + 2.9 * a - 0.65 * w1 - 0.25 * w2 - 0.2 * w3 - 0.75 * w4 - 0.25 * a * w2 - 0.35 * a * w3))
M4 <- M + node("y", distr = "rbern", prob = 1 - plogis(4.25 + 2.9 * a - 0.65 * w1 - 0.25 * w2 - 0.2 * w3 - 0.75 * w4 - 0.35 * a * w2 - 0.45 * a * w3))

# Simulate  sample (n = 5000) in line with the data generating process and outcome generation under M2
Mset <- set.DAG(M2)
Odat <- simcausal::sim(DAG = Mset, n = 5000, rndseed = 345, verbose = F)
Odat1 <- Odat
Odat1$w1 <- factor(Odat$w1, levels = c(0:1), labels = c("young", "old"))
Odat1$w2 <- factor(Odat$w2, levels = c(1:5), labels = c("1_high", "2_high_middle", "3_middle", "4_low_middle", "5_low"))
Odat1$w3 <- factor(Odat$w3, levels = c(1:6), labels = c("none", "hypertension", "diabetes", "other", "alzheimer", "heart attack"))
Odat1$w4 <- factor(Odat$w4, levels = c(1:4), labels = c("stage 1", "stage 2", "stage 3", "stage 4"))
Odat$a   <- factor(Odat$a,  levels = c(0:1), labels = c("chemo", "chemo+radio"))
summary(Odat1)

# Estimate TRUE marginal causal odds ratio  for all 4 data setups: 
# M1 under no interaction between treatment and covariates = homogeneous treatment effect
# M2 mild interaction between treatment and covariates = mild heterogeneous treatment effect
# M3 interaction between treatment and covariates = heterogeneous treatment effect
# M4 strong interaction between treatment and covariates = stronger heterogeneous treatment effect

#############################
# Under M1 (no interaction) #
#############################
Mset <- set.DAG(M1)
a1 <- node("a", distr = "rbern", prob = 1)
Mset <- Mset + action("a1", nodes = a1)
a0 <- node("a", distr = "rbern", prob = 0)
Mset <- Mset + action("a0", nodes = a0)

dat <- simcausal::sim(DAG = Mset, actions = c("a1", "a0"), n = 10000000, rndseed = 666, verbose = FALSE)

Mset <- set.targetE(Mset, outcome = "y", param = "a1")
M1y1 <- eval.target(Mset, data = dat)$res              # P(death with chemo and radio) = 0.0352066 
Mset <- set.targetE(Mset, outcome = "y", param = "a0")  
M1y0 <- eval.target(Mset, data = dat)$res              # P(death with chemo) = 0.3253395
Mset <- set.targetE(Mset, outcome = "y", param = "a1-a0") 
True_ATE_M1 <- eval.target(Mset, data = dat)$res       # Evaluate true ATE

# True Marginal Odds Ratio of Treatment M1
# True_mor_M1 <- (0.0352066 /(1-0.0352066)) / (0.3253395 / (1 - 0.3253395)); True_mor_M1
True_mor_M1 <- M1y1 * (1 - M1y0) / ((1 - M1y1) * M1y0)
cat("MOR_M1:", True_mor_M1[[1]])

###############################
# Under M2 (mild interaction) #
###############################
Mset <- set.DAG(M2)
a1 <- node("a", distr = "rbern", prob = 1)
Mset <- Mset + action("a1", nodes = a1)
a0 <- node("a", distr = "rbern", prob = 0)
Mset <- Mset + action("a0", nodes = a0)

dat <- simcausal::sim(DAG = Mset, actions = c("a1", "a0"), n = 10000000, rndseed = 666, verbose = FALSE)

Mset <- set.targetE(Mset, outcome = "y", param = "a1")
M2y1 <- eval.target(Mset, data = dat)$res               # P(death with chemo and radio) = 0.0869532  
Mset <- set.targetE(Mset, outcome = "y", param = "a0")   
M2y0 <- eval.target(Mset, data = dat)$res               # P(death with chemo) = 0.3253395 
Mset <- set.targetE(Mset, outcome = "y", param = "a1-a0")
True_ATE_M2 <- eval.target(Mset, data = dat)$res        # Evaluate true ATE

# True Marginal Odds Ratio of Treatment M2
# True_mor_M2  <- (0.0869532  /(1-0.0869532))/(0.3253395/(1-0.3253395)); True_mor_M2
True_mor_M2 <- M2y1 * (1 - M2y0) / ((1 - M2y1) * M2y0)
cat("MOR_M2:", True_mor_M2[[1]])

##########################
# Under M3 (interaction) #
##########################
Mset <- set.DAG(M3)
a1 <- node("a", distr = "rbern", prob = 1)
Mset <- Mset + action("a1", nodes = a1)
a0 <- node("a", distr = "rbern", prob = 0)
Mset <- Mset + action("a0", nodes = a0)

dat <- simcausal::sim(DAG = Mset, actions = c("a1", "a0"), n = 10000000, rndseed = 666, verbose = FALSE)

Mset <- set.targetE(Mset, outcome = "y", param = "a1")
M3y1 <- eval.target(Mset, data = dat)$res               # P(death with chemo and radio) = 0.1830828  
Mset <- set.targetE(Mset, outcome = "y", param = "a0")    
M3y0 <- eval.target(Mset, data = dat)$res               # P(death with chemo) = 0.3253395
Mset <- set.targetE(Mset, outcome = "y", param = "a1-a0")
True_ATE_M3 <- eval.target(Mset, data = dat)$res        # Evaluate true ATE

# True Marginal Odds Ratio of Treatment M3
# True_mor_M3  <- (0.1830828 /(1-0.1830828))/(0.3253395/(1-0.3253395))
True_mor_M3 <- M3y1 * (1 - M3y0) / ((1 - M3y1) * M3y0)
cat("MOR_M3:", True_mor_M3[[1]])

##################################
# Under M4 (stronger interaction)#
##################################
Mset <- set.DAG(M4)
a1 <- node("a", distr = "rbern", prob = 1)
Mset <- Mset + action("a1", nodes = a1)
a0 <- node("a", distr = "rbern", prob = 0)
Mset <- Mset + action("a0", nodes = a0)

dat <- simcausal::sim(DAG = Mset, actions = c("a1", "a0"), n = 10000000, rndseed = 666, verbose = FALSE)

Mset <- set.targetE(Mset, outcome = "y", param = "a1")
M4y1 <- eval.target(Mset, data = dat)$res               # P(death with chemo and radio) = 0.2687067
Mset <- set.targetE(Mset, outcome = "y", param = "a0")  
M4y0 <- eval.target(Mset, data = dat)$res               # P(death with chemo) = 0.3253395
Mset <- set.targetE(Mset, outcome = "y", param = "a1-a0")
True_ATE_M4 <- eval.target(Mset, data = dat)$res        # Evaluate true ATE

# True Marginal Odds Ratio of Treatment M3
# True_mor_M4  <- (0.2687067   /(1-0.2687067))/(0.3253395/(1-0.3253395))
True_mor_M4 <- M4y1 * (1 - M4y0) / ((1 - M4y1) * M4y0)
cat("MOR_M4:", True_mor_M4[[1]])

# Plot for the Directed Acyclic Graph implied by the structural causal model
# pdf("plot-dag.pdf", width = 10, height = 6)
  plotDAG(Mset, xjitter = 0.3, yjitter = 0.04,
          edge_attrs = list(width = 0.5, arrow.width = 0.2, arrow.size = 0.3),
          vertex_attrs = list(size = 12, label.cex = 0.8))
# dev.off()

##########################
# Monte Carlo Simulation #
##########################

# Number of simulation runs
R <- 10000

results_reg <- matrix(NA, ncol = 4, nrow = R)
results_g   <- matrix(NA, ncol = 4, nrow = R)

for(r in 1:R) try({
    if (r%%10 == 0) cat(paste("This is simulation run number", r, "\n"))
    # M1
    ###########################
    Mset <- suppressWarnings(set.DAG(M1, verbose = FALSE))
    simdat <- suppressWarnings(simcausal::sim(DAG = Mset, n = 5000, verbose = FALSE))
    #
    results_reg[r, 1] <- exp(coefficients(glm(y ~ a + w1 + w2 + w3 + w4, data = simdat, family = binomial)))[2]
    #
    is1 <- simdat
    is1$a <- rep(1, dim(simdat)[1])
    is2 <- simdat
    is2$a <- rep(0, dim(simdat)[1])
    ms1 <- glm(y ~ a + w1 + w2 + w3 + w4, data = simdat, family = binomial)
    results_g[r, 1] <- ((mean(predict(ms1, type = "response", newdata = is1))) / (1 - mean(predict(ms1, type = "response", newdata = is1)))) / ((mean(predict(ms1, type = "response", newdata = is2))) / (1 - mean(predict(ms1, type = "response", newdata = is2))))
    # M2
    ##########################
    Mset <- suppressWarnings(set.DAG(M2, verbose = FALSE))
    simdat <- suppressWarnings(simcausal::sim(DAG = Mset, n = 5000, verbose = FALSE))
    #
    results_reg[r, 2] <- exp(coefficients(glm(y ~ a + w1 + w2 + w3 + w4, data = simdat, family = binomial)))[2]
    #
    is1 <- simdat
    is1$a <- rep(1, dim(simdat)[1])
    is2 <- simdat
    is2$a <- rep(0, dim(simdat)[1])
    ms1 <- glm(y ~ a + w1 + w2 + w3 + w4 + a * w2 + a * w3, data = simdat, family = binomial)
    results_g[r, 2] <- ((mean(predict(ms1, type = "response", newdata = is1))) / (1 - mean(predict(ms1, type = "response", newdata = is1))))/((mean(predict(ms1, type = "response", newdata = is2))) / (1 - mean(predict(ms1, type = "response", newdata = is2))))
    # M3
    ##########################
    Mset <- suppressWarnings(set.DAG(M3, verbose = FALSE))
    simdat <- suppressWarnings(simcausal::sim(DAG = Mset, n = 5000, verbose = FALSE))
    #
    results_reg[r, 3] <- exp(coefficients(glm(y ~ a + w1 + w2 + w3 + w4, data = simdat, family = binomial)))[2]
    #
    is1 <- simdat
    is1$a <- rep(1, dim(simdat)[1])
    is2 <- simdat
    is2$a <- rep(0, dim(simdat)[1])
    ms1 <- glm(y ~ a + w1 + w2 + w3 + w4 + a*w2 + a * w3, data = simdat, family = binomial)
    results_g[r, 3] <- ((mean(predict(ms1, type = "response", newdata = is1))) / (1 - mean(predict(ms1, type = "response", newdata = is1)))) / ((mean(predict(ms1, type = "response", newdata = is2))) / (1 - mean(predict(ms1, type = "response", newdata = is2))))
    # M4
    ##########################
    Mset <- suppressWarnings(set.DAG(M4, verbose = F))
    simdat <- suppressWarnings(simcausal::sim(DAG = Mset, n = 5000, verbose = FALSE))
    #
    results_reg[r, 4] <- exp(coefficients(glm(y ~ a + w1 + w2 + w3 + w4, data = simdat, family = binomial)))[2]
    #
    is1 <- simdat
    is1$a <- rep(1, dim(simdat)[1])
    is2 <- simdat
    is2$a <- rep(0, dim(simdat)[1])
    ms1 <- glm(y ~ a + w1 + w2 + w3 + w4 + a*w2 + a*w3, data = simdat, family = binomial)
    results_g[r, 4] <- ((mean(predict(ms1, type = "response", newdata = is1)))/(1 - mean(predict(ms1, type = "response", newdata = is1)))) / ((mean(predict(ms1, type = "response", newdata = is2))) / (1 - mean(predict(ms1, type = "response", newdata = is2))))
    ##########################
})

# Estimating BIAS for i) naive regression analysis  and ii) g-formula
BIAS_reg <- abs(apply(results_reg, 2, mean) - c(True_mor_M1, True_mor_M2, True_mor_M3, True_mor_M4))
BIAS_g   <- abs(apply(results_g,   2, mean) - c(True_mor_M1, True_mor_M2, True_mor_M3, True_mor_M4))

# Plot Bias.pdf
# change path if required
pdf("Figure2.pdf", width = 8) 
  plot(NA, xlab = "Simulation setting", ylab = "Absolute bias", cex.lab = 1.3, lwd = 2, axes = FALSE, ylim = c(-0.001, 0.12), xlim = c(0.8, 4.2))
  axis(side = 1, at = 1:4, labels = c("no EM*", "some EM", "EM", "strong EM"), cex.axis = 1.2)
  axis(side = 2, at = seq(0, 0.12, 0.02), las = 1, cex.axis = 1.2)   # ),labels = c("", "2%", "4%", "6%", "8%", "10%", "12%"))
  abline(h = 0, lwd = 2, lty = 2, col = "black")  
  lines(1:4, BIAS_reg, type = "p", cex = 2, lwd = 2, col = "red",  pch = 17)
  lines(1:4, BIAS_g,   type = "p", cex = 2, lwd = 2, col = "blue", pch = 16)
  legend("topleft", c("Multivariable Logistic regression", "G-Formula"), col = c("red", "blue"), pch = c(17, 16), cex = 1.3, bty = "n")
  mtext("*EM = Effect modification", side = 1, adj = 0.05, padj = -1.8, cex = 1.2, bty = "n", outer = TRUE)
dev.off()

# Plot Bias.tiff
tiff("Figure2.tiff", width = 520)
  plot(NA, xlab = "Simulation setting", ylab = "Absolute bias", cex.lab = 1.3, lwd = 2, axes = FALSE, ylim = c(-0.001, 0.12), xlim = c(0.8, 4.2))
  axis(side = 1, at = 1:4, labels = c("no EM*", "some EM", "EM", "strong EM"), cex.axis = 1.2)
  axis(side = 2, at = seq(0, 0.12, 0.02), las = 1, cex.axis = 1.2)   # ),labels = c("", "2%", "4%", "6%", "8%", "10%", "12%"))
  abline(h = 0, lwd = 2, lty = 2, col = "black")  
  lines(1:4, BIAS_reg, type = "p", cex = 2, lwd = 2, col = "red",  pch = 17)
  lines(1:4, BIAS_g,   type = "p", cex = 2, lwd = 2, col = "blue", pch = 16)
  legend("topleft", c("Multivariable Logistic regression", "G-Formula"), col = c("red", "blue"), pch = c(17, 16), cex = 1.3, bty = "n")
  mtext("*EM = Effect modification", side = 1, adj = 0.05, padj = -1.8, cex = 1.2, bty = "n", outer = TRUE)
dev.off()

# Cite as: Luque-Fernandez MA, Redondo-Sanchez D, Schomaker M. 
# Effect Modification and Collapsibility in Evaluations of Public Health Interventions. 
# American Journal of Public Health. 2019;109(3):e12-e3.