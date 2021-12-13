###################################################################################################################################################
# Supplementary Material to Schomaker M, Hogger S, Johnson LF, Hoffmann CJ, Barnighausen T, Heumann C.                                            #
# Simultaneous Treatment of Missing Data and Measurement Error in HIV Research Using Multiple Overimputation.                                     #
# Epidemiology. 2015;26(5):628-36.                                                                                                                #
#                                                                                                                                                 #
# Simulations from the paper                                                                                                                      #
# Designed by Sara Hogger, updated by Michael Schomaker                                                                                           # 
#                                                                                                                                                 #
###################################################################################################################################################


# Set working directory
#setwd('C:/Users/schomakm/Dropbox/Documents/Research_Statistics/Code_Reproduce/Overimputation')
setwd("//home//mschomaker//Code_Reproduce//Overimputation")

library(Amelia)       # amelia
library(survival)     # coxph
library(mvtnorm)      # rmvnorm
library(copula)

################################
# 1) Function: Simulate Data   #
#    (only Cox model,lognormal)#
################################

sim_data <- function(model, n, S, covariates, var_me=F, dev_me=F, conf_me=F, prob_mv, lp, wrong){

  # List to return simulated data
  data <- list()

  # Draw lognormal variables via copula
  if (covariates=="lognormal"){
  mycopula <-  mvdc(claytonCopula(1, dim=2),c("lnorm","lnorm"),list(list(meanlog=4.286, sdlog=1.086),list(meanlog=10.760, sdlog=1.808607)))
  if(lp=="+var"){mycopula <-  mvdc(claytonCopula(1, dim=6),c("lnorm","lnorm","binom","weibull","exp","gamma"),list(list(meanlog=4.286, sdlog=1.086),list(meanlog=10.760, sdlog=1.808607),list(size=1, prob=0.65),list(shape=1.75,scale=1.9),list(rate=1),list(shape=0.25,scale=2)))}
  mycovdata <- rMvdc(n,mycopula)
  if(lp!="+var"){  colnames(mycovdata)<- c("x1","x2")}
  if(lp=="+var"){  colnames(mycovdata)<- c("x1","x2","x3","x4","x5","x6")}
  mycovdata <- as.data.frame(mycovdata)
  x_true <- mycovdata
  if(lp!="+var"){x_true <- as.data.frame(cbind(log(x_true[,1]),log10(x_true[,2])))}
  if(lp=="+var"){x_true <- as.data.frame(cbind(log(x_true[,1]),log10(x_true[,2]),x_true[,3:6]))}
  if(lp!="+var"){colnames(x_true)<- c("x1","x2")}
  if(lp=="+var"){colnames(x_true)<- c("x1","x2","x3","x4","x5","x6")}
  y_true <- 0.5 + 0.75*x_true[,1] + 1.25*x_true[,2]
  } 
  
  # Draw design-matrix for normal covariates 
  if (covariates=="normal"){
    x_true <- rmvnorm(n, mean = rep(0, 2), sigma = matrix(c(1,0.2,0.2,1), nrow=2)) 
    colnames(x_true)<- c("x1","x2")         
  }


  # Draw measurement error, missing values and error-term epsilon
  for (s in 1:S){

    # Observed covariates
    x <- x_true

    # homoscedastic measurement error
    if(dev_me[1]==F & conf_me[1]==F){
      x[,1] <- x[,1] + rnorm(n, sd=sqrt(var_me[1]))
      x[,2] <- x[,2] + rnorm(n, sd=sqrt(var_me[2]))
    } else{stop("not homoscedastic ME")}

    if (model=="lm") {
      trueb <- c(0.75,1.25)
      y_true <- 0.5 + 0.75*x_true[,1] + 1.25*x_true[,2]
      y <- y_true + 0.505*rnorm(n)      
      data[[s]] <- data.frame(x=x, y=y)
    }
    
    if (model=="logit") {
      trueb <- c(0.2,-0.7)
      pi <- exp(0.5 + 0.2*x_true[,1] - 0.7*x_true[,2]) / (1 + exp(0.5 + 0.2*x_true[,1] - 0.7*x_true[,2]))
      y <- rbinom(n, 1, prob=pi)
      data[[s]] <- data.frame(x=x, y=y)
    }

     if (model=="cox") {
      trueb <- c(-0.3,0.3)
      if(lp=="small"){trueb <- c(-0.1,0.1)} 
      # Survival-time T, Censoring Time C
      D <- -log(runif(n)) / (0.1 * exp(trueb[1]*x_true[,1] + trueb[2]*x_true[,2]))
      C <- -log(runif(n)) / 0.2
      y <- apply(cbind(D,C), 1, min)
    } 

    # missing values
    if(prob_mv[1]=="MARs"){
    ind_mv2 <- as.logical(rbinom(n, 1, prob= (1-(1/(1+exp(1-4*(y)))))))
    x[ind_mv2,2] <- NA
    x[ind_mv2,1] <- NA
    }
    
    if(prob_mv[1]=="MARl"){
    ind_mv1 <- as.logical(rbinom(n, 1, prob=(1-(1/(0.02*(y)^2+1))) ))
    x[ind_mv1,2] <- NA
    x[ind_mv1,1] <- NA
    }
    
    if(prob_mv[1]=="MCARs"){
    ind_mv3 <- as.logical(rbinom(n, 1, prob=0.1))
    x[ind_mv3,2] <- NA
    x[ind_mv3,1] <- NA
    }
    
    if(prob_mv[1]=="MCARl"){
    ind_mv4 <- as.logical(rbinom(n, 1, prob=0.4))
    x[ind_mv4,2] <- NA
    x[ind_mv4,1] <- NA
    }
    
    
    
    # Calculate response Variable

     if (model=="cox") {
      c <- 1 * (D < C)
      data[[s]] <- data.frame(x=x, y=y, c=c)
      if(lp!="+var"){colnames(data[[s]])[1:2]<- c("x.1","x.2")}
      if(lp=="+var"){colnames(data[[s]])[1:6]<- c("x.1","x.2","x.3","x.4","x.5","x.6")}
    } 
     if(model!="cox"){
     data[[s]] <- data.frame(x=x, y=y)
     colnames(data[[s]])[1:2]<- c("x.1","x.2")
     }  


  }

  return(list("data"=data, "x_true"=x_true))
}


################################
# 2) Function: Imputation      #
################################

impute_data <- function(data, var_me=F, dev_me=F, conf_me=F, amelia_m, lp, model, wrong, ...){


# Lists to return (over-)imputed datasets
data_imp <- list()
data_overimp <- list()


# Loop s through all simulated datasets
for (s in 1:length(data)){

  print(paste("(Over)impute data nr.", s ))

  ### Number of observations
  n <- nrow(data[[s]])

  ### Indexes of columns with covariates
  index_x <- c(1,2)#grep("x", names(data[[s]]))

  ### Check if there are missing values in the dataset
  if(sum(is.na(data[[s]]))!=0) {

  ### Just impute missing values, ignore measurement error
  if(model=="cox"){
  if(lp!="+var"){
  imp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  incheck=F,
                  overimp=NULL,
                  logs="y"
                  )}
                  
  
  if(lp=="+var"){
  imp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  incheck=F,
                  overimp=NULL,
                  noms="x.3",
                  logs=c("y","x.5","x.6")
                  )}
  }
  
  if(model=="lm"){
   imp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  incheck=F,
                  overimp=NULL,
                  )}
                  
   if(model=="logit"){
   imp <- amelia(x=data[[s]],
                  m=amelia_m,
                  noms="y",
                  p2s=0,
                  incheck=F,
                  overimp=NULL,
                  )}
  
  

  data_imp[[s]] <- do.call(rbind, imp$imputations)
  data_imp[[s]]$id <- as.factor(rep(1:n, times=imp$m))
  data_imp[[s]]$imp <- rep(1:imp$m, each=n)
  


  ### Calculate prior and overimp
  priors <- list()
  overimp <- list()

  for(j in index_x) {

    # Prior for homoscedastic measurement error
    if (dev_me[1]==F & conf_me[1]==F){
        priors[[j]] <- na.omit(cbind("row" = 1:n,         # Michael: changed to define priors only for mismeasured, not missing data
                          "column"= rep(j, times=n),
                          "mean" = data[[s]][,j],
                          "sd" = sqrt(var_me[j])))
    }  else{stop("not homoscedastic ME")}

    # Cells to overimpute
    overimp[[j]] <- cbind("row" = which(!is.na(data[[s]][,j])),
                         "column"= j)
  }

  priors <- do.call(rbind, priors)
  overimp <- do.call(rbind, overimp)
  }


  ### Calculate priors and overimp if there are no missing values
  if(sum(is.na(data[[s]]))==0) {

  priors <- list()
  overimp <- list()

  for(j in index_x) {

      # Prior for homoscedastic measurement error
      if (dev_me[1]==F & conf_me[1]==F){
        priors[[j]] <- na.omit(cbind("row" = 1:n,
                             "column"= rep(j, times=n),
                             "mean" = data[[s]][,j],
                             "sd" = sqrt(var_me[j])+wrong))
      }  else{stop("not homoscedastic ME")}

      # Cells to overimpute
      overimp[[j]] <- cbind("row" = which(!is.na(data[[s]][,j])),
                            "column"= j)
    }

    priors <- do.call(rbind, priors)
    overimp <- do.call(rbind, overimp)
  }


  ### Run amelia
  if(model=="cox"){
  if(lp!="+var"){
  overimp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  priors=priors,
                  overimp=overimp,
                  incheck=F,
                  logs="y"
                  )}
                  
  
  if(lp=="+var"){
  overimp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  incheck=F,
                  priors=priors,
                  overimp=overimp,
                  noms=c("x.3"),
                  logs=c("y","x.5","x.6")
                  )}
  }
  
  if(model=="lm"){
  overimp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  priors=priors,
                  overimp=overimp,
                  incheck=F
                  )}
  
  
  if(model=="logit"){
  overimp <- amelia(x=data[[s]],
                  m=amelia_m,
                  p2s=0,
                  priors=priors,
                  overimp=overimp,
                  noms="y",
                  incheck=F
                  )}
  

  ### Save (over-) imputed datasets
  data_overimp[[s]] <- do.call(rbind, overimp$imputations)
  data_overimp[[s]]$id <- as.factor(rep(1:n, times=overimp$m))
  data_overimp[[s]]$imp <- rep(1:overimp$m, each=n)


} # End loop s

return(list("data_imp"=data_imp, "data_overimp"=data_overimp))

} # End function impute_data


### Function to aggregate the standard deviation
aggr_std <- function(x){
  sqrt(mean(x^2, na.rm=T) + (1+1/length(x)) * var(x, na.rm=T))
}


################################
# 3) Function: Check density   #
################################

check_density <- function(density_true, density_org, density_imp, density_overimp, main="", max=0){

  plot(density_true, lwd = 2, cex.lab=1.5, cex.axis=1.5,
       main=main, ylab="Relative Dichte",
       xlab="x.1",
       ylim=c(0, max(density_true$y, density_org$y, density_overimp[[1]]$y, max)),
       col="forestgreen"
  )

  lines(density_org, lwd = 2)

  try(
    for (i in 1:length(density_imp)){
      lines(density_imp[[i]], col="royalblue1")
    },
    silent=T
  )

  for (i in 1:length(density_overimp)){
    lines(density_overimp[[i]], col="orangered1")
  }

} # End function check_density


################################
# 4) Function: Simulation      #
################################

simulate <- function(model, n=n, S=S, covariates, var_me, dev_me, conf_me, prob_mv, lp, amelia_m=5, wrong) {

  ###-------------------------------------------------
  ### Create data
  if (model=="cox"){
    set.seed(230312)
    data_sim <- sim_data(model="cox", n=n, S=S, covariates=covariates, var_me=var_me, dev_me=dev_me, conf_me=conf_me, prob_mv=prob_mv, lp=lp)
    data_org <- data_sim$data
    formula <- as.formula(Surv(y,c)~x.1+x.2)
    if(lp=="+var"){formula <- as.formula(Surv(y,c)~x.1+x.2+x.3+x.4+x.5+x.6)}
    beta_true <- c(-0.3,0.3)
    if(lp=="small"){beta_true <- c(-0.1,0.1)} 
    if(lp=="+var"){beta_true <- c(-0.3,0.3,0,0,0,0)}
  }  

  if (model=="lm"){
    set.seed(230312)
    data_sim <- sim_data(model="lm", n=n, S=S, covariates="normal", var_me=var_me, dev_me=dev_me, conf_me=conf_me, prob_mv=prob_mv, lp=lp)
    data_org <- data_sim$data
    formula <- as.formula(y~x.1+x.2)
    beta_true <- c(0.5, 0.75,1.25)
  } 
  if (model=="logit"){
    set.seed(230312)
    data_sim <- sim_data (model="logit", n=n, S=S, covariates="normal", var_me=var_me, dev_me=dev_me, conf_me=conf_me, prob_mv=prob_mv)
    data_org <- data_sim$data
    formula <- as.formula(y~x.1+x.2)
    beta_true <- c(0.5,0.2,-0.7)
  } 

  ###-------------------------------------------------
  ### (Over-) Impute data and check densities (only for x.1, first simulated dataset)
  data <- impute_data(data=data_org, var_me=var_me, dev_me=dev_me, conf_me=conf_me, amelia_m=amelia_m, lp=lp, wrong=wrong, model=model)

  # density with true values
  density_true <- density(data_sim$x_true[,1])

  # Density with observed values
  density_org <- density(data_org[[1]][["x.1"]], na.rm=T, bw=density_true$bw)

  # Density with imputed values
  if (prob_mv[1]!=0){
    density_imp <- list()
    for (i in 1:amelia_m){
      density_imp[[i]] <- density(data$data_imp[[1]][["x.1"]][data$data_imp[[1]]$imp==i], bw=density_true$bw)
    }
  }   else density_imp <- NULL

  # Density with overimputed values
  density_overimp <- list()
  for (i in 1:amelia_m){
    density_overimp[[i]] <- density(data$data_overimp[[1]][["x.1"]][data$data_overimp[[1]]$imp==i], bw=density_true$bw)
  }


  ###-------------------------------------------------
  ### Calculate models

  ### Matrices to save results
  if (model=="cox"){
    coefs <- coxph(formula=formula, data=data_org[[1]])
  } else {
    coefs <- lm(formula=formula, data=data_org[[1]])
  }

  m <- length(coefs$coefficients)

  estimates_org <- data.frame(matrix(data=NA, nrow=S, ncol=2*m))
  names(estimates_org) <- c(names(coefs$coefficients), paste(names(coefs$coefficients), "_std", sep=""))
  predicted_org <- data.frame(matrix(data=NA, nrow=S, ncol=n))

  estimates_imp <- data.frame(matrix(data=NA, nrow=S*amelia_m, ncol=2*m+2))
  names(estimates_imp) <- c(names(coefs$coefficients), paste(names(coefs$coefficients), "_std", sep=""), "data", "imp")

  estimates_overimp <- estimates_imp


  ### Loop through the datasets
  for (k in 1:S) {

    ### Calculate models for the original datasets
    if (model=="lm") {
      model_org <- lm(formula=formula, data=data_org[[k]], x=T)
      estimates_org[k, 1:(2*m)] <- c(summary(model_org)$coefficients[,1], summary(model_org)$coefficients[,2])
      predicted_org[k, which(!is.na(rowSums(data_org[[k]])))] <- predict(model_org)
    }
    if (model=="logit") {
      model_org <- glm(formula=formula, family=binomial, data=data_org[[k]], maxit=100, x=T)
      estimates_org[k, 1:(2*m)] <- c(summary(model_org)$coefficients[,1], summary(model_org)$coefficients[,2])
      predicted_org[k, which(!is.na(rowSums(data_org[[k]])))] <- predict(model_org, type="response")
    }
    if (model=="cox") {
      model_org <- coxph(formula=formula, data=data_org[[k]])
      estimates_org[k, 1:(2*m)] <- c(summary(model_org)$coefficients[,1], summary(model_org)$coefficients[,3])
    }


    ### Loop through the imputations
    
    for (l in 1:amelia_m) {

      ### In which row should the results be saved
      index <- (k-1)*amelia_m + l

      ### Calculate models for the (over-)imputed datasets
      if (model=="lm") {          
        try(model_imp <- lm(formula=formula, data=data$data_imp[[k]][which(data$data_imp[[k]]$imp==l),], x=T), silent=T)
        try(estimates_imp[index, 1:(2*m)] <- c(summary(model_imp)$coefficients[,1], summary(model_imp)$coefficients[,2]), silent=T)
        
        model_overimp <- lm(formula=formula, data=data$data_overimp[[k]][which(data$data_overimp[[k]]$imp==l),], x=T)
        estimates_overimp[index, 1:(2*m)] <- c(summary(model_overimp)$coefficients[,1], summary(model_overimp)$coefficients[,2])
      }
      
      if (model=="logit") { 
        try(model_imp <- glm(formula=formula, family=binomial, data=data$data_imp[[k]][which(data$data_imp[[k]]$imp==l),], maxit=100, x=T), silent=T)
        try(estimates_imp[index, 1:(2*m)] <- c(summary(model_imp)$coefficients[,1], summary(model_imp)$coefficients[,2]), silent=T)
        
        model_overimp <- glm(formula=formula, family=binomial, data=data$data_overimp[[k]][which(data$data_overimp[[k]]$imp==l),], maxit=100, x=T)
        estimates_overimp[index, 1:(2*m)] <- c(summary(model_overimp)$coefficients[,1], summary(model_overimp)$coefficients[,2])
      }      
      if (model=="cox") {
        try(model_imp <- coxph(formula=formula, data=data$data_imp[[k]][which(data$data_imp[[k]]$imp==l),]), silent=T)
        try(estimates_imp[index, 1:(2*m)] <- c(summary(model_imp)$coefficients[,1], summary(model_imp)$coefficients[,3]), silent=T)

        model_overimp <- coxph(formula=formula, data=data$data_overimp[[k]][which(data$data_overimp[[k]]$imp==l),])
        estimates_overimp[index, 1:(2*m)] <- c(summary(model_overimp)$coefficients[,1], summary(model_overimp)$coefficients[,3])
      }

      estimates_imp[index, which(colnames(estimates_imp)=="data")] <- k
      estimates_imp[index, which(colnames(estimates_imp)=="imp")] <- l

      estimates_overimp[index, which(colnames(estimates_overimp)=="data")] <- k
      estimates_overimp[index, which(colnames(estimates_overimp)=="imp")] <- l

    } # End loop through all imputations

  } # End loop through the datasets


  # Aggregate estimates and predicted values of the imputations
  try(estimates_imp_mean <- aggregate(estimates_imp[,1:m],
                                    by=list(estimates_imp[,which(colnames(estimates_imp)=="data")]),
                                    mean), silent=T)
  try(estimates_imp_std <- aggregate(estimates_imp[,(m+1):(2*m)],
                                   by=list(estimates_imp[,which(colnames(estimates_imp)=="data")]),
                                   aggr_std), silent=T)
  try(estimates_imp_aggr  <- merge(estimates_imp_mean, estimates_imp_std)[,-1], silent=T)

  estimates_overimp_mean <- aggregate(estimates_overimp[,1:m],
                                  by=list(estimates_overimp[,which(colnames(estimates_overimp)=="data")]),
                                  mean)
  estimates_overimp_std <- aggregate(estimates_overimp[,(m+1):(2*m)],
                                 by=list(estimates_overimp[,which(colnames(estimates_overimp)=="data")]),
                                 aggr_std)
  estimates_overimp_aggr <- merge(estimates_overimp_mean, estimates_overimp_std)[,-1]


  ### Calculate MSE_beta and quality of prediction MSE_y
  MSE_beta_org <- apply(estimates_org[,1:m], 1, function(x, beta_true){1/m*(sum((x - beta_true)^2))}, beta_true=beta_true)
  MSE_beta_imp <- apply(estimates_imp_aggr[,1:m], 1, function(x, beta_true){1/m*(sum((x - beta_true)^2))}, beta_true=beta_true)
  MSE_beta_overimp <- apply(estimates_overimp_aggr[,1:m], 1, function(x, beta_true){1/m*(sum((x - beta_true)^2))}, beta_true=beta_true)


  ### Calculate differences between MSE_beta_org-MSE_beta_imp, MSE_beta_org-MSE_beta_imp bzw. MSE_y
  MSE_beta_imp_diff <- MSE_beta_org - MSE_beta_imp
  MSE_beta_overimp_diff <- MSE_beta_org - MSE_beta_overimp

  ### Output MSE and density
  return(list("MSE_beta_org" = MSE_beta_org,
              "MSE_beta_imp" = MSE_beta_imp,
              "MSE_beta_overimp" = MSE_beta_overimp,
              "MSE_beta_imp_diff"=MSE_beta_imp_diff,
              "MSE_beta_overimp_diff"=MSE_beta_overimp_diff,
              "density_true"=density_true,
              "density_org"=density_org,
              "density_imp"=density_imp,
              "density_overimp"=density_overimp,
              "estimates_org"=estimates_org,
              "estimates_imp"=estimates_imp_aggr,
              "estimates_overimp"=estimates_overimp_aggr
              ))

} # End function simulate

##########################
# 5) Perform Simulations #
##########################

# number of simulations
sims = 1000

# Main simulation
cox_nm  <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F, wrong=0)
cox_m   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F, wrong=0)

# Variation of sample size - needed for Figure 1
cox_nm_s1  <- simulate(model="cox", n=500, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F,wrong=0)
cox_m_s1   <- simulate(model="cox", n=500, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F,wrong=0)

cox_nm_s2  <- simulate(model="cox", n=1000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F,wrong=0)
cox_m_s2   <- simulate(model="cox", n=1000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F,wrong=0)

cox_nm_s3  <- simulate(model="cox", n=2500, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F,wrong=0)
cox_m_s3   <- simulate(model="cox", n=2500, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F,wrong=0)

cox_nm_s4  <- simulate(model="cox", n=7500, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F,wrong=0)
cox_m_s4   <- simulate(model="cox", n=7500, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F,wrong=0)

cox_nm_s5  <- simulate(model="cox", n=10000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F,wrong=0)
cox_m_s5   <- simulate(model="cox", n=10000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F,wrong=0)

# Several variations, for supplementary material

# Variation: Measurement Error
cox_nm_me1  <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.04, 0.0225), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F, wrong=0)
cox_m_me1   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.04, 0.0225), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F, wrong=0)

cox_nm_me2  <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.09, 0.0961), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F, wrong=0)
cox_m_me2   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.09, 0.0961), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F, wrong=0)

# Variation: Missing data
cox_m_m1   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARl","MARl"), dev_me=F, conf_me=F, wrong=0)
cox_m_m2   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MCARs","MCARs"), dev_me=F, conf_me=F, wrong=0)
cox_m_m3   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MCARl","MCARl"), dev_me=F, conf_me=F, wrong=0)

# Variation: linear predictor
cox_nm_lp1  <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="small", prob_mv=c(0, 0)         , dev_me=F, conf_me=F, wrong=0)
cox_m_lp1   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="small", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F, wrong=0)

cox_nm_lp2  <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="+var", prob_mv=c(0, 0)         , dev_me=F, conf_me=F, wrong=0)
cox_m_lp2   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="+var", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F, wrong=0)

# Variation: assumed ME variance
cox_nm_w1  <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c(0, 0)         , dev_me=F, conf_me=F, wrong=0.1)
cox_m_w1   <- simulate(model="cox", n=5000, S=sims, covariates="lognormal", var_me=c(0.0676, 0.065025), lp="default", prob_mv=c("MARs","MARs"), dev_me=F, conf_me=F, wrong=0.1)


#####################################################
# 6) Functions for plotting results of simulations  #
#####################################################

plot_estimates <- function(x_org, x_imp, x_overimp, xlab, main="", min, max, model,para=list(-0.3, 0.3)) {

  if (model=="cox") {

    names <- c(expression(X[1]),expression(X[2]))

    boxplot(para, at=c(1:length(x_org)), boxcol="forestgreen", medcol="forestgreen", ylim=c(min, max), xaxt="n", yaxt="n")

    if(is.na(x_imp[[1]][1])){
      boxplot(x_org, boxwex=0.15, col="gray70", at=c(1:ncol(x_org))-0.1,
              names=F, ylab="Estimates", xlab=xlab, main=main, cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
      suppressWarnings(boxplot(x_imp, boxwex=0.15, col="royalblue1", at=c(1:length(x_org)), names=names, cex.lab=1.5, cex.axis=1.5, add=T))
      boxplot(x_overimp, boxwex=0.15, col="orangered1", at=c(1:length(x_org)+0.1), names=F, xaxt="n", yaxt="n", add=T)
    }
    if(!is.na(x_imp[[1]][1])){
      boxplot(x_org, boxwex=0.15, col="gray70", at=c(1:length(x_org))-0.2,
              names=F, ylab="Estimates", xlab=xlab, main=main, cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
      boxplot(x_imp, boxwex=0.15, col="royalblue1", at=c(1:length(x_org)), names=names, cex.lab=1.5, cex.axis=1.5, add=T)
      boxplot(x_overimp, boxwex=0.15, col="orangered1", at=c(1:length(x_org)+0.2), names=F, xaxt="n", yaxt="n", add=T)
    }
  }
  
  if (model=="lm") {
    
    names <- c(expression(X[0]),expression(X[1]),expression(X[2]))
    
    boxplot(para, at=c(1:length(x_org)), boxcol="forestgreen", medcol="forestgreen", ylim=c(min, max), xaxt="n", yaxt="n")
    
    if(is.na(x_imp[[1]][1])){      
      boxplot(x_org, boxwex=0.15, col="gray70", at=c(1:ncol(x_org))-0.1, 
              names=F, ylab="Estimates", xlab=xlab, main=main, cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
      suppressWarnings(boxplot(x_imp, boxwex=0.15, col="royalblue1", at=c(1:length(x_org)), names=names, cex.lab=1.5, cex.axis=1.5, add=T))
      boxplot(x_overimp, boxwex=0.15, col="orangered1", at=c(1:length(x_org)+0.1), names=F, xaxt="n", yaxt="n", add=T) 
    }
    if(!is.na(x_imp[[1]][1])){
      boxplot(x_org, boxwex=0.15, col="gray70", at=c(1:length(x_org))-0.2,  
              names=F, ylab="Estimates", xlab=xlab, main=main, cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
      boxplot(x_imp, boxwex=0.15, col="royalblue1", at=c(1:length(x_org)), names=names, cex.lab=1.5, cex.axis=1.5, add=T)
      boxplot(x_overimp, boxwex=0.15, col="orangered1", at=c(1:length(x_org)+0.2), names=F, xaxt="n", yaxt="n", add=T)
    }
  }
  
  if (model=="logit") {
    
    names <- c(expression(X[0]),expression(X[1]),expression(X[2]))
    
    boxplot(para, at=c(1:length(x_org)), boxcol="forestgreen", medcol="forestgreen", ylim=c(min, max), xaxt="n", yaxt="n")
    
    if(is.na(x_imp[[1]][1])){      
      boxplot(x_org, boxwex=0.15, col="gray70", at=c(1:ncol(x_org))-0.1, 
              names=F, ylab="Estimates", xlab=xlab, main=main, cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
      suppressWarnings(boxplot(x_imp, boxwex=0.15, col="royalblue1", at=c(1:length(x_org)), names=names, cex.lab=1.5, cex.axis=1.5, add=T))
      boxplot(x_overimp, boxwex=0.15, col="orangered1", at=c(1:length(x_org)+0.1), names=F, xaxt="n", yaxt="n", add=T) 
    }
    if(!is.na(x_imp[[1]][1])){
      boxplot(x_org, boxwex=0.15, col="gray70", at=c(1:length(x_org))-0.2,  
              names=F, ylab="Estimates", xlab=xlab, main=main, cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
      boxplot(x_imp, boxwex=0.15, col="royalblue1", at=c(1:length(x_org)), names=names, cex.lab=1.5, cex.axis=1.5, add=T)
      boxplot(x_overimp, boxwex=0.15, col="orangered1", at=c(1:length(x_org)+0.2), names=F, xaxt="n", yaxt="n", add=T)
    }
  }

} # End function plot_estimates



########################################################
# 7) PRODUCE FINAL RESULTS (not all reported in paper) #
########################################################

path_results <- getwd()

###################
# MAIN SIMULATION #
###################

### a) point estimates ###

pdf(file=paste(path_results, "/main_sim_point_estimates.pdf", sep=""))
boxplot(list(-0.3, 0.3), at=c(1:length(cox_nm$estimates_org[,1:2])), ylim=c(-0.5,0.75), boxcol="forestgreen", medcol="forestgreen", xaxt="n", yaxt="n")
boxplot(cox_nm$estimates_org[,1:2], boxwex=0.15, col="gray70", at=c(1:ncol(cox_nm$estimates_org[,1:2]))-0.1,names = F, ylab="Estimates", cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
boxplot(list(-10.3, 10.3), boxwex=0.15, col="royalblue1", at=c(1:ncol(cox_nm$estimates_org[,1:2])), names=c(expression(X[1]),expression(X[2])), cex.lab=1.5, cex.axis=1.5, add=T)
boxplot(cox_nm$estimates_overimp[,1:2], boxwex=0.15, col="orangered1", at=c(1:length(cox_nm$estimates_org[,1:2])+0.1), names=F, xaxt="n", yaxt="n", add=T)
legend("topleft", legend=c("true parameter","complete cases","multiple overimputation"), col=c("forestgreen","gray70","orangered1"), lwd=2, cex=1.2)
dev.off()

pdf(file=paste(path_results, "/main_sim_point_estimates_2.pdf", sep=""))
plot_estimates(x_org=cox_m$estimates_org[,1:2],
               x_imp=cox_m$estimates_imp[,1:2],
               x_overimp=cox_m$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.75)
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

### b) MSE ###
mse_0_x1_org         <- c(mean((cox_nm$estimates_org[,1]+0.3)^2),(mean(cox_nm$estimates_org[,1]+0.3)))
mse_0_x2_org         <- c(mean((cox_nm$estimates_org[,2]-0.3)^2),(mean(cox_nm$estimates_org[,2]-0.3)))

mse_0_x1_overimp     <- c(mean((cox_nm$estimates_overimp[,1]+0.3)^2),(mean(cox_nm$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp     <- c(mean((cox_nm$estimates_overimp[,2]-0.3)^2),(mean(cox_nm$estimates_overimp[,2]-0.3)))

mse_1_x1_org         <- c(mean((cox_m$estimates_org[,1]+0.3)^2),(mean(cox_m$estimates_org[,1]+0.3)))
mse_1_x2_org         <- c(mean((cox_m$estimates_org[,2]-0.3)^2),(mean(cox_m$estimates_org[,2]-0.3)))

mse_1_x1_imp         <- c(mean((cox_m$estimates_imp[,1]+0.3)^2),(mean(cox_m$estimates_imp[,1]+0.3)))
mse_1_x2_imp         <- c(mean((cox_m$estimates_imp[,2]-0.3)^2),(mean(cox_m$estimates_imp[,2]-0.3)))

mse_1_x1_overimp     <- c(mean((cox_m$estimates_overimp[,1]+0.3)^2),(mean(cox_m$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp     <- c(mean((cox_m$estimates_overimp[,2]-0.3)^2),(mean(cox_m$estimates_overimp[,2]-0.3)))

# Summarize MSE and Bias (beta1, X1)
pdf(file=paste(path_results, "/main_sim_MSE_Bias_X1.pdf", sep=""))
plot(0 , xlab="" , ylab = "" , main = "" , sub="MSE and squared Bias for different methods (X1)", type = "n" , axes=F, xlim=c(0.5,2.5), ylim=c(0,0.011))
lines(0.8,mse_0_x1_org[1],type="p",col="black",pch=16,cex=2.5)
lines(0.8,mse_0_x1_org[2],type="p",col="black",pch=18,cex=2.5)
#lines(1,  mse_0_x1_imp[1],type="p",col="forestgreen",pch=16,cex=2.5)
#lines(1,  mse_0_x1_imp[2],type="p",col="forestgreen",pch=18,cex=2.5)
lines(1.2,mse_0_x1_overimp[1],type="p",col="red",pch=16,cex=2.5)
lines(1.2,mse_0_x1_overimp[2],type="p",col="red",pch=18,cex=2.5)
lines(1.8,mse_1_x1_org[1],type="p",col="black",pch=16,cex=2.5)
lines(1.8,mse_1_x1_org[2],type="p",col="black",pch=18,cex=2.5)
lines(2,  mse_1_x1_imp[1],type="p",col="forestgreen",pch=16,cex=2.5)
lines(2,  mse_1_x1_imp[2],type="p",col="forestgreen",pch=18,cex=2.5)
lines(2.2,mse_1_x1_overimp[1],type="p",col="red",pch=16,cex=2.5)
lines(2.2,mse_1_x1_overimp[2],type="p",col="red",pch=18,cex=2.5)
axis(side = 1, at = c(0.5,1,2,2.5), tick=T, las=1, labels = c(" ","no missing data","missing data"," "))
axis(side=2, at = c(0,0.005,0.010), labels=c("0","0.005","0.010"))
legend("topleft",bty="n",col=c("black","forestgreen","red"),legend=c("Complete Cases","Multiple Imputation","Multiple Overimputation"), lwd=2, cex=1.25)
legend("left",bty="n",col=c("grey","grey"),legend=c("MSE","sq. Bias"), pch=c(16,18), cex=1.5)
dev.off()

# Summarize MSE and Bias (beta2, X2)
pdf(file=paste(path_results, "/main_sim_point_MSE_Bias_X2.pdf", sep=""))
plot(0 , xlab="" , ylab = "" , main = "" , sub="MSE and squared Bias for different methods (X2)", type = "n" , axes=F, xlim=c(0.5,2.5), ylim=c(0,0.025))
lines(0.8,mse_0_x2_org[1],type="p",col="black",pch=16,cex=2.5)
lines(0.8,mse_0_x2_org[2],type="p",col="black",pch=18,cex=2.5)
#lines(1,  mse_0_x2_imp[1],type="p",col="forestgreen",pch=16,cex=2.5)
#lines(1,  mse_0_x2_imp[2],type="p",col="forestgreen",pch=18,cex=2.5)
lines(1.2,mse_0_x2_overimp[1],type="p",col="red",pch=16,cex=2.5)
lines(1.2,mse_0_x2_overimp[2],type="p",col="red",pch=18,cex=2.5)
lines(1.8,mse_1_x2_org[1],type="p",col="black",pch=16,cex=2.5)
lines(1.8,mse_1_x2_org[2],type="p",col="black",pch=18,cex=2.5)
lines(2,  mse_1_x2_imp[1],type="p",col="forestgreen",pch=16,cex=2.5)
lines(2,  mse_1_x2_imp[2],type="p",col="forestgreen",pch=18,cex=2.5)
lines(2.2,mse_1_x2_overimp[1],type="p",col="red",pch=16,cex=2.5)
lines(2.2,mse_1_x2_overimp[2],type="p",col="red",pch=18,cex=2.5)
axis(side = 1, at = c(0.5,1,2,2.5), tick=T, las=1, labels = c(" ","no missing data","missing data"," "))
axis(side=2)
legend("topleft",bty="n",col=c("black","forestgreen","red"),legend=c("Complete Cases","Multiple Imputation","Multiple Overimputation"), lwd=2, cex=1.25)
legend("left",bty="n",col=c("grey","grey"),legend=c("MSE","sq. Bias"), pch=c(16,18), cex=1.5)
dev.off()

### c) Standard error ####

se_nm_org_1 <- c(sd(cox_nm$estimates_org[,1]),mean(cox_nm$estimates_org[,3]))
se_nm_org_2 <- c(sd(cox_nm$estimates_org[,2]),mean(cox_nm$estimates_org[,4]))

se_nm_overimp_1 <- c(sd(cox_nm$estimates_overimp[,1]),mean(cox_nm$estimates_overimp[,3]))
se_nm_overimp_2 <- c(sd(cox_nm$estimates_overimp[,2]),mean(cox_nm$estimates_overimp[,4]))


se_m_org_1 <- c(sd(cox_m$estimates_org[,1]),mean(cox_m$estimates_org[,3]))
se_m_org_2 <- c(sd(cox_m$estimates_org[,2]),mean(cox_m$estimates_org[,4]))

se_m_imp_1 <- c(sd(cox_m$estimates_imp[,1]),mean(cox_m$estimates_imp[,3]))
se_m_imp_2 <- c(sd(cox_m$estimates_imp[,2]),mean(cox_m$estimates_imp[,4]))

se_m_overimp_1 <- c(sd(cox_m$estimates_overimp[,1]),mean(cox_m$estimates_overimp[,3]))
se_m_overimp_2 <- c(sd(cox_m$estimates_overimp[,2]),mean(cox_m$estimates_overimp[,4]))

summaryse <- cbind(rbind(se_nm_org_1,se_nm_overimp_1,se_m_org_1,se_m_imp_1,se_m_overimp_1),rbind(se_nm_org_2,se_nm_overimp_2,se_m_org_2,se_m_imp_2,se_m_overimp_2))
colnames(summaryse) <- c("Sim SE, X1","Est SE, X1","Sim SE, X2","Est SE, X2")
rownames(summaryse) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryse,file=paste(path_results, "/main_sim_summmary_SE.csv", sep=""))

### d) MSE & Bias Table ###

summaryMSE <- cbind(rbind(mse_0_x1_org,mse_0_x1_overimp,mse_1_x1_org,mse_1_x1_imp,mse_1_x1_overimp),
                    rbind(mse_0_x2_org,mse_0_x2_overimp,mse_1_x2_org,mse_1_x2_imp,mse_1_x2_overimp))
colnames(summaryMSE) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryMSE,file=paste(path_results, "/main_sim_summmary_MSE_Bias.csv", sep=""))


################################
# Variation: Measurement Error #
################################


pdf(file=paste(path_results, "/supp_fig_f.pdf", sep=""))
boxplot(list(-0.3, 0.3), at=c(1:length(cox_nm_me1$estimates_org[,1:2])), ylim=c(-0.5,0.5), boxcol="forestgreen", medcol="forestgreen", xaxt="n", yaxt="n")
boxplot(cox_nm_me1$estimates_org[,1:2], boxwex=0.15, col="gray70", at=c(1:ncol(cox_nm_me1$estimates_org[,1:2]))-0.1,names = F, ylab="Estimates", cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
boxplot(list(-10.3, 10.3), boxwex=0.15, col="royalblue1", at=c(1:ncol(cox_nm_me1$estimates_org[,1:2])), names=c(expression(X[1]),expression(X[2])), cex.lab=1.5, cex.axis=1.5, add=T)
boxplot(cox_nm_me1$estimates_overimp[,1:2], boxwex=0.15, col="orangered1", at=c(1:length(cox_nm_me1$estimates_org[,1:2])+0.1), names=F, xaxt="n", yaxt="n", add=T)
legend("topleft", legend=c("true parameter","complete cases","multiple overimputation"), col=c("forestgreen","gray70","orangered1"), lwd=2, cex=1.2)
dev.off()

pdf(file=paste(path_results, "/supp_fig_g.pdf", sep=""))
plot_estimates(x_org=cox_m_me1$estimates_org[,1:2],
               x_imp=cox_m_me1$estimates_imp[,1:2],
               x_overimp=cox_m_me1$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5)
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

mse_0_x1_org_me1         <- c(mean((cox_nm_me1$estimates_org[,1]+0.3)^2),(mean(cox_nm_me1$estimates_org[,1]+0.3)))
mse_0_x2_org_me1         <- c(mean((cox_nm_me1$estimates_org[,2]-0.3)^2),(mean(cox_nm_me1$estimates_org[,2]-0.3)))

mse_0_x1_overimp_me1     <- c(mean((cox_nm_me1$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_me1$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_me1     <- c(mean((cox_nm_me1$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_me1$estimates_overimp[,2]-0.3)))

mse_1_x1_org_me1         <- c(mean((cox_m_me1$estimates_org[,1]+0.3)^2),(mean(cox_m_me1$estimates_org[,1]+0.3)))
mse_1_x2_org_me1         <- c(mean((cox_m_me1$estimates_org[,2]-0.3)^2),(mean(cox_m_me1$estimates_org[,2]-0.3)))

mse_1_x1_imp_me1         <- c(mean((cox_m_me1$estimates_imp[,1]+0.3)^2),(mean(cox_m_me1$estimates_imp[,1]+0.3)))
mse_1_x2_imp_me1         <- c(mean((cox_m_me1$estimates_imp[,2]-0.3)^2),(mean(cox_m_me1$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_me1     <- c(mean((cox_m_me1$estimates_overimp[,1]+0.3)^2),(mean(cox_m_me1$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_me1     <- c(mean((cox_m_me1$estimates_overimp[,2]-0.3)^2),(mean(cox_m_me1$estimates_overimp[,2]-0.3)))


summaryMSE_me1 <- cbind(rbind(mse_0_x1_org_me1,mse_0_x1_overimp_me1,mse_1_x1_org_me1,mse_1_x1_imp_me1,mse_1_x1_overimp_me1),
                    rbind(mse_0_x2_org_me1,mse_0_x2_overimp_me1,mse_1_x2_org_me1,mse_1_x2_imp_me1,mse_1_x2_overimp_me1))
colnames(summaryMSE_me1) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_me1) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryMSE_me1,file=paste(path_results, "/supp_table_f_g.csv", sep=""))


pdf(file=paste(path_results, "/supp_fig_d.pdf", sep=""))
boxplot(list(-0.3, 0.3), at=c(1:length(cox_nm_me2$estimates_org[,1:2])), ylim=c(-0.5,0.5), boxcol="forestgreen", medcol="forestgreen", xaxt="n", yaxt="n")
boxplot(cox_nm_me2$estimates_org[,1:2], boxwex=0.15, col="gray70", at=c(1:ncol(cox_nm_me2$estimates_org[,1:2]))-0.1,names = F, ylab="Estimates", cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
boxplot(list(-10.3, 10.3), boxwex=0.15, col="royalblue1", at=c(1:ncol(cox_nm_me2$estimates_org[,1:2])), names=c(expression(X[1]),expression(X[2])), cex.lab=1.5, cex.axis=1.5, add=T)
boxplot(cox_nm_me2$estimates_overimp[,1:2], boxwex=0.15, col="orangered1", at=c(1:length(cox_nm_me2$estimates_org[,1:2])+0.1), names=F, xaxt="n", yaxt="n", add=T)
legend("topleft", legend=c("true parameter","complete cases","multiple overimputation"), col=c("forestgreen","gray70","orangered1"), lwd=2, cex=1.2)
dev.off()


pdf(file=paste(path_results, "/supp_fig_e.pdf", sep=""))
plot_estimates(x_org=cox_m_me2$estimates_org[,1:2],
               x_imp=cox_m_me2$estimates_imp[,1:2],
               x_overimp=cox_m_me2$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5)
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

mse_0_x1_org_me2         <- c(mean((cox_nm_me2$estimates_org[,1]+0.3)^2),(mean(cox_nm_me2$estimates_org[,1]+0.3)))
mse_0_x2_org_me2         <- c(mean((cox_nm_me2$estimates_org[,2]-0.3)^2),(mean(cox_nm_me2$estimates_org[,2]-0.3)))

mse_0_x1_overimp_me2     <- c(mean((cox_nm_me2$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_me2$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_me2     <- c(mean((cox_nm_me2$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_me2$estimates_overimp[,2]-0.3)))

mse_1_x1_org_me2         <- c(mean((cox_m_me2$estimates_org[,1]+0.3)^2),(mean(cox_m_me2$estimates_org[,1]+0.3)))
mse_1_x2_org_me2         <- c(mean((cox_m_me2$estimates_org[,2]-0.3)^2),(mean(cox_m_me2$estimates_org[,2]-0.3)))

mse_1_x1_imp_me2         <- c(mean((cox_m_me2$estimates_imp[,1]+0.3)^2),(mean(cox_m_me2$estimates_imp[,1]+0.3)))
mse_1_x2_imp_me2         <- c(mean((cox_m_me2$estimates_imp[,2]-0.3)^2),(mean(cox_m_me2$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_me2     <- c(mean((cox_m_me2$estimates_overimp[,1]+0.3)^2),(mean(cox_m_me2$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_me2     <- c(mean((cox_m_me2$estimates_overimp[,2]-0.3)^2),(mean(cox_m_me2$estimates_overimp[,2]-0.3)))



summaryMSE_me2 <- cbind(rbind(mse_0_x1_org_me2,mse_0_x1_overimp_me2,mse_1_x1_org_me2,mse_1_x1_imp_me2,mse_1_x1_overimp_me2),
                    rbind(mse_0_x2_org_me2,mse_0_x2_overimp_me2,mse_1_x2_org_me2,mse_1_x2_imp_me2,mse_1_x2_overimp_me2))
colnames(summaryMSE_me2) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_me2) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryMSE_me2,file=paste(path_results, "/supp_table_d_e.csv", sep=""))


################################
# Variation: Missing data      #
################################

pdf(file=paste(path_results, "/supp_fig_c.pdf", sep=""))
plot_estimates(x_org=cox_m_m1$estimates_org[,1:2],
               x_imp=cox_m_m1$estimates_imp[,1:2],
               x_overimp=cox_m_m1$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5)
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()


mse_1_x1_org_m1         <- c(mean((cox_m_m1$estimates_org[,1]+0.3)^2),(mean(cox_m_m1$estimates_org[,1]+0.3)))
mse_1_x2_org_m1         <- c(mean((cox_m_m1$estimates_org[,2]-0.3)^2),(mean(cox_m_m1$estimates_org[,2]-0.3)))

mse_1_x1_imp_m1         <- c(mean((cox_m_m1$estimates_imp[,1]+0.3)^2),(mean(cox_m_m1$estimates_imp[,1]+0.3)))
mse_1_x2_imp_m1         <- c(mean((cox_m_m1$estimates_imp[,2]-0.3)^2),(mean(cox_m_m1$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_m1     <- c(mean((cox_m_m1$estimates_overimp[,1]+0.3)^2),(mean(cox_m_m1$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_m1     <- c(mean((cox_m_m1$estimates_overimp[,2]-0.3)^2),(mean(cox_m_m1$estimates_overimp[,2]-0.3)))



summaryMSE_m1 <- cbind(rbind(mse_1_x1_org_m1,mse_1_x1_imp_m1,mse_1_x1_overimp_m1),
                    rbind(mse_1_x2_org_m1,mse_1_x2_imp_m1,mse_1_x2_overimp_m1))
colnames(summaryMSE_m1) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_m1) <- c("m: CC","m: MI","m: MO")
write.csv(summaryMSE_m1,file=paste(path_results, "/supp_table_c.csv", sep=""))




pdf(file=paste(path_results, "/supp_fig_a.pdf", sep=""))
plot_estimates(x_org=cox_m_m2$estimates_org[,1:2],
               x_imp=cox_m_m2$estimates_imp[,1:2],
               x_overimp=cox_m_m2$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5)
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

mse_1_x1_org_m2         <- c(mean((cox_m_m2$estimates_org[,1]+0.3)^2),(mean(cox_m_m2$estimates_org[,1]+0.3)))
mse_1_x2_org_m2         <- c(mean((cox_m_m2$estimates_org[,2]-0.3)^2),(mean(cox_m_m2$estimates_org[,2]-0.3)))

mse_1_x1_imp_m2         <- c(mean((cox_m_m2$estimates_imp[,1]+0.3)^2),(mean(cox_m_m2$estimates_imp[,1]+0.3)))
mse_1_x2_imp_m2         <- c(mean((cox_m_m2$estimates_imp[,2]-0.3)^2),(mean(cox_m_m2$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_m2     <- c(mean((cox_m_m2$estimates_overimp[,1]+0.3)^2),(mean(cox_m_m2$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_m2     <- c(mean((cox_m_m2$estimates_overimp[,2]-0.3)^2),(mean(cox_m_m2$estimates_overimp[,2]-0.3)))



summaryMSE_m2 <- cbind(rbind(mse_1_x1_org_m2,mse_1_x1_imp_m2,mse_1_x1_overimp_m2),
                    rbind(mse_1_x2_org_m2,mse_1_x2_imp_m2,mse_1_x2_overimp_m2))
colnames(summaryMSE_m2) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_m2) <- c("m: CC","m: MI","m: MO")
write.csv(summaryMSE_m2,file=paste(path_results, "/supp_table_a.csv", sep=""))



pdf(file=paste(path_results, "/supp_fig_b.pdf", sep=""))
plot_estimates(x_org=cox_m_m3$estimates_org[,1:2],
               x_imp=cox_m_m3$estimates_imp[,1:2],
               x_overimp=cox_m_m3$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5)
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

mse_1_x1_org_m3         <- c(mean((cox_m_m3$estimates_org[,1]+0.3)^2),(mean(cox_m_m3$estimates_org[,1]+0.3)))
mse_1_x2_org_m3         <- c(mean((cox_m_m3$estimates_org[,2]-0.3)^2),(mean(cox_m_m3$estimates_org[,2]-0.3)))

mse_1_x1_imp_m3         <- c(mean((cox_m_m3$estimates_imp[,1]+0.3)^2),(mean(cox_m_m3$estimates_imp[,1]+0.3)))
mse_1_x2_imp_m3         <- c(mean((cox_m_m3$estimates_imp[,2]-0.3)^2),(mean(cox_m_m3$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_m3     <- c(mean((cox_m_m3$estimates_overimp[,1]+0.3)^2),(mean(cox_m_m3$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_m3     <- c(mean((cox_m_m3$estimates_overimp[,2]-0.3)^2),(mean(cox_m_m3$estimates_overimp[,2]-0.3)))



summaryMSE_m3 <- cbind(rbind(mse_1_x1_org_m3,mse_1_x1_imp_m3,mse_1_x1_overimp_m3),
                    rbind(mse_1_x2_org_m3,mse_1_x2_imp_m3,mse_1_x2_overimp_m3))
colnames(summaryMSE_m3) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_m3) <- c("m: CC","m: MI","m: MO")
write.csv(summaryMSE_m3,file=paste(path_results, "/supp_table_b.csv", sep=""))

################################
# Variation: Linear Predictor  #
################################

pdf(file=paste(path_results, "/supp_fig_h.pdf", sep=""))
boxplot(list(-0.1, 0.1), at=c(1:length(cox_nm_lp1$estimates_org[,1:2])), ylim=c(-0.5,0.5), boxcol="forestgreen", medcol="forestgreen", xaxt="n", yaxt="n")
boxplot(cox_nm_lp1$estimates_org[,1:2], boxwex=0.15, col="gray70", at=c(1:ncol(cox_nm_lp1$estimates_org[,1:2]))-0.1,names = F, ylab="Estimates", cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
boxplot(list(-10.3, 10.3), boxwex=0.15, col="royalblue1", at=c(1:ncol(cox_nm_lp1$estimates_org[,1:2])), names=c(expression(X[1]),expression(X[2])), cex.lab=1.5, cex.axis=1.5, add=T)
boxplot(cox_nm_lp1$estimates_overimp[,1:2], boxwex=0.15, col="orangered1", at=c(1:length(cox_nm_lp1$estimates_org[,1:2])+0.1), names=F, xaxt="n", yaxt="n", add=T)
legend("topleft", legend=c("true parameter","complete cases","multiple overimputation"), col=c("forestgreen","gray70","orangered1"), lwd=2, cex=1.2)
dev.off()

pdf(file=paste(path_results, "/supp_fig_i.pdf", sep=""))
plot_estimates(x_org=cox_m_lp1$estimates_org[,1:2],
               x_imp=cox_m_lp1$estimates_imp[,1:2],
               x_overimp=cox_m_lp1$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5,para=list(-0.1, 0.1))
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

mse_0_x1_org_lp1         <- c(mean((cox_nm_lp1$estimates_org[,1]+0.1)^2),(mean(cox_nm_lp1$estimates_org[,1]+0.1)))
mse_0_x2_org_lp1         <- c(mean((cox_nm_lp1$estimates_org[,2]-0.1)^2),(mean(cox_nm_lp1$estimates_org[,2]-0.1)))

mse_0_x1_overimp_lp1     <- c(mean((cox_nm_lp1$estimates_overimp[,1]+0.1)^2),(mean(cox_nm_lp1$estimates_overimp[,1]+0.1)))
mse_0_x2_overimp_lp1     <- c(mean((cox_nm_lp1$estimates_overimp[,2]-0.1)^2),(mean(cox_nm_lp1$estimates_overimp[,2]-0.1)))

mse_1_x1_org_lp1         <- c(mean((cox_m_lp1$estimates_org[,1]+0.1)^2),(mean(cox_m_lp1$estimates_org[,1]+0.1)))
mse_1_x2_org_lp1         <- c(mean((cox_m_lp1$estimates_org[,2]-0.1)^2),(mean(cox_m_lp1$estimates_org[,2]-0.1)))

mse_1_x1_imp_lp1         <- c(mean((cox_m_lp1$estimates_imp[,1]+0.1)^2),(mean(cox_m_lp1$estimates_imp[,1]+0.1)))
mse_1_x2_imp_lp1         <- c(mean((cox_m_lp1$estimates_imp[,2]-0.1)^2),(mean(cox_m_lp1$estimates_imp[,2]-0.1)))

mse_1_x1_overimp_lp1     <- c(mean((cox_m_lp1$estimates_overimp[,1]+0.1)^2),(mean(cox_m_lp1$estimates_overimp[,1]+0.1)))
mse_1_x2_overimp_lp1     <- c(mean((cox_m_lp1$estimates_overimp[,2]-0.1)^2),(mean(cox_m_lp1$estimates_overimp[,2]-0.1)))



summaryMSE_lp1 <- cbind(rbind(mse_0_x1_org_lp1,mse_0_x1_overimp_lp1,mse_1_x1_org_lp1,mse_1_x1_imp_lp1,mse_1_x1_overimp_lp1),
                    rbind(mse_0_x2_org_lp1,mse_0_x2_overimp_lp1,mse_1_x2_org_lp1,mse_1_x2_imp_lp1,mse_1_x2_overimp_lp1))
colnames(summaryMSE_lp1) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_lp1) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryMSE_lp1,file=paste(path_results, "/supp_table_h_i.csv", sep=""))



pdf(file=paste(path_results, "/supp_fig_j.pdf", sep=""))
boxplot(list(-0.3, 0.3), at=c(1:length(cox_nm_lp2$estimates_org[,1:2])), ylim=c(-0.5,0.5), boxcol="forestgreen", medcol="forestgreen", xaxt="n", yaxt="n")
boxplot(cox_nm_lp2$estimates_org[,1:2], boxwex=0.15, col="gray70", at=c(1:ncol(cox_nm_lp2$estimates_org[,1:2]))-0.1,names = F, ylab="Estimates", cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
boxplot(list(-10.3, 10.3), boxwex=0.15, col="royalblue1", at=c(1:ncol(cox_nm_lp2$estimates_org[,1:2])), names=c(expression(X[1]),expression(X[2])), cex.lab=1.5, cex.axis=1.5, add=T)
boxplot(cox_nm_lp2$estimates_overimp[,1:2], boxwex=0.15, col="orangered1", at=c(1:length(cox_nm_lp2$estimates_org[,1:2])+0.1), names=F, xaxt="n", yaxt="n", add=T)
legend("topleft", legend=c("true parameter","complete cases","multiple overimputation"), col=c("forestgreen","gray70","orangered1"), lwd=2, cex=1.2)
dev.off()


pdf(file=paste(path_results, "/supp_fig_k.pdf", sep=""))
plot_estimates(x_org=cox_m_lp2$estimates_org[,1:2],
               x_imp=cox_m_lp2$estimates_imp[,1:2],
               x_overimp=cox_m_lp2$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5,para=list(-0.3, 0.3))
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()


mse_0_x1_org_lp2         <- c(mean((cox_nm_lp2$estimates_org[,1]+0.3)^2),(mean(cox_nm_lp2$estimates_org[,1]+0.3)))
mse_0_x2_org_lp2         <- c(mean((cox_nm_lp2$estimates_org[,2]-0.3)^2),(mean(cox_nm_lp2$estimates_org[,2]-0.3)))

mse_0_x1_overimp_lp2     <- c(mean((cox_nm_lp2$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_lp2$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_lp2     <- c(mean((cox_nm_lp2$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_lp2$estimates_overimp[,2]-0.3)))

mse_1_x1_org_lp2         <- c(mean((cox_m_lp2$estimates_org[,1]+0.3)^2),(mean(cox_m_lp2$estimates_org[,1]+0.3)))
mse_1_x2_org_lp2         <- c(mean((cox_m_lp2$estimates_org[,2]-0.3)^2),(mean(cox_m_lp2$estimates_org[,2]-0.3)))

mse_1_x1_imp_lp2         <- c(mean((cox_m_lp2$estimates_imp[,1]+0.3)^2),(mean(cox_m_lp2$estimates_imp[,1]+0.3)))
mse_1_x2_imp_lp2         <- c(mean((cox_m_lp2$estimates_imp[,2]-0.3)^2),(mean(cox_m_lp2$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_lp2     <- c(mean((cox_m_lp2$estimates_overimp[,1]+0.3)^2),(mean(cox_m_lp2$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_lp2     <- c(mean((cox_m_lp2$estimates_overimp[,2]-0.3)^2),(mean(cox_m_lp2$estimates_overimp[,2]-0.3)))


summaryMSE_lp2 <- cbind(rbind(mse_0_x1_org_lp2,mse_0_x1_overimp_lp2,mse_1_x1_org_lp2,mse_1_x1_imp_lp2,mse_1_x1_overimp_lp2),
                    rbind(mse_0_x2_org_lp2,mse_0_x2_overimp_lp2,mse_1_x2_org_lp2,mse_1_x2_imp_lp2,mse_1_x2_overimp_lp2))
colnames(summaryMSE_lp2) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_lp2) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryMSE_lp2,file=paste(path_results, "/supp_table_j_k.csv", sep=""))

#####################################################
# Variation: Assumption Measurement Error Variance  #
#####################################################


pdf(file=paste(path_results, "/supp_fig_l.pdf", sep=""))
boxplot(list(-0.3, 0.3), at=c(1:length(cox_nm_w1$estimates_org[,1:2])), ylim=c(-0.5,0.5), boxcol="forestgreen", medcol="forestgreen", xaxt="n", yaxt="n")
boxplot(cox_nm_w1$estimates_org[,1:2], boxwex=0.15, col="gray70", at=c(1:ncol(cox_nm_w1$estimates_org[,1:2]))-0.1,names = F, ylab="Estimates", cex.lab=1.5, cex.axis=1.5, xaxt="n", add=T)
boxplot(list(-10.3, 10.3), boxwex=0.15, col="royalblue1", at=c(1:ncol(cox_nm_w1$estimates_org[,1:2])), names=c(expression(X[1]),expression(X[2])), cex.lab=1.5, cex.axis=1.5, add=T)
boxplot(cox_nm_w1$estimates_overimp[,1:2], boxwex=0.15, col="orangered1", at=c(1:length(cox_nm_w1$estimates_org[,1:2])+0.1), names=F, xaxt="n", yaxt="n", add=T)
legend("topleft", legend=c("true parameter","complete cases","multiple overimputation"), col=c("forestgreen","gray70","orangered1"), lwd=2, cex=1.2)
dev.off()


pdf(file=paste(path_results, "/supp_fig_m.pdf", sep=""))
plot_estimates(x_org=cox_m_w1$estimates_org[,1:2],
               x_imp=cox_m_w1$estimates_imp[,1:2],
               x_overimp=cox_m_w1$estimates_overimp[,1:2],
               xlab="", model="cox", min=-0.5, max=0.5,para=list(-0.3, 0.3))
legend("topleft", legend=c("true parameter","complete cases","multiple imputation","multiple overimputation"), col=c("forestgreen","gray70","royalblue1","orangered1"), lwd=2, cex=1.2)
dev.off()

mse_0_x1_org_w1         <- c(mean((cox_nm_w1$estimates_org[,1]+0.3)^2),(mean(cox_nm_w1$estimates_org[,1]+0.3)))
mse_0_x2_org_w1         <- c(mean((cox_nm_w1$estimates_org[,2]-0.3)^2),(mean(cox_nm_w1$estimates_org[,2]-0.3)))

mse_0_x1_overimp_w1     <- c(mean((cox_nm_w1$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_w1$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_w1     <- c(mean((cox_nm_w1$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_w1$estimates_overimp[,2]-0.3)))

mse_1_x1_org_w1         <- c(mean((cox_m_w1$estimates_org[,1]+0.3)^2),(mean(cox_m_w1$estimates_org[,1]+0.3)))
mse_1_x2_org_w1         <- c(mean((cox_m_w1$estimates_org[,2]-0.3)^2),(mean(cox_m_w1$estimates_org[,2]-0.3)))

mse_1_x1_imp_w1         <- c(mean((cox_m_w1$estimates_imp[,1]+0.3)^2),(mean(cox_m_w1$estimates_imp[,1]+0.3)))
mse_1_x2_imp_w1         <- c(mean((cox_m_w1$estimates_imp[,2]-0.3)^2),(mean(cox_m_w1$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_w1     <- c(mean((cox_m_w1$estimates_overimp[,1]+0.3)^2),(mean(cox_m_w1$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_w1     <- c(mean((cox_m_w1$estimates_overimp[,2]-0.3)^2),(mean(cox_m_w1$estimates_overimp[,2]-0.3)))

summaryMSE_w1 <- cbind(rbind(mse_0_x1_org_w1,mse_0_x1_overimp_w1,mse_1_x1_org_w1,mse_1_x1_imp_w1,mse_1_x1_overimp_w1),
                    rbind(mse_0_x2_org_w1,mse_0_x2_overimp_w1,mse_1_x2_org_w1,mse_1_x2_imp_w1,mse_1_x2_overimp_w1))
colnames(summaryMSE_w1) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_w1) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")
write.csv(summaryMSE_w1,file=paste(path_results, "/supp_table_l_m.csv", sep=""))



#########################################################################################################################

# Wilcoxon Tests, referred to in the paper, p.632
wilcox.test(cox_nm$estimates_overimp[,1],cox_nm$estimates_org[,1],paired=TRUE)
wilcox.test(cox_nm$estimates_overimp[,2],cox_nm$estimates_org[,2],paired=TRUE)

wilcox.test(cox_m$estimates_overimp[,1],cox_m$estimates_org[,1],paired=TRUE)
wilcox.test(cox_m$estimates_overimp[,2],cox_m$estimates_org[,2],paired=TRUE)
wilcox.test(cox_m$estimates_overimp[,1],cox_m$estimates_imp[,1],paired=TRUE)
wilcox.test(cox_m$estimates_overimp[,2],cox_m$estimates_imp[,2],paired=TRUE)

########################################################################################################################

##########################
# RESULTS BY SAMPLE SIZE #
##########################

# 500
mse_0_x1_org_s1         <- c(mean((cox_nm_s1$estimates_org[,1]+0.3)^2),(mean(cox_nm_s1$estimates_org[,1]+0.3)))
mse_0_x2_org_s1         <- c(mean((cox_nm_s1$estimates_org[,2]-0.3)^2),(mean(cox_nm_s1$estimates_org[,2]-0.3)))

mse_0_x1_overimp_s1     <- c(mean((cox_nm_s1$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_s1$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_s1     <- c(mean((cox_nm_s1$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_s1$estimates_overimp[,2]-0.3)))

mse_1_x1_org_s1         <- c(mean((cox_m_s1$estimates_org[,1]+0.3)^2),(mean(cox_m_s1$estimates_org[,1]+0.3)))
mse_1_x2_org_s1         <- c(mean((cox_m_s1$estimates_org[,2]-0.3)^2),(mean(cox_m_s1$estimates_org[,2]-0.3)))

mse_1_x1_imp_s1         <- c(mean((cox_m_s1$estimates_imp[,1]+0.3)^2),(mean(cox_m_s1$estimates_imp[,1]+0.3)))
mse_1_x2_imp_s1         <- c(mean((cox_m_s1$estimates_imp[,2]-0.3)^2),(mean(cox_m_s1$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_s1     <- c(mean((cox_m_s1$estimates_overimp[,1]+0.3)^2),(mean(cox_m_s1$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_s1     <- c(mean((cox_m_s1$estimates_overimp[,2]-0.3)^2),(mean(cox_m_s1$estimates_overimp[,2]-0.3)))


summaryMSE_s1 <- cbind(rbind(mse_0_x1_org_s1,mse_0_x1_overimp_s1,mse_1_x1_org_s1,mse_1_x1_imp_s1,mse_1_x1_overimp_s1),
                    rbind(mse_0_x2_org_s1,mse_0_x2_overimp_s1,mse_1_x2_org_s1,mse_1_x2_imp_s1,mse_1_x2_overimp_s1))
colnames(summaryMSE_s1) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_s1) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")

# 1000
mse_0_x1_org_s2         <- c(mean((cox_nm_s2$estimates_org[,1]+0.3)^2),(mean(cox_nm_s2$estimates_org[,1]+0.3)))
mse_0_x2_org_s2         <- c(mean((cox_nm_s2$estimates_org[,2]-0.3)^2),(mean(cox_nm_s2$estimates_org[,2]-0.3)))

mse_0_x1_overimp_s2     <- c(mean((cox_nm_s2$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_s2$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_s2     <- c(mean((cox_nm_s2$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_s2$estimates_overimp[,2]-0.3)))

mse_1_x1_org_s2         <- c(mean((cox_m_s2$estimates_org[,1]+0.3)^2),(mean(cox_m_s2$estimates_org[,1]+0.3)))
mse_1_x2_org_s2         <- c(mean((cox_m_s2$estimates_org[,2]-0.3)^2),(mean(cox_m_s2$estimates_org[,2]-0.3)))

mse_1_x1_imp_s2         <- c(mean((cox_m_s2$estimates_imp[,1]+0.3)^2),(mean(cox_m_s2$estimates_imp[,1]+0.3)))
mse_1_x2_imp_s2         <- c(mean((cox_m_s2$estimates_imp[,2]-0.3)^2),(mean(cox_m_s2$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_s2     <- c(mean((cox_m_s2$estimates_overimp[,1]+0.3)^2),(mean(cox_m_s2$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_s2     <- c(mean((cox_m_s2$estimates_overimp[,2]-0.3)^2),(mean(cox_m_s2$estimates_overimp[,2]-0.3)))


summaryMSE_s2 <- cbind(rbind(mse_0_x1_org_s2,mse_0_x1_overimp_s2,mse_1_x1_org_s2,mse_1_x1_imp_s2,mse_1_x1_overimp_s2),
                    rbind(mse_0_x2_org_s2,mse_0_x2_overimp_s2,mse_1_x2_org_s2,mse_1_x2_imp_s2,mse_1_x2_overimp_s2))
colnames(summaryMSE_s2) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_s2) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")

# 2500
mse_0_x1_org_s3         <- c(mean((cox_nm_s3$estimates_org[,1]+0.3)^2),(mean(cox_nm_s3$estimates_org[,1]+0.3)))
mse_0_x2_org_s3         <- c(mean((cox_nm_s3$estimates_org[,2]-0.3)^2),(mean(cox_nm_s3$estimates_org[,2]-0.3)))

mse_0_x1_overimp_s3     <- c(mean((cox_nm_s3$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_s3$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_s3     <- c(mean((cox_nm_s3$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_s3$estimates_overimp[,2]-0.3)))

mse_1_x1_org_s3         <- c(mean((cox_m_s3$estimates_org[,1]+0.3)^2),(mean(cox_m_s3$estimates_org[,1]+0.3)))
mse_1_x2_org_s3         <- c(mean((cox_m_s3$estimates_org[,2]-0.3)^2),(mean(cox_m_s3$estimates_org[,2]-0.3)))

mse_1_x1_imp_s3         <- c(mean((cox_m_s3$estimates_imp[,1]+0.3)^2),(mean(cox_m_s3$estimates_imp[,1]+0.3)))
mse_1_x2_imp_s3         <- c(mean((cox_m_s3$estimates_imp[,2]-0.3)^2),(mean(cox_m_s3$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_s3     <- c(mean((cox_m_s3$estimates_overimp[,1]+0.3)^2),(mean(cox_m_s3$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_s3     <- c(mean((cox_m_s3$estimates_overimp[,2]-0.3)^2),(mean(cox_m_s3$estimates_overimp[,2]-0.3)))


summaryMSE_s3 <- cbind(rbind(mse_0_x1_org_s3,mse_0_x1_overimp_s3,mse_1_x1_org_s3,mse_1_x1_imp_s3,mse_1_x1_overimp_s3),
                    rbind(mse_0_x2_org_s3,mse_0_x2_overimp_s3,mse_1_x2_org_s3,mse_1_x2_imp_s3,mse_1_x2_overimp_s3))
colnames(summaryMSE_s3) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_s3) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")

# 7500
mse_0_x1_org_s4         <- c(mean((cox_nm_s4$estimates_org[,1]+0.3)^2),(mean(cox_nm_s4$estimates_org[,1]+0.3)))
mse_0_x2_org_s4         <- c(mean((cox_nm_s4$estimates_org[,2]-0.3)^2),(mean(cox_nm_s4$estimates_org[,2]-0.3)))

mse_0_x1_overimp_s4     <- c(mean((cox_nm_s4$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_s4$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_s4     <- c(mean((cox_nm_s4$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_s4$estimates_overimp[,2]-0.3)))

mse_1_x1_org_s4         <- c(mean((cox_m_s4$estimates_org[,1]+0.3)^2),(mean(cox_m_s4$estimates_org[,1]+0.3)))
mse_1_x2_org_s4         <- c(mean((cox_m_s4$estimates_org[,2]-0.3)^2),(mean(cox_m_s4$estimates_org[,2]-0.3)))

mse_1_x1_imp_s4         <- c(mean((cox_m_s4$estimates_imp[,1]+0.3)^2),(mean(cox_m_s4$estimates_imp[,1]+0.3)))
mse_1_x2_imp_s4         <- c(mean((cox_m_s4$estimates_imp[,2]-0.3)^2),(mean(cox_m_s4$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_s4     <- c(mean((cox_m_s4$estimates_overimp[,1]+0.3)^2),(mean(cox_m_s4$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_s4     <- c(mean((cox_m_s4$estimates_overimp[,2]-0.3)^2),(mean(cox_m_s4$estimates_overimp[,2]-0.3)))


summaryMSE_s4 <- cbind(rbind(mse_0_x1_org_s4,mse_0_x1_overimp_s4,mse_1_x1_org_s4,mse_1_x1_imp_s4,mse_1_x1_overimp_s4),
                    rbind(mse_0_x2_org_s4,mse_0_x2_overimp_s4,mse_1_x2_org_s4,mse_1_x2_imp_s4,mse_1_x2_overimp_s4))
colnames(summaryMSE_s4) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_s4) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")

# 10000
mse_0_x1_org_s5         <- c(mean((cox_nm_s5$estimates_org[,1]+0.3)^2),(mean(cox_nm_s5$estimates_org[,1]+0.3)))
mse_0_x2_org_s5         <- c(mean((cox_nm_s5$estimates_org[,2]-0.3)^2),(mean(cox_nm_s5$estimates_org[,2]-0.3)))

mse_0_x1_overimp_s5     <- c(mean((cox_nm_s5$estimates_overimp[,1]+0.3)^2),(mean(cox_nm_s5$estimates_overimp[,1]+0.3)))
mse_0_x2_overimp_s5     <- c(mean((cox_nm_s5$estimates_overimp[,2]-0.3)^2),(mean(cox_nm_s5$estimates_overimp[,2]-0.3)))

mse_1_x1_org_s5         <- c(mean((cox_m_s5$estimates_org[,1]+0.3)^2),(mean(cox_m_s5$estimates_org[,1]+0.3)))
mse_1_x2_org_s5         <- c(mean((cox_m_s5$estimates_org[,2]-0.3)^2),(mean(cox_m_s5$estimates_org[,2]-0.3)))

mse_1_x1_imp_s5         <- c(mean((cox_m_s5$estimates_imp[,1]+0.3)^2),(mean(cox_m_s5$estimates_imp[,1]+0.3)))
mse_1_x2_imp_s5         <- c(mean((cox_m_s5$estimates_imp[,2]-0.3)^2),(mean(cox_m_s5$estimates_imp[,2]-0.3)))

mse_1_x1_overimp_s5     <- c(mean((cox_m_s5$estimates_overimp[,1]+0.3)^2),(mean(cox_m_s5$estimates_overimp[,1]+0.3)))
mse_1_x2_overimp_s5     <- c(mean((cox_m_s5$estimates_overimp[,2]-0.3)^2),(mean(cox_m_s5$estimates_overimp[,2]-0.3)))


summaryMSE_s5 <- cbind(rbind(mse_0_x1_org_s5,mse_0_x1_overimp_s5,mse_1_x1_org_s5,mse_1_x1_imp_s5,mse_1_x1_overimp_s5),
                    rbind(mse_0_x2_org_s5,mse_0_x2_overimp_s5,mse_1_x2_org_s5,mse_1_x2_imp_s5,mse_1_x2_overimp_s5))
colnames(summaryMSE_s5) <- c("MSE X1","Bias X1","MSE X2","Bias X2")
rownames(summaryMSE_s5) <- c("nm: CC","nm: MO","m: CC","m: MI","m: MO")

summaryMSE_s1
summaryMSE_s2
summaryMSE_s3
summaryMSE
summaryMSE_s4
summaryMSE_s5

sample_size <- c(500,1000,2500,5000,7500,10000)

Bias_X1_MO_nm <- c(summaryMSE_s1[2,2],summaryMSE_s2[2,2],summaryMSE_s3[2,2],summaryMSE[2,2],summaryMSE_s4[2,2],summaryMSE_s5[2,2])
Bias_X1_CC_nm <- c(summaryMSE_s1[1,2],summaryMSE_s2[1,2],summaryMSE_s3[1,2],summaryMSE[1,2],summaryMSE_s4[1,2],summaryMSE_s5[1,2])

Bias_X2_MO_nm <- c(summaryMSE_s1[2,4],summaryMSE_s2[2,4],summaryMSE_s3[2,4],summaryMSE[2,4],summaryMSE_s4[2,4],summaryMSE_s5[2,4])
Bias_X2_CC_nm <- c(summaryMSE_s1[1,4],summaryMSE_s2[1,4],summaryMSE_s3[1,4],summaryMSE[1,4],summaryMSE_s4[1,4],summaryMSE_s5[1,4])

MSE_X1_MO_nm <- c(summaryMSE_s1[2,1],summaryMSE_s2[2,1],summaryMSE_s3[2,1],summaryMSE[2,1],summaryMSE_s4[2,1],summaryMSE_s5[2,1])
MSE_X1_CC_nm <- c(summaryMSE_s1[1,1],summaryMSE_s2[1,1],summaryMSE_s3[1,1],summaryMSE[1,1],summaryMSE_s4[1,1],summaryMSE_s5[1,1])

MSE_X2_MO_nm <- c(summaryMSE_s1[2,3],summaryMSE_s2[2,3],summaryMSE_s3[2,3],summaryMSE[2,3],summaryMSE_s4[2,3],summaryMSE_s5[2,3])
MSE_X2_CC_nm <- c(summaryMSE_s1[1,3],summaryMSE_s2[1,3],summaryMSE_s3[1,3],summaryMSE[1,3],summaryMSE_s4[1,3],summaryMSE_s5[1,3])


Bias_X1_MO_m <- c(summaryMSE_s1[5,2],summaryMSE_s2[5,2],summaryMSE_s3[5,2],summaryMSE[5,2],summaryMSE_s4[5,2],summaryMSE_s5[5,2])
Bias_X1_MI_m <- c(summaryMSE_s1[4,2],summaryMSE_s2[4,2],summaryMSE_s3[4,2],summaryMSE[4,2],summaryMSE_s4[4,2],summaryMSE_s5[4,2])
Bias_X1_CC_m <- c(summaryMSE_s1[3,2],summaryMSE_s2[3,2],summaryMSE_s3[3,2],summaryMSE[3,2],summaryMSE_s4[3,2],summaryMSE_s5[3,2])

MSE_X1_MO_m <- c(summaryMSE_s1[5,1],summaryMSE_s2[5,1],summaryMSE_s3[5,1],summaryMSE[5,1],summaryMSE_s4[5,1],summaryMSE_s5[5,1])
MSE_X1_MI_m <- c(summaryMSE_s1[4,1],summaryMSE_s2[4,1],summaryMSE_s3[4,1],summaryMSE[4,1],summaryMSE_s4[4,1],summaryMSE_s5[4,1])
MSE_X1_CC_m <- c(summaryMSE_s1[3,1],summaryMSE_s2[3,1],summaryMSE_s3[3,1],summaryMSE[3,1],summaryMSE_s4[3,1],summaryMSE_s5[3,1])


Bias_X2_MO_m <- c(summaryMSE_s1[5,4],summaryMSE_s2[5,4],summaryMSE_s3[5,4],summaryMSE[5,4],summaryMSE_s4[5,4],summaryMSE_s5[5,4])
Bias_X2_MI_m <- c(summaryMSE_s1[4,4],summaryMSE_s2[4,4],summaryMSE_s3[4,4],summaryMSE[4,4],summaryMSE_s4[4,4],summaryMSE_s5[4,4])
Bias_X2_CC_m <- c(summaryMSE_s1[3,4],summaryMSE_s2[3,4],summaryMSE_s3[3,4],summaryMSE[3,4],summaryMSE_s4[3,4],summaryMSE_s5[3,4])

MSE_X2_MO_m <- c(summaryMSE_s1[5,3],summaryMSE_s2[5,3],summaryMSE_s3[5,3],summaryMSE[5,3],summaryMSE_s4[5,3],summaryMSE_s5[5,3])
MSE_X2_MI_m <- c(summaryMSE_s1[4,3],summaryMSE_s2[4,3],summaryMSE_s3[4,3],summaryMSE[4,3],summaryMSE_s4[4,3],summaryMSE_s5[4,3])
MSE_X2_CC_m <- c(summaryMSE_s1[3,3],summaryMSE_s2[3,3],summaryMSE_s3[3,3],summaryMSE[3,3],summaryMSE_s4[3,3],summaryMSE_s5[3,3])

#  X1, no missing data
pdf(file=paste(path_results, "/main_sim_sample_size_beta1_nomiss.pdf", sep=""),width=12)
par(mfrow=c(1,2))
plot(sample_size[-1] ,Bias_X1_MO_nm[-1],type="l",col="red",lty=3,lwd=3,xlim=c(1000,10000),ylim=c(0,0.08),xlab="Sample Size",ylab="Bias")
lines(sample_size[-1],Bias_X1_CC_nm[-1],type="l",lwd=3,col="black")
legend("topright",bty="n",col=c("black","red"),legend=c("Complete Cases","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,3))
plot(sample_size[-1] ,MSE_X1_MO_nm[-1],type="l",col="red",lwd=3,lty=3,xlim=c(1000,10000),ylim=c(0,0.01),xlab="Sample Size",ylab="MSE")
lines(sample_size[-1],MSE_X1_CC_nm[-1],type="l",col="black",lwd=3)
legend("topright",bty="n",col=c("black","red"),legend=c("Complete Cases","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,3))
dev.off()

#  X2, no missing data
pdf(file=paste(path_results, "/main_sim_sample_size_beta2_nomiss.pdf", sep=""),width=12)
par(mfrow=c(1,2))
plot(sample_size[-1] ,Bias_X2_MO_nm[-1],type="l",col="red",lwd=3,lty=3,xlim=c(1000,10000),ylim=c(-0.075,0),xlab="Sample Size",ylab="Bias")
lines(sample_size[-1],Bias_X2_CC_nm[-1],type="l",col="black",lwd=3)
legend("bottomright",bty="n",col=c("black","red"),legend=c("Complete Cases","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,3))
plot(sample_size[-1] ,MSE_X2_MO_nm[-1],type="l",col="red",lwd=3,lty=3,xlim=c(1000,10000),ylim=c(0,0.015),xlab="Sample Size",ylab="MSE")
lines(sample_size[-1],MSE_X2_CC_nm[-1],type="l",col="black",lwd=3)
legend("topright",bty="n",col=c("black","red"),legend=c("Complete Cases","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,3))
dev.off()

#  X1, missing data   - REPORTED IN PAPER
pdf(file=paste(path_results, "/_FIGURE1.pdf", sep=""),width=12)
par(mfrow=c(1,2))
plot(sample_size[-1] ,Bias_X1_MO_m[-1],type="l",col="red",lty=3,lwd=3,xlim=c(1000,10000),ylim=c(0,0.06),xlab="Sample Size",ylab="Bias")
lines(sample_size[-1],Bias_X1_CC_m[-1],type="l",lwd=3,col="black")
lines(sample_size[-1],Bias_X1_MI_m[-1],type="l",lwd=3,lty=2,col="royalblue1")  
legend("topright",bty="n",col=c("black","royalblue1","red"),legend=c("Complete Cases","Multiple Imputation","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,2,3))
plot(sample_size[-1] ,MSE_X1_MO_m[-1],type="l",col="red",lwd=3,lty=3,xlim=c(1000,10000),ylim=c(0,0.006),xlab="Sample Size",ylab="MSE")
lines(sample_size[-1],MSE_X1_CC_m[-1],type="l",col="black",lwd=3)
lines(sample_size[-1],MSE_X1_MI_m[-1],type="l",col="royalblue1",lwd=3,lty=2)
legend("topright",bty="n",col=c("black","royalblue1","red"),legend=c("Complete Cases","Multiple Imputation","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,2,3))
dev.off()

#  X2, missing data
pdf(file=paste(path_results, "/main_sim_sample_size_beta2_miss.pdf", sep=""),width=12)
par(mfrow=c(1,2))
plot(sample_size[-1] ,Bias_X2_MO_m[-1],type="l",col="red",lty=3,lwd=3,xlim=c(1000,10000),ylim=c(-0.08,0),xlab="Sample Size",ylab="Bias")
lines(sample_size[-1],Bias_X2_CC_m[-1],type="l",lwd=3,col="black")
lines(sample_size[-1],Bias_X2_MI_m[-1],type="l",lwd=3,lty=2,col="royalblue1")  
legend("bottomright",bty="n",col=c("black","royalblue1","red"),legend=c("Complete Cases","Multiple Imputation","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,2,3))
plot(sample_size[-1] ,MSE_X2_MO_m[-1],type="l",col="red",lwd=3,lty=3,xlim=c(1000,10000),ylim=c(0,0.015),xlab="Sample Size",ylab="MSE")
lines(sample_size[-1],MSE_X2_CC_m[-1],type="l",col="black",lwd=3)
lines(sample_size[-1],MSE_X2_MI_m[-1],type="l",col="royalblue1",lwd=3,lty=2)
legend("topright",bty="n",col=c("black","royalblue1","red"),legend=c("Complete Cases","Multiple Imputation","Multiple Overimputation"), lwd=3, cex=1.25,lty=c(1,2,3))
dev.off()

####################
# TABLE FROM PAPER #
####################

rbind(
t(summaryMSE)[c(2,4),],
t(summaryMSE)[c(1,3),]
)



#########################################################################################################################
save.image(paste(path_results, "/results.Rdata", sep=""))
#########################################################################################################################




