#############################################################################
# Supplementary Material to Schomaker M, Davies MA, Cornell M, Ford N.      #
# Assessing the risk of dolutegravir for women of childbearing potential.   #
# Lancet Global Health. 2018;6(9):e958-e9.                                  #                                             
#                                                                           #                                         
# R-Code for reproduction of figure                                         #
#############################################################################


# P(X>=4)
1-pbinom(3,426,0.001)

#
pfunc <- function(xv){1-pbinom(3,xv,0.001)}
pfunc.1 <- function(xv){1-pbinom(4,xv,0.001)}
pfunc.2 <- function(xv){1-pbinom(5,xv,0.001)}
pfunc.3 <- function(xv){1-pbinom(6,xv,0.001)}
pfunc.4 <- function(xv){1-pbinom(7,xv,0.001)}
pfunc.5 <- function(xv){1-pbinom(8,xv,0.001)}
pfunc.6 <- function(xv){1-pbinom(9,xv,0.001)}

# not reported in paper
plot(c(400:4000),pfunc(c(400:4000)),type="l",xlab="number of patients",ylab="P(4 or more events)",cex.lab=1.5,lwd=2)
title(main = "Probability of 4 or more events\n depending on patients number")
abline(h=0.05,col="red",lwd=2)


# not reported in paper
plot(c(400:4000),pfunc(c(400:4000)),type="l",xlab="number of pregnant women",ylab="P(x or more events)",cex.lab=1.5,lwd=2)
lines(c(400:4000),pfunc.1(c(400:4000)),type="l",cex.lab=1.5,lwd=2,col=3,lty=2)
lines(c(400:4000),pfunc.2(c(400:4000)),type="l",cex.lab=1.5,lwd=2,col=4,lty=3)
lines(c(400:4000),pfunc.3(c(400:4000)),type="l",cex.lab=1.5,lwd=2,col=5,lty=4)
lines(c(400:4000),pfunc.4(c(400:4000)),type="l",cex.lab=1.5,lwd=2,col=6,lty=5)
lines(c(400:4000),pfunc.5(c(400:4000)),type="l",cex.lab=1.5,lwd=2,col=7,lty=6)
lines(c(400:4000),pfunc.6(c(400:4000)),type="l",cex.lab=1.5,lwd=2,col=8,lty=7)
title(main = "Probability of x or more adverse events\n depending on number of patients")
abline(h=0.05,col="red",lwd=2)
legend("topleft",c("4 events","5 events","6 events","7 events","8 events","9 events","10 events"),col=c(1,3:8),lty=1:7)


# not reported in paper
plot(c(250:2500),pfunc(c(250:2500)),type="l",xlab="number of pregnant women",ylab="P(x or more events)",cex.lab=1.5,lwd=2,axes=F,ylim=c(0,0.12))
axis(side = 1, at = c(250,426,1000,1426,2000,2500))
axis(side = 2, at = c(0,0.02,0.04,0.06,0.08,0.1,0.12),labels=c("","2%","4%","6%","8%","10%","12%"))
lines(c(250:2500),pfunc.1(c(250:2500)),type="l",cex.lab=1.5,lwd=2,col=3,lty=2)
lines(c(250:2500),pfunc.2(c(250:2500)),type="l",cex.lab=1.5,lwd=2,col=4,lty=3)
lines(c(250:2500),pfunc.3(c(250:2500)),type="l",cex.lab=1.5,lwd=2,col=5,lty=4)
lines(c(250:2500),pfunc.4(c(250:2500)),type="l",cex.lab=1.5,lwd=2,col=6,lty=5)
points(426,pfunc(426),cex=2.5,pch=19)
points(1426,pfunc(1426),cex=2.5,pch=17,col="red")
points(1426,pfunc.2(1426),cex=2.5,pch=18,col="blue")  
abline(h=0.05,col="red",lwd=2,lty=2)     
legend("left",c("current situation","+1000 patients / +no event","+1000 patients /+2 events"),col=c("black","red","blue"),pch=c(19,17,18),cex=1.1,bty="n")
legend("topleft",c("4 or more events","5 or more events","6 or more events","7 or more events","8 or more events"),col=c(1,3:6),lty=1:5,bty="n")


# binom.test
binom.test(c(4,422),p=0.001,alternative="greater")
pfunc2 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(4,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}
pfunc2.1 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(5,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}
pfunc2.2 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(6,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}
pfunc2.3 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(7,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}
pfunc2.4 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(8,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}
pfunc2.5 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(9,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}
pfunc2.6 <- function(xv){
results <- rep(NA,length(xv))
for(i in 1:length(xv)){results[i] <- binom.test(c(10,xv[i]),p=0.001,alternative="greater")$conf.int[1]}
return(results)
}


plot(c(400:4004),pfunc2(c(396:4000)),type="l",xlab="number of patients",ylab="lower conf. limit",cex.lab=1.5,lwd=2)
title(main = "Lower confidence limit \n depending on patients number")
abline(h=0.001,col="red",lwd=2)


#
plot(c(400:4004),pfunc2(c(396:4000)),type="l",xlab="number of patients",ylab="lower conf. limit",cex.lab=1.5,lwd=2)
lines(c(400:4004),pfunc2.1(c(396:4000)),type="l",cex.lab=1.5,lwd=2,col=3,lty=2)
lines(c(400:4004),pfunc2.2(c(396:4000)),type="l",cex.lab=1.5,lwd=2,col=4,lty=3)
lines(c(400:4004),pfunc2.3(c(396:4000)),type="l",cex.lab=1.5,lwd=2,col=5,lty=4)
lines(c(400:4004),pfunc2.4(c(396:4000)),type="l",cex.lab=1.5,lwd=2,col=6,lty=5)
lines(c(400:4004),pfunc2.5(c(396:4000)),type="l",cex.lab=1.5,lwd=2,col=7,lty=6)
lines(c(400:4004),pfunc2.6(c(396:4000)),type="l",cex.lab=1.5,lwd=2,col=8,lty=7)
title(main = "Lower 95% confidence limit \n depending on number of patients and events")
abline(h=0.001,col="red",lwd=2)
legend("topright",c("4 events","5 events","6 events","7 events","8 events","9 events","10 events"),col=c(1,3:8),lty=1:7)


# FIGURE FROM PAPER
plot(c(250:2504),pfunc2(c(246:2500)),type="l",xlab="number of pregnant women",ylab="lower 95% confidence limit",cex.lab=1.5,lwd=2,axes=F,ylim=c(0,0.007))
lines(c(250:2504),pfunc2.1(c(246:2500)),type="l",cex.lab=1.5,lwd=2,col=3,lty=2)
lines(c(250:2504),pfunc2.2(c(246:2500)),type="l",cex.lab=1.5,lwd=2,col=4,lty=3)
lines(c(250:2504),pfunc2.3(c(246:2500)),type="l",cex.lab=1.5,lwd=2,col=5,lty=4)
lines(c(250:2504),pfunc2.4(c(246:2500)),type="l",cex.lab=1.5,lwd=2,col=6,lty=5)
abline(h=0.001,col="red",lwd=2)
legend(2000,0.0045,c("4 events","5 events","6 events","7 events","8 events"),col=c(1,3:6),lty=1:5,lwd=2)
axis(side = 1, at = c(250,426,1000,1426,2000,2500))
axis(side = 2, at = c(0,0.001,0.002,0.003,0.004,0.005,0.006,0.007),labels=c("","0.1%","0.2%","0.3%","0.4%","0.5%","0.6%","0.7%"))
points(426,pfunc2(422),cex=2.5,pch=19)
points(1426,pfunc2(1422),cex=2.5,pch=17,col="red")
points(1426,pfunc2.2(1422),cex=2.5,pch=18,col="blue")       
legend("topright",c("current situation","+1000 patients / +no event","+1000 patients /+2 events"),col=c("black","red","blue"),pch=c(19,17,18),cex=1.1,bty="n")





