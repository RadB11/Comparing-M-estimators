mad <- function(error)
{
  MAD <- median(abs(error-median(error)))
  (1.483*MAD)^2 
}

######### Huber psi function #################
H.psi <- function(error, tune) {
  error*(abs(error) <= tune) + sign(error)*tune*(abs(error) > tune)  
}

######### Huber weight function #################

H.wt <- function(error, tune) {
  1*(abs(error) <= tune)  
}

######### derivative of Huber's function

H.deri <- function(error, tune) {
  1*(abs(error) <= tune)
}

##################################
tau.hat.huber <- function(s, tune) {
  n <- length(s)
  tau <- sum(abs(s) <= tune)^2/(n*sum((abs(s) <= tune)*s^2 + tune^2 *(abs(s) > tune)))
}

#-------------------------Tukey's biweight function----

######### Tukey's psi function #################

T.psi <- function(error, tune) {
  error*(1-(error/tune)^2)^2*(abs(error) <= tune)
}

######### Tukey's weight function #################

T.wt <- function(error, tune) {
  (1-(error/tune)^2)^2*(abs(error) <= tune)
}

######### derivative of Tukey's function

T.deri <- function(error, tune) {
  (1-(error/tune)^2)*(1-5*(error/tune)^2)*(abs(error)<=tune)
}

##################################################
tau.hat.tukey <- function(s, tune) {
  tau <- mean(T.deri(s, tune))^2/(mean((abs(s) <= tune)*T.psi(s,tune)^2))
}
#--------------- Andrew's function -------------

A.psi <- function(error, tune) {
  tune*sin(error/tune)*(abs(error) <= tune*pi)
}


A.wt <- function(error, tune) {
  tune*(sin(error/tune)*(abs(error) <= tune*pi))/error
}

######### derivative of Andrew's function

A.deri <- function(error, tune) {
  cos(error/tune)*(abs(error) <= tune*pi)
}

##################################################
tau.hat.Andrew <- function(s, tune) {
  tau <- mean(A.deri(s, tune))^2/(mean((abs(s) <= tune*pi)*A.psi(s,tune)^2))
}

################# Alamgir's function

# psi function
Alam.psi <- function(error, tune) {
  ((16*error*exp(-2*(error/tune)^2))/(1 + exp(-(error/tune)^2))^4)*(abs(error) <= tune)
}

# weight function
Alam.wt <- function(error, tune) {
  ((16*exp(-2*(error/tune)^2))/(1 + exp(-(error/tune)^2))^4)*(abs(error) <= tune)
}

# 2nd derivative 
Alam.deri <- function(error, tune){
  ((16*exp(-2*(error/tune)^2)*(1- 4*(error/tune)^2*(1+2*exp(-(error/tune)^2))))/(1 + exp(-(error/tune)^2))^5)*(abs(error) <= tune)
}
################################################ tau hat
tau.hat.Alam <- function(s, tune) {
  tau <- mean(Alam.deri(s, tune))^2/(mean((abs(s) <= tune)*Alam.psi(s,tune)^2))
}


###################### Insha's function
Insha.psi <- function(error, tune) {
  (error/(1+(error/tune)^4)^2)*(abs(error) >= 0)
}

# weight function
Insha.wt <- function(error, tune) {
  (1/(1+(error/tune)^4)^2)*(abs(error) >= 0)
}

# 2nd derivative 
Insha.deri <- function(error, tune){
  ((((1+(error/tune)^2) - 8*(error/tune)^4))/((1+(error/tune)^4)^3)) *(abs(error) >= 0)
}
##################################

tau.hat.Insha <- function(s, tune) {
  tau <- mean(Insha.deri(s, tune))^2/(mean((abs(s) >= 0)*Insha.psi(s,tune)^2))
}

#################### Asad's function

Asad.psi <- function(error, tune){
  ((2*error/3) * (1-(error/tune)^4)^2)*(abs(error) <= tune)
}

# weight function
Asad.wt <- function(error, tune){
  ((2/3)*(1-(error/tune)^4)^2)*(abs(error) <= tune)
}

# derivative

Asad.deri <- function(error, tune){
  (0.6666667*((1-(error/tune)^4)^2 - 8*(error/tune)^4*(1-(error/tune)^4)))*(abs(error) <= tune)
}

#####################################

tau.hat.Asad <- function(s, tune) {
  tau <- mean(Asad.deri(s, tune))^2/(mean((abs(s) <= tune)*Asad.psi(s,tune)^2))
}

######################## Qader beta function
Qader.beta.psi <- function(error, tune){
  ((error/16)*(1-(error/tune)^2)^2 ) *(abs(error) <= tune)
}

# weight function
Qader.beta.wt <- function(error, tune){
  ((1/16)*(1-(error/tune)^2)^2 ) *(abs(error) <= tune)
}

# derivative

Qader.beta.deri <- function(error, tune){
  (0.0625*(1-(error/tune)^2)*(-4*(error/tune)^2 + (1-(error/tune)^2)))*(abs(error) <= tune)
}

#####################################

tau.hat.beta.Qader <- function(s, tune) {
  tau <- mean(Qader.beta.deri(s, tune))^2/(mean((abs(s) <= tune)*Qader.beta.psi(s,tune)^2))
}






n <- 10000

B0 <- 1
B1 <- 1


#tune = 2.1


################ Tua hat
x_1 <- rnorm(n)
e <- rt(n, df=5)
#e <- rnorm(n)
y = B0 + B1*x_1 + e
wts <- rep(1,n)

tune_range <- seq(0.4,10,.1)
tau.h <- vector()
for (indx in 1:length(tune_range)){
  for (itr in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- H.wt(err,tune_range[indx])
  }
  
  tau.h[indx] <- tau.hat.huber(err,tune_range[indx])
  
}

pdf("TauH.pdf")
par(mar=c(5, 5, 2, 1))
plot(tune_range,tau.h,type = "l",main = "Huber" , mgp = c(2.4, 1, 0), xlab ="", xlim=c(0,10),
     las=1, lwd=2, cex.lab=2, cex.axis=2, ylab="")
title(xlab="c", cex.lab=2,line = 2.5)
title(ylab=expression(hat(tau)), cex.lab=2,line = 3)
dev.off()
#--------------- for Tukey ----------
wts <- rep(1,n)
tau.T <- vector()
for (indx in 1:length(tune_range)){
  for (itr in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- T.wt(err,tune_range[indx])
  }
  
  tau.T[indx] <- tau.hat.tukey(err,tune_range[indx])
  
}

pdf("TauT1.pdf")
par(mar=c(5, 5, 2, 1))
plot(tune_range,tau.T,type = "l", main = "Tukey", mgp = c(2.4, 1, 0), xlab ="", xlim=c(0,10),
     las=1, lwd=2, cex.lab=2, cex.axis=2, ylab="")
title(xlab="c", cex.lab=2,line = 2.5)
title(ylab=expression(hat(tau)), cex.lab=2,line = 3)
dev.off()

#-------------- for Andrew's ----------
wts <- rep(1,n)
tau.A <- vector()
for (indx in 1:length(tune_range)){
  for (itr in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- A.wt(err,tune_range[indx])
  }
  
  tau.A[indx] <- tau.hat.Andrew(err,tune_range[indx])
  
}

pdf("TauA.pdf")
par(mar=c(5, 5, 2, 1))
plot(tune_range,tau.A,type = "l",main = "Andrew", mgp = c(2.4, 1, 0), xlab ="", xlim=c(0,10),
     las=1, lwd=2, cex.lab=2, cex.axis=2, ylab="")
title(xlab="c", cex.lab=2,line = 2.5)
title(ylab=expression(hat(tau)), cex.lab=2,line = 3)
dev.off()

# pdf("TauHTA.pdf")
# par(mar=c(5, 5, 2, 1))
# plot(tune_range,tau.A,type = "l", col = "red", ylim = c(0,1.5),main = "Huber, Tukey and Andrew efficiency factors for c", mgp = c(2.4, 1, 0), xlab ="", xlim=c(0,5),
#      las=1, lwd=2, cex.lab=2, cex.axis=2, ylab="")
# points(tune_range,tau.T,type = "l", col = "green" , main = "Tukey", mgp = c(2.4, 1, 0), xlab ="", xlim=c(0,5),
#        las=1, lwd=2, cex.lab=2, cex.axis=2, ylab="")
# points(tune_range,tau.h,type = "l", col = "blue", main = "Tukey", mgp = c(2.4, 1, 0), xlab ="", xlim=c(0,5),
#        las=1, lwd=2, cex.lab=2, cex.axis=2, ylab="")
# title(xlab="c", cex.lab=2,line = 2.5)
# title(ylab=expression(hat(tau)), cex.lab=2,line = 3)
# 
# legend("topright", legend = c("Huber", "Tukey", "Andrew"),
#        cex = 2,box.lwd = 1,box.lty = 1,
#        text.font = 1.5,col = c("red","green", "blue"),
#        lty = 1, pch = 16)
# 
# 
# dev.off()


#########################################
# Find maximum tau for c vlaues and 
#then boxplot for several simulations
########################################
n= 5000
tune_range <- seq(0.5,5,.2)
res <- replicate(5,{
  
  x_1 <- rnorm(n)
  e <- rt(n, df=2)
  #e <- rnorm(n)
  y = B0 + B1*x_1 + e
  wts <- rep(1,n)
  
  tau.h <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- H.wt(err,tune_range[indx])
    }
    
    tau.h[indx] <- tau.hat.huber(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.h,type = "l")
  tau.hat.h <- max(tau.h)
  
  #--------------- for Tukey ----------
  wts <- rep(1,n)
  tau.T <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- T.wt(err,tune_range[indx])
    }
    
    tau.T[indx] <- tau.hat.tukey(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.T,type = "l")
  tau.hat.T <- max(tau.T)
  
  #-------------- for Andrew's ----------
  wts <- rep(1,n)
  tau.A <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- A.wt(err,tune_range[indx])
    }
    
    tau.A[indx] <- tau.hat.Andrew(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.A,type = "l")
  tau.hat.A <- max(tau.A)
  
  
  #################### For Alamgir
  
  wts <- rep(1,n)
  tau.Alam <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- Alam.wt(err,tune_range[indx])
    }
    
    tau.Alam[indx] <- tau.hat.Alam(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.A,type = "l")
  tau.hat.Alam <- max(tau.Alam)
  
  
  wts <- rep(1,n)
  tau.Insha <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- Insha.wt(err,tune_range[indx])
    }
    
    tau.Insha[indx] <- tau.hat.Insha(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.A,type = "l")
  tau.hat.Insha <- max(tau.Insha)
  
  
  #------------------------- For Asad
  
  wts <- rep(1,n)
  tau.Asad <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- Asad.wt(err,tune_range[indx])
    }
    
    tau.Asad[indx] <- tau.hat.Asad(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.Asad,type = "l")
  tau.hat.Asad <- max(tau.Asad)
  
  #------------------------------ For Qader beta function
  wts <- rep(1,n)
  tau.beta.Qader <- vector()
  for (indx in 1:length(tune_range)){
    for (itr in 1:10){
      mod <- lm(y~x_1,weights=wts)
      SS <- mad(mod$resid)
      err <- mod$resid/sqrt(SS)
      wts <- Qader.beta.wt(err,tune_range[indx])
    }
    
    tau.beta.Qader[indx] <- tau.hat.beta.Qader(err,tune_range[indx])
    
  }
  
  #plot(tune_range,tau.Qader,type = "l")
  tau.hat.beta.Qader <- max(tau.beta.Qader)
  
  ##############################################
  # Calculating Coefficients for tune with highest tau
  ##############################################
  
  
  #------------- Huber method --------------
  
  
  coeff_H <- rep(NA,2)
  
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- H.wt(err,tune_range[min(which(tau.hat.h==tau.h))])
    coeff_H <- mod$coefficients
  }
  
  # ---------------- Tukey method -------------
  
  coeff_T <- rep(NA, 2)
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- T.wt(err,tune_range[min(which(tau.hat.T==tau.T))])
    coeff_T <- mod$coefficients
  }
  
  #--------------------- Alamgir mwthod-------------
  
  coeff_Alam <- rep(NA, 2)
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- Alam.wt(err,tune_range[min(which(tau.hat.Alam==tau.Alam))])
    coeff_Alam <- mod$coefficients
  }
  
  #------------------------ Insha method
  coeff_Insha <- rep(NA, 2)
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- Insha.wt(err,tune_range[min(which(tau.hat.Insha==tau.Insha))])
    coeff_Insha <- mod$coefficients
  }
  
  #------------------------------ Asad method
  
  coeff_Asad <- rep(NA, 2)
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- Asad.wt(err,tune_range[min(which(tau.hat.Asad==tau.Asad))])
    coeff_Asad <- mod$coefficients
  }
  
  #------------------------- Qader beta method
  
  coeff_beta_Qader <- rep(NA, 2)
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- Qader.beta.wt(err,tune_range[min(which(tau.hat.beta.Qader==tau.beta.Qader))])
    coeff_beta_Qader <- mod$coefficients
  }
  
  #-----------------OLS method-----------
  
  mod <- lm(y~x_1)
  coeff_OLS <- mod$coefficients
  
  #----------- Andrew's method ---------
  
  coeff_Ad <- rep(NA, 2)
  for (i in 1:10){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    wts <- A.wt(err,tune_range[min(which(tau.hat.A==tau.A))])
    #abline(mod)
    coeff_Ad <- mod$coefficients
  }
  
  c(coeff_H[1],coeff_T[1],coeff_Ad[1],coeff_Alam[1],coeff_Insha[1], coeff_Asad[1],coeff_beta_Qader[1],coeff_OLS[1],
    coeff_H[2],coeff_T[2],coeff_Ad[2],coeff_Alam[2], coeff_Insha[2], coeff_Asad[2],coeff_beta_Qader[2],coeff_OLS[2])
})


#pdf(file="box.pdf")
#par(mar=c(5, 5, 2, 1))
boxplot(res~row(res), outline=FALSE, las=2, xaxt="n",at=c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17),ylim=c(.85,1.12),
        cex.axis=2, cex.lab=2, col=c("red","green", "blue", "gray","orange","green3", "coral1", "black"))
abline(h=1)
#points(apply(res,1,mean),pch=16)
mtext(expression(beta[0]), at=2.5, line = 0)
mtext(expression(beta[1]), at=7.5, line = 0)
legend("bottomright", legend = c("Huber", "Tukey", "Andrew","Alamgir", "Insha", "Asad", "Qader-beta",  "OLS"),
       cex = 1.5,box.lwd = 1,box.lty = 1,
       text.font = 2,col = c("red","green", "blue", "gray", "orange", "green3","coral1", "black"),
       pch = 2,horiz = T)
#dev.off()






#############################################
# Obtaining c values on 10 simulations for Huber, Tukey, Andrew's
#############################################


n <- c(50, 100, 500, 1000)

B0 <- 1
B1 <- 1



res <- replicate(20,{
  eff_c.H <- vector()
  for(samp in 1:length(n)){
    x_1 <- rnorm(n[samp])
    e <- rt(n[samp], df=1)
    
    y = B0 + B1*x_1 + e
    wts <- rep(1,n[samp])
    
    tune_range <- seq(0.2,3,.1)
    tau.h <- vector()
    for (indx in 1:length(tune_range)){
      for (itr in 1:10){
        mod <- lm(y~x_1,weights=wts)
        err <- scale(mod$resid)
        wts <- H.wt(err,tune_range[indx])
      }
      
      tau.h[indx] <- tau.hat.huber(err,tune_range[indx])
      
    }
    
    tau.hat.h <- max(tau.h)
    
    eff_c.H[samp] <- tune_range[min(which(tau.hat.h==tau.h))]
  }
  eff_c.H
})

eff_c <- as.vector(res) 
arr <- as.vector(array(c(1:4),dim = c(length(n),20)))


boxplot(eff_c~arr, outline=FALSE, ylim = c(0,1), las=2, xaxt="n",at=c(1,2,3,4),
        cex.axis=2, cex.lab=2, col=c("red","green"))


#--------------- for Tukey ----------

res <- replicate(10,{
  eff_c.T <- vector()
  for(samp in 1:length(n)){
    x_1 <- rnorm(n[samp])
    e <- rt(n[samp], df=1)
    
    y = B0 + B1*x_1 + e
    wts <- rep(1,n[samp])
    
    tau.T <- vector()
    for (indx in 1:length(tune_range)){
      for (itr in 1:10){
        mod <- lm(y~x_1,weights=wts)
        err <- scale(mod$resid)
        wts <- T.wt(err,tune_range[indx])
      }
      
      tau.T[indx] <- tau.hat.tukey(err,tune_range[indx])
      
    }
    
    tau.hat.T <- max(tau.T)
    eff_c.T[samp] <- tune_range[min(which(tau.hat.T==tau.T))]
  }
  eff_c.T
})

eff_c <- as.vector(res) 
arr <- as.vector(array(c(1:4),dim = c(length(n),10)))


boxplot(eff_c~arr, outline=FALSE, ylim = c(0,3), las=2, xaxt="n",at=c(1,2,3,4),
        cex.axis=2, cex.lab=2, col=c("red","green"))



#-------------- for Andrew's ----------

res <- replicate(20,{
  eff_c.A <- vector()
  for(samp in 1:length(n)){
    x_1 <- rnorm(n[samp])
    e <- rt(n[samp], df=1)
    
    y = B0 + B1*x_1 + e
    wts <- rep(1,n[samp])
    
    tau.A <- vector()
    for (indx in 1:length(tune_range)){
      for (itr in 1:10){
        mod <- lm(y~x_1,weights=wts)
        err <- scale(mod$resid)
        wts <- T.wt(err,tune_range[indx])
      }
      
      tau.A[indx] <- tau.hat.Andrew(err,tune_range[indx])
      
    }
    
    tau.hat.A <- max(tau.A)
    
    eff_c.A[samp] <- tune_range[min(which(tau.hat.A==tau.A))]
  }
  eff_c.A
})

eff_c <- as.vector(res) 
arr <- as.vector(array(c(1:4),dim = c(length(n),20)))


boxplot(eff_c~arr, outline=FALSE, ylim = c(0,1), las=2, xaxt="n",at=c(1,2,3,4),
        cex.axis=2, cex.lab=2, col=c("red","green"))


####################################################
# 1000 simulatons, cosidering different contimination and then calculating 
# mean and sd of B1, B2 over thousand simulatons
####################################################

# lamda percent of Contimination
cont <- function(lambda, n){
  (lambda*n)/100
  
}

contimination <- c(cont(5,n), cont(10,n), cont(15,n), cont(25,n))


#---------------------- Huber

coeff_H <- matrix(NA,length(contimination),2)


res <- replicate(5,{
  
  x_1 <- rnorm(n)
  
  for(j in 1:length(contimination)){
    e <- c(rnorm(n-contimination[j]),rt(contimination[j], df=1))
    
    y = B0 + B1*x_1 + e
    wts <- rep(1,n)
    
    tune_range <- seq(0.2,3,.1)
    tau.h <- vector()
    for (indx in 1:length(tune_range)){
      for (itr in 1:10){
        mod <- lm(y~x_1,weights=wts)
        err <- scale(mod$resid)
        wts <- H.wt(err,tune_range[indx])
      }
      
      tau.h[indx] <- tau.hat.huber(err,tune_range[indx])
      
    }
    
    #plot(tune_range,tau.h,type = "l")
    tau.hat.h <- max(tau.h)
    
    
    for (i in 1:10){
      mod <- lm(y~x_1,weights=wts)
      err <- scale(mod$resid)
      wts <- H.wt(err,tune_range[min(which(tau.hat.h==tau.h))])
      coeff_H[j,] <- mod$coefficients
    }
  }
  coeff_H
})

apply(res, MARGIN=c(1,2),mean)
apply(res, MARGIN=c(1,2),sd)


#------------------------------- Tukey

coeff_T <- matrix(NA, length(contimination),2)

res <- replicate(5,{
  
  x_1 <- rnorm(n)
  
  for(j in 1:length(contimination)){
    e <- c(rnorm(n-contimination[j]),rt(contimination[j], df=1))
    
    y = B0 + B1*x_1 + e
    wts <- rep(1,n)
    
    tune_range <- seq(0.2,3,.1)
    wts <- rep(1,n)
    tau.T <- vector()
    for (indx in 1:length(tune_range)){
      for (itr in 1:10){
        mod <- lm(y~x_1,weights=wts)
        err <- scale(mod$resid)
        wts <- T.wt(err,tune_range[indx])
      }
      
      tau.T[indx] <- tau.hat.tukey(err,tune_range[indx])
      
    }
    
    tau.hat.T <- max(tau.T)
    
    for (i in 1:10){
      mod <- lm(y~x_1,weights=wts)
      err <- scale(mod$resid)
      wts <- T.wt(err,tune_range[min(which(tau.hat.T==tau.T))])
      coeff_T[j,] <- mod$coefficients
    }
  }
  coeff_T
})

apply(res, MARGIN=c(1,2), mean)
apply(res, MARGIN=c(1,2), sd)


#------------------------------- Andrew

coeff_Ad <- matrix(NA, length(contimination), 2)
res <- replicate(5,{
  
  x_1 <- rnorm(n)
  
  for(j in 1:length(contimination)){
    e <- c(rnorm(n-contimination[j]),rt(contimination[j], df=1))
    
    y = B0 + B1*x_1 + e
    wts <- rep(1,n)
    
    tune_range <- seq(0.2,3,.1)
    
    
    wts <- rep(1,n)
    tau.A <- vector()
    for (indx in 1:length(tune_range)){
      for (itr in 1:10){
        mod <- lm(y~x_1,weights=wts)
        err <- scale(mod$resid)
        wts <- A.wt(err,tune_range[indx])
      }
      
      tau.A[indx] <- tau.hat.Andrew(err,tune_range[indx])
      
    }
    
    tau.hat.A <- max(tau.A)
    for (i in 1:10){
      mod <- lm(y~x_1,weights=wts)
      err <- scale(mod$resid)
      wts <- A.wt(err,tune_range[min(which(tau.hat.A==tau.A))])
      #abline(mod)
      coeff_Ad[j,] <- mod$coefficients
    }
  }
  coeff_Ad
})
apply(res, MARGIN=c(1,2), mean)
apply(res, MARGIN=c(1,2), sd)


#------------------------------ OLS


coeff_OLS <- matrix(NA, length(contimination), 2)

res <- replicate(5,{
  
  x_1 <- rnorm(n)
  
  for(j in 1:length(contimination)){
    e <- c(rnorm(n-contimination[j]),rt(contimination[j], df=1))
    
    y = B0 + B1*x_1 + e
    
    mod <- lm(y~x_1)
    coeff_OLS[j,] <- mod$coefficients
  }
  coeff_OLS
})

apply(res, MARGIN=c(1,2), mean)
apply(res, MARGIN=c(1,2), sd)





#################################################
#   ploting B1 considering c with high effiency factor and other c values 
#   between 0.2 & 10 with difference 0.1
#################################################

n <- 1000

x_1 <- rnorm(n)
e <- rt(n, df=1)

y = B0 + B1*x_1 + e
wts <- rep(1,n)

tune_range <- seq(0.2,10,.1)

#--------------- For Huber -----------

coeff <- matrix(NA, length(tune_range),2)
tau.h <- vector()
for (indx in 1:length(tune_range)){
  for (itr in 1:10){
    mod <- lm(y~x_1,weights=wts)
    err <- scale(mod$resid)
    wts <- H.wt(err,tune_range[indx])
    coeff[indx,] <- mod$coefficients
  }
  
  tau.h[indx] <- tau.hat.huber(err,tune_range[indx])
  
}

plot(tune_range,tau.h,type = "l")
tau.hat.h <- max(tau.h)

plot(tune_range, coeff[,1], type = "l")
abline(v=tune_range[which(tau.hat.h==tau.h)])
abline(h=1)


#--------------- for Tukey ----------
coeff <- matrix(NA, length(tune_range),2)
tau.T <- vector()
for (indx in 1:length(tune_range)){
  for (itr in 1:10){
    mod <- lm(y~x_1,weights=wts)
    err <- scale(mod$resid)
    wts <- T.wt(err,tune_range[indx])
    coeff[indx,] <- mod$coefficients
  }
  
  tau.T[indx] <- tau.hat.tukey(err,tune_range[indx])
  
}

#plot(tune_range,tau.T,type = "l")
tau.hat.T <- max(tau.T)

plot(tune_range, coeff[,1], type = "l", ylim = c(0.5,3))
abline(v=tune_range[which(tau.hat.T==tau.T)])
abline(h=1)


#-------------- for Andrew's ----------
coeff <- matrix(NA, length(tune_range),2)
tau.A <- vector()
for (indx in 1:length(tune_range)){
  for (itr in 1:10){
    mod <- lm(y~x_1,weights=wts)
    err <- scale(mod$resid)
    wts <- A.wt(err,tune_range[indx])
    coeff[indx,] <- mod$coefficients
  }
  
  tau.A[indx] <- tau.hat.Andrew(err,tune_range[indx])
  
}

tau.hat.A <- max(tau.A)

plot(tune_range, coeff[,1], type = "l", ylim = c(0.5,3))
abline(v=tune_range[which(tau.hat.A==tau.A)])
abline(h=1)
