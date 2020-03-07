# We are comapring the M-estimators based on the effciency factor

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
tau.hat.huber <- function(x, par) {
  n <- length(x)
  tau <- -sum(abs(x) <= par)^2/(n*sum((abs(x) <= par)*x^2 + par^2 *(abs(x) > par)))
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
tau.hat.tukey <- function(x, par) {
  tau <- -mean(T.deri(x, par))^2/(mean((abs(x) <= par)*T.psi(x,par)^2))
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
tau.hat.Andrew <- function(x, par) {
  tau <- -mean(A.deri(x, par))^2/(mean((abs(x) <= par*pi)*A.psi(x,par)^2))
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
tau.hat.Alam <- function(x, par) {
  tau <- -mean(Alam.deri(x, par))^2/(mean((abs(x) <= par)*Alam.psi(x,par)^2))
}


###################### Insha's function
Insha.psi <- function(error, tune) {
  (error/(1+(error/tune)^4)^2)
}

# weight function
Insha.wt <- function(error, tune) {
  (1/(1+(error/tune)^4)^2)
}

# 2nd derivative 
Insha.deri <- function(error, tune){
  ((((1+(error/tune)^2) - 8*(error/tune)^4))/((1+(error/tune)^4)^3))
}
# 2nd derivative 
# Insha.deri <- function(error, tune){
#   (tune^12-7*tune^8*error^4)/(tune^4+error^4)^3
# }
##################################

tau.hat.Insha <- function(x, par) {
  tau <- -mean(Insha.deri(x, par))^2/mean(Insha.psi(x,par)^2)
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
  ((2/3)*((1-(error/tune)^4)^2 - 8*(error/tune)^4*(1-(error/tune)^4)))*(abs(error) <= tune)
}

#####################################

tau.hat.Asad <- function(x, par) {
  tau <- -mean(Asad.deri(x, par))^2/(mean((abs(x) <= par)*Asad.psi(x, par)^2))
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

tau.hat.beta.Qader <- function(x, par) {
  tau <- -mean(Qader.beta.deri(x, par))^2/(mean((abs(x) <= par)*Qader.beta.psi(x, par)^2))
}



B0 <- 1
B1 <- 1


n = 1000


res <- replicate(500,{

  x_1 <- rnorm(n)
  
  ########## error distributions
  e <- rt(n, df=1)
  
  #e <- rnorm(n)
  
  #e <- rchisq(n, df = 4)
  #e <- rcauchy(n, 0,100)
  
  
  
  y = B0 + B1*x_1 + e
  

  #------------- Huber method --------------
  
  wts <- rep(1,n)
  for (i in 1:20){
   mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.h <- (optim(par =  3, fn = tau.hat.huber, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- H.wt(err,tau.hat.h)
    coeff_H <- mod$coefficients
  }
  
  # ---------------- Tukey method -------------
  
  wts <- rep(1,n)
  for (i in 1:20){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.T <- (optim(par =  4.685, fn = tau.hat.tukey, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- T.wt(err,tau.hat.T)
    coeff_T <- mod$coefficients
  }
  
  #--------------------- Alamgir mwthod-------------
  
  wts <- rep(1,n)
  for (i in 1:20){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.Alamv <- (optim(par =  3, fn = tau.hat.Alam, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- Alam.wt(err,tau.hat.Alamv)
    coeff_Alam <- mod$coefficients
  }
  
  #------------------------ Insha method
  
  wts <- rep(1,n)
  for (i in 1:20){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.I <- (optim(par =  3, fn = tau.hat.Insha, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- Insha.wt(err,tau.hat.I)
    coeff_Insha <- mod$coefficients
  }
  
  #---------------------------- Asad method
  
  wts <- rep(1,n)
  for (i in 1:20){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.As <- (optim(par =  3, fn = tau.hat.Asad, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- Asad.wt(err,tau.hat.As)
    coeff_Asad <- mod$coefficients
  }
  
  #---------------------------- Qader beta function method
  
  wts <- rep(1,n)
  for (i in 1:20){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.Q <- (optim(par =  3, fn = tau.hat.beta.Qader, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- Qader.beta.wt(err,tau.hat.Q)
    coeff_beta_Qader <- mod$coefficients
  }
  
  #-----------------OLS method----------- It can also be included for comparison
  
  mod <- lm(y~x_1)
  coeff_OLS <- mod$coefficients
  
  #----------- Andrew's method ---------
  
  wts <- rep(1,n)
  for (i in 1:20){
    mod <- lm(y~x_1,weights=wts)
    SS <- mad(mod$resid)
    err <- mod$resid/sqrt(SS)
    tau.hat.A <- (optim(par =  3, fn = tau.hat.Andrew, x = err,lower = .3, upper=20,method = "Brent"))$par
    wts <- A.wt(err,tau.hat.A)
    #abline(mod)
    coeff_Ad <- mod$coefficients
  }
  
  c(coeff_T[1],coeff_Ad[1],coeff_Alam[1],coeff_Insha[1], coeff_Asad[1],coeff_beta_Qader[1],
    coeff_T[2],coeff_Ad[2],coeff_Alam[2], coeff_Insha[2], coeff_Asad[2],coeff_beta_Qader[2])
})


pdf(file="boxplot.pdf")
par (mar=c ( 6 , 7 , 2 , 1 ) , mgp = c ( 4 , 1 , 0 ) )
boxplot(res~row(res), outline=F, las=2, xaxt="n",at=c(1,2,3,4,5,6,8,9,10,11,12,13),
        cex.axis=1.3, cex.lab=1.5, ylim = c(0.5,1.5),
        main = "M-estimators comparison based on efficiency factor" , xlab = "" ,
        ylab = expression(beta),
        col=c("green", "blue", "red","orange","green3", "coral1"))
abline(h=1)
#points(apply(res,1,mean),pch=16)
mtext(expression(hat(beta)[0]), at=4, line = -29.5, cex = 1.5)
mtext(expression(hat(beta)[1]), at=12, line = -29.5, cex = 1.5)


legend("topright", legend = c("Tukey", "Andrew","Alamgir", "Insha", "Asad", "Qader-beta"),
       cex = .7,box.lwd = 1,box.lty = 1,
       text.font = 2,col = c("green", "blue", "red", "orange", "green3","coral1"),
       pch = 15,horiz = F)
dev.off()


