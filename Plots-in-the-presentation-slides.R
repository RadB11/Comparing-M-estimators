# R codes for Plots in the presentation slides

######### Huber psi function #################
H.psi <- function(error, tune) {
  error*(abs(error) <= tune) + sign(error)*tune*(abs(error) > tune)  
}

######### Huber weight function #################

H.wt <- function(error, tune) {
  1*(abs(error) <= tune)  
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


#--------------- Andrew's function -------------

A.psi <- function(error, tune) {
  tune*sin(error/tune)*(abs(error) <= tune*pi)
}


A.wt <- function(error, tune) {
  tune*(sin(error/tune)*(abs(error) <= tune*pi))/error
}

############################################################



pdf("Hpsi.pdf")
par(mar=c(4, 5, 2, 1))
x <-sort(rnorm(10000))
Psi <- H.psi(x,2)
plot(x,Psi, type = "l", xlab = expression(epsilon), ylim=c(-3,3),xlim=c(-3,3),
     ylab = expression(Psi(epsilon)),las=1, lwd=2, cex.lab=2, cex.axis=2)
abline(h=0,lwd=2)
abline(0,1,lwd=2)
dev.off()

pdf("Hwt.pdf")
par(mar=c(4, 5, 2, 1))
w <- H.wt(x,2)
plot(x,w, type = "l", xlab = expression(epsilon), ylim=c(-0,1), xlim=c(-3,3),
     ylab = expression(w(epsilon)),las=1, lwd=2, cex.lab=2, cex.axis=2)
abline(h=1,lwd=2)
dev.off()

pdf("Tpsi.pdf")
par(mar=c(4, 5, 2, 1))
Psi <- T.psi(x,2)
plot(x,Psi, type = "l", xlab = expression(epsilon), xlim=c(-3,3),
     ylab = expression(Psi(epsilon)),las=1, lwd=2, cex.lab=2, cex.axis=2)
abline(h=0,lwd=2)
abline(0,1,lwd=2)
dev.off()

pdf("Twt.pdf")
par(mar=c(4, 5, 2, 1))
w <- T.wt(x,2)
plot(x,w, type = "l", xlab = expression(epsilon), ylim=c(-0,1), xlim=c(-3,3),
     ylab = expression(w(epsilon)),las=1, lwd=2, cex.lab=2, cex.axis=2)
abline(h=1,lwd=2)
dev.off()

pdf("Apsi.pdf")
par(mar=c(4, 5, 2, 1))
Psi <- A.psi(x,.9)
plot(x,Psi, type = "l", xlab = expression(epsilon), xlim=c(-3,3), 
     ylab = expression(Psi(epsilon)),las=1, lwd=2, cex.lab=2, cex.axis=2)
abline(h=0,lwd=2)
abline(0,1,lwd=2)
dev.off()

pdf("Awt.pdf")
par(mar=c(4, 5, 2, 1))
w <- A.wt(x,.9)
plot(x,w, type = "l", xlab = expression(epsilon), ylim=c(-0,1), xlim=c(-3,3),
     ylab = expression(w(epsilon)),las=1, lwd=2, cex.lab=2, cex.axis=2)
abline(h=1,lwd=2)
dev.off()




# Comparing all methods' efficiency

n <- 1000


B_0 <- 1
B_1 <- 1

res <- replicate(100,{
  
  X_1 <- rnorm(n)
  #
  
  
  for(i in 1:length(part_of_n)){
    
    e <- c(rnorm(n - part_of_n[i]), rnorm(part_of_n[i]))
    y <- B_0 + B_1 * X_1 + e
    coeff[i,] <- lm(y~X_1)$coefficients
  }
  coeff
})

coeff <- as.vector(res) 
arr <- as.vector(array(c(1:8),dim = c(length(part_of_n),2,100)))

boxplot(coeff~arr, outline=FALSE, las=2, xaxt="n",at=c(1,2,3,4,5,6,7,8),
        cex.axis=2, cex.lab=2, col=c("red","green"))