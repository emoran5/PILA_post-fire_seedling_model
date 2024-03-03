#Survival model
#S6 - Intercept, Site, Year/age, Ht, Clim  

load("PILASdl_Geno_Data.RData")

library(mvtnorm)


###functions
b.update <- function(LM,X1,sinv,vinvert,bprior){  ##Update betas
    sx <- crossprod(X1)*sinv   #Problem here because some heights are NA
    sy <- crossprod(X1,LM)*sinv #There are some NAs here too bc of Xs
  bigv <- solve(sx+vinvert)  #inverse of (n/sig + 1/priorvar) ###all NAs
  smallv <- sy+vinvert%*%bprior  #(n*meany/sig + priormean/priorvar) ###all NAs
  b <- t(rmvnorm(1,bigv%*%smallv,bigv))
  return(b)
}

v.update <- function(Lm,X1,beta,N.surv){  #update variance
   sx <- crossprod((Lm-X1%*%beta))
    u1 <- S1 + 0.5*N.surv +1
    u2 <- S2 + 0.5*sx
  return(1/rgamma(1,u1,u2))
}

lm.calc <- function(X1,beta,sig,N.surv){
	q <- X1%*%beta
	e <- rnorm(N.surv,0,sqrt(sig))
	lm <- q+e
	return(lm)
}

l1func <- function (Lm,Y.PILA.surv.all,N.surv){   #Calculates likelihood of Y given current parameters
	theta <- inv.logit(Lm)
	l1 <- rep(1,N.surv)
	for (i in 1:N.surv){
		l1[i] <- (theta[i]^Y.PILA.surv.all[i])*((1-theta[i])^(1-Y.PILA.surv.all[i]))
	}
	return(l1)
}

logit <- function(x){
	a <- log(x/(1-x))
	return(a)
}
inv.logit <- function(x){
	a <- exp(x)/(1+exp(x))
}


###########################Setup
##### Y.PILA.surv.all - vector of survival (1) vs. death (0) observations; N.surv is length

####### Setup design matrix
#X1 <- cbind(rep(1,N.surv),X.Sites.S,X.Yr.S,X.Ht.S,X.Sno.S) 
#X1 <- cbind(rep(1,N.surv),X.Sites.S,X.Yr.S,X.Ht.S,X.P.S) 
#X1 <- cbind(rep(1,N.surv),X.Sites.S,X.Yr.S,X.Ht.S,X.JMax.S) 
#X1 <- cbind(rep(1,N.surv),X.Sites.S,X.Yr.S,X.Ht.S,X.JMin.S) 
X1 <- cbind(rep(1,N.surv),X.Sites.S,X.Yr.S,X.Ht.S,X.CWD.S) 

I <- ncol(X1)  #columns of design matrix

###### Priors
Bpm <- c(2.5,rep(0,3),rep(-0.3,3),0.1,0)
#Intercept = 0.92% survival; with years becomes 65%, 85%, 3 x 90%; with height, 10 cm +intercept+ 1st yr is 66.8%, 20 cm is 69%
#Bps <- c(3,rep(3,3),rep(3,3),3,3)   #variances
Bps <- c(6,rep(6,3),rep(6,3),6,6)   #variances
bprior <- matrix(Bpm,I,1)  #beta mean matrix
vinvert <- solve(diag(Bps))  # 1/beta variances

 ###IG prior for error parameter
#S1 <- 2.01 #must be >2.  Find appropriate as follows...
  IGM <- (0.3)^2  #desired variance mean based on sd
#S2 <- (S1-1)*IGM
# IGV <- (S2^2)/(((S1-1)^2)*(S1-2)); IGV  #check if has high enough variance

S1 <- 2.01
S2 <- 0.35


###Starting values
Beta <- c(runif(1,1,5),runif(3,-0.1,0.1),runif(3,-0.6,0.2),runif(1,-0.1,0.2),runif(1,-0.2,0.2))
Sig <- runif(1,0.05,0.5)  

LM <- lm.calc(X1,Beta,Sig,N.surv)

ncyc <- 25000  #Number of Gibbs Steps
bgibbs <- matrix(NA,I,ncyc+1)
sgibbs <- rep(NA,ncyc+1)

bgibbs[,1] <- Beta
sgibbs[1] <- Sig

thin <- seq(12000,(ncyc+1),by=40)  #Burn-in and thinning

for (g in 1:ncyc){
	#update beta and sigma
	sinv <- 1/Sig
	Beta <- b.update(LM,X1,sinv,vinvert,bprior)  
	Sig <- v.update(LM,X1,Beta,N.surv)
	bgibbs[,g+1] <- Beta
	sgibbs[g+1] <- Sig

	Lm.new <- lm.calc(X1,Beta,Sig,N.surv)

	l1.now <- l1func(LM,Y.PILA.surv.all,N.surv)
	l1.new <- l1func(Lm.new,Y.PILA.surv.all,N.surv)

	pnow <- log(l1.now); pnew <- log(l1.new)
	r <- exp(pnew-pnow)

	z <- runif(N.surv,0,1)

	for(i in 1:N.surv){
		if (r[i]>z[i]) LM[i] <- Lm.new[i]
	}
	print(g)
}

par(mfrow=c(4,2))
  plot(seq(1,(ncyc+1)),bgibbs[1,],type='l')
   abline(h=Bpm[1],col='red')
  plot(seq(1,(ncyc+1)),bgibbs[2,],type='l')
   abline(h=Bpm[2],col='red')    
  plot(seq(1,(ncyc+1)),bgibbs[4,],type='l')  
   abline(h=Bpm[4],col='red')
  plot(seq(1,(ncyc+1)),bgibbs[6,],type='l')
   abline(h=Bpm[6],col='red')  
  plot(seq(1,(ncyc+1)),bgibbs[7,],type='l')
   abline(h=Bpm[7],col='red')   
  plot(seq(1,(ncyc+1)),bgibbs[10,],type='l')
   abline(h=Bpm[10],col='red')  
  plot(seq(1,(ncyc+1)),bgibbs[11,],type='l')
   abline(h=Bpm[11],col='red')   
  plot(seq(1,(ncyc+1)),sgibbs,type='l')
   abline(h=IGM,col='red')

  
beta.thin <- bgibbs[,thin]
sig.thin <- sgibbs[thin]

#### Posterior parameter output
outpar <- matrix(NA,I+1,5)
outpar[,1]<- c(Bpm,IGM)
outpar[,2]<- c(apply(beta.thin,1,mean),mean(sig.thin))
outpar[,3]<- c(apply(beta.thin,1,sd),sd(sig.thin))
outpar[,4]<- c(apply(beta.thin,1, quantile, probs= 0.025),quantile(sig.thin,0.025))
outpar[,5]<- c(apply(beta.thin,1, quantile, probs= 0.975),quantile(sig.thin,0.975))
colnames(outpar) <- c('prior','estimate','se','.025','.975')
rownames(outpar) <- c('Int','site17','site18','site19','Yr3','Yr4','Yr5','Ht','Clim','sig')
outpar


######### Predictive Loss Calculation
ng <- length(thin)
ypred <- matrix(0,ng,N.surv)
for (j in 1:ng){
	bg <- beta.thin[,j]; sd <- sqrt(sig.thin[j])

	q <- X1%*%bg
	e <- rnorm(N.surv,0,sd)
	lam <- q+e
	theta <- inv.logit(lam)
	for(i in 1:N.surv){
	ypred[j,i] <- rbinom(1,1,theta[i])
	}
}

ym <- apply(ypred,2,mean)
yv <- apply(ypred,2,var)

gm <- sum((ym-Y.PILA.surv.all)^2)
pm <- sum(yv)
Dm <- gm +pm
Dm
