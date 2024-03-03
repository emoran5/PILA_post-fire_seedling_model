#Growth model
#G6 - Intercept, Site, Year/age, Ht, Clim  

load("PILASdl_Geno_Data.RData")

library(mvtnorm)


###functions
tnorm <- function(n,lo,hi,mu,sig){ 
z <- runif(n,0,1)
qlo <- pnorm((lo-mu)/sig)
qhi <- pnorm((hi-mu)/sig)
rr <- mu + sig*qnorm(qlo+z*(qhi-qlo))
return(rr)
}

b.update <- function(LM,X1,sinv,vinvert,bprior){  ##Update betas
    sx <- crossprod(X1)*sinv   
    sy <- crossprod(X1,LM)*sinv 
  bigv <- solve(sx+vinvert)  
  smallv <- sy+vinvert%*%bprior  
  b <- t(rmvnorm(1,bigv%*%smallv,bigv))
  return(b)
}

v.update <- function(Lm,X1,beta,N.surv){  #update variance
   sx <- crossprod((Lm-X1%*%beta))
    u1 <- S1 + 0.5*N.surv +1
    u2 <- S2 + 0.5*sx
  return(1/rgamma(1,u1,u2))
}


###########################Setup
##### Y.PILA.grow.all - vector of observations; N.grow is length

####### Setup design matrix
#X1 <- cbind(rep(1,N.grow),X.Sites.G,X.Yr.G,X.Ht.G,X.Sno.G) 
#X1 <- cbind(rep(1,N.grow),X.Sites.G,X.Yr.G,X.Ht.G,X.P.G) 
#X1 <- cbind(rep(1,N.grow),X.Sites.G,X.Yr.G,X.Ht.G,X.JMax.G) 
#X1 <- cbind(rep(1,N.grow),X.Sites.G,X.Yr.G,X.Ht.G,X.JMin.G) 
X1 <- cbind(rep(1,N.grow),X.Sites.G,X.Yr.G,X.Ht.G,X.CWD.G) 

I <- ncol(X1)  #columns of design matrix

###### Priors
Bpm <- c(5,rep(0,3),rep(0,4),0.15,0)
#Intercept = 5 cm annual growth; 0 site-level and year means; 
Bps <- c(6,rep(6,3),rep(6,4),6,6)   #variances
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
Beta <- c(runif(1,2,10),runif(3,-2,2),runif(4,-2,2),runif(1,-0.1,0.5),runif(1,-1,1))
Sig <- runif(1,0.05,0.5) 

ncyc <- 25000  #Number of Gibbs Steps
bgibbs <- matrix(NA,I,ncyc+1)
sgibbs <- rep(NA,ncyc+1)

bgibbs[,1] <- Beta
sgibbs[1] <- Sig

thin <- seq(12000,(ncyc+1),by=40)  #Burn-in and thinning

for (g in 1:ncyc){
	#update beta and sigma
	sinv <- 1/Sig
	Beta <- b.update(Y.PILA.grow.all,X1,sinv,vinvert,bprior)  
	Sig <- v.update(Y.PILA.grow.all,X1,Beta,N.surv)
	bgibbs[,g+1] <- Beta
	sgibbs[g+1] <- Sig

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
rownames(outpar) <- c('Int','site17','site18','site19','Yr1','Yr2','Yr3','Yr4','Ht','Clim','sig')
outpar


######### Predictive Loss Calculation
ng <- length(thin)
ypred <- matrix(0,ng,N.grow)
for (j in 1:ng){
	bg <- beta.thin[,j]; sd <- sqrt(sig.thin[j])

	q <- X1%*%bg
	for(i in 1:N.grow){
	ypred[j,i] <- tnorm(1,0,10,q[i],sd)		
	}
}

ym <- apply(ypred,2,mean)
yv <- apply(ypred,2,var)

gm <- sum((ym-Y.PILA.grow.all)^2)
pm <- sum(yv)
Dm <- gm +pm
Dm
