##Basing this on TagSdlDataOrg5 in **

bframe1 <-data.frame(read.table("Peavine_PILA_dec2_trim.txt",header=T,fill=T,row.names=NULL))
bframe2 <-data.frame(read.table("Stumpy2018_PILAdec2_trim.txt",header=T,fill=T,row.names=NULL)) 
bframe3 <-data.frame(read.table("Stumpy2017_PILAdec2_trim.txt",header=T,fill=T,row.names=NULL))

CWD <- data.frame(read.table("CWD_91_21.txt",header=T,fill=T,row.names=NULL))
Snow <- data.frame(read.table("Sno_91_21.txt",header=T,fill=T,row.names=NULL))
Precip <- data.frame(read.table("Precip_91_21.txt",header=T,fill=T,row.names=NULL))
JMin <- data.frame(read.table("JMin_91_21.txt",header=T,fill=T,row.names=NULL))
JMax <- data.frame(read.table("JMax_91_21.txt",header=T,fill=T,row.names=NULL))

#### PILA seedlings 

PILA.19.Sdl <- cbind(bframe1[,4],bframe1[,7],bframe1[,9],bframe1[,3],bframe1[,12:16],bframe1[,18],bframe1[,23])
PILA.18.Sdl <- cbind(bframe2[,3],bframe2[,6],bframe2[,8],bframe2[,2],bframe2[,12:15],bframe2[,18],bframe2[,21],bframe2[,25])
PILA.17.Sdl <- cbind(bframe3[,3:4],bframe3[,6],bframe3[,2],bframe3[,11:13],bframe3[,16],bframe3[,19],bframe3[,22],bframe3[,26])

colnames(PILA.19.Sdl)<-c('Species','SZ','ElevB','Lot_Name','PlantYr','MortYr','Ht17','Ht18','Ht19','Ht20','Ht21')
colnames(PILA.18.Sdl)<-c('Species','SZ','ElevB','Lot_Name','PlantYr','MortYr','Ht17','Ht18','Ht19','Ht20','Ht21')
colnames(PILA.17.Sdl)<-c('Species','SZ','ElevB','Lot_Name','PlantYr','MortYr','Ht17','Ht18','Ht19','Ht20','Ht21')

#seedling numbers
N.PILA.19 <- nrow(PILA.19.Sdl)
N.PILA.18 <- nrow(PILA.18.Sdl)
N.PILA.17 <- nrow(PILA.17.Sdl)

PILA.19.survmat <- matrix(1,N.PILA.19,3);PILA.19.Htmat <- matrix(NA,N.PILA.19,3)
PILA.18.survmat <- matrix(1,N.PILA.18,4);PILA.18.Htmat <- matrix(NA,N.PILA.18,4)
PILA.17.survmat <- matrix(1,N.PILA.17,5);PILA.17.Htmat <- matrix(NA,N.PILA.17,5)

mn.Ht.19 <- 0; tot.m.19 <- 0
for(i in 1:N.PILA.19){
	if(PILA.19.Sdl[i,9]>0) {
		mn.Ht.19 <- mn.Ht.19 + PILA.19.Sdl[i,9]
		tot.m.19 <- tot.m.19 + 1
	}
}
I.mn.Ht.19 <- mn.Ht.19/tot.m.19 #15.04

mn.Ht.18 <- 0; tot.m.18 <- 0
for(i in 1:N.PILA.18){
	if(PILA.18.Sdl[i,8]>0) {
		mn.Ht.18 <- mn.Ht.18 + PILA.18.Sdl[i,8]
		tot.m.18 <- tot.m.18 + 1
	}
}
I.mn.Ht.18 <- mn.Ht.18/tot.m.18 #12.89

mn.Ht.17 <- 0; tot.m.17 <- 0
for(i in 1:N.PILA.17){
	if(PILA.17.Sdl[i,7]>0) {
		mn.Ht.17 <- mn.Ht.17 + PILA.17.Sdl[i,7]
		tot.m.17 <- tot.m.17 + 1
	}
}
I.mn.Ht.17 <- mn.Ht.17/tot.m.17 #14.57


for(i in 1:N.PILA.19){
	if(PILA.19.Sdl[i,'MortYr']>0){
		if(PILA.19.Sdl[i,'MortYr']==2019) PILA.19.survmat[i,]<-0
		if(PILA.19.Sdl[i,'MortYr']==2020) PILA.19.survmat[i,c(2,3)]<-0
		if(PILA.19.Sdl[i,'MortYr']==2021) PILA.19.survmat[i,3]<-0
	}
	#Year 1 ht
	if(PILA.19.Sdl[i,'Ht19']>0) PILA.19.Htmat[i,1] <- PILA.19.Sdl[i,'Ht19']
		if(PILA.19.Sdl[i,'Ht19']<0) PILA.19.Htmat[i,1] <- 15.04 #average initial height
	#Year 2 ht
	if(PILA.19.Sdl[i,'Ht20']>0) PILA.19.Htmat[i,2] <- PILA.19.Sdl[i,'Ht20']
		if(PILA.19.Sdl[i,'Ht20']<0){
			if(PILA.19.Sdl[i,'Ht21']>0) PILA.19.Htmat[i,2] <- sum(PILA.19.Htmat[i,1],PILA.19.Sdl[i,'Ht21'])/2
			if(PILA.19.Sdl[i,'Ht21']<0) PILA.19.Htmat[i,2] <- PILA.19.Htmat[i,1]
		}
	#Year 3 ht
	if(PILA.19.Sdl[i,'Ht21']>0) PILA.19.Htmat[i,3] <- PILA.19.Sdl[i,'Ht21']
		if(PILA.19.Sdl[i,'Ht21']<0 & PILA.19.survmat[i,3]==1) PILA.19.Htmat[i,3] <- 
			PILA.19.Htmat[i,2]
}

for(i in 1:N.PILA.18){    
	if(PILA.18.Sdl[i,'MortYr']>0){
		if(PILA.18.Sdl[i,'MortYr']==2018) PILA.18.survmat[i,]<-0
		if(PILA.18.Sdl[i,'MortYr']==2019) PILA.18.survmat[i,c(2,3,4)]<-0
		if(PILA.18.Sdl[i,'MortYr']==2020) PILA.18.survmat[i,c(3,4)]<-0
		if(PILA.18.Sdl[i,'MortYr']==2021) PILA.18.survmat[i,4]<-0
	}
	#Year 1 ht
	if(PILA.18.Sdl[i,'Ht18']>0) PILA.18.Htmat[i,1] <- PILA.18.Sdl[i,'Ht18']
		if(PILA.18.Sdl[i,'Ht18']<0) PILA.18.Htmat[i,1] <- 12.89 #average initial height
	#Year 2 ht
	if(PILA.18.Sdl[i,'Ht19']>0) PILA.18.Htmat[i,2] <- PILA.18.Sdl[i,'Ht19']
		if(PILA.18.Sdl[i,'Ht19']<0){
			if(PILA.18.Sdl[i,'Ht20']>0) PILA.18.Htmat[i,2] <- sum(PILA.18.Htmat[i,1],PILA.18.Sdl[i,'Ht20'])/2
			if(PILA.18.Sdl[i,'Ht20']<0 & PILA.18.Sdl[i,'Ht21']>0) PILA.18.Htmat[i,2] <- 
				(sum(PILA.18.Htmat[i,1],PILA.18.Sdl[i,'Ht21'])/3)
			if(PILA.18.Sdl[i,'Ht20']<0 & PILA.18.Sdl[i,'Ht21']<0 & PILA.18.survmat[i,2]==1) PILA.18.Htmat[i,2] <- 
				PILA.18.Htmat[i,1]
		}
	#Year 3 ht
	if(PILA.18.Sdl[i,'Ht20']>0) PILA.18.Htmat[i,3] <- PILA.18.Sdl[i,'Ht20']
		if(PILA.18.Sdl[i,'Ht20']<0){
			if(PILA.18.Sdl[i,'Ht21']>0) PILA.18.Htmat[i,3] <- sum(PILA.18.Htmat[i,2],PILA.18.Sdl[i,'Ht21'])/2
			if(PILA.18.Sdl[i,'Ht21']<0 & PILA.18.survmat[i,3]==1) PILA.18.Htmat[i,3] <- PILA.18.Htmat[i,2]
		}
	#Year 4 ht
	if(PILA.18.Sdl[i,'Ht21']>0) PILA.18.Htmat[i,4] <- PILA.18.Sdl[i,'Ht21']
		if(PILA.18.Sdl[i,'Ht21']<0 & PILA.18.survmat[i,4]==1) PILA.18.Htmat[i,4] <- PILA.18.Htmat[i,3]
}


for(i in 1:N.PILA.17){
	if(PILA.17.Sdl[i,'MortYr']>0){
		if(PILA.17.Sdl[i,'MortYr']==2017) PILA.17.survmat[i,]<-0
		if(PILA.17.Sdl[i,'MortYr']==2018) PILA.17.survmat[i,c(2,3,4,5)]<-0
		if(PILA.17.Sdl[i,'MortYr']==2019) PILA.17.survmat[i,c(3,4,5)]<-0
		if(PILA.17.Sdl[i,'MortYr']==2020) PILA.17.survmat[i,c(4,5)]<-0
		if(PILA.17.Sdl[i,'MortYr']==2021) PILA.17.survmat[i,5]<-0
	}
	#Year 1 ht
	if(PILA.17.Sdl[i,'Ht17']>0) PILA.17.Htmat[i,1] <- PILA.17.Sdl[i,'Ht17']
		if(PILA.17.Sdl[i,'Ht17']<0) PILA.17.Htmat[i,1] <- 14.57 #average initial height
	#Year 2 ht	
	if(PILA.17.Sdl[i,'Ht18']>0) PILA.17.Htmat[i,2] <- PILA.17.Sdl[i,'Ht18']
		if(PILA.17.Sdl[i,'Ht18']<0){
			if(PILA.17.Sdl[i,'Ht19']>0) PILA.17.Htmat[i,2] <- sum(PILA.17.Htmat[i,1],PILA.17.Sdl[i,'Ht19'])/2
		  	if(PILA.17.Sdl[i,'Ht19']<0 & PILA.17.Sdl[i,'Ht20']>0) PILA.17.Htmat[i,2] <- 
		  		(sum(PILA.17.Htmat[i,1],PILA.17.Sdl[i,'Ht20'])/3)
			if(PILA.17.Sdl[i,'Ht19']<0 & PILA.17.Sdl[i,'Ht20']<0 & PILA.17.survmat[i,2]==1) PILA.17.Htmat[i,2] <- 
				PILA.17.Htmat[i,1]		  	
		}
	#Year 3 ht
	if(PILA.17.Sdl[i,'Ht19']>0) PILA.17.Htmat[i,3] <- PILA.17.Sdl[i,'Ht19']
		if(PILA.17.Sdl[i,'Ht19']<0) {
			if(PILA.17.Sdl[i,'Ht20']>0) PILA.17.Htmat[i,3] <- sum(PILA.17.Htmat[i,2],PILA.17.Sdl[i,'Ht20'])/2
		   	if(PILA.17.Sdl[i,'Ht20']<0 & PILA.17.Sdl[i,'Ht21']>0) PILA.17.Htmat[i,3] <- 
		   		(sum(PILA.17.Htmat[i,1],PILA.17.Sdl[i,'Ht21'])/3)
			if(PILA.17.Sdl[i,'Ht20']<0 & PILA.17.Sdl[i,'Ht21']<0 & PILA.17.survmat[i,3]==1) PILA.17.Htmat[i,3] <- 	
				PILA.17.Htmat[i,2]		  	
		}
    #Year 4 ht
	if(PILA.17.Sdl[i,'Ht20']>0) PILA.17.Htmat[i,4] <- PILA.17.Sdl[i,'Ht20']
		if(PILA.17.Sdl[i,'Ht20']<0){
			if(PILA.17.Sdl[i,'Ht21']>0) PILA.17.Htmat[i,4] <- sum(PILA.17.Htmat[i,3],PILA.17.Sdl[i,'Ht21'])/2
			if(PILA.17.Sdl[i,'Ht21']<0 & PILA.17.survmat[i,4]==1) PILA.17.Htmat[i,4] <- PILA.17.Htmat[i,3]
		}
	#Year 5 ht
	if(PILA.17.Sdl[i,'Ht21']>0) PILA.17.Htmat[i,5] <- PILA.17.Sdl[i,'Ht21']
		if(PILA.17.Sdl[i,'Ht21']<0 & PILA.17.survmat[i,5]==1) PILA.17.Htmat[i,5] <- PILA.17.Htmat[i,4]
}



#########Graphs of survival and growth
PILA.17.start <- nrow(PILA.17.survmat); PILA.18.start <- nrow(PILA.18.survmat)
PILA.19.start <- nrow(PILA.19.survmat)

PILA.17.total <- c(PILA.17.start,apply(PILA.17.survmat,2,sum))
PILA.17.rel_surv <- PILA.17.total/PILA.17.start

PILA.18.total <- c(PILA.18.start,apply(PILA.18.survmat,2,sum))
PILA.18.rel_surv <- PILA.18.total/PILA.18.start

PILA.19.total <- c(PILA.19.start,apply(PILA.19.survmat,2,sum))
PILA.19.rel_surv <- PILA.19.total/PILA.19.start

X1 <- c(0,0.5,1,2,3,4) #years from planting


Sdl.elev.17 <- c(3000,4000,4500,5000,5500,6000); Sdl.elev.18 <- c(4500,5000,5250,5500,6000); Sdl.elev.19 <- c(4500,5000,5500,6000,7500)

PILA.30.17 <- which(PILA.17.Sdl$ElevB == 3000); PILA.40.17 <- which(PILA.17.Sdl$ElevB == 4000)
PILA.45.17 <- which(PILA.17.Sdl$ElevB == 4500); PILA.50.17 <- which(PILA.17.Sdl$ElevB == 5000)
PILA.55.17 <- which(PILA.17.Sdl$ElevB == 5500); PILA.60.17 <- which(PILA.17.Sdl$ElevB == 6000)

PILA.30.17.total <-c(length(PILA.30.17),apply(PILA.17.survmat[PILA.30.17,],2,sum))
PILA.30.17.rel_surv <- PILA.30.17.total/length(PILA.30.17)
PILA.40.17.total <-c(length(PILA.40.17),apply(PILA.17.survmat[PILA.40.17,],2,sum))
PILA.40.17.rel_surv <- PILA.40.17.total/length(PILA.40.17)
PILA.45.17.total <-c(length(PILA.45.17),apply(PILA.17.survmat[PILA.45.17,],2,sum))
PILA.45.17.rel_surv <- PILA.45.17.total/length(PILA.45.17)
PILA.50.17.total <-c(length(PILA.50.17),apply(PILA.17.survmat[PILA.50.17,],2,sum))
PILA.50.17.rel_surv <- PILA.50.17.total/length(PILA.50.17)
PILA.55.17.total <-c(length(PILA.55.17),apply(PILA.17.survmat[PILA.55.17,],2,sum))
PILA.55.17.rel_surv <- PILA.55.17.total/length(PILA.55.17)
PILA.60.17.total <-c(length(PILA.60.17),apply(PILA.17.survmat[PILA.60.17,],2,sum))
PILA.60.17.rel_surv <- PILA.60.17.total/length(PILA.60.17)

PILA.45.18 <- which(PILA.18.Sdl$ElevB == 4500); PILA.50.18 <- which(PILA.18.Sdl$ElevB == 5000)
PILA.52.18 <- which(PILA.18.Sdl$ElevB == 5250)
PILA.55.18 <- which(PILA.18.Sdl$ElevB == 5500); PILA.60.18 <- which(PILA.18.Sdl$ElevB == 6000)

PILA.45.18.total <-c(length(PILA.45.18),apply(PILA.18.survmat[PILA.45.18,],2,sum))
PILA.45.18.rel_surv <- PILA.45.18.total/length(PILA.45.18)
PILA.50.18.total <-c(length(PILA.50.18),apply(PILA.18.survmat[PILA.50.18,],2,sum))
PILA.50.18.rel_surv <- PILA.50.18.total/length(PILA.50.18)
PILA.52.18.total <-c(length(PILA.52.18),apply(PILA.18.survmat[PILA.52.18,],2,sum))
PILA.52.18.rel_surv <- PILA.52.18.total/length(PILA.52.18)
PILA.55.18.total <-c(length(PILA.55.18),apply(PILA.18.survmat[PILA.55.18,],2,sum))
PILA.55.18.rel_surv <- PILA.55.18.total/length(PILA.55.18)
PILA.60.18.total <-c(length(PILA.60.18),apply(PILA.18.survmat[PILA.60.18,],2,sum))
PILA.60.18.rel_surv <- PILA.60.18.total/length(PILA.60.18)

PILA.45.19 <- which(PILA.19.Sdl$ElevB == 4500); PILA.50.19 <- which(PILA.19.Sdl$ElevB == 5000)
PILA.55.19 <- which(PILA.19.Sdl$ElevB == 5500); PILA.60.19 <- which(PILA.19.Sdl$ElevB == 6000)
PILA.75.19 <- which(PILA.19.Sdl$ElevB == 7500)

PILA.45.19.total <-c(length(PILA.45.19),apply(PILA.19.survmat[PILA.45.19,],2,sum))
PILA.45.19.rel_surv <- PILA.45.19.total/length(PILA.45.19)
PILA.50.19.total <-c(length(PILA.50.19),apply(PILA.19.survmat[PILA.50.19,],2,sum))
PILA.50.19.rel_surv <- PILA.50.19.total/length(PILA.50.19)
PILA.55.19.total <-c(length(PILA.55.19),apply(PILA.19.survmat[PILA.55.19,],2,sum))
PILA.55.19.rel_surv <- PILA.55.19.total/length(PILA.55.19)
PILA.60.19.total <-c(length(PILA.60.19),apply(PILA.19.survmat[PILA.60.19,],2,sum))
PILA.60.19.rel_surv <- PILA.60.19.total/length(PILA.60.19)
PILA.75.19.total <-c(length(PILA.75.19),apply(PILA.19.survmat[PILA.75.19,],2,sum))
PILA.75.19.rel_surv <- PILA.75.19.total/length(PILA.75.19)

jpeg(filename="PILA_Surv_all.jpg",width=10,height=10,units="in",res=500)

par(mfrow=c(2,2))

plot(X1,PILA.17.rel_surv,type='l',lty=1,xlab='Year Since Planting',main='Overall PILA survival',ylab='Relative survival',col='red',ylim=c(0,1),lwd=2)
points(X1[1:5],PILA.18.rel_surv,type='l',lty=2,col='blue',lwd=2)
points(X1[1:4],PILA.19.rel_surv,type='l',lty=3,col='purple',lwd=2)
legend(2.5,0.85,legend=c("2017 site","2018 site","2019 site"),col=c('red','blue','purple'),lty=c(1,2,3),lwd=c(2,2,2),cex=0.7)

plot(X1,PILA.17.rel_surv,type='l',lty=1,xlab='Year Since Planting',main='2017 PILA by elev',ylab='Relative survival',ylim=c(0,1),lwd=2)
points(X1,PILA.30.17.rel_surv,type='l',col='red',lwd=2,lty=2)
points(X1,PILA.40.17.rel_surv,type='l',col='orange',lwd=2,lty=2)
points(X1,PILA.45.17.rel_surv,type='l',col='green',lwd=2)
points(X1,PILA.50.17.rel_surv,type='l',col='light blue',lwd=2,lty=2)
points(X1,PILA.55.17.rel_surv,type='l',col='blue',lwd=2,lty=2)
points(X1,PILA.60.17.rel_surv,type='l',col='purple',lwd=2,lty=2)
legend(3,0.9,legend=c('All',Sdl.elev.17),col=c('black','red','orange','green','light blue','blue','purple'),lty=c(1,2,2,1,2,2,2),lwd=c(2,2,2,2,2,2,2),cex=0.7)

plot(X1[1:5],PILA.18.rel_surv,type='l',lty=1,xlab='Year Since Planting',main='2018 PILA by elev',ylab='Relative survival',ylim=c(0,1),xlim=c(0,4),lwd=2)
points(X1[1:5],PILA.45.18.rel_surv,type='l',col='green',lwd=2)
points(X1[1:5],PILA.50.18.rel_surv,type='l',col='light blue',lwd=2,lty=2)
points(X1[1:5],PILA.52.18.rel_surv,type='l',col='cornflowerblue',lwd=2,lty=2)
points(X1[1:5],PILA.55.18.rel_surv,type='l',col='blue',lwd=2,lty=2)
points(X1[1:5],PILA.60.18.rel_surv,type='l',col='purple',lwd=2,lty=2)
legend(3,0.9,legend=c('All',Sdl.elev.18),col=c('black','green','light blue','cornflowerblue','blue','purple'),lty=c(1,1,2,2,2,2),lwd=c(2,2,2,2,2,2),cex=0.7)

plot(X1[1:4],PILA.19.rel_surv,type='l',lty=1,xlab='Year Since Planting',main='2019 PILA by elev',ylab='Relative survival',ylim=c(0,1),lwd=2,xlim=c(0,4))
points(X1[1:4],PILA.45.19.rel_surv,type='l',col='green',lwd=2)
points(X1[1:4],PILA.50.19.rel_surv,type='l',col='light blue',lwd=2,lty=2)
points(X1[1:4],PILA.55.19.rel_surv,type='l',col='blue',lwd=2,lty=2)
points(X1[1:4],PILA.60.19.rel_surv,type='l',col='purple',lwd=2,lty=2)
points(X1[1:4],PILA.75.19.rel_surv,type='l',col='darkorchid4',lwd=2,lty=2)
legend(3,0.9,legend=c('All',Sdl.elev.19),col=c('black','green','light blue','blue','purple','darkorchid4'),lty=c(1,1,2,2,2,2),lwd=c(2,2,2,2,2,2),cex=0.7)

dev.off()

#Heights

PILA.17.Ht <- apply(PILA.17.Htmat,2,mean,na.rm=T);PILA.18.Ht <- apply(PILA.18.Htmat,2,mean,na.rm=T)
PILA.19.Ht <- apply(PILA.19.Htmat,2,mean,na.rm=T)

X2 <- c(0,1,2,3,4)

PILA.30.17.Ht <- apply(PILA.17.Htmat[PILA.30.17,],2,mean,na.rm=T)
PILA.40.17.Ht <- apply(PILA.17.Htmat[PILA.40.17,],2,mean,na.rm=T)
PILA.45.17.Ht <- apply(PILA.17.Htmat[PILA.45.17,],2,mean,na.rm=T)
PILA.50.17.Ht <- apply(PILA.17.Htmat[PILA.50.17,],2,mean,na.rm=T)
PILA.55.17.Ht <- apply(PILA.17.Htmat[PILA.55.17,],2,mean,na.rm=T)
PILA.60.17.Ht <- apply(PILA.17.Htmat[PILA.60.17,],2,mean,na.rm=T)

PILA.45.18.Ht <- apply(PILA.18.Htmat[PILA.45.18,],2,mean,na.rm=T)
PILA.50.18.Ht <- apply(PILA.18.Htmat[PILA.50.18,],2,mean,na.rm=T)
PILA.52.18.Ht <- apply(PILA.18.Htmat[PILA.52.18,],2,mean,na.rm=T)
PILA.55.18.Ht <- apply(PILA.18.Htmat[PILA.55.18,],2,mean,na.rm=T)
PILA.60.18.Ht <- apply(PILA.18.Htmat[PILA.60.18,],2,mean,na.rm=T)

PILA.45.19.Ht <- apply(PILA.19.Htmat[PILA.45.19,],2,mean,na.rm=T)
PILA.50.19.Ht <- apply(PILA.19.Htmat[PILA.50.19,],2,mean,na.rm=T)
PILA.55.19.Ht <- apply(PILA.19.Htmat[PILA.55.19,],2,mean,na.rm=T)
PILA.60.19.Ht <- apply(PILA.19.Htmat[PILA.60.19,],2,mean,na.rm=T)
PILA.75.19.Ht <- apply(PILA.19.Htmat[PILA.75.19,],2,mean,na.rm=T)


jpeg(filename="PILA_Ht_all.jpg",width=10,height=10,units="in",res=500)

par(mfrow=c(2,2))

plot(X2,PILA.17.Ht,type='l',lty=1,xlab='Year Since Planting',main='Overall PILA height',ylab='Height (cm)',col='red',ylim=c(0,60),lwd=2)
points(X2[1:4],PILA.18.Ht,type='l',lty=2,col='blue',lwd=2)
points(X2[1:3],PILA.19.Ht,type='l',lty=3,col='purple',lwd=2)
legend(2.8,30,legend=c("2017 site","2018 site","2019 site"),col=c('red','blue','purple'),lty=c(1,2,3),lwd=c(2,2,2),cex=0.7)

plot(X2,PILA.17.Ht,type='l',lty=1,xlab='Year Since Planting',main='2017 PILA by elev',ylab='Height (cm)',lwd=2,ylim=c(0,60))
points(X2,PILA.30.17.Ht,type='l',col='red',lwd=2,lty=2)
points(X2,PILA.40.17.Ht,type='l',col='orange',lwd=2,lty=2)
points(X2,PILA.45.17.Ht,type='l',col='green',lwd=2)
points(X2,PILA.50.17.Ht,type='l',col='light blue',lwd=2,lty=2)
points(X2,PILA.55.17.Ht,type='l',col='blue',lwd=2,lty=2)
points(X2,PILA.60.17.Ht,type='l',col='purple',lwd=2,lty=2)
legend(2.8,30,legend=c('All',Sdl.elev.17),col=c('black','red','orange','green','light blue','blue','purple'),lty=c(1,2,2,1,2,2,2),lwd=c(2,2,2,2,2,2,2),cex=0.7)

plot(X2[1:4],PILA.18.Ht,type='l',lty=1,xlab='Year Since Planting',main='2018 PILA by elev',ylab='Height (cm)',xlim=c(0,4),lwd=2,ylim=c(0,60))
points(X2[1:4],PILA.45.18.Ht,type='l',col='green',lwd=2)
points(X2[1:4],PILA.50.18.Ht,type='l',col='light blue',lwd=2,lty=2)
points(X2[1:4],PILA.52.18.Ht,type='l',col='cornflowerblue',lwd=2,lty=2)
points(X2[1:4],PILA.55.18.Ht,type='l',col='blue',lwd=2,lty=2)
points(X2[1:4],PILA.60.18.Ht,type='l',col='purple',lwd=2,lty=2)
legend(2.8,30,legend=c('All',Sdl.elev.18),col=c('black','green','light blue','cornflowerblue','blue','purple'),lty=c(1,1,2,2,2,2),lwd=c(2,2,2,2,2,2),cex=0.7)

plot(X2[1:3],PILA.19.Ht,type='l',lty=1,xlab='Year Since Planting',main='2019 PILA by elev',ylab='Height (cm)',ylim=c(0,60),lwd=2,xlim=c(0,4))
points(X2[1:3],PILA.45.19.Ht,type='l',col='green',lwd=2)
points(X2[1:3],PILA.50.19.Ht,type='l',col='light blue',lwd=2,lty=2)
points(X2[1:3],PILA.55.19.Ht,type='l',col='blue',lwd=2,lty=2)
points(X2[1:3],PILA.60.19.Ht,type='l',col='purple',lwd=2,lty=2)
points(X2[1:3],PILA.75.19.Ht,type='l',col='darkorchid4',lwd=2,lty=2)
legend(2.8,30,legend=c('All',Sdl.elev.19),col=c('black','green','light blue','blue','purple','darkorchid4'),lty=c(1,1,2,2,2,2),lwd=c(2,2,2,2,2,2),cex=0.7)

dev.off()



#########Growth matrix

PILA.19.Grow <- cbind((PILA.19.Htmat[,2]-PILA.19.Htmat[,1]),(PILA.19.Htmat[,3]-PILA.19.Htmat[,2]))
PILA.18.Grow <- cbind((PILA.18.Htmat[,2]-PILA.18.Htmat[,1]),(PILA.18.Htmat[,3]-PILA.18.Htmat[,2]),(PILA.18.Htmat[,4]-PILA.18.Htmat[,3]))
PILA.17.Grow <- cbind((PILA.17.Htmat[,2]-PILA.17.Htmat[,1]),(PILA.17.Htmat[,3]-PILA.17.Htmat[,2]),(PILA.17.Htmat[,4]-PILA.17.Htmat[,3]),(PILA.17.Htmat[,5]-PILA.17.Htmat[,4]))

#######Elevation & seed zones

PILA.17.Elev <- matrix(0,N.PILA.17,8); PILA.18.Elev <- matrix(0,N.PILA.18,8)
PILA.19.Elev <- matrix(0,N.PILA.19,8)

PILA.17.Elev[PILA.30.17,1] <- 1 ; PILA.17.Elev[PILA.40.17,2] <- 1
PILA.17.Elev[PILA.45.17,3] <- 1 ; PILA.17.Elev[PILA.50.17,4] <- 1
PILA.17.Elev[PILA.55.17,6] <- 1 ; PILA.17.Elev[PILA.60.17,7] <- 1
PILA.18.Elev[PILA.45.18,3] <- 1 ; PILA.18.Elev[PILA.50.18,4] <- 1
PILA.18.Elev[PILA.52.18,5] <- 1
PILA.18.Elev[PILA.55.18,6] <- 1 ; PILA.18.Elev[PILA.60.18,7] <- 1
PILA.19.Elev[PILA.45.19,3] <- 1 ; PILA.19.Elev[PILA.50.19,4] <- 1
PILA.19.Elev[PILA.55.19,6] <- 1 ; PILA.19.Elev[PILA.60.19,7] <- 1
PILA.19.Elev[PILA.75.19,7] <- 1

PILA.526.17 <- which(PILA.17.Sdl$SZ == 526); PILA.533.17 <- which(PILA.17.Sdl$SZ == 533)
PILA.531.17 <- which(PILA.17.Sdl$SZ == 531); PILA.540.17 <- which(PILA.17.Sdl$SZ == 540)
PILA.532.17 <- which(PILA.17.Sdl$SZ == 532) 
PILA.526.18 <- which(PILA.18.Sdl$SZ == 526); PILA.533.18 <- which(PILA.18.Sdl$SZ == 533)
PILA.531.18 <- which(PILA.18.Sdl$SZ == 531); PILA.540.18 <- which(PILA.18.Sdl$SZ == 540)
PILA.532.18 <- which(PILA.18.Sdl$SZ == 532) 
PILA.526.19 <- which(PILA.19.Sdl$SZ == 526); PILA.533.19 <- which(PILA.19.Sdl$SZ == 533)
PILA.531.19 <- which(PILA.19.Sdl$SZ == 531); PILA.540.19 <- which(PILA.19.Sdl$SZ == 540)
PILA.532.19 <- which(PILA.19.Sdl$SZ == 532); PILA.997.19 <- which(PILA.19.Sdl$SZ == 997)

PILA.17.SZ <- matrix(0,N.PILA.17,6); PILA.18.SZ <- matrix(0,N.PILA.18,6)
PILA.19.SZ <- matrix(0,N.PILA.19,6)
PILA.17.SZ[PILA.526.17,1] <- 1 ; PILA.17.SZ[PILA.531.17,2] <- 1
PILA.17.SZ[PILA.532.17,3] <- 1 ; PILA.17.SZ[PILA.533.17,4] <- 1
PILA.17.SZ[PILA.540.17,5] <- 1 
PILA.18.SZ[PILA.526.18,1] <- 1 ; PILA.18.SZ[PILA.531.18,2] <- 1
PILA.18.SZ[PILA.532.18,3] <- 1 ; PILA.18.SZ[PILA.533.18,4] <- 1
PILA.18.SZ[PILA.540.18,5] <- 1 
PILA.19.SZ[PILA.526.19,1] <- 1 ; PILA.19.SZ[PILA.531.19,2] <- 1
PILA.19.SZ[PILA.532.19,3] <- 1 ; PILA.19.SZ[PILA.533.19,4] <- 1
PILA.19.SZ[PILA.540.19,5] <- 1 ; PILA.19.SZ[PILA.997.19,6] <- 1



#########Setting up x's and y's

Y.PILA.surv.all <- numeric(0)  #survival vector over all years
Y.PILA.grow.all <- numeric(0)  #growth vector over all years

X.Sites.S <- numeric(0) #Which planting site (survival)?
X.Yr.S <- numeric(0) #Which year/age (survival)?
X.Ht.S <- numeric(0) #Seedling height (survival)
X.Elev.S <- numeric(0) #Elevation bands (survival)
X.SZ.S <- numeric(0) #Seed Zone (survival)
X.CWD.S <- numeric(0) #CWD (survival)
X.P.S <- numeric(0) #Precip (survival)
X.JMin.S <- numeric(0) #JMin (survival)
X.JMax.S <- numeric(0) #JMax (survival)
X.Sno.S <- numeric(0) #Snow (survival)

X.Sites.G <- numeric(0) #Which planting site (growth)?
X.Yr.G <- numeric(0) #Which year/age (growth)?
X.Ht.G <- numeric(0) #Seedling height (growth)
X.Elev.G <- numeric(0) #Elevation bands (growth)
X.SZ.G <- numeric(0) #Seed Zone (survival)
X.CWD.G <- numeric(0) #CWD (growth)
X.P.G <- numeric(0) #Precip (growth)
X.JMin.G <- numeric(0) #JMin (growth)
X.JMax.G <- numeric(0) #JMax (growth)
X.Sno.G <- numeric(0) #Snow (growth)

for(t in 1:5){ ###looping over years
	
	if(t==1){    ###Year 1, 2017 site only
		#Survival.all
		surv.temp <- PILA.17.survmat[,t]
		N.sdl.17 <- length(surv.temp)

		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp)
	
		#planting site
		site.temp.S <- matrix(0,N.sdl.17,3); site.temp.S[,1] <- 1   #(1,0,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.17,5); Yr.temp.S[,1] <- 1   #(1,0,0,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.17.Htmat[,1])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.17.Elev)
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.17.SZ)
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[26,2],N.PILA.17)); X.P.S <- c(X.P.S,rep(Precip[26,2],N.PILA.17))
		X.Sno.S <- c(X.Sno.S,rep(Snow[26,2],N.PILA.17))
		X.JMin.S <- c(X.JMin.S,rep(JMin[26,2],N.PILA.17)); X.JMax.S <- c(X.JMax.S,rep(JMax[26,2],N.PILA.17))
	}
	
	if(t==2){   ###Year 2, 2017 and 2018 sites
		site.17.live.2 <- which(PILA.17.survmat[,1]==1) #which seedlings alive in year 1, 2017 site			
		surv.temp.1 <- PILA.17.survmat[site.17.live.2,t]
		surv.temp.2 <- PILA.18.survmat[,(t-1)]
		N.sdl.17 <- length(surv.temp.1); N.sdl.18 <- length(surv.temp.2)
		
		site.17.live.2G <- which(PILA.17.survmat[,1]==1 & PILA.17.survmat[,2]==1) #which seedlings alive in year 1 &2, 2017 site	
		gro.temp.1 <- PILA.17.Grow[site.17.live.2G,(t-1)]
		N.sdl.17G <- length(gro.temp.1)

		###2017 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.1)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.17,3); site.temp.S[,1] <- 1   #(1,0,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.17,5); Yr.temp.S[,2] <- 1   #(0,1,0,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.17G,5); Yr.temp.G[,1] <- 1   #(1,0,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.17.Htmat[site.17.live.2,1])
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[site.17.live.2G,1])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.17.Elev[site.17.live.2,])
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev[site.17.live.2G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.17.SZ[site.17.live.2,])
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ[site.17.live.2G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[27,2],N.sdl.17)); X.P.S <- c(X.P.S,rep(Precip[27,2],N.sdl.17))
		X.Sno.S <- c(X.Sno.S,rep(Snow[27,2],N.sdl.17))
		X.JMin.S <- c(X.JMin.S,rep(JMin[27,2],N.sdl.17)); X.JMax.S <- c(X.JMax.S,rep(JMax[27,2],N.sdl.17))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[27,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[27,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[27,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[27,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[27,2],N.sdl.17G))
		
		###2018 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.2)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.18,3); site.temp.S[,2] <- 1   #(0,1,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.18,5); Yr.temp.S[,1] <- 1   #(1,0,0,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.18.Htmat[,1])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.18.Elev)
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.18.SZ)
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[27,3],N.PILA.18)); X.P.S <- c(X.P.S,rep(Precip[27,3],N.PILA.18))
		X.Sno.S <- c(X.Sno.S,rep(Snow[27,3],N.PILA.18))
		X.JMin.S <- c(X.JMin.S,rep(JMin[27,3],N.PILA.18)); X.JMax.S <- c(X.JMax.S,rep(JMax[27,3],N.PILA.18))
	}
	
	if(t==3){   ###Year 3, all sites
		site.17.live.3 <- which(PILA.17.survmat[,2]==1) #which seedlings alive in year 2, 2017 site	
		site.18.live.3 <- which(PILA.18.survmat[,1]==1) #which seedlings alive in year 2, 2018 site			
		surv.temp.1 <- PILA.17.survmat[site.17.live.3,t]
		surv.temp.2 <- PILA.18.survmat[site.18.live.3,(t-1)]
		surv.temp.3 <- PILA.19.survmat[,(t-2)]
		N.sdl.17 <- length(surv.temp.1); N.sdl.18 <- length(surv.temp.2); N.sdl.19 <- length(surv.temp.3)
		
		site.17.live.3G <- which(PILA.17.survmat[,2]==1 & PILA.17.survmat[,3]==1) #which seedlings alive in year 2&3, 2017 site	
		gro.temp.1 <- PILA.17.Grow[site.17.live.3G,(t-1)]; N.sdl.17G <- length(gro.temp.1)
		site.18.live.3G <- which(PILA.18.survmat[,1]==1 & PILA.18.survmat[,2]==1) #which seedlings alive in year 2&3, 2018 site	
		gro.temp.2 <- PILA.18.Grow[site.18.live.3G,(t-2)]; N.sdl.18G <- length(gro.temp.2)

		###2017 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.1)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.17,3); site.temp.S[,1] <- 1   #(1,0,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.17,5); Yr.temp.S[,3] <- 1   #(0,0,1,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.17G,5); Yr.temp.G[,2] <- 1   #(0,1,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.17.Htmat[site.17.live.3,2])
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[site.17.live.3G,2])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.17.Elev[site.17.live.3,])
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev[site.17.live.3G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.17.SZ[site.17.live.3,])
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ[site.17.live.3G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[28,2],N.sdl.17)); X.P.S <- c(X.P.S,rep(Precip[28,2],N.sdl.17))
		X.Sno.S <- c(X.Sno.S,rep(Snow[28,2],N.sdl.17))
		X.JMin.S <- c(X.JMin.S,rep(JMin[28,2],N.sdl.17)); X.JMax.S <- c(X.JMax.S,rep(JMax[28,2],N.sdl.17))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[28,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[28,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[28,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[28,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[28,2],N.sdl.17G))

		
		###2018 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.2)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.2)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.18,3); site.temp.S[,2] <- 1   #(0,1,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.18G,3); site.temp.G[,2] <- 1   #(0,1,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.18,5); Yr.temp.S[,2] <- 1   #(0,1,0,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.18G,5); Yr.temp.G[,1] <- 1   #(1,0,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.18.Htmat[site.18.live.3,1])
		X.Ht.G <- c(X.Ht.G,PILA.18.Htmat[site.18.live.3G,1])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.18.Elev[site.18.live.3,])
		X.Elev.G <- rbind(X.Elev.G,PILA.18.Elev[site.18.live.3G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.18.SZ[site.18.live.3,])
		X.SZ.G <- rbind(X.SZ.G,PILA.18.SZ[site.18.live.3G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[28,3],N.sdl.18)); X.P.S <- c(X.P.S,rep(Precip[28,3],N.sdl.18))
		X.Sno.S <- c(X.Sno.S,rep(Snow[28,3],N.sdl.18))
		X.JMin.S <- c(X.JMin.S,rep(JMin[28,3],N.sdl.18)); X.JMax.S <- c(X.JMax.S,rep(JMax[28,3],N.sdl.18))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[28,3],N.sdl.18G)); X.P.G <- c(X.P.G,rep(Precip[28,3],N.sdl.18G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[28,3],N.sdl.18G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[28,3],N.sdl.18G)); X.JMax.G <- c(X.JMax.G,rep(JMax[28,3],N.sdl.18G))
		
		###2019 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.3)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.19,3); site.temp.S[,3] <- 1   #(0,0,1)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.19,5); Yr.temp.S[,1] <- 1   #(1,0,0,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.19.Htmat[,1])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.19.Elev)
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.19.SZ)
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[28,4],N.PILA.19)); X.P.S <- c(X.P.S,rep(Precip[28,4],N.PILA.19))
		X.Sno.S <- c(X.Sno.S,rep(Snow[28,4],N.PILA.19))
		X.JMin.S <- c(X.JMin.S,rep(JMin[28,4],N.PILA.19)); X.JMax.S <- c(X.JMax.S,rep(JMax[28,4],N.PILA.19))
		
	}	

	if(t==4){   ###Year 4, all sites
		site.17.live.4 <- which(PILA.17.survmat[,3]==1) #which seedlings alive in year 3, 2017 site	
		site.18.live.4 <- which(PILA.18.survmat[,2]==1) #which seedlings alive in year 3, 2018 site
		site.19.live.4 <- which(PILA.19.survmat[,1]==1) #which seedlings alive in year 3, 2019 site				
		surv.temp.1 <- PILA.17.survmat[site.17.live.4,t]
		surv.temp.2 <- PILA.18.survmat[site.18.live.4,(t-1)]
		surv.temp.3 <- PILA.19.survmat[site.19.live.4,(t-2)]
		N.sdl.17 <- length(surv.temp.1); N.sdl.18 <- length(surv.temp.2); N.sdl.19 <- length(surv.temp.3)
		
		site.17.live.4G <- which(PILA.17.survmat[,3]==1 & PILA.17.survmat[,4]==1) #which seedlings alive in year 3&4, 2017 site	
		gro.temp.1 <- PILA.17.Grow[site.17.live.4G,(t-1)]; N.sdl.17G <- length(gro.temp.1)
		site.18.live.4G <- which(PILA.18.survmat[,2]==1 & PILA.18.survmat[,3]==1) #which seedlings alive in year 3&4, 2018 site	
		gro.temp.2 <- PILA.18.Grow[site.18.live.4G,(t-2)]; N.sdl.18G <- length(gro.temp.2)
		site.19.live.4G <- which(PILA.19.survmat[,1]==1 & PILA.19.survmat[,2]==1) #which seedlings alive in year 3&4, 2019 site	
		gro.temp.3 <- PILA.19.Grow[site.19.live.4G,(t-3)]; N.sdl.19G <- length(gro.temp.3)

		###2017 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.1)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.17,3); site.temp.S[,1] <- 1   #(1,0,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
				
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.17,5); Yr.temp.S[,4] <- 1   #(0,0,0,1,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.17G,5); Yr.temp.G[,3] <- 1   #(0,0,1,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.17.Htmat[site.17.live.4,3])
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[site.17.live.4G,3])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.17.Elev[site.17.live.4,])
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev[site.17.live.4G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.17.SZ[site.17.live.4,])
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ[site.17.live.4G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[29,2],N.sdl.17)); X.P.S <- c(X.P.S,rep(Precip[29,2],N.sdl.17))
		X.Sno.S <- c(X.Sno.S,rep(Snow[29,2],N.sdl.17))
		X.JMin.S <- c(X.JMin.S,rep(JMin[29,2],N.sdl.17)); X.JMax.S <- c(X.JMax.S,rep(JMax[29,2],N.sdl.17))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[29,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[29,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[29,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[29,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[29,2],N.sdl.17G))
		
		###2018 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.2)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.2)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.18,3); site.temp.S[,2] <- 1   #(0,1,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.18G,3); site.temp.G[,2] <- 1   #(0,1,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.18,5); Yr.temp.S[,3] <- 1   #(0,0,1,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.18G,5); Yr.temp.G[,2] <- 1   #(0,1,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.18.Htmat[site.18.live.4,2])
		X.Ht.G <- c(X.Ht.G,PILA.18.Htmat[site.18.live.4G,2])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.18.Elev[site.18.live.4,])
		X.Elev.G <- rbind(X.Elev.G,PILA.18.Elev[site.18.live.4G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.18.SZ[site.18.live.4,])
		X.SZ.G <- rbind(X.SZ.G,PILA.18.SZ[site.18.live.4G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[29,3],N.sdl.18)); X.P.S <- c(X.P.S,rep(Precip[29,3],N.sdl.18))
		X.Sno.S <- c(X.Sno.S,rep(Snow[29,3],N.sdl.18))
		X.JMin.S <- c(X.JMin.S,rep(JMin[29,3],N.sdl.18)); X.JMax.S <- c(X.JMax.S,rep(JMax[29,3],N.sdl.18))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[29,3],N.sdl.18G)); X.P.G <- c(X.P.G,rep(Precip[29,3],N.sdl.18G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[29,3],N.sdl.18G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[29,3],N.sdl.18G)); X.JMax.G <- c(X.JMax.G,rep(JMax[29,3],N.sdl.18G))
		
		###2019 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.3)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.3)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.19,3); site.temp.S[,3] <- 1   #(0,0,1)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.19G,3); site.temp.G[,3] <- 1   #(0,0,1)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.19,5); Yr.temp.S[,2] <- 1   #(0,1,0,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.19G,5); Yr.temp.G[,1] <- 1   #(1,0,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.19.Htmat[site.19.live.4,1])
		X.Ht.G <- c(X.Ht.G,PILA.19.Htmat[site.19.live.4G,1])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.19.Elev[site.19.live.4,])
		X.Elev.G <- rbind(X.Elev.G,PILA.19.Elev[site.19.live.4G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.19.SZ[site.19.live.4,])
		X.SZ.G <- rbind(X.SZ.G,PILA.19.SZ[site.19.live.4G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[29,4],N.sdl.19)); X.P.S <- c(X.P.S,rep(Precip[29,4],N.sdl.19))
		X.Sno.S <- c(X.Sno.S,rep(Snow[29,4],N.sdl.19))
		X.JMin.S <- c(X.JMin.S,rep(JMin[29,4],N.sdl.19)); X.JMax.S <- c(X.JMax.S,rep(JMax[29,4],N.sdl.19))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[29,4],N.sdl.19G)); X.P.G <- c(X.P.G,rep(Precip[29,4],N.sdl.19G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[29,4],N.sdl.19G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[29,4],N.sdl.19G)); X.JMax.G <- c(X.JMax.G,rep(JMax[29,4],N.sdl.19G))
	}
	
	if(t==5){   ###Year 5, all sites
		site.17.live.5 <- which(PILA.17.survmat[,4]==1) #which seedlings alive in year 4, 2017 site	
		site.18.live.5 <- which(PILA.18.survmat[,3]==1) #which seedlings alive in year 4, 2018 site
		site.19.live.5 <- which(PILA.19.survmat[,2]==1) #which seedlings alive in year 4, 2019 site				
		surv.temp.1 <- PILA.17.survmat[site.17.live.5,t]
		surv.temp.2 <- PILA.18.survmat[site.18.live.5,(t-1)]
		surv.temp.3 <- PILA.19.survmat[site.19.live.5,(t-2)]
		N.sdl.17 <- length(surv.temp.1); N.sdl.18 <- length(surv.temp.2); N.sdl.19 <- length(surv.temp.3)
		
		site.17.live.5G <- which(PILA.17.survmat[,4]==1 & PILA.17.survmat[,5]==1) #which seedlings alive in year 4&5, 2017 site	
		gro.temp.1 <- PILA.17.Grow[site.17.live.5G,(t-1)]; N.sdl.17G <- length(gro.temp.1)
		site.18.live.5G <- which(PILA.18.survmat[,3]==1 & PILA.18.survmat[,4]==1) #which seedlings alive in year 4&5, 2018 site	
		gro.temp.2 <- PILA.18.Grow[site.18.live.5G,(t-2)]; N.sdl.18G <- length(gro.temp.2)
		site.19.live.5G <- which(PILA.19.survmat[,2]==1 & PILA.19.survmat[,3]==1) #which seedlings alive in year 4&5, 2019 site	
		gro.temp.3 <- PILA.19.Grow[site.19.live.5G,(t-3)]; N.sdl.19G <- length(gro.temp.3)

		###2017 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.1)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.17,3); site.temp.S[,1] <- 1   #(1,0,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.17,5); Yr.temp.S[,5] <- 1   #(0,0,0,0,1)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.17G,5); Yr.temp.G[,4] <- 1   #(0,0,0,1,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.17.Htmat[site.17.live.5,4])
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[site.17.live.5G,4])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.17.Elev[site.17.live.5,])
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev[site.17.live.5G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.17.SZ[site.17.live.5,])
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ[site.17.live.5G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[30,2],N.sdl.17)); X.P.S <- c(X.P.S,rep(Precip[30,2],N.sdl.17))
		X.Sno.S <- c(X.Sno.S,rep(Snow[30,2],N.sdl.17))
		X.JMin.S <- c(X.JMin.S,rep(JMin[30,2],N.sdl.17)); X.JMax.S <- c(X.JMax.S,rep(JMax[30,2],N.sdl.17))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[30,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[30,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[30,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[30,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[30,2],N.sdl.17G))
		
		###2018 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.2)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.2)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.18,3); site.temp.S[,2] <- 1   #(0,1,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.18G,3); site.temp.G[,2] <- 1   #(0,1,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.18,5); Yr.temp.S[,4] <- 1   #(0,0,0,1,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.18G,5); Yr.temp.G[,3] <- 1   #(0,0,1,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.18.Htmat[site.18.live.5,3])
		X.Ht.G <- c(X.Ht.G,PILA.18.Htmat[site.18.live.5G,3])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.18.Elev[site.18.live.5,])
		X.Elev.G <- rbind(X.Elev.G,PILA.18.Elev[site.18.live.5G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.18.SZ[site.18.live.5,])
		X.SZ.G <- rbind(X.SZ.G,PILA.18.SZ[site.18.live.5G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[30,3],N.sdl.18)); X.P.S <- c(X.P.S,rep(Precip[30,3],N.sdl.18))
		X.Sno.S <- c(X.Sno.S,rep(Snow[30,3],N.sdl.18))
		X.JMin.S <- c(X.JMin.S,rep(JMin[30,3],N.sdl.18)); X.JMax.S <- c(X.JMax.S,rep(JMax[30,3],N.sdl.18))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[30,3],N.sdl.18G)); X.P.G <- c(X.P.G,rep(Precip[30,3],N.sdl.18G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[30,3],N.sdl.18G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[30,3],N.sdl.18G)); X.JMax.G <- c(X.JMax.G,rep(JMax[30,3],N.sdl.18G))
		
		###2019 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.3)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.3)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.19,3); site.temp.S[,3] <- 1   #(0,0,1)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.19G,3); site.temp.G[,3] <- 1   #(0,0,1)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.19,5); Yr.temp.S[,3] <- 1   #(0,0,1,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.19G,5); Yr.temp.G[,2] <- 1   #(0,1,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.S <- c(X.Ht.S,PILA.19.Htmat[site.19.live.5,2])
		X.Ht.G <- c(X.Ht.G,PILA.19.Htmat[site.19.live.5G,2])
		
		#Elev
		X.Elev.S <- rbind(X.Elev.S,PILA.19.Elev[site.19.live.5,])
		X.Elev.G <- rbind(X.Elev.G,PILA.19.Elev[site.19.live.5G,])
		
		#SZ
		X.SZ.S <- rbind(X.SZ.S,PILA.19.SZ[site.19.live.5,])
		X.SZ.G <- rbind(X.SZ.G,PILA.19.SZ[site.19.live.5G,])
		
		#Climate
		X.CWD.S <- c(X.CWD.S,rep(CWD[30,4],N.sdl.19)); X.P.S <- c(X.P.S,rep(Precip[30,4],N.sdl.19))
		X.Sno.S <- c(X.Sno.S,rep(Snow[30,4],N.sdl.19))
		X.JMin.S <- c(X.JMin.S,rep(JMin[30,4],N.sdl.19)); X.JMax.S <- c(X.JMax.S,rep(JMax[30,4],N.sdl.19))
		
		X.CWD.G <- c(X.CWD.G,rep(CWD[30,4],N.sdl.19G)); X.P.G <- c(X.P.G,rep(Precip[30,4],N.sdl.19G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[30,4],N.sdl.19G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[30,4],N.sdl.19G)); X.JMax.G <- c(X.JMax.G,rep(JMax[30,4],N.sdl.19G))
	}
} #year loop end



N.surv <- length(Y.PILA.surv.all)
N.grow <- length(Y.PILA.grow.all)

X.Site.Elev.S <- matrix(0,N.surv,24); X.Site.Elev.G <- matrix(0,N.grow,24)

for(i in 1:N.surv){
	if(X.Elev.S[i,1]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,1] <-1
	if(X.Elev.S[i,1]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,2] <-1
	if(X.Elev.S[i,1]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,3] <-1
	if(X.Elev.S[i,2]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,4] <-1
	if(X.Elev.S[i,2]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,5] <-1
	if(X.Elev.S[i,2]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,6] <-1
	if(X.Elev.S[i,3]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,7] <-1
	if(X.Elev.S[i,3]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,8] <-1
	if(X.Elev.S[i,3]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,9] <-1
	if(X.Elev.S[i,4]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,10] <-1
	if(X.Elev.S[i,4]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,11] <-1
	if(X.Elev.S[i,4]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,12] <-1
	if(X.Elev.S[i,5]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,13] <-1
	if(X.Elev.S[i,5]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,14] <-1
	if(X.Elev.S[i,5]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,15] <-1
	if(X.Elev.S[i,6]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,16] <-1
	if(X.Elev.S[i,6]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,17] <-1
	if(X.Elev.S[i,6]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,18] <-1
	if(X.Elev.S[i,7]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,19] <-1
	if(X.Elev.S[i,7]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,20] <-1
	if(X.Elev.S[i,7]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,21] <-1
	if(X.Elev.S[i,8]==1 & X.Sites.S[i,1]==1) X.Site.Elev.S[i,22] <-1
	if(X.Elev.S[i,8]==1 & X.Sites.S[i,2]==1) X.Site.Elev.S[i,23] <-1
	if(X.Elev.S[i,8]==1 & X.Sites.S[i,3]==1) X.Site.Elev.S[i,24] <-1
}

for(i in 1:N.grow){
	if(X.Elev.G[i,1]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,1] <-1
	if(X.Elev.G[i,1]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,2] <-1
	if(X.Elev.G[i,1]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,3] <-1
	if(X.Elev.G[i,2]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,4] <-1
	if(X.Elev.G[i,2]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,5] <-1
	if(X.Elev.G[i,2]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,6] <-1
	if(X.Elev.G[i,3]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,7] <-1
	if(X.Elev.G[i,3]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,8] <-1
	if(X.Elev.G[i,3]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,9] <-1
	if(X.Elev.G[i,4]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,10] <-1
	if(X.Elev.G[i,4]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,11] <-1
	if(X.Elev.G[i,4]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,12] <-1
	if(X.Elev.G[i,5]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,13] <-1
	if(X.Elev.G[i,5]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,14] <-1
	if(X.Elev.G[i,5]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,15] <-1
	if(X.Elev.G[i,6]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,16] <-1
	if(X.Elev.G[i,6]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,17] <-1
	if(X.Elev.G[i,6]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,18] <-1
	if(X.Elev.G[i,7]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,19] <-1
	if(X.Elev.G[i,7]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,20] <-1
	if(X.Elev.G[i,7]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,21] <-1
	if(X.Elev.G[i,8]==1 & X.Sites.G[i,1]==1) X.Site.Elev.G[i,22] <-1
	if(X.Elev.G[i,8]==1 & X.Sites.G[i,2]==1) X.Site.Elev.G[i,23] <-1
	if(X.Elev.G[i,8]==1 & X.Sites.G[i,3]==1) X.Site.Elev.G[i,24] <-1
}


X.SZ.Elev.S <- matrix(0,N.surv,16); X.SZ.Elev.G <- matrix(0,N.grow,16)

for(i in 1:N.surv){
	if(X.SZ.S[i,1]==1 & X.Elev.S[i,3]==1) X.SZ.Elev.S[i,1] <-1 #526_4500
    if(X.SZ.S[i,1]==1 & X.Elev.S[i,4]==1) X.SZ.Elev.S[i,2] <-1 #526_5000
    if(X.SZ.S[i,1]==1 & X.Elev.S[i,5]==1) X.SZ.Elev.S[i,3] <-1 #526_5250
    if(X.SZ.S[i,1]==1 & X.Elev.S[i,6]==1) X.SZ.Elev.S[i,4] <-1 #526_5500
    if(X.SZ.S[i,2]==1 & X.Elev.S[i,1]==1) X.SZ.Elev.S[i,5] <-1 #531_3000
    if(X.SZ.S[i,2]==1 & X.Elev.S[i,2]==1) X.SZ.Elev.S[i,6] <-1 #531_4000
    if(X.SZ.S[i,2]==1 & X.Elev.S[i,3]==1) X.SZ.Elev.S[i,7] <-1 #531_4500
    if(X.SZ.S[i,2]==1 & X.Elev.S[i,4]==1) X.SZ.Elev.S[i,8] <-1 #531_5000
    if(X.SZ.S[i,2]==1 & X.Elev.S[i,6]==1) X.SZ.Elev.S[i,9] <-1 #531_5500
    if(X.SZ.S[i,2]==1 & X.Elev.S[i,7]==1) X.SZ.Elev.S[i,10] <-1 #531_6000
    if(X.SZ.S[i,3]==1 & X.Elev.S[i,4]==1) X.SZ.Elev.S[i,11] <-1 #532_5000
    if(X.SZ.S[i,3]==1 & X.Elev.S[i,6]==1) X.SZ.Elev.S[i,12] <-1 #532_5500
    if(X.SZ.S[i,3]==1 & X.Elev.S[i,7]==1) X.SZ.Elev.S[i,13] <-1 #532_6000
    if(X.SZ.S[i,4]==1 & X.Elev.S[i,7]==1) X.SZ.Elev.S[i,14] <-1 #533_6000
    if(X.SZ.S[i,5]==1 & X.Elev.S[i,7]==1) X.SZ.Elev.S[i,15] <-1 #540_6000
	if(X.SZ.S[i,6]==1 & X.Elev.S[i,8]==1) X.SZ.Elev.S[i,16] <-1 #540_7500
	}
	

save.image(file="PILASdl_Data.RData")