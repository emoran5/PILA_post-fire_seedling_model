##Basing this on TagSdlDataOrg5 in **

bframe1 <-data.frame(read.table("Peavine_PILA_dec2_trim.txt",header=T,fill=T,row.names=NULL))
bframe2 <-data.frame(read.table("Stumpy2018_PILAdec2_trim.txt",header=T,fill=T,row.names=NULL)) 
bframe3 <-data.frame(read.table("Stumpy2017_PILAdec2_trim.txt",header=T,fill=T,row.names=NULL))

CWD <- data.frame(read.table("CWD_91_21.txt",header=T,fill=T,row.names=NULL))
Snow <- data.frame(read.table("Sno_91_21.txt",header=T,fill=T,row.names=NULL))
Precip <- data.frame(read.table("Precip_91_21.txt",header=T,fill=T,row.names=NULL))
JMin <- data.frame(read.table("JMin_91_21.txt",header=T,fill=T,row.names=NULL))
JMax <- data.frame(read.table("JMax_91_21.txt",header=T,fill=T,row.names=NULL))

Sno.Geno <-data.frame(read.table("favor_allele_Sno2.txt",header=T,fill=T,row.names=NULL)) #185 loci
JMx.Geno <-data.frame(read.table("favor_allele_Tmax2.txt",header=T,fill=T,row.names=NULL)) #82 loci
JMn.Geno <-data.frame(read.table("favor_allele_Tmn2.txt",header=T,fill=T,row.names=NULL)) #15 loci
CWD.Geno <-data.frame(read.table("favor_allele_Cwd2.txt",header=T,fill=T,row.names=NULL)) #20 loci
P.Geno <-data.frame(read.table("favor_allele_Ppt2.txt",header=T,fill=T,row.names=NULL)) #3  loci

#### Genotyped PILA seedlings 

Geno.19 <- which(bframe1$Genotype==1); Geno.18 <- which(bframe2$Genotype==1)
Geno.17 <- which(bframe3$Genotype==1)
#65+55+41 =161

NSdl <- nrow(Sno.Geno)


PILA.19.Sdl <- cbind(bframe1[Geno.19,4],bframe1[Geno.19,7],bframe1[Geno.19,9],bframe1[Geno.19,3],bframe1[Geno.19,12:16],bframe1[Geno.19,18],bframe1[Geno.19,23],bframe1[Geno.19,10])
PILA.18.Sdl <- cbind(bframe2[Geno.18,3],bframe2[Geno.18,6],bframe2[Geno.18,8],bframe2[Geno.18,2],bframe2[Geno.18,12:15],bframe2[Geno.18,18],bframe2[Geno.18,21],bframe2[Geno.18,25],bframe2[Geno.18,10])
PILA.17.Sdl <- cbind(bframe3[Geno.17,3:4],bframe3[Geno.17,6],bframe3[Geno.17,2],bframe3[Geno.17,11:13],bframe3[Geno.17,16],bframe3[Geno.17,19],bframe3[Geno.17,22],bframe3[Geno.17,26],bframe3[Geno.17,9])

colnames(PILA.19.Sdl)<-c('Species','SZ','ElevB','Lot_Name','PlantYr','MortYr','Ht17','Ht18','Ht19','Ht20','Ht21','Tag')
colnames(PILA.18.Sdl)<-c('Species','SZ','ElevB','Lot_Name','PlantYr','MortYr','Ht17','Ht18','Ht19','Ht20','Ht21','Tag')
colnames(PILA.17.Sdl)<-c('Species','SZ','ElevB','Lot_Name','PlantYr','MortYr','Ht17','Ht18','Ht19','Ht20','Ht21','Tag')

#seedling numbers
N.PILA.19 <- nrow(PILA.19.Sdl) #65
N.PILA.18 <- nrow(PILA.18.Sdl) #55
N.PILA.17 <- nrow(PILA.17.Sdl) #41

PILA.19.survmat <- matrix(1,N.PILA.19,3);PILA.19.Htmat <- matrix(NA,N.PILA.19,3)
PILA.18.survmat <- matrix(1,N.PILA.18,4);PILA.18.Htmat <- matrix(NA,N.PILA.18,4)
PILA.17.survmat <- matrix(1,N.PILA.17,5);PILA.17.Htmat <- matrix(NA,N.PILA.17,5)


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
			if(PILA.19.Sdl[i,'Ht21']<0) PILA.19.Htmat[i,2] <- PILA.19.Htmat[i,1]  ##update
		}
	#Year 3 ht
	if(PILA.19.Sdl[i,'Ht21']>0) PILA.19.Htmat[i,3] <- PILA.19.Sdl[i,'Ht21']
		if(PILA.19.Sdl[i,'Ht21']<0 & PILA.19.survmat[i,3]==1) PILA.19.Htmat[i,3] <- 
			PILA.19.Htmat[i,2] ##update
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


PILA.30.17 <- which(PILA.17.Sdl$ElevB == 3000) # 1 sdl
PILA.40.17 <- which(PILA.17.Sdl$ElevB == 4000)  # 1 sdl
PILA.45.17 <- which(PILA.17.Sdl$ElevB == 4500) #23 sdl
PILA.50.17 <- which(PILA.17.Sdl$ElevB == 5000) #4 sdl
PILA.55.17 <- which(PILA.17.Sdl$ElevB == 5500) #2 sdl
PILA.60.17 <- which(PILA.17.Sdl$ElevB == 6000) #10 sdl

PILA.45.18 <- which(PILA.18.Sdl$ElevB == 4500) #3 sdl
PILA.50.18 <- which(PILA.18.Sdl$ElevB == 5000) #14 sdl
PILA.52.18 <- which(PILA.18.Sdl$ElevB == 5250) #2 sdl
PILA.55.18 <- which(PILA.18.Sdl$ElevB == 5500) #14 sdl
PILA.60.18 <- which(PILA.18.Sdl$ElevB == 6000) #22 sdl

PILA.45.19 <- which(PILA.19.Sdl$ElevB == 4500) #5 sdl
PILA.50.19 <- which(PILA.19.Sdl$ElevB == 5000) #12 sdl
PILA.55.19 <- which(PILA.19.Sdl$ElevB == 5500) #18 sdl
PILA.60.19 <- which(PILA.19.Sdl$ElevB == 6000) #26 sdl
PILA.75.19 <- which(PILA.19.Sdl$ElevB == 7500) #4 sdl


#########Growth matrix

PILA.19.Grow <- cbind((PILA.19.Htmat[,2]-PILA.19.Htmat[,1]),(PILA.19.Htmat[,3]-PILA.19.Htmat[,2]))
PILA.18.Grow <- cbind((PILA.18.Htmat[,2]-PILA.18.Htmat[,1]),(PILA.18.Htmat[,3]-PILA.18.Htmat[,2]),(PILA.18.Htmat[,4]-PILA.18.Htmat[,3]))
PILA.17.Grow <- cbind((PILA.17.Htmat[,2]-PILA.17.Htmat[,1]),(PILA.17.Htmat[,3]-PILA.17.Htmat[,2]),(PILA.17.Htmat[,4]-PILA.17.Htmat[,3]),(PILA.17.Htmat[,5]-PILA.17.Htmat[,4]))


#######Elevation

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

#########Genotypes
##Adjusting tags
PILA.17.Sdl[4,12] <- 219; PILA.17.Sdl[34,12] <- 810
PILA.17.Sdl[35,12] <- 831

#Matrices to hold proportion favored genotypes by climate gradient
PILA.17.Geno <- matrix(0,N.PILA.17,5) #Cols Sno,JMax,CWD,JMin,Precip
PILA.18.Geno <- matrix(0,N.PILA.18,5)
PILA.19.Geno <- matrix(0,N.PILA.19,5)

#Matrices to hold proportion of loci genotyped
PILA.17.PG <- matrix(0,N.PILA.17,5) 
PILA.18.PG <- matrix(0,N.PILA.18,5)
PILA.19.PG <- matrix(0,N.PILA.19,5)

NLoc.Sno <- (ncol(Sno.Geno)-1) #185 loci
NLoc.JMx <- (ncol(JMx.Geno)-1) #82 loci
NLoc.CWD <- (ncol(CWD.Geno)-1) #20 loci
NLoc.JMn <- (ncol(JMn.Geno)-1) #15 loci
NLoc.P <- (ncol(P.Geno)-1) #3 loci

###alternate version
#Matrices to hold proportion favored genotypes by climate gradient
PILA.17.Geno2 <- matrix(0,N.PILA.17,5) #Cols Sno,JMax,CWD,JMin,Precip
PILA.18.Geno2 <- matrix(0,N.PILA.18,5)
PILA.19.Geno2 <- matrix(0,N.PILA.19,5)


for (i in 1:N.PILA.19){
		TS <- which(Sno.Geno$ID==PILA.19.Sdl[i,12]) #Find tag 
	
	#Snow-associated loci 
		he.S <- length(which(Sno.Geno[TS,]==1)) #heterozy for favored
		ho.S <- length(which(Sno.Geno[TS,]==2)) #homozy for favored
		none.S <- length(which(Sno.Geno[TS,]==0)) #no favored
		T.g.S <- he.S + ho.S + none.S  #number of loci genotyped
		FG.S <- ((2*ho.S)+he.S)/(2*T.g.S) #Proportion favored alleles
		FG.S2 <- (ho.S+he.S)/T.g.S #Proportion favored alleles
	PILA.19.Geno[i,1] <- FG.S;  PILA.19.Geno2[i,1] <- FG.S2
	PILA.19.PG[i,1] <- T.g.S/NLoc.Sno 
	
	#JMx-associated loci 
		he.JMx <- length(which(JMx.Geno[TS,]==1)) #heterozy for favored
		ho.JMx <- length(which(JMx.Geno[TS,]==2)) #homozy for favored
		none.JMx <- length(which(JMx.Geno[TS,]==0)) #no favored
		T.g.JMx <- he.JMx + ho.JMx + none.JMx  #number of loci genotyped
		FG.JMx <- ((2*ho.JMx)+he.JMx)/(2*T.g.JMx) #Proportion favored alleles
		FG.JMx2 <- (ho.JMx+he.JMx)/T.g.JMx #Proportion favored alleles
	PILA.19.Geno[i,2] <- FG.JMx;  PILA.19.Geno2[i,2] <- FG.JMx2
	PILA.19.PG[i,2] <- T.g.JMx/NLoc.JMx 
	
	#CWD-associated loci 
		he.CWD <- length(which(CWD.Geno[TS,]==1)) #heterozy for favored
		ho.CWD <- length(which(CWD.Geno[TS,]==2)) #homozy for favored
		none.CWD <- length(which(CWD.Geno[TS,]==0)) #no favored
		T.g.CWD <- he.CWD + ho.CWD + none.CWD  #number of loci genotyped
		FG.CWD <- ((2*ho.CWD)+he.CWD)/(2*T.g.CWD) #Proportion favored alleles
		FG.CWD2 <- (ho.CWD+he.CWD)/T.g.CWD #Proportion favored alleles
	PILA.19.Geno[i,3] <- FG.CWD;  PILA.19.Geno2[i,3] <- FG.CWD2
	PILA.19.PG[i,3] <- T.g.CWD/NLoc.CWD 
	
	#JMn-associated loci 
		he.JMn <- length(which(JMn.Geno[TS,]==1)) #heterozy for favored
		ho.JMn <- length(which(JMn.Geno[TS,]==2)) #homozy for favored
		none.JMn <- length(which(JMn.Geno[TS,]==0)) #no favored
		T.g.JMn <- he.JMn + ho.JMn + none.JMn  #number of loci genotyped
		FG.JMn <- ((2*ho.JMn)+he.JMn)/(2*T.g.JMn) #Proportion favored alleles
		FG.JMn2 <- (ho.JMn+he.JMn)/T.g.JMn #Proportion favored alleles
	PILA.19.Geno[i,4] <- FG.JMn;    PILA.19.Geno2[i,4] <- FG.JMn2
	PILA.19.PG[i,4] <- T.g.JMn/NLoc.JMn 	
	
	#Precip-associated loci 
		he.P <- length(which(P.Geno[TS,]==1)) #heterozy for favored
		ho.P <- length(which(P.Geno[TS,]==2)) #homozy for favored
		none.P <- length(which(P.Geno[TS,]==0)) #no favored
		T.g.P <- he.P + ho.P + none.P  #number of loci genotyped
		FG.P <- ((2*ho.P)+he.P)/(2*T.g.P) #Proportion favored alleles
		FG.P2 <- (ho.P+he.P)/T.g.P #Proportion favored alleles
	PILA.19.Geno[i,5] <- FG.P;  PILA.19.Geno2[i,5] <- FG.P2
	PILA.19.PG[i,5] <- T.g.P/NLoc.P
}

for (i in 1:N.PILA.18){
	TS <- which(Sno.Geno$ID==PILA.18.Sdl[i,12]) #Find tag 
	
	#Snow-associated loci 
		he.S <- length(which(Sno.Geno[TS,]==1)) #heterozy for favored
		ho.S <- length(which(Sno.Geno[TS,]==2)) #homozy for favored
		none.S <- length(which(Sno.Geno[TS,]==0)) #no favored
		T.g.S <- he.S + ho.S + none.S  #number of loci genotyped
		FG.S <- ((2*ho.S)+he.S)/(2*T.g.S) #Proportion favored alleles
		FG.S2 <- (ho.S+he.S)/T.g.S #Proportion favored alleles
	PILA.18.Geno[i,1] <- FG.S;   PILA.18.Geno2[i,1] <- FG.S2
	PILA.18.PG[i,1] <- T.g.S/NLoc.Sno 
	
	#JMx-associated loci 
		he.JMx <- length(which(JMx.Geno[TS,]==1)) #heterozy for favored
		ho.JMx <- length(which(JMx.Geno[TS,]==2)) #homozy for favored
		none.JMx <- length(which(JMx.Geno[TS,]==0)) #no favored
		T.g.JMx <- he.JMx + ho.JMx + none.JMx  #number of loci genotyped
		FG.JMx <- ((2*ho.JMx)+he.JMx)/(2*T.g.JMx) #Proportion favored alleles
		FG.JMx2 <- (ho.JMx+he.JMx)/T.g.JMx #Proportion favored alleles
	PILA.18.Geno[i,2] <- FG.JMx;   PILA.18.Geno2[i,2] <- FG.JMx2
	PILA.18.PG[i,2] <- T.g.JMx/NLoc.JMx 
	
	#CWD-associated loci 
		he.CWD <- length(which(CWD.Geno[TS,]==1)) #heterozy for favored
		ho.CWD <- length(which(CWD.Geno[TS,]==2)) #homozy for favored
		none.CWD <- length(which(CWD.Geno[TS,]==0)) #no favored
		T.g.CWD <- he.CWD + ho.CWD + none.CWD  #number of loci genotyped
		FG.CWD <- ((2*ho.CWD)+he.CWD)/(2*T.g.CWD) #Proportion favored alleles
		FG.CWD2 <- (ho.CWD+he.CWD)/T.g.CWD #Proportion favored alleles
	PILA.18.Geno[i,3] <- FG.CWD;   PILA.18.Geno2[i,3] <- FG.CWD2
	PILA.18.PG[i,3] <- T.g.CWD/NLoc.CWD 
	
	#JMn-associated loci 
		he.JMn <- length(which(JMn.Geno[TS,]==1)) #heterozy for favored
		ho.JMn <- length(which(JMn.Geno[TS,]==2)) #homozy for favored
		none.JMn <- length(which(JMn.Geno[TS,]==0)) #no favored
		T.g.JMn <- he.JMn + ho.JMn + none.JMn  #number of loci genotyped
		FG.JMn <- ((2*ho.JMn)+he.JMn)/(2*T.g.JMn) #Proportion favored alleles
		FG.JMn2 <- (ho.JMn+he.JMn)/T.g.JMn #Proportion favored alleles
	PILA.18.Geno[i,4] <- FG.JMn;        PILA.18.Geno2[i,4] <- FG.JMn2
	PILA.18.PG[i,4] <- T.g.JMn/NLoc.JMn 	
	
	#Precip-associated loci 
		he.P <- length(which(P.Geno[TS,]==1)) #heterozy for favored
		ho.P <- length(which(P.Geno[TS,]==2)) #homozy for favored
		none.P <- length(which(P.Geno[TS,]==0)) #no favored
		T.g.P <- he.P + ho.P + none.P  #number of loci genotyped
		FG.P <- ((2*ho.P)+he.P)/(2*T.g.P) #Proportion favored alleles
		FG.P2 <- (ho.P+he.P)/T.g.P #Proportion favored alleles
	PILA.18.Geno[i,5] <- FG.P;   PILA.18.Geno2[i,5] <- FG.P2
	PILA.18.PG[i,5] <- T.g.P/NLoc.P	
}

for (i in 1:N.PILA.17){
	TS <- which(Sno.Geno$ID==PILA.17.Sdl[i,12]) #Find tag 
	
	#Snow-associated loci 
		he.S <- length(which(Sno.Geno[TS,]==1)) #heterozy for favored
		ho.S <- length(which(Sno.Geno[TS,]==2)) #homozy for favored
		none.S <- length(which(Sno.Geno[TS,]==0)) #no favored
		T.g.S <- he.S + ho.S + none.S  #number of loci genotyped
		FG.S <- ((2*ho.S)+he.S)/(2*T.g.S) #Proportion favored alleles
		FG.S2 <- (ho.S+he.S)/T.g.S #Proportion favored alleles
	PILA.17.Geno[i,1] <- FG.S;    PILA.17.Geno2[i,1] <- FG.S2
	PILA.17.PG[i,1] <- T.g.S/NLoc.Sno 
	
	#JMx-associated loci 
		he.JMx <- length(which(JMx.Geno[TS,]==1)) #heterozy for favored
		ho.JMx <- length(which(JMx.Geno[TS,]==2)) #homozy for favored
		none.JMx <- length(which(JMx.Geno[TS,]==0)) #no favored
		T.g.JMx <- he.JMx + ho.JMx + none.JMx  #number of loci genotyped
		FG.JMx <- ((2*ho.JMx)+he.JMx)/(2*T.g.JMx) #Proportion favored alleles
		FG.JMx2 <- (ho.JMx+he.JMx)/T.g.JMx #Proportion favored alleles
	PILA.17.Geno[i,2] <- FG.JMx;   PILA.17.Geno2[i,2] <- FG.JMx2
	PILA.17.PG[i,2] <- T.g.JMx/NLoc.JMx 
	
	#CWD-associated loci 
		he.CWD <- length(which(CWD.Geno[TS,]==1)) #heterozy for favored
		ho.CWD <- length(which(CWD.Geno[TS,]==2)) #homozy for favored
		none.CWD <- length(which(CWD.Geno[TS,]==0)) #no favored
		T.g.CWD <- he.CWD + ho.CWD + none.CWD  #number of loci genotyped
		FG.CWD <- ((2*ho.CWD)+he.CWD)/(2*T.g.CWD) #Proportion favored alleles
		FG.CWD2 <- (ho.CWD+he.CWD)/T.g.CWD #Proportion favored alleles
	PILA.17.Geno[i,3] <- FG.CWD;   PILA.17.Geno2[i,3] <- FG.CWD2
	PILA.17.PG[i,3] <- T.g.CWD/NLoc.CWD 
	
	#JMn-associated loci 
		he.JMn <- length(which(JMn.Geno[TS,]==1)) #heterozy for favored
		ho.JMn <- length(which(JMn.Geno[TS,]==2)) #homozy for favored
		none.JMn <- length(which(JMn.Geno[TS,]==0)) #no favored
		T.g.JMn <- he.JMn + ho.JMn + none.JMn  #number of loci genotyped
		FG.JMn <- ((2*ho.JMn)+he.JMn)/(2*T.g.JMn) #Proportion favored alleles
		FG.JMn2 <- (ho.JMn+he.JMn)/T.g.JMn #Proportion favored alleles
	PILA.17.Geno[i,4] <- FG.JMn;    PILA.17.Geno2[i,4] <- FG.JMn2
	PILA.17.PG[i,4] <- T.g.JMn/NLoc.JMn 	
	
	#Precip-associated loci 
		he.P <- length(which(P.Geno[TS,]==1)) #heterozy for favored
		ho.P <- length(which(P.Geno[TS,]==2)) #homozy for favored
		none.P <- length(which(P.Geno[TS,]==0)) #no favored
		T.g.P <- he.P + ho.P + none.P  #number of loci genotyped
		FG.P <- ((2*ho.P)+he.P)/(2*T.g.P) #Proportion favored alleles
		FG.P2 <- (ho.P+he.P)/T.g.P #Proportion favored alleles
	PILA.17.Geno[i,5] <- FG.P;   PILA.17.Geno2[i,5] <- FG.P2
	PILA.17.PG[i,5] <- T.g.P/NLoc.P
	
}




#########Setting up x's and y's

Y.PILA.surv.all <- numeric(0)  #survival vector over last two years
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
X.Geno.S <- numeric(0) #Genotypes (survival)
X.Geno.S2 <- numeric(0) #Genotypes (survival)

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
X.Geno.G <- numeric(0) #Genotypes (survival)
X.Geno.G2 <- numeric(0) #Genotypes (survival)

for(t in 2:5){ ###looping over years
	
	if(t==2){   ###Year 2, 2017 and 2018 sites
		gro.temp.1 <- PILA.17.Grow[,(t-1)]
		N.sdl.17G <- length(gro.temp.1)

		###2017 site seedlings
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.G <- matrix(0,N.sdl.17G,4); Yr.temp.G[,1] <- 1   #(1,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[,1])
		
		#Elev
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev)
		
		#SZ
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ)
		
		#Climate
		X.CWD.G <- c(X.CWD.G,rep(CWD[27,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[27,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[27,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[27,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[27,2],N.sdl.17G))
		
		#Genotype
		X.Geno.G <- rbind(X.Geno.G,PILA.17.Geno)	
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.17.Geno2)	
	}
	
	if(t==3){   ###Year 3, all sites	
		gro.temp.1 <- PILA.17.Grow[,(t-1)]; N.sdl.17G <- length(gro.temp.1)
		gro.temp.2 <- PILA.18.Grow[,(t-2)]; N.sdl.18G <- length(gro.temp.2)

		###2017 site seedlings
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.G <- matrix(0,N.sdl.17G,4); Yr.temp.G[,2] <- 1   #(0,1,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[,2])
		
		#Elev
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev)
		
		#SZ
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ)
		
		#Climate
		X.CWD.G <- c(X.CWD.G,rep(CWD[28,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[28,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[28,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[28,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[28,2],N.sdl.17G))

		#Genotype
		X.Geno.G <- rbind(X.Geno.G,PILA.17.Geno)
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.17.Geno2)
		
		###2018 site seedlings
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.2)
		
		#planting site
		site.temp.G <- matrix(0,N.sdl.18G,3); site.temp.G[,2] <- 1   #(0,1,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.G <- matrix(0,N.sdl.18G,4); Yr.temp.G[,1] <- 1   #(1,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.G <- c(X.Ht.G,PILA.18.Htmat[,1])
		
		#Elev
		X.Elev.G <- rbind(X.Elev.G,PILA.18.Elev)
		
		#SZ
		X.SZ.G <- rbind(X.SZ.G,PILA.18.SZ)
		
		#Climate		
		X.CWD.G <- c(X.CWD.G,rep(CWD[28,3],N.sdl.18G)); X.P.G <- c(X.P.G,rep(Precip[28,3],N.sdl.18G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[28,3],N.sdl.18G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[28,3],N.sdl.18G)); X.JMax.G <- c(X.JMax.G,rep(JMax[28,3],N.sdl.18G))
		
		#Genotype
		X.Geno.G <- rbind(X.Geno.G,PILA.18.Geno)
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.18.Geno2)
	}	

	if(t==4){   ###Year 4, all sites
		site.17.live.4G <- which(PILA.17.survmat[,3]==1 & PILA.17.survmat[,4]==1) #which seedlings alive in year 3&4, 2017 site	
		site.18.live.4G <- which(PILA.18.survmat[,2]==1 & PILA.18.survmat[,3]==1) #which seedlings alive in year 3&4, 2018 site
		site.19.live.4G <- which(PILA.19.survmat[,1]==1 & PILA.19.survmat[,2]==1) #which seedlings alive in year 3&4, 2019 site	
		gro.temp.1 <- PILA.17.Grow[site.17.live.4G,(t-1)]; N.sdl.17G <- length(gro.temp.1)
		gro.temp.2 <- PILA.18.Grow[site.18.live.4G,(t-2)]; N.sdl.18G <- length(gro.temp.2)
		gro.temp.3 <- PILA.19.Grow[site.19.live.4G,(t-3)]; N.sdl.19G <- length(gro.temp.3)

		###2017 site seedlings
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.1)
		
		#planting site
		site.temp.G <- matrix(0,N.sdl.17G,3); site.temp.G[,1] <- 1   #(1,0,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
				
		#Year/age
		Yr.temp.G <- matrix(0,N.sdl.17G,4); Yr.temp.G[,3] <- 1   #(0,0,1,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.G <- c(X.Ht.G,PILA.17.Htmat[site.17.live.4G,3])
		
		#Elev
		X.Elev.G <- rbind(X.Elev.G,PILA.17.Elev[site.17.live.4G,])
		
		#SZ
		X.SZ.G <- rbind(X.SZ.G,PILA.17.SZ[site.17.live.4G,])
		
		#Climate	
		X.CWD.G <- c(X.CWD.G,rep(CWD[29,2],N.sdl.17G)); X.P.G <- c(X.P.G,rep(Precip[29,2],N.sdl.17G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[29,2],N.sdl.17G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[29,2],N.sdl.17G)); X.JMax.G <- c(X.JMax.G,rep(JMax[29,2],N.sdl.17G))
		
		#Genotype
		X.Geno.G <- rbind(X.Geno.G,PILA.17.Geno[site.17.live.4G,])
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.17.Geno2[site.17.live.4G,])
		
		###2018 site seedlings
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.2)
		
		#planting site
		site.temp.G <- matrix(0,N.sdl.18G,3); site.temp.G[,2] <- 1   #(0,1,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.G <- matrix(0,N.sdl.18G,4); Yr.temp.G[,2] <- 1   #(0,1,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.G <- c(X.Ht.G,PILA.18.Htmat[site.18.live.4G,2])
		
		#Elev
		X.Elev.G <- rbind(X.Elev.G,PILA.18.Elev[site.18.live.4G,])
		
		#SZ
		X.SZ.G <- rbind(X.SZ.G,PILA.18.SZ[site.18.live.4G,])
		
		#Climate		
		X.CWD.G <- c(X.CWD.G,rep(CWD[29,3],N.sdl.18G)); X.P.G <- c(X.P.G,rep(Precip[29,3],N.sdl.18G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[29,3],N.sdl.18G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[29,3],N.sdl.18G)); X.JMax.G <- c(X.JMax.G,rep(JMax[29,3],N.sdl.18G))
		
		#Genotype
		X.Geno.G <- rbind(X.Geno.G,PILA.18.Geno[site.18.live.4G,])
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.18.Geno2[site.18.live.4G,])
		
		###2019 site seedlings
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.3)
		
		#planting site
		site.temp.G <- matrix(0,N.sdl.19G,3); site.temp.G[,3] <- 1   #(0,0,1)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.G <- matrix(0,N.sdl.19G,4); Yr.temp.G[,1] <- 1   #(1,0,0,0)
		X.Yr.G <- rbind(X.Yr.G,Yr.temp.G)
		
		#Height
		X.Ht.G <- c(X.Ht.G,PILA.19.Htmat[site.19.live.4G,1])
		
		#Elev
		X.Elev.G <- rbind(X.Elev.G,PILA.19.Elev[site.19.live.4G,])
		
		#SZ
		X.SZ.G <- rbind(X.SZ.G,PILA.19.SZ[site.19.live.4G,])
		
		#Climate		
		X.CWD.G <- c(X.CWD.G,rep(CWD[29,4],N.sdl.19G)); X.P.G <- c(X.P.G,rep(Precip[29,4],N.sdl.19G))
		X.Sno.G <- c(X.Sno.G,rep(Snow[29,4],N.sdl.19G))
		X.JMin.G <- c(X.JMin.G,rep(JMin[29,4],N.sdl.19G)); X.JMax.G <- c(X.JMax.G,rep(JMax[29,4],N.sdl.19G))
		
		#Genotype
		X.Geno.G <- rbind(X.Geno.G,PILA.19.Geno[site.19.live.4G,])
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.19.Geno2[site.19.live.4G,])
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
		Yr.temp.S <- matrix(0,N.sdl.17,3); Yr.temp.S[,3] <- 1   #(0,0,1)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.17G,4); Yr.temp.G[,4] <- 1   #(0,0,0,1)
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
		
		#Genotype
		X.Geno.S <- rbind(X.Geno.S,PILA.17.Geno[site.17.live.5,])
		X.Geno.G <- rbind(X.Geno.G,PILA.17.Geno[site.17.live.5G,])
		X.Geno.S2 <- rbind(X.Geno.S2,PILA.17.Geno2[site.17.live.5,])
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.17.Geno2[site.17.live.5G,])
		
		###2018 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.2)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.2)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.18,3); site.temp.S[,2] <- 1   #(0,1,0)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.18G,3); site.temp.G[,2] <- 1   #(0,1,0)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.18,3); Yr.temp.S[,2] <- 1   #(0,1,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.18G,4); Yr.temp.G[,3] <- 1   #(0,0,1,0)
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
		
		#Genotype
		X.Geno.S <- rbind(X.Geno.S,PILA.18.Geno[site.18.live.5,])
		X.Geno.G <- rbind(X.Geno.G,PILA.18.Geno[site.18.live.5G,])
		X.Geno.S2 <- rbind(X.Geno.S2,PILA.18.Geno2[site.18.live.5,])
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.18.Geno2[site.18.live.5G,])
		
		###2019 site seedlings
		Y.PILA.surv.all <- c(Y.PILA.surv.all,surv.temp.3)
		Y.PILA.grow.all <- c(Y.PILA.grow.all,gro.temp.3)
		
		#planting site
		site.temp.S <- matrix(0,N.sdl.19,3); site.temp.S[,3] <- 1   #(0,0,1)
		X.Sites.S <- rbind(X.Sites.S,site.temp.S)
		site.temp.G <- matrix(0,N.sdl.19G,3); site.temp.G[,3] <- 1   #(0,0,1)
		X.Sites.G <- rbind(X.Sites.G,site.temp.G)
		
		#Year/age
		Yr.temp.S <- matrix(0,N.sdl.19,3); Yr.temp.S[,1] <- 1   #(1,0,0)
		X.Yr.S <- rbind(X.Yr.S,Yr.temp.S)
		Yr.temp.G <- matrix(0,N.sdl.19G,4); Yr.temp.G[,2] <- 1   #(0,1,0,0)
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
		
		#Genotype
		X.Geno.S <- rbind(X.Geno.S,PILA.19.Geno[site.19.live.5,])
		X.Geno.G <- rbind(X.Geno.G,PILA.19.Geno[site.19.live.5G,])
		X.Geno.S2 <- rbind(X.Geno.S2,PILA.19.Geno2[site.19.live.5,])
		X.Geno.G2 <- rbind(X.Geno.G2,PILA.19.Geno2[site.19.live.5G,])
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

for(i in 1:N.grow){
	if(X.SZ.G[i,1]==1 & X.Elev.G[i,3]==1) X.SZ.Elev.G[i,1] <-1 #526_4500
    if(X.SZ.G[i,1]==1 & X.Elev.G[i,4]==1) X.SZ.Elev.G[i,2] <-1 #526_5000
    if(X.SZ.G[i,1]==1 & X.Elev.G[i,5]==1) X.SZ.Elev.G[i,3] <-1 #526_5250
    if(X.SZ.G[i,1]==1 & X.Elev.G[i,6]==1) X.SZ.Elev.G[i,4] <-1 #526_5500
    if(X.SZ.G[i,2]==1 & X.Elev.G[i,1]==1) X.SZ.Elev.G[i,5] <-1 #531_3000
    if(X.SZ.G[i,2]==1 & X.Elev.G[i,2]==1) X.SZ.Elev.G[i,6] <-1 #531_4000
    if(X.SZ.G[i,2]==1 & X.Elev.G[i,3]==1) X.SZ.Elev.G[i,7] <-1 #531_4500
    if(X.SZ.G[i,2]==1 & X.Elev.G[i,4]==1) X.SZ.Elev.G[i,8] <-1 #531_5000
    if(X.SZ.G[i,2]==1 & X.Elev.G[i,6]==1) X.SZ.Elev.G[i,9] <-1 #531_5500
    if(X.SZ.G[i,2]==1 & X.Elev.G[i,7]==1) X.SZ.Elev.G[i,10] <-1 #531_6000
    if(X.SZ.G[i,3]==1 & X.Elev.G[i,4]==1) X.SZ.Elev.G[i,11] <-1 #532_5000
    if(X.SZ.G[i,3]==1 & X.Elev.G[i,6]==1) X.SZ.Elev.G[i,12] <-1 #532_5500
    if(X.SZ.G[i,3]==1 & X.Elev.G[i,7]==1) X.SZ.Elev.G[i,13] <-1 #532_6000
    if(X.SZ.G[i,4]==1 & X.Elev.G[i,7]==1) X.SZ.Elev.G[i,14] <-1 #533_6000
    if(X.SZ.G[i,5]==1 & X.Elev.G[i,7]==1) X.SZ.Elev.G[i,15] <-1 #540_6000
	if(X.SZ.G[i,6]==1 & X.Elev.G[i,8]==1) X.SZ.Elev.G[i,16] <-1 #540_7500
	}


save.image(file="PILASdl_Geno_Data.RData")