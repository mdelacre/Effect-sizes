##############################################################################################
#####                          Distr d'échantillonnage de mu1-mu2                        #####
##############################################################################################

mudiffDE <- function(mu1,sd){

  matr <- matrix(0,1000000,3)

  for (i in 1:1000000){
    y1=rnorm(20,mu1,sd=sd)
    y2=rnorm(20,0,sd=sd)
    matr[i,1]=mean(y1)-mean(y2)  # diff mean
    matr[i,2]=sqrt((var(y1)+var(y2))/2) # pooled sd
    matr[i,3]=(mean(y1)-mean(y2))/sqrt((var(y1)+var(y2))/2) # pooled sd
  }
return(matr)  
}

A=mudiffDE(mu1=1,sd=1)
B=mudiffDE(mu1=5,sd=1)
C=mudiffDE(mu1=14,sd=1)

D=mudiffDE(mu1=20,sd=1)

# mu1-mu2
## La moyenne de cette distribution = mu1-mu2
mean(A[,1]) # sqrt(1/20+1/20) = .316
mean(B[,1]) # sqrt(1/20+1/20) = .316
mean(C[,1]) # sqrt(1/20+1/20) = .316
## mais la variance de cette distribution ne dépend PAS de mu1-mu2
sd(A[,1]) # sqrt(1/20+1/20) = .316
sd(B[,1]) # sqrt(1/20+1/20) = .316
sd(C[,1]) # sqrt(1/20+1/20) = .316

# pooled sd
## La moyenne de cette distribution ne dépend pas de mu1-mu2
mean(A[,2]) 
mean(B[,2]) 
mean(C[,2]) 
# la variance de cette distribution ne dépend PAS non plus de mu1-mu2
sd(A[,2]) 
sd(B[,2]) 
sd(C[,2]) 

# Cohen's d
## La moyenne de cette distribution ne dépend pas de mu1-mu2
mean(A[,3]) 
mean(B[,3]) 
mean(C[,3]) 
# la variance de cette distribution ne dépend PAS non plus de mu1-mu2
sd(A[,3]) 
sd(B[,3])



Cohen <- function(sd,nSims){
  
  SD <- rep(0,50)

  for (k in 1:50){

    cohen=rep(0,nSims)
    for (i in 1:nSims){
      y1=rnorm(20,k,sd=sd)
      y2=rnorm(20,0,sd=sd)
      cohen[i]=(mean(y1)-mean(y2))/sqrt((var(y1)+var(y2))/2) # pooled sd
    }
    
  SD[k]=sd(cohen)
    
  }
  return(SD)  
}

A=Cohen(sd=1,nSims=100000)
B=Cohen(sd=5,nSims=100000)

data=data.frame(mudiff=1:50,A)
linearMod <- lm(A ~ mudiff, data=data)
plot(1:50,A) # quand sd=1, relation linéaire entre meandiff et variance de la distribution d'échantillonnage. 
abline(linearMod[[1]][1], linearMod[[1]][2])
abline(v=1)

data=data.frame(mudiff=1:50,B)
linearMod <- lm(B ~ mudiff, data=data)
plot(1:50,B) # quand sd=1, relation linéaire entre meandiff et variance de la distribution d'échantillonnage. 
abline(linearMod[[1]][1], linearMod[[1]][2])
abline(v=1)




write.table(A,"testA.txt",sep=";",dec=",")
lines(1:50,B) # quand sd != 1, relation non linéaire


plot(1:50,B,xlim=c(1,50))


##############################################################################################
#####                                 REMARQUE GLASS DS                                  #####
##############################################################################################

##### REMARQUE RELATIVE AU FAIT QUE LE BIAIS A L AIR PIRE QD ON PREND SD2 DANS GLASS DS  #####
##############################################################################################

## La distribution d'échantillonnage de l'écart-type
# Plus sigma est grand, plus l'écart-type de la distribution d'échantillonnage est grand
# Mais en termes proportionnel (par rapport à la taille de sigma), ça n'est QUE fonction du n!

test=function(sd,nSims,n){
  
  Moy<-rep(0,nSims)
  Ecartype<-rep(0,nSims)
  Vari<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y=rnorm(n,sd=sd)  

    Moy[i] <- mean(y)
    Ecartype[i] <- sd(y)
    Vari[i] <- var(y)
    
  }
  
return(c(sd(Moy),sd(Ecartype),sd(Vari),var(Moy),var(Ecartype),var(Vari)))
}

A=test(sd=.1,nSims=1000000,n=100)
B=test(sd=1,nSims=1000000,n=100)
C=test(sd=10,nSims=1000000,n=100)
D=test(sd=100,nSims=1000000,n=100)

##### Distribution d'échantillonnage de la moyenne #####

## SD de la distribution
plot(1:4,c(A[1],B[1],C[1],D[1]),pch=19,xaxt="n",bty="n",xlab="population SD",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[1],B[1],C[1],D[1]),col="black",lty=1)
text(c(A[1]+.1,B[1]-.1,C[1]+.1,D[1]+.1), labels=round(c(A[1],B[1],C[1],D[1]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.1,1,10,100))

points(1:4,c(A[1]/.1,B[1]/1,C[1]/10,D[1]/100),col="darkgrey",pch=19)
lines(1:4,c(A[1]/.1,B[1]/1,C[1]/10,D[1]/100),col="darkgrey",lty=2)
text(c(A[1]/.1+.1,B[1]/1+.1,C[1]/10+.1,D[1]/100+.1), labels=round(c(A[1]/.1,B[1]/1,C[1]/10,D[1]/100),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.1,1,10,100))
legend("top",c("SD(sdest)","SD(sdest)/sd"),col=c("black","darkgrey"),pch=19)

A[1]/.1
B[1]/1
C[1]/10
D[1]/100

## Variance de la distribution
plot(1:4,c(A[4],B[4],C[4],D[4]),pch=19,xaxt="n",bty="n",xlab="population SD",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[4],B[4],C[4],D[4]),col="black",lty=1)
text(c(A[4]+.1,B[4]-.1,C[4]+.1,D[4]+.1), labels=round(c(A[4],B[4],C[4],D[4]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,10,100,10000))

points(1:4,c(A[4]/.01,B[4]/10,C[4]/100,D[4]/10000),col="darkgrey",pch=19)
lines(1:4,c(A[4]/.01,B[4]/10,C[4]/100,D[4]/10000),col="darkgrey",lty=2)
text(c(A[4]/.01+.1,B[4]/10+.1,C[4]/100+.1,D[4]/10000+.1), labels=round(c(A[4]/.01,B[4]/10,C[4]/100,D[4]/10000),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,10,100,10000))
legend("top",c("sd(varest)","sd(varest)/var"),col=c("black","darkgrey"),pch=19)

A[5]/.1^2  # A[5]/.1^2-(A[2]/.1)^2=0
B[5]/1^2 # B[5]/1^2-(B[2]/1)^2=0
C[5]/10^2 # C[5]/10^2-(C[2]/10)^2=0
D[5]/100^2 # (D[2]/100)^2-D[5]/100^2=0

##### Distribution d'échantillonnage de l'écart-type #####

## SD de la distribution
plot(1:4,c(A[2],B[2],C[2],D[2]),pch=19,xaxt="n",bty="n",xlab="population SD",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[2],B[2],C[2],D[2]),col="black",lty=1)
text(c(A[2]+.1,B[2]-.1,C[2]+.1,D[2]+.1), labels=round(c(A[2],B[2],C[2],D[2]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.1,1,10,100))

points(1:4,c(A[2]/.1,B[2]/1,C[2]/10,D[2]/100),col="darkgrey",pch=19)
lines(1:4,c(A[2]/.1,B[2]/1,C[2]/10,D[2]/100),col="darkgrey",lty=2)
text(c(A[2]/.1+.1,B[2]/1+.1,C[2]/10+.1,D[2]/100+.1), labels=round(c(A[2]/.1,B[2]/1,C[2]/10,D[2]/100),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.1,1,10,100))
legend("top",c("SD(sdest)","SD(sdest)/sd"),col=c("black","darkgrey"),pch=19)

A[2]/.1
B[2]/1
C[2]/10
D[2]/100

## Variance de la distribution
plot(1:4,c(A[5],B[5],C[5],D[5]),pch=19,xaxt="n",bty="n",xlab="population SD",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[5],B[5],C[5],D[5]),col="black",lty=1)
text(c(A[5]+.1,B[5]-.1,C[5]+.1,D[5]+.1), labels=round(c(A[5],B[5],C[5],D[5]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,10,100,10000))

points(1:4,c(A[5]/.01,B[5]/10,C[5]/100,D[5]/10000),col="darkgrey",pch=19)
lines(1:4,c(A[5]/.01,B[5]/10,C[5]/100,D[5]/10000),col="darkgrey",lty=2)
text(c(A[5]/.01+.1,B[5]/10+.1,C[5]/100+.1,D[5]/10000+.1), labels=round(c(A[5]/.01,B[5]/10,C[5]/100,D[5]/10000),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,10,100,10000))
legend("top",c("sd(varest)","sd(varest)/var"),col=c("black","darkgrey"),pch=19)

A[5]/.1^2  # A[5]/.1^2-(A[2]/.1)^2=0
B[5]/1^2 # B[5]/1^2-(B[2]/1)^2=0
C[5]/10^2 # C[5]/10^2-(C[2]/10)^2=0
D[5]/100^2 # (D[2]/100)^2-D[5]/100^2=0

##### Distribution d'échantillonnage de la variance #####

## SD de la distribution
plot(1:4,c(A[3],B[3],C[3],D[3]),pch=19,xaxt="n",bty="n",xlab="population SD",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[3],B[3],C[3],D[3]),col="black",lty=1)
text(c(A[3]+.1,B[3]-.1,C[3]+.1,D[3]+.1), labels=round(c(A[3],B[3],C[3],D[3]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,1,100,10000))

points(1:4,c(A[3]/.01,B[3]/1,C[3]/100,D[3]/10000),col="darkgrey",pch=19)
lines(1:4,c(A[3]/.01,B[3]/1,C[3]/100,D[3]/10000),col="darkgrey",lty=2)
text(c(A[3]/.01+100,B[3]/1+100,C[3]/100+100,D[3]/10000+100), labels=round(c(A[3]/.01,B[3]/1,C[3]/100,D[3]/10000),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.1,1,10,100))
legend("top",c("var(sdest)","var(sdest)/var"),col=c("black","darkgrey"),pch=19)

A[3]/.1^2
B[3]/1
C[3]/10^2
D[3]/100^2

## Variance de la distribution
plot(1:4,c(A[6],B[6],C[6],D[6]),pch=19,xaxt="n",bty="n",xlab="population SD",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[6],B[6],C[6],D[6]),col="black",lty=1)
text(c(A[6]+.1,B[6]-.1,C[6]+.1,D[6]+.1), labels=round(c(A[6],B[6],C[6],D[6]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,10,100,10000))

points(1:4,c(A[6]/.01^2,B[6]/10^2,C[6]/100^2,D[6]/10000^2),col="darkgrey",pch=19)
lines(1:4,c(A[6]/.01^2,B[6]/10^2,C[6]/100^2,D[6]/10000^2),col="darkgrey",lty=2)
text(c(A[6]/.01^2+.1,B[6]/10^2+.1,C[6]/100^2+.1,D[6]/10000^2+.1), labels=round(c(A[6]/.01,B[6]/10,C[6]/100,D[6]/10000),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.01,10,100,10000))
legend("top",c("sd(varest)","sd(varest)/var"),col=c("black","darkgrey"),pch=19)

A[6]/.1^4 # A[6]/.1^4-(A[3]/.1^2)^2=0
B[6]/1^4 # B[6]/1^4 -(B[3]/1)^2=0
C[6]/10^4 # C[6]/10^4-(C[3]/10^2)^2=0
D[6]/100^4 # D[6]/100^4-(D[3]/100^2)^2=0

##### POURQUOI LE BIAIS RELATIF EST PIRE QD ON PREND SD1 QD G1 < 0, ET SD2  QUAND G1 > 0? ####
##############################################################################################

# Application du script "Simulations" (auquel j'ai juste modifié le "setwd" pour stocker les résultats dans le dossier "Réflexion interprétation simulation")
get_simu(nSims=10000,n1=20,n2=20,sd1=2,sd2=2,m1=0,m2=1,skew=-6.32,kurt=95.75)
get_simu(nSims=10000,n1=20,n2=20,sd1=2,sd2=2,m1=1,m2=0,skew=-6.32,kurt=95.75)  
get_simu(nSims=10000,n1=20,n2=20,sd1=2,sd2=2,m1=0,m2=1,skew=6.32,kurt=95.75)
get_simu(nSims=10000,n1=20,n2=20,sd1=2,sd2=2,m1=1,m2=0,skew=6.32,kurt=95.75)

Folder="C:/Users/Marie/Documents/Github_projects/Effect-sizes/Réflexion interprétation simulations"
File=readRDS(list.files(Folder,pattern=".*G1=")[i])
param <- str_extract_all(list.files(Folder,pattern=".*G1=")[i], "[[:digit:]]+\\.*[[:digit:]]*")
n1 <- as.numeric(param[[1]][5])
n2 <- as.numeric(param[[1]][6])
m1 <- as.numeric(param[[1]][7])
m2 <- as.numeric(param[[1]][8])
sd1 <- as.numeric(param[[1]][9])
sd2 <- as.numeric(param[[1]][10])
glass_delta1 <- (m1-m2)/sd1
glass_delta2 <- (m1-m2)/sd2

bias_glass1 <- mean(File[,11]) - glass_delta1 
bias_glass2 <- mean(File[,12]) - glass_delta2

relbias_glass1 <- (mean(File[,11]) - glass_delta1)/glass_delta1 
relbias_glass2 <- (mean(File[,12]) - glass_delta2)/glass_delta2
param

bias_glass1  
bias_glass2

relbias_glass1  
relbias_glass2

# On voit très nettement qu'en cas d'asymétrie négative, toujours pire de prendre le sd du groupe ayant la moyenne la plus élevée comme standardiseur.
# Pourtant, même écart-type, même taille, même G1, même G2 (seule la moyenne change).

## La distribution d'échantillonnage de l'écart-type quand asymétrie négative?

# Ecart-type de la distribution?
test=function(sd,nSims=100000,m1,m2,n,skew1,skew2,kurt=95.75){
  
  Ecartype1<-rep(0,nSims)
  Ecartype2<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew1*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew2*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    Ecartype1[i] <- sd(y1)
    Ecartype2[i] <- sd(y2)

  }
  
  plot(density(Ecartype1))
  lines(density(Ecartype2),col="red",lty=2)
  #return(c(sd(Ecartype1),sd(Ecartype2)))
}

A=test(sd=1,m1=0,m2=0,skew1=-6.32,skew2=6.32,n=20)


test=function(sd,nSims=100000,m1,m2,n,skew1,skew2,kurt=95.75){
  
  diffmoyneg<-rep(0,nSims)
  diffmoypos<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew1*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew1*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  

    y3=rpearson(n,moments=c(m1,sd^2,skewness=skew2*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y4=rpearson(n,moments=c(m2,sd^2,skewness=skew2*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    diffmoyneg[i] <- mean(y1)-mean(y2)
    diffmoypos[i] <- mean(y3)-mean(y4)
    
  }
  
  plot(density(diffmoyneg))
  lines(density(diffmoypos),col="red",lty=2)
  #return(c(sd(Ecartype1),sd(Ecartype2)))
}

A=test(sd=1,m1=0,m2=0,skew1=-6.32,skew2=6.32,n=20)
B=test(sd=1,m1=1,m2=0,skew1=-6.32,skew2=6.32,n=20)





test=function(sd,nSims=100000,m1,m2,n,skew,kurt=95.75){
  
  Glass1<-rep(0,nSims)
  Glass2<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    Glass1[i] <- (mean(y1)-mean(y2))/sd(y1)
    Glass2[i] <- (mean(y1)-mean(y2))/sd(y2)
    
  }
  
  plot(density(Glass1))
  lines(density(Glass2),col="red",lty=2)
  #return(c(sd(Ecartype1),sd(Ecartype2)))
}

A=test(sd=1,m1=1,m2=0,skew=6.32,n=20)
B=test(sd=1,m1=1,m2=0,skew=-6.32,n=20)



A=test(sd=1,m1=0,m2=0,skew=-6.32,n=20)
B=test(sd=1,m1=0,m2=0,skew=6.32,n=20)

B=test(sd=1,m1=0,m2=0,skew=6.32,n=20)




A=test(sd=1,m1=0,m2=0,skew=-6.32,n=20)
B=test(sd=1,m1=0,m2=0,skew=6.32,n=20)



test=function(sd,nSims=100000,m1,m2,m3,n,skew,kurt=95.75){
  
  Ecartype1<-rep(0,nSims)
  Ecartype2<-rep(0,nSims)
  Ecartype3<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y3=rpearson(n,moments=c(m3,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    Ecartype1[i] <- sd(y1)
    Ecartype2[i] <- sd(y2)
    Ecartype3[i] <- sd(y3)
    
  }
  
  return(c(sd(Ecartype1),sd(Ecartype2),sd(Ecartype3)))
}

test(sd=5,m1=0,m2=0,m3=0,n=20,skew=6.32)

B # Ne semble absolument pas dépendre de la moyenne!

# Moyenne de la distribution?
test=function(sd,nSims=100000,m1,m2,n,skew,kurt=95.75){
  
  Ecartype1<-rep(0,nSims)
  Ecartype2<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    Ecartype1[i] <- sd(y1)
    Ecartype2[i] <- sd(y2)
    
  }
  
  return(c(mean(Ecartype1),mean(Ecartype2)))
}

A=test(sd=1,m1=5,m2=0,skew=-6.32,n=20)
B=test(sd=1,m1=0,m2=5,skew=-6.32,n=20)

A
B #  ne semble pas non plus dépendre de la moyenne...

# Allure générale de la distribution?
test=function(sd,nSims=100000,m1,m2,n,skew,kurt=95.75){
  
  Ecartype1<-rep(0,nSims)
  Ecartype2<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    Ecartype1[i] <- sd(y1)
    Ecartype2[i] <- sd(y2)
    
  }

  par(mfrow=c(2,1))    
  plot(0,xlim=c(0,20),col="white")
  lines(density(Ecartype1),xlim=c(0,20),col="pink")
  lines(density(Ecartype2),xlim=c(0,20),lty=2,col="black")
  
  return(c(min(Ecartype1),min(Ecartype2),max(Ecartype1),max(Ecartype2),skewness(Ecartype1),skewness(Ecartype2)))
  
}

A=test(sd=1,m1=5,m2=0,skew=-6.32,n=20)
B=test(sd=1,m1=0,m2=5,skew=-6.32,n=20)
A
B

test=function(sd,nSims=100000,m1,m2,n,skew,kurt=95.75){
  
  Glass1<-rep(0,nSims)
  Glass2<-rep(0,nSims)
  
  for (i in 1:nSims){
    
    y1=rpearson(n,moments=c(m1,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    y2=rpearson(n,moments=c(m2,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
    
    
    Glass1[i] <- (mean(y1)-mean(y2))/sd(y1)
    Glass2[i] <- (mean(y1)-mean(y2))/sd(y2)
    
  }
  
  return(c(mean(Glass1),mean(Glass2)))
}

A=test(sd=1,m1=1,m2=0,skew=-6.32,n=20) # Biais plus important quand on prend sd1 comme standardiseur
A
A-((1-0)/1)
B=test(sd=1,m1=0,m2=1,skew=-6.32,n=20) # Biais plus important quand on prend sd2 comme standardiseur
B
B-((0-1)/1)

plot(density(y1))
plot(density(y2))

library(moments)
skewness(y1)
skewness(y2)
##### Distribution d'échantillonnage de l'écart-type #####

## SD de la distribution
plot(1:4,c(A[2],B[2],C[2],D[2]),pch=19,xaxt="n",bty="n",xlab="sample mean",ylab="SD",cex.lab=1.5,main="distribution d'échantillonnage de l'écart-type")
lines(1:4,c(A[2],B[2],C[2],D[2]),col="black",lty=1)
text(c(A[2]+.1,B[2]-.1,C[2]+.1,D[2]+.1), labels=round(c(A[2],B[2],C[2],D[2]),2), col="black",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(-5,0,5,10))

points(1:4,c(A[2]/.1,B[2]/1,C[2]/10,D[2]/100),col="darkgrey",pch=19)
lines(1:4,c(A[2]/.1,B[2]/1,C[2]/10,D[2]/100),col="darkgrey",lty=2)
text(c(A[2]/.1+.1,B[2]/1+.1,C[2]/10+.1,D[2]/100+.1), labels=round(c(A[2]/.1,B[2]/1,C[2]/10,D[2]/100),2), col="darkgrey",cex= .8, offset = 10)
axis(1,at=1:4,labels=c(.1,1,10,100))
legend("top",c("SD(sdest)","SD(sdest)/sd"),col=c("black","darkgrey"),pch=19)

A[2]/.1
B[2]/1
C[2]/10
D[2]/100

#####                 POURQUOI GLASS CRAINT TOUJOURS PLUS QUE LES AUTRES?                 ####
##############################################################################################

# Car l'écart-type est estimé sur base d'un n plus petit.

test=function(sd,nSims,n){
  
  Ecartype<-rep(0,nSims)

  for (i in 1:nSims){
    
    y=rnorm(n,sd=sd)  
    Ecartype[i] <- sd(y)

  }

  return(c(mean(Ecartype),sd(Ecartype),skewness(Ecartype),kurtosis(Ecartype)))
}

A=test(sd=5,nSims=1000000,n=20)
B=test(sd=5,nSims=1000000,n=40)
C=test(sd=5,nSims=1000000,n=80)

round(A,2)
round(B,2)
round(C,2)

# On voit que quand n augmente, l'écart-type de la distribution d'échantillonnage de SD diminue, ainsi que l'asymétrie (on se rapproche de plus en plus d'un truc symétrique)





#####           QD SIGMA1=SIGMA2 et VAR INEGALE. QUID DE LA VARIANCE DE SHIEH?            ####
##############################################################################################

