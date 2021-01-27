####################################################################################################
### INTERACTION BETWEEN SKEWNESS (RIGHT VS. LEFT) AND WORST CHOICE BETWEEN SD1 AND SD2 FOR GLASS ###
####################################################################################################

library(PearsonDS)
library(moments)

test=function(sd,nSims=100000,m1,m2,n1,n2,skew,kurt=95.75){
   
   glasspos1<-rep(0,nSims)
   glasspos2<-rep(0,nSims)
   
   sd1<-rep(0,nSims)
   sd2<-rep(0,nSims)
   meandiff<-rep(0,nSims)
   mean1<-rep(0,nSims)
   mean2<-rep(0,nSims)
   
   y=rpearson(1000000,moments=c(m1,sd^2,skewness=skew*(n1-2)/sqrt(n1*(n1-1)),kurtosis=(kurt*(n1-2)*(n1-3)-6*(n1-1))/(n1^2-1)+3))  
   
   for (i in 1:nSims){

      y2 <- sample(y,size=n1,replace=TRUE)
      y1 <- sample(y,size=n2,replace=TRUE)
      sd1[i] <- sd(y1)
      sd2[i] <- sd(y2)
      meandiff[i] <- mean(y1)-mean(y2)
      mean1[i] <- mean(y1)
      mean2[i] <- mean(y2)
      
      glasspos1[i] <- (mean(y1)-mean(y2))/sd(y1)
      glasspos2[i] <- (mean(y1)-mean(y2))/sd(y2)
      
   }
   
   plot(density(glasspos1),col="blue",lty=2)
   lines(density(glasspos2),col="green",lty=2)
   
   return(cbind(sd1,sd2,meandiff,glasspos1,glasspos2,mean1,mean2))
}

######################################################################################################################
##              When distributions are not skewed, whatever one choses sd1 or sd2, results are identical            ##
######################################################################################################################

# When n1=n2

n1 <- 20
n2 <- 20
N <- n1+n2
B=test(sd=1,m1=0,m2=0,skew=0,n1=n1,n2=n2) # symmetric distribution

cor(B[,1],B[,3]) # no correlation between sd1 and meandiff
cor(B[,2],B[,3]) # no correlation between sd2 and meandiff
cor(B[,1],B[,2]) # no correlation between sd1 and sd2

stdz_shieh <- sqrt(B[,1]^2/(n1/N)+B[,2]^2/(n2/N))
stdz_cohendprime <- sqrt((B[,1]^2+B[,2]^2)/2)

cor(stdz_shieh,B[,3]) # no correlation between stdz_shieh and the mean difference
cor(stdz_cohendprime,B[,3]) # no correlation between stdz_cohendprime and the mean difference

mean(B[,4]) # = 0
mean(B[,5]) # = 0

cor(B[,3],B[,6]) # positive correlation between meandiff and m1 (the larger m1, the larger meandiff)
cor(B[,3],B[,7]) # negative correlation between meandiff and m2 (the lower m2, the larger meandiff) --> of course, as we compute m1-m2
cor(B[,6],B[,7]) # no correlation between m1 and m2

cor(B[,1],B[,6])  # no correlation between m and sd
cor(B[,2],B[,7])  # because both sampling distributions of m and sd are symmetric

table(B[,4]<B[,5]) # true, about half the time

table(B[,4]<0)
table(B[,5]<0)

median(B[,4]) # median = 0
median(B[,5]) 

# When n1!=n2

n1 <- 20
n2 <- 40
N <- n1+n2

C=test(sd=1,m1=0,m2=0,skew=0,n1=n1,n2=n2) # symmetric distribution

cor(C[,1],C[,3]) # no correlation between sd1 and meandiff
cor(C[,2],C[,3]) # no correlation between sd2 and meandiff
cor(C[,1],C[,2]) # no correlation between sd1 and sd2

stdz_shieh <- sqrt(C[,1]^2/(n1/N)+C[,2]^2/(n2/N))
stdz_cohendprime <- sqrt((C[,1]^2+C[,2]^2)/2)

cor(stdz_shieh,C[,3]) # no correlation between stdz_shieh and the mean difference
cor(stdz_cohendprime,C[,3]) # no correlation between stdz_cohendprime and the mean difference

######################################################################################################################
##                When distributions are skewed, results vary as a function of the chosen standardizer              ##
######################################################################################################################

# When n1=n2

n1 <- 20
n2 <- 20
N <- n1+n2
D=test(sd=1,m1=0,m2=0,skew=6.32,n1=n1,n2=n2) # symmetric distribution

stdz_cohen <- sqrt(((n1-1)*D[,1]^2+(n2-1)*D[,2]^2)/(n1+n2-2))
stdz_shieh <- sqrt(D[,1]^2/(n1/N)+D[,2]^2/(n2/N))
stdz_cohendprime <- sqrt((D[,1]^2+D[,2]^2)/2)

plot(D[,1],D[,3])
cor(D[,1],D[,3]) # .628; positive correlation between sd1 and meandiff
cor(D[,2],D[,3]) # -.628; negative correlation between sd2 and meandiff
cor(D[,1],D[,2]) # no correlation between sd1 and sd2

cor(stdz_cohen,D[,3]) # no correlation between stdz_shieh and the mean difference
cor(stdz_shieh,D[,3]) # no correlation between stdz_shieh and the mean difference
cor(stdz_cohendprime,D[,3]) # no correlation between stdz_cohendprime and the mean difference

# When n1!=n2

n1 <- 20
n2 <- 100
N <- n1+n2

E=test(sd=1,m1=0,m2=0,skew=6.32,n1=n1,n2=n2) # symmetric distribution

cor(E[,1],E[,3]) # positive correlation between sd1 and meandiff
cor(E[,2],E[,3]) # negative correlation between sd2 and meandiff
# Note: the largest correlation is the one between meandiff and the sd extracted from the largest group.
cor(E[,1],E[,2]) # no correlation between sd1 and sd2

stdz_cohen <- sqrt(((n1-1)*E[,1]^2+(n2-1)*E[,2]^2)/(n1+n2-2))
stdz_shieh <- sqrt(E[,1]^2/(n1/N)+E[,2]^2/(n2/N))
stdz_cohendprime <- sqrt((E[,1]^2+E[,2]^2)/2)

cor(stdz_cohen,E[,3]) # correlation: -.735
cor(stdz_shieh,E[,3]) # no correlation
cor(stdz_cohendprime,E[,3]) # correlation: -.492

##################################################################
##                   studying right-skewed case                 ##
##################################################################
# Note: when y1 and y2 are left-skewed, we find opposite results

##### 1) m1=m2

A=test(sd=1,m1=0,m2=0,skew=6.32,n=20)
skewness(A[,4])
skewness(A[,5])

cor(A[,1],A[,3]) # positive correlation between meandiff and sd1 
cor(A[,2],A[,3]) # negative correlation between meandiff and sd2
cor(A[,1],A[,2]) # no correlation between sd1 and sd2

mean(A[,4]) # below 0 (i.e. below the real delta)
mean(A[,5]) # above 0 (i.e. above the real delta)

cor(A[,3],A[,6]) # positive correlation between meandiff and m1 (the larger m1, the larger meandiff)
cor(A[,3],A[,7]) # negative correlation between meandiff and m2 (the lower m2, the larger meandiff) --> of course, as we compute m1-m2
cor(A[,6],A[,7]) # no correlation between m1 and m2

cor(A[,1],A[,6])  # positive correlation between m and sd (always true when y1 and y2 are right-skewed, because the sampling distribution of the mean 
cor(A[,2],A[,7])  # will be right-skewed - and as in all configurations, the sampling distribution of sd is of course also right-skewed).

table(A[,4]<A[,5]) # true in mosts cases (i.e. 86%)

table(A[,4]<0)
table(A[,5]<0)

median(A[,4]) # about 50% of glass's measures are > 0, and other half is below 0 (therefore, 0 is the median)
median(A[,5]) # but the sampling distribution of glass is either positive, or negative, median != mean

# These results are not dependent on the 'rpearson' function (I've reproduced the test with rchisq and rsnorm, and I found similar results)

##### 2) m1>m2

C=test(sd=1,m1=0,m2=0,skew=6.32,n=20)
skewness(C[,4])
skewness(C[,5])

cor(C[,1],C[,3]) # positive correlation between meandiff and sd1 
cor(C[,2],C[,3]) # negative correlation between meandiff and sd2
cor(C[,1],C[,2]) # no correlation between sd1 and sd2

mean(C[,4]) # below 0 (i.e. below the real delta)
mean(C[,5]) # above 0 (i.e. above the real delta)

cor(C[,3],C[,6]) # positive correlation between meandiff and m1 (the larger m1, the larger meandiff)
cor(C[,3],C[,7]) # negative correlation between meandiff and m2 (the lower m2, the larger meandiff) --> of course, as we compute m1-m2
cor(C[,6],C[,7]) # no correlation between m1 and m2

cor(C[,1],C[,6])  # positive correlation between m and sd (always true when y1 and y2 are right-skewed, because the sampling distribution of the mean 
cor(C[,2],C[,7])  # will be right-skewed - and as in all configurations, the sampling distribution of sd is of course also right-skewed).

table(C[,4]<C[,5]) # true in mosts cases (i.e. 86%)

table(C[,4]<0)
table(C[,5]<0)

median(C[,4]) # about 50% of glass's measures are > 0, and other half is below 0 (therefore, 0 is the median)
median(C[,5]) # but the sampling distribution of glass is either positive, or negative, median != mean


######################################################################################################################
##                              When distributions are skewed, n1=n2 but sigma1!=sigma2                             ##
######################################################################################################################

library(PearsonDS)
library(moments)

test=function(sd1,sd2,nSims=100000,m1,m2,n,skew,kurt=95.75){
   
   glasspos1<-rep(0,nSims)
   glasspos2<-rep(0,nSims)
   
   SD1<-rep(0,nSims)
   SD2<-rep(0,nSims)
   meandiff<-rep(0,nSims)
   mean1<-rep(0,nSims)
   mean2<-rep(0,nSims)
   
   Y1=rpearson(1000000,moments=c(m1,sd1^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
   Y2=rpearson(1000000,moments=c(m1,sd2^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
   
   for (i in 1:nSims){
      
      y1 <- sample(Y1,size=n,replace=TRUE)
      y2 <- sample(Y2,size=n,replace=TRUE)
      SD1[i] <- sd(y1)
      SD2[i] <- sd(y2)
      meandiff[i] <- mean(y1)-mean(y2)
      mean1[i] <- mean(y1)
      mean2[i] <- mean(y2)
      
      glasspos1[i] <- (mean(y1)-mean(y2))/sd(y1)
      glasspos2[i] <- (mean(y1)-mean(y2))/sd(y2)
      
   }
   
   plot(density(glasspos1),col="blue",lty=2)
   lines(density(glasspos2),col="green",lty=2)
   
   return(cbind(SD1,SD2,meandiff,glasspos1,glasspos2,mean1,mean2))
}

F=test(sd1=1,sd2=3,m1=0,m2=0,skew=6.32,n=20) # symmetric distribution

cor(F[,1],F[,3]) # positive correlation between sd1 and meandiff
cor(F[,2],F[,3]) # negative correlation between sd2 and meandiff
# Note: the largest correlation is the one between meandiff and the sd extracted from the largest group.
cor(F[,1],F[,2]) # no correlation between sd1 and sd2

stdz_cohen <- sqrt(((n1-1)*F[,1]^2+(n2-1)*F[,2]^2)/(n1+n2-2))
stdz_shieh <- sqrt(F[,1]^2/(n1/N)+F[,2]^2/(n2/N))
stdz_cohendprime <- sqrt((F[,1]^2+F[,2]^2)/2)

cor(stdz_cohen,F[,3]) # correlation: -.822
cor(stdz_cohendprime,F[,3]) # correlation: -.768
cor(stdz_shieh,F[,3]) # -.54
# Normal que Shieh soit un peu meilleur que Cohendprime, étant donné que corrélation moins forte entre meandiff et le