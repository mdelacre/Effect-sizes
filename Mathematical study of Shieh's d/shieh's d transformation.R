# Shieh's d depends on the n-ratio.
# For the same amount of differences between two means, same SD and sd-ratio
# Shieh's d will vary as a function of the n-ratio.

# Case 1: when sd1 = sd2

n=500
shieh<-function(n1,mu1=1,mu2=0,sd1=2,sd2=2)
{
n1 <- n1
n2 <- n-n1
#pooled_sd<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2))
q1<-n1/n
q2<-n2/n
#cohen_d<-(mu1-mu2)/pooled_sd
shieh_d<-(mu1-mu2)/sqrt(sd1^2/q1+sd2^2/q2)
return(cbind(n1,n2,nratio=n1/n2,shieh=shieh_d))
}

# Let's make vary n-ratio, keeping N constant.

rep <-data.frame(t(sapply(seq(1,n-1,1),shieh)))
colnames(rep)=c("n1","n2","nratio","shieh")
rep

rep$nratio[rep$shieh==max(rep$shieh)] # As one can see, Shieh's d achieves its maximum value when n1 = n2. It's always true (not only in this example)

# when n1=n2, shieh's d = .25
## Note: shieh's d = cohen's d/2 (whatever sd1=sd2 or not)
# see the word file "Relation between Shieh and Cohen when n1=n2"

##############################################################################################################################################
###  QUESTION: WHAT IS THE RELATION BETWEEN the n-ratio and the multiplier required in order to obtain same shieh's d value as when n1=n1  ###
##############################################################################################################################################

## I would like to find a mathematical transformation of Shieh's d measure, in order to be able to answer the
## following question: "what would be the shieh's d effect size is design was balanced?"
## The correction would work for any cases where there is homoscedasticity

# It seems very improbable to find a transformation that fits all configurations of sd, nratio, N... but it would be great to do so!


## WHAT IS THE RELATION BETWEEN SHIEH's D and the n-ratio?
plot(rep$nratio,rep$shieh,ylim=c(0,.3),pch=19,cex=.3,xlab="n-ratio",ylab="shieh's d")
abline(h=max(rep$shieh),lty=2,col="lightgrey")
abline(v=rep$nratio[rep$shieh==max(rep$shieh)],lty=2,col="lightgrey")
# Zoom on the plot when n-ratio < 1
plot(rep$nratio[rep$nratio < 1],rep$shieh[rep$nratio < 1],ylim=c(0,.3),xlim=c(min(rep$nratio),1),pch=19,cex=.3,xlab="n-ratio",ylab="shieh's d")
obj<-loess(rep$shieh[rep$nratio < 1] ~ rep$nratio[rep$nratio < 1])
lines(obj,lty=2,col="red")
# Zoom on the plot when n-ratio > 1
plot(rep$nratio[rep$nratio > 1],rep$shieh[rep$nratio > 1],ylim=c(0,.3),xlim=c(1,max(rep$nratio)),pch=19,cex=.3,xlab="n-ratio",ylab="shieh's d")
obj<-loess(rep$shieh[rep$nratio > 1] ~ rep$nratio[rep$nratio > 1])
lines(obj,lty=2,col="red")

# Compute the required correction in order to achieve same value as when n1=n2:
coeff <-max(rep$shieh)/rep$shieh # multiplier in order to obtain the same value than when n1=n2...
# Exemple: rep[180,]    when n1=180 and n2=20, one should multiply shieh's d by coeff[180] = 1.666667 in order to have .25

# Is there a simple relation between the n-ratio and the correction to perform in shieh's d?
plot(rep$nratio,coeff,pch=19,cex=.3,xlab="n-ratio",ylab="multiplier")
abline(h=max(rep$shieh),lty=2,col="lightgrey")
abline(v=rep$nratio[rep$shieh==max(rep$shieh)],lty=2,col="lightgrey")
# Zoom on the plot when n-ratio < 1
plot(rep$nratio[rep$nratio < 1],coeff[rep$nratio < 1],xlim=c(min(rep$nratio),1),pch=19,cex=.3,xlab="n-ratio",ylab="multiplier")
obj<-loess(coeff[rep$nratio < 1] ~ rep$nratio[rep$nratio < 1])
lines(obj,lty=2,col="red")
# Zoom on the plot when n-ratio > 1
plot(rep$nratio[rep$nratio > 1],coeff[rep$nratio > 1],xlim=c(1,max(rep$nratio)),pch=19,cex=.3,xlab="n-ratio",ylab="multiplier")
obj<-loess(coeff[rep$nratio > 1] ~ rep$nratio[rep$nratio > 1])
lines(obj,lty=2,col="red")

write.table(cbind(rep,multiplier=coeff),"multipliervsnratio.txt",sep=";",dec=",")
# trying to figure out a law underlying this relation!
