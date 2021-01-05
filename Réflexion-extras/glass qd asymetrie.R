
######################################################################################################################
##              When distributions are not skewed, whatever one choses sd1 or sd2, results are identical            ##
######################################################################################################################

B=test(sd=1,m1=0,m2=0,skew=0,n=20) # symmetric distribution

cor(B[,1],B[,3]) # no correlation between sd1 and meandiff
cor(B[,2],B[,3]) # no correlation between sd2 and meandiff
cor(B[,1],B[,2]) # no correlation between sd1 and sd2

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

# However, when distributions are skewed, one observes strange things

##################################################################
##                   studying right-skewed case                 ##
##################################################################

return(cbind(sd1,sd2,meandiff,glasspos1,glasspos2,mean1,mean2))

# Note: when y1 and y2 are left-skewed, we find opposite results

##### 1) m1=m2

A=test(nSims=100000,sd=1,m1=0,m2=0,skew=6.32,n=20)

# Glass is m1-m2, divided either by sd1 or by sd2.

# While the sampling distribution of both m1 and m2 are right-skewed
# The sampling distribution of m1-m2 is symmetric
plot(density(A[,3]))
skewness(A[,3]) # approximately 0

# The sampling distribution of s1 and s2 are right-skewed, and identical
plot(density(A[,1]))
lines(density(A[,2]),col="red")

## So why aren't the sampling distribution of both measures of glass identical?

# Because there is a correlation between mean and sd.
cor(A[,1],A[,6])  # positive correlation between m and sd (always true when y1 and y2 are right-skewed, because the sampling distribution of the mean 
cor(A[,2],A[,7])  # will be right-skewed - and as in all configurations, the sampling distribution of sd is of course also right-skewed).

                  # ATTENTION: we are talking about the SAMPLE mean and the SAMPLE sd! the estimates!


cor(A[,1],A[,3]) # positive correlation between meandiff and s1 
cor(A[,2],A[,3]) # negative correlation between meandiff and s2
cor(A[,1],A[,2]) # no correlation between sd1 and sd2

plot(density(A[,3]))
abline(v=0)

# when m1-m2=0 (vertical line on the plot) <--> m1 = m2 <--> s1=s2.  Whatever one divides M.diff by sd1 or sd2, it results in the same Glass's measure.

# when m1 > m2 (right part under the curve): m1-m2 > 0 and s1 > s2
# when m1 < m2 (left part under the curve): m1-m2 < 0 and s1 < s2.

### Glass delta when m1-m2 is divided by s1 ##

# when m1-m2>  0 (right part under the curve): s1 > s2
# s1 is the largest sd estimate, and it can be very large (because the sampling dist of sd is right-skewed)
# a positive number divided by a very large s --> we stay close of 0

# when m1-m2<  0 (right part under the curve): s1 < s2
# s1 is the smallest sd estimate, but it cannot decrease as much as it can increase, because the sd cannot be < 0.
# For these reason, negative estimates of m1-m2 will be divided by smaller numbers than positive numbers, explaining why glass's d 
# divided by sd1 is left-skewed.

### Glass delta when m1-m2 is divided by s2 ##

# that's the opposite

# when m1-m2>  0 (right part under the curve): s1 > s2
# s2 is the smallest sd estimate

# when m1-m2<  0 (right part under the curve): s1 < s2
# s2 is the largest sd estimate




# M.diff is divided by the larger s, which can be very large as the sampling distribution of sd is right-skewed


# when m1 < m2 (left part under the curve): s1 < s2.

### Glass delta when m1-m2 is divided by s2 ##
mean(A[,4]) # below 0 (i.e. below the real delta)
mean(A[,5]) # above 0 (i.e. above the real delta)

cor(A[,3],A[,6]) # positive correlation between meandiff and m1 (the larger m1, the larger meandiff)
cor(A[,3],A[,7]) # negative correlation between meandiff and m2 (the lower m2, the larger meandiff) --> of course, as we compute m1-m2
cor(A[,6],A[,7]) # no correlation between m1 and m2

cor(A[,1],A[,2]) # no correlation between sd1 and sd2


table(A[,4]<A[,5]) # true in mosts cases (i.e. 86%)

table(A[,4]<0)
table(A[,5]<0)

median(A[,4]) # about 50% of glass's measures are > 0, and other half is below 0 (therefore, 0 is the median)
median(A[,5]) # but the sampling distribution of glass is either positive, or negative, median != mean

plot(density(A[,3])) # when mu1-mu2=0, the mean diff is symmetrical
skewness(A[,3])


# These results are not dependent on the 'rpearson' function (I've reproduced the test with rchisq and rsnorm, and I found similar results)
## Note: there is a correlation betwen SAMPLE mean and sample sd (but no correlation between POPULATION MEAN and sample sd).
# Important to realize it!

##### 2) m1!=m2

C=test(nSims=100000,sd=1,m1=1,m2=0,skew=6.32,n=20) # m1 > m2 --> les 2 glass sont right-skewed
#D=test(sd=1,m1=0,m2=1,skew=6.32,n=20) # m1 < m2 --> les 2 glass sont left-skewed

skewness(C[,4])
skewness(C[,5])

cor(C[,1],C[,3]) # 
cor(C[,2],C[,3]) # 
cor(C[,1],C[,2]) # 

mean(C[,4]) # 
mean(C[,5]) # 

cor(C[,3],C[,6]) # 
cor(C[,3],C[,7]) # 
cor(C[,6],C[,7]) # 

cor(C[,1],C[,6])  # 
cor(C[,2],C[,7])  # 

table(C[,4]<C[,5]) # 
table(C[,4]<0)
table(C[,5]<0)

median(C[,4]) # 
median(C[,5]) # 

plot(density(C[,3])) # when mu1-mu2=0, the mean diff is symmetrical
skewness(C[,3])


C=test(nSims=100000,sd=1,m1=5,m2=0,skew=6.32,n=20) # m1 > m2 --> les 2 glass sont right-skewed



##################################################################
##                  correlation between n and sd                ##
##################################################################

##### 1) does it depend on the mean?

test=function(sd,nSims=1000000,m1,m2,m3,m4,n,skew,kurt=95.75){
   
   sd1<-rep(0,nSims)
   sd2<-rep(0,nSims)
   sd3<-rep(0,nSims)
   sd4<-rep(0,nSims)
   
   mean1<-rep(0,nSims)
   mean2<-rep(0,nSims)
   mean3<-rep(0,nSims)
   mean4<-rep(0,nSims)
   
   for (i in 1:nSims){
      
      y1 <- rpearson(n,moments=c(m1,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
      y2 <- rpearson(n,moments=c(m2,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
      y3 <- rpearson(n,moments=c(m3,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
      y4 <- rpearson(n,moments=c(m4,sd^2,skewness=skew*(n-2)/sqrt(n*(n-1)),kurtosis=(kurt*(n-2)*(n-3)-6*(n-1))/(n^2-1)+3))  
      
      sd1[i] <- sd(y1)
      sd2[i] <- sd(y2)
      sd3[i] <- sd(y3)
      sd4[i] <- sd(y4)
      
      mean1[i] <- mean(y1)
      mean2[i] <- mean(y2)
      mean3[i] <- mean(y3)
      mean4[i] <- mean(y4)
      
   }
   
   return(cbind(mean1,mean2,mean3,mean4,sd1,sd2,sd3,sd4))
}

F=test(sd=1,m1=1,m2=2,m3=3,m4=4,skew=6.32,n=20)

cor(F[,1],F[,5])
cor(F[,2],F[,6])
cor(F[,3],F[,7])
cor(F[,4],F[,8]))

# No it does not. 