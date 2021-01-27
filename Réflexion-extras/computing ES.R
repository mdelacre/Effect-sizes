##### D familiy

meandiff <- mu1-mu2

### Assuming normality AND homoscedasticity
# 1) Cohen's d
pooled_sd<-sqrt(((n1-1)*sd1^2+(n2-1)*sd2^2)/(n1+n2-2)) # pooled_sd
cohen_d <- meandiff/pooled_sd

### corrected version = hedge's g correction
hedges_cohen <- cohen_d*(1-3/(4*(n1+n2-2)-1))

### Assuming only normality
# 1) Keselman
average_sd <- sqrt((sd1^2+sd2^2)/2) # if first group = control group
average_d <- meandiff/average_sd # attention: c'est critiqué, car estimateur correspond à une population artificielle

# 2) Glass's d 
ctrl_sd <- sd1 # if first group = control group
glass_d <- meandiff/sd1
  
### corrected version = hedge's g correction
hedges_glass <- glass_d*(1-3/(4*(n1-1)-1))

# 3) Shieh's d
q1<-n1/n
q2<-n2/n
Shieh_stdz <- sqrt(sd1^2/q1+sd2^2/q2)

shieh_d <- meandiff/shieh_stdz



