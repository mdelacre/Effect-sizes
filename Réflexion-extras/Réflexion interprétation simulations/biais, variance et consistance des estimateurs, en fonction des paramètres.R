### COHEN ###

N=seq(4,1000)
cohen=0

# Biais
#Rem.: on ne peut calculer le biais que si N >=4

Nsize=NULL
bias=NULL
for (i in seq_len(200)){
  N=i+3
  Nsize[i]=N
  bias[i] = ((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))-1
}

plot(Nsize,bias)

# Le terme par lequel est multiplié delta, pour calculer le biais, est TOUJOURS positif. C'est pourquoi on dit que le d de 
# Cohen surestime le biais.

# Il semblerait que le biais de Cohen ne dépende pas du nratio. Ca sera à checker via les simulations!

# Variance
#Rem.: on ne peut calculer le biais que si N >=5

cohen=0
N=5
n1=2
n2=N-n1

(((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*cohen^2))-cohen^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2

# Effet du nratio, pour valeurs de N et de Cohen constantes
var_cohen <- function(cohen,N){

  variance=NULL
  nratio=NULL
  for(i in seq_len(N-2)){
    n1=i+1
    n2=N-n1
    nratio[i] <- n1/n2
    variance[i] <- (((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*cohen^2))-cohen^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2  
  }
  
  plot(log(nratio),variance,main=paste("best n-ratio = ", nratio[variance==min(variance)]))  
  }

var_cohen(cohen=-.52,N=300)


# Effet du Cohen, pour un N et nratio constant
var_cohen2 <- function(nratio,N){
  n1=(N*nratio)/(1+nratio)
  n2=N-n1

  coh=seq(-1,1,.1)
  
  variance=NULL
  cohen=NULL
   for(i in seq_len(length(coh))){
    delta=coh[i]
    cohen[i] = delta
    variance[i] <- (((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*delta^2))-delta^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2  
  }
  
  plot(cohen,variance)  
}

var_cohen2(nratio=.1,N=100)


# Effet du N, pour un nratio et une valeur de d de Cohen constante
var_cohen3 <- function(nratio,cohen){
  
  Nsize=seq(4,300,2)
  
  variance=NULL 
  Ns=NULL
  for(i in seq_len(length(Nsize))){
    N=Nsize[i]
    Ns[i]=N
    n1=(N*nratio)/(1+nratio)
    n2=N-n1
    
    variance[i] <- (((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*cohen^2))-cohen^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2  
  }
  
  plot(Ns,variance)  
}

var_cohen3(nratio=13,cohen=-4)
#nratio[variance==min(variance)]
plot(N,bias,pch=19,cex=.3)
abline(h=0)

### HEDGES ###
coeff=NULL
nsize=NULL
for (i in seq_len(1000)){
N=i+3
nsize[i] <- N
coeff[i] <- (gamma((N-2)/2)/(sqrt((N-2)/2)*gamma((N-3)/2)))^2
}

min(coeff)
max(coeff)
nsize
plot(nsize,coeff)

### GLASS ###
# rem.: on ne peut calculer le biais que si n_c >=3
n=30 #n = n_c

mult <- NULL
Nsize <- NULL
for (i in seq_len(100)){
  n=i+2
  Nsize[i]=n
  mult[i]=((sqrt((n-1)/2)*gamma((n-2)/2))/gamma((n-1)/2))-1  
}

plot(mult,Nsize)

# Rem.: quand je compare Glass et Cohen, on voit clairement que le biais ne dépend absolument pas de nexp et N

# Variance

#1 ) Impact de delta Glass

sd_e=7
sd_c=1
glass=NULL
ne=20
nc=40
var <- NULL
term1 <- NULL
term2 <- NULL

tailleeff<-seq(1,10,.1)

for (i in seq_len(length(tailleeff))){
  glass_delta1=tailleeff[i]
  glass[i]=glass_delta1
  var[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)-glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
  term1[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)
  term2[i] <- glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
}

plot(glass,var) # Plus Nc grandit, moindre est la variance (logique)
abline(h=0)

# Plus la taille d'effet est grande, plus la variance est grande
# On peut démontrer que (nc-1)/(nc-3) est tj plus grand que ((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
# Voici une démo graphique (non mathématique)

A <- NULL
B <- NULL
Nc <- NULL
for (i in 5:200){
nc = i
Nc[i] <- nc
A[i] <-  (nc-1)/(nc-3)
B[i] <-  ((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
}

plot(Nc,A)
lines(Nc,B,col="red")

# Impact de l'augmentation de N, en gardant soit Nc soit Ne constant

test = function(sd_e,sd_c){

      # Plot 1: impact du nc grandissant, pour un ne constant  
      glass_delta1=1
      ne=100
      var <- NULL
      term1 <- NULL
      term2 <- NULL
    
      Ntot=NULL
      for (i in 1:200){
        nc=i+5
        Ntot[i]=nc+ne
        var[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)-glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
        term1[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)
        term2[i] <- glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
      }
    
      par(mfrow=c(1,2))
      limsup<-max(term1,term2,var)
      plot(Ntot,var,ylim=c(0,limsup),main="Nc variable") # Plus Nc grandit, moindre est la variance (logique)
      abline(h=0)
      # Mais elle ne tendra jamais vers 0, car le 1er terme est tj plus grand que le 2ème
      lines(Ntot,term1,col="red")
      lines(Ntot,term2,col="blue")
      #term1>term2 --> always TRUE
      abline(v=200)
      
      
      # Plot 2: impact du ne grandissant, pour un nc constant  
      glass_delta1=1
      nc=100
      var <- NULL
      term1 <- NULL
      term2 <- NULL
    
      Ntot=NULL
      for (i in 1:200){
        ne=i+5
        Ntot[i]=nc+ne
        var[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)-glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
        term1[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)
        term2[i] <- glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
      }
      
      limsup<-max(term1,term2,var)
      plot(Ntot,var,ylim=c(0,limsup),main="Ne variable") # Plus Nc grandit, moindre est la variance (logique)
      abline(h=0)
      # Mais elle ne tendra jamais vers 0, car le 1er terme est tj plus grand que le 2ème
      lines(Ntot,term1,col="red")
      lines(Ntot,term2,col="blue")
      #term1>term2 --> always TRUE
      abline(v=200)
}

# Pour un N constant, impact du n-ratio? 

test2 = function(sd_e,sd_c){
  
  glass_delta1=1
  N=100
  var <- NULL
  term1 <- NULL
  term2 <- NULL
  
  Nratio=NULL
  for (i in 1:89){
    nc=i+5
    ne=N-nc
    Nratio[i]=nc/ne
    var[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)-glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
    term1[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)
    term2[i] <- glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
  }
  
  limsup<-max(term1,term2,var)
  plot(log(Nratio),var,ylim=c(0,limsup),main="SDR positif quand sdc > sde") # Plus Nc grandit, moindre est la variance (logique)
  abline(h=0)
  # Mais elle ne tendra jamais vers 0, car le 1er terme est tj plus grand que le 2ème
  lines(log(Nratio),term1,col="red")
  lines(log(Nratio),term2,col="blue")
  abline(h=max(var),lty=2)
  #term1>term2 --> always TRUE
  
}

# garder soit Ne soit Nc constant, et faire varier l'autre...
test(sd_e=5,sd_c=5) # quand SDR=1 # Là, il faut mieux ajouter des gens dans le groupe Nc pour faire diminuer la variance
test(sd_e=1,sd_c=10) # quand sdC >> sdE, augmenter ne ne sert pratiquement à rien
test(sd_e=10,sd_c=1) # quand sdC << sdE mieux vaut augmenter Ne que Nc

Nratio[var==min(var)]
# garder N constant, et faire varier le ratio entre Nc et Ne...
test2(sd_e=5,sd_c=5) # quand SDR=1 # Là, il faut mieux ajouter des gensn dans le groupe Nc pour faire diminuer la variance
test2(sd_e=1,sd_c=10) # quand sdC >> sdE, augmenter ne ne sert pratiquement à rien
test2(sd_e=10,sd_c=1) # quand sdC << sdE mieux vaut augmenter Ne que Nc


# CCL
# Plus N est grand mieux c'est, mais pour un N constant, mieux vaut avoir plus de gens dans le groupe control que dans le groupe experimental (tant pour le biais que pour la variance)
Nratio[var==min(var)]


# QUand var1!=var2 et n1=n2

glass_delta1=1
sd_e=1
ne=20
nc=20
var <- NULL
term1 <- NULL
term2 <- NULL

sd_ratio=c(.1,.25,.5,2,4,10)
SDR=NULL
for (i in seq_len(length(sd_ratio))){
  sd_c=sd_e*sd_ratio[i]
  SDR[i]=sd_e/sd_c
  var[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)-glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
  term1[i] <- (nc-1)/(nc-3)*(1/nc+sd_e^2/(ne*sd_c^2)+glass_delta1^2)
  term2[i] <- glass_delta1^2*((sqrt((nc-1)/2)*gamma((nc-2)/2))/gamma((nc-1)/2))^2
}

maxi=max(var,term1,term2)
plot(log(SDR),var,ylim=c(0,maxi))
lines(log(SDR),term1,col="red")
lines(log(SDR),term2,col="blue")
abline(v=0)
### Shieh ###


# Est-ce que le biais dépend du nratio? pour un N constant

s1=2
s2=2
shieh = 2
N <- 120

DF <- NULL
BIAS <- NULL
NRATIO <- NULL
for (i in 1:100){
  n1 <- i+10
  n2 <- N-n1
  NRATIO[i] <- n1/n2
  df <- (s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
  DF[i] <- df
  BIAS[i] <- shieh*(sqrt(df/2)*gamma((df-1)/2))/(gamma(df/2))-1
  
}
  
plot(log(NRATIO),BIAS) # tant qu'il y a homoscédasticité, le biais est minimum quand le n ratio = 1.
NRATIO[BIAS==min(BIAS)]

# Qu'est-ce qui compte le plus? le N ou le nratio dans le calcul du biais?

s1=2
s2=2
shieh = 2
n1 <- 20

DF <- NULL
BIAS <- NULL
N <- NULL   
for (i in 1:100){
  n2 <- n1+i-1
  N[i] <- n1+n2
  df <- (s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
  DF[i] <- df
  BIAS[i] <- shieh*(sqrt(df/2)*gamma((df-1)/2))/(gamma(df/2))-1
  
}

plot(N,BIAS) # tant qu'il y a homoscédasticité, très vite, le nratio différent de 1 va plus jouer que l'augmentation du n
             # ce qui fait que si on rajoute tous les sujets dans un seul groupe, ça va augmenter le biais

N[BIAS==min(BIAS)]
# Qu'est-ce qui compte le plus? le N ou le nratio dans le calcul de la variance?

s1=2
s2=2
shieh = 2
n1 <- 20

DF <- NULL
VAR <- NULL
NTOT <- NULL   
for (i in 1:100){
  n2 <- n1+i-1
  N <- n1+n2
  NTOT[i] <- N
  df <- (s1^2/n1+s2^2/n2)^2/((s1^2/n1)^2/(n1-1)+(s2^2/n2)^2/(n2-1))
  DF[i] <- df
  VAR[i] <- df/((df-2)*N)*(1+N*shieh^2)-shieh^2*((sqrt(df/2)*gamma((df-1)/2))/gamma(df/2))^2
}

plot(NTOT,VAR) # tant qu'il y a homoscédasticité, très vite, le nratio différent de 1 va plus jouer que l'augmentation du n
# ce qui fait que si on rajoute tous les sujets dans un seul groupe, ça va augmenter le biais


























### COHEN D' ###

# quand N augmente, le biais diminue
Nsize=NULL
bias=NULL
DF=NULL

sd1=2
sd2=10

for (i in seq_len(200)){
   n1=15*i
  n2=i
  N = n1+n2
  Nsize[i]=N
  df = ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4) 
  DF[i] = df
  bias[i] = sqrt(df/2)*gamma((df-1)/2)/gamma(df/2)-1
}

plot(Nsize,bias)
plot(DF,bias)

plot(Nsize,bias)

# Le terme par lequel est multiplié delta, pour calculer le biais, est TOUJOURS positif. C'est pourquoi on dit que le d de 
# Cohen surestime le biais.
Nratio=NULL
bias=NULL
DF=NULL

sd1=2
sd2=10


for (i in seq_len(200)){
  N = 206
  n2=i+1
  n1=N-(i+1)
  Nratio[i]=n1/n2
  df = ((n1-1)*(n2-1)*(sd1^2+sd2^2)^2)/((n2-1)*sd1^4+(n1-1)*sd2^4) 
  DF[i] = df
  bias[i] = sqrt(df/2)*gamma((df-1)/2)/gamma(df/2)-1
}

plot(Nratio,bias)
plot(log(Nratio),bias)

Nratio[bias==min(bias)]


plot(DF,bias)

plot(Nsize,bias)



# quand Nratio augmente, pour un n constant?




# Il semblerait que le biais de Cohen ne dépende pas du nratio. Ca sera à checker via les simulations!

# Variance
#Rem.: on ne peut calculer le biais que si N >=5

cohen=0
N=5
n1=2
n2=N-n1

(((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*cohen^2))-cohen^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2

# Effet du nratio, pour valeurs de N et de Cohen constantes
var_cohen <- function(cohen,N){
  
  variance=NULL
  nratio=NULL
  for(i in seq_len(N-2)){
    n1=i+1
    n2=N-n1
    nratio[i] <- n1/n2
    variance[i] <- (((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*cohen^2))-cohen^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2  
  }
  
  plot(log(nratio),variance,main=paste("best n-ratio = ", nratio[variance==min(variance)]))  
}

var_cohen(cohen=-.52,N=300)


# Effet du Cohen, pour un N et nratio constant
var_cohen2 <- function(nratio,N){
  n1=(N*nratio)/(1+nratio)
  n2=N-n1
  
  coh=seq(-1,1,.1)
  
  variance=NULL
  cohen=NULL
  for(i in seq_len(length(coh))){
    delta=coh[i]
    cohen[i] = delta
    variance[i] <- (((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*delta^2))-delta^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2  
  }
  
  plot(cohen,variance)  
}

var_cohen2(nratio=.1,N=100)


# Effet du N, pour un nratio et une valeur de d de Cohen constante
var_cohen3 <- function(nratio,cohen){
  
  Nsize=seq(4,300,2)
  
  variance=NULL 
  Ns=NULL
  for(i in seq_len(length(Nsize))){
    N=Nsize[i]
    Ns[i]=N
    n1=(N*nratio)/(1+nratio)
    n2=N-n1
    
    variance[i] <- (((N-2)/((N-4)*n1*n2/N))*(1+n1*n2/N*cohen^2))-cohen^2*((sqrt((N-2)/2)*gamma((N-3)/2))/gamma((N-2)/2))^2  
  }
  
  plot(Ns,variance)  
}

var_cohen3(nratio=13,cohen=-4)
#nratio[variance==min(variance)]
plot(N,bias,pch=19,cex=.3)
abline(h=0)
