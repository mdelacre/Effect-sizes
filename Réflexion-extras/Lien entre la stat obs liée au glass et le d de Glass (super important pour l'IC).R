
Folder="D:/ES MEASURES/G1=0,G2=0"

j=75

filepath = paste0(Folder,"/",list.files(Folder)[j])
file=readRDS(filepath)

# Le titre des fichiers contenant les valeurs de tous les paramètres de la population,
# j'extrais ces différentes valeurs via str_extract_all
    param <- str_extract_all(list.files(Folder)[j], "[[:digit:]]+\\.*[[:digit:]]*")
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])
    N <- n1+n2
    nratio <- n1/n2

#distribution théorique attendue pour le d de glass, utilisant sd1 comme standardiseur

      mean1 <- file[,1]
      mean2 <- file[,2]
      s1 <- file[,3]
      s2 <- file[,4]

      t_obs <- ((mean1-mean2)/(s1*sqrt(1/n1+s2^2/(n2*s1^2))))

    glass1<-file[,11] #file[,11] = la colonne contenant les valeurs de d de glass utilisant sd1 comme standardiseur, pour chaque itération de la simulation
    glass1bis <- t_obs*sqrt(1/n1+s2^2/(n2*s1^2))

    round(glass1-glass1bis,6)

    mean1 <- file[,1]
    mean2 <- file[,2]
    s1 <- file[,3]
    s2 <- file[,4]

    t_obs <- ((mean1-mean2)/(s2*sqrt(1/n2+s1^2/(n1*s2^2))))
    glass2<-file[,12]
    glass2bis <- t_obs*sqrt(1/n2+s1^2/(n1*s2^2))

    round(glass2-glass2bis,6)

    mean1 <- file[,1]
    mean2 <- file[,2]
    s1 <- file[,3]
    s2 <- file[,4]

    t_obs <- (sqrt(n1*n2)*(mean1-mean2))/sqrt(n2*s1^2+n1*s2^2)
    cohenprime<-file[,17]
    cohenprimebis <- t_obs*sqrt((2*(n2*s1^2+n1*s2^2))/(n1*n2*(s1^2+s2^2)))
    round(cohenprime-cohenprimebis,6)


  plot(density(file[,17]))
  df=((n1-1)*(n2-1)*(s1^2+s2^2)^2)/((n2-1)*s1^4+(n1-1)*s2^4)
  ncp=((m1-m2)/sqrt((sd1^2+sd2^2)/2))*sqrt((n1*n2*(sd1^2+sd2^2))/(2*(n2*sd1^2+n1*sd2^2)))
  test <- sqrt((2*(n2*s1^2+n1*s2^2))/(n1*n2*(s1^2+s2^2)))*rt(1000000,df,ncp)
  lines(density(test),col="red")



#    lines(density()





    s_pooled <- sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2))
    t_obs2 <- (mean1-mean2)/sqrt(s_pooled^2*(1/n1+1/n2))
    cohen1 <- file[,9]
    cohen1bis <- t_obs2*sqrt((n1+n2)/(n1*n2))

    round(cohen1-cohen1bis,4)

