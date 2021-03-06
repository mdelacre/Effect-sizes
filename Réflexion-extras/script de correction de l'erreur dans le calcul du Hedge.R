library(stringr)

Mainfolder="C:/Users/Marie/Documents/ES MEASURES/"
subfolder=list.files(Mainfolder)
Folder=paste0(Mainfolder,subfolder)

for (i in 4){

  for (j in seq_len(length(list.files(Folder[i])))){ 


    filepath = paste0(Folder[i],"/",list.files(Folder[i])[j])
    file=readRDS(filepath)

    param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*")
    if(as.numeric(param[[1]][2])==2.08){
      G1=-2.08
    } else {G1=as.numeric(param[[1]][2])}
    G2 <- as.numeric(param[[1]][4])
    n1 <- as.numeric(param[[1]][5])
    n2 <- as.numeric(param[[1]][6])
    N <- n1+n2
    m1 <- as.numeric(param[[1]][7])
    m2 <- as.numeric(param[[1]][8])
    sd1 <- as.numeric(param[[1]][9])
    sd2 <- as.numeric(param[[1]][10])

    # ligne contenant le g de Hedges:        
    file[,10]=file[,9]*(1-3/(4*(n1+n2)-9))  

    # Extraction of the ES matrix
    chem <- "C:/Users/Marie/Documents/ES MEASURES/"
    sschem<- paste0("G1=",G1,",G2=",G2) 
    setwd(dir=paste0(chem,sschem)) # destination file  
    
    fname=paste0("G1=",G1,", G2=",G2,", n=[",n1,",",n2,"], means=[",m1,",",m2,"], sds=[",sd1,",",sd2,"].rds")
    saveRDS(file, file = fname)
    
  }
  
}
