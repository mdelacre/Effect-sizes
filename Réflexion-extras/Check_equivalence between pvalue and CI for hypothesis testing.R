library(stringr)

#############  Equivalence between hypothesis testing and proportion of true effect sizes out of the CI: checking #############

### Check in order to be sure that I computed confidence limits correcly (especially for transformed shieh's d)

# Check 1: 
# [P(Student's p-value > alpha) = percentage of Cohen's d confidence intervals covering 0] = TRUE?
# [P(Welch's p-value > alpha) = percentage of Shieh's d (and transformed Shieh's d) confidence intervals covering 0] = TRUE?

# Difference between proportion of pvalue > alpha and percentage of real effect size within confidence limits is computed for Cohen's d, Shieh's and transformed Shieh's d
# And stored in .txt files
# Expectations: difference should be null in any cases

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI")

Mainfolder="C:/Users/Marie/Documents/CI measures/check1/"
subfolder=list.files(Mainfolder)

for (k in seq_len(length(subfolder))){
  subfolder2=list.files(paste0(Mainfolder,subfolder[k],"/"))  
  Folder=paste0(Mainfolder,subfolder[k],"/",subfolder2)
  
  for (i in seq_len(length(Folder))){
    
    # set up empty container for all estimated parameters
    equiv <-matrix(0,length(list.files(Folder[i])),5+6)
    # 5 columns for n1,n2,m1,m2, sd1 and sd2
    # 2 columns for the proportion of Cohen's CI covering 0, and comparison between this frequency and the proportion of Student's p-value > alpha
    # 2 columns for the proportion of Shieh's CI covering 0, and comparison between this frequency and the proportion of Welch's p-value > alpha
    # 2 columns for the proportion of transformed Shieh's CI covering 0, and comparison between this frequency and the proportion of Welch's p-value > alpha
    
    colnames(equiv) <- c("conf.level","n1","n2","m1-m2","sd1/sd2",
                         "p_Cohen_NH0","p_Shieh_NH0","p_transfoshieh_NH0",
                         "CICohen_vs_Student","CIShieh_vs_Welch","CItransfoShieh_vs_Welch")
    
    for (j in seq_len(length(list.files(Folder[i])))){
      filename=list.files(Folder[i])[j]
      filepath = paste0(Folder[i],"/",filename)
      file=readRDS(filepath)   
      param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*") 
      conf.level <- as.numeric(param[[1]][1])
      n1 <- as.numeric(param[[1]][6])   
      n2 <- as.numeric(param[[1]][7])
      m1 <- as.numeric(param[[1]][8])
      m2 <- as.numeric(param[[1]][9])
      sd1 <- as.numeric(param[[1]][10])
      sd2 <- as.numeric(param[[1]][11])
      alternative=strsplit(filename,split=",")[[1]][2]
      
      Student_NRH0 <- sum(file[,1]>= (1-conf.level))/length(file[,1])
      Welch_NRH0 <- sum(file[,6]>= (1-conf.level))/length(file[,6])
      
      if (alternative=="two.sided"){ # If I compute a two sided confidence interval
        
        Cohen_NH0 <- sum(file[,4]<= 0 & file[,5] >= 0)/length(file[,4]) # Frequency of Cohen's d confidence intervals covering 0
        Shieh_NH0 <-  sum(file[,8]<= 0 & file[,9] >= 0)/length(file[,8]) # Frequency of Shieh's d confidence intervals covering 0
        transfoshieh_NH0 <- sum(file[,11]<= 0 & file[,12] >= 0)/length(file[,11]) # Frequency of transformed Shieh's d confidence intervals covering 0
        
      } else if (alternative=="greater"){
        
        Cohen_NH0 <- sum(file[,4]<= 0)/length(file[,4]) # Frequency of Cohen's d confidence intervals covering 0
        Shieh_NH0 <-  sum(file[,8]<= 0)/length(file[,8]) # Frequency of Shieh's d confidence intervals covering 0
        transfoshieh_NH0 <- sum(file[,11]<= 0)/length(file[,11]) # Frequency of transformed Shieh's d confidence intervals covering 0
        
        
      } else if (alternative=="less"){
        
        Cohen_NH0 <- sum(file[,5] >= 0)/length(file[,5]) # Frequency of Cohen's d confidence intervals covering 0
        Shieh_NH0 <-  sum(file[,9] >= 0)/length(file[,9]) # Frequency of Shieh's d confidence intervals covering 0
        transfoshieh_NH0 <- sum(file[,12] >= 0)/length(file[,12]) # Frequency of transformed Shieh's d confidence intervals covering 0
        
      }
      
      
      equiv[j,1] <- conf.level
      equiv[j,2] <- n1
      equiv[j,3] <- n2
      equiv[j,4] <- m1-m2
      equiv[j,5] <- sd1/sd2
      
      equiv[j,6] <- Cohen_NH0
      equiv[j,7] <- Shieh_NH0
      equiv[j,8] <- transfoshieh_NH0
      equiv[j,9] <- Cohen_NH0-Student_NRH0
      equiv[j,10] <-Shieh_NH0-Welch_NRH0  
      equiv[j,11] <-transfoshieh_NH0-Welch_NRH0
      
    }
    
    setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI")
    shapeparam<-str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
    sub<-paste0(subfolder[k],"; G1=",shapeparam[[1]][2],",G2=",shapeparam[[1]][4]) # Sign is not extracted but it's not a big deal
    write.table(equiv,paste0(sub,",equivalence.txt"),sep=";",dec=",")
    
  }
}

# Are there non null differences detected?

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI/check1/")
folders=list.files("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI/check1/")

discrep <- matrix(0,21,2)
colnames(discrep) <- c("min discrepancy","max discrepancy")
for (i in seq_len(length(folders))){
  file<-read.table(folders[i],sep=";",dec=",")  
  discrep[i,]<-c(min(file[,9:11]),max(file[,9:11]))
}

discrep
min(discrep)
max(discrep) # The biggest difference is 6e-05.


# Check 2: 
# If there is equivalence between hypothesis testing with pvalue and confidence limits, 
# the number of rows for which p-value > (or <) alpha should be exactly the same as the 
# number of rows for which p-value > (or <) alpha AND simultaneously
# real effect size is within (or out of) confidence limits of this effect size

setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI/check2")

Mainfolder="C:/Users/Marie/Documents/CI measures/"
subfolder=list.files(Mainfolder)

for (k in seq_len(length(subfolder))){
  subfolder2=list.files(paste0(Mainfolder,subfolder[k],"/"))  
  Folder=paste0(Mainfolder,subfolder[k],"/",subfolder2)
  
  for (i in seq_len(length(Folder))){
    
    # set up empty container for all estimated parameters
    equiv <-matrix(0,length(list.files(Folder[i])),5+6)
    # 5 columns for n1,n2,m1,m2, sd1 and sd2
    # 2 columns for cohen's d
    # 2 columns for shieh's d
    # 2 columns for transformed shieh's d
    
    colnames(equiv) <- c("conf.level","n1","n2","m1-m2","sd1/sd2",
                         "Cohen_TRUE","cohen_FALSE",
                         "shieh_TRUE","shieh_FALSE",
                         "trshieh_TRUE","trshieh_FALSE")
    
    for (j in seq_len(length(list.files(Folder[i])))){
      filename=list.files(Folder[i])[j]
      filepath = paste0(Folder[i],"/",filename)
      file=readRDS(filepath)   
      param <- str_extract_all(list.files(Folder[i])[j], "[[:digit:]]+\\.*[[:digit:]]*") 
      conf.level <- as.numeric(param[[1]][1])
      n1 <- as.numeric(param[[1]][6])   
      n2 <- as.numeric(param[[1]][7])
      m1 <- as.numeric(param[[1]][8])
      m2 <- as.numeric(param[[1]][9])
      sd1 <- as.numeric(param[[1]][10])
      sd2 <- as.numeric(param[[1]][11])
      alternative=strsplit(filename,split=",")[[1]][2]
      
      
      # If there is equivalence between hypothesis testing with pvalue and confidence limits, 
      # difference between 
      ## the number of rows for which p-value > (or <) alpha, and   
      #  the number of rows for which simultaneously p-value > (or <) alpha AND real effect size is within (or out of) confidence limits of this effect size
      
      Cohen_equiv <-  table(file[,1]>.05)-table(file[,1]>.05&(file[,2]>file[,4]&file[,2]<file[,5]))         #file[,1]= Welch's p-value
                                                                                                            #file[,2]= Cohen's d-value
                                                                                                            #file[,4&5]= lower and upper cohen's d limits
      
      Shieh_equiv <-   table(file[,6]>.05)-table(file[,6]>.05&(file[,7]>file[,8]&file[,7]<file[,9]))        #file[,6]= Welch's p-value   
                                                                                                            #file[,7]= shieh's d-value
                                                                                                            #file[,8&9]= lower and upper shieh's d limits
      
      trshieh_equiv <-  table(file[,6]>.05)-table(file[,6]>.05&(file[,10]>file[,11]&file[,10]<file[,12]))   #file[,6]= Welch's p-value
                                                                                                            #file[,10]= transformed shieh's d-value
                                                                                                            #file[,11&12]= lower and upper transformed shieh's d limits
      
      equiv[j,1] <- conf.level
      equiv[j,2] <- n1
      equiv[j,3] <- n2
      equiv[j,4] <- m1-m2
      equiv[j,5] <- sd1/sd2
      
      equiv[j,6:7] <- Cohen_equiv
      equiv[j,8:9] <- Shieh_equiv
      equiv[j,10:11] <- trshieh_equiv 

    }
    
    setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI/check2/")
    shapeparam<-str_extract_all(Folder[i], "[[:digit:]]+\\.*[[:digit:]]*")
    sub<-paste0(subfolder[k],"; G1=",shapeparam[[1]][2],",G2=",shapeparam[[1]][4]) # Sign is not extracted but it's not a big deal
    write.table(equiv,paste0(sub,",equivalence2.txt"),sep=";",dec=",")
    
  }
}


# Do we strictly have 0 everywhere?
setwd("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI/check2/")
folders=list.files("C:/Users/Marie/Documents/Github_projects/Effect-sizes/Scripts outputs/Confidence intervals/Equivalence pval CI/check2/")

diff <- data.frame(0)
for (i in seq_len(length(folders))){
  file<-read.table(folders[i],sep=";",dec=",")
  rep<-as.numeric(c(file[,6],file[,7],file[,8],file[,9],file[,10],file[,11]))
  j= length(names(table(rep)))
  diff[i,1]<-length(rep)
  diff[i,2:(j+1)]<-names(table(rep))
}


colnames(diff) <- c("nb",paste(expand.grid("rep",1:j)[,1],expand.grid("rep",1:j)[,2]))
diff # Only one value everywhere = strictly 0 --> Strict equivalence between both methods of hypothesis testing!
