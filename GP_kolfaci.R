## Genomic Prediction Kolfaci Project
### Colombian corporation for agricultural research AGROSAVIA and EAFIT University
### Felipe López-Hernández et al. 2023

#### Uploading R libraries
library(BGLR)
library(ggplot2)
library(readr)
library(rcompanion)
library(readxl)
library(kableExtra)
library(GAPIT3)

#### Uploading the Allelic Variant File
# Using GAPIT for convert the Hapmap file format to numeric format

setwd("D:/AGROSAVIA_ 2020 _BUHO/PAPERS_2020/Korea Case/Paper GP/Genotype/Yield")

myG <- as.data.frame(read_delim("Yield_KKNimp.hmp.txt",
                                delim = "\t", escape_double = FALSE,
                                col_names = FALSE, trim_ws = TRUE))
myGAPIT <- GAPIT(G=myG, output.numerical=TRUE,
                 Geno.View.output=FALSE)

# Reading the numerical genotype file created by GAPIT
X <- read.delim("GAPIT.Genotype.Numerical.txt", row.names=1)

X[X==2] <- -1
X[X==1] <- 7
X[X==0] <- 1
X[X==7] <- 0

X[1:10,1:10] %>%
  kbl(caption = "Allelic Variant ") %>%
  kable_classic(full_width = F, html_font = "Cambria")
dim(X)

# Uploading the phenotypic data related to Yield
YLP1 <- read_delim("Index_YLP_Kolfaci.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
names <- as.data.frame(rownames(X))
colnames(names) <- c("Taxa")
YLP <- merge(YLP1, names, by="Taxa")
YLP <- YLP[,1:5]
YLP %>%
  kbl(caption = "Phenotypic Table Yield") %>%
  kable_classic(full_width = F, html_font = "Cambria")

# Defining the environment (e.g. 4 environments)
Environment <- 2
# 1, 2 or 3
y<-YLP[,Environment]

# Cross Validation
# Establishing seed
set.seed(0000)

# Cross-validation 
sets <- sample(c(1,2,3,4,5), size = 83, replace = TRUE)
for(i in 1:5){
  yNa=y
  whichNa=(sets==i)
  yNa[whichNa]=NA
  
  ETA=list(list(X=X,model="BayesC"))
  fmR<-BGLR(y=y,ETA=ETA,nIter=1000,burnIn=1000,thin=10)
  
  # Recording variances
  varU=scan("ETA_1_varU.dat")
  varE=scan("varE.dat")
  
  # Recording correlation in training dataset
  COR.trt<-cor(fmR$yHat,y, use="pairwise.complete.obs")
  if(i==1){CORst<-COR.trt}
  if(i!=1){CORst<-cbind(CORst,COR.trt)}
  
  # Recording correlation in testing dataset
  COR.tst<-cor(fmR$yHat[whichNa],y[whichNa], use="pairwise.complete.obs")
  if(i==1){CORs<-COR.tst}
  if(i!=1){CORs<-cbind(CORs,COR.tst)}
  
  # Recording heritabilities
  h2<-varU/(varU+varE)
  if(i==1){hs<-mean(h2,na.rm=T)}
  if(i!=1){hs<-cbind(hs,mean(h2,na.rm=T))}
  
  # Recording mean square error in training dataset
  MSE.trn<-mean(na.omit((fmR$yHat[-whichNa]-y[-whichNa])^2))
  if(i==1){MSE.trn_s<-mean(MSE.trn,na.rm=T)}
  if(i!=1){MSE.trn_s<-cbind(MSE.trn_s,mean(MSE.trn,na.rm=T))}
  
  # Recording mean square error in testing dataset
  MSE.tst<-mean(na.omit((fmR$yHat[whichNa]-y[whichNa])^2))
  if(i==1){MSE.tst_s<-mean(MSE.tst,na.rm=T)}
  if(i!=1){MSE.tst_s<-cbind(MSE.tst_s,mean(MSE.tst,na.rm=T))}
  
  # Genomic Estimated Breeding Values GEBV 
  GEBVs <- fmR$yHat-fmR$mu
  if(i==1){GEBVst_s<-GEBVs}
  if(i!=1){GEBVst_s<-cbind(GEBVst_s,GEBVs)}
  
  if(i==3){ ##Please select de best fold in the raw table from paper
    beta_hat<-fmR$ETA[[1]]$b
    # plot(beta_hat,col="royalblue1",cex=.7,pch=1)
    Betas <- as.data.frame(sort(beta_hat,decreasing=T))
    colnames(Betas) <- "beta"
    Betas <- cbind(rownames(Betas),Betas$beta)
    colnames(Betas) <- c("Molecuar Variant","beta")
    write.table(Betas,"Betas.csv", row.names = FALSE, sep = "\t")
  }
  
}

