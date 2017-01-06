##############################################################
#######        Microarray Data Preprocessing             #####
########       Alison G Paquette, January 3rd 2017  ##########
##############################################################


#Part 0A: Load All packages & Functions used in this analysis####

library(oligo) #For Affymetrix Data
library(lumi) #For Illumina Data

#We use Biomart to convert probes to Genes
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

#Functions Written for this Analysis
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

RemoveDupProbes<-function(IDs, NormData){
  Report<-as.data.frame(matrix(NA,nrow = 4,ncol=1))
  rownames(Report)<-c("Initial Probes","Probes w/ no Genes","Duplicated Probes","Unique Probes")

  probeids=row.names(NormData)
  Report[1,1]<-length(probeids)

  ArrayProbes<-unique(IDs[,1])
  Report[2,1]<-(length(probeids)-length(ArrayProbes))

  #Isolate Duplicated Probes
  Duplicated<-subset(IDs,duplicated(IDs[,1])|duplicated(IDs[,1],fromLast=T)==T)

  DupTable<-as.data.frame(table(Duplicated[,1]))
  rownames(DupTable)<-as.character(DupTable$Var1)

  #Create list of Data that is Duplicated and subset this Data
  Duplications<-rownames(DupTable)
  Report[3,1]<-length(Duplications)

  #Isolate Unique Probes
  UniqueProbes<-outersect(ArrayProbes,Duplications)
  Unique<-subset(IDs,duplicated(IDs[,1])==F)
  rownames(Unique)<-as.character(Unique[,1])

  if (dim(Unique)[1]!=length(ArrayProbes)){
    stop("Number of Unique Probes Not Equal To Array Probes")
  }

  Unique<-(Unique[UniqueProbes,])
  rownames(Unique)<-Unique[,1]
  Report[4,1]<-length(rownames(Unique))

  print(Report)
  #Create Dataset with these unique probes
  Unique.Data<-merge(Unique,(NormData),by='row.names',all=F)
  Unique.Data<-Unique.Data[,-1]
  if (dim(Unique)[1]+(dim(DupTable)[1])!=length(ArrayProbes)){
    stop("Unique and Non Unique Probes do NOT equal Probe IDs")
  }
  (Unique.Data)
}

CondenseProbes<-function(ProbeData1){

  rownames(ProbeData1)<-ProbeData1[,1]
  Ensembl_IDs<-list(as.character(ProbeData1[,2]))
  TempData<-ProbeData1[,-c(1:3)]
  CollapsedData<-aggregate(TempData,by=Ensembl_IDs,FUN=mean,na.rm=TRUE, na.action=NULL)
  rownames(CollapsedData)<-CollapsedData[,1]
  CollapsedData<-CollapsedData[,-1]

}

#Part OB: Load covariate data and identify what samples to use####
#Here, I am loading all covariate data available on GEO from publicly available Datasets, which was manually curated

Covar<-read.csv("~/Documents/Project 1. Preterm Birth/PaperPreprocessingJan/Preprocessing/1262016UpdatedCovar.csv")

rownames(Covar)<-as.character(Covar$SampleIdentifier)

#Eliminating PE infants, Infants with known Chorioamniotis, and Infants 36-38 weeks
NoPE<-subset(Covar,!PETstatus=="Y"|is.na(PETstatus))
NoCA<-subset(NoPE,!CA.Diagnosis=="Yes"|is.na(CA.Diagnosis))
NoInbetween<-subset(NoCA,!is.na(birth.type))
Covar2$StudyIdentifier<-as.character(droplevels(Covar2$StudyIdentifier))
table(Covar2$StudyIdentifier)

#Generate table of term and preterm infants/Sutdy
Preterm<-subset(Covar2,birth.type=="PRETERM")
Term<-subset(Covar2,birth.type=="TERM")
rbind(table(Preterm$StudyIdentifier),table(Term$StudyIdentifier))


##Example 1: GSE73374: Processong w/ Affymetrix Hugene 2.0ST array (Gene Version) ####
**Example 1: : HuGene-2_0-st (Gene Version)**

#1. Subset Raw Data files BEFORE preprocessing
SampleID<-as.character(rownames(subset(Covar2,StudyIdentifier=="GSE73374")))

#Data munging: Getting .CEL messy filenames to The actual filenames
setwd("~/Documents/Project 1. Preterm Birth/GSE73374 Raw Data/GSE73374_RAW/")
filenames <- list.celfiles()
filenames.tmp<-as.data.frame(strsplit(filenames, split="_"))
filenames.tmp2<-cbind(t(filenames.tmp),filenames)
rownames(filenames.tmp2)<-filenames.tmp2[,1]
filenames.tmp2<-filenames.tmp2[SampleID,]
filenames<-filenames.tmp2[,3]

rawData <- read.celfiles(filenames)


#2. Preprocess/Prepare Datasets

NormData <- rma(rawData)#Note: The output of RMA using Oligo is already log2 normalized (unlike xps, previous packages, as described in https://www.bioconductor.org/packages/release/bioc/vignettes/oligo/inst/doc/oug.pdf)
NormData<-exprs(NormData)

par(mfrow=c(2,1))
boxplot(exprs(rawData),main=c(X,"Raw"))
boxplot(exprs(NormData),main=c(X,"Normalized"))


#3: Convert Probes to Genes
#Obtain Genes for each probe from BIOMART
IDs<-getBM(filters="affy_hugene_2_0_st_v1",attributes= c("affy_hugene_2_0_st_v1","ensembl_gene_id","hgnc_symbol"),values=row.names(NormData), mart=ensembl)

#Remove Duplicated Probes and print output of Probe processing
ProbeData1<-RemoveDupProbes(IDs,NormData)

#Condense Probes by mean
GSE73374<-CondenseProbes(ProbeData1)

#save(GSE73374,file="~/Documents/Project 1. Preterm Birth/PaperPreprocessingJan/Preprocessing/GSE733374_12282016.RData")

