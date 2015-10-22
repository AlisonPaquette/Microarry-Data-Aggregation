######################################################################
##############CODE TO PROCESS AND MERGE GEO DATASETS##################


##Merging affy, illumina excetera datasets that were non normalized
#Normalization done through lumi
#Merging through inSilicoDB

#for each group, step 1: load file, 2: concatenate redundant probes to genes 3: Insert Mod file into expressionset object 4: Save!


library(GEOquery)
library(lumi)
######################################################
#I have hand curated a meta covariate file (attached)
metacovar<-read.csv("IlluminaMetaCovariate2.csv")
PT<-subset(metacovar,PretermGroup=="PT")
rownames(PT)<-PT$SampleIdentifier
####################GSE25906##########################
SOI<-subset(PT,StudyIdentifier=="GSE25906B")
SAMPLEIDS<-rownames(SOI)
#This is the GEO non normalized file non normalized from the internet, or from Authors
fileName = "~/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/GSE25906_non-normalized_BatchB.txt"
mydata.lumi = lumiR(fileName, inputAnnotation = FALSE)
mydata.lumi<-mydata.lumi[ ,SAMPLEIDS] #subset by the samples you want
dim(exprs(mydata.lumi)) #Make sure this is the length of Sample IDS!
lumi.T<-lumiT(mydata.lumi,method='log2')
lumi.N<-lumiN(lumi.T,method="quantile")

#GSE25906: Step 2####
lumi.nuId <- addNuID2lumi(lumi.N, lib.mapping='lumiHumanIDMapping',verbose=T)
nuIDs<-rownames(exprs(lumi.nuId))
EntrezIDs<-as.data.frame(nuID2EntrezID(rownames(exprs(lumi.nuId)), lib.mapping='lumiHumanIDMapping'))
colnames(EntrezIDs)<-"EntrezID"
temp<-merge(EntrezIDs,(exprs(lumi.nuId)),by='row.names')
dim(temp)
##Remove Probes that Do not track to a gene
temp$EntrezID<-as.numeric(as.character(temp$EntrezID))
temp2<-subset(temp,!is.na(EntrezID))
temp2<-temp2[,-1]
dim(temp2)-dim(temp) #Lost -17946 Probes

#Merge probes by genes
x<-list(as.character(temp2$EntrezID))
temp3<-aggregate(temp2,by=x,FUN=mean)
rownames(temp3)<-temp3$EntrezID
temp3<-temp3[,-c(1,2)]
temp3<-as.matrix(temp3)

GEO<-getGEO('GSE25906')
gse <- GEO$GSE25906_series_matrix.txt.gz
gse <- gse[ ,SAMPLEIDS]
dim(exprs(gse)) ##Make sure this is equivalent to the number of sample ids you want!
head(exprs(gse))
exprs(gse)<-temp3##Insert your final file here!
head(exprs(gse))
#GSE25906: Step 4####
fDataOriginal<-fData(gse)
MyIlluminaData<-rownames(exprs(lumi.nuId))
fData2<-fDataOriginal[MyIlluminaData,]
fData(gse)<-fData2

GSE25906B<-gse
dim(exprs(GSE25906B))
dim(exprs(lumi.nuId))

save(GSE25906B,file="GSE25906BProcessed10212015.RData")


####################GSE54618##########################
#GSE54618: Step 1####
SOI<-subset(PT,StudyIdentifier=="GSE54618")
SAMPLEIDS<-rownames(SOI)
fileName = "~/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/GSE54618_nonnormalized.txt"
mydata.lumi = lumiR(fileName, inputAnnotation = FALSE)
mydata.lumi<-mydata.lumi[ ,SAMPLEIDS]
dim(exprs(mydata.lumi)) #Make sure this is the length of Sample IDS!
lumi.T<-lumiT(mydata.lumi,method='log2')
lumi.N<-lumiN(lumi.T,method="quantile")





#GSE54618: Step 2####
lumi.nuId <- addNuID2lumi(lumi.N, lib.mapping='lumiHumanIDMapping',verbose=T)
nuIDs<-rownames(exprs(lumi.nuId))
EntrezIDs<-as.data.frame(nuID2EntrezID(rownames(exprs(lumi.nuId)), lib.mapping='lumiHumanIDMapping'))
colnames(EntrezIDs)<-"EntrezID"
temp<-merge(EntrezIDs,(exprs(lumi.nuId)),by='row.names')
dim(temp)
##Remove Probes that Do not track to a gene
temp$EntrezID<-as.numeric(as.character(temp$EntrezID))
temp2<-subset(temp,!is.na(EntrezID))
temp2<-temp2[,-1]
dim(temp2)-dim(temp) #Lost -3271

#Merge probes by genes
x<-list(as.character(temp2$EntrezID))
temp3<-aggregate(temp2,by=x,FUN=mean)
rownames(temp3)<-temp3$EntrezID
temp3<-temp3[,-c(1,2)]
temp3<-as.matrix(temp3)
dim(temp3)-dim(temp2) #lost 13460

#GSE54618: Step 3####
GEO<-getGEO('GSE54618')
gse <- GEO$GSE54618_series_matrix.txt.gz
gse <- gse[ ,SAMPLEIDS]
dim(exprs(gse)) ##Make sure this is equivalent to the number of sample ids you want!
head(exprs(gse))
exprs(gse)<-temp3##Insert your final file here!
head(exprs(gse))
#GSE54618: Step 4####
fDataOriginal<-fData(gse)
MyIlluminaData<-rownames(exprs(lumi.nuId))
fData2<-fDataOriginal[MyIlluminaData,]
fData(gse)<-fData2

GSE54618<-gse
dim(exprs(GSE54618))
dim(exprs(lumi.nuId))

save(GSE54618,file="GSE54618Processed10212015.RData")

####################GSE44711##########################
#GSE44711: Step 1####
SOI<-subset(PT,StudyIdentifier=="GSE44711")
SAMPLEIDS<-rownames(SOI)
fileName = "~/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/GSE44711_non-normalized_data.txt"
mydata.lumi = lumiR(fileName, inputAnnotation = FALSE)
mydata.lumi<-mydata.lumi[ ,SAMPLEIDS]
dim(exprs(mydata.lumi)) #Make sure this is the length of Sample IDS!
lumi.T<-lumiT(mydata.lumi,method='log2')
lumi.N<-lumiN(lumi.T,method="quantile")



#GSE44711: Step 2####
lumi.nuId <- addNuID2lumi(lumi.N, lib.mapping='lumiHumanIDMapping',verbose=T)
nuIDs<-rownames(exprs(lumi.nuId))
EntrezIDs<-as.data.frame(nuID2EntrezID(rownames(exprs(lumi.nuId)), lib.mapping='lumiHumanIDMapping'))
colnames(EntrezIDs)<-"EntrezID"
temp<-merge(EntrezIDs,(exprs(lumi.nuId)),by='row.names')
dim(temp)
##Remove Probes that Do not track to a gene
temp$EntrezID<-as.numeric(as.character(temp$EntrezID))
temp2<-subset(temp,!is.na(EntrezID))
temp2<-temp2[,-1]
dim(temp2)-dim(temp) #Lost -3271

#Merge probes by genes
x<-list(as.character(temp2$EntrezID))
temp3<-aggregate(temp2,by=x,FUN=mean)
rownames(temp3)<-temp3$EntrezID
temp3<-temp3[,-c(1,2)]
temp3<-as.matrix(temp3)
dim(temp3)-dim(temp2) #lost 13460

#GSE44711: Step 3####
GEO<-getGEO('GSE44711')
gse <- GEO$GSE44711_series_matrix.txt.gz
gse <- gse[ ,SAMPLEIDS]
dim(exprs(gse)) ##Make sure this is equivalent to the number of sample ids you want!
head(exprs(gse))
exprs(gse)<-temp3##Insert your final file here!
head(exprs(gse))
#GSE44711: Step 4####
fDataOriginal<-fData(gse)
MyIlluminaData<-rownames(exprs(lumi.nuId))
fData2<-fDataOriginal[MyIlluminaData,]
fData(gse)<-fData2

GSE44711<-gse
dim(exprs(GSE44711))
dim(exprs(lumi.nuId))

save(GSE44711,file="GSE44711Processed10212015.RData")

###########MERGING DATA##################
library(inSilicoDb)
InSilicoLogin("YOUR LOGIN NAME","YOUR LOGIN PASSWORD")
library(inSilicoMerging)

##Load Files####


load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/GSE25906BProcessed10212015.RData")
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/GSE44711Processed10212015.RData")
load("/Users/alisonpaquette/Documents/Project 1. Preterm Birth/1021PTBStratifiedAnalysis/GSE54618Processed10212015.RData")
#
esets = list(GSE44711,GSE25906B,GSE54618)

eset_notnorm = merge(esets) #No batch correction
eset_COMBAT = merge(esets, method="COMBAT") #Emperical Bayes Batch Correction
eset_GENENORM= merge(esets, method="GENENORM") #Mean Centered Normalization


save(eset_notnorm,file="illuminamergedwnobatchcorrection.Rdata")
save(eset_COMBAT,file="illuminamergedwcombat.Rdata")
save(eset_GENENORM,file="illuminamergedwmeannorm.Rdata")


