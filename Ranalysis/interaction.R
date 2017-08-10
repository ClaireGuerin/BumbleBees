#---- Import dataset ----
rm(list = ls())

setwd("H:/Academia/BumbleBees2016/Behav_Ovaries/Behav/Odyssey/allFiles")

compiled = read.table("compiledDataForR.csv", sep=",", na.strings="NA", dec=".", strip.white=TRUE, header=FALSE)

dataHeaders = read.table("compiledDataForR_headers.txt", header = FALSE, sep = '\t')
colnames(compiled) = t(dataHeaders)

compiled$id1 = as.factor(compiled$id1) # to transform back into original char form ASCII code, use function intToUtf8()
compiled$id2 = as.factor(compiled$id2)
compiled$colony = as.factor(compiled$colony)
compiled$chamber = as.factor(compiled$chamber) 

summary(compiled)

factorLevels = levels(compiled$colony)
nLevels = length(factorLevels)

for (factor in 1:nLevels){
  subsetName = paste(factorName,intToUtf8(factorLevels[factor]), sep = "")
  dataSubset = compiled[compiled[,factorInd] == factorLevels[factor],]
  assign(subsetName, dataSubset)
}

#---- work in progress ----
ids = unique(c(as.character(colonyA$id1),as.character(colonyA$id2)))

for(id in 1:length(ids)){
  findBin1 = which(colonyA$id1 == ids[id])
  findBin2 = which(colonyA$id2 == ids[id])
  
  if(length(findBin2) > 0){
    Bcharacteristics = colonyA[findBin2[1],2:4]
  }else{}
  
  meanInterE = mean(colonyA$int.ellps[c(findBin1,findBin2)], na.rm = TRUE)
  meanInterP = mean(colonyA$int.proba[c(findBin1,findBin2)], na.rm = TRUE)
}

##---- Description data ----

plot(colonyA$ov1,colonyA$int.ellps)
which(colonyA$int.ellps>0)
plot(jitter(colonyA$ov1,0.2),colonyA$int.proba)
which(colonyA$int.proba>0)
