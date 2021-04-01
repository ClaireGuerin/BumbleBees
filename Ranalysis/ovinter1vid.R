#---- Import data -----

rm(list = ls())
#library(latex2exp)

setwd("H:/Academia/BumbleBees2016/Behav_Ovaries/Behav/Odyssey/allFiles/track")

compiled = read.table("compiledDataForR.csv", sep=",", na.strings="NA", dec=".", strip.white=TRUE, header=FALSE)

dataHeaders = read.table("compiledDataForR_headers.txt", header = FALSE, sep = '\t')
colnames(compiled) = t(dataHeaders)

compiled$id1 = as.factor(compiled$id1) # to transform back into original char form ASCII code, use function intToUtf8()
compiled$id2 = as.factor(compiled$id2)
compiled$colony = as.factor(compiled$colony)
compiled$chamber = as.factor(compiled$chamber) 

summary(compiled)

#---- Re-organize data ----

compiled$idcol1 = apply(compiled[,c('id1', 'colony')], 1, function(x) paste(x[1], x[2], collapse=':'))
compiled$idcol2 = apply(compiled[,c('id2', 'colony')], 1, function(x) paste(x[1], x[2], collapse=':'))
idcol = with(compiled, c(idcol1, idcol2))

allindiv = levels(factor(idcol))

meaninter = sapply(allindiv, function(currindiv) {
  idx1 = compiled$idcol1==currindiv
  idx2 = compiled$idcol2==currindiv
  idx = c(idx1, idx2)
  meaninter = mean(compiled$int.proba[idx], na.rm=T)
  
  return(meaninter)
})

ovscore = rep(NA,length(allindiv))
for(indiv in 1:length(allindiv)){
  idx1 = compiled$idcol1==allindiv[indiv]
  
  if(length(idx1) >= 1){
    ovscore[indiv] = compiled$ov1[idx1][1]
  }else{
    idx2 = compiled$idcol1==allindiv[indiv]
    ovscore[indiv] = compiled$ov1[idx2][1]
}}

