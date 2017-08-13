#---- Import data -----

rm(list = ls())
library(latex2exp)

setwd("H:/Academia/BumbleBees2016/Behav_Ovaries/Behav")

compiled = read.table("compiledEllipseDataForR.csv", sep=",", na.strings="NA", dec=".", strip.white=TRUE, header=FALSE)

dataHeaders = read.table("compiledEllipseDataForR_headers.txt", header = FALSE, sep = '\t')
colnames(compiled) = t(dataHeaders)

compiled$id1 = as.factor(compiled$id1) # to transform back into original char form ASCII code, use function intToUtf8()
compiled$id2 = as.factor(compiled$id2)
compiled$colony = as.factor(compiled$colony)
compiled$chamber = as.factor(compiled$chamber) 

summary(compiled)

#---- Stable behavior across days? ----

compiled$idcol1 = apply(compiled[,c('id1', 'colony')], 1, function(x) paste(x[1], x[2], collapse=':'))
compiled$idcol2 = apply(compiled[,c('id2', 'colony')], 1, function(x) paste(x[1], x[2], collapse=':'))
idcol = with(compiled, c(idcol1, idcol2))
allindiv = levels(factor(idcol))

days = unique(compiled$day)

daySubsets = lapply(days, function(x) subset(compiled, compiled$day==x))

daysind = sapply(daySubsets, function(currsubset) {
  meaninter = sapply(allindiv, function(currindiv) {
    idx1 = currsubset$idcol1==currindiv
    idx2 = currsubset$idcol2==currindiv
    idx = c(idx1, idx2)
    meaninter = mean(currsubset$int.ellps[idx], na.rm=T)
    return(meaninter)
  })
  return(meaninter)
})


combs = combn(1:4, 2)

cortest = t(apply(combs, 2, function(comb) {
  cortest = cor.test(daysind[,comb[1]], daysind[,comb[2]], method='spearman', use='na.or.complete')
  pval = cortest$p.value
  corstat = cortest$estimate
  return(c(corstat,pval))
}))
0.001 = cbind(t(combs), cortest)
colnames(cortest) = c('day1', 'day2', 'cor', 'pval') 

colnames(daysind) = c('Day 1','Day 2','Day 3','Day 6')

panel.cor <- function(x, y, cex=2, pch=20)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  cortest = cor.test(x, y, method='spearman', use='na.or.complete')
  pval = round(cortest$p.value, 4)
  if(pval<0.001) pval = 'p-value < 0.001' else pval = paste0('p-value = ', pval)
  corstat = round(cortest$estimate, 2)
  txt = paste0('r = ', corstat, ', \n', pval)
  if(missing(cex)) cex <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex)
}

p = pairs(log(daysind), cex = 1.5, pch=20, row1attop = F, upper.panel = panel.cor, lower.panel = points)
mtext('Mean interaction (log)',1, line = 3.5)
mtext('Mean interaction (log)',2, line = 2.5)
