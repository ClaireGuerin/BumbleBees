#---- Import data -----

rm(list = ls())

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

nIndiv = length(levels(ovData$ID))
tags = levels(ovData$ID)
indivCols = heat.colors(nIndiv)

plot(Interaction~ExpDay, data = ovData, col = indivCols[ID], pch = 20)

plot(NULL, xlim = c(min(ovData$ExpDay), max(ovData$ExpDay)), ylim = c(min(ovData$Interaction), max(ovData$Interaction)), xlab = "Day", ylab = "Interaction frequency")

for (i in 1:nIndiv){
  
  ind = tags[i]
  keepRows = which(ovData$ID == ind)
  tmpDat = ovData[keepRows,c(4,8)]
  points(Interaction~ExpDay, data = tmpDat, col = indivCols[i])
}

ovDvpt = ovData$OvScoreCeil
scores = sort(unique(ovDvpt))
nScores = length(scores)
ovCol = brewer.pal(n = nScores, name = "Dark2")

lastDay = which(ovData$ExpDay == 6)
firtData = ovData[-lastDay,]
dat = ovData
boxplot(Interaction~ID, data = dat, col = ovCol[ovDvpt+1])
legend("topright", legend = as.character(scores), fill=ovCol)

stabModMult = lm(Interaction~ID*ExpDay, data = ovData)
stabModAdd = lm(Interaction~ID+ExpDay, data = ovData)
stabModSin = lm(Interaction~ExpDay, data = ovData)
summary(stabModMult)
summary(stabModAdd)
summary(stabModSin)
plot(residuals(stabModSin))
qqnorm(residuals(stabModSin))

pcData = ovData
row.names(pcData) = pcData$ID

#res.pca = PCA()

crossTab = xtabs(Interaction ~ ID + OvScoreCeil, data = ovData)
res.chisq = chisq.test(crossTab)

# does interaction vary (mean) over time?
meanInteractions = tapply(X = ovData$Interaction, INDEX = as.factor(ovData$ExpDay), FUN = mean)
plot(sort(unique(ovData$ExpDay)), meanInteractions)

plot(Interaction~ExpDay, data = ovData)

day1 = ovData[which(ovData$ExpDay == 1),]
day1[which(day1$Interaction == 0),] = NA
day2 = ovData[which(ovData$ExpDay == 2),]
day2[which(day2$Interaction == 0),] = NA
day3 = ovData[which(ovData$ExpDay == 3),]
day3[which(day3$Interaction == 0),] = NA
day6 = ovData[which(ovData$ExpDay == 6),]
day6[which(day6$Interaction == 0),] = NA

par(mfrow=c(3,2))
plot(day2$Interaction~day1$Interaction, pch = 20, xlab = 'Day 1', ylab = 'Day 2', cex.lab = 1.5)
legend('topright', legend = as.character(round(cor(day2$Interaction,day1$Interaction, use = "pairwise.complete.obs", method = "spearman"),2)), bty = "n")
plot(day3$Interaction~day1$Interaction, pch = 20, xlab = 'Day 1', ylab = 'Day 3', cex.lab = 1.5)
legend('topright', legend = as.character(round(cor(day3$Interaction,day1$Interaction, use = "pairwise.complete.obs", method = "spearman"),2)), bty = "n")
plot(day6$Interaction~day1$Interaction, pch = 20, xlab = 'Day 1', ylab = 'Day 6', cex.lab = 1.5)
legend('topright', legend = as.character(round(cor(day6$Interaction,day1$Interaction, use = "pairwise.complete.obs", method = "spearman"),2)), bty = "n")
plot(day3$Interaction~day2$Interaction, pch = 20, xlab = 'Day 2', ylab = 'Day 3', cex.lab = 1.5)
legend('topright', legend = as.character(round(cor(day3$Interaction,day2$Interaction, use = "pairwise.complete.obs", method = "spearman"),2)), bty = "n")
plot(day6$Interaction~day2$Interaction, pch = 20, xlab = 'Day 2', ylab = 'Day 6', cex.lab = 1.5)
legend('topright', legend = as.character(round(cor(day6$Interaction,day2$Interaction, use = "pairwise.complete.obs", method = "spearman"),2)), bty = "n")
plot(day6$Interaction~day3$Interaction, pch = 20, xlab = 'Day 3', ylab = 'Day 6', cex.lab = 1.5)
legend('topright', legend = as.character(round(cor(day6$Interaction,day3$Interaction, use = "pairwise.complete.obs", method = "spearman"),2)), bty = "n")

indivMean = aggregate(Interaction~ID, data = ovData, FUN = mean)
indivMeanSort = indivMean[with(indivMean, order(Interaction)),]

indivMeanSort$ID
ovData$ID = factor(ovData$ID, levels = indivMeanSort$ID)
with(ovData, boxplot(Interaction~ID))
