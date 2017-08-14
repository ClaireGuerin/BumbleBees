#---- Import dataset ----
rm(list = ls())

setwd("H:/Academia/BumbleBees2016/Behav_Ovaries/Behav/Odyssey/allFiles/")

ovData = read.table("ovExp.csv", sep=",", na.strings="NA", dec=".", strip.white=TRUE, header=FALSE)

colnames(ovData) = c("ID","Colony","Chamber","ExpDay","DayTime","Length","Width","Interaction","OvScore")

ovData$Chamber[which(ovData$Chamber == 1)] = "Nest"
ovData$Chamber[which(ovData$Chamber == 2)] = "Foraging"

ovData$ID = as.factor(ovData$ID)
ovData$Colony = as.factor(ovData$Colony)
ovData$Chamber = as.factor(ovData$Chamber)
ovData$OvScoreFloor = floor(ovData$OvScore)
ovData$OvScoreCeil = ceiling(ovData$OvScore)
ovData$Counts = as.integer(ovData$Interaction*1000000)

summary(ovData)

queenLess = ovData[ovData$Width<0.5,]

plot(Interaction~Length, data=queenLess, pch = 20, col = "gray")
plot(Interaction~Width, data=queenLess, pch = 20, col = "gray")
qls.m = lm(Interaction~Width*Length, data = queenLess)
summary(qls.m)

#---- Load libraries ----

# install.packages("pscl")
# hurdle and zero-inflated models
library(pscl) 
#install.packages("RColorBrewer") 
# nice colors for scientific plots
library(RColorBrewer) 
#install.packages("extrafont")
# CM font... and others
library(extrafont) 
loadfonts()
#font_install('fontcm') 
# violin-plots
library(vioplot) 
#install.packages("lmtest")
# Likelihood Ratio Test (LRT) when method anova(object, test = "lrt") not callable on object
library(lmtest) 
#install.packages("ResourceSelection")
# Hosmer-Lemeshow Goodness of Fit (GOF) Test
library(ResourceSelection)
#install.packages("pgirmess")
# Permutation test for lm, lme and glm (binomial and Poisson) objects
library(pgirmess)
library(FactoMineR)
library(scatterplot3d)

#---- Set up plotting parameters ----

#display.brewer.all() # to see available color palettes in brewer

old.par <- par(pch = 20)
on.exit(par(old.par))
pchCol = c(20,17,18)
colCol = c("deeppink4","turquoise4","firebrick3")
color_transparent30 <- adjustcolor(colCol, alpha.f = 0.3) 
color_transparent20 <- adjustcolor(colCol, alpha.f = 0.2)
brewCol = brewer.pal(n = 3, name = "Paired")
whichCol = brewCol

#---- Explore Data ----


plot(Interaction~OvScoreCeil, data = ovData)
outlier = ovData$ID[which(ovData$Interaction == max(ovData$Interaction[which(ovData$OvScoreCeil == 2)]))]
outDay = ovData$ExpDay[which(ovData$Interaction == max(ovData$Interaction[which(ovData$OvScoreCeil == 2)]))]
outChbr = ovData$Chamber[which(ovData$Interaction == max(ovData$Interaction[which(ovData$OvScoreCeil == 2)]))]
outDat = ovData[which(ovData$ID == outlier),]
except = which(ovData$ID == outlier & ovData$ExpDay == outDay & ovData$Chamber == outChbr)
ovData[except,]
ovDataExcept = ovData[-except,]

colPalette = brewer.pal(n = 8, name = "Dark2")

# Distribution Histogram
hist(ovData$Interaction, col = colPalette[2], breaks = 1000) 
# clear overrepresentation of zeros
NoZeroInteraction = ovData$Interaction[-which(ovData$Interaction == 0)]
hist(NoZeroInteraction, col = colPalette[3], breaks = 1000)

# Estimate frequency density
plot(density(ovData$Interaction), lwd = 2, col = colPalette[1], main = "Density estimate of interaction rates - Zeros included")
plot(density(NoZeroInteraction), lwd = 2, col = colPalette[1], main = "Density estimate of interaction rates - Zeros Excluded")
# Clear skewness of the data

# Empirical cumulative distribution function
plot(ecdf(ovData$Interaction), main = "Empirical cumulative distribution function - Zeros included", col = colPalette[6])
plot(ecdf(NoZeroInteraction), main = "Empirical cumulative distribution function - Zeros excluded", col = colPalette[6])

# Distributions that might fit:
# Chi-squared, exponential, F, gamma, geometric, log-normale, negative binomial, poisson, weibull

qqline(y, distribution = qchisq)
qqline(y, distribution = qchisq)
qqline(y, distribution = qchisq)
qqline(y, distribution = qchisq)

#---- Linear Models ----

# 1- Simple regression
ov.lm = lm(OvScore~Interaction, data=ovData)
summary(ov.lm)
# 2 - Log transform response variable ovary score

cst = 1e-100*range(ovData$Interaction)[2] # log(0)=-Inf can't be evaluated in lm
# solution: add a very small constant value
ovData$LogScore = log(ovData$OvScore+cst)
ov.lm.ln = lm(LogScore~Interaction, data=ovData)

# 3 - Without social interaction values of zero

ov.lm.nz = lm(LogScore~Interaction, data=ovDataNoZero) # even worse...
summary(ov.lm.nz)
# Extract parameters of the fit

coef  <- coefficients(fit)       # coefficients
resid <- residuals(fit)          # residuals
pred  <- predict(fit)            # fitted values
rsq   <- summary(fit)$r.squared  # R-sq for the fit
se    <- summary(fit)$sigma      # se of the fit

stat.coef  <- summary(fit)$coefficients
coef    <- stat.coef[,1]    # 1st column: coefficients (same as above)
se.coef <- stat.coef[,2]    # 2nd column: se for each coef
t.coef  <- stat.coef[,3]    # 3rd column: t-value for each coef
p.coef  <- stat.coef[,4]    # 4th column: p-value for each coefficient

#---- Visualize ----

plot(ovData$OvScore~ovData$Interaction, col=whichCol[ovData$Colony], cex=1.5, main = "Are cheaters asocial?", xlab = "Frequency of social interactions (per hour)", ylab = "Ovary development score")
abline(ov.lm, col='grey')

plot(ovData$LogScore~ovData$Interaction, col=whichCol[ovData$Colony], cex=1.5)
abline(ov.lm.ln, col='mediumseagreen')

# remove zero interactions
ovDataNoZero = ovData[-which(ovData$Interaction==0),]
plot(ovDataNoZero$OvScore~ovDataNoZero$Interaction, col=whichCol[ovData$Colony], cex=1.5)

#---- Qualitative response ----

# As categorical data

ovData$OvScoreCat = as.factor(ceiling(ovData$OvScore))
boxplot(Interaction~OvScoreCat, data=ovData)
levels(ovData$OvScoreCat)

# As categorical yes / no data

dvpOv = which(ovData$OvScore > 1)

ovData$OvCat = "Sterile"
ovData$OvCat[dvpOv] = "Fertile"
ovData$OvCat = as.factor(ovData$OvCat)

boxplot(ovData$Interaction~ovData$OvCat)

# As binary data

ovData$OvBin = 0
ovData$OvBin[dvpOv] = 1

summary(ovData)

with(ovData, plot(Interaction,OvBin))

#---- Poisson distribution ----

pois.log = glm(Counts ~ OvScore, data=ovData, family=poisson(link = "log"))
pois.ident = glm(Counts ~ OvScore, data=ovData, family=poisson(link = "identity"))
pois.ident.aug = glm(Counts + 1 ~ OvScore, data=ovData, family=poisson(link = "identity"))

summary(pois.log)
summary(pois.ident) # overdispersion in both: deviance much higher than df

plot(Counts ~ OvScore, data=ovData, col = whichCol[Colony], cex = 2)
points(ovData$OvScore, fitted(pois.log), col="red", pch = "*", cex = 2)
points(ovData$OvScore, fitted(pois.ident), col="green", pch = "*", cex = 2)

# Solution 1: Quasi-Poisson

quasipois.log = glm(Counts ~ OvScore, data=ovData, family=quasipoisson(link = "log"))
quasipois.ident = glm(Counts ~ OvScore, data=ovData, family=quasipoisson(link = "identity"))
quasipois.ident.aug = glm(Counts + 1 ~ OvScore, data=ovData, family=quasipoisson(link = "identity"))

#---- Hurdle ----

# Solution 2: Hurdle model
# Corrects for both overdispersion and overrepresentation of zero counts
# combines a count data model fcount(y; x, beta) (that is left-truncated at y = 1) and a zero hurdle model fzero(y; z, gamma) (right-censored at y = 1)

hurdle.pois = hurdle(Counts ~ OvScoreCeil, data = ovData, dist = "poisson", zero.dist = "binomial", link = "logit")
summary(hurdle.pois)

# Test significance of explanatory variable
# Agresti (1990) suggests that you should use the LRT instead of the Wald test for small sample sizes or if the parameters are large.

hurdletest(hurdle.pois) # can only be applied if dist and zer.dist are the same
linearHypothesis(hurdle.pois)

points(ovData$OvScore, predict(hurdle.pois), pch = "*", col = "blue", cex = 2)


AIC(hurdle.pois,pois.log, quasipois.log)
lrtest(hurdle.pois,pois.log, quasipois.log)

coef = coefficients(hurdle.pois)
beta = coef[2]
gamma = coef[4]

# Which distri: Poisson, Negative binomial or Geometric distribution?

#hurdle.nbin = hurdle(Counts ~ OvScore, data = ovData, dist = "negbin")
hurdle.nbin = hurdle(Counts ~ OvScoreCeil, data = ovData, dist = "negbin")
hurdle.geom = hurdle(Counts ~ OvScoreCeil, data = ovData, dist = "geometric")

AIC(hurdle.geom,hurdle.pois, hurdle.nbin)
lrtest(hurdle.geom,hurdle.pois, hurdle.nbin)

# nbin dispersion seems to explain best the data

qqline(residuals(hurdle.pois, type = "pearson"), distribution = qnbinom)

hoslem.test(ovData$Counts, fitted(hurdle.nbin))

#---- Compare scores as float, floor and ceiling ----
# 1 is float, 2 is floor and 3 is ceiling

hurdle.pois1 = hurdle(Counts ~ OvScore, data = ovData, dist = "poisson")
hurdle.nbin1 = hurdle(Counts ~ OvScore, data = ovData, dist = "negbin")
hurdle.geom1 = hurdle(Counts ~ OvScore, data = ovData, dist = "geometric")

hurdle.pois2 = hurdle(Counts ~ OvScoreFloor, data = ovData, dist = "poisson")
hurdle.nbin2 = hurdle(Counts ~ OvScoreFloor, data = ovData, dist = "negbin")
hurdle.geom2 = hurdle(Counts ~ OvScoreFloor, data = ovData, dist = "geometric")

hurdle.pois3 = hurdle(Counts ~ OvScoreCeil, data = ovData, dist = "poisson")
hurdle.nbin3 = hurdle(Counts ~ OvScoreCeil, data = ovData, dist = "negbin")
hurdle.geom3 = hurdle(Counts ~ OvScoreCeil, data = ovData, dist = "geometric")

aic.res = AIC(hurdle.pois1, hurdle.pois2, hurdle.pois3, hurdle.nbin1, hurdle.nbin2, hurdle.nbin3, hurdle.geom1, hurdle.geom2, hurdle.geom3)
ordered.aic = aic.res[order(aic.res$AIC),]

lrt.res = lrtest(hurdle.pois1, hurdle.pois2, hurdle.pois3, hurdle.nbin1, hurdle.nbin2, hurdle.nbin3, hurdle.geom1, hurdle.geom2, hurdle.geom3)
lrt.res[order(lrt.res[2]),]

wald.res = waldtest(hurdle.pois1, hurdle.pois2, hurdle.pois3, hurdle.nbin1, hurdle.nbin2, hurdle.nbin3, hurdle.geom1, hurdle.geom2, hurdle.geom3, test = "F") # ! Only compares nested models

#---- Goodness of fit ----
# Pearson Chi squared

ovDataExcept = ovData[-which(ovData$Counts>60000),]

dataSet = ovData
defMod = hurdle(Counts ~ OvScoreFloor, data = dataSet, dist = "negbin")
expVar = dataSet$OvScoreFloor


# Deviance goodness of fit test
pchisq(abs(defMod$loglik), df=defMod$df.residual, lower.tail=FALSE) # crap!!!

#---- Visual Diagnosis of the model ----

### Generate distribution and use violin plots 

imPath = "H:/Academia/MEME/S3 - HARVARD/Report/Figs/EffectOvSocInteraction.pdf"
pdf(imPath, family="CM Roman", width=14, height=14)

par(mfrow = c(1,1), mar = c(5.1, 5.1, 4.1, 2.1))

## Observed data
# probs = theta/(theta+mu)
scoresList = sort(unique(expVar))
theta = defMod$theta
meanObs = rep(NA,length(scoresList))
MUs = rep(NA,length(scoresList))
lenObs = rep(NA,length(scoresList))
obs = "black"
obs.trans = adjustcolor(obs, alpha.f = 0.6) 

## Simulated data
randomInt = c()
randomScores = c()
nGen = 10000
viol = "grey"
viol.trans = adjustcolor(viol, alpha.f = 0.6) 

jitter.a = 0.1

plot(NULL, xlim=c(-1,5), ylim=c(0,32500), ylab="Total interaction time (in number of frames)", xlab="Ovary development score", cex.lab = 2, cex.axis = 2)

for (i in 1:length(scoresList)){
  scoreI = which(expVar == i-1)
  meanObs[i] = mean(dataSet$Counts[scoreI]) 
  MUs[i] = defMod$fitted.values[scoreI][1]
  lenObs[i] = length(scoreI)
  rdNB = rnbinom(n = nGen, size = theta, mu = MUs[i])
  vioplot(rdNB, col=viol.trans, horizontal=FALSE, at=i-1, add=TRUE, lty = "blank", colMed = "tomato", drawRect = FALSE)
}
points(Counts~expVar, data = dataSet, pch = 20, col = obs.trans)
points(fitted(defMod)~expVar, pch = 18, col = "coral", cex = 2)
points(jitter(scoresList,amount = jitter.a), meanObs, pch = 18, col = "mediumseagreen", cex = 2)

legend("topright", legend = c("Individual Observations",expression(paste("Observed ", mu["i"])), expression(paste("Predicted ", mu["i"])),"Kernel density (NB)"), col = c(obs.trans,"mediumseagreen", "coral", NA), fill = c(NA, NA, NA, viol.trans), border = rep(NA, 4), pch = c(20,18,18,23), bty = "n", cex = 2, pt.cex = 2)

dev.off(dev.cur())

### Permutations
# Compare actual results to results randomly generated with same distribution but no effect of ovary development score 

# Step 1: Build an approximate permutation distribution
#PermTest(defMod, B = 1000) # Not working on hurdle class object

# Step 1: Create approximate permutation distribution
theo.permut = factorial(dim(ovData)[1])

nPermut = 19
modList <- vector("list", nPermut)

for (i in 1:nPermut){
  # Single permutation without replacement
  sPerm = sample(ovData$OvScoreCeil, length(ovData$OvScoreCeil), FALSE)
  modList[[i]] = hurdle(ovData$Counts~sPerm, dist = "negbin")
  plot()
}

aic.comp = AIC(modList[[1]],modList[[2]],modList[[3]],modList[[4]],modList[[5]],modList[[6]],modList[[7]],modList[[8]],modList[[9]],modList[[10]],modList[[11]],modList[[12]],modList[[13]],modList[[14]],modList[[15]],modList[[16]],modList[[17]],modList[[18]],modList[[19]],defMod)
aic.comp.ord = aic.comp[order(aic.comp$AIC),]

lrt.comp = lrtest(modList[[1]],modList[[2]],modList[[3]],modList[[4]],modList[[5]],modList[[6]],modList[[7]],modList[[8]],modList[[9]],modList[[10]],modList[[11]],modList[[12]],modList[[13]],modList[[14]],modList[[15]],modList[[16]],modList[[17]],modList[[18]],modList[[19]],defMod)
lrt.comp.ord = lrt.comp[order(lrt.comp[2]),]

par(mfrow = c(4,5))
j = 1
maxJ = 19

while (j <= 19){
  randomInt = c()
  randomScores = c()
  for (i in 1:5){
    intTemp = rnbinom(n = ns[i], size = theta, mu = mu[i])
    randomInt = append(randomInt,intTemp)
    randomScores = append(randomScores, rep(i-1, ns[i]))
  }
  plot(randomInt~randomScores)
  j = j + 1
}
plot(Counts~OvScore, data = ovData)

plot(hurdle.nbin$residuals~predict(hurdle.nbin))



dist <- replicate(2000, diff(by(ovData$Counts, sample(ovData$OvScoreCeil, length(ovData$OvScoreCeil), FALSE), mean)))
hist(dist, col = "black", breaks = 100)
abline(v = diff(by(ovData$Counts, ovData$OvScoreCeil, mean)), col = "blue", lwd = 2)



#---- Stable behaviour across days? ----

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

#---- Effect of each factor individually ----

# 1 - Effect of colony

boxplot(Interaction~Colony, data = ovData)
col.m = lm(Interaction~Colony, data = ovData)
colEffect = anova(col.m)
colEffect2 = aov(Interaction~Colony, data = ovData)

tukey = TukeyHSD(colEffect2) # colony '65' is significantly different from the others

# 2 - Effect of chamber

boxplot(Interaction~Chamber, data = ovData)
chb.m = lm(Interaction~Chamber, data = ovData)
chbEffect = anova(col.m)
chbEffect2 = aov(Interaction~Chamber, data = ovData)
summary(chbEffect2)

tukey = TukeyHSD(chbEffect2) # Interaction level significantly higher in FC

# 3 - Effect of day

plot(Interaction~ExpDay, data = ovData)
day.m = lm(Interaction~ExpDay, data = ovData)
summary(day.m) # no significant effect... in appearance. BUT non-independance of samples!!! --> should check effect of each factor, for each day

# 4 - Effect of Time of the day
plot(Interaction~DayTime, data = ovData)
tim.m = lm(Interaction~DayTime, data = ovData)
summary(tim.m) # there seems to be an effect, eventhough the time range is quite limited... Explains very little of the variation
abline(tim.m)

# 5 - Effect of bee body size (length and width)

len = aggregate(Length~ID, data = ovData, FUN = mean)
wid = aggregate(Width~ID, data = ovData, FUN = mean)
siz = merge(len,wid)

cor(siz$Length,siz$Width) # of course, highly correlated (0.8354607)
body.m = lm(Length~Width, data = siz)
summary(body.m) # p < 2e-16 ***; Rsq = 0.7
plot(Length~Width, data = siz, pch = 20, col = "grey")
abline(body.m, lwd = 2)

scatterplot3d(ovData$Width, ovData$Length, ovData$Interaction, pch = 20, color = "grey")

ovData$meanSize = rowMeans(subset(ovData, select = c(Width, Length)), na.rm = TRUE)
mlw.m = lm(Interaction~meanSize, data = ovData)
siz.m = lm(Interaction~Width*Length, data = ovData)
len.m = lm(Interaction~Length, data = ovData)
wid.m = lm(Interaction~Width, data = ovData)
all.m = lm(Interaction~Colony*ExpDay*DayTime*Length*Width*OvScoreCeil, data = ovData)
col.m = lm(Interaction~Colony, data = ovData)
ids.m = lm(Interaction~ID, data = ovData)
day.m = lm(Interaction~ExpDay, data = ovData) # but not independant
plot(Interaction~ExpDay, data = ovData, pch = 20, col = "gray")
abline(day.m)


# maybe a more relevant measure would be the approx. ground surface covered by the body? area for an ellipse is calculated by pi*a*b

ovData$Surf = pi*(ovData$Length/2)*(ovData$Width/2)
plot(Interaction~Surf, data = ovData, pch = 20, col = "grey")
srf.m = lm(Interaction~Surf, data = ovData) # no significant effect
plot(residuals(srf.m))
qqnorm(residuals(srf.m))

#---- Repeated measures Anova? ----
# PB: data not normally distributed...



