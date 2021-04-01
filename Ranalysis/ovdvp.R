# 10/08/2017
# Claire Guerin
# Analyzing ovarian development data from bumble bee dissections (3 colonies A, B and C)

library(Rcmdr)
library(RColorBrewer) 
library(latex2exp)
library(grDevices)
#install.packages('xtable')
library(xtable)

A <- 
  read.table("H:/Academia/BumbleBees2016/Behav_Ovaries/Dissection/Data/AforR.txt",
             header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
B <- 
  read.table("H:/Academia/BumbleBees2016/Behav_Ovaries/Dissection/Data/BforR.txt",
             header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
C <- 
  read.table("H:/Academia/BumbleBees2016/Behav_Ovaries/Dissection/Data/CforR.txt",
             header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

summary(A)

brewCol = brewer.pal(n = 9, name = "YlOrRd")

#layout(mat = matrix(C(1,4,4,2,4,4,3,4,4), ncol=3, byrow = TRUE), widths = c(1,1,1), heights = c(1,1,1))
par(mfrow = c(3,1))
hist(A$ovscore, col=brewCol, main=TeX("Colony $A$"), xlab=TeX("Ovary Score"), cex.lab = 1.5, cex.main = 2)
hist(B$ovscore, col=brewCol, main=TeX("Colony $B$"), xlab=TeX("Ovary Score"), cex.lab = 1.5, cex.main = 2)
hist(C$ovscore, col=brewCol, main=TeX("Colony $C$"), xlab=TeX("Ovary Score"), cex.lab = 1.5, cex.main = 2)

ABC = list(A,B,C)

colony = do.call('c', mapply(function(x,y) rep(x, y), c('A', 'B', 'C'), sapply(ABC, nrow)))

ABC = do.call('rbind', ABC)

ABC$colony = as.factor(colony)

ABCQL = ABC[-which(ABC$ovscore == 4),]

ncolors = nlevels(ABC$colony)
colonycols <- brewer.pal(ncolors,'Dark2')
colonycols_alpha = adjustcolor(colonycols, alpha.f = 0.5)
colonyNum = apply(cbind(as.character(ABC$colony)), 1, utf8ToInt) - utf8ToInt('A') + 1

plot(ABC$ovscore~ABC$weight.mg, col = 'black', bg=colonycols_alpha[colonyNum], pch = 21, cex = 1.5, xlab = "Weight (mg)", ylab = "Ovarian development score")
legend('bottomright', legend = c('Colony A','Colony B','Colony C'), col = 'black', pch = 21, pt.bg = colonycols_alpha, bty = 'n', pt.cex = 1.5)

qlreg = lm(ovscore~weight.mg, data = ABCQL)
summary(qlreg)
abline(qlreg, lwd = 1.5)
text(250,4,TeX('$R^2 = 0.15$, $p-value = 0.0018$'))

Boxplot(weight.mg~colony, data=ABC, id.method="y")
Boxplot(weight.mg~colony, data=ABCQL, id.method="y")

AnovaModel.1 <- aov(weight.mg ~ colony, data=ABCQL)
summary(AnovaModel.1)
with(ABCQL, numSummary(weight.mg, groups=colony, statistics=c("mean", 
                                                              "sd")))
local({
  .Pairs <- glht(AnovaModel.1, linfct = mcp(colony = "Tukey"))
  print(summary(.Pairs)) # pairwise tests
  print(confint(.Pairs)) # confidence intervals
  print(cld(.Pairs)) # compact letter display
  old.oma <- par(oma=c(0,5,0,0))
  plot(confint(.Pairs))
  par(old.oma)
})

AnovaModel.2 <- aov(ovscore ~ colony, data=ABCQL)
summary(AnovaModel.2)
with(ABCQL, numSummary(ovscore, groups=colony, statistics=c("mean", "sd")))
local({
  .Pairs <- glht(AnovaModel.2, linfct = mcp(colony = "Tukey"))
  print(summary(.Pairs)) # pairwise tests
  print(confint(.Pairs)) # confidence intervals
  print(cld(.Pairs)) # compact letter display
  old.oma <- par(oma=c(0,5,0,0))
  plot(confint(.Pairs))
  par(old.oma)
})


contrasts(ABCQL$colony)=contr.poly(nlevels(ABCQL$colony)) 
model.ancova=aov(ovscore~weight.mg*colony, data=ABCQL)
Anova(model.ancova, type="III") 
summary.lm(model.ancova)
model.ancova.corrected = aov(ovscore~weight.mg, data=ABCQL)
Anova(model.ancova.corrected, type="III") 
summary.lm(model.ancova.corrected)


xtable(a, digits=0)