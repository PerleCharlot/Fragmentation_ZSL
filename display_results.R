setwd("C:/ZSL/Madagascar")
library(data.table)
data <-fread('madagascar_loss_1992.csv')
#data <-fread('madagascar_loss_2015.csv')
data$percentAFS <- (data$loss_AFS*100)/data$area_BSS
data$percentSS <- (data$loss_SS*100)/data$area_BSS
data$percentAF <- (data$loss_AF*100)/data$area_BSS
data$percentFS <- (data$loss_FS*100)/data$area_BSS

a = hist(data$loss_AFS, labels=TRUE, nclass = 10,main='Habitat loss in Madagascar - AFS',
     xlab='km2', ylab = 'number of species')#km2
a$counts
b = hist(data$loss_SS, labels=TRUE, nclass = 10,main='Habitat loss in Madagascar - SS',
     xlab='km2', ylab = 'number of species')#km2
b$counts
c = hist(data$loss_AF, labels=TRUE, nclass = 10,main='Habitat loss in Madagascar - AF',
     xlab='km2', ylab = 'number of species')#km2
c$counts

a=hist(data$percentAFS, labels=TRUE, nclass=10, main='% habitat loss in Madagascar - AFS',
     xlab='%',ylab = 'number of species')#% inital suitable habitat
a$counts
b = hist(data$percentAF, labels=TRUE, nclass=10, main='% habitat loss in Madagascar - AF',
     xlab='%',ylab = 'number of species')#% inital suitable habitat
b$counts
c=hist(data$percentSS, labels=TRUE, nclass=2, main='% habitat loss in Madagascar - SS',
     xlab='%',ylab = 'number of species', xlim=c(0,60))#% inital suitable habitat
c$counts

### comparison 1992/2015
data1992 <-fread('madagascar_loss_1992.csv')#7626 in 1992
data2015 <-fread('madagascar_loss_2015.csv')
com <- data2015

#difference 1992-2015 in initial habitat area_BSS, km2 and %
data1992[,1] = order(data1992$idspecies)
data2015[,1] = order(data2015$idspecies)

com$diff_BSSkm <- data1992[,2] - data2015[,2] #km2
com$diff_BSS <- ((data1992[,2] - data2015[,2])/data1992[,2]) *100 #%
com$diff_AFS <- (data1992$area_AFS - data2015$area_AFS)/data1992$area_AFS *100#%
com$diff_AF <- (data1992$area_AF - data2015$area_AF)/data1992$area_AF *100
com$diff_SS <- (data1992$area_ASS - data2015$area_ASS)/data1992$area_ASS *100
com <- subset(com, select = c(idspecies,diff_BSSkm, diff_BSS, diff_AFS,diff_AF, diff_SS))

write.csv(com, file = 'Madagascar_92_15_comparison.csv')
