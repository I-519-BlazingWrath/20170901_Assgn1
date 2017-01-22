#setwd("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro")

plant<-read.csv("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\HF_plants.csv")
print(dim(plant))
#There are 96 rows. This means there are 96 "sites"" overall.
#However, sites were resampled yearly each year for 4 years. Too, there were four different treatments.
#Thus there are 96/4/4 = 6 replicate sites for each treatment for each year.
#There are 43 columns.
print(names(plant))
#Three columns correspond to year, plot, and treatment. The other 40 are spp.
#install.packages("dplyr")
require(dplyr)

C_2009_plant<-filter(plant, treatment==1, year==2009)
C_2009_plant<-C_2009_plant[4:43]
print(C_2009_plant)  
#class(plant.matrix.C)
#So, for each year & treatment combination, there are 6 sites.
print(S.obs(C_2009_plant))
#richness can vary a lot between plots: ranges 11 to 20 in these six plots.

#I will plot a histogram for the first site
plotme2<-as.numeric(C_2009_plant[1,])
print(plotme2)
hist(plotme2, xlab="Abundance", ylab="Number of Taxa", main="Histogram of Taxon Abundances", breaks=30)
#once again, many spp are rare, few are abundant

C09_Chao2<-c(1:6)
for (i in c(1:6)){
  temp<-(S.chao2(i, C_2009_plant))
  print(temp)
}
#estimated richness varies almost as much among these sites, here ranging from 12 to 21.
#estimated richness values are not much larger than observed richness,
#which indicated that the sites were sampled thoroughly.
C09_H<-c(1:6)
for (i in c(1:6)) {
  temp<-ShanH(C_2009_plant[i,])
  C09_H[i]<-temp
}
print(C09_H)
#Note that the plots with highest richness don't correspond very well to those with high Shannon's diversity

#Now I would like to draw a rank-abundance curve for the first site
C09_1_RAC_vector<-radfit(C_2009_plant[1,])
print(C09_1_RAC_vector)
plot(C09_1_RAC_vector, las = 1, cex.lab = 1.4, cex.axis = 1.25)
#Lowest AIC model is preemption, although for Mandelbrot DeltaAIC=1.568,
#which may not be enough difference to distinguish between the two models definitively.