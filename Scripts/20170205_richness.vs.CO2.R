plant<-read.csv("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\HF_plants.csv",header=TRUE,sep=",")
#print(dim(plant))
#print(names(plant))
S.obs <- function(x=""   ){
  rowSums(x>0) *1 
}
#install.packages("dplyr")
require(dplyr)
plants.sortedT<-arrange(plant, treatment)
plants.sortedY<-arrange(plants.sortedT,year)
plants.sortedP<-arrange(plants.sortedY,plot)
my2009plants.sortedY<-filter(plants.sortedY, year==2009, plot!=16 & plot!=12 & plot!=19)
#dim(my2009plants.sortedY)
sbs2009plants.sortedYt<-my2009plants.sortedY[4:43]
sbs2009plants.sortedYtt<-as.data.frame(sbs2009plants.sortedYt)

'''
my2008plants.sortedY<-filter(plants.sortedY, year==2008)
#dim(my2008plants.sortedY)
sbs2008plants.sortedYt<-my2008plants.sortedY[4:43]
sbs2008plants.sortedYtt<-as.data.frame(sbs2008plants.sortedYt)

my2007plants.sortedY<-filter(plants.sortedY, year==2007)
#dim(my2007plants.sortedY)
sbs2007plants.sortedYt<-my2007plants.sortedY[4:43]
sbs2007plants.sortedYtt<-as.data.frame(sbs2007plants.sortedYt)

my2006plants.sortedY<-filter(plants.sortedY, year==2006)
#dim(my2006plants.sortedY)
sbs2006plants.sortedYt<-my2006plants.sortedY[4:43]
sbs2006plants.sortedYtt<-as.data.frame(sbs2006plants.sortedYt)
'''
soilenv<-read.csv("C:\\Users\\rmoge\\GitHub\\QB2017_DivPro\\Data\\HF_soilresp.csv")
mysoilenv_holes<-na.omit(soilenv,year==2009,month==1)
my2009soilenv_1<-filter(soilenv,year==2009,month==1,plot!=12 & plot!=16 & plot!=19)
co2_1<-c(my2009soilenv_1[8])
co2_1<-
print(co2_1)
tester<-c(S.obs(sbs2009plants.sortedYtt))
sbs2009_1.sobs<-c(tester,tester)
toyr<-cbind(data.frame(sbs2009_1.sobs, co2_1))
lmfit<-lm(co2flux~sbs2009_1.sobs,data=toyr)
summary(lmfit)
plot(co2flux~sbs2009_1.sobs,data=toyr)
abline(lmfit)

temporary_soil<-my2009soilenv[-9]
sbs2009soilenv<-na.omit(as.matrix(temporary_soil[8:11]))
mycolnames<-colnames(sbs2009plants.sortedYtt, do.NULL = TRUE, prefix = "col")