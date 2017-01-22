require(vegan)
require(dplyr)

#observed richness
S.obs <- function(x=""   ){
  rowSums(x>0) *1 
}


#Good's coverage
C<-function(x=""){
  1- (rowSums(x==1)/rowSums(x))
}
C(BCI)

xC<-function(x=""){
  (rowSums(x==1))
}
C(BCI)
xC(BCI)
yC<-function(x=""){
  (rowSums(x))
}


#Richness estimators
S.chao1<- function(x=""){
  S.obs(x)+(sum(x==1)^2 / (2*sum(x==2)))
}

S.chao2<-function(site="",SbyS=""){
  SbyS=as.data.frame(SbyS)
  x=SbyS[site, ]
  SbyS.pa<-(SbyS>0)*1
  Q1=sum(colSums(SbyS.pa)==1)
  Q2=sum(colSums(SbyS.pa)==2)
  S.chao2 = S.obs(x) + (Q1^2)/(2 * Q2)
  return(S.chao2)
}
#go back through the handout and make sure you understand this

S.ace<-function(x="",thresh=""){
  x<-x[x>0]#exclude zero abundance 
  S.abund<-length(which(x>thresh))
  S.rare <-length(which(x<=thresh))
  singlt <-length(which(x==1))
  N.rare <-sum(x[which(x<=thresh)])
  C.ace  <-1-(singlt/N.rare)
  i      <-c(1:thresh)
  count<-function(i,y){length(y[y==i])}
  a.1    <-sapply(i, count, x)
  f.1    <-(i*(i-1))*a.1
  G.ace  <-(S.rare/C.ace)*sum(f.1)/(N.rare*(N.rare-1))
  S.ace  <- S.abund + (S.rare/C.ace) + (singlt/C.ace)*max(G.ace,0)
  return(S.ace)
}

#Rank-abundance curve, or can use radfit from vegan
RAC <- function(x = ""){#begin RA curve function
  x = as.vector(x)#force a varb type
  x.ab = x[x > 0]#ignore spp for which no obsvns
  x.ab.ranked = x.ab[order(x.ab, decreasing = TRUE)]#rank ordered
  return(x.ab.ranked)#return the ranked vector so that it can be stored as an outside variable
}



#Simpson's Evenness
SimpE <- function(x = ""){
  S <- S.obs(x)#obsvd richness
  x = as.data.frame(x)
  D <- diversity(x, "inv")
  E <- (D)/S
  return(E)
}


#Evar
Evar <- function(x){#define
  x <- as.vector(x[x > 0])#take nonzero values of a vector
  1 - (2/pi)*atan(var(log(x)))#look at variance of the vector, then standarize it
}

#Shannon's H (diversity index)
ShanH <- function(x = ""){
  H = 0
  for (n_i in x){
    if(n_i > 0) {
      p = n_i / sum(x)
      H = H - p*log(p)
    }
  }
  return(H)
}


#Simpson's Dominance (use 1/D or 1-D for diversity)
SimpD <- function(x = ""){
  D = 0
  N = sum(x)
  for (n_i in x){
    D = D + (n_i^2)/(N^2)
  }
  return(D)
}

#D.inv <- 1/SimpD(site1)
#D.sub <- 1-SimpD(site1)

