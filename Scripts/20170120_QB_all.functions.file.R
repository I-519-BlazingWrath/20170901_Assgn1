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




#Beta diversity functions
# Beta diversity function that calcuates Whittaker's Beta diversity
beta.w1 = function(site.by.species = ""){
  SbyS.pa = decostand(site.by.species, method = "pa")
  S = ncol(SbyS.pa[,which(colSums(SbyS.pa) > 0)])
  a.bar = mean(specnumber(SbyS.pa))
  b.w = round(S/a.bar, 3)
  return(b.w)
}

# Beta diversity function that calcuates Whittaker's Beta diversity- the inpt is a sbys matrix with optional arguments 
beta.w2 = function(site.by.species = "", sitenum1 = "", sitenum2 = "", pairwise = FALSE){
  if (pairwise == TRUE){
    if (sitenum1 == "" | sitenum2 == "") {
      print("Error: please specify sites to compare") 
      return(NA)}
    site1 = site.by.species[sitenum1,]
    site2 = site.by.species[sitenum2,]
    site1 = subset(site1, select = site1 > 0)
    site2 = subset(site2, select = site2 > 0)
    gamma = union(colnames(site1), colnames(site2))
    s = length(gamma)
    a.bar = mean(c(specnumber(site1), specnumber(site2)))
    b.w = round(s/a.bar - 1, 3)
    return(b.w)
  }
  else{
    SbyS.pa = decostand(site.by.species, method = "pa")
    S = ncol(SbyS.pa[,which(colSums(SbyS.pa) > 0)])
    a.bar = mean(specnumber(SbyS.pa))
    b.w = round(S/a.bar, 3)
    SbyS.pa <- decostand(site.by.species, method = "pa")
    S = ncol(SbyS.pa[,which(colSums(SbyS.pa) > 0)])
    a.bar = mean(specnumber(SbyS.pa))
    b.w = round(S/a.bar, 3)
    return(c(b.w,a.bar,S))
  }
}
# The output is beta, alpha, gamma- in that order

# Add Species scores function...The add.spec.scores() function in the BiodiversityR package produces species coordinates in ordination space with reflect the strength and direction each species has on ordination at the different sites.
`add.spec.scores` <-function(ordi,comm,method="cor.scores",multi=1,Rscale=F,scaling="1") {
  ordiscores <- scores(ordi,display="sites")
  n <- ncol(comm)
  p <- ncol(ordiscores)
  specscores <- array(NA,dim=c(n,p))
  rownames(specscores) <- colnames(comm)
  colnames(specscores) <- colnames(ordiscores)
  if (method == "cor.scores") {
    for (i in 1:n) {
      for (j in 1:p) {specscores[i,j] <- cor(comm[,i],ordiscores[,j],method="pearson")}
    }
  }
  if (method == "wa.scores") {specscores <- wascores(ordiscores,comm)}
  if (method == "pcoa.scores") {
    rownames(ordiscores) <- rownames(comm)
    eigenv <- ordi$eig
    accounted <- sum(eigenv)
    tot <- 2*(accounted/ordi$GOF[2])-(accounted/ordi$GOF[1])
    eigen.var <- eigenv/(nrow(comm)-1)
    neg <- length(eigenv[eigenv<0])
    pos <- length(eigenv[eigenv>0])
    tot <- tot/(nrow(comm)-1)
    eigen.percen <- 100*eigen.var/tot
    eigen.cumpercen <- cumsum(eigen.percen)
    constant <- ((nrow(comm)-1)*tot)^0.25
    ordiscores <- ordiscores * (nrow(comm)-1)^-0.5 * tot^-0.5 * constant
    p1 <- min(p, pos)
    for (i in 1:n) {
      for (j in 1:p1) {
        specscores[i,j] <- cor(comm[,i],ordiscores[,j])*sd(comm[,i])/sd(ordiscores[,j])
        if(is.na(specscores[i,j])) {specscores[i,j]<-0}
      }
    }
    if (Rscale==T && scaling=="2") {
      percen <- eigen.var/tot
      percen <- percen^0.5
      ordiscores <- sweep(ordiscores,2,percen,"/")   
      specscores <- sweep(specscores,2,percen,"*")
    }
    if (Rscale==F) {
      specscores <- specscores / constant
      ordiscores <- ordi$points
    }        
    ordi$points <- ordiscores
    ordi$eig <- eigen.var
    ordi$eig.percen <- eigen.percen
    ordi$eig.cumpercen <- eigen.cumpercen
    ordi$eigen.total <- tot
    ordi$R.constant <- constant
    ordi$Rscale <- Rscale
    ordi$scaling <- scaling
  }
  specscores <- specscores * multi    
  ordi$cproj <- specscores
  return(ordi)
}



#######Species spatial abundance distribution calculator
ssad <- function(x){
  ad <- c(2, 3)
  ad <- OTUs[, otu]
  ad = as.vector(t(x = ad))
  ad = ad[ad > 0]
}