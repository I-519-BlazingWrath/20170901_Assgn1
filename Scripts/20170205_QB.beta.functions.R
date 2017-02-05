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

