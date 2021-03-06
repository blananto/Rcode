library(ncdf4) # NetCDF
library(fBasics) #timPalette
library(ismev) #gp.fit
library(splancs) #inpip
library(fields) #Tps
library(KernSmooth) #bkde2D
#library(MASS) #kde2d
library(plotrix) #draw.circle
library(mev) #egp2.fit
require(akima) #interp: bilinear interpolation
library(zoo) #rollapply
library(RANN) #nn2
library(scoringRules) #crps
library(evd) #fpot
require(mvtnorm)


compare.loglik.gamma.nei<-function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05"){
  descr<-rbind(
    c("sing","sing"),
    c("cel","sing"),
    #c("pers","sing"),
    c("pnei","sing"),
    c("snei","sing"),
    c("rsing","rsing"),
    c("cel","rsing"),
    #c("pers","rsing"),
    c("pnei","rsing"),
    c("snei","rsing"),
    #c("pers","rlsing")
    c("sing","rsing"))
  
  ndesc<-nrow(descr)
  
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE)
  
  strs<-c("05","1","2","5")
  
  loglik.mat<-NULL
  coln<-NULL
  for (str in strs){
    descriptors<-descr
    descriptors[,1]<-paste0(descriptors[,1],str)
    descriptors[,2]<-paste0(descriptors[,2],str)
    descriptors[substr(descriptors,1,3)=="cel"]<-"cel"
    print(descriptors)
    
    for (i in 1:nrow(descriptors)){
      print(descriptors[i,])
      load(file=paste0(get.dirstr(k),"fit.loglik.gamma/",descriptors[i,1],"_",descriptors[i,2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
      loglik.mat<-cbind(loglik.mat,fit$loglik.i)
      coln<-c(coln,paste0(descr[i,1],"-",descr[i,2]))
    }
  }
  print(loglik.mat[1:10,])
  colnames(loglik.mat)<-coln
  
  idx<-which(apply(loglik.mat,1,function(v) all(!is.na(v))))
  print(length(idx))
  #  idx<-intersect(idx,which(!is.na(fit$loglik.i)))
  #  print(length(idx))
  
  sumloglik<-apply(loglik.mat[idx,],2,sum)
  print(sumloglik)
  
  idx.1<-intersect(idx,soso$ix[1:62*12])
  sumloglik.1<-apply(loglik.mat[idx.1,],2,sum)
  print(sumloglik.1)
  print(rbind(precip[idx.1],fit$param[idx.1,1]*fit$param[idx.1,2]))
  
  idx.2<-intersect(idx,soso$ix[1:100])
  sumloglik.2<-apply(loglik.mat[idx.2,],2,sum)
  print(sumloglik.2)
  print(precip[idx.2])
  
  print(c(length(idx),length(idx.1),length(idx.2)))
  
  load(file=paste0("2_Travail/20CR/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  #normalize<-length(idx)
  #normalize.1<-length(idx.1)
  #normalize.2<-length(idx.2)
  
  normalize<-sum(score$loglik.gamma[idx])
  normalize.1<-sum(score$loglik.gamma[idx.1])
  normalize.2<-sum(score$loglik.gamma[idx.2])
  
  if (seasonal) seastr<-"seasonal"
  else seastr<-"overall"
  
  png(file=paste0(get.dirstr(k),"compare.loglik.gamma/loglik.gamma.nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".png"),width=14,height=8,units="in",res=72)
  par(mar=c(9,4,4,2),mfrow=c(1,3))
  plot(c(1,ndesc),range(sumloglik-normalize),axes=FALSE,xlab="",ylab="",type="n")
  for (k in 1:length(strs)) {
    points(1:ndesc,sumloglik[(k-1)*ndesc+(1:ndesc)]-normalize,col=k)
    lines(1:ndesc,sumloglik[(k-1)*ndesc+(1:ndesc)]-normalize,col=k)
  }
  #legend(legend=1:3,col=1:3)
  axis(2)
  axis(1,labels=coln,at=1:length(sumloglik),las=3)
  box()
  
  title(paste0("Gamma distribution, all precip>0 ",seastr))
  legend("bottomleft",legend=strs,col=1:length(strs),pch=1)
  abline(v=1:ndesc,lty=2,col=gray(0.5))
  
  plot(c(1,ndesc),range(sumloglik.1-normalize.1),axes=FALSE,xlab="",ylab="",type="n")
  for (k in 1:length(strs)) {
    points(1:ndesc,sumloglik.1[(k-1)*ndesc+(1:ndesc)]-normalize.1,col=k)
    lines(1:ndesc,sumloglik.1[(k-1)*ndesc+(1:ndesc)]-normalize.1,col=k)
  }
  axis(2)
  axis(1,labels=coln,at=1:length(sumloglik.1),las=3)
  box()
  title("Gamma distribution, monthly max (62*12 largest)")
  abline(v=1:ndesc,lty=2,col=gray(0.5))
  
  plot(c(1,ndesc),range(sumloglik.2-normalize.2),axes=FALSE,xlab="",ylab="",type="n")
  for (k in 1:length(strs)) {
    points(1:ndesc,sumloglik.2[(k-1)*ndesc+(1:ndesc)]-normalize.2,col=k)
    lines(1:ndesc,sumloglik.2[(k-1)*ndesc+(1:ndesc)]-normalize.2,col=k)
  }
  axis(2)
  axis(1,labels=coln,at=1:length(sumloglik.2),las=3)
  box()
  title("Gamma distribution, 100 largest precip")
  abline(v=1:ndesc,lty=2,col=gray(0.5))
  
  dev.off()
  
}

compare.loglik.p0.nei<-function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05"){
  descr<-rbind(
    c("sing","sing"),
    c("cel","sing"),
    #c("pers","sing"),
    c("pnei","sing"),
    c("snei","sing"),
    c("rsing","rsing"),
    c("cel","rsing"),
    #c("pers","rsing"),
    c("pnei","rsing"),
    c("snei","rsing"),
    #c("pers","rlsing")
    c("sing","rsing"))
  
  ndesc<-nrow(descr)
  
  precip<-get.precip(nbdays,start,end)
  
  strs<-c("05","1","2","5")
  
  loglik.mat<-NULL
  coln<-NULL
  for (str in strs){
    descriptors<-descr
    descriptors[,1]<-paste0(descriptors[,1],str)
    descriptors[,2]<-paste0(descriptors[,2],str)
    descriptors[substr(descriptors,1,3)=="cel"]<-"cel"
    print(descriptors)
    
    for (i in 1:nrow(descriptors)){
      print(descriptors[i,])
      load(file=paste0(get.dirstr(k),"fit.loglik.p0/",descriptors[i,1],"_",descriptors[i,2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
      loglik.mat<-cbind(loglik.mat,fit$loglik)
      coln<-c(coln,paste0(descr[i,1],"-",descr[i,2]))
    }
  }
  colnames(loglik.mat)<-coln
  
  idx<-which(apply(loglik.mat,1,function(v) all(!is.na(v))))
  print(length(idx))
  sumloglik<-apply(loglik.mat[idx,],2,sum)
  
  load(file=paste0("2_Travail/20CR/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
  print(score$loglik.p0)
  print(summary(score$loglik.p0))
  
  normalize<-sum(score$loglik.p0)
  
  png(file=paste0(get.dirstr(k),"compare.loglik.p0/loglik.p0.nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".png"),width=7,height=8,units="in",res=72)
  par(mar=c(5.5,4,4,2))
  plot(c(1,ndesc),range(sumloglik-normalize),axes=FALSE,xlab="",ylab="",type="n")
  for (k in 1:length(strs)) {
    points(1:ndesc,sumloglik[(k-1)*ndesc+(1:ndesc)]-normalize,col=k)
    lines(1:ndesc,sumloglik[(k-1)*ndesc+(1:ndesc)]-normalize,col=k)
  }
  axis(2)
  axis(1,labels=coln,at=1:length(sumloglik),las=3)
  box()
  if (seasonal) seastr<-"seasonal"
  else seastr<-"overall"
  title(paste0("p0 ",seastr))
  legend("bottomleft",legend=strs,col=1:5,pch=1)
  abline(v=1:ndesc,lty=2,col=gray(0.5))
  dev.off()
  
}

# Trace les parametres des lois empiriques obtenues selon les differents couples d'indicateurs (un parametre de la loi par png)
compare.empir<-function(case,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=FALSE){
  
  if (case==0){ # seulement k=1
    descr<-rbind(
      c("cel","sing"),
      c("snei","sing"),
      c("cel","rsing"),
      c("snei","rsing"),
      c("sing","rsing"))
    compstr<-"compare0"
    kvec<-1
    kstr.save<-kvec
  }
  if (case==1){ # k = 1 ? 3 et enregistrement en amont dans l'arborescence
    descr<-rbind(
      c("cel","sing"),
      c("snei","sing"),
      c("cel","rsing"),
      c("snei","rsing"),
      c("sing","rsing"))
    
    neiseq<-c("05")
    compstr<-"compare1"
    kvec<-1:3
    kstr.save<-NULL
  }
  
  
  print(descr)
  
  param.list<-list()
  for (i in 1:nrow(descr)){ # pour chaque couple d'indicateurs
    param.i.list<-list()
    param.i.list$param<-list()
    param.i.list$descr<-list()
    
    for (nn in 1:length(kvec)){ # pour chaque k
      descriptor1<-descr[i,1]
      descriptor2<-descr[i,2]
      k<-kvec[nn]
      if (descriptor1!="cel") descriptor1<-paste0(descriptor1,"05") # si on n'est pas sur la celerite, on ajoute un "05" aux indicateurs
      if (descriptor2!="cel") descriptor2<-paste0(descriptor2,"05")
      print(c(descriptor1,descriptor2))
      print(paste0(get.dirstr(k),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
      load(file=paste0(get.dirstr(k),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata")) # Import des parametres de la loi empirique
      
      descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE) # Import des descripteurs
      descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE)
      
      param.i.list$param[[nn]]<-param # remplissage des listes
      param.i.list$descr[[nn]]<-cbind(descr1,descr2)
    }
    param.list[[i]]<-param.i.list # liste contenant les deux precedentes listes pour le couple d'indicateur donne
  }
  
  coln<-colnames(param.list[[1]]$param[[1]])
  print(coln)
  
  zlims<-matrix(NA,ncol=2,nrow=length(coln)) # deux colonnes, une pour min une pour max
  for (p in 1:length(coln)){ # pour chaque parametre, i.e. chaque ligne de la matrice zlims
    tmp<-NULL
    for (i in 1:nrow(descr)){ # pour chaque couple d'indicateurs
      for (nn in 1:length(kvec)){ # pour chaque k
        tmp<-range(c(tmp,range(param.list[[i]]$param[[nn]][,p],na.rm=TRUE))) # plage de valeur du parametre de la loi empirique
      }
    }
    zlims[p,]<-tmp # on remplit la ligne de la matrice zlims
  }
  print(zlims)
  
  
  for (p in 1:length(coln)){ # pour chaque parametre de la loi empirique
    #graphics.off()
    print(paste0(get.dirstr(kstr.save),"compare.empir",get.CVstr(CV),"/",compstr,"_",coln[p],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".png"))
    png(file=paste0(get.dirstr(kstr.save),"compare.empir",get.CVstr(CV),"/",compstr,"_",coln[p],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".png"),width=nrow(descr)*3.5,height=4*length(kvec),units="in",res=72)
    par(mfrow=c(length(kvec),nrow(descr))) # on divise la fenetre graphique par le nombre de couples d'indicateurs
    
    zlim<-zlims[p,]
    
    for (nn in 1:length(kvec)){ # pour chaque k
      for (i in 1:nrow(descr)){ # pour chaque couple d'indicateurs
        
        descriptor1<-descr[i,1]
        descriptor2<-descr[i,2]
        
        param.pik<-param.list[[i]]$param[[nn]][,p]
        plot(param.list[[i]]$descr[[nn]][,1],param.list[[i]]$descr[[nn]][,2],col=getcol(param.pik,range=zlim),xlab=descriptor1,ylab=descriptor2) # plot du parametre colorie en fonction des deux indicateurs en x et en y
        addscale(zlim,r=2) # legende de la gamme de couleur
        title(paste0(coln[p],": ",descriptor2," vs. ",descriptor1))
        legend("topright",legend=paste0(round(min(param.pik,na.rm=TRUE),2),"-",round(max(param.pik,na.rm=TRUE),2)),bty="n") # legende de l'amplitude du parametre dans le plan
      }
    }
    dev.off()
  }
  
}

compute_density<-function(descriptor1,descriptor2,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE){
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize)
  
  idx.notna<-unique(c(!is.na(descr1),!is.na(descr2)))
  descr1<-descr1[idx.notna]
  descr2<-descr2[idx.notna]
  
  print(sd(descr1))
  print(sd(descr2))
  
  bandwidth<-c(sd(descr1)/5,sd(descr2)/5)
  gridsize<-c(200,200)
  dens<-bkde2D(cbind(descr1,descr2),bandwidth=bandwidth,gridsize = gridsize,range.x = list(c(min(descr1)-sd(descr1),max(descr1)+sd(descr1)),c(min(descr2)-sd(descr2),max(descr2)+sd(descr2))))
  dens$prob<-dens$fhat*(dens$x1[2]-dens$x1[1])*(dens$x2[2]-dens$x2[1])
  
  graphics.off()
  png(file=paste0(get.dirstr(k),"compute_density/image_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),".png"),width=600,heigh=600)
  par(mar=c(5.1,5.1,4.1,2.1))
  l<-max(c(max(dens$x1)-min(dens$x1),max(dens$x2)-min(dens$x2)))
  limx<-c(min(dens$x1),min(dens$x1)+l)
  limy<-c(min(dens$x2),min(dens$x2)+l)
  print("lims")
  print(lims)
  image.plot(x=dens$x1,y=dens$x2,z=dens$fhat,main=paste0("Density bkde2D, ",getZstr(k),", ",dist,sep=""),xlab=descriptor1,ylab=descriptor2,cex.main=1.7,cex.lab=1.7,cex.axis=1.5,legend.width=2,legend.cex = 2,xlim=limx,ylim=limy)
  #points(descr1,descr2,cex=0.1)
  print(dim(dens$fhat))
  dev.off()
  
  png(file=paste0(get.dirstr(k),"compute_density/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),".png"),width=600,heigh=600)
  par(mar=c(5.1,5.1,4.1,2.1))
  plot(descr1,descr2,main=paste0("Scatter plot, ",getZstr(k),", ",dist,sep=""),xlab=descriptor1,ylab=descriptor2,cex.main=1.7,cex.lab=1.7,cex.axis=1.5,xlim=limx,ylim=limy,cex=0.1)
  #points(descr1,descr2,cex=0.1)
  dev.off()
  
  save(dens,file=paste0(get.dirstr(k),"compute_density/density_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),".Rdata"))
  
  #case with kde2d: tres bruite
  #dens<-kde2d(descr1, descr2, h=bandwidth, n = gridsize, lims = c(min(descr1)-1.5*bandwidth[1], max(descr1)+1.5*bandwidth[1],min(descr2)-1.5*bandwidth[2], max(descr2)+1.5*bandwidth[2]))
  #quartz()
  #image.plot(x=dens$x,y=dens$y,z=dens$z,main=paste0("Density kde2d, ",getZstr(k),", ",dist,sep=""),xlab=descriptor1,ylab=descriptor2,cex.main=1.7,cex.lab=1.7,cex.axis=1.5,legend.width=2,legend.cex = 2)
  #print(dim(dens$z))
}

# Calcul du nouveau score RMSE_new_bon (formule de base, donc ok mais long)
compute_RMSE_new_bon<-function(k,start="1950-01-01",end="2011-12-31",rean){
  
  # Import geopotentiel et score RMSE
  load.nc()
  dat<-getdata(k,start,end,rean)
  N<-dim(dat)[3]
  
  # Calcul nouveau score
  moy_geo<-apply(dat,3,mean)
  moy_geo <- rep(moy_geo,each=dim(dat)[1]*dim(dat)[2])
  dim(moy_geo) <- dim(dat)
  norm <- dat - moy_geo
  dist.list<-list()
  
  for(i in 1:(N-1)){
    if (i%%50==0) print(i)
    normj<- norm[,,(i+1):N]
    normi<-array(norm[,,i],dim(normj))
    dist.list[[i]] <- apply(normi-normj,3,function(v) sqrt(mean(v^2)))
  }
  dist.list
}

# Ancienne fonction fit.empir
fit.empir.old<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean){
  
  print(paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  precip<-get.precip(nbdays,start,end)
  
  descr<-NULL
  for (i in 1:length(descriptors)) descr<-cbind(descr,get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize,rean)) # recupere les indicateurs
  
  if (nrow(descr) !=length(precip)) stop("Probleme taille des donnees")
  
  N<-length(precip)
  
  param<-matrix(NA,ncol=7,nrow=N)
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){
      tmp<-get.closest(i,descr,precip,CV,nbdays,radtype)
      pi<-precip[tmp$idx]
      ni<-length(pi)
      ni1<-sum(pi>0)
      parami<-c(tmp$nb,1-ni1/ni,mean(pi[pi>0]),sd(pi[pi>0]),skewness(pi[pi>0]),kurtosis(pi[pi>0]),tmp$radius) # remplissage de la ligne: nbre de voisins, proba jour sec, moyenne et ecart-type des pluies non nulles, coefficient d'assymetrie et d'applatissement, rayon du cercle
    }
    else parami<-rep(NA,7) # si un des deux indicateurs est nul, ligne du tableau parametre vide
    param[i,]<-parami
  }
  colnames(param)<-c("nbnei","p0","mean","sd","skewness","kurtosis","radius")
  
  save(param,file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
}

# Ajustement d'une loi de Pareto selon le voisinage au sens des indicateurs
fit.loglik.egp<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean){
  
  print(paste0(get.dirstr(k,rean),"fit.loglik.egp",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  precip<-get.precip(nbdays,start,end)
  
  descr<-NULL
  for (i in 1:length(descriptors)) descr<-cbind(descr,get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize,rean))
  
  if (nrow(descr) !=length(precip)) stop("Probleme taille des donnees")
  
  N<-length(precip)
  
  param<-matrix(NA,ncol=3,nrow=N)
  loglik.i<-rep(NA,N)
  opt<-list()
  
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){ # si les indicateurs sont non nuls
      tmp<-get.closest(i,descr,precip,CV,nbdays,radtype)
      pi<-precip[tmp$idx]
      ni<-length(pi)
      ni0<-sum(pi>10^-10)
      pi0<-pi[pi>10^-10]
      init=c(0.9, fpot(pi0,quantile(pi0,0.8),shape=0.05)$param) # initialisation de la loi de Pareto
      opt[[i]]<-egp2.fit(pi0, model=1, method="mle", init=init,rounded=0,CI=FALSE,plots=FALSE) # Optimisation des parametres kappa, sigma et xi
      param[i,]<-opt[[i]]$fit$mle
    }
  }
  
  colnames(param)<-c("kappa","sigma","xi") 
  
  for (i in 1:N){
    if (precip[i]>0 & !is.na(descr[i,1]) & !is.na(descr[i,2])) loglik.i[i]<- degp2(precip[i], kappa=param[i,"kappa"], sigma=param[i,"sigma"], xi=param[i,"xi"], type=1, log=TRUE) # log likelihood de la valeur de pluie si non nulle => probabilite d'avoir cette pluie selon la loi de Pareto ajustee
  }  
  fit<-list(param=param,opt=opt,loglik.i=loglik.i)
  
  save(fit,file=paste0(get.dirstr(k,rean),"fit.loglik.egp",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
}

# Ajustement d'une loi de Pareto selon le voisinage au sens de l'analogie classique
fit.loglik.egp.A<-function(rad,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  precip<-get.precip(nbdays,start,end)
  
  N<-length(precip)
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[[rad]]
  
  param<-matrix(NA,ncol=3,nrow=N)
  loglik.i<-rep(NA,N)
  opt<-list()
  for (i in 1:N){
    if (i %%50==0) print(i)
    pi<-precip[setdiff(nei[[i]],(i-nbdays+1):(i+nbdays-1))]
    pi<-pi[which(!is.na(pi))]
    ni<-length(pi)
    ni0<-sum(pi>10^-10)
    pi0<-pi[pi>10^-10]
    init=c(0.9, fpot(pi0,quantile(pi0,0.8),shape=0.05)$param)
    opt[[i]]<-egp2.fit(pi0, model=1, method="mle", init=init,rounded=0,CI=FALSE,plots=FALSE)
    param[i,]<-opt[[i]]$fit$mle
  }
  
  colnames(param)<-c("kappa","sigma","xi") 
  
  for (i in 1:N){
    if (precip[i]>0) loglik.i[i]<- degp2(precip[i], kappa=param[i,"kappa"], sigma=param[i,"sigma"], xi=param[i,"xi"], type=1, log=TRUE)
  }  
  fit<-list(param=param,opt=opt,loglik.i=loglik.i)
  
  save(fit,file=paste0(get.dirstr(k,rean),"fit.loglik.egp-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Ajustement d'une loi gamma selon le voisinage au sens des indicateurs
fit.loglik.gamma<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean){
  
  print(paste0(get.dirstr(k,rean),"fit.loglik.gamma",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  precip<-get.precip(nbdays,start,end)
  
  descr<-NULL
  for (i in 1:length(descriptors)) descr<-cbind(descr,get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize,rean))
  
  if (nrow(descr) !=length(precip)) stop("Probleme taille des donnees")
  
  N<-length(precip)
  
  nloglik.fct<-function(param,x){
    -sum(dgamma(x, shape=param[2], scale = param[1], log = TRUE) )
  }
  
  param<-matrix(NA,ncol=2,nrow=N)
  loglik.i<-rep(NA,N)
  opt<-list()
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){ # si les indicateurs sont non nuls
      tmp<-get.closest(i,descr,precip,CV,nbdays,radtype)
      pi<-precip[tmp$idx]
      ni<-length(pi)
      ni0<-sum(pi>0)
      pi0<-pi[pi>0]
      mean0<-mean(pi0)
      var0<-var(pi0)
      scale<-var0/mean0
      shape<-mean0^2/var0
      opt[[i]]<-optim(c(scale,shape),nloglik.fct,method="L-BFGS-B",x=pi0,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) # Optimisation de scale and shape par maximum de vraisemblance
      names(opt[[i]]$par)<-c("scale","shape")
      param[i,]<-opt[[i]]$par
    }
  }
  
  colnames(param)<-c("scale","shape")
  
  for (i in 1:N){
    if (precip[i]>0) loglik.i[i]<- -nloglik.fct(param[i,],precip[i]) # log likelihood de la valeur de pluie si non nulle => probabilite d'avoir cette pluie selon la loi gamma ajustee
  }  
  
  fit<-list(param=param,opt=opt,loglik.i=loglik.i)
  
  save(fit,file=paste0(get.dirstr(k,rean),"fit.loglik.gamma",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
}

# Ajustement d'une loi gamma selon le voisinage au sens de l'analogie classique
fit.loglik.gamma.A<-function(rad,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  precip<-get.precip(nbdays,start,end)
  
  N<-length(precip)
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[[rad]]
  
  nloglik.fct<-function(param,x){
    -sum(dgamma(x, shape=param[2], scale = param[1], log = TRUE) )
  }
  
  param<-matrix(NA,ncol=2,nrow=N)
  loglik.i<-rep(NA,N)
  opt<-list()
  for (i in 1:N){
    if (i %%50==0) print(i)
    pi<-precip[setdiff(nei[[i]],(i-nbdays+1):(i+nbdays-1))]
    pi<-pi[which(!is.na(pi))]
    ni<-length(pi)
    ni0<-sum(pi>0)
    pi0<-pi[pi>0]
    mean0<-mean(pi0)
    var0<-var(pi0)
    scale<-var0/mean0
    shape<-mean0^2/var0
    #print(pi)
    #print(pi0)
    #print(precip[nei[[i]]])
    opt[[i]]<-optim(c(scale,shape),nloglik.fct,method="L-BFGS-B",x=pi0,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) 
    names(opt[[i]]$par)<-c("scale","shape")
    param[i,]<-opt[[i]]$par
  }
  
  colnames(param)<-c("scale","shape")
  
  for (i in 1:N){
    if (precip[i]>0) loglik.i[i]<- -nloglik.fct(param[i,],precip[i])
  }  
  
  fit<-list(param=param,opt=opt,loglik.i=loglik.i)
  
  save(fit,file=paste0(get.dirstr(k,rean),"fit.loglik.gamma-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Taille de la fenetre spatiale d'analogie (fonction Simon)
getinfo_window=function(k){ 
  #lat et lon correspondent aux degrés de longitude et latitude pour lesquels on a les données des géopot 500 et 1000
  lon<-nc[[1]]$dim$lon$vals # de -30 à 50 deg E tranches de 2 deg -> 41 valeurs
  lat<-nc[[1]]$dim$lat$vals #de 24 - 72 deg N tranches de 2 deg ->25 valeurs
  #time=data500$dim$time$vals # 58804 valeurs par tranches de 24h (du 01-01-1851 au 31-12-2011), valeurs à 9h chaque jour
  #paramètre fenêtre d'analogie
  
  c_lon<-6 #centre de la fenêtre longitude
  c_lat<-44 #centre de la fenêtre latitude
  
  if(k==1) {#500hPA
    d_lon<-32
    d_lat<-16
  }
  if (k==2){ #1000 hPA
    d_lon<-16
    d_lat<-8
  }
  infolon<-c(which(lon==c_lon-d_lon/2),which(lon==c_lon+d_lon/2)-which(lon==c_lon-d_lon/2)+1) #1er point de grille et nbre de points de grille
  infolat<-c(which(lat==c_lat-d_lat/2),which(lat==c_lat+d_lat/2)-which(lat==c_lat-d_lat/2)+1)
  
  return(rbind(infolon,infolat))#matrice 2x2 
  
} 

# Si rayon fixe, calcule le rayon du cercle voisin en moyennant les ecart-types des deux descripteurs
get.radius<-function(descr,norm){
  sdvec<-apply(descr,2,sd,na.rm=TRUE)
  radius<-mean(sdvec)/norm
}

# plot du copule colorie par la valeur de alpha, comme dans Faranda (2020)
plot.alpha.faranda<-function(rho,qua){
  
  # Copule gaussienne
  z<-rmvnorm(2*10^3,mean=c(0,0),sigma=matrix(c(1,rho,rho,1),ncol=2))
  xvec<-pnorm(z[,1])
  yvec<-pnorm(z[,2])
  plot(xvec,yvec,xlim=c(0,1),ylim=c(0,1))
  
  # Alpha
  alpha <- NULL
  for(i in 1:2000){
    print(i)
    dist.x <- -log(abs(z[i,1]-z[,1]))
    qua.x <- quantile(dist.x,qua)
    
    dist.y <- -log(abs(z[i,2]-z[,2]))
    qua.y <- quantile(dist.y,qua)
    
    ana.x <- which(dist.x>qua.x & dist.x!=Inf)
    ana.y <- which(dist.y>qua.y & dist.y!=Inf)
    ana <- which(dist.x>qua.x & dist.y>qua.y & dist.x!=Inf & dist.y!=Inf)
    
    alpha[i] <- length(ana)/length(ana.x)
  }
  
  # Plot
  plot(xvec,yvec,xlim=c(0,1),ylim=c(0,1),col=getcol(alpha))
  addscale(alpha)
}

# Ancienne fonction plot.empir
plot.empir.old<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,rean,empir=TRUE,obs=FALSE){
  
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  if(!empir) precip <- get.precip(nbdays,start,end)
  
  print(paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata"))
  load(file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata")) # importation des parametres de la loi empirique
  param<-param[,colnames(param)!="kurtosis"]
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE,rean)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE,rean)
  
  if(empir){
    print(paste0(get.dirstr(k,rean),"plot.empir",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"))
    png(file=paste0(get.dirstr(k,rean),"plot.empir",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"),width=10,height=7,units="in",res=72)
    par(mfrow=c(2,3))
    for (i in 1:ncol(param)){
      param.i<-param[,i]
      #q1<-quantile(param.i,0.005,na.rm=TRUE)
      #q2<-quantile(param.i,0.995,na.rm=TRUE)
      #param.i[param.i<q1]<-q1
      #param.i[param.i>q2]<-q2
      plot(descr1,descr2,col=getcol(param.i),xlab=descriptor1,ylab=descriptor2)
      addscale(param.i)
      ind.extr <- get.ind.extr(nbre = 62,nbdays = 3)
      ind.extr <- as.vector(ind.extr)
      points(descr1[ind.extr],descr2[ind.extr],pch=19)
      title(colnames(param)[i])
    }
  }
  
  if(obs){ 
    print(paste0(get.dirstr(k,rean),"plot.empir",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_obs_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"))
    png(file=paste0(get.dirstr(k,rean),"plot.empir",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_obs_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"),width=10,height=7,units="in",res=72)#,width=10,height=7,units="in",res=72)
    layout(matrix(1:2,1,2))
    
    plot(descr1,descr2,col=getcol(precip*nbdays),xlab=descriptor1,ylab=descriptor2)#.,pch=19,cex=0.6)
    addscale(precip*nbdays)
    title("Observed precipitation")
    
    plot(descr1,descr2,col=as.numeric(precip==0)+1,xlab=descriptor1,ylab=descriptor2)#,pch=19,cex=0.6)
    legend("topleft",c("Rain","No rain"),col=c(1,2),bty="n",pch=19)
    title("Rain - No rain")
  }
  
  dev.off()
}

# Trace les boxplot des interquantiles des scores pour les sequences seches et les sequences des pluies fortes
plot.interquantile <- function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  # Importation des scores et des precip
  load(file=paste0(get.dirstr(k,rean),"save.score.A/score_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  precip <- get.precip(nbdays,start,end)
  
  # Calcul des quantiles des scores
  quant <- c(0.001,0.002,0.005,0.01,0.02,0.05,0.1)
  
  qua <- lapply(score,function(v) quantile(v/nbdays,probs = quant)) # on divise par nbdays pour ramener le score entre 0 et 1
  #inter.qua <- lapply(qua, function(v) c(v[-1]-v[-length(v)]))
  
  #quant.norm <- 1/c(quant[-1]-quant[-length(quant)]) # on ramene l'espace interquantile a une valeur normalisee sur 100% des donnees
  #inter.qua.norm <- lapply(inter.qua, function(v) v*quant.norm)
  
  # Pluies qui nous interessent
  soso <- sort(precip,index.return=TRUE,decreasing=TRUE)
  ind <- vector("list",3)
  names(ind) <- c("Toutes sequences","Sequences seches","62 sequences pluie forte")
  ind[[1]] <- 1:length(precip)
  ind[[2]] <- which(precip == 0)
  ind[[3]] <- soso$ix[1:62]
  
  # Boxplot
  tab <- do.call(rbind,qua)
  #colnames(tab) <- paste0(colnames(tab)," - ",names(qua[[1]])[1:ncol(tab)])
  ylim <- c(0,0.6)
  
  for(i in 1:length(ind)){
    tab.plot <- tab[ind[[i]],]
    
    filename <- paste0(get.dirstr(k,rean),"plot.interquantile/",i,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png")
    png(file=filename,width=18,height=8,units="in",res=72)
    boxplot(tab.plot, ylim = ylim, main = paste0("Interquantiles ",dist," normalises - ",names(ind)[[i]]))
    abline(v=1:ncol(tab),lty=2,col=gray(0.5))
    
    graphics.off()
  }
  
}

# Trace la repartition des parametres de la loi de Pareto dans le plan des indicateurs
plot.loglik.egp<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=FALSE){
  
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  print(c(descriptor1,descriptor2))
  
  load(file=paste0(get.dirstr(k),"fit.loglik.egp",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  load(file=paste0(get.dirstr(k),"save.stat.egp",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  png(file=paste0(get.dirstr(k),"plot.loglik.egp",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".png"),width=24,height=4,units="in",res=72)
  par(mfrow=c(1,6))
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE)
  
  for (i in 1:6){
    if (i==1) {
      nam<-"kappa"
      val<-fit$param[,"kappa"]
    }
    if (i==2) {
      nam<-"sigma"
      val<-fit$param[,"sigma"]
    }
    if (i==3){
      nam<-"xi"
      val<-fit$param[,"xi"]
    }
    if (i==4){
      nam<-"mean"
      val<-statmat[,"mean"]
    }
    if (i==5){
      nam<-"sd"
      val<-statmat[,"sd"]
    }
    if (i==6){
      nam<-"skewness"
      val<-statmat[,"skewness"]
    }
    plot(descr1,descr2,col=getcol(val),xlab=descriptor1,ylab=descriptor2)
    addscale(val)
    title(nam)
  }
  dev.off()
  
}

# Trace la repartition des parametres de la loi gamma dans le plan des indicateurs
plot.loglik.gamma<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=FALSE){
  
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  print(c(descriptor1,descriptor2))
  
  load(file=paste0(get.dirstr(k),"fit.loglik.gamma",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  png(file=paste0(get.dirstr(k),"plot.loglik.gamma",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".png"),width=20,height=4,units="in",res=72)
  par(mfrow=c(1,5))
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE)
  
  for (i in 1:5){
    if (i==1) {
      nam<-"scale"
      val<-fit$param[,"scale"]
    }
    if (i==2) {
      nam<-"shape"
      val<-fit$param[,"shape"]
    }
    if (i==3){
      nam<-"mean"
      val<-fit$param[,"scale"]*fit$param[,"shape"]
    }
    if (i==4){
      nam<-"sd"
      val<-fit$param[,"scale"]*sqrt(fit$param[,"shape"])
    }
    if (i==5){
      nam<-"skewness"
      val<-2/sqrt(fit$param[,"shape"])
    }
    plot(descr1,descr2,col=getcol(val),xlab=descriptor1,ylab=descriptor2)
    addscale(val)
    title(nam)
  }
  dev.off()
  
}

# Enregistrement du fichier precip pour cumul sur une ou plusieurs journees (sur plusieurs journees: mean et pas sum mais ca revient au meme)
save.precip<-function(nbdays,start="1950-01-01",end="2011-12-31"){
  if (nbdays==1) precip<-make.precip1(start,end)
  else {
    precip<-get.precip(1,start,end)
    precip<-rollapply(precip,width=nbdays,FUN=mean)
    precip[precip<0]<-0 #sometimes -10^-15 instead of 0
  }
  write.csv(precip,file=paste0("/Users/juliette/DataLTHE/METEO-FRANCE/DAILY/Isere@Grenoble/Isere@Grenoble_cum",nbdays,"day_",start,"_",end,".csv"),row.names=FALSE)
}

# Enregistrement des dates dans le membre utilise
savedates<-function(){
  dates<-substr(as.POSIXct(nc[[1]]$dim$time$vals*3600,origin='1800-01-01 00:00',tz="GMT"),1,10) #format "YYYY-MM-DD"
  save(dates,file=paste0("2_Travail/20CR/Membre_",member,"/dates.Rdata"))
}

#savedates<-function(start,end){
#  dates<-substr(as.POSIXct(nc[[1]]$dim$time$vals*3600,origin='1800-01-01 00:00',tz="GMT"),1,10) #format "YYYY-MM-DD"
#  dates[which(d==start):which(d==end)]
#  save(dates,file=paste0("compute_dist/",dist,"_member",member,"_Z500+1000_",start,"_",end,".Rdata"))
#}

scatterplot.criteria<-function(neistr,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=FALSE){
  
  dd<-dir(paste0(get.dirstr(k),"fit.empir/"))
  dd<-dd[grep(paste0("ing",neistr,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),".Rdata"),dd)]
  print(dd)
  ldd<-strsplit(dd,"_")
  
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE)
  id0<-which(precip==0)
  
  for (i in 1:length(ldd)){
    descriptor1<-ldd[[i]][1]
    descriptor2<-ldd[[i]][2]
    
    load(file=paste0(get.dirstr(k),"fit.empir/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),".Rdata"))
    
    descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize)
    descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize)
    descr<-cbind(descr1,descr2)
    N<-nrow(descr)
    
    for (case in 0:1){
      png(file=paste0(get.dirstr(k),"scatterplot.criteria/",case,"_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),".png"),width=10,height=10,units="in",res=72)
      
      if (case==0) {
        id<-id0
        tid<-rep(0,length(id0))
      }
      if (case==1) {
        id<-soso$ix[1:20]
        tid<-1:20
      }
      
      plot(descr1,descr2,col=getcol(param[,"nbnei"]),xlab=descriptor1,ylab=descriptor2)
      title(paste0(descriptor2," vs. ",descriptor1))
      shadowtext(descr1[id],descr2[id],labels=tid,col="black",bg="white",font=2)
      addscale(param[,"nbnei"])
      
      dev.off()
    }
  }
  
}

# Amelioration possible de l'analogie classique pour modeliser les precipitations
test.ana <- function(){
  
  # Import precip et jour d'interet
  precip=get.precip(1)
  precip[2302]
  
  # Distribution des analogues
  ana=get.ana(date=getdates()[2302],rank = 1:224,k = 1,dist = "TWS",ref = start,nbdays = 1,start,end,rean)
  plot.ecdf(precip[ana$ind[1:112]])
  abline(v=precip[2302])
  
  # Analogues les plus proches selon la moyenne du geopotentiel
  geo=getdata(k = 1,start,end,rean)
  tmp=apply(geo,3,mean)
  ran=sort(abs(tmp[ana$ind]-tmp[2302]),index.return=T)$ix
  ana.new=ana$ind[ran[1:112]]
  plot.ecdf(precip[ana.new])
  abline(v=precip[2302])
  
  # crps avant apres
  crps_sample(precip[2302],precip[ana$ind[1:112]])
  crps_sample(precip[2302],precip[ana.new])
  
  # distance moyenne des analogues à la moyenne du geopotentiel cible pour les pluies croissantes
  for(i in 1:22645){
    tmp1[i]=mean(tmp[i]-tmp[nei[[i]]])
  }
  
  ind=sort(precip,index.return=T)$ix
  plot(sort(precip),abs(tmp1[ind]))
  points(sort(precip[ind1]),abs(tmp1[ind])[ind1],col="red")
  # => plus on va vers les extrêmes, plus les écarts à la moyenne diminuent
}

# Test sur la reconstitution des indices analogues de chaque date a partir de la liste
test.getdist4i<-function(N,len=NULL){
  l<-list();
  for (i in 1:(N-1)) l[[i]]<-paste(i,(i+1):N)
  ul<-unlist(l)
  for (i in 1:N) print(ul[selind_season(i,len,ul)])
}

# Test sur la reconstitution des indices analogues de chaque date a partir de la liste et saisonnalite
test.getdist4i_season<-function(N){
  l<-list();
  for (i in 1:(N-1)) l[[i]]<-paste(i,(i+1):N)
  ul<-unlist(l)
  for (i in 1:N) print(getdist4i_season(i,ul,N,3,10))
}














