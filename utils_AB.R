library(ncdf4) # NetCDF
library(fBasics) #timPalette
library(ismev) #gp.fit
library(splancs) #inpip
library(fields) #Tps
library(KernSmooth) #bkde2D
library(plotrix) #draw.circle
library(mev) #egp2.fit
library(akima) #interp: bilinear interpolation
library(zoo) #rollapply
library(RANN) #nn2
library(scoringRules) #crps
library(evd) #fpot
library(gpclib) #gpc.poly
library(rgdal) #readOGR
library(tidyr) # gather
library(ggplot2) #ggplot
library(ggpubr) #ggarrange
library(deSolve) #ode
library(plot3D) #lines3D
library(PACBO) # runiform_ball
library(rgeos) # pour maptools
library(maptools) # wrld_simpl
library(RColorBrewer) # brewer.pal

# Fonctions graphiques
addcircle<-function(radius){
  usr<-par("usr")
  insetx <- 0.1 * (usr[2L] - usr[1L])
  insety <- 0.1 * (usr[4L] - usr[3L])
  draw.circle(x=usr[2L] - insetx, y=usr[3L] + insety, radius=radius, lty=2)
}
addscale<-function(vec,r=2,legend="",centered=FALSE,rev=FALSE){
  
  usr<-par("usr")
  insetx <- 0.05 * (usr[2L] - usr[1L])
  insety <- 0.1 * (usr[4L] - usr[3L])
  xmin<-usr[1L]*1.05 + insetx
  ymin<-usr[4L] - insety
  width <- 0.4 * (usr[2L] - usr[1L])
  heigh <- 0.05 * (usr[4L] - usr[3L])
  shift <- 0.05 * (usr[4L] - usr[3L])
  seg <- 0.02 * (usr[4L] - usr[3L])
  
  col=colpalette(rev)
  
  if (centered) {
    m<-max(abs(vec),na.rm=TRUE)
    minv<- -m
    maxv<- m
  }
  else{
    minv<-min(vec,na.rm=TRUE)
    maxv<-max(vec,na.rm=TRUE)
  }
  
  rect(xmin+seq(0,(0.99*width),length=100),rep(ymin,100),xmin+width/100+seq(0,(0.99*width),length=100),rep(ymin+heigh,100),col=col,border=col)
  text(xmin,ymin-shift,round(minv,r),cex=1)
  text(xmin+width,ymin-shift,round(maxv,r),cex=1)
  text(xmin+width/2,ymin-shift,round((minv+maxv)/2,r),cex=1)
  rect(xmin,ymin,xmin+width,ymin+heigh)
  segments(xmin,ymin,xmin,ymin-seg)
  segments(xmin+width/2,ymin,xmin+width/2,ymin-seg)
  segments(xmin+width,ymin,xmin+width,ymin-seg)
  text(xmin+width/2,ymin+heigh,labels=legend,pos=3)
}
colpalette<-function(rev=FALSE){
  col<-timPalette(113)[10:110]
  if (rev) col<-rev(col)
  col
}
getcol<-function(vec=NULL,range=NULL,transparent=FALSE,centered=FALSE,rev=FALSE) {
  if (centered) {
    m<-max(abs(vec),na.rm=TRUE)
    minv<- -m
    maxv<- m
  }
  else{
    minv<-min(vec,na.rm=TRUE)
    maxv<-max(vec,na.rm=TRUE)
  }
  colvec<-colpalette(rev)
  if (is.null(vec)) col<-colvec
  else if (is.null(range)) col<-colvec[99*(vec-minv)/(maxv-minv)+1] # Attribution d'un ratio entre 0 et 100 pour attribution d'une couleur a chaque valeur du vecteur
  else {
    vec[vec>range[2]]<-range[2]
    vec[vec<range[1]]<-range[1]
    col<-colvec[99*(vec-min(range,na.rm=TRUE))/(max(range,na.rm=TRUE)-min(range,na.rm=TRUE))+1]
    #heat.colors(100)[99*(vec-min(vec))/(max(vec)-min(vec))+1]
  }
  
  if (transparent) {
    col<-apply(col2rgb(col),2,function(v) rgb(v[1],v[2],v[3],alpha=127, maxColorValue=255))
  }
  col[is.na(col)]<-gray(0.5)
  col
}

# Compare notre BV au BV Isere@StEgreve EDF
compare.bv.edf<-function(){
  
  # Import nouveau BV
  poly <- readOGR("2_Travail/Rresults/compare.bv.edf", "DTG_bv_isere_stegreve")
  bv_edf <- poly@polygons[[1]]@Polygons[[1]]@coords/1000
  
  # Export nouveau BV
  write.csv(bv_edf,file="2_Travail/Rresults/compare.bv.edf/border_Isere@Grenoble_EDF.csv",row.names = F)
  
  # Carte des deux BVs
  pdf(file=paste0("2_Travail/Rresults/compare.bv.edf/map_compare_bv.pdf"),width = 7.5,height = 7.5)
  load(file=paste("2_Travail/Data/Carto/griddata_1x1_IsereSavoieHautesAlpes.Rdata",sep=""))
  Fx<-griddata$Fx
  Fy<-griddata$Fy
  Fz<-griddata$Fz*1000
  image.plot(Fx,Fy,Fz,col=gray(seq(0.1,0.99,length=100)),xlab="X (km) - Lambert II extended",ylab="Y (km) - Lambert II extended",legend.line=-2.3, cex.axis=1.3, cex.lab=1.3)
  bord<-read.csv("2_Travail/Data/Carto/border_Isere@Grenoble.csv",sep=",")
  ui<-unique(bord[,1])
  a1<-as(bord[bord[,1]==unique(bord[,1])[1],2:3],"gpc.poly")
  a2<-as(bord[bord[,1]==unique(bord[,1])[2],2:3],"gpc.poly")
  a<-union(a1,a2)
  for (i in ui[-(1:2)]) {
    a1<-as(bord[bord[,1]==i,2:3],"gpc.poly")
    a<-union(a1,a)
  } 
  lines(a@pts[[1]]$x,a@pts[[1]]$y,col="white",lwd=3)
  lines(bv_edf,col="red",lwd=3)
  graphics.off()
  
}

# Calcule et trace les CRPSS pour analogie indicateurs pour les differents couples d'indicateurs et differents rayons 
compare.crps<-function(which="",k=NULL,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="",standardize=TRUE,CV=TRUE,rean,comp=F){
  
  # Indicateurs a comparer
  descr<-list(
    c("persnei","celnei"),
    c("celnei","singnei"),
    c("persnei","singnei"),
    c("celnei","rsingnei"),
    c("persnei","rsingnei"),
    c("singnei","rsingnei")
    #c("TWSgeo","rsingnei"),
    #c("A05","A05")
  )
  
  threeday <- rep(F,6)
  ndesc <- ifelse(comp,length(descr)-1,length(descr))
  namdescr <- lapply(descr,nam2str)
  
  # Pluies
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE) # Classement des pluies de la plus forte a la plus faible avec indices
  
  # Graphique pour un rayon donne (base)
  if (which=="") { 
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<-k
    kstr.save<-paste0("_k",k)
    distance <-dist
    diststr.save<-dist
    legs<-"05"
    style <- 1:length(legs)
    type  <- rep(1,length(legs))
    dirstr.save<-get.dirstr(k,rean)
  }
  # Graphique pour les differents rayons
  if (which=="rad"){ 
    rads<-paste0("nrn",c("02","05"))
    radstr.save<-""
    kvec<-k
    kstr.save<-paste0("_k",k)
    distance <-dist
    diststr.save<-dist
    legs<-c("02","05")
    style <- 1:length(legs)
    type  <- rep(1,length(legs))
    dirstr.save<-get.dirstr(k,rean)
  }
  # Graphique pour les differentes distances (pour un rayon donne)
  if (which=="dist") { 
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<-k
    kstr.save<-paste0("_k",k)
    distance <-c("TWS","RMSE")
    diststr.save<-dist # pour stocker dist et l'utiliser pour le chemin du graphique, meme apres la boucle
    legs<-c("TWS","RMSE")
    style <- 1:length(legs)
    type  <- rep(1,length(legs))
    dirstr.save<-get.dirstr(k,rean)
  }
  # Graphique pour les differents geopotentiels
  if (which=="k") {
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<- c(1,2,4)
    kstr.save<-""
    distance <-dist
    diststr.save<-dist
    legs<-kvec
    style <- 1:length(legs)
    type  <- rep(1,length(legs))
    dirstr.save<-get.dirstr(rean = rean)
  }
  # Graphique pour les differents geopotentiels et distances
  if (which=="k_dist") {
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<- c(1,2)
    kstr.save<-""
    distance <-c("TWS","RMSE")
    diststr.save<-dist
    legs<-c("500hPa_TWS","1000hPa_TWS","500hPa_RMSE","1000hPa_RMSE")
    style <- c(rep("blue",2),rep("red",2))
    type  <- c(1,2,1,2)
    dirstr.save<-get.dirstr(rean = rean)
  }
  
  colnam<-c("all ecdf","pos ecdf","p0 binom")
  
  # pour les lois de probabilite ajustees
  for (nam in colnam){
    print(nam)
    crps.mat<-NULL
    coln<-NULL
    # pour les distances selectionnees
    for(dist in distance){
      # pour les k selectionnes
      for (k in kvec){ 
        # pour les rayons selectionnes
        for (rad in rads){ 

          descriptors<-lapply(descr,paste.descr,"05")
          
          for (i in 1:length(descriptors)){ # Import et recuperation des scores selon les differents couples d'indicateurs pour une distribution
            if (substr(descriptors[[i]][1],1,1)=="A") load(file=paste0(get.dirstr(3,rean),"compute_crps-CV.A/02_",dist,"_member",member,"_k",3,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
            else load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,ifelse(threeday[i],"_threeday",""),".Rdata"))
            crps.mat<-cbind(crps.mat,crps[,nam])
            coln<-c(coln,paste0(namdescr[[i]],collapse="-"))
          }
        }
      }
    }
    
    # indice des lignes completes sans NA
    idx<-which(apply(crps.mat,1,function(v) all(!is.na(v)))) 
    print(length(idx))
    
    # CRPS moyen par indicateur
    meancrps<-apply(crps.mat[idx,],2,mean)
    
    if(substr(nam,1,3)=="pos"){
      # 3% fortes pluies ete
      ete <- substr(getdates(start,as.Date(end)-nbdays+1),6,7) %in% c("05","06","07")
      prec_ete <- precip[ete]
      qua <- quantile(prec_ete[prec_ete>0],probs=c(0.97))
      idx.ete <- intersect(idx,which(ete & precip>qua))
      meancrps.ete<-apply(crps.mat[idx.ete,],2,mean)
      
      # quantile 0.94 pluies positives
      qua <- quantile(precip[precip>0],probs=c(0.94))
      idx.0 <- intersect(idx,which(precip>qua))
      meancrps.0<-apply(crps.mat[idx.0,],2,mean)
      
      # 62*12 pluies fortes
      idx.1<-intersect(idx,soso$ix[1:(62*12)])
      meancrps.1<-apply(crps.mat[idx.1,],2,mean)
      
      # 62 pluies fortes
      idx.2<-intersect(idx,soso$ix[1:62])
      meancrps.2<-apply(crps.mat[idx.2,],2,mean)
      
      # par tranche de 10% de pluie positive
      sosobis <- rev(soso$ix[soso$x>0]) # indice des pluies positives croissantes
      pos <- round(quantile(1:length(sosobis),probs=seq(0,1,0.1),names=FALSE),0)
      idx.quant <- list()
      meancrps.quant <- list()
      for(i in 1:(length(pos)-1)){
        idx.quant[[i]] <- intersect(idx,sosobis[pos[i]:pos[i+1]])
        meancrps.quant[[i]]<-apply(crps.mat[idx.quant[[i]],],2,mean)
      }
      
    }
    
    # Score climato moyen
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
    
    normalize<-mean(score$crps[idx,nam])
    if(substr(nam,1,3)=="pos"){
      normalize.ete<-mean(score$crps[idx.ete,nam])
      normalize.0<-mean(score$crps[idx.0,nam])
      normalize.1<-mean(score$crps[idx.1,nam])
      normalize.2<-mean(score$crps[idx.2,nam])
      normalize.quant <- NULL
      for(i in 1:(length(pos)-1)) normalize.quant[i] <- mean(score$crps[idx.quant[[i]],nam])
    }
    
    # Graphiques
    
    if (seasonal) seastr<-"seasonal"
    else seastr<-"overall"
    
    filename<-paste0(dirstr.save,"compare.crps/",which,match(nam,colnam),"_",diststr.save,"_member",member,kstr.save,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),radstr.save,".png")
    
    # Fenêtre graphique
    if (substr(nam,1,3)=="pos"){
      png(file=filename,width=24,height=8,units="in",res=72)
      par(mar=c(max(nchar(coln))*2.3/3,4,4,2),mfrow=c(1,4))
    }else{
      png(file=filename,width=7,height=8,units="in",res=72)
      par(mar=c(max(nchar(coln))*2/3,4,3,2))
    }
    
    # plot
    plot(c(1,ndesc),c(0,max(1-meancrps/normalize)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45)) # Definition des proprietes du graphique
    for (j in 1:length(legs)) { # pour chaque legs, on ajoute une courbe
      points(1:ndesc,1-meancrps[(j-1)*ndesc+(1:ndesc)]/normalize,col=style[j],lty=type[j],pch=19)
      lines(1:ndesc,1-meancrps[(j-1)*ndesc+(1:ndesc)]/normalize,col=style[j],lty=type[j])
    }
    
    # Mise en forme
    if(comp) abline(h = 1-meancrps[(j-1)*(ndesc+1)+(ndesc+1)]/normalize,lty=3,lwd=2,col="red")
    box()
    axis(1,labels=coln,at=1:length(meancrps),las=3,cex.axis=ifelse(substr(nam,1,3)=="pos",1.8,1.3))
    axis(2)
    grid()
    abline(v=1:ndesc,lty=2,col=gray(0.5))
    title(paste0(nam," ",seastr))
    legend("topleft",legend=legs,pch=19,bty ="n",lty = type,col=style,
           cex=ifelse(substr(nam,1,3)=="pos",1.8,1.3))
    
    # Pluies positives
    if (substr(nam,1,3)=="pos"){
      
      # 3% fortes pluies ete
      plot(c(1,ndesc),c(0,max(1-meancrps.ete/normalize.ete)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45))
      for (j in 1:length(legs)) {
        points(1:ndesc,1-meancrps.ete[(j-1)*ndesc+(1:ndesc)]/normalize.ete,col=ifelse(which=="k_dist",style[j],j),lty=ifelse(which=="k_dist",type[j],1),pch=19)
        lines(1:ndesc,1-meancrps.ete[(j-1)*ndesc+(1:ndesc)]/normalize.ete,col=ifelse(which=="k_dist",style[j],j),lty=ifelse(which=="k_dist",type[j],1))
      }
      if(comp) abline(h = 1-meancrps.ete[(j-1)*(ndesc+1)+(ndesc+1)]/normalize.ete,lty=3,lwd=2,col="red")
      box()
      axis(1,labels=coln,at=1:length(meancrps.ete),las=3, cex.axis=1.8)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," summer max (",length(idx.ete)," largest - 3%)"))
      legend("topleft",legend=legs,pch=19,bty ="n",lty = type,col=style,cex=1.8)
      
      # 62*12
      plot(c(1,ndesc),c(0,max(1-meancrps.1/normalize.1)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45))
      for (j in 1:length(legs)) {
        points(1:ndesc,1-meancrps.1[(j-1)*ndesc+(1:ndesc)]/normalize.1,col=ifelse(which=="k_dist",style[j],j),lty=ifelse(which=="k_dist",type[j],1),pch=19)
        lines(1:ndesc,1-meancrps.1[(j-1)*ndesc+(1:ndesc)]/normalize.1,col=ifelse(which=="k_dist",style[j],j),lty=ifelse(which=="k_dist",type[j],1))
      }
      if(comp) abline(h = 1-meancrps.1[(j-1)*(ndesc+1)+(ndesc+1)]/normalize.1,lty=3,lwd=2,col="red")
      box()
      axis(1,labels=coln,at=1:length(meancrps.1),las=3, cex.axis=1.8)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," monthly max (62*12 largest)"))
      legend("topleft",legend=legs,pch=19,bty ="n",lty = type,col=style,cex=1.8)
      
      # 62
      plot(c(1,ndesc),c(0,max(1-meancrps.2/normalize.2)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45))
      for (j in 1:length(legs)) {
        points(1:ndesc,1-meancrps.2[(j-1)*ndesc+(1:ndesc)]/normalize.2,col=ifelse(which=="k_dist",style[j],j),lty=ifelse(which=="k_dist",type[j],1),pch=19)
        lines(1:ndesc,1-meancrps.2[(j-1)*ndesc+(1:ndesc)]/normalize.2,col=ifelse(which=="k_dist",style[j],j),lty=ifelse(which=="k_dist",type[j],1))
      }
      if(comp) abline(h = 1-meancrps.2[(j-1)*(ndesc+1)+(ndesc+1)]/normalize.2,lty=3,lwd=2,col="red")
      box()
      axis(1,labels=coln,at=1:length(meancrps.2),las=3, cex.axis=1.8)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," yearly max (62 largest)"))
      legend("topleft",legend=legs,pch=19,bty ="n",lty = type,col=style,cex=1.8)
      
      # par tranche de 10% de pluie positive
      #plot(c(1,ndesc),c(0,max(1-meancrps.1/normalize.1)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(-0.4,0.8))
      #for (j in 1:(length(pos)-1)) {
      #  points(1:ndesc,1-meancrps.quant[[j]]/normalize.quant[j],col=getcol(j,range = c(1,length(pos)-1)),lty=1,pch=19)
      #  lines(1:ndesc,1-meancrps.quant[[j]]/normalize.quant[j],col=getcol(j,range = c(1,length(pos)-1)),lty=1)
      #}
      #box()
      #axis(1,labels=coln,at=1:length(meancrps.1),las=3, cex.axis=1.8)
      #axis(2)
      #grid()
      #abline(v=1:ndesc,lty=2,col=gray(0.5))
      #title(paste0(nam," quantile"))
      #nleg <- matrix(seq(0.1,1,0.1),2,5,byrow = T)
      #legend("topleft",legend=nleg,pch=19,lty = type,bty="n",col=getcol(nleg),cex=1.6,ncol = 5,x.intersp=0.2)
      
      # Au dessus du quantile 0.94 de pluie positive
      graphics.off()
      filename<-paste0(dirstr.save,"compare.crps/",which,"4_",diststr.save,"_member",member,kstr.save,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),radstr.save,".png")
      png(file=filename,width=7,height=8,units="in",res=72)
      par(mar=c(max(nchar(coln))*2/3,4,3,2))
      
      plot(c(1,ndesc),c(0,max(1-meancrps.0/normalize.0)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.5)) # Definition des proprietes du graphique
      for (j in 1:length(legs)) { # pour chaque legs, on ajoute une courbe
        points(1:ndesc,1-meancrps.0[(j-1)*ndesc+(1:ndesc)]/normalize.0,col=style[j],lty=type[j],pch=19)
        lines(1:ndesc,1-meancrps.0[(j-1)*ndesc+(1:ndesc)]/normalize.0,col=style[j],lty=type[j])
      }
      
      # Mise en forme
      if(comp) abline(h = 1-meancrps.0[(j-1)*(ndesc+1)+(ndesc+1)]/normalize.0,lty=3,lwd=2,col="red")
      box()
      axis(1,labels=coln,at=1:length(meancrps),las=3,cex.axis=1.3)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," 0.94 quantile"))
      legend("topleft",legend=legs,pch=19,bty ="n",lty = type,col=style,cex=1.3)
      
    }
    
    
    
    dev.off()
  }
  
}

# Calcule et trace les CRPSS pour analogie classique pour differents rayons
compare.crps.A<-function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE) # classement des pluies avec index
  
  if (seasonal) strs<-c("1","2","5","10") # si seasonal, on prend des pourcentages d'analogues plus eleves car population plus petite
  else strs<-c("01","02","05","1","2","5","10")
  
  ns<-length(strs)
  colnam<-c("all ecdf","pos ecdf","p0 binom")#,"pos gamma")
  
  for (nam in colnam){ # pour chaque modele
    crps.mat<-NULL
    for (str in strs){ # pour chaque rayon, import des scores CRPS et remplissage de la matrice des scores
      load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/",str,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
      crps.mat<-cbind(crps.mat,crps[,nam])
    }
    
    colnames(crps.mat)<-strs
    
    idx<-which(apply(crps.mat,1,function(v) all(!is.na(v)))) # indices sans aucun NA
    print(length(idx))
    #  idx<-intersect(idx,which(!is.na(fit$loglik.i)))
    #  print(length(idx))
    
    meancrps<-apply(crps.mat[idx,],2,mean) # moyenne pour toutes les dates sans NA
    
    idx.1<-intersect(idx,soso$ix[1:62*12]) # moyenne pour les 62*12 plus fortes pluies
    meancrps.1<-apply(crps.mat[idx.1,],2,mean)
    
    idx.2<-intersect(idx,soso$ix[1:100]) # moyenne pour les 62 plus fortes pluies
    meancrps.2<-apply(crps.mat[idx.2,],2,mean)
    
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata")) # Import du score climatologique
    
    normalize<-mean(score$crps[idx,nam]) # moyenne des scores pour toutes les pluies
    normalize.1<-mean(score$crps[idx.1,nam]) # les 62*12 plus fortes
    normalize.2<-mean(score$crps[idx.2,nam]) # les 62 plus fortes
    
    if (seasonal) seastr<-"seasonal"
    else seastr<-"overall"
    
    filename<-paste0(get.dirstr(k,rean),"compare.crps.A/",match(nam,colnam),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png")
    if (substr(nam,1,3)=="pos"){ # si "pos", trois fenetres graphiques
      png(file=filename,width=14,height=8,units="in",res=72)
      par(mar=c(9,4,4,2),mfrow=c(1,3))
    }
    else{ # sinon une seule
      png(file=filename,width=7,height=8,units="in",res=72)
      par(mar=c(5.5,4,4,2))
    }
    plot(c(1,ns),c(0,max(1-meancrps/normalize)),axes=FALSE,xlab="",ylab="",type="n") # plot des CRPSS toutes pluies
    points(1:ns,1-meancrps/normalize)
    lines(1:ns,1-meancrps/normalize)
    
    #legend(legend=1:3,col=1:3)
    axis(2)
    axis(1,labels=strs,at=1:length(meancrps),las=3)
    box()
    
    title(paste0(nam," ",seastr))
    abline(v=1:ns,lty=2,col=gray(0.5))
    grid()
    
    if (substr(nam,1,3)=="pos"){ # si "pos"
      plot(c(1,ns),c(0,max(1-meancrps.1/normalize.1)),axes=FALSE,xlab="",ylab="",type="n") # plot des CRPSS des 62*12 plus fortes pluies
      points(1:ns,1-meancrps.1/normalize.1)
      lines(1:ns,1-meancrps.1/normalize.1)
      axis(2)
      axis(1,labels=strs,at=1:length(meancrps.1),las=3)
      box()
      title(paste0(nam," monthly max (62*12 largest)"))
      abline(v=1:ns,lty=2,col=gray(0.5))
      grid()
      
      plot(c(1,ns),c(0,max(1-meancrps.2/normalize.2)),axes=FALSE,xlab="",ylab="",type="n") # plot des CRPSS des 62 plus fortes pluies
      points(1:ns,1-meancrps.2/normalize.2)
      lines(1:ns,1-meancrps.2/normalize.2)
      axis(2)
      axis(1,labels=strs,at=1:length(meancrps.2),las=3)
      box()
      title(paste0(nam," 100 largest precip"))
      abline(v=1:ns,lty=2,col=gray(0.5))
      grid()
    }
    
    dev.off()
  }
}

# Calcule et trace les CPRSS et comparaison avec meilleure analogie classique
compare.crps.ana <- function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,CV=TRUE,rean){
  
  rad <- "nrn05"
  
  descr<-list(
    c("celnei","singnei"),
    c("singnei","rsingnei")
  )
  
  ndesc <- length(descr)
  namdescr <- lapply(descr,nam2str)
  
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE) # Classement des pluies de la plus forte a la plus faible avec indices

  colnam<-c("all ecdf","pos ecdf","p0 binom")

  for (nam in colnam){
    
    print(nam)
    crps.mat<-NULL
    coln<-NULL
    descriptors<-lapply(descr,paste.descr,"05")
    
    for (i in 1:length(descriptors)){ # Import et recuperation des scores selon les differents couples d'indicateurs pour une distribution
      load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".Rdata"))
      crps.mat<-cbind(crps.mat,crps[,nam])
      coln<-c(coln,paste0(namdescr[[i]],collapse="-"))
    }
    
    idx<-which(apply(crps.mat,1,function(v) all(!is.na(v)))) 
    print(length(idx))
    
    # CRPS moyen
    meancrps<-apply(crps.mat[idx,],2,mean)
    
    if(nam =="pos ecdf"){
      # quantile 0.94 pluies positives
      qua <- quantile(precip[precip>0],probs=c(0.94))
      idx.0 <- intersect(idx,which(precip>qua))
      meancrps.0<-apply(crps.mat[idx.0,],2,mean)
      
      # 62*12 pluies fortes
      idx.1<-intersect(idx,soso$ix[1:(62*12)])
      meancrps.1<-apply(crps.mat[idx.1,],2,mean)
      
      # 62 pluies fortes
      idx.2<-intersect(idx,soso$ix[1:62])
      meancrps.2<-apply(crps.mat[idx.2,],2,mean)
    }
    
    # Analogie classique
    load(file=paste0(get.dirstr(3,rean),"compute_crps-CV.A/02_TWS_member",member,"_k",3,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    meancrps.A<-mean(crps[idx,nam])
    if(nam =="pos ecdf"){
    meancrps.A.0<-mean(crps[idx.0,nam])
    meancrps.A.1<-mean(crps[idx.1,nam])
    meancrps.A.2<-mean(crps[idx.2,nam])
    }
    
    # CRPS climato moyen
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
    
    normalize<-mean(score$crps[idx,nam])
    if(nam == "pos ecdf"){
      normalize.0<-mean(score$crps[idx.0,nam])
      normalize.1<-mean(score$crps[idx.1,nam])
      normalize.2<-mean(score$crps[idx.2,nam])
    }
    
    # Graphiques
    filename<-paste0(get.dirstr(k,rean),"compare.crps.ana/",nam,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
    pdf(file = filename, width = 6, height = 6)
    par(mar=c(11,5,3,3))
    
    plot(c(1,ndesc),c(0,max(c(1-meancrps/normalize,1-meancrps.A/normalize))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
    points(1:ndesc,1-meancrps/normalize,pch=19)
    lines(1:ndesc,1-meancrps/normalize)
    abline(h=1-meancrps.A/normalize,col="red",lty=2)
    box()
    axis(1,labels=coln,at=1:length(meancrps),las=3,cex.axis=1.2)
    axis(2)
    grid()
    abline(v=1:ndesc,lty=2,col=gray(0.5))
    if (nam=="all ecdf") title("Overall")
    if (nam=="pos ecdf") title("Non-zero rainfall")
    if (nam=="p0 binom") title("Rain/No rain")
    graphics.off()
    
    if (nam =="pos ecdf"){
      # 0.94
      filename<-paste0(get.dirstr(k,rean),"compare.crps.ana/",nam,"_0.94_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps.0/normalize.0,1-meancrps.A.0/normalize.0))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps.0/normalize.0,pch=19)
      lines(1:ndesc,1-meancrps.0/normalize.0)
      abline(h=1-meancrps.A.0/normalize.0,col="red",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps.0),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," > 0.94 quantile"))
      graphics.off()

      # 62*12
      filename<-paste0(get.dirstr(k,rean),"compare.crps.ana/",nam,"_monthly_max_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps.1/normalize.1,1-meancrps.A.1/normalize.1))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps.1/normalize.1,pch=19)
      lines(1:ndesc,1-meancrps.1/normalize.1)
      abline(h=1-meancrps.A.1/normalize.1,col="red",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps.1),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," 12*62 largest rainfalls"))
      graphics.off()
      
      # 62
      filename<-paste0(get.dirstr(k,rean),"compare.crps.ana/",nam,"_yearly_max_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps.2/normalize.2,1-meancrps.A.2/normalize.2))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps.2/normalize.2,pch=19)
      lines(1:ndesc,1-meancrps.2/normalize.2)
      abline(h=1-meancrps.A.2/normalize.2,col="red",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps.2),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," 62 largest rainfalls"))
      graphics.off()
    }
  }
}

# Calcule et trace les CRPSS pour la double analogie combinee
compare.crps.comb<-function(which = 1,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",p=1,rad="05",standardize=TRUE,CV=TRUE,rean){
  
  descr <- list(
    c("sing","sing"),
    c("cel","sing"),
    c("snei","sing"),
    c("rsing","rsing"),
    c("cel","rsing"),
    c("snei","rsing"),
    c("sing","rsing"),
    c("cel","sing","rsing"),
    c("snei","sing","rsing"))
  
  descr2  <- lapply(descr,paste.descr,"05")
  ndescr  <- length(descr)
  
  precip <- get.precip(nbdays,start,end)
  soso   <- sort(precip,index.return=TRUE,decreasing=TRUE) # Classement des pluies de la plus forte a la plus faible avec indices
  
  colnam<-c("all ecdf","pos ecdf","p0 binom")
  
  for(nam in colnam){
    
    print(nam)
    
    # Import des scores analogie simple, indicateurs simples, double analogie
    load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    crps.ana     <- crps[,nam]
    
    crps.ind     <- matrix(NA,length(precip),ndescr)
    crps.ana.ind <- matrix(NA,length(precip),ndescr)
    
    coln <- vector(length = ndescr)
    
    for(i in 1:ndescr){
      
      load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descr2[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_nrn",rad,".Rdata"))
      crps.ind[,i] <- crps[,nam]
      
      load(paste0(get.dirstr(k,rean),"compute_crps_comb",get.CVstr(CV),"/",paste0(descr2[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_p",p,"_rad",rad,".Rdata"))
      crps.ana.ind[,i] <- crps[,nam]
      
      coln[i] <- paste0(descr[[i]],collapse="-")
    }
    
    colnames(crps.ind) <- colnames(crps.ana.ind) <- coln
    
    # Lignes sans NA
    tmp <- cbind(crps.ind, crps.ana.ind, crps.ana)
    idx<-which(apply(tmp,1,function(v) all(!is.na(v)))) # indice des lignes completes sans NA
    print(length(idx))
    
    # Moyenne des crps sur toute les sequences
    meancrps.ana     <- mean(crps.ana[idx])
    meancrps.ind     <- apply(crps.ind[idx,],2,mean)
    meancrps.ana.ind <- apply(crps.ana.ind[idx,],2,mean)
    
    # Si pluies positives, moyenne aussi sur les plus grosses pluies
    if(nam == "pos ecdf"){
      idx.1<-intersect(idx,soso$ix[1:(62*12)])
      meancrps.ana.1     <- mean(crps.ana[idx.1])
      meancrps.ind.1     <- apply(crps.ind[idx.1,],2,mean)
      meancrps.ana.ind.1 <- apply(crps.ana.ind[idx.1,],2,mean)
      
      idx.2<-intersect(idx,soso$ix[1:62])
      meancrps.ana.2     <- mean(crps.ana[idx.2])
      meancrps.ind.2     <- apply(crps.ind[idx.2,],2,mean)
      meancrps.ana.ind.2 <- apply(crps.ana.ind[idx.2,],2,mean)
    }
    
    # Score climato
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata")) # Import du score climatologique
    
    normalize<-mean(score$crps[idx,nam])
    if(nam == "pos ecdf"){
      normalize.1<-mean(score$crps[idx.1,nam])
      normalize.2<-mean(score$crps[idx.2,nam])
    }
    
    # Graphique
    if (seasonal) seastr<-"seasonal"
    else seastr<-"overall"
    filename<-paste0(get.dirstr(k,rean),"compare.crps.comb/",match(nam,colnam),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_p",p,"_rad",rad,".png")
    
    if (substr(nam,1,3)=="pos"){ # si "pos", trois fenetres graphiques
      png(file=filename,width=18,height=8,units="in",res=72)
      par(mar=c(9,4,4,2),mfrow=c(1,3))
    }
    else{ # sinon, une seule fenetre graphique
      png(file=filename,width=7,height=8,units="in",res=72)
      par(mar=c(5.5,4,4,2))
    }
    
    plot(c(1,ndescr),c(0,0.5),axes=FALSE,xlab="",ylab="",type="n") # Definition des proprietes du graphique
    points(1:ndescr,1-meancrps.ind/normalize,col="black",pch=19,type = "o")
    points(1:ndescr,1-meancrps.ana.ind/normalize,col="red",pch=19,type = "o")
    abline(h = 1-meancrps.ana/normalize,lty = 2, lwd = 2)
    
    axis(2)
    axis(1,labels=coln,at=1:ndescr,las=3)
    box()
    
    title(paste0(nam," ",seastr))
    grid()
    abline(v=1:ndescr,lty=2,col=gray(0.5))
    legend("topleft",legend = c("Ana","Ana + Ind", "Ind"), col = c("black","red","black"),
           lty = c(2,1,1), pch = c(NA,19,19), bty = "n", y.intersp = 1)
    
    
    if (substr(nam,1,3)=="pos"){ # si "pos" dans nam, on ajoute deux autres graphiques au png avec les 62*12 puis 12 plus fortes pluies
      # graphique 62*12
      plot(c(1,ndescr),c(0,0.5),axes=FALSE,xlab="",ylab="",type="n") # Definition des proprietes du graphique
      points(1:ndescr,1-meancrps.ind.1/normalize.1,col="black",pch=19,type = "o")
      points(1:ndescr,1-meancrps.ana.ind.1/normalize.1,col="red",pch=19,type = "o")
      abline(h = 1-meancrps.ana.1/normalize.1,lty = 2, lwd = 2)
      
      axis(2)
      axis(1,labels=coln,at=1:ndescr,las=3)
      box()
      
      title(paste0(nam," monthly max (62*12 largest)"))
      abline(v=1:ndescr,lty=2,col=gray(0.5))
      grid()
      
      # graphique 62
      plot(c(1,ndescr),c(0,0.5),axes=FALSE,xlab="",ylab="",type="n") # Definition des proprietes du graphique
      points(1:ndescr,1-meancrps.ind.2/normalize.2,col="black",pch=19,type = "o")
      points(1:ndescr,1-meancrps.ana.ind.2/normalize.2,col="red",pch=19,type = "o")
      abline(h = 1-meancrps.ana.2/normalize.2,lty = 2, lwd = 2)
      
      axis(2)
      axis(1,labels=coln,at=1:ndescr,las=3)
      box()
      
      title(paste0(nam," yearly max (62 largest)"))
      abline(v=1:ndescr,lty=2,col=gray(0.5))
      grid()
    }
    
    dev.off()
  }
  
}

# Comparaison des indicateurs par saison
compare.crps.sais<-function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",standardize=TRUE,CV=TRUE,rean){
  
  # Indicateurs a comparer
  descr<-list(
    c("sing","rsing"),
    c("singnei","rsingnei"),
    c("A","A"))
  ndesc <- length(descr)
  namdescr <- lapply(descr,nam2str)
  
  # Pluies
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE)
  
  # Boucle sur les trois scores: all,p>0,p0
  colnam<-c("all ecdf","pos ecdf","p0 binom")
  
  for (nam in colnam){
    print(nam)
    crps.mat<-NULL
    coln<-NULL
    
    descriptors<-lapply(descr,paste.descr,"05")
    
    for (i in 1:length(descriptors)){
      if (substr(descriptors[[i]][1],1,1)=="A") load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/05_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
      else load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
      crps.mat<-cbind(crps.mat,crps[,nam])
      coln<-c(coln,paste0(namdescr[[i]],collapse="-"))
    }
    
    # Score climato
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
    crpss.mat <- apply(crps.mat,2,function(x) 1-x/score$crps[,nam])
    
    # indice des lignes completes sans NA
    idx<-which(apply(crpss.mat,1,function(v) all(!is.na(v)))) 
    print(length(idx))
    
    # saisonnalite
    sais <- substr(getdates(start,as.character(as.Date(end)-nbdays+1)),6,7)
    sais[sais %in% c("12","01","02")] <- "DJF";sais[sais %in% c("03","04","05")] <- "MAM"
    sais[sais %in% c("06","07","08")] <- "JJA";sais[sais %in% c("09","10","11")] <- "SON"
    crpss.mat <- cbind(as.data.frame(crpss.mat),sais)
    colnames(crpss.mat) <- c(coln,"Saison")
    
    crpss.all <- gather(data = crpss.mat[idx,],key = "Indicateurs",value = "CRPSS",1:ndesc)
    crpss.all$Indicateurs <- factor(crpss.all$Indicateurs, levels = coln)
    crpss.all$Saison <- factor(crpss.all$Saison, levels = c("DJF","MAM","JJA","SON"))
    
    if(substr(nam,1,3)=="pos"){
      crpss.mmax <- gather(data = crpss.mat[intersect(idx,soso$ix[1:(62*12)]),],key = "Indicateurs",value = "CRPSS",1:ndesc)
      crpss.mmax$Indicateurs <- factor(crpss.mmax$Indicateurs, levels = coln)
      crpss.mmax$Saison <- factor(crpss.mmax$Saison, levels = c("DJF","MAM","JJA","SON"))
      
      crpss.ymax <- gather(data = crpss.mat[intersect(idx,soso$ix[1:62]),],key = "Indicateurs",value = "CRPSS",1:ndesc)
      crpss.ymax$Indicateurs <- factor(crpss.ymax$Indicateurs, levels = coln)
      crpss.ymax$Saison <- factor(crpss.ymax$Saison, levels = c("DJF","MAM","JJA","SON"))
    }
    
    # boxplot
    filename<-paste0(get.dirstr(k,rean),"compare.crps.sais/",match(nam,colnam),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".pdf")
    
    a <- ggplot(crpss.all, aes(x=Saison, y=CRPSS, fill=Indicateurs)) + 
      theme_bw()+
      geom_boxplot(outlier.shape = NA)+
      ylim(0,1)+
      stat_summary(fun.y = mean, geom="point",colour="black", size=1,position=position_dodge(width=0.75))+
      ggtitle(paste0(nam," overall"))+
      theme(plot.title = element_text(hjust = 0.5),legend.position = ifelse(substr(nam,1,3)=="pos","none","right"))+
      scale_fill_manual(values=c(getcol(vec = 1:(ndesc-1)),"red"))
      
    if(substr(nam,1,3)=="pos"){
      # 62*12
      b <- ggplot(crpss.mmax, aes(x=Saison, y=CRPSS, fill=Indicateurs)) + 
        theme_bw()+
        geom_boxplot(outlier.shape = NA)+
        ylim(0,1)+
        stat_summary(fun.y = mean, geom="point",colour="black", size=1,position=position_dodge(width=0.75))+
        ggtitle(paste0(nam," monthly max (62*12 largest)"))+
        theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+
        scale_fill_manual(values=c(getcol(vec = 1:(ndesc-1)),"red"))
      #62
      c <- ggplot(crpss.ymax, aes(x=Saison, y=CRPSS, fill=Indicateurs)) + 
        theme_bw()+
        geom_boxplot(outlier.shape = NA)+
        ylim(0,1)+
        stat_summary(fun.y = mean, geom="point",colour="black", size=1,position=position_dodge(width=0.75))+
        ggtitle(paste0(nam," yearly max (62 largest)"))+
        theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+
        scale_fill_manual(values=c(getcol(vec = 1:(ndesc-1)),"red"))
      
      a <- ggarrange(a, b, c, ncol = 3, nrow = 1)
    }
    
    ggsave(filename = filename,plot = a,width = ifelse(substr(nam,1,3)=="pos",12,8),height = 5)
    graphics.off()
    
  }
}

# Calcule et trace les CRPSS pour la double analogie: classique puis indicateurs 
compare.crps.TL<-function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radAna ="10",radInd="05",standardize=TRUE,CV=TRUE,rean){
  
  descr <- list(
    c("sing","sing"),
    c("cel","sing"),
    c("snei","sing"),
    c("rsing","rsing"),
    c("cel","rsing"),
    c("snei","rsing"),
    c("sing","rsing"),
    c("cel","sing","rsing"),
    c("snei","sing","rsing"))
  
  descr2  <- lapply(descr,paste.descr,"05")
  ndescr  <- length(descr)
  
  precip <- get.precip(nbdays,start,end)
  soso   <- sort(precip,index.return=TRUE,decreasing=TRUE) # Classement des pluies de la plus forte a la plus faible avec indices
  
  colnam<-c("all ecdf","pos ecdf","p0 binom")
  
  for(nam in colnam){
    
    print(nam)
    
    # Import des scores analogie simple, indicateurs simples, double analogie
    load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/",radInd,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    crps.ana     <- crps[,nam]
    
    crps.ind     <- matrix(NA,length(precip),ndescr)
    crps.ana.ind <- matrix(NA,length(precip),ndescr)
    
    coln <- vector(length = ndescr)
    
    for(i in 1:ndescr){
      
      load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descr2[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_nrn",radInd,".Rdata"))
      crps.ind[,i] <- crps[,nam]
      
      load(paste0(get.dirstr(k,rean),"compute_crps_TL",get.CVstr(CV),"/",paste0(descr2[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radAna,"_",radInd,".Rdata"))
      crps.ana.ind[,i] <- crps[,nam]
      
      coln[i] <- paste0(descr[[i]],collapse="-")
    }

    coln[coln!="A-A"] <- nam2str(coln[coln!="A-A"])
    colnames(crps.ind) <- colnames(crps.ana.ind) <- coln
    
    # Lignes sans NA
    tmp <- cbind(crps.ind, crps.ana.ind, crps.ana)
    idx<-which(apply(tmp,1,function(v) all(!is.na(v)))) # indice des lignes completes sans NA
    print(length(idx))
   
    # Moyenne des crps sur toute les sequences
    meancrps.ana     <- mean(crps.ana[idx])
    meancrps.ind     <- apply(crps.ind[idx,],2,mean)
    meancrps.ana.ind <- apply(crps.ana.ind[idx,],2,mean)
    
    # Si pluies positives, moyenne aussi sur les plus grosses pluies
    if(nam == "pos ecdf"){
    idx.1<-intersect(idx,soso$ix[1:(62*12)])
    meancrps.ana.1     <- mean(crps.ana[idx.1])
    meancrps.ind.1     <- apply(crps.ind[idx.1,],2,mean)
    meancrps.ana.ind.1 <- apply(crps.ana.ind[idx.1,],2,mean)
    
    idx.2<-intersect(idx,soso$ix[1:62])
    meancrps.ana.2     <- mean(crps.ana[idx.2])
    meancrps.ind.2     <- apply(crps.ind[idx.2,],2,mean)
    meancrps.ana.ind.2 <- apply(crps.ana.ind[idx.2,],2,mean)
    }
    
    # Score climato
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata")) # Import du score climatologique
    
    normalize<-mean(score$crps[idx,nam])
    if(nam == "pos ecdf"){
      normalize.1<-mean(score$crps[idx.1,nam])
      normalize.2<-mean(score$crps[idx.2,nam])
    }
    
    # Graphique
    if (seasonal) seastr<-"seasonal"
    else seastr<-"overall"
    filename<-paste0(get.dirstr(k,rean),"compare.crps.TL/",match(nam,colnam),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radAna,"_",radInd,".png")
    
    if (substr(nam,1,3)=="pos"){ # si "pos", trois fenetres graphiques
      png(file=filename,width=18,height=8,units="in",res=72)
      par(mar=c(11,4,4,2),mfrow=c(1,3))
    }
    else{ # sinon, une seule fenetre graphique
      png(file=filename,width=7,height=8,units="in",res=72)
      par(mar=c(8.5,4,4,2))
    }
    
    plot(c(1,ndescr),c(0,0.5),axes=FALSE,xlab="",ylab="",type="n") # Definition des proprietes du graphique
    points(1:ndescr,1-meancrps.ind/normalize,col="black",pch=19,type = "o")
    points(1:ndescr,1-meancrps.ana.ind/normalize,col="red",pch=19,type = "o")
    abline(h = 1-meancrps.ana/normalize,lty = 2, lwd = 2)
    
    axis(2)
    axis(1,labels=coln,at=1:ndescr,las=3,cex.axis=1.3)
    box()
    
    title(paste0(nam," ",seastr))
    grid()
    abline(v=1:ndescr,lty=2,col=gray(0.5))
    legend("topleft",legend = c("Analog","Analog + Descriptors", "Descriptors"), col = c("black","red","black"),
           lty = c(2,1,1), pch = c(NA,19,19), bty = "n", y.intersp = 1,cex=1.3)
    
    
    if (substr(nam,1,3)=="pos"){ # si "pos" dans nam, on ajoute deux autres graphiques au png avec les 62*12 puis 12 plus fortes pluies
      # graphique 62*12
      plot(c(1,ndescr),c(0,0.5),axes=FALSE,xlab="",ylab="",type="n") # Definition des proprietes du graphique
      points(1:ndescr,1-meancrps.ind.1/normalize.1,col="black",pch=19,type = "o")
      points(1:ndescr,1-meancrps.ana.ind.1/normalize.1,col="red",pch=19,type = "o")
      abline(h = 1-meancrps.ana.1/normalize.1,lty = 2, lwd = 2)
      
      axis(2)
      axis(1,labels=coln,at=1:ndescr,las=3,cex.axis=1.3)
      box()
      
      title(paste0(nam," monthly max (62*12 largest)"))
      abline(v=1:ndescr,lty=2,col=gray(0.5))
      grid()
      
      # graphique 62
      plot(c(1,ndescr),c(0,0.5),axes=FALSE,xlab="",ylab="",type="n") # Definition des proprietes du graphique
      points(1:ndescr,1-meancrps.ind.2/normalize.2,col="black",pch=19,type = "o")
      points(1:ndescr,1-meancrps.ana.ind.2/normalize.2,col="red",pch=19,type = "o")
      abline(h = 1-meancrps.ana.2/normalize.2,lty = 2, lwd = 2)
      
      axis(2)
      axis(1,labels=coln,at=1:ndescr,las=3,cex.axis=1.3)
      box()
      
      title(paste0(nam," yearly max (62 largest)"))
      abline(v=1:ndescr,lty=2,col=gray(0.5))
      grid()
    }
    
    dev.off()
  }
  
}

# Calcule et trace les CPRSS et comparaison avec les Weather Patterns d'EDF
compare.crps.wp <- function(k,dist,nbdays=1,start="1950-01-01",end="2011-12-31",standardize=TRUE,CV=TRUE,rean,ana=F){
  
  rad <- "nrn05"
  
  descr<-list(
    c("celnei","singnei"),
    c("singnei","rsingnei")
  )
  
  #descr <-list(
  #  c("singnei","rsingnei"),
  #  c("celnei_2","rsingnei"),
  #  c("sing05_2nei","rsingnei")
  #)
  
  ndesc <- length(descr)
  namdescr <- lapply(descr,nam2str)
  
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE) # Classement des pluies de la plus forte a la plus faible avec indices
  
  colnam<-c("all ecdf","pos ecdf","p0 binom")
  
  for (nam in colnam){
    
    print(nam)
    crps.mat<-NULL
    coln<-NULL
    descriptors<-lapply(descr,paste.descr,"05")
    
    for (i in 1:length(descriptors)){ # Import et recuperation des scores selon les differents couples d'indicateurs pour une distribution
      load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".Rdata"))
      crps.mat<-cbind(crps.mat,crps[,nam])
      coln<-c(coln,paste0(namdescr[[i]],collapse="-"))
    }
    
    idx<-which(apply(crps.mat,1,function(v) all(!is.na(v)))) 
    print(length(idx))
    
    # CRPS moyen
    meancrps<-apply(crps.mat[idx,],2,mean)
    
    if(nam =="pos ecdf"){
      # quantile 0.94 pluies positives
      qua <- quantile(precip[precip>0],probs=c(0.94))
      idx.0 <- intersect(idx,which(precip>qua))
      meancrps.0<-apply(crps.mat[idx.0,],2,mean)
      
      # 62*12 pluies fortes
      idx.1<-intersect(idx,soso$ix[1:(62*12)])
      meancrps.1<-apply(crps.mat[idx.1,],2,mean)
      
      # 62 pluies fortes
      idx.2<-intersect(idx,soso$ix[1:62])
      meancrps.2<-apply(crps.mat[idx.2,],2,mean)
    }
    
    # WP EDF
    load(file="2_Travail/Rresults/compute.crps.wp/crps_wp.Rdata")
    meancrps.wp<-mean(crps[idx,nam])
    if(nam =="pos ecdf"){
      meancrps.wp.0<-mean(crps[idx.0,nam])
      meancrps.wp.1<-mean(crps[idx.1,nam])
      meancrps.wp.2<-mean(crps[idx.2,nam])
    }
    
    # CRPS climato moyen
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
    
    normalize<-mean(score$crps[idx,nam])
    if(nam == "pos ecdf"){
      normalize.0<-mean(score$crps[idx.0,nam])
      normalize.1<-mean(score$crps[idx.1,nam])
      normalize.2<-mean(score$crps[idx.2,nam])
    }
    
    # Graphiques
    filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
    pdf(file = filename, width = 6, height = 6)
    par(mar=c(11,5,3,3))
    
    plot(c(1,ndesc),c(0,max(c(1-meancrps/normalize,1-meancrps.wp/normalize))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
    points(1:ndesc,1-meancrps/normalize,pch=19)
    lines(1:ndesc,1-meancrps/normalize)
    abline(h=1-meancrps.wp/normalize,col="royalblue",lty=2)
    box()
    axis(1,labels=coln,at=1:length(meancrps),las=3,cex.axis=1.2)
    axis(2)
    grid()
    abline(v=1:ndesc,lty=2,col=gray(0.5))
    if (nam=="all ecdf") title("Overall")
    if (nam=="pos ecdf") title("Non-zero rainfall")
    if (nam=="p0 binom") title("Rain/No rain")
    graphics.off()
    
    if (nam =="pos ecdf"){
      # 0.94
      filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_0.94_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps.0/normalize.0,1-meancrps.wp.0/normalize.0))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps.0/normalize.0,pch=19)
      lines(1:ndesc,1-meancrps.0/normalize.0)
      abline(h=1-meancrps.wp.0/normalize.0,col="royalblue",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps.0),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," > 0.94 quantile"))
      graphics.off()
      
      # 62*12
      filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_monthly_max_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps.1/normalize.1,1-meancrps.wp.1/normalize.1))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps.1/normalize.1,pch=19)
      lines(1:ndesc,1-meancrps.1/normalize.1)
      abline(h=1-meancrps.wp.1/normalize.1,col="royalblue",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps.1),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," 12*62 largest rainfalls"))
      graphics.off()
      
      # 62
      filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_yearly_max_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps.2/normalize.2,1-meancrps.wp.2/normalize.2))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps.2/normalize.2,pch=19)
      lines(1:ndesc,1-meancrps.2/normalize.2)
      abline(h=1-meancrps.wp.2/normalize.2,col="royalblue",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps.2),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      title(paste0(nam," 62 largest rainfalls"))
      graphics.off()
    }
    
    # Analogie classique
    if(ana){
      load(file=paste0(get.dirstr(3,rean),"compute_crps-CV.A/02_TWS_member",member,"_k",3,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
      meancrps.A<-mean(crps[idx,nam])
      if(nam =="pos ecdf"){
        meancrps.A.0<-mean(crps[idx.0,nam])
        meancrps.A.1<-mean(crps[idx.1,nam])
        meancrps.A.2<-mean(crps[idx.2,nam])
      }
      
      # Graphiques
      filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_ana_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
      pdf(file = filename, width = 6, height = 6)
      par(mar=c(11,5,3,3))
      
      plot(c(1,ndesc),c(0,max(c(1-meancrps/normalize,1-meancrps.wp/normalize,1-meancrps.A/normalize))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
      points(1:ndesc,1-meancrps/normalize,pch=19)
      lines(1:ndesc,1-meancrps/normalize)
      abline(h=1-meancrps.wp/normalize,col="royalblue",lty=2)
      abline(h=1-meancrps.A/normalize,col="red",lty=2)
      box()
      axis(1,labels=coln,at=1:length(meancrps),las=3,cex.axis=1.2)
      axis(2)
      grid()
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      if (nam=="all ecdf") title("Overall")
      if (nam=="pos ecdf") title("Non-zero rainfall")
      if (nam=="p0 binom") title("Rain/No rain")
      graphics.off()
      
      if (nam =="pos ecdf"){
        # 0.94
        filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_ana_0.94_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
        pdf(file = filename, width = 6, height = 6)
        par(mar=c(11,5,3,3))
        
        plot(c(1,ndesc),c(0,max(c(1-meancrps.0/normalize.0,1-meancrps.wp.0/normalize.0,1-meancrps.A.0/normalize.0))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
        points(1:ndesc,1-meancrps.0/normalize.0,pch=19)
        lines(1:ndesc,1-meancrps.0/normalize.0)
        abline(h=1-meancrps.wp.0/normalize.0,col="royalblue",lty=2)
        abline(h=1-meancrps.A.0/normalize.0,col="red",lty=2)
        box()
        axis(1,labels=coln,at=1:length(meancrps.0),las=3,cex.axis=1.2)
        axis(2)
        grid()
        abline(v=1:ndesc,lty=2,col=gray(0.5))
        title(paste0(nam," > 0.94 quantile"))
        graphics.off()
        
        # 62*12
        filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_ana_monthly_max_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
        pdf(file = filename, width = 6, height = 6)
        par(mar=c(11,5,3,3))
        
        plot(c(1,ndesc),c(0,max(c(1-meancrps.1/normalize.1,1-meancrps.wp.1/normalize.1,1-meancrps.A.1/normalize.1))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
        points(1:ndesc,1-meancrps.1/normalize.1,pch=19)
        lines(1:ndesc,1-meancrps.1/normalize.1)
        abline(h=1-meancrps.wp.1/normalize.1,col="royalblue",lty=2)
        abline(h=1-meancrps.A.1/normalize.1,col="red",lty=2)
        box()
        axis(1,labels=coln,at=1:length(meancrps.1),las=3,cex.axis=1.2)
        axis(2)
        grid()
        abline(v=1:ndesc,lty=2,col=gray(0.5))
        title(paste0(nam," 12*62 largest rainfalls"))
        graphics.off()
        
        # 62
        filename<-paste0(get.dirstr(k,rean),"compare.crps.wp/",nam,"_ana_yearly_max_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".pdf")
        pdf(file = filename, width = 6, height = 6)
        par(mar=c(11,5,3,3))
        
        plot(c(1,ndesc),c(0,max(c(1-meancrps.2/normalize.2,1-meancrps.wp.2/normalize.2,1-meancrps.A.2/normalize.2))*1.2),axes=FALSE,xlab="",ylab="CRPSS",type="n")
        points(1:ndesc,1-meancrps.2/normalize.2,pch=19)
        lines(1:ndesc,1-meancrps.2/normalize.2)
        abline(h=1-meancrps.wp.2/normalize.2,col="royalblue",lty=2)
        abline(h=1-meancrps.A.2/normalize.2,col="red",lty=2)
        box()
        axis(1,labels=coln,at=1:length(meancrps.2),las=3,cex.axis=1.2)
        axis(2)
        grid()
        abline(v=1:ndesc,lty=2,col=gray(0.5))
        title(paste0(nam," 62 largest rainfalls"))
        graphics.off()
      }
    }
  }
}

# Calcul de la relative sinfularite pour tous les rayons possibles
compute.rsing<- function(vec,sing=F){
  som <- cumsum(vec)
  pos <- 1:length(vec)
  if(sing) {res <- som/pos
  }else{res <- som/pos/vec}
  res
}

# Calcul des indicateurs
compute_criteria<-function(k,dist,start="1950-01-01",end="2011-12-31",update=FALSE,rean,threeday=FALSE){
  gc()
  
  print(paste0(get.dirstr(k,rean),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
  
  if (update) {
    load(file=paste0(file=paste0(get.dirstr(k,rean),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
    coln<-colnames(criteria)
  }
  
  dist.vec<-getdist(k,dist,start,end,rean,threeday)
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  dates<-getdates(start,ifelse(threeday,as.character(as.Date(end)-3+1),end))
  N<-length(dates)
  gc()
  
  U<-c(0,(N-1):1); # U = 0, 22644, 22643, 22642, ...
  sU<-sapply(1:(N-1),function(x) sum(U[1:x])) # on fait la somme de U[1], U[1:2], etc pour obtenir la position de la derniere distance qui separe chaque date
  gc()
  
  if (!update) {
    #coln.new<-c("cel","mind","sing05","sing1","sing2","sing5","lsing05","lsing1","lsing2","lsing5","pers05","pers1","pers2","pers5","q05","q1","q2","q5","pcel","pnei05","pnei1","pnei2","pnei5","snei05","snei1","snei2","snei5")
    coln.new<-c("q05","cel","snei05","sing05")
  }
  if (update) {
    coln.new<-c("celnei","persnei","singnei","rsingnei","accnei")
    #coln.new <- c("rsingnei")
  }
  
  criteria.new<-matrix(NA,ncol=length(coln.new),nrow=N)
  colnames(criteria.new)<-coln.new
  
  for (i in 1:N){
    ddi<-getdist4i(i,dist.vec,N,sU)
    #if (seasonal){
    #  j<-selind_season(i,len=30,dates)
    #  di<-ddi[j]
    #}
    #else 
      
    di<-ddi
    n<-length(di)
    
    gc()
    
    soso<-sort(di,index.return=TRUE) # classement par plus petit score, et donne les positions
    qi05<-di[soso$ix[(0.005*n)]]     # quantile 0.5%
    #qi1<-di[soso$ix[(0.01*n)]]       # 1%
    #qi2<-di[soso$ix[(0.02*n)]]       # 2%
    #qi5<-di[soso$ix[(0.05*n)]]       # 5%
    #qi10<-di[soso$ix[(0.1*n)]]       # 10%
    #Donne le rayon du cercle (les 0.5% les plus proches, les 1% les plus proches...)
    
    #gc()
    
    #if(i!=1) idi05v <- idi05
    idi05<-soso$ix[2:(0.005*n)] # recupere la position des 0.5% les plus proches
    #idi1<-soso$ix[2:(0.01*n)]   # des 1% les plus proches
    #idi2<-soso$ix[2:(0.02*n)]   # des 2% les plus proches
    #idi5<-soso$ix[2:(0.05*n)]   # des 5% les plus proches
    #idi10<-soso$ix[2:(0.1*n)]   # des 10% les plus proches
    
    #gc()
    #
    #ri05<-rle(as.integer(di<=qi05)) # calcule le temps que reste le score a l'interieur ou l'exterieur du quantile 0.5%
    #ri1<-rle(as.integer(di<=qi1))   # du quantile 1%
    #ri2<-rle(as.integer(di<=qi2))   # du quantile 2%
    #ri5<-rle(as.integer(di<=qi5))   # du quantile 5%
    #ri10<-rle(as.integer(di<=qi10)) # du quantile 10%
    
    tmp<-NULL
    
    for (cc in coln.new){
      
      # Celerite
      #if (cc=="cel") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,ddi[i-1])} # score de la journee avec la veille
      #if (cc=="cel+1") {if (i==N) tmp<-c(tmp,NA) else tmp<-c(tmp,ddi[i+1])} # score de la journee avec le lendemain
      
      #if(cc=="celAv"){if(1 %in% idi05) idi05 <- idi05[idi05!=1]; tmp <- c(tmp,mean(di[idi05-1]))} # score entre jour J et veille des voisins
      #if(cc=="celAvNorm"){if(1 %in% idi05) idi05 <- idi05[idi05!=1]; tmp <- c(tmp,mean(di[idi05-1])/mean(di[idi05]))}
      #if(cc=="celAvR"){if(1 %in% idi05) idi05 <- idi05[idi05!=1]; tmp <- c(tmp,mean(di[idi05-1])/qi05)}
      #
      #if(cc=="celAv1"){if(N %in% idi05) idi05 <- idi05[idi05!=N]; tmp <- c(tmp,mean(di[idi05+1]))} # score entre jour J et le lendemain des voisins
      #if(cc=="celAv1Norm"){if(N %in% idi05) idi05 <- idi05[idi05!=N]; tmp <- c(tmp,mean(di[idi05+1]/mean(di[idi05])))}
      #if(cc=="celAv1R"){if(N %in% idi05) idi05 <- idi05[idi05!=N]; tmp <- c(tmp,mean(di[idi05+1]/qi05))}
      #
      #if (cc=="celAp") {if (i==1) tmp<-c(tmp,NA) else {if(N %in% idi05v) idi05v <- idi05v[idi05v!=N]; tmp<-c(tmp,mean(di[idi05v+1]))}} # score entre jour J et lendemains des voisins de J-1
      #if (cc=="celApNorm") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean(di[idi05v+1])/mean(di[idi05]))}
      #if (cc=="celApR") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean(di[idi05v+1])/qi05)}
      #
      #if (cc=="celV") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean(di[idi05v]))} # score entre jour J et voisins de J-1
      #if (cc=="celVNorm") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean(di[idi05v])/mean(di[idi05]))}
      #if (cc=="celVR") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean(di[idi05v])/qi05)}
      
      # Acceleration
      #if (cc=="acc") {if (i==1) tmp<-c(tmp,NA) else tmp <- c(tmp,criteria[i,"cel"]/criteria[i-1,"cel"])}
      #if (cc=="acc10") {if (i %in% 1:10) tmp<-c(tmp,NA) else tmp <- c(tmp,criteria[i,"cel"]/mean(criteria[(i-10):(i-1),"cel"],na.rm=TRUE))}
      #if (cc=="accR") {if (i==1) tmp<-c(tmp,NA) else tmp <- c(tmp,criteria[i,"cel"]/criteria[i-1,"cel"]/qi05)}
      #if (cc=="acc10R") {if (i %in% 1:10) tmp<-c(tmp,NA) else tmp <- c(tmp,criteria[i,"cel"]/mean(criteria[(i-10):(i-1),"cel"],na.rm=TRUE)/qi05)}
      
      if (cc=="accnei") tmp <- c(tmp,mean(criteria[idi05,"acc"],na.rm=TRUE))
      #if (cc=="accneiR") tmp <- c(tmp,mean(criteria[idi05,"acc"],na.rm=TRUE)/qi05)
      # Minimum distance
      #if (cc=="mind") tmp<-c(tmp,di[soso$ix[2]]) # score minimum obtenu avec la meilleure journee analogue
      #
      ## Singularite
      #if (cc=="sing05") tmp<-c(tmp,mean(di[idi05])) # moyenne des distances (scores) des 0.5% les plus proches
      #if (cc=="sing1") tmp<-c(tmp,mean(di[idi1]))   # des 1% les plus proches
      #if (cc=="sing2") tmp<-c(tmp,mean(di[idi2]))   # des 2% les plus proches
      #if (cc=="sing5") tmp<-c(tmp,mean(di[idi5]))   # des 5% les plus proches
      #if (cc=="sing10") tmp<-c(tmp,mean(di[idi10])) # des 10% les plus proches
      if (cc=="singnei") tmp <- c(tmp,mean(criteria[idi05,"sing05"],na.rm=TRUE))
      if (cc=="rsingnei") tmp <- c(tmp,mean(criteria[idi05,"rsing05"],na.rm=TRUE))
      #
      ## Log Singularite
      #if (cc=="lsing05") tmp<-c(tmp,mean(log(di[idi05]))) # moyenne du logarithme des distances des 0.5% les plus proches
      #if (cc=="lsing1") tmp<-c(tmp,mean(log(di[idi1])))   # des 1% les plus proches
      #if (cc=="lsing2") tmp<-c(tmp,mean(log(di[idi2])))   # des 2% les plus proches
      #if (cc=="lsing5") tmp<-c(tmp,mean(log(di[idi5])))   # des 5% les plus proches
      #if (cc=="lsing10") tmp<-c(tmp,mean(log(di[idi10]))) # des 10% les plus proches
      #
      ## Persistence
      #if (cc=="pers05") tmp<-c(tmp,mean(ri05$lengths[ri05$values==1])) # moyenne du temps passe a l'interieur du quantile 0.5%
      #if (cc=="pers1") tmp<-c(tmp,mean(ri1$lengths[ri1$values==1]))    # du quantile 1%
      #if (cc=="pers2") tmp<-c(tmp,mean(ri2$lengths[ri2$values==1]))    # du quantile 2%
      #if (cc=="pers5") tmp<-c(tmp,mean(ri5$lengths[ri5$values==1]))    # du quantile 5%
      #if (cc=="pers10") tmp<-c(tmp,mean(ri10$lengths[ri10$values==1])) # du quantile 10%
      #
      ## Quantiles
      #if (cc=="q05") tmp<-c(tmp,qi05) # Quantile 0.5%
      #if (cc=="q1") tmp<-c(tmp,qi1)   # Quantile 1%
      #if (cc=="q2") tmp<-c(tmp,qi2)   # Quantile 2%
      #if (cc=="q5") tmp<-c(tmp,qi5)   # Quantile 5%
      #if (cc=="q10") tmp<-c(tmp,qi10) # Quantile 10%
      #
      ## Probabilite celerite: probabilite d'avoir un score en dessous de celui de la veille
      #if (cc=="pcel"){if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,ecdf(dim1)(ddi[i-1]))} # calcul de la fonction de repartition empirique (F=i/n, fonction creneau) et probabilite au non depassement de la veille (donne une indication sur son rang dans les analogues)
      #
      ## Persistence neighbour: probabilite qu'un jour du voisinage soit dans le voisinage de la veille le jour precedent (qu'il suive la meme trajectoire)
      #if (cc=="pnei05") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi05-1) %in% idi05m1))} # moyenne de la condition: les veilles des journees analogues sont dans les 0.5% les plus proches de la veille de la journee cible
      #if (cc=="pnei1") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi1-1) %in% idi1m1))}    # dans les 1% les plus proches
      #if (cc=="pnei2") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi2-1) %in% idi2m1))}    # dans les 2% les plus proches
      #if (cc=="pnei5") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi5-1) %in% idi5m1))}    # dans les 5% les plus proches
      #if (cc=="pnei10") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi10-1) %in% idi10m1))} # dans les 10% les plus proches
      #
      ## Persistence neighbour: probabilite qu'un jour du voisinage soit deja dans le voisinage le jour precedent
      #if (cc=="snei05") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi05-1) %in% idi05))} # moyenne de la condition: les veilles des journees analogues sont deja dans les 0.5% les plus proches de la journee cible
      #if (cc=="snei1") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi1-1) %in% idi1))}    # dans les 1% les plus proches
      #if (cc=="snei2") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi2-1) %in% idi2))}    # dans les 2% les plus proches
      #if (cc=="snei5") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi5-1) %in% idi5))}    # dans les 5% les plus proches
      #if (cc=="snei10") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi10-1) %in% idi10))} # dans les 10% les plus proches
      if (cc=="persnei") tmp <- c(tmp,mean(criteria[idi05,"snei05"],na.rm=TRUE))
      
      ## Celerite neighbour: on moyenne la celerite des plus proches voisins
      if (cc=="celnei") tmp<-c(tmp,mean(criteria[idi05,"cel"],na.rm=TRUE)) # moyenne des celerites des 0.5% les plus proches
      #if (cc=="celnei1") tmp<-c(tmp,mean(criteria[idi1,"cel"],na.rm=TRUE))   # des 1% les plus proches
      #if (cc=="celnei2") tmp<-c(tmp,mean(criteria[idi2,"cel"],na.rm=TRUE))   # des 2% les plus proches
      #if (cc=="celnei5") tmp<-c(tmp,mean(criteria[idi5,"cel"],na.rm=TRUE))   # des 5% les plus proches
      #if (cc=="celnei10") tmp<-c(tmp,mean(criteria[idi10,"cel"],na.rm=TRUE)) # des 10% les plus proches
      
      ## TWSgeo: score TWS entre les geopotentiels 500 et 1000 normalises par leur sd respectifs
      #if (cc=="TWSgeonei") tmp<-c(tmp,mean(criteria[idi05,"TWSgeo"],na.rm=TRUE))
      
      # La singularite relative rsing (ou local dimension) est calculee dans la fonction get.descriptor
      
    }
    
    #dim1<-di
    #idi05m1<-idi05
    #idi1m1<-idi1
    #idi2m1<-idi2
    #idi5m1<-idi5
    #idi10m1<-idi10
    
    criteria.new[i,]<-tmp
    
    if (i %% 50==0) {
      print(i)
      #print(criteria.new[(i-49):i,])
    }
    
    #dim1<-di
    gc()
    
  }
  if (update) {
    criteria<-cbind(criteria,criteria.new)
    colnames(criteria)<-c(coln,coln.new)
  }
  else {
    criteria<-criteria.new
    colnames(criteria)<-coln.new
  }
  
  save(criteria,file=paste0(file=paste0(get.dirstr(k,rean),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
}

# Calcul des CRPS apres analogie au sens des indicateurs
compute_crps<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean,best=FALSE,threeday=FALSE){
  
  print(paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,ifelse(threeday,"_threeday",""),".Rdata"))
  
  precip<-get.precip(nbdays,start,end)
  N<-length(precip)
  
  descr <- matrix(NA,N,length(descriptors))
  for (i in 1:length(descriptors)) descr[,i]<-get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize,rean,threeday) # on importe les descripteurs
  
  #filegamma<-paste0(get.dirstr(k,rean),"fit.loglik.gamma",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata")
  #fileegp<-paste0(get.dirstr(k,rean),"fit.loglik.egp",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata")
  #print(filegamma)
  #if (file.exists(filegamma))  { # on importe les loi gamma
  #  print("exists!")
  #  load(file=filegamma)
  #  fit.gamma<-fit
  #}
  #print(fileegp)
  #if (file.exists(fileegp))  { # on importe les lois de Pareto
  #  print("exists!")
  #  load(file=fileegp)
  #  fit.egp<-fit
  #}
  if(best){
    load(file="2_Travail/20CR/Rresults/overall/k2/fit.empir-CV/persnei_singnei_RMSE_member1_k2_mean3day_1950-01-01_2011-12-31_std_nrn05.Rdata")
    paramp0 <- param[,"p0"]
  }
  
  load(file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,ifelse(threeday,"_threeday",""),".Rdata")) # on importe la loi empirique
  
  crps<-matrix(NA,ncol=5,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom","pos gamma","pos egp") # memes colonnes que pour score climato
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){ # si les indicateurs du jour ne sont pas nuls
      tmp<-get.closest(i,descr,precip,CV,nbdays,radtype) # on prend les voisins selon un certain rayon
      pi<-precip[tmp$idx]
      
      if(best){
        nbsec <- paramp0[i]*length(pi) # nombre de 0 a avoir
        pi <- pi[pi != 0] # on enleve tous les 0
        pi <- c(pi,rep(0,nbsec)) # on ajoute le nombre souhaite
        crps[i,"all ecdf"]<-crps_sample(precip[i],pi) # crps all ecdf
        crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,paramp0[i]) # crps p0 binom
      } else{
        crps[i,"all ecdf"]<-crps_sample(precip[i],pi) # crps all ecdf
        crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,param[i,"p0"]) # crps p0 binom
      }
      if (precip[i]>0) {
        crps[i,"pos ecdf"]<-crps_sample(precip[i],pi[pi>0]) # crps pos ecdf
        #if (file.exists(filegamma)) crps[i,"pos gamma"]<-crps_gamma(precip[i],shape=fit.gamma$param[i,"shape"],scale=fit.gamma$param[i,"scale"]) # crps pos gamma
        #if (file.exists(fileegp)) {
        #  simu<-regp2(1000,kappa=fit.egp$param[i,"kappa"],sigma=fit.egp$param[i,"sigma"],xi=fit.egp$param[i,"xi"])
        #  crps[i,"pos egp"]<-crps_sample(precip[i],simu) # crps pos egp
        #  gc()
        #}
      }
    }
  }
  save(crps,file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),ifelse(best,"_best",""),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,ifelse(threeday,"_threeday",""),".Rdata"))
}

# Calcul des CRPS apres analogie au sens de l'analogie classique et des indicateurs combines
compute_crps_comb<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,p=1,rad="05",CV=TRUE,rean){
  
  precip<-get.precip(nbdays,start,end)
  N<-length(precip)
  
  load(file=paste0(get.dirstr(k,rean),"fit.empir.comb",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_p",p,"_rad",rad,".Rdata"))
  crps<-matrix(NA,ncol=5,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom","pos gamma","pos egp") # memes colonnes que pour score climato
  
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (length(fit.emp$ind[[i]]) != 0){
      pi<-precip[fit.emp$ind[[i]]]
      crps[i,"all ecdf"]<-crps_sample(precip[i],pi) # crps all ecdf
      crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,fit.emp$param[i,"p0"]) # crps p0 binom
      if (precip[i]>0) {
        crps[i,"pos ecdf"]<-crps_sample(precip[i],pi[pi>0]) # crps pos ecdf
      }
    }
  }
  
  save(crps,file=paste0(get.dirstr(k,rean),"compute_crps_comb",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_p",p,"_rad",rad,".Rdata"))
}

# Calcul des CRPS apres analogie au sens de l'analogie classique ET des indicateurs
compute_crps_TL<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radAna ="10",radInd="05",CV=TRUE,rean){
  
  precip<-get.precip(nbdays,start,end)
  N<-length(precip)
  
  load(file=paste0(get.dirstr(k,rean),"fit.empir.TL",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radAna,"_",radInd,".Rdata")) # on importe la loi empirique
  
  crps<-matrix(NA,ncol=5,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom","pos gamma","pos egp") # memes colonnes que pour score climato
  
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (length(fit.emp$ind[[i]]) != 0){
      pi<-precip[fit.emp$ind[[i]]]
      crps[i,"all ecdf"]<-crps_sample(precip[i],pi) # crps all ecdf
      crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,fit.emp$param[i,"p0"]) # crps p0 binom
      if (precip[i]>0) {
        crps[i,"pos ecdf"]<-crps_sample(precip[i],pi[pi>0]) # crps pos ecdf
      }
    }
  }
  
  save(crps,file=paste0(get.dirstr(k,rean),"compute_crps_TL",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radAna,"_",radInd,".Rdata"))
}

# Calcul des CRPS apres analogie au sens de l'analogie classique
compute_crps.A<-function(rad,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  print(paste0(get.dirstr(k,rean),"compute_crps-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  precip<-get.precip(nbdays,start,end)
  
  N<-length(precip)
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[[rad]]
  
  #load(file=paste0(get.dirstr(k,rean),"fit.loglik.gamma-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  #fit.gamma<-fit
  
  #load(file=paste0(get.dirstr(k,rean),"fit.loglik.egp-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  #fit.egp<-fit
  
  load(file=paste0(get.dirstr(k,rean),"fit.loglik.p0-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  fit.p0<-fit
  
  crps<-matrix(NA,ncol=5,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom","pos gamma","pos egp")
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(precip[i])){
      pi<-precip[setdiff(nei[[i]],(i-nbdays+1):(i+nbdays-1))]
      pi<-pi[which(!is.na(pi))]
      crps[i,"all ecdf"]<-crps_sample(precip[i],pi)
      crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,fit.p0$p0[i])
      if (precip[i]>0 & sum(pi)!=0) {
        crps[i,"pos ecdf"]<-crps_sample(precip[i],pi[pi>0])
        #crps[i,"pos gamma"]<-crps_gamma(precip[i],shape=fit.gamma$param[i,"shape"],scale=fit.gamma$param[i,"scale"])
        #simu<-regp2(1000,kappa=fit.egp$param[i,"kappa"],sigma=fit.egp$param[i,"sigma"],xi=fit.egp$param[i,"xi"])
        #crps[i,"pos egp"]<-crps_sample(precip[i],simu)
        #gc()
      }
    }
    
  }
  
  save(crps,file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Calcul des CRPS avec les WP et saison a risque
compute.crps.wp <- function(start="1950-01-01",end="2011-12-31"){
  
  # Import
  precip <- get.precip(1,start,end)
  wp <- get.wp(start,end,risk=T)
  load("2_Travail/Rresults/save.wp.ind/WP_ind.Rdata")
  
  # Calcul
  N <- length(precip)
  crps<-matrix(NA,ncol=3,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom")
  
  for(i in 1:N){
    if (i %%50==0) print(i)
    pi <- precip[wp.ind[[wp[i]]]]
    p0 <- sum(pi==0)/length(pi)
    
    crps[i,"all ecdf"]<-crps_sample(precip[i],pi) # crps all ecdf
    crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,p0) # crps p0 binom
    if (precip[i]>0) crps[i,"pos ecdf"]<-crps_sample(precip[i],pi[pi>0]) # crps pos ecdf
  }
  
  # Sauvegarde
  save(crps,file="2_Travail/Rresults/compute.crps.wp/crps_wp.Rdata")
  
}

# Calcul des distances de maniere generique
compute_dist.gen<-function(k,dist,start="1950-01-01",end="2011-12-31",rean){
  if (k %in% 1:2) {
    if (dist %in% c("TWS","RMSE","RMSE_I","RMSE_II")){
      if (dist=="TWS") dist.list<-compute_TWS(k,start,end,rean)
      if (dist=="RMSE") dist.list<-compute_RMSE(k,start,end,rean)
      if (dist=="RMSE_I") dist.list<-compute_RMSE_I(k,start,end,rean)
      if (dist=="RMSE_II") dist.list<-compute_RMSE_II(k,start,end,rean)
      save(dist.list,file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
    }
    else{ # si on demande un nTWS, sTWS, nRMSE ou sRMSE, les distances sont normalisees par la moyenne ou l'ecart type des distances
      dist0<-substr(dist,2,nchar(dist))
      load(file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist0,"_member",member,"_Z",Z,"_",start,"_",end,".Rdata"))
      type<-substr(dist,1,1)
      print(type)
      if (type=="n") norm<-mean(unlist(dist.list))
      if (type=="s") norm<-sd(unlist(dist.list))
      print(norm)
      for (i in 1:length(dist.list)) dist.list[[i]]<-dist.list[[i]]/norm
      
      save(dist.list,file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
    }
  }
  else { # si k=3, les scores sont deja censes etre calcules pour 500 et 1000, donc on fait juste la moyenne des deux
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist,"_member",member,"_Z500_",start,"_",end,".Rdata"))
    print("1")
    dist1.list<-dist.list
    rm(dist.list);gc(TRUE)
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist,"_member",member,"_Z1000_",start,"_",end,".Rdata"))
    print("2")
    for (i in 1:length(dist.list)) dist.list[[i]]<-(dist.list[[i]]+dist1.list[[i]])/2
    rm(dist1.list);gc(TRUE)
    save(dist.list,file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
  }
}

# Calcul des CRPS climato: performance des differentes distributions sachant les autres jours (pas de selection par analogie)
compute_score_climato<-function(nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  precip<-get.precip(nbdays,start,end)
  N<-length(precip)
  
  nloglik.fct<-function(param,x){
    -sum(dgamma(x, shape=param[2], scale = param[1], log = TRUE) ) # calcul de la densite de probabilite de la loi gamma, et somme des loglikelihood
  }
  
  precip0<-precip[precip>0]
  mean0<-mean(precip0)
  var0<-var(precip0)
  scale<-var0/mean0 # premiere estimation du parametre d'echelle avec la methode des moments: var/mean
  shape<-mean0^2/var0 # premiere estimation du parametre de forme avec la methode des moments: mean^2/var
  print("1")
  opt<-optim(c(scale,shape),nloglik.fct,method="L-BFGS-B",x=precip0,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) # optimisation des parametres de forme et d'echelle pour maximiser la fonction du maximum de vraisemblance
  print("2")
  names(opt$par)<-c("scale","shape") # opt contient les nouveaux parametres de la loi gamma optimisee par max de vraisemblance
  print(opt)
  
  precip0<-precip[precip>10^-10]
  init=c(0.9, fpot(precip0,quantile(precip0,0.8),shape=0.05)$param) # Peak Over Threshold: maximum de vraisemblance utilisant la representation generale de Pareto (seuil = q80%)
  opt.egp<-egp2.fit(precip0, model=1, method="mle", init=init,rounded=0,CI=FALSE,plots=FALSE) # optimisation des parametres de la loi precedente
  
  #F2<-function(x) pegp2(x,kappa=opt.egp$fit$mle["kappa"],sigma=opt.egp$fit$mle["sigma"],xi=opt.egp$fit$mle["xi"])^2
  #F3<-function(x) (pegp2(x,kappa=opt.egp$fit$mle["kappa"],sigma=opt.egp$fit$mle["sigma"],xi=opt.egp$fit$mle["xi"])-1)^2
  
  israin<-(precip>0)
  p0<-1-sum(israin)/sum(!is.na(precip)) # probabilite empirique d'avoir un jour sec (ou un sequence seche)
  
  crps<-matrix(NA,ncol=5,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom","pos gamma","pos egp") # ecdf = Empirical Cumulative Distribution Function, egp = Extended Generelised Pareto
  crps[,"p0 binom"]<-crps_binom(israin,1,1-p0) # loi binomiale
  loglik.gamma<-rep(NA,N)
  loglik.egp<-rep(NA,N)
  for (i in 1:N){
    if (i %%50==0) print(i)
    crps[i,"all ecdf"]<-crps_sample(precip[i],precip) # crps entre distribution empirique du jour et distribution empirique de tous les jours
    if (precip[i]>0) {
      crps[i,"pos ecdf"]<-crps_sample(precip[i],precip0) # crps entre distribution empirique du jour et distribution empirique des jours pluvieux
      crps[i,"pos gamma"]<-crps_gamma(precip[i],shape=opt$par["shape"],scale=opt$par["scale"]) # crps entre distribution empirique du jour et distribution loi gamma des jours pluvieux
      loglik.gamma[i]<- dgamma(precip[i], shape=opt$par["shape"], scale = opt$par["scale"], log = TRUE)
      simu<-regp2(1000,kappa=opt.egp$fit$mle["kappa"],sigma=opt.egp$fit$mle["sigma"],xi=opt.egp$fit$mle["xi"])
      crps[i,"pos egp"]<-crps_sample(precip[i],simu) # crps entre distribution empirique du jour et distribution Pareto des jours pluvieux
      #crps[i,"pos egp"]<-integrate(F2,lower=0,upper=precip[i])$value+integrate(F3,lower=precip[i],upper=Inf)$value
      loglik.egp[i]<- degp2(precip[i], kappa=opt.egp$fit$mle["kappa"], sigma=opt.egp$fit$mle["sigma"], xi=opt.egp$fit$mle["xi"], type=1, log=TRUE)
      gc()
    }
  }
  
  score<-list(crps=crps,p0=p0,par.gamma=opt$par,opt.gamma=opt,loglik.p0=(1-israin)*log(p0)+israin*log(1-p0),loglik.gamma=loglik.gamma,par.egp=opt.egp$fit$mle,opt.egp=opt.egp,loglik.egp=loglik.egp)
  
  save(score,file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Calcul des scores RMSE
compute_RMSE<-function(k,start="1950-01-01",end="2011-12-31",rean){
  dat<-getdata(k,start,end,rean)
  N<-dim(dat)[3]
  dist.list<-list()
  #tt<-Sys.time()
  for (i in 1:(N-1)){#1:
    if (i%%10==0) print(i)
    j<-(i+1):N
    
    datj<-dat[,,j]
    if (is.matrix(datj)) datj<-array(datj,c(dim(datj),1)) # si derniere date, on force le format de la matrice en array
    dati<-array(dat[,,i],dim(datj)) # on repete le geopotentiel du jour le meme nombre de fois que le nombre de dates ulterieures
    
    dist.list[[i]]<-apply(dati-datj,3,function(v) sqrt(mean(v^2))) # difference d'altitude au carre entre chaque point des deux jours, puis moyenne de ces differences et racine carree
  }
  dist.list
}

# Calcul de RMSE_I: delta des pts normalises par la moyenne du jour
compute_RMSE_I<-function(k,start="1950-01-01",end="2011-12-31",rean){
  
  # Import geopotentiel et score RMSE
  load.nc(rean)
  dat<-getdata(k,start,end,rean)
  N<-dim(dat)[3]
  RMSE<-getdist(k,"RMSE",start,end,rean)
  
  # Calcul nouveau score
  moy_geo<-apply(dat,3,mean)
  dist.list<-RMSE
  
  for(i in 1:(N-1)){
    if (i%%50==0) print(i)
    dist.list[[i]] <- sqrt(dist.list[[i]]^2 - (moy_geo[i]-moy_geo[(i+1):length(moy_geo)])^2)
  }
  dist.list
}

# Calcul de RMSE_II: delta des moyennes du jour
compute_RMSE_II<-function(k,start="1950-01-01",end="2011-12-31",rean){
  
  # Import scores RMSE
  load(paste0("2_Travail/",rean,"/Rresults/compute_dist/RMSE_I_member1_k",k,"_",start,"_",end,".Rdata"))
  dist.list.rmseI <- dist.list
  
  load(paste0("2_Travail/",rean,"/Rresults/compute_dist/RMSE_member1_k",k,"_",start,"_",end,".Rdata"))
  dist.list.rmse <- dist.list
  rm(dist.list)
  
  # Calcul de RMSE II
  dist.list <- lapply(mapply("-",lapply(dist.list.rmse,function(v) v^2),lapply(dist.list.rmseI,function(v) v^2),SIMPLIFY = FALSE),sqrt)
  dist.list
}

# Calcul des scores TWS
compute_TWS<-function(k,start="1950-01-01",end="2011-12-31",rean){
  gradlist<-grad(k,start,end,rean)
  attach(gradlist)#gradlon,gradlat
  N<-dim(gradlon)[3]
  dist.list<-list()
  #tt<-Sys.time()
  for (i in 1:(N-1)){#1:
    if (i%%10==0) print(i)
    j<-(i+1):N
    gradlonj<-gradlon[,,j]
    gradlatj<-gradlat[,,j]
    if (is.matrix(gradlonj)) gradlonj<-array(gradlonj,c(dim(gradlonj),1)) # si on arrive a la derniere journee, on force la matrice en array pour derouler le code qui suit sans probleme
    if (is.matrix(gradlatj)) gradlatj<-array(gradlatj,c(dim(gradlatj),1))
    
    gradloni<-array(gradlon[,,i],dim(gradlonj)) # on repete j fois le champs de pression de la journee i
    gradlati<-array(gradlat[,,i],dim(gradlatj))
    
    #### numerateur du TWS
    fct<-function(v) sum(abs(v))
    dif_grad<-apply(gradloni-gradlonj,3,fct)+apply(gradlati-gradlatj,3,fct) # on calcule la somme des differences entre les gradients de chaque point de grille du jour i et le gradient de chaque point de grille de tous les jours j
    
    #### denominateur du TWS
    max_grad<-apply(pmax(abs(gradloni),abs(gradlonj)),3,sum)+apply(pmax(abs(gradlati),abs(gradlatj)),3,sum) # on calcule la somme des gradients maximums pour chaque point de grille entre ceux du jour i et ceux de tous les jours j 
    
    ####calcul du TWS
    dist.list[[i]]<-dif_grad/max_grad/2
    
    #if (i==10){
    #  print(dist.mat[1:10,1:10])
    #  print(Sys.time()-tt)
    #  stop()
    #}
    gc()
  }
  dist.list
}

# Calcul de l'indicateur D500_1000 pour les 4 distances
compute_D500_1000<-function(start="1950-01-01",end="2011-12-31",rean){
  
  # Import des champs de geopotentiels
  load.nc(rean)
  dat500<-getdata(1,start,end,rean)
  dat1000<-getdata(2,start,end,rean,large_win = T)
  
  # Calcul TWS
  print("TWS")
  gradlist500 <- grad(k=1,start,end,rean)
  gradlist1000 <- grad(k=2,start,end,rean,large_win = T)
  fct <- function(v) sum(abs(v))
  num <- apply(gradlist500$gradlon - gradlist1000$gradlon,3,fct) + apply(gradlist500$gradlat - gradlist1000$gradlat,3,fct)
  den <- apply(pmax(abs(gradlist500$gradlon),abs(gradlist1000$gradlon)),3,sum)+apply(pmax(abs(gradlist500$gradlat),abs(gradlist1000$gradlat)),3,sum)
  D500_1000_TWS <- num/den/2
  
  load("2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_TWS_member1_k1_1950-01-01_2011-12-31.Rdata")
  criteria <- cbind(criteria,D500_1000_TWS)
  colnames(criteria)[ncol(criteria)] <- "D_500_1000"
  save(criteria, file="2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_TWS_member1_k1_1950-01-01_2011-12-31.Rdata")
  
  # RMSE
  print("RMSE")
  D500_1000_RMSE <- apply(dat500-dat1000,3,function(v) sqrt(mean(v^2)))
  
  load("2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_RMSE_member1_k1_1950-01-01_2011-12-31.Rdata")
  criteria <- cbind(criteria,D500_1000_RMSE)
  colnames(criteria)[ncol(criteria)] <- "D_500_1000"
  save(criteria, file="2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_RMSE_member1_k1_1950-01-01_2011-12-31.Rdata")
  
  # RMSE_II
  print("RMSE_II")
  D500_1000_RMSE_II <- apply(dat500,3,mean) - apply(dat1000,3,mean)
  
  load("2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_RMSE_II_member1_k1_1950-01-01_2011-12-31.Rdata")
  criteria <- cbind(criteria,D500_1000_RMSE_II)
  colnames(criteria)[ncol(criteria)] <- "D_500_1000"
  save(criteria, file="2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_RMSE_II_member1_k1_1950-01-01_2011-12-31.Rdata")
  
  # RMSE_I
  print("RMSE_I")
  D500_1000_RMSE_I <- sqrt(D500_1000_RMSE^2 - D500_1000_RMSE_II^2)
  
  load("2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_RMSE_I_member1_k1_1950-01-01_2011-12-31.Rdata")
  criteria <- cbind(criteria,D500_1000_RMSE_I)
  colnames(criteria)[ncol(criteria)] <- "D_500_1000"
  save(criteria, file="2_Travail/20CR/Rresults/overall/k1/compute_criteria/criteria_RMSE_I_member1_k1_1950-01-01_2011-12-31.Rdata")
  
}

# Trouve le jour de NetCDF associe a la date
date_to_number<-function(nc,day,rean){
  which(getdays(nc,rean)==day)
}

# Calcul des parametres de la loi empirique selon le voisinage au sens des indicateurs
fit.empir<-function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,threeday=c(F,F)){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype)
  }                     
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype)
  }
  
  # Import des precip et des indicateurs
  precip<-get.precip(nbdays,start,end)
  
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize,rean[2],threeday[2])
  descr<- cbind(descr1,descr2)
  
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=FALSE,rean[1],threeday[1]) # pas de standardisation pour la densite de pts dans le plan
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=FALSE,rean[2],threeday[2])
  descrBis<- cbind(descr1,descr2)
  rad <- mean(c(sd(descrBis[,1],na.rm=TRUE),sd(descrBis[,2],na.rm=TRUE)))/2
  
  if (nrow(descr) != length(precip)) stop("Probleme taille des donnees")
  
  N<-length(precip)
  
  param<-matrix(NA,ncol=7,nrow=N)
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){
      tmp<-get.closest(i,descr,precip,CV,nbdays,radtype)
      pi<-precip[tmp$idx]
      ni<-length(pi)
      ni1<-sum(pi>0)
      count <- nn2(data = descrBis[tmp$idef,],query = t(descrBis[i,]),searchtype = "radius",radius =  rad,k = length(tmp$idef)) # nombre de voisins dans le rayon (tmp$idef: indices ou les deux descr sont non NA)
      nb <- rowSums(count$nn.idx>0)-1
      parami<-c(nb,1-ni1/ni,mean(pi[pi>0]),sd(pi[pi>0]),skewness(pi[pi>0]),kurtosis(pi[pi>0]),tmp$radius) # remplissage de la ligne: nbre de voisins, proba jour sec, moyenne et ecart-type des pluies non nulles, coefficient d'assymetrie et d'applatissement, rayon du cercle
    }
    else parami<-rep(NA,7) # si un des deux indicateurs est nul, ligne du tableau parametre vide
    param[i,]<-parami
  }
  colnames(param)<-c("nbnei","p0","mean","sd","skewness","kurtosis","radius")
  
  if(TRUE %in% threeday) {comp <- paste0("_threeday",which(threeday==T))
  } else {comp <- ""}
  save(param,file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,comp,".Rdata"))
}

# Calcul des parametres de la loi empirique selon le voisinage au sens de l'analogie classique
fit.empir.A<-function(rad,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  precip<-get.precip(nbdays,start,end)
  
  N<-length(precip)
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[[rad]]
  
  param<-matrix(NA,ncol=6,nrow=N)
  for (i in 1:N){
    if (i %%50==0) print(i)
    pi<-precip[setdiff(nei[[i]],(i-nbdays+1):(i+nbdays-1))] # on retire les indices des sequences dans lesquelles la journee i intervient
    ni<-length(pi)
    ni1<-sum(pi>0)
    parami<-c(ni,1-ni1/ni,mean(pi[pi>0]),sd(pi[pi>0]),skewness(pi[pi>0]),kurtosis(pi[pi>0])) # remplissage de la ligne: nbre de voisins, proba jour sec, moyenne et ecart-type des pluies non nulles, coefficient d'assymetrie et d'applatissement, rayon du cercle
    param[i,]<-parami
  }
  colnames(param)<-c("nbnei","p0","mean","sd","skewness","kurtosis")
  
  save(param,file=paste0(get.dirstr(k,rean),"fit.empir-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# fit une cdf empirique avec deux niveaux d'analogie: analogie classique puis indicateurs
fit.empir.TL <- function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radAna ="10",radInd="05",CV=TRUE,rean){
  
  # Import des pluies
  precip<-get.precip(nbdays,start,end)
  
  # Import des voisins analogie classique
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[[radAna]]
  
  # Import des indicateurs
  descr<-NULL
  for (i in 1:length(descriptors)) descr<-cbind(descr,get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize,rean)) # recupere les indicateurs
  
  if (nrow(descr) !=length(precip)) stop("Probleme taille des donnees")
  if (nrow(descr) !=length(nei)) stop("Probleme taille des donnees")
  
  N <- length(precip)
  rad1 <- str2prop(radAna)
  rad2 <- str2prop(radInd)
  radtype <- paste0("nrn",rad2/rad1*100)
  
  param<-matrix(NA,ncol=7,nrow=N)
  ind <- vector("list",N)
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){
      tmp<-get.closest(1,descr[nei[[i]],],precip[nei[[i]]],CV=FALSE,nbdays,radtype) # le premier argument est la position du jour dans descr => toujours 1 dans nei
      tmp$idx <- nei[[i]][tmp$idx] # on remet les indices dans le bon referentiel
      if (CV) tmp$idx<-setdiff(tmp$idx,(i-nbdays+1):(i+nbdays-1)) # post cross-validation, impossible en l'etat dans get.closest
      ind[[i]] <- tmp$idx
      
      pi<-precip[tmp$idx]
      ni<-length(pi)
      ni1<-sum(pi>0)
      parami<-c(ni,1-ni1/ni,mean(pi[pi>0]),sd(pi[pi>0]),skewness(pi[pi>0]),kurtosis(pi[pi>0]),tmp$radius) # remplissage de la ligne: nbre de voisins, proba jour sec, moyenne et ecart-type des pluies non nulles, coefficient d'assymetrie et d'applatissement, rayon du cercle
    }
    else parami<-rep(NA,7) # si un des deux indicateurs est nul, ligne du tableau parametre vide
    param[i,]<-parami
  }
  colnames(param)<-c("nbnei","p0","mean","sd","skewness","kurtosis","radius")
  
  fit.emp <- list(param = param,ind = ind)
  save(fit.emp,file=paste0(get.dirstr(k,rean),"fit.empir.TL",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radAna,"_",radInd,".Rdata"))
}

# fit une cdf empirique avec deux analogies (classique et indicateurs) combinees (p=1: analogie classique; p=0: analogie indicateurs)
fit.empir.comb <- function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,p=1,rad="05",CV=TRUE,rean){
  
  # Import des pluies
  precip<-get.precip(nbdays,start,end)
  N <- length(precip)
  
  # Import des scores analogie classique
  load(paste0(get.dirstr(k,rean),"save.score.A/score_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  # Import des indicateurs et calcul des distances euclidiennes
  descr<-NULL
  for (i in 1:length(descriptors)) descr<-cbind(descr,get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize = TRUE,rean)) # Import des indicateurs standardises
  mat.dist <- as.matrix(dist(descr, method = "euclidean"))
  pos.na <- apply(descr,1,function(x) any(is.na(x)))
  mat.dist[pos.na,] <- mat.dist[,pos.na] <- NA
  gc()
  
  # Verification taille
  if (nrow(descr) != length(precip)) stop("Probleme taille des donnees")
  if (nrow(descr) != nrow(mat.dist)) stop("Probleme taille des donnees")
  
  # Calcul de la distance combinee et ajustement empirique
  param<-matrix(NA,ncol=7,nrow=N)
  ind <- vector("list",N)
  
  for (i in 1:N){
    if (i %%50==0) print(i)
    if(!any(is.na(descr[i,]))){
    
      dist_ana <- score[[i]]/sd(score[[i]],na.rm=TRUE) # on normalise les deux distances pour qu'elles aient le meme poids dans dist_tot
      dist_ind <- mat.dist[i,]/sd(mat.dist[i,],na.rm=TRUE)
    
      dist_tot <- p * dist_ana + (1-p) * dist_ind
      pos <- sort(dist_tot, index.return = TRUE, na.last = TRUE) # NA a la fin (on ne peux pas prendre la premiere journee)
      idx <- pos$ix[1:(str2prop(rad)*N)]
      radius <- pos$x[str2prop(rad)*N]
      if (CV) idx <- setdiff(idx,(i-nbdays+1):(i+nbdays-1)) 
      ind[[i]]  <- idx
      
      pi<-precip[idx]
      ni<-length(pi)
      ni1<-sum(pi>0)
      parami<-c(ni,1-ni1/ni,mean(pi[pi>0]),sd(pi[pi>0]),skewness(pi[pi>0]),kurtosis(pi[pi>0]),radius) # remplissage de la ligne: nbre de voisins, proba jour sec, moyenne et ecart-type des pluies non nulles, coefficient d'assymetrie et d'applatissement, rayon du cercle
    }
    else parami<-rep(NA,7) # si un des deux indicateurs est nul, ligne du tableau parametre vide
    param[i,] <- parami
  }
  
  colnames(param)<-c("nbnei","p0","mean","sd","skewness","kurtosis","radius")
  fit.emp <- list(param = param,ind = ind)
  save(fit.emp,file=paste0(get.dirstr(k,rean),"fit.empir.comb",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_p",p,"_rad",rad,".Rdata"))
  
}

# Ajustement d'une loi en p0, a partir de p0 calculee dans fit.empir
fit.loglik.p0<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean){
  
  load(file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  p0<-param[,"p0"]
  
  precip<-get.precip(nbdays,start,end)
  israin<-(precip>0)
  
  loglik.i<-(1-israin)*log(p0) + israin * log(1-p0)
  loglik.i[p0==0 & precip>0] <- 0 #proba de 1 d'avoir de la pluie. Avec le log, fait des NA car (1-israin)*log(p0)=0*-Inf.
  
  fit<-list(p0=p0,loglik.i=loglik.i)
  
  save(fit,file=paste0(get.dirstr(k,rean),"fit.loglik.p0",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),".Rdata"))
}

# Ajustement d'une loi en p0, a partir de p0 calculee dans fit.empir.A
fit.loglik.p0.A<-function(rad,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  print(paste0(get.dirstr(k,rean),"fit.empir-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  load(file=paste0(get.dirstr(k,rean),"fit.empir-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  p0<-param[,"p0"]
  
  precip<-get.precip(nbdays,start,end)
  israin<-(precip>0)
  
  loglik.i<-(1-israin)*log(p0) + israin * log(1-p0)
  loglik.i[p0==0 & precip>0] <- 0 #proba de 1 d'avoir de la pluie. Avec le log ça fait des NA car (1-israin)*log(p0)=0*-Inf.
  
  fit<-list(p0=p0,loglik.i=loglik.i)
  
  save(fit,file=paste0(get.dirstr(k,rean),"fit.loglik.p0-CV.A/",rad,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Sort la matrice de donnees de dat500 (k=1) ou dat1000 (k=2) dans notre fenetre d'analogie pour le jour j (de 1=01/01/1851 a 58804=31/12/2011)
getdata<-function(k,day0,day1=day0,rean,large_win=F){
  num0<-date_to_number(nc[[k]],day0,rean)
  num1<-date_to_number(nc[[k]],day1,rean)
  N<-length(num0:num1)
  infowind<-getinfo_window(k,large_win)
  ncvar_get(nc[[k]],varid="hgt",start=c(infowind[1,1],infowind[2,1],num0),count=c(infowind[1,2],infowind[2,2],N))
}

# Part de la serie de dates complete et ressort seulement le subset voulu
getdates<-function(start="1950-01-01",end="2011-12-31"){
  load(file=paste0("2_Travail/20CR/Data/Membre_",member,"/dates.Rdata"))
  dates[which(dates==start):which(dates==end)]
}

# Date associee a chaque pas de temps dans le fichier NetCDF
getdays<-function(nc,rean){
  if(rean == "20CR") orig <- "1800"
  if(rean == "ERA20C" | rean =="ERA20C_18") orig <- "1900"
  substr(as.POSIXct(nc$dim$time$vals*3600,origin=paste0(orig,'-01-01 00:00'),tz="GMT"),1,10) #format "YYYY-MM-DD"
}

# Importe la liste contenant la distance souhaitee
getdist<-function(k,dist,start="1950-01-01",end="2011-12-31",rean,threeday=FALSE){
  load(file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",ifelse(threeday,"3day_",""),dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
  gc()
  return(dist.list)
}

# Recupere le vecteur distance correspondant au jour voulu
getdist4i<-function(i,dist.vec,N,sU){
  if (i==N) di<-dist.vec[sU[1:(i-1)]+(i-1):1] # si i=N, on prend comme reference le changement de journee defini par sU auquel on ajoute le nombre de jours pour timber sur l'analogie avec N
  else if (i==1) di<-dist.vec[1:(N-i)] # si i=1, on prend le debut du vecteur distance
  else {
    di<-dist.vec[c(sU[1:(i-1)]+(i-1):1,sU[i]+1:(N-i))] # si 1<i<N, on va chercher les analogies sur les autres dates, auxquelles on concatene les analogies de la date i
  }
  if (i==1) di<-c(0,di) # si i=1, le score nul est en premiere position
  else if (i==N) di<-c(di,0) # si i=N, le score nul est en derniere position
  else di<-c(di[1:(i-1)],0,di[i:(N-1)]) # si 1<i<N, le score nul est a la position i (la taille du precedent di etait de N-1, maintenant de N)
  gc()
  
  return(di)
}

# Definit la fenetre spatiale d'analogie pour extraire des donnees de data500 et data1000 
getinfo_window<-function(k,large_win = F){ 
  # lat et lon correspondent aux degres de longitude et latitude pour lesquels on a les donnees des geopot 500 et 1000
  lon<-nc[[1]]$dim$lon$vals # de -30 a 50 deg E tranches de 2 deg -> 41 valeurs
  lat<-nc[[1]]$dim$lat$vals #de 24 a 72 deg N tranches de 2 deg -> 25 valeurs
  # time=data500$dim$time$vals # 58804 valeurs par tranches de 24h (du 01-01-1851 au 31-12-2011), valeurs a 9h chaque jour
  # parametre fenetre d'analogie
  
  c_lon<-6 # centre de la fenetre longitude
  c_lat<-44 # centre de la fenetre latitude
  
  if(k==1 | large_win==T) {# 500hPA
    d_lon<-32 # on prend 32 en longitude pour 500hPa
    d_lat<-16 # on prend 16 en latitude pour 500hPa
  }
  if (k==2 & large_win==F){ # 1000 hPA
    d_lon<-16 # on prend 16 en longitude pour 1000hPa
    d_lat<-8  # on prend 8 en latitude pour 1000hPa
  }
  
  deb_lon <- which.min(abs(lon - (c_lon-d_lon/2)))
  fin_lon <- which.min(abs(lon - (c_lon+d_lon/2)))
  len_lon <- fin_lon - deb_lon + 1
  
  deb_lat <- which.min(abs(lat - (c_lat-d_lat/2)))
  fin_lat <- which.min(abs(lat - (c_lat+d_lat/2)))
  len_lat <- fin_lat - deb_lat + 1
  
  infolon <- c(deb_lon,len_lon) # indice du 1er point de grille et nbre de points de grille
  infolat <- c(deb_lat,len_lat)
  
  return(rbind(infolon,infolat))# matrice 2x2 
  
} 

# Ressort le nom du geopotentiel en fonction de k
getZstr<-function(k){
  if (k==1) return("Z500")
  if (k==2) return("Z1000")
  if (k==3) return("Z500+1000")
}

# Renvoie les indices des analogue souhaites par rapport a ref ainsi que leur score
get.ana<-function(date,rank=c(1,50,100),ref="1851-01-01",k,dist,nbdays=1,start="1950-01-01",end="2011-12-31",rean){
  
  # Indice des voisins
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[["05"]]
  
  pos <- which(getdates()==date)
  res <- vector("list",length = 2)
  names(res) <- c("ind","score")
  res[[1]] <- nei[[pos]][rank+1] + get.delta(ref = ref,start = start) # rank +1 pour enlever le premier qui est la date cible
  
  # Score des voisins
  dist.list <- score.to.mat(k,dist,nbdays,start,end,rean)
  vec <- dist.list[,pos]
  res[[2]] <- sort(vec)[rank+1]
  
  res
}

# Ressort les indices des journees dans le voisinage du jour cible (dans le plan des indicateurs)
get.closest<-function(i,descr,precip,CV=TRUE,nbdays=3,radtype){
  radstr<-substr(radtype,4,nchar(radtype))
  if (substr(radtype,1,3)=="fix"){ # si le rayon demande est fixe (une valeur)
    radius<-get.radius(descr,as.numeric(radstr))
    m0<-matrix(rep(descr[i,],nrow(descr)),ncol=2,byrow=TRUE) # repetition des descripteurs du jour cible
    d<-apply((m0-descr)^2,1,sum)
    idx<-which(d<=radius^2) # idx = indice des jours dans le rayon du cercle
  }
  
  if (substr(radtype,1,3)=="nrn"){ # si le rayon demande n'est pas fixe (un pourcentage de voisins)
    idef<-which(apply(descr,1,function(v) all(!is.na(v))))
    n<-nn2(descr[idef,],t(descr[i,]),str2prop(radstr)*length(precip)) # Nearest Neighbour Search
    radius<-max(n$nn.dists) # le rayon est la distance maximum selectionnee dans n
    idx<-sort(idef[n$nn.idx]) # idx = indice des jours dans le rayon du cercle
  }
  
  if (CV) idx<-setdiff(idx,(i-nbdays+1):(i+nbdays-1)) # si Cross-Validation, on retire les indices des sequences de nbdays jours dans lesquelles le jour i intervient
  pi<-precip[idx]
  
  
  #ni0<-sum(pi>0) # nombre de jours pluvieux dans le voisinage selectionne
  #if (ni0<20) {
  #  #we increase the radius to embrace >=20 positive rainfall
  #  ipos<-which(precip>0 & !is.na(descr[,1]) & !is.na(descr[,2])) # selection des pluies positives avec indicateurs non nuls
  #  if (CV) ipos<-setdiff(ipos,(i-nbdays+1):(i+nbdays-1)) # si cross-validation, on enleve les sequences incluent le jour i
  #  n<-nn2(descr[ipos,],t(descr[i,]),20) # recherche de 20 voisins parmi les jours pluvieux
  #  radius<-max(n$nn.dists)
  #  #idx<-which(d<=radius^2)
  #  idx<-sort(ipos[n$nn.idx]) # nouveau voisinage
  #}
  
  lclose<-list(idx=idx,radius=radius,idef=idef)
  return(lclose)
}

# Retourne une chaine de caractere "-CV" pour le chemin des dossiers
get.CVstr<-function(CV){
  if (CV) return("-CV")
  else return("")
}

# Calcule le delta pdt a appliquer pour ramener les indices a une reference anterieure
get.delta <- function(ref = "1949-12-31", start = "1950-01-01"){
  
  res <- length(seq(as.Date(ref),as.Date(start),by="days")) - 1
  res
  
}

# Import et mise en forme des indicateurs: calculs supplementaires, moyenne glissante sur trois jours
get.descriptor<-function(descriptor,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,rean,threeday=FALSE){
  
  if(rean == "20CR" && end != "2011-12-31"){
    load(file=paste0(get.dirstr(k,rean),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_2011-12-31.Rdata"))
    end_diff <- length(seq(as.Date(end),as.Date("2011-12-31"),"days"))
    criteria <- criteria[1:(nrow(criteria) - end_diff + 1),]
    
  } else {load(file=paste0(get.dirstr(k,rean),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))}
  
  if (descriptor=="celcelnei") descr<-criteria[,"cel"]/criteria[,"celnei"] # division de la celerite du jour par la celerite moyenne des voisins
  else if (substr(descriptor,1,5)=="rsing" & substr(descriptor,nchar(descriptor)-2,nchar(descriptor))!="nei") {
    nb<-substr(descriptor,6,nchar(descriptor))
    descr<-criteria[,paste0("sing",nb)]/criteria[,paste0("q",nb)] # calcul de rsing (local dimension), en divisant la singularite par le quantile souhaite
  }
  #else if (substr(descriptor,1,6)=="rlsing") {
  #  nb<-substr(descriptor,7,nchar(descriptor))
  #  descr<-criteria[,paste0("lsing",nb)]-log(criteria[,paste0("q",nb)]) # calcul de rlsing (log local dimension)
  #}
  else descr<-criteria[,descriptor] # sinon, on selectionne tout simplement la colonne specifiee dans descriptor
  
  if (nbdays>1 && !threeday) descr<-rollapply(descr,width=nbdays,FUN=mean) # moyenne glissante de l'indicateur sur nbdays jours
  
  if (standardize) return(descr/sd(descr,na.rm=TRUE)) # si standardise, on divise l'indicateur par l'ecart type de sa serie
  else return(descr)
}

# Modifie le repertoire d'ecriture des resultats selon k et overall/seasonal 
get.dirstr<-function(k=NULL,rean){
  if (seasonal) {dirstr<-paste0("2_Travail/",rean,"/Rresults/seasonal/")
  } else dirstr<-paste0("2_Travail/",rean,"/Rresults/overall/")
  if (!is.null(k)) dirstr<-paste0(dirstr,"k",k,"/")
  dirstr
}

# Renvoie les indices des nbre events extremes de precip a partir d'une certaine reference
get.ind.extr <- function(nbre, ref = "1950-01-01", nbdays=3, start="1950-01-01", end="2011-12-31", nei=FALSE, bv="all"){
  
  # Classement
  precip <- get.precip(nbdays, start, end, bv)
  tri <- sort(precip, decreasing = TRUE, index.return = TRUE)
  ind <- tri$ix[1:nbre]
  
  # Si chevauchement de deux jours (dt=1), on enleve les voisins et on prend d'autres candidats
  if(nei){
    cond = 0.1 # init
    tmp  = 0.1 # init
    while(cond>0){
      if(cond!=0.1) ind <- c(ind,tri$ix[(nbre+1+tmp):(nbre+tmp+cond)])
      pos <- NULL
      for(i in 1:length(ind)) pos <- c(pos,ifelse(match(ind+1,ind)[i]>i,match(ind+1,ind)[i],i))    
      pos <- sort(na.omit(pos))
      if(length(pos!=0)) ind <- ind[-pos]
      tmp <- tmp+cond
      cond <- nbre-length(ind)
    }
  }
  
  # Reference
  if(ref != start){
    dt <- get.delta(ref, start)
    ind <- ind + dt
  }
  return(ind)
}

# Renvoie les indices des max annuels ou mensuels de precip
get.ind.max <- function(type="year",nbdays=3,start="1950-01-01", end="2011-12-31",bv="all"){
  
  # Import
  precip <- get.precip(nbdays, start, end, bv)
  dates <- as.Date(getdates(start,end))
  length(dates) <- length(precip)
  
  # Traitement
  ag <- as.Date(format(dates,ifelse(type=="year","%Y-01-01","%Y-%m-01")))
  pos <- aggregate(precip,by=list(ag),which.max)
  pos <- pos[,2] + as.numeric(pos[,1]-pos[1,1])
  pos
}

# Renvoie les indices des min et max d'un descripteur a partir d'une certaine reference
get.ind.min.max.descr <- function(descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31"){
  
  descr <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,
                          start = start,end = end,standardize = F,rean = rean)
  res <- list()
  res$x <- c(min(descr,na.rm=T),max(descr,na.rm=T))
  res$idx <- c(which.min(descr),which.max(descr))
  res$xdate <- getdates(start,end)[res$idx]
  res
}

# Importation des donnees de precipitation
get.precip<-function(nbdays,start="1950-01-01",end="2011-12-31",bv="all"){
  
  precip <- read.csv(file=paste0("2_Travail/Data/Precip/",ifelse(bv=="all","Isere@Grenoble",bv),"_cum",nbdays,"day_1950-01-01_2011-12-31.csv"))
  precip <- precip[,1]
  
  if(start != "1950-01-01"){
    start_diff <- length(seq(as.Date("1950-01-01"),as.Date(start),"days"))
    precip <- precip[start_diff:length(precip)]
  }
  if(end != "2011-12-31"){
    end_diff <- length(seq(as.Date(end),as.Date("2011-12-31"),"days"))
    precip <- precip[1:(length(precip) - end_diff + 1)]
  }
  return(precip)
}

# Retourne une chaine de caractere "_std" pour le chemin des dossiers
get.stdstr<-function(standardize){
  if (standardize) return("_std")
  else return("")
}

# Importation des WP EDF
get.wp <- function(start="1950-01-01",end="2011-12-31",risk=F){
  
  wp <- read.table("2_Travail/Data/WP/type_temps_jour.txt", quote="\"", comment.char="")
  wp[,1] <- as.character(format(as.Date(wp[,1],"%d/%m/%Y"),"%Y-%m-%d"))
  wp <- wp[which(wp[,1]==start):which(wp[,1]==end),2]
  
  if(risk){
    dates <- getdates(start,end)
    wp[substr(dates,6,7) %in% c("03","04","05","06","07","08")] <- wp[substr(dates,6,7) %in% c("03","04","05","06","07","08")]+8
  }
  wp
}

# Calcule les gradients pour 500 (k=1) ou 1000 (k=2) pour tous les jours
grad<-function(k,day0,day1,rean,large_win=F){
  x=getdata(k,day0,day1,rean,large_win)  #data500 ou data1000, 3 dims
  gradlon=(makegrad(x,1)) #selon lon
  gradlat=(makegrad(x,2)) #selon lat
  return(list(gradlon=gradlon,gradlat=gradlat))
}

# Comparaison des sequences de differents nbdays jours
illustr.precip.seq<-function(nbdays,start="1950-01-01",end="2011-12-31"){
  
  precip<-get.precip(nbdays,start,end)
  cdf<-ecdf(precip)
  precip0<-get.precip(1,start,end)
  cdf0<-ecdf(precip0)
  
  graphics.off()
  
  tmp<-NULL
  n<-length(precip)
  for (i in 1:nbdays) tmp<-cbind(tmp,precip0[(i:(n+i-1))]) 
  tmp<-t(apply(tmp,1,sort))
  
  rep <- round(table(tmp[sort(precip,decreasing=T,index.return=T)$ix[1:62],1] > quantile(precip0,0.9))/62*100,0) # points noirs dessus/dessous q90
  
  png(file=paste0("2_Travail/Rresults/illustr.precip.seq/plotseq_",nbdays,"_",start,"_",end,".png"),width=350,heigh=404,units="px",res=72)
  plot(range(precip0),range(precip0),type="n",xlab=paste0(nbdays,"-day precip (mm/day)"),ylab="daily precip (mm/day)")
  for (i in 1:ncol(tmp)) points(precip,tmp[,i],col=i,pch=19,cex=0.8)
  abline(v=quantile(precip,1-62/length(precip)),lty=2)
  text(x=quantile(precip,1-62/length(precip))*1.1,y=67,"q = 0.9973",font = 3,srt = 90)
  abline(h=quantile(precip0,0.9),lty=2)
  text(x=55,y=quantile(precip0,0.9)*1.2,"q = 0.9",font = 3)
  text(x=65,y=5,paste0(rep[1],"%"),font = 2,cex=1.2)
  text(x=65,y=20,paste0(rep[2],"%"),font = 2,cex=1.2)
  dev.off()
}

# Trace la carte du BV avec pluvios, rivieres, villes
image.region<-function(pluvios = TRUE, sousBV = F){
  
  pdf(file=paste0("2_Travail/Rresults/image.region/map_region",ifelse(pluvios,"_pluvios",""),ifelse(sousBV,"_sousBV",""),".pdf"),width = 7.5,height = 7.5)
  
  # Fond de carte
  load(file=paste("2_Travail/Data/Carto/griddata_1x1_IsereSavoieHautesAlpes.Rdata",sep=""))
  Fx<-griddata$Fx
  Fy<-griddata$Fy
  Fz<-griddata$Fz*1000
  #windows(width=7.5,heigh=6.5)
  image.plot(Fx,Fy,Fz,col=gray(seq(0.1,0.99,length=100)),xlab="X (km) - Lambert II extended",ylab="Y (km) - Lambert II extended",legend.line=-2.3, cex.axis=1.3, cex.lab=1.3)
  
  # Bordure du BV et des sous-BV
  bord<-read.csv("2_Travail/Data/Carto/borders-lambert93-lambertII.csv",sep=";")
  bord<-cbind(bord$XLII.m,bord$YLII.m)/1000
  lines(bord,col="white",lwd=3)
  #geometry::polyarea(bord[,1],bord[,2]) pour calculer la surface!
  #bord <- bord[bord[,1] %in% c(670,19,35,124,248,176,596,104),] # Drac
  #bord <- bord[bord[,1] %in% c(87,576,706,577,797,621,136,180,1092,116,792),] # Isere
  #ui<-unique(bord[,1])
  #for (i in ui) lines(bord[bord[,1]==i,2],bord[bord[,1]==i,3],col="red")
  
  #a1<-as(bord[bord[,1]==unique(bord[,1])[1],2:3],"gpc.poly")
  #a2<-as(bord[bord[,1]==unique(bord[,1])[2],2:3],"gpc.poly")
  #a<-union(a1,a2) # Initialisation puis boucle

  # Tout le BV
  #for (i in ui[-(1:2)]) {
  #  a1<-as(bord[bord[,1]==i,2:3],"gpc.poly")
  #  a<-union(a1,a)
  #} 
  #lines(a@pts[[1]]$x,a@pts[[1]]$y,col="white",lwd=3)
  
  if(sousBV){ # On retrace le Drac dessus
    bord <- bord[bord[,1] %in% c(670,19,35,124,248,176,596,104),] # Drac
    ui<-unique(bord[,1])
    a1<-as(bord[bord[,1]==unique(bord[,1])[1],2:3],"gpc.poly")
    a2<-as(bord[bord[,1]==unique(bord[,1])[2],2:3],"gpc.poly")
    a<-union(a1,a2) # Initialisation puis boucle
    for (i in ui[-(1:2)]) {
      a1<-as(bord[bord[,1]==i,2:3],"gpc.poly")
      a<-union(a1,a)
    } 
    lines(a@pts[[1]]$x,a@pts[[1]]$y,col="white",lwd=3)
  }
  # Rivieres
  riv<-readOGR("2_Travail/Data/Carto/Principales/Principales/riviere2_3.shp")
  id<-c(552,554,553,381)#isere,Drac,Arc,Romanche
  
  for (i in id){
    tmp<-riv@lines[[i]]@Lines[[1]]@coords/1000
    lines(tmp[,1],tmp[,2],col="cyan")
  }
  
  # Pluvios
  if(pluvios){
    load("2_Travail/Data/Carto/Isere@Grenoble_24h_minlength50.Rdata")
    points(l$coord[,1],l$coord[,2],pch=21,bg="red")
  }
  
  # Noms des rivieres
  if(!sousBV){
  shadowtext(938,2074+1,"Isère",col="cyan",srt=28,cex=1.2,bg="darkblue")
  shadowtext(880,2042+1,"Isère",col="cyan",srt=60,cex=1.2,bg="darkblue")
  shadowtext(846,2034+1,"Isère",col="cyan",srt=64,cex=1.2,bg="darkblue")
  shadowtext(934,2035,"Arc",col="cyan",srt=-5,cex=1.2,bg="darkblue")
  shadowtext(907-2.5,2010-2,"Romanche",col="cyan",srt=5,cex=1.2,bg="darkblue")
  shadowtext(904.5,1972+1,"Drac",col="cyan",srt=20,cex=1.2,bg="darkblue")
  
  # Noms des villes
  r<-rbind(
    c(867.137,2026.491,"Grenoble"),
    c(894.062,2107.205,"Annecy"),
    c(897.82,1959.14,"Gap"),
    c(915,2082,"Albertville"),
    c(940.146,2031.530,"Modane")#,
    #c(945.384,2078.211,"Bourg St Maurice")
  )
  shadowtext(r[,1],r[,2],r[,3],pos=c(1,1,1,1,1,3),cex=1.2,col="white",bg="darkblue",adj=c(0,0),r=0.09,font=3)#font=2
  points(r[,1],r[,2],pch=22,col="darkblue",bg="white",cex=.8)
  }
  graphics.off()
}

# Trace la carte avec taille des pluvios par cumul
image.cumul<-function(crue=FALSE){
  
  # Import du fond de carte
  load(file=paste("2_Travail/Data/Carto/griddata_1x1_IsereSavoieHautesAlpes.Rdata",sep=""))
  Fx<-griddata$Fx
  Fy<-griddata$Fy
  Fz<-griddata$Fz*1000
  graphics.off()
  
  # Import des bordures
  bord<-read.csv("2_Travail/Data/Carto/borders-lambert93-lambertII.csv",sep=",")
  bord<-cbind(bord$XLII.m,bord$YLII.m)/1000
  #ui<-unique(bord[,1])
  #
  #a1<-as(bord[bord[,1]==unique(bord[,1])[1],2:3],"gpc.poly")
  #a2<-as(bord[bord[,1]==unique(bord[,1])[2],2:3],"gpc.poly")
  #a<-union(a1,a2) # Initialisation puis boucle
  #for (i in ui[-(1:2)]) {
  #  a1<-as(bord[bord[,1]==i,2:3],"gpc.poly")
  #  a<-union(a1,a)
  #}
  
  # Import des rivieres
  riv<-readOGR("2_Travail/Data/Carto/Principales/Principales/riviere2_3.shp")
  id<-c(552,554,553,381)#isere,Drac,Arc,Romanche
  
  # Import des pluvios et de la pluie de BV
  load("2_Travail/Data/Carto/Isere@Grenoble_24h_minlength50.Rdata")
  precip <- get.precip(nbdays = 3)
  
  # Events extremes de pluie
  if(!crue){
    ind    <- get.ind.extr(10,nei=TRUE)
    event  <- getdates()[ind]
  } else {# Events extremes de crues
    load("2_Travail/Rresults/image.cumul/crues_1969_2011_Isere_StGervais.Rdata")
    crues <- as.data.frame(crues[crues[,4]>0.75,]) # on selectionne les events quadriennaux
    crues <- crues[order(crues[,2],decreasing = T),]
    event <- as.character(crues[,1]-3) # on prend les pluies de j-3,j-2,j-1 de la crue (car 8hj a 7hj+1)
    ind   <- match(event,getdates())
    
    # Quantiles de pluie associes
    distrib <- ecdf(precip[precip>0])
    qua <- round(distrib(precip[ind]),4)
    
    png(file = "2_Travail/Rresults/image.cumul/Quant_rain_flood.png",width = 517,height = 369,units = "px",res = 72)
    plot(qua*100,ylim=c(60,100),col="royalblue",pch=19,xlab="10 largest flow rates 1969-2011",ylab="3-day positive rainfall quantile (%)")
    abline(h=94,col="red",lty=2)
    graphics.off()
    
    prec_ete <- precip[substr(getdates(),6,7) %in% c("05","06","07")]
    distrib <- ecdf(prec_ete[prec_ete>0])
    qua_ete <- round(distrib(precip[ind]),4)
    
    # Histogramme des crues par mois
    png(file = "2_Travail/Rresults/image.cumul/Histogram_flood_month.png",width = 419,height = 322,units = "px",res = 72)
    tmp <- as.numeric(substr(crues[,1],6,7))
    hist(tmp,0:12,axes = FALSE,xlab="Month",col = "royalblue",main="10 largest flow rates 1969 -2011",border = "white")
    lines(c(0,12),c(0,0))
    axis(2)
    axis(1,at = 0.5:11.5,labels = 1:12,tick = FALSE,padj = -1)
    graphics.off()
  }
  
  for(i in 1:length(event)){
    print(paste0("Carte ",i,"/",length(event)))
    pdf(file=paste0("2_Travail/Rresults/image.cumul/Map_",ifelse(crue,"flood_","rainfall_"),"event_",i,".pdf"),width=7.5,heigh=7.5)
    
    # Carte de base
    image.plot(Fx,Fy,Fz,col=gray(seq(0.1,0.99,length=100)),xlab="X (km) - Lambert II extended",ylab="Y (km) - Lambert II extended",
               legend.line=-2.3, cex.axis=1.3, cex.lab=1.3, legend.shrink = 0.5) # fond
    lines(bord,col="white",lwd=3)
    #lines(a@pts[[1]]$x,a@pts[[1]]$y,col="white",lwd=3) # contour
    for (j in id){ # rivieres
      tmp<-riv@lines[[j]]@Lines[[1]]@coords/1000
      lines(tmp[,1],tmp[,2],col="cyan")
    }
    
    title(main=paste0(event[i]," - ",as.Date(event[i])+2," (BV: ",round(precip[ind[i]],1)," mm/day)",
          ifelse(crue,paste0("  -  Qm = ",crues[i,2]," m3/s - TR = ",round(1/(1-crues[i,4]),0)," ans"),"")),
          cex.main=ifelse(crue,0.9,1))

    if(crue){ # Affichage des quantiles
      par(xpd=TRUE)
      if(substr(getdates()[ind[i]],6,7) %in% c("05","06","07")){
        text(x = 915,y = 2125,labels = paste0("q(pos) = ",qua[i],"     q(pos;5,6,7) = ",qua_ete[i]),cex = 0.9,font=3,col="gray20")
      } else{text(x = 915,y = 2125,labels = paste0("q(pos) = ",qua[i]),cex = 0.9,font=3,col="gray20")}
      
    }
    
    # Pluvios
    pos <- length(seq(as.Date("1913-01-01"),as.Date(event[i]),by="day"))
    cum <- unname(cbind(l$coord[,1],l$coord[,2],l$ymat[pos,],l$ymat[pos+1,],l$ymat[pos+2,]))
    cum <- cbind(cum,round(apply(cum[,3:5]*0.1,1,mean),1)) # cumul moyen sur trois jours
    
    points(cum[,1],cum[,2],pch=21,bg="red",cex=cum[,6]/50*2)
    par(xpd=TRUE)
    legend(993,2116,c("25","50","75","100"),col="black",pch=21,pt.cex=c(1,2,3,4),
           pt.bg="red",bty="n",cex=0.8,title = "Cumul (mm/day)",x.intersp = c(1.1,1.2,1.3,1.5),
           y.intersp = c(1.05,1.08,1.15,1.3))
    dev.off()
  }
}

# Charge les fichiers NetCDF 500 hPa et 1000 hPa d'un membre donne
load.nc<-function(rean = "20CR"){
  if(rean == "20CR"){
    nc500<-nc_open(paste0("2_Travail/20CR/Data/Membre_",member,"/20Crv2c_Membre_",member,"_HGT500_1851-2011_daily.nc"))
    nc1000<-nc_open(paste0("2_Travail/20CR/Data/Membre_",member,"/20Crv2c_Membre_",member,"_HGT1000_1851-2011_daily.nc"))
  }
  if(rean == "ERA20C"){
    nc500<-nc_open("2_Travail/ERA20C/Data/ERA20C_HGT500_1900_2010_daily.nc")
    nc1000<-nc_open("2_Travail/ERA20C/Data/ERA20C_HGT1000_1900_2010_daily.nc")
  }
  if(rean == "ERA20C_18"){
    nc500<-nc_open("2_Travail/ERA20C_18/Data/ERA20C_HGT500_1900_2010_18h.nc")
    nc1000<-nc_open("2_Travail/ERA20C_18/Data/ERA20C_HGT1000_1900_2010_18h.nc")
  }
  
  nc<<-list(nc500=nc500,nc1000=nc1000)
}

# Sort la matrice des gradients pour input= data500 ou data1000 sur notre fenetre d'analogie, l=1 longitude, l=2 latitude
makegrad<-function(mat,l){
  dims<-dim(mat)
  if (l==1) gradmat<-mat[2:dims[[1]],,]-mat[1:(dims[[1]]-1),,] # l=1 1ere ligne= difference entre 2eme et 1ere ligne, 2eme ligne=difference entre la 3eme et le 2eme ligne --> une ligne de moins que input
  if (l==2) gradmat<-mat[,2:dims[[2]],]-mat[,1:(dims[[2]]-1),] # l=2 1ere colonne= difference entre 2eme et 1ere colonne, 2eme colonne=difference entre la 3eme et la 2eme colonne  --> une colonne de moins que input
  return(gradmat)
}

# Genere la pluie pour un sous BV
make.precip1.isere<-function(start="1950-01-01",end="2011-12-31") {
  
  bord<-read.csv("2_Travail/Data/Carto/border_Isere@Grenoble.csv",sep=",") #coordonnées des bassins
  
  load("2_Travail/Data/Carto/Isere@Grenoble_24h_minlength50.Rdata") #precip journalières
  
  #bassin<-c(87,567,706,577,797,621,136,180,1092,116,792) #ensemble des bassins de l'Isère
  bassin<-c(670,19,35,124,248,176,596,104) #ensemble des bassins du Drac
  
  #indice des 1er et derniers jours dans le Rdata des precip 
  id1<-which(rownames(l$ymat)==paste0(substr(start,1,4),substr(start,6,7),substr(start,9,10)))
  id2<-which(rownames(l$ymat)==paste0(substr(end,1,4),substr(end,6,7),substr(end,9,10)))
  
  grids<-merge((seq(min(bord[,2]),max(bord[,2]),by=1)),(seq(min(bord[,3]),max(bord[,3]),by=1)))  #on quadrille (pas 1km), pour pouvoir ensuite calculer à l'échelle du domaine en ne prenant que les points du quadrillage qui sont dans le domaine
  
  idx<-NULL
  for (i in 1:length(bassin)) { 
    idx<-c(idx,inpip(grids,bord[bord[,1]==bassin[i],2:3])) #indice des points du quadrillage qui sont dans mon domaine
  }
  #plot(grids[,1],grids[,2])
  #points(grids[idx,1],grids[idx,2],col="red")
  
  ivec<-id1:id2
  precip<-rep(NA,length(ivec))
  for (i in ivec) {
    if ((i-ivec[1]+1) %%50==0) print(i-ivec[1]+1)
    fit<- Tps(l$coord[,1:2],l$ymat[i,]/10,give.warnings=F)
    tmp<-predict(fit,cbind(grids[idx,1],grids[idx,2]))  #contient toutes les estimations des points du quadrillage qui sont dans le domaine
    tmp[tmp<0]=0  #on remet à 0 toutes les valeurs estimées de pluies qui sont négatives
    precip[i-ivec[1]+1]<-mean(tmp,na.rm = T) # precip contient la valeur moyenne du bassin pour chaque jour entre debut et fin
  }
  
  precip
}

# Carte de geopotentiel d'un jour donne
map.geo <- function(date,rean,k,save=F){
  
  # Import des donnees
  load.nc(rean)
  if(k==1){nc <- nc$nc500
  }else{nc <- nc$nc1000}
  
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  num0<-date_to_number(nc,date,rean)
  
  geo <- ncvar_get(nc,varid="hgt",start=c(1,1,num0),count=c(length(lon),length(lat),1))
  title <- ifelse(k==1,"500 hPa","1000 hPa")
  zlim <- ifelse(rep(k==1,2),c(4800,6100),c(-400,500))
  
  # Carte
  if(save) png(filename = paste0("2_Travail/20CR/Rresults/overall/k",k,"/map.geo/",date,"_k",k,".png"))
  image.plot(lon,lat,geo,xlim=c(-20,25),ylim=c(25,70),asp=1,zlim=zlim,
             col=rev(brewer.pal(n = 11, name = "RdBu")),
             xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(date," - ",title),
             legend.line=-2.3, cex.axis=1.3, cex.lab=1.3, cex.main=1.3)
  
  data(wrld_simpl)
  plot(wrld_simpl, add = TRUE)
  box()
  if(save) graphics.off()
  
}

# Carte de geopotentiels des sequences de plus fortes precipitations, avec distribution des analogues
map.extr <- function(k,N="02",start,end,rean){
  
  if(N=="02") n <- "0.2"
  if(N=="05") n <- "0.5"
  
  # Import des donnees
  precip <- get.precip(3,start,end)
  dates <- as.Date(getdates(start,end))
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_TWS_member",member,"_k",k,"_mean",3,"day_",start,"_",end,".Rdata"))
  neiTWS<-nei[[N]]
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_RMSE_member",member,"_k",k,"_mean",3,"day_",start,"_",end,".Rdata"))
  neiRMSE<-nei[[N]]
  
  # Cartes
  config <- c("62_max","annual_max")
  q90 <- quantile(precip[precip>0],probs=0.9)
  
  for(j in 1:length(config)){
    
    if(j==1) ind <- sort(precip,decreasing = T,index.return=T)$ix[1:62]
    if(j==2) ind <- get.ind.max(type = "year",start = start,end = end)
    dates.extr <- dates[ind]
    
    pdf(file = paste0("2_Travail/20CR/Rresults/overall/k",k,"/map.extr/map_",config[j],"_k",k,"_",N,"_",start,"_",end,".pdf"),width = 13,height = 8)
    layout(matrix(1:12,3,4,byrow = T))
    par(mar=c(4.5,5,4,6.5))
    print(paste0(length(dates.extr)," max"))
    
    for(i in 1:length(dates.extr)){
      
      print(paste0(i,"/",length(dates.extr)))
      
      # les 3 cartes
      map.geo(dates.extr[i],rean,k); map.geo(dates.extr[i]+1,rean,k); map.geo(dates.extr[i]+2,rean,k)
      
      # la distribution des analogues
      precip.i <- precip[ind[i]]
      neiTWS.i <- neiTWS[[ind[i]]][-1]
      neiRMSE.i <- neiRMSE[[ind[i]]][-1]
      
      plot(ecdf(precip[neiTWS.i]),xlim=c(0,60),lwd=2,xlab="Precipitation (mm/d)",
           main=paste0(length(neiTWS.i)," analog sequences (",n,"%)"))
      lines(ecdf(precip[neiRMSE.i]),col="purple")
      abline(v=precip.i,col="darkblue",lwd=2)
      text(precip.i-3,0.8,"obs",col="darkblue",font=2,srt=90)
      abline(v=q90,col="red",lwd=2)
      text(q90-3,0.8,"q90",col="red",font=2,srt=90)
      legend(ifelse(precip.i<40,40,20),0.2,c("TWS","RMSE"),col=c("black","purple"),lwd=2,pch=19,bty="n",cex=0.8)
      
    }
    graphics.off()
  }
}

# Noms propres des couples d'indicateurs
nam2str<-function(nams,cloud=FALSE){

  for(i in 1:length(nams)){
    if(nams[i] == "snei") nams[i] <- "pers"
    if(nams[i] == "sing05_2") nams[i] <- "singNew"
    if(nams[i] == "sing05_2_500") nams[i] <- "singNew500"
    if(nams[i] == "sing05_2nei") nams[i] <- "singneiNew"
    if(nams[i] == "celnei_2") nams[i] <- "celneiNew"
    if(nams[i] == "rsingnei_best") nams[i] <- "rsingnei-best"
    if(nams[i] == "TWSgeo") nams[i] <- "TWS_500_1000"
  }
  
  if(!cloud && length(nams)==2 && nams != c("A","A")){
    if(nams[1]==nams[2]) nams <- nams[1]
  }
  nams
}

# Ajout d'une chaine de caractere a l'indicateur (05,1,2,5,10 a tous les indicateurs sauf cel)
paste.descr<-function(descr,str){
  
  cond <- c(substr(descr,1,3) == "cel",
            substr(descr,1,3) == "acc",
            substr(descr,1,8) == "sing05_2",
            substr(descr,1,8) == "rsingnei",
            substr(descr,1,6) =="celnei",
            descr == "singnei",
            descr == "persnei",
            descr == "TWSgeo")
  
  cond <- matrix(data = cond,nrow = length(cond)/2,ncol = 2,byrow = TRUE)
  descr[apply(cond,2,sum)==0] <- paste0(descr[apply(cond,2,sum)==0],str)
  descr
}

# plot les analogues du plus/moins attracteur pour expliquer difference sing et rsing
plot.ana<-function(){
  pos <- get.ind.min.max("rsing05")
  ana_min <- get.ana(date = getdates()[pos$idx[1]],rank=1:112,ref = "1950-01-01",
                     k = 1,dist = "TWS",nbdays = 1,start,end,rean)$score
  ana_max <- get.ana(date = getdates()[pos$idx[2]],rank=1:112,ref = "1950-01-01",
                     k = 1,dist = "TWS",nbdays = 1,start,end,rean)$score
  
  for(i in 1:3){
    png(filename = paste0("2_Travail/20CR/Rresults/overall/k1/plot.ana/plot_ana_",i,".png"),width = 610,height = 390,units = "px")
    par(bg="grey")
    plot(c(0,0.5),c(0,0),type="n",axes=FALSE,xlab="",ylab="",ylim=c(-40,100))
    points(0,0,pch=4)
    points(ana_min,rep(0,112),col="red",pch=4)
    points(ana_max,rep(0,112),col="royalblue",pch=4)
    axis(1,pos = -10)
    text(0.25,-40,"TWS score")
    if(i!=1){
      abline(v=mean(ana_min),lty=2,lwd=2)
      abline(v=mean(ana_max),lty=2,lwd=2)
      arrows(mean(ana_min),80,mean(ana_max),80,code = 3,length = 0.1)
      text(0.28,90,expression(paste(Delta,"sing")))}
    if(i==3){
      lines(density(ana_min))
      lines(density(ana_max))
      arrows(mean(ana_min)+0.022,15,mean(ana_max)-0.01,15,code = 3,length = 0.1)
      text(0.28,25,expression(paste(Delta,"rsing")))
    }
    graphics.off()
  }
}

# plot le bilan des valeurs des indicateurs pluies fortes & sequences seches QUALITATIF
plot.bilan<-function(k,dist,val=c(5,5,5,5),type="fort"){
  
  descr <- c("celnei","persnei","singnei","rsingnei")
  
  png(filename = paste0("2_Travail/Rresults/plot.bilan/plot_bilan_k",k,"_",dist,"_",type,".png"),width = 400,height = 400,units = "px")
  par(pty="s")
  
  plot(c(1,4),c(1,10),type="n",ylim=c(1,10),pch=19,xaxt="n",yaxt="n",xlab="",ylab="")
  axis(side=1, at=1:4, labels = descr)
  axis(side=2, at=1:10, labels = 1:10)
  grid()
  points(1:4,val,pch=19,cex=1.5)
  title(paste0(ifelse(k==1,"500","1000"),"hPa - ",dist))
  graphics.off()
}

# plot le bilan des valeurs des indicateurs pluies fortes & sequences seches QUANTITATIF
plot.bilan.quant <- function(nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  comb  <- c("500hPa - TWS","500hPa - RMSE","500hPa - RMSE_I","500hPa - RMSE_II",
             "1000hPa - TWS","1000hPa - RMSE","1000hPa - RMSE_I","1000hPa - RMSE_II")
  descr <- c("celnei","persnei","singnei","rsingnei","accnei")
  
  precip <- get.precip(nbdays,start,end)
  
  ind_min <- which(precip==0)
  ind_max <- sort(precip,decreasing=T,index.return=T)$ix[1:62]
  #ind_max <- which(precip>quantile(precip[precip>0],probs=0.94))
  ind=list(ind_min,ind_max)
  
  for(i in 1:length(ind)){
    png(filename = paste0("2_Travail/20CR/Rresults/overall/plot.bilan.quant/boxplot_",ifelse(i==1,"p0","pmax"),"_member",member,"_",nbdays,"day.png"),width = 1200,height = 800,units = "px")
    layout(matrix(1:length(comb),2,length(comb)/2,byrow = T))
    
    for(j in comb){
      k <- ifelse(substr(j,1,3)=="500",1,2)
      if(substr(j,nchar(j)-2,nchar(j))=="TWS") dist <- "TWS"
      if(substr(j,nchar(j)-3,nchar(j))=="RMSE") dist <- "RMSE"
      if(substr(j,nchar(j)-5,nchar(j))=="RMSE_I") dist <- "RMSE_I"
      if(substr(j,nchar(j)-6,nchar(j))=="RMSE_II") dist <- "RMSE_II"
      
      tab.descr <- matrix(NA,length(getdates(start,end))-nbdays+1,length(descr))
      colnames(tab.descr) <- descr
      ind.descr <- NULL
      
      for(l in descr){
        tab.descr[,l] <- get.descriptor(descriptor = l,k = k,dist = dist,nbdays = nbdays,start = start,end = end,
                                        standardize = F,rean = rean)
      }
      tab.qua <- apply(tab.descr,2,function(x) y=(rank(x)-1)/(length(x)-1)*100)
      
      #png(filename = paste0("2_Travail/20CR/Rresults/overall/plot.bilan.quant/boxplot_",ifelse(i==1,"p0_","pmax_"),dist,"_member",member,"_k",k,"_mean",nbdays,"day.png"),width = 400,height = 400,units = "px")
      boxplot(tab.qua[ind[[i]],],ylim=c(0,100),outline=F,col=ifelse(i==1,"mistyrose","lightsteelblue1"),ylab="Quantile (%)",type="n")
      grid()
      title(j)
    }
    graphics.off()
  }
}
  
# plot la comparaison des p0 et mean(p>0) par mois entre la climato, l'analogie classique et les indicateurs
plot.clim<-function(rean,k,descriptorsPos,descriptorsp0,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",str="05",CV=TRUE){
  
  # Import climato
  precip <- get.precip(nbdays,start,end)
  precip[precip<0.01] <- 0
  
  # Import de l'analogie classique
  load(paste0(get.dirstr(k[1],rean),"fit.empir",get.CVstr(CV),".A/",str,"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  paramAna <- param[,c("p0","mean")]
  
  # Import pluies positives de descriptorsPos
  load(paste0(get.dirstr(k[1],rean),"fit.empir",get.CVstr(CV),"/",descriptorsPos[1],"_",descriptorsPos[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata"))
  paramPos <- param[,"mean"]
  
  # Import p0 de descriptorsp0
  load(paste0(get.dirstr(k[2],rean),"fit.empir",get.CVstr(CV),"/",descriptorsp0[1],"_",descriptorsp0[2],"_",dist[2],"_member",member,"_k",k[2],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata"))
  paramp0 <- param[,"p0"]
  rm(param)
  
  # Boxplot pluies positives
  tab <- data.frame(Date = getdates(start,as.Date(end)-nbdays+1),
                    Obs  = precip,
                    Ana  = paramAna[,"mean"],
                    Ind  = paramPos)
  
  tab$Month <- substr(tab$Date,6,7)
  barplot(height = aggregate(tab$Obs,by=list(tab$Month),FUN=function(x) sum(x)/62)[,2],
          col="royalblue",names.arg = month.abb,ylim=c(0,120),xlab="Month",ylab="Cumul (mm/day")
  
  tab[tab$Obs==0,] <- NA # on ne prend que les jours pluvieux
  tab <- gather(data = tab,key = "Method",value = "Rain",2:4)
  tab$Method <- factor(tab$Method, levels = c("Obs","Ana","Ind"))
  tab <- na.omit(tab)
  
  ggplot(tab, aes(x=Month, y=Rain, fill=Method)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,3,1.5,1.5),"cm"),legend.key.size = unit(2,"cm"),plot.title = element_text(hjust = 0.5,vjust = 5))+
    geom_boxplot()+
    stat_summary(fun.y = mean, geom="point",colour="black", size=2,position=position_dodge(width=0.75)) +
    ylim(0,10)+
    ylab("Rain (mm/day)")+
    scale_x_discrete(labels=month.abb)+
    scale_fill_manual(values=c("royalblue","orange","red"))+
    geom_hline(yintercept = 0,linetype="dashed")+
    ggtitle("Pluies positives observées et modelisées")+
  
  ggsave(paste0("2_Travail/20CR/Rresults/overall/plot.clim/boxplot_pos_",descriptorsPos[1],"_",descriptorsPos[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"), width = 16, height = 9, dpi = 100)
  graphics.off()
  
  # Boxplot p0
  tab <- data.frame(Date = getdates(start,as.Date(end)-nbdays+1),
                    Obs  = precip,
                    Ana  = paramAna[,"p0"],
                    Ind  = paramp0)
  
  tab$Month <- substr(tab$Date,6,7)
  clim <- aggregate(tab$Obs,by=list(tab$Month),FUN=function(x) sum(x==0)/length(x))
  colnames(clim) <- c("Month","p0")
  tab <- gather(data = tab,key = "Method",value = "p0",3:4)
  tab <- na.omit(tab)
  
  ggplot(tab, aes(x=Month, y=p0)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,3,1.5,1.5),"cm"),legend.key.size = unit(2,"cm"),plot.title = element_text(hjust = 0.5,vjust = 5))+
    geom_boxplot(aes(fill=Method))+
    stat_summary(aes(group=Method),fun.y = mean, geom="point",colour="black", size=2,position=position_dodge(width=0.75)) +
    ylim(0,0.5)+
    ylab("p0")+
    scale_x_discrete(labels=month.abb)+
    scale_fill_manual(values=c("orange","red"))+
    geom_hline(yintercept = 0,linetype="dashed")+
    ggtitle("Probabilité de séquence sèche")+
    geom_line(data = clim,aes(x = Month, y = p0,group = 2),col="royalblue",size=2)
    
    ggsave(paste0("2_Travail/20CR/Rresults/overall/plot.clim/boxplot_p0_",descriptorsp0[1],"_",descriptorsp0[2],"_",dist[2],"_member",member,"_k",k[2],"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"), width = 16, height = 9, dpi = 100)
  graphics.off()
  
  
}

# plot 3D d'un nuage de points representant le plan des geopotentiels
plot.cloud<-function(){
  
  # Non Uniforme
  x <- rnorm(2264,1,1/3) # quasi 100% des valeurs entre mu-3sigma et mu+3sigma donc diametre boule de 2 (rayon de 1)
  y <- rnorm(2264,1,1/3)
  z <- rnorm(2264,1,1/3)
  mat <- cbind(x,y,z)
  sing <- rsing <- vector(length=nrow(mat))
  
  for(i in 1:nrow(mat)){
  count <- nn2(data = mat[-i,],query = t(mat[i,]),k = nrow(mat)-1)$nn.dists
  sing[i] <- mean(count)
  rsing[i] <- mean(count)/count[length(count)]
  }
  
  ampl_sing <- range(sing)
  ampl_rsing <- range(c(rsing,0.7))
  
  pdf("2_Travail/Rresults/plot.cloud/plot_sing.pdf",width = 7.5,height = 7.5)
  points3D(mat[,1],mat[,2],mat[,3],cex=0.6,pch=19,bty="g",colvar=sing,main="all sing")
  graphics.off()
  pdf("2_Travail/Rresults/plot.cloud/plot_rsing.pdf",width = 7.5,height = 7.5)
  points3D(mat[,1],mat[,2],mat[,3],cex=0.6,pch=19,bty="g",colvar=rsing,clim=ampl_rsing,main="all rsing")
  graphics.off()
  
  # Uniforme
  mat <- runiform_ball(2264,3,1)
  sing <- rsing <- vector(length=nrow(mat))
  
  for(i in 1:nrow(mat)){
    count <- nn2(data = mat[-i,],query = t(mat[i,]),k = nrow(mat)-1)$nn.dists
    sing[i] <- mean(count)
    rsing[i] <- mean(count)/count[length(count)]
  }
  
  pdf("2_Travail/Rresults/plot.cloud/plot_sing_uni.pdf",width = 7.5,height = 7.5)
  points3D(mat[,1],mat[,2],mat[,3],cex=0.6,pch=19,bty="g",colvar=sing,clim=ampl_sing,main="all sing")
  graphics.off()
  pdf("2_Travail/Rresults/plot.cloud/plot_rsing_uni.pdf",width = 7.5,height = 7.5)
  points3D(mat[,1],mat[,2],mat[,3],cex=0.6,pch=19,bty="g",colvar=rsing,clim=ampl_rsing,main="all rsing")
  graphics.off()
  
}

# plot la comparaison des CRPS indicateurs, analogie classique et climato pour differents quantiles de pluie
plot.crps<-function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",str="05",CV=TRUE){
  
  # Import des scores ind, ana, climato
  load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",descriptors[1],"_",descriptors[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata"))
  crps_ind <- crps[,"pos ecdf"]
  
  load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/",str,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  #load(file=paste0("2_Travail/",rean,"/Rresults/overall/k3/compute_crps-CV.A/05_TWS_member1_k3_mean3day_1950-01-01_2011-12-31.Rdata"))
  crps_ana <- crps[,"pos ecdf"]
  
  load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
  crps_clim <- score$crps[,"pos ecdf"]
  
  rm(crps,score)
  
  # Quantiles de pluie
  precip <- get.precip(nbdays,start,end)
  pos <- sort(precip,index.return=T)$ix[-c(1:length(precip[precip==0]))] # ordre croissant des pluies positives
  quant <- unname(round(quantile(1:length(pos),probs=seq(0,1,0.1)),0)) # position des quantiles
  
  # Graphique ind VS climato
  png(file=paste0(get.dirstr(k,rean),"plot.crps",get.CVstr(CV),"/plot_ind_clim_",descriptors[1],"_",descriptors[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"),width = 14,height = 7,units="in",res = 200)
  par(mfrow=c(2,5))
  for(i in 1:(length(quant)-1)){
    indic <- pos[quant[i]:quant[i+1]]
    xlim  <- range(crps_clim[indic],crps_ind[indic],na.rm=T)
    par(pty="s")
    plot(crps_clim[indic],crps_ind[indic],col=getcol(i,range=c(1,10),rev=T),xlim=xlim,ylim=xlim,
         xlab="CRPS climato",ylab="CRPS ind",
         main=paste0("Quantile ",(i-1)*10,"% - ",i*10,"% (",round(min(precip[indic],na.rm=T),1),"-",round(max(precip[indic],na.rm=T),1)," mm/day)"))
    abline(0,1,lwd=2)
  }
  dev.off()
  
  # Graphique ana VS climato
  png(file=paste0(get.dirstr(k,rean),"plot.crps",get.CVstr(CV),"/plot_ana_clim_",descriptors[1],"_",descriptors[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"),width = 14,height = 7,units="in",res = 200)
  par(mfrow=c(2,5))
  for(i in 1:(length(quant)-1)){
    indic <- pos[quant[i]:quant[i+1]]
    xlim  <- range(crps_clim[indic],crps_ana[indic],na.rm=T)
    par(pty="s")
    plot(crps_clim[indic],crps_ana[indic],col=getcol(i,range=c(1,10),rev=T),xlim=xlim,ylim=xlim,
         xlab="CRPS climato",ylab="CRPS ana",
         main=paste0("Quantile ",(i-1)*10,"% - ",i*10,"% (",round(min(precip[indic],na.rm=T),1),"-",round(max(precip[indic],na.rm=T),1)," mm/day)"))
    abline(0,1,lwd=2)
  }
  dev.off()
  
  # Boxplot 0-90%
  df <- data.frame(Climatologie=crps_clim[pos[1:quant[10]]],
                   Analogues=crps_ana[pos[1:quant[10]]],
                   Indicateurs=crps_ind[pos[1:quant[10]]])
                   
  df$Quantile <- as.character(c(rep(1:9,each=nrow(df)/9),rep(9,4)))
  df <- gather(data = df,key = "Methode",value = "CRPS",1:3)
  df$Methode <- factor(df$Methode, levels = c("Climatologie","Analogues","Indicateurs"))
  
  ggplot(df, aes(x=Quantile, y=CRPS, fill=Methode)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,3,1.5,1.5),"cm"),legend.key.size = unit(2,"cm"))+
    geom_boxplot(outlier.shape = NA)+
    ylim(0,5)+
    scale_x_discrete(breaks=as.character(1:9),labels=paste0(seq(10,90,10),"%"))+
    scale_fill_manual(values=c("chartreuse3","royalblue","red"))+
    geom_hline(yintercept = 0,linetype="dashed")+
    ggtitle("CRPS sur la distribution des differents quantiles de pluie")+
    theme(plot.title = element_text(hjust = 0.5))
    
  ggsave(paste0(get.dirstr(k,rean),"plot.crps",get.CVstr(CV),"/boxplot_0_90_",descriptors[1],"_",descriptors[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"), width = 16, height = 9, dpi = 100)
  graphics.off()
  
  # Boxplot 90-100%
  df <- data.frame(Climatologie=crps_clim[pos[quant[10]:quant[11]]],
                   Analogues=crps_ana[pos[quant[10]:quant[11]]],
                   Indicateurs=crps_ind[pos[quant[10]:quant[11]]],
                   Quantile=as.character(10))
  
  df <- gather(data = df,key = "Methode",value = "CRPS",1:3)
  df$Methode <- factor(df$Methode, levels = c("Climatologie","Analogues","Indicateurs"))
  
  ggplot(df, aes(x=Quantile, y=CRPS, fill=Methode)) + 
    theme_bw()+
    theme(plot.margin = unit(c(2,3,2,2),"cm"),legend.key.size = unit(2,"cm"))+
    geom_boxplot()+
    ylim(0,30)+
    scale_x_discrete(breaks=as.character(1),labels=paste0(10,"%"))+
    scale_fill_manual(values=c("chartreuse3","royalblue","red"))+
    geom_hline(yintercept = 0,linetype="dashed")+
    ggtitle("CRPS sur la distribution des differents quantiles de pluie")+
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0(get.dirstr(k,rean),"plot.crps",get.CVstr(CV),"/boxplot_90_100_",descriptors[1],"_",descriptors[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"), width = 8, height = 5, dpi = 100)
  graphics.off()
  
}

# plot la comparaison des CRPS analogie classique et climato pour differents quantiles de pluie
plot.crps.ana<-function(rean,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",str="05",CV=TRUE){
  
  # Import des scores ana & climato
  load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/",str,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  crps_ana <- crps[,"pos ecdf"]
  
  load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata"))
  crps_clim <- score$crps[,"pos ecdf"]
  
  rm(crps,score)
  
  # Quantiles de pluie
  precip <- get.precip(nbdays,start,end)
  pos <- sort(precip,index.return=T)$ix[-c(1:length(precip[precip==0]))] # ordre croissant des pluies positives
  quant <- unname(round(quantile(1:length(pos),probs=seq(0,1,0.1)),0)) # position des quantiles
  
  # Boxplot
  df <- data.frame(Climatology=crps_clim[pos], Analog=crps_ana[pos])
  
  df$Quantile <- c(rep(1:10,each=nrow(df)/10),rep(10,4))
  df$Quantile <- paste0((df$Quantile-1)*10,"% - ",df$Quantile*10,"%")
  rain.qua <- quantile(precip[precip>0],probs=seq(0,1,0.1))
  rain <- NULL
  for(i in 1:10) rain[i] <- paste0("[",round(rain.qua[i],2),"-",round(rain.qua[i+1],2),"]")
  rain <- c(rep(rain,each=nrow(df)/10),rep(rain[10],4))
  df$Quantile <- paste0(df$Quantile,"\n ",rain)
  
  df <- gather(data = df,key = "Method",value = "CRPS",1:2)
  df$Method <- factor(df$Method, levels = c("Climatology","Analog"))
  
  ggplot(df, aes(x=Quantile, y=CRPS, fill=Method)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,3,1.5,1.5),"cm"),legend.key.size = unit(2,"cm"))+
    geom_boxplot(outlier.shape = NA)+
    xlab("Quantile (mm/d)")+
    theme(axis.text.x = element_text(color="black", face="bold",size=9,angle=45,vjust=0.5))+
    ylim(0,5)+
    scale_fill_manual(values=c("royalblue","red"))+
    geom_hline(yintercept = 0,linetype="dashed")+
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(paste0(get.dirstr(k,rean),"plot.crps.ana",get.CVstr(CV),"/boxplot_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"), width = 10, height = 5, dpi = 100)
  graphics.off()
}

# plot la distribution modelisee et observee de certaines sequences de pluie
plot.distrib<-function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE){
  
  # Import des precip et des indicateurs
  precip<-get.precip(nbdays,start,end)
  
  descr1<-get.descriptor(descriptors[1],k,dist,nbdays,start,end,standardize=TRUE,rean)
  descr2<-get.descriptor(descriptors[2],k,dist,nbdays,start,end,standardize=TRUE,rean)
  descr <- cbind(descr1,descr2)
  
  load(file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei<-nei[["05"]]
  
  # Sequences qu'on veut ploter
  ind <- sample(which(precip >0.1 & precip<1),size = 9)
  dates <- paste0(getdates(start,end)[ind]," - ",getdates(start,end)[ind+2])
  
  png(file=paste0(get.dirstr(k,rean),"plot.distrib",get.CVstr(CV),"/plot_",descriptors[1],"_",descriptors[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"),width=17,height=13,units="in",res=200)
  par(mfrow=c(3,3))
  
  for(i in ind){
    # climato
    pi<-precip[setdiff(1:length(precip),(i-nbdays+1):(i+nbdays-1))]
    plot(ecdf(pi),col="green",xlab="Precipitations (mm/day)",main=paste0(dates[i == ind]," (",round(precip[i],1)," mm/day)"),lwd=3)
    
    # obs
    segments(0,0,precip[i],0,lwd = 2)
    segments(precip[i],0,precip[i],1,lwd = 2)
    segments(precip[i],1,100,1,lwd = 2)
    
    # analogie classique
    pi<-precip[setdiff(nei[[i]],(i-nbdays+1):(i+nbdays-1))]
    lines(ecdf(pi),col="royalblue",cex=0.8)
    
    # indicateurs
    tmp<-get.closest(i,descr,precip,CV,nbdays,radtype)
    pi<-precip[tmp$idx]
    lines(ecdf(pi),col="red",cex=0.8)

    legend("bottomright",inset=.05,legend = c("Obs","Climato","Analog","Indicateurs"),cex=1.2,
           lty=1,lwd = c(2,2,1,1),pch = c(NA,NA,19,19),bty="n",col=c("black","green","royalblue","red"))
    
  }
  graphics.off()
 
  # Analogie seule
  png(file=paste0(get.dirstr(k,rean),"plot.distrib",get.CVstr(CV),"/plot_ana_seule_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".png"),width=450,height=350,units="px")
  pos <- 7970
  date <- paste0(getdates(start,end)[pos]," - ",getdates(start,end)[pos+2])
  pi<-precip[setdiff(nei[[pos]],(pos-nbdays+1):(pos+nbdays-1))]
  plot(ecdf(pi),col="red",xlab="Precipitations (mm/day)",main=paste0(date," (",round(precip[pos],1)," mm/day)"),lwd=3)
  segments(0,0,precip[pos],0,lwd = 2)
  segments(precip[pos],0,precip[pos],1,lwd = 2)
  segments(precip[pos],1,100,1,lwd = 2)
  legend("bottomright",inset=.05,legend = c("Obs","Analog"),cex=1.2,lty=1,lwd = c(2,1),pch = c(NA,19),bty="n",col=c("black","red"))
  graphics.off()
  
}

# Trace la repartition des parametres empiriques dans le plan des indicateurs (un couple d'indicateur par png)
plot.empir<-function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,obs=FALSE,threeday=c(F,F)){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import des precip, du fit.empir, des indicateurs
  precip <- get.precip(nbdays,start,end)
  
  if(TRUE %in% threeday) {comp <- paste0("_threeday",which(threeday==T))
  } else {comp <- ""}
  load(file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,comp,".Rdata"))
  param<-param[,colnames(param)!="kurtosis"]
  
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=FALSE,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=FALSE,rean[2],threeday[2])
  
  # Nom propre des indicateurs
  cond <- substr(descriptors,nchar(descriptors)-1,nchar(descriptors)) == substr(radtype,nchar(radtype)-1,nchar(radtype))
  descriptors[cond] <- substr(descriptors[cond],1,nchar(descriptors[cond])-2)
  namdescr <- nam2str(descriptors, cloud = TRUE)
  
  # Creation du plot
  png(file=paste0(path1,"plot.empir",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,".png"),width=640,height=450,units="px",res=72) # manip substr pour enlever le std}
  par(mfrow=c(2,3))
  gamme <- list(6)
  gamme[[1]] <- c(0,5380)
  gamme[[2]] <- c(0,0.7)
  gamme[[3]] <- c(0.02,16.32)
  
  for (i in 1:ncol(param)){
    param.i<-param[,i]
    plot(descr1,descr2,
         col=getcol(param.i,range = ifelse(rep(i %in% c(1:3),2),gamme[[i]],range(param.i,na.rm=T))),
         xlab=namdescr[1],
         ylab=namdescr[2],
         ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
    if(i %in% c(1:3)) {addscale(vec = c(param.i,gamme[[i]]))
    } else {addscale(vec = param.i)}
    text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.9,
         y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
         paste0(round(min(param.i,na.rm=T),2),"-",round(max(param.i,na.rm=T),2)))
    #ind.extr <- get.ind.extr(nbre = 3,nbdays = 3)
    #points(descr1[ind.extr],descr2[ind.extr],pch=19)
    title(colnames(param)[i])
  }
  dev.off()
  
  if(obs){
    # Indices
    ind <- seq(1,length(precip),250)[seq(1,90,length.out = 12)]
    ord <- sort(precip,decreasing = T,index.return=T)
    
    # plot
    png(file=paste0(path1,"plot.empir",get.CVstr(CV),"/plot_obs_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,".png"),width=14,height=11,units="in",res=72) # manip substr pour enlever le std
    par(mfrow=c(3,4))
    
    for(i in ind){
      param.i<-param[,1]
      plot(descr1,descr2,
           col=getcol(param.i,range = gamme[[1]]),
           xlab=namdescr[1],
           ylab=namdescr[2],
           ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
      addscale(vec = c(param.i,gamme[[1]]))
      text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.9,
           y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
           paste0(round(min(param.i,na.rm=T),2),"-",round(max(param.i,na.rm=T),2)))
      points(descr1[ord$ix[i:(i+249)]],descr2[ord$ix[i:(i+249)]],pch=19)
      title(paste0("Rank ",i,"-",i+249,", ",round(ord$x[i+249],1),"-",round(ord$x[i],1)," mm/day"))
    }
    dev.off()
  }
  
}

# plot.empir pour presentations
plot.empir.clean<-function(sel=c("p0","mean"),rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,obs=FALSE,threeday=c(F,F)){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import des precip, du fit.empir, des indicateurs
  precip <- get.precip(nbdays,start,end)
  
  if(TRUE %in% threeday) {comp <- paste0("_threeday",which(threeday==T))
  } else {comp <- ""}
  load(file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,comp,".Rdata"))
  param<-param[,colnames(param)!="kurtosis"]
  
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=FALSE,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=FALSE,rean[2],threeday[2])
  
  # Nom propre des indicateurs
  cond <- substr(descriptors,nchar(descriptors)-1,nchar(descriptors)) == substr(radtype,nchar(radtype)-1,nchar(radtype))
  descriptors[cond] <- substr(descriptors[cond],1,nchar(descriptors[cond])-2)
  namdescr <- nam2str(descriptors, cloud = TRUE)
  
  # Creation du plot
  pdf(file=paste0(path1,"plot.empir.clean",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,".pdf"),width=8,height=3) # manip substr pour enlever le std
  par(mfrow=c(1,3))
  
  gamme <- list()
  if(all(sel=="mean")){
    gamme[[1]] <- c(0,5380)
    gamme[[2]] <- c(0.02,16.32)
  }else{
    gamme[[1]] <- c(0,5380)
    gamme[[2]] <- c(0,0.7)
    gamme[[3]] <- c(0.02,16.32)
  }
    
  param <- param[,c("nbnei",sel)]
  
  for (i in 1:ncol(param)){
    param.i<-param[,i]
    plot(descr1,descr2,
         col=getcol(param.i,range = ifelse(rep(i %in% c(1:3),2),gamme[[i]],range(param.i,na.rm=T))),
         #xlab=namdescr[1],
         #ylab=namdescr[2],
         xlab="Célérité",
         ylab="Singularité",
         ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
    if(i %in% c(1:3)) {addscale(vec = c(param.i,gamme[[i]]))
    } else {addscale(vec = param.i)}
    text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.8,
         y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
         paste0(round(min(param.i,na.rm=T),2),"-",round(max(param.i,na.rm=T),2)))
    #ind.extr <- get.ind.extr(nbre = 3,nbdays = 3)
    #points(descr1[ind.extr],descr2[ind.extr],pch=19)
    title(colnames(param)[i])
  }
  dev.off()
  
  if(obs){
    # Indices
    ind <- seq(1,length(precip),250)[seq(1,90,length.out = 12)]
    ord <- sort(precip,decreasing = T,index.return=T)
    
    # plot
    png(file=paste0(path1,"plot.empir",get.CVstr(CV),"/plot_obs_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,".png"),width=14,height=11,units="in",res=72) # manip substr pour enlever le std
    par(mfrow=c(3,4))
    
    for(i in ind){
      param.i<-param[,1]
      plot(descr1,descr2,
           col=getcol(param.i,range = gamme[[1]]),
           xlab=namdescr[1],
           ylab=namdescr[2],
           ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
      addscale(vec = c(param.i,gamme[[1]]))
      text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.9,
           y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
           paste0(round(min(param.i,na.rm=T),2),"-",round(max(param.i,na.rm=T),2)))
      points(descr1[ord$ix[i:(i+249)]],descr2[ord$ix[i:(i+249)]],pch=19)
      title(paste0("Rank ",i,"-",i+249,", ",round(ord$x[i+249],1),"-",round(ord$x[i],1)," mm/day"))
    }
    dev.off()
  }
  
}

# plot.empir pour presentations, avec obs
plot.empir.clean.obs<-function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F)){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import des precip, du fit.empir, des indicateurs
  precip <- get.precip(nbdays,start,end)
  
  if(TRUE %in% threeday) {comp <- paste0("_threeday",which(threeday==T))
  } else {comp <- ""}
  load(file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,comp,".Rdata"))
  param<-param[,colnames(param)!="kurtosis"]
  
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=FALSE,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=FALSE,rean[2],threeday[2])
  
  # Nom propre des indicateurs
  cond <- substr(descriptors,nchar(descriptors)-1,nchar(descriptors)) == substr(radtype,nchar(radtype)-1,nchar(radtype))
  descriptors[cond] <- substr(descriptors[cond],1,nchar(descriptors[cond])-2)
  namdescr <- nam2str(descriptors, cloud = TRUE)
  
  # Creation du plot
  pdf(file=paste0(path1,"plot.empir.clean.obs",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,"_obs.pdf"),width=8,height=3) # manip substr pour enlever le std
  par(mfrow=c(1,3))
  
  gamme <- c(0,5380)
  param <- param[,c("nbnei")]
  
  plot(descr1,descr2,
       col=getcol(param,range = gamme),
       #xlab=namdescr[1],
       #ylab=namdescr[2],
       xlab="Célérité",
       ylab="Singularité",
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
  addscale(vec = c(param,gamme))
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.8,
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
       paste0(round(min(param,na.rm=T),2),"-",round(max(param,na.rm=T),2)))
  title(names(param))
  
  plot(descr1,descr2,
       col=getcol(param,range = gamme),
       #xlab=namdescr[1],
       #ylab=namdescr[2],
       xlab="Célérité",
       ylab="Singularité",
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
  addscale(vec = c(param,gamme))
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.8,
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
       paste0(round(min(param,na.rm=T),2),"-",round(max(param,na.rm=T),2)))
  ind.extr <- get.ind.extr(nbre = 250,nbdays = nbdays)
  points(descr1[ind.extr],descr2[ind.extr],pch=19,cex=0.5)
  title(names(param))
  
graphics.off()

}

# plot.empir seulement pour la moyenne, avec evenements extremes et rectangles de selections
plot.empir.mean<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,rean,extr = TRUE,rect = TRUE,ref = "1950-01-01"){
  
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  print(paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata"))
  load(file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata")) # importation des parametres de la loi empirique
  param<-param[,colnames(param)=="mean"]
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE,rean)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE,rean)
  
  print(paste0(get.dirstr(k,rean),"plot.empir.mean",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".pdf"))
  pdf(file=paste0(get.dirstr(k,rean),"plot.empir.mean",get.CVstr(CV),"/plot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(FALSE),"_",radtype,".pdf"),width=7,height=7)
  
  par(pty="s")
  plot(descr1,descr2,col=getcol(param),xlab=descriptor1,ylab=descriptor2)
  addscale(param)
  title("mean")
  
  if(extr){
    ind.extr <- get.ind.extr(nbre = 3,nbdays = 3)
    points(descr1[ind.extr],descr2[ind.extr],pch=19,cex=1.5)
  }
  
  if(rect){
    # Definition des boites auxquelles on s'interesse
    boxes <- list(4)
    boxes[[1]] <- c(0.14,0.88,0.18,0.895)
    boxes[[2]] <- c(0.17,0.92,0.21,0.935)
    boxes[[3]] <- c(0.24,0.91,0.28,0.925)
    boxes[[4]] <- c(0.30,0.94,0.34,0.955)
    
    
    for(i in 1:length(boxes)){
      ind <- boxes[[i]]
      # Trace des boites
      rect(xleft = ind[1],
           ybottom = ind[2],
           xright = ind[3],
           ytop = ind[4],
           lwd = 2)
      
      # Indices des journees concernees
      pos1 <- which(descr1>ind[1] & descr1<ind[3])
      pos2 <- which(descr2>ind[2] & descr2<ind[4])
      
      # Position des sequences de 3 jours qui sont dans le rectangle
      pos  <- pos1[pos1 %in% pos2] 
      
      # Position de tous les jours appartenant a ces sequences
      if(nbdays!=1){
        tmp <- list(nbdays)
        for(j in 1:nbdays) tmp[[j]] <- pos+j-1
        pos <- unique(sort(unlist(tmp)))
      } 
      
      # Nombre de journees dans le rectangles
      print(paste0("Rectangle ",i,": ",length(pos)," journees considerees"))
      
      # Passage en indices selon la reference souhaitee
      delta <- get.delta(ref = ref)
      pos <- pos + delta
      write(x = pos, file = paste0(get.dirstr(k,rean),"plot.empir.mean",get.CVstr(CV),"/Indices_rect_",i,".txt"), ncolumns = length(pos),sep = ",")
    }
  }
  
  dev.off()
}

# Trace la repartition des points dans le plan et colorie la selection par indicateurs
plot.empir.sel<-function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import des precip, du fit.empir, des indicateurs
  precip <- get.precip(nbdays,start,end)
  
  load(file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,".Rdata"))
  param<-param[,colnames(param)!="kurtosis"]
  
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=T,rean[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=T,rean[2])
  
  # Nom propre des indicateurs
  cond <- substr(descriptors,nchar(descriptors)-1,nchar(descriptors)) == substr(radtype,nchar(radtype)-1,nchar(radtype))
  descriptors[cond] <- substr(descriptors[cond],1,nchar(descriptors[cond])-2)
  namdescr <- nam2str(descriptors, cloud = TRUE)
  
  # Creation du plot
  png(file=paste0(path1,"plot.empir.sel",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),".png"),width=520,height=400,units="px",res=72) # manip substr pour enlever le std
  plot(descr1,descr2,
       asp=1,
       col="grey",
       xlab=namdescr[1],
       ylab=namdescr[2])#,
       #ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+max(descr2,na.rm=T)-min(descr2,na.rm=T))
    
  ind <- get.closest(22643,cbind(descr1,descr2),precip,CV,nbdays,radtype)
  points(descr1[ind$idx],descr2[ind$idx],col="black",pch=21,bg="red")
  points(descr1[22643],descr2[22643],pch=24,cex=1.1,col="white",bg="black")
  
  
  dev.off()
  
}

# plot le papillon de Lorenz en 3d colorie par indicateur
plot.lorenz <- function(){
  
  # papillon
  parameters <- c(a = -8/3,b = -10,c = 28)
  state <- c(X = 1,Y = 1,Z = 1)
  
  Lorenz<-function(t, state, parameters) {
    with(as.list(c(state, parameters)),{
      # rate of change
      dX <- a*X + Y*Z
      dY <- b * (Y-Z)
      dZ <- -X*Y + c*Y - Z
      
      # return the rate of change
      list(c(dX, dY, dZ))
    }) # end with(as.list ...
  }
  
  times <- seq(0, 100, by = 0.01)
  out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
  
  # calcul de sing et rsing
  sing <- rsing <- vector(length = nrow(out))
  
  for(i in 1:nrow(out)){
    if (i %% 50 == 0) print(i)
    nei <-nn2(data = out[-i,-1],query = t(out[i,-1]),k = 0.005*nrow(out))$nn.dists
    sing[i] <- mean(nei)
    rsing[i] <- sing[i]/max(nei)
  }
  
  pdf("2_Travail/Rresults/plot.lorenz/plot_sing_rsing.pdf",width = 7.5,height = 7.5)
  lines3D(x = out[,2],y = out[,3],z=out[,4],colvar=sing,bty="g",theta=-60,phi=40,main="sing")
  lines3D(x = out[,2],y = out[,3],z=out[,4],colvar=rsing,bty="g",theta=-60,phi=40,main="rsing")
  
  descr <- cbind(sing,rsing)
  rad <- mean(c(sd(descr[,1],na.rm=TRUE),sd(descr[,2],na.rm=TRUE)))/2
  nb <- vector(length = nrow(descr))
  
  for(i in 1:nrow(descr)){
    if (i %% 50 == 0) print(i)
    count <- nn2(data = descr,query = t(descr[i,]),searchtype = "radius",radius =  rad,k=nrow(descr)) # nombre de voisins dans le rayon
    nb[i] <- rowSums(count$nn.idx>0)-1
  }
  plot(sing,rsing,ylim=c(min(rsing),min(rsing)+(max(rsing)-min(rsing))*1.2),col=getcol(nb))
  addscale(nb)
  graphics.off()
}

# plot p0 en fonction du WP et de la saison
plot.p0.wp <- function(){
  
  # Import
  load("2_Travail/Rresults/save.wp.ind/WP_ind.Rdata")
  precip <- get.precip(1)
  
  # Calcul
  p0 <- lapply(wp.ind,function(x) {p=precip[x];res=table(p==0)[2]/length(p);return(res)})
  p0 <- unlist(p0,use.names = F)
  dim(p0) <- c(8,2)
  
  # Graphique
  png(filename = "2_Travail/Rresults/plot.p0.wp/plot_p0_wp.png",width = 580,height = 440,units = "px")
  plot(p0[,1],pch=19,ylim=c(0,0.8),col="red",xlab="Weather Pattern",ylab="p0")
  lines(p0[,1],col="red")
  points(p0[,2],pch=19)
  lines(p0[,2])
  grid()
  legend("topleft",c("Season-at-risk (Sept-Fev)","Low-risk-season (Mar-Aug)"),col=c("red","black"),pch=19,lty=1,bty="n")
  graphics.off()
  
}

# Trace le boxplot des quantiles d'un indicateur, pour 500 & 1000, pour les 4 distances
plot.quant.descr <- function(descr,nbdays,start="1950-01-01",end="2011-12-31",rean){
  
  # Import de l'indicateurs pour k1, k2, et les 4 distances
  k <- c(1,2)
  dist <- c("TWS","RMSE","RMSE_I","RMSE_II")
  mat <- NULL
  
  for(i in k){
    for(j in dist){
      des <- get.descriptor(descriptor = descr,k = i,dist = j,nbdays = nbdays,
                              start = start,end = end,standardize = F,rean = rean)
      mat <- cbind(mat,des)
      colnames(mat)[ncol(mat)] <- paste0(j,"_",i)
    }
  }
  
  # Traitement
  ind.year <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  ind.mon <- get.ind.max(type = "month",nbdays = nbdays,start = start,end = end)
  ind <- list(ind.year,ind.mon)
  
  qua <- apply(mat,2,function(v) ecdf(v)(v)*100)
  l <- vector("list",2)
  
  for(i in 1:length(ind)){
  tmp <- qua[ind[[i]],]
  tmp <- pivot_longer(data = as.data.frame(tmp),1:8,names_to = "dist",values_to = "quantile")
  tmp$geo <- substr(tmp$dist,nchar(tmp$dist),nchar(tmp$dist))
  tmp$dist <- substr(tmp$dist,1,nchar(tmp$dist)-2)
  tmp$dist <- factor(tmp$dist,levels = dist)
  l[[i]] <- tmp
  }
  
  # Boxplot
  ggplot(l[[1]], aes(x=dist, y=quantile, fill=geo)) + 
    theme_bw()+
    theme(plot.margin = unit(c(2,3,2,2),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))+
    stat_boxplot(geom = "errorbar", width = 0.2,col="darkblue",position = ) +
    geom_boxplot(outlier.shape = NA)+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")
  
  ggplot(l[[2]], aes(x=dist, y=quantile, fill=geo)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,3,1.5,1.5),"cm"),legend.key.size = unit(2,"cm"))+
    geom_boxplot(outlier.shape = NA)+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")
  
}

# Trace les hyétorgammes des 62 plus grosses séquences de pluie
plot.rain <- function(nbdays=3,start="1950-01-01",end="2011-12-31"){
  
  precip <- get.precip(nbdays=1,start,end)
  dates <- seq(as.Date(start),as.Date(end),by="day")
  
  ind <- get.ind.extr(nbre = 62, ref = start, nbdays, start, end, nei = TRUE)
  pdf(file = paste0("2_Travail/Rresults/plot.rain/Hietograph","_",nbdays,"day.pdf"))
  layout(matrix(1:20,4,5,byrow=TRUE))
  for(i in 1:length(ind)){
    barplot(precip[ind[i]:(ind[i]+nbdays-1)],space=0,col = "royalblue",main = dates[ind[i]],names.arg = c("D1","D2","D3"),ylab="Precip (mm)",cex.names = 0.8)
  }
  graphics.off()
  
  # Frequence d'occurence des sequences de pluie fortes dans l'annee
  png(file = "2_Travail/Rresults/plot.rain/Histogram_month.png",width = 419,height = 322,units = "px",res = 72)
  tmp <- as.numeric(substr(dates[ind],6,7))
  hist(tmp,0:12,axes = FALSE,xlab="Month",col = "royalblue",border = "white",main="62 largest 3-day precip 1950-2011")
  lines(c(0,12),c(0,0))
  axis(2)
  axis(1,at = 0.5:11.5,labels = 1:12,tick = FALSE,padj = -1)
  graphics.off()
  
  # Hyetogramme des cumuls de pluie sur l'annee
  png(file = "2_Travail/Rresults/plot.rain/Hyetograph_month.png",width = 460,height = 322,units = "px",res = 72)
  precip <- get.precip(1)
  dates <- getdates(start,end)
  month <- substr(dates,6,7)
  barplot(aggregate(precip,by=list(month),FUN=function(x) sum(x)/62)[,2],ylim=c(0,130),cex.names = 0.8,
          names.arg = month.abb,col="royalblue",xlab="Month",ylab="precip (mm)",main="Monthly precip accumulation")
  graphics.off()          

}

# plot le ratio de pluie entre les sous BV sur les sequences
plot.ratio.BV <- function(nbdays=3,start="1950-01-01",end="2011-12-31",rain=FALSE,flood=FALSE,drac=FALSE,isere=FALSE){
  
  precip_BV    <- get.precip(3,start,end)
  precip_Isere <- get.precip(nbdays,start,end,bv="Isere")
  precip_Drac  <- get.precip(nbdays,start,end,bv="Drac")
  
  ratio <- precip_Isere/precip_Drac
  ratio[ratio==Inf] <- NA
  tmp <- NULL
  supp0 <- "BV"
  supp <- ""
  
  if(nbdays==1){
    n<-length(precip_BV)
    for (i in 1:3) tmp<-cbind(tmp,ratio[(i:(n+i-1))]) 
  }
  
  # Evenement fort sur un sous BV (Drac ou Isere en abscisses)
  if(drac) {precip_BV <- get.precip(3,start,end,bv="Drac");supp0 <- "Drac"}
  if(isere) {precip_BV <- get.precip(3,start,end,bv="Isere");supp0 <- "Isere"}
  
  # Plus fortes pluies
  if(rain){
    pos <- sort(precip_BV,decreasing = T,index.return = T)$ix[1:62]
    
    precip_BV <- precip_BV[pos]
    ratio <- ratio[pos]
    tmp <- tmp[pos,]
    supp <- "rain_"
  }
  
  # Plus fortes crues
  if(flood){
    load("2_Travail/Rresults/image.cumul/crues_1969_2011_Isere_StGervais.Rdata")
    crues <- as.data.frame(crues[crues[,4]>0.75,]) # on selectionne les events quadriennaux
    event <- as.character(crues[,1]-3) # on prend les pluies de j-3,j-2,j-1 de la crue (car 8hj a 7hj+1)
    pos   <- match(event,getdates())
    
    precip_BV <- precip_BV[pos]
    ratio <- ratio[pos]
    tmp <- tmp[pos,]
    supp <- "flood_"
  }
  
  
  png(file=paste0("2_Travail/Rresults/plot.ratio.BV/plotratio_",supp0,"_",supp,nbdays,"days_",start,"_",end,".png"),width=520,heigh=480,units="px",res=72)
  plot(range(precip_BV),range(ratio,na.rm=T),log="y",type="n",xlab=paste0(supp0," ","3-day precip (mm/day)"),ylab="precip Isere / precip Drac")
  if(nbdays==1) {
    for (i in 1:ncol(tmp)) points(precip_BV,tmp[,i],col=i,pch=19,cex=0.9)
    legend("topright",c("J1","J2","J3"),col=c(1,2,3),pch=19,pt.cex = 0.9,bty="n")
  }else{points(precip_BV,ratio,pch=19,cex=0.9)}
  abline(h=1,col="red",lwd=2)
  if(nbdays==1){
    pt <- apply(tmp,2,mean,na.rm=T)
    for(i in 1:length(pt)) abline(h=pt[i],lty=2,col=i)
  }
  dev.off()
}

# plot le ratio de pluvios concernes pour les differents quantiles de precip sur le BV
plot.ratio.pluvios <- function(nbdays,seuil=0,start="1950-01-01",end="2011-12-31"){
  
  # Import des donnees station et BV
  precip <- get.precip(nbdays,start,end)
  ind <- which(precip>seuil)
  precip <- precip[ind]
  
  load("2_Travail/Data/Precip/data_precip.Rdata")
  data_sta <- as.data.frame(data$precip/10)
  rownames(data_sta) <- as.Date(rownames(data_sta),"%Y%m%d")
  data_sta <- subset(data_sta,rownames(data_sta)>=start & rownames(data_sta)<=end)
  
  # Traitement des donnees
  if(nbdays>1) data_sta <- rollapply(data_sta,width=nbdays,FUN=mean)
  ratio <- apply(data_sta,1,function(v) sum(v>seuil,na.rm=T)/sum(!is.na(v))*100)
  ratio <- ratio[ind]
  qua <- ecdf(precip)(precip)*10
  qua[which.max(qua)] <- max(qua)-0.001 # -0.0001 car ecdf donne quantile 100% pour valeur max
  qua <- as.integer(qua)+1
  df <- as.data.frame(cbind(qua,ratio))
  
  xlabel <- NULL
  for(i in 1:10){
    q1 <- round(quantile(precip,probs=(i-1)*0.1),1)
    q2 <- round(quantile(precip,probs=i*0.1),1)
    xlabel[i] <- paste0(i*10-10,"%-",i*10,"%\n",q1,"-",q2)
  }
  
  # Graphique
  ggplot(df, aes(x=qua, y=ratio, group=qua)) + 
    theme_bw()+
    theme(plot.margin = unit(c(2,3,2,2),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10))+
    stat_boxplot(geom = "errorbar", width = 0.2,col="darkblue") +
    geom_boxplot(outlier.shape = NA,fill="cornflowerblue",color="darkblue")+
    xlab(paste0("Isere catchment basin ",nbdays,"-day positive precipitation (mm/d)"))+
    ylab(paste0("% of stations with positive \n",nbdays,"-day precipitation"))+
    scale_x_discrete(limits = 1:10,labels=xlabel)
  
  ggsave(filename = paste0("2_Travail/Rresults/plot.ratio.pluvios/ratio_",nbdays,"day_seuil",seuil,"_",start,"_",end,".png"),
         width = 28,height = 15,units = "cm")
  graphics.off()
}
  
# Trace la fonction de repartition des scores
plot.score <- function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean,rsing=FALSE,sing=FALSE){
  
  # Importation et mise en forme des scores
  dist.list <- score.to.mat(k,dist,nbdays,start,end,rean)
  dist1 <- dist
  
  # Pluies qui nous interessent
  precip <- get.precip(nbdays,start,end)
  soso <- sort(precip,index.return=TRUE,decreasing=TRUE)
  ind <- vector("list",3)
  ind[[1]] <- 1:(length(precip)-1)
  ind[[2]] <- which(precip == 0)
  ind[[3]] <- soso$ix[1:1247]
  
  leg <- ifelse(nbdays>1,"Sequences","Days")
  comp.leg <- c(nrow(dist.list),length(precip[precip==0]),1247)
  names(ind) <- c(paste0("All ",leg),paste0("Dry ",leg),paste0("Heavy Rain ",leg))
  
  # Mise en forme et calculs
  score <- apply(dist.list,2,sort)
  score <- t(score)
  score <- score[,-1]
  
  if(rsing){
    score <- t(apply(score,1,compute.rsing))
    dist1 <- "Relative singularity"
  }  
  if(sing){
    score <- t(apply(score,1,compute.rsing,sing=TRUE))
    dist1 <- "Singularity" 
  }  
  
  legy <- paste0(dist1," ",nbdays,ifelse(nbdays==1," day"," days"))
  legm <- paste0(dist1," Repartition - ",nbdays,ifelse(nbdays==1," day"," days"))
  
  moy <- med <- q5 <- q95 <- matrix(NA,3,ncol(score))
  
  for(i in 1:length(ind)){
    score.ind <- score[ind[[i]],]
    moy[i,] <- colMeans(score.ind)
    med[i,] <- apply(score.ind,2,median)
    q5[i,]  <- apply(score.ind,2,quantile,probs=0.05)
    q95[i,] <- apply(score.ind,2,quantile,probs=0.95)
  }
  
  # Graphique complet
  png(file=paste0(get.dirstr(k,rean),"plot.score/plot_repartition_",ifelse(rsing,"rsing_",""),ifelse(sing,"sing_",""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width=930,height=520,units="px",res=72)
  plot(med[1,],xlab=leg,ylab=legy,ylim = c(min(q5),max(q95)),type="l",main=legm,lwd=3)
  lines(med[2,],col = "red",lwd=2)
  lines(q5[2,],col = "red",lty=2)
  lines(q95[2,],col = "red",lty=2)
  lines(med[3,],col = "royalblue",lwd=2)
  lines(q5[3,],col = "royalblue",lty=2)
  lines(q95[3,],col = "royalblue",lty=2)
  grid()
  legend(ifelse(rsing,"bottomleft","topleft"),inset=.02,c(paste0("Q50% ",names(ind)," (",comp.leg,")"),"Q05%,Q95%"),col=c("black","red","royalblue","black"),cex=1.2,lty=c(1,1,1,2),lwd = c(3,2,2,1),bty="n")
  graphics.off()
  
  # Chevelu
  png(file=paste0(get.dirstr(k,rean),"plot.score/plot_chevelu_",ifelse(rsing,"rsing_",""),ifelse(sing,"sing_",""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width=930,height=520,units="px",res=72)
  plot(moy[1,],xlab=leg,ylab=legy,ylim = c(min(score[ind[[3]],]),max(score[ind[[3]],])),type="n",main=legm,lwd=3)
  for(i in 1:62) lines(score[ind[[3]][i],],col = "royalblue")
  lines(moy[1,],lwd=3)
  grid()
  legend(ifelse(rsing,"bottomleft","topleft"),inset=.02,c(paste0(names(ind)[1]," Mean (",comp.leg[1],")"),paste0(names(ind)[3]," (",comp.leg[3],")")),col=c("black","royalblue"),cex=1.2,lty=1,lwd = c(3,1),bty="n")
  graphics.off()
  
  # Graphique faibles quantiles
  #png(file=paste0(get.dirstr(k,rean),"plot.score/member",member,"_k",k,"_",dist,"_",nbdays,"day_Repartition_TWS_0.01.png"),width=18,height=8,units="in",res=72)
  #plot(moy[1,1:(ncol(moy)*0.01)],xlab="Days",ylab="TWS",type="l",ylim = c(0.15,0.31),main="Repartition des TWS journaliers moyennes entre séquences - 1%")
  #lines(moy[2,1:(ncol(moy)*0.01)],col = "red")
  #lines(moy[3,1:(ncol(moy)*0.01)],col = "royalblue")
  #abline(v = c(0.01,0.005,0.002,0.001)*ncol(moy))
  #graphics.off()
  
}

plot.TWSgeo<-function(k,dist,nbdays,start,end,rean){
  
  # Import
  precip <- get.precip(nbdays,start,end)
  geo <- get.descriptor(descriptor = "TWSgeonei",k = k,dist = dist,nbdays = nbdays,start = start,end = end,
                 standardize = F,rean = rean)
  # Traitement
  ind <- list(2)
  ind[[1]] <- sort(precip,decreasing=T,index.return=T)$ix[1:(62*12)]
  ind[[2]] <- sort(precip,decreasing=T,index.return=T)$ix[1:62]
  
  # Boxplot
  png(filename = paste0("2_Travail/20CR/Rresults/overall/plot.TWSgeo/boxplot_",dist,"_member",member,"_k",k,"_mean",nbdays,"day.png"),width = 400,height = 400,units = "px")
  plot(1,1,type="n",xlim=c(0,3),ylim=range(geo),xlab="",ylab="",xaxt="n")
  grid()
  boxplot(geo,at=0.5,xlim=c(0,3),col="lightsteelblue1",add=T)
  boxplot(geo[ind[[1]]],at=1.5,add=T,col="lightsteelblue1")
  boxplot(geo[ind[[2]]],at=2.5,add=T,col="lightsteelblue1")
  axis(side = 1,at = c(0.5,1.5,2.5),labels = c("all","monthly max","yearly max"))
  graphics.off()
}

# Concatenation et aggregation de netcdf sous R
reshape.netcdf <- function(z = "1000"){
  
  stock <- list()
  for(i in 1:4){ # Importation et stockage des arrays
    print(paste0("Importation ",i,"/4"))
    tmp        <- nc_open(paste0("2_Travail/ERA-20C/Data/Z",z,"_era20c",i,".nc"))
    stock[[i]] <- ncvar_get(nc = tmp, varid = "z", start = c(1,1,1), count = c(tmp$dim$longitude$len,tmp$dim$latitude$len,tmp$dim$time$len))
  }
  
  # Groupement en une seule array
  print("Concatenation")
  arr <- do.call(abind,stock)
  
  # Passage en hauteur
  print("Convertion en hauteur")
  arr <- round(arr/9.81,2)
  
  # Agregation au pdt journalier
  print("Agregation")
  fourth_dim <- dim(arr)[3]/4
  dim(arr)   <- c(dim(arr)[1:2],4,fourth_dim)
  arr_final  <- apply(arr,c(1,2,4),mean)
  
  # Export au format NetCDF
  time.seq <- seq(from = 9, to = fourth_dim*24, by = 24)
  dim1 <- ncdim_def(name = "longitude", units = tmp$dim$longitude$units, vals = tmp$dim$longitude$vals)
  dim2 <- ncdim_def(name = "latitude", units = tmp$dim$latitude$units, vals = tmp$dim$latitude$vals)
  dim3 <- ncdim_def(name = "time", units = tmp$dim$time$units, vals = as.integer(time.seq), calendar = "gregorian")
  
  res.form   <- ncvar_def(name = "hgt", units = "gpm", dim = list(dim1,dim2,dim3),prec = "double",longname = "Geopotential")
  res.create <- nc_create(filename = paste0("2_Travail/ERA-20C/Data/ERA20C_HGT",z,"_1900_2010_daily.nc"), vars = res.form)
  ncvar_put(nc = res.create, varid = res.form,vals = arr_final)
  
}

# Lance les fonctions souhaitees
run<-function(k,dist,nbdays,str,radtype,start,end,rean){
  
  descr<-list(
    c("persnei","celnei"),
    c("celnei","singnei"),
    c("persnei","singnei"),
    c("celnei","rsingnei"),
    c("persnei","rsingnei"),
    c("singnei","rsingnei")
  )
  
  descr<-lapply(descr,paste.descr,str)
  print(descr)
  
  
  for (i in 1:length(descr)){
    fit.empir(rean = rean, k = k, descriptors = descr[[i]], dist = dist,
              nbdays = nbdays, start = start, end = end, radtype = M)
    
    plot.empir(rean = rean, k = k, descriptors = descr[[i]], dist = dist,
              nbdays = nbdays, start = start, end = end, radtype = M)
    
    compute_crps(descriptors = descr[[i]], k = unique(k), dist = unique(dist), nbdays = nbdays,
                 start = start, end = end, radtype = M, rean = unique(rean))
  }
  
}

# Lance les fonctions souhaitees
run.TL<-function(k,dist,nbdays,str,radAna,radInd,start,end,rean){
  
  descr<-list(
    c("cel","sing"),
    c("snei","sing"),
    c("cel","rsing"),
    c("snei","rsing"),
    c("sing","rsing"),
    c("sing","sing"),
    c("rsing","rsing"),
    c("cel","sing","rsing"),
    c("snei","sing","rsing")
  )
  
  descr<-lapply(descr,paste.descr,str)
  print(descr)
  
  for (i in 1:length(descr)){
    fit.empir.TL(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, start = start, end = end, radAna = radAna, radInd = radInd, rean = rean)
    compute_crps_TL(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, start = start, end = end, radAna = radAna, radInd = radInd, rean = rean)
  }
  
  compare.crps.TL(k = k, dist = dist, nbdays = nbdays, start = start, end = end, radAna = radAna, radInd = radInd, rean = rean)
}

# Lance les fonctions souhaitees
run.comb<-function(k,dist,nbdays,str,start,end,p,rad,rean){
  
  descr<-list(
    c("cel","sing"),
    c("snei","sing"),
    c("cel","rsing"),
    c("snei","rsing"),
    c("sing","rsing"),
    c("sing","sing"),
    c("rsing","rsing"),
    c("cel","sing","rsing"),
    c("snei","sing","rsing")
  )
  
  descr<-lapply(descr,paste.descr,str)
  print(descr)
  
  for (i in 1:length(descr)){
    fit.empir.comb(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, start = start, end = end, p = p, rad = rad, rean = rean)
    compute_crps_comb(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, start = start, end = end, p = p, rad = rad, rean = rean)
  }
  
  #compare.crps.TL(k = k, dist = dist, nbdays = nbdays, start = start, end = end, radAna = radAna, radInd = radInd, rean = rean)
}

# Renvoie les indices des X% des journees/sequences les plus analogues (pour les sequences, somme des scores de chaque ieme journee des deux sequences)
save.nei.A<-function(k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean){
  gc()
  
  print(paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  dist.vec<-getdist(k,dist,start,end,rean)
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  dates<-getdates(start,end)
  N<-length(dates)
  gc()
  
  U<-c(0,(N-1):1);
  sU<-sapply(1:(N-1),function(x) sum(U[1:x]))
  gc()
  
  nei<-list()
  namnei<-c("01","02","05","1","2","5","10")
  
  for (i in 1:length(namnei)) nei[[i]]<-list()
  names(nei)<-namnei
  
  ddi.mat<-NULL
  for (i in 1:(N-nbdays+1)){
    if (nbdays>1){
      if (i==1){
        for (i2 in 1:nbdays) ddi.mat<-rbind(ddi.mat,getdist4i(i2,dist.vec,N,sU))
      }
      else {
        ddi.mat<-ddi.mat[-1,]
        ddi.mat<-rbind(ddi.mat,getdist4i(i+nbdays-1,dist.vec,N,sU))
      }
      #ddi<-apply(ddi.mat,2,sum)
      tmp.mat<-matrix(NA,ncol=N-nbdays+1,nrow=nbdays)
      tmp.mat[1,]<-ddi.mat[1,1:(N-nbdays+1)]
      for (i2 in (2:nbdays)) tmp.mat[i2,]<-ddi.mat[i2,i2:(N-nbdays+i2)]
      ddi<-apply(tmp.mat,2,sum)
    }
    else{
      ddi<-getdist4i(i,dist.vec,N,sU)
    }
    
    if (seasonal){
      j<-selind_season(i,len=30,dates)
      di<-ddi[j]
    }
    else {
      di<-ddi
      j<-1:length(ddi)
    }
    n<-length(di)
    
    gc()
    
    soso<-sort(di,index.return=TRUE)
    for (nn in names(nei)) nei[[nn]][[i]]<-j[soso$ix[1:(str2prop(nn)*n)]]
    
    if (i %% 50==0) {
      print(i)
    }
  }
  save(nei,file=paste0(get.dirstr(k,rean),"save.nei.A/nei_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
}

# Genere les memes fichiers que compute.dist mais pour des sequences de nbdays jours
save.score.A <- function(k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean){
  
  gc()
  
  dist.vec<-getdist(k,dist,start,end,rean,nbdays=1)
  gc()
  
  dist.vec<-unlist(dist.vec)
  gc()
  
  dates<-getdates(start,end)
  N<-length(dates)
  gc()
  
  U<-c(0,(N-1):1)
  sU<-sapply(1:(N-1),function(x) sum(U[1:x]))
  gc()
  
  dist.list <-vector("list",N-nbdays)
  ddi.mat<-NULL
  
  for (i in 1:(N-nbdays)){
    if (i %% 50==0) {
      print(i)
    }
    
    if (nbdays>1){
      if (i==1){
        for (i2 in 1:nbdays) ddi.mat<-rbind(ddi.mat,getdist4i(i2,dist.vec,N,sU))
      }
      else {
        ddi.mat<-ddi.mat[-1,]
        ddi.mat<-rbind(ddi.mat,getdist4i(i+nbdays-1,dist.vec,N,sU))
      }
      
      tmp.mat<-matrix(NA,ncol=N-nbdays+1,nrow=nbdays)
      tmp.mat[1,]<-ddi.mat[1,1:(N-nbdays+1)]
      for (i2 in (2:nbdays)) tmp.mat[i2,]<-ddi.mat[i2,i2:(N-nbdays+i2)]
      ddi<-apply(tmp.mat,2,sum)
    }
    else{
      ddi<-getdist4i(i,dist.vec,N,sU)
    }
    
    if (seasonal){
      j<-selind_season(i,len=30,dates)
      di<-ddi[j]
    }
    else {
      di<-ddi
    }
    
    dist.list[[i]] <- di[(i+1):length(di)]
    gc()
    
  }
  
  save(dist.list, file = paste0("2_Travail/20CR/Rresults/compute_dist/",ifelse(nbdays>1,paste0(nbdays,"day_"),""),dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
}

# Indices des 16 differentes classes de wp (saison a risque de 1 a 8, non risque de 9 a 16)
save.wp.ind <- function(start="1950-01-01",end="2011-12-31"){
  
  wp <- get.wp(start,end,risk=T)
  delta <- get.delta(ref="1851-01-01",start = start)
  wp.ind <- list()
  
  for(i in 1:16){
    wp.ind[[i]] <- which(wp==i)
    write(x = wp.ind[[i]]+delta, file = paste0("2_Travail/Rresults/save.wp.ind/WP",i,".txt"), ncolumns = length(wp.ind[[i]]),sep = ",")
  }
  
  save(wp.ind,file="2_Travail/Rresults/save.wp.ind/WP_ind.Rdata")
  
}

# Tranforme la liste de score en matrice diagonale
score.to.mat <- function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  load(file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",ifelse(nbdays>1,paste0(nbdays,"day_"),""),dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
  for(i in 1:length(dist.list)) dist.list[[i]] <- c(rep(0,length(dist.list)+1-length(dist.list[[i]])),dist.list[[i]])
  dist.list <- do.call(rbind,dist.list)
  last <- rep(0,ncol(dist.list))
  dist.list <- rbind(dist.list,last) # ajout de la derniere ligne
  dist.list2 <- t(dist.list)
  dist.list <- dist.list + dist.list2
}

# Trace le boxplot des max max annuels par mois de precipitation journaliere
season.at.risk <- function(start="1950-01-01",end="2011-12-31"){
  
  precip <- get.precip(1,start,end)
  mon_max <- aggregate(precip,by=list(substr(getdates(start,end),1,7)),FUN=max)
  
  png("2_Travail/Rresults/season.at.risk/boxplot_daily_rainfall.png",width = 670,height = 500,units = "px")
  boxplot(mon_max[,2]~as.factor(substr(mon_max[,1],6,7)),outline=F,names=month.abb,xlab="Month",ylab="Precipitation (mm)")
  points(aggregate(mon_max[,2],by=list(substr(mon_max[,1],6,7)),FUN=quantile,probs=0.9),col="red",pch=19) # on ajoute le quantile 90%
  abline(v = 2.5,col="red",lty=2)
  abline(v = 8.5,col="red",lty=2)
  legend("topleft","90% quantile",col="red",pch=19,bty="n")
  title("Annual maxima")
  graphics.off()
  
}

# Renvoie les indices des dates dans les + ou - len jours de chaque annee
selind_season<-function(i,len,dates){
  #print(dates[i])
  i0<-which(substr(dates,6,10)==substr(dates[i],6,10))
  
  ind<-NULL
  N<-length(dates)
  #print(N)
  for (i1 in i0) {
    if (i1-len>=1 & i1+len<=N) ind<-c(ind,(i1-len):(i1+len)) # on ajoute plus ou moins len jours a la date (longueur de 2*len+1)
    else if (i1-len<1) ind<-c(ind,1:(2*len+1)) # si la fenetre bloque sur le debut du vecteur, on prend les premiers 2*len+1 jours
    else if (i1+len>N) ind<-c(ind,(N-2*len+1):N) # si la fenetre bloque sur la fin du vecteur, on prend les derniers 2*len+1 jours
  }
  
  #dates[ind]
  ind
}

# Ecrit un mot sur une carte a la position souhaitee
shadowtext <- function(x, y=NULL, labels, col='white', bg='black',
                       theta= seq(pi/4, 2*pi, length.out=8), r=0.1, ... ) {
  
  xy <- xy.coords(x,y)
  xo <- r*strwidth('x')
  yo <- r*strheight('x')
  
  for (i in theta) {
    text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
  }
  text(xy$x, xy$y, labels, col=col, ... )
}

# Renvoie un numerique a partir d'une chaine de caracteres (0.005, 0.01)
str2prop<-function(str){
  
  if(substr(str,1,1) == "0") {
    prop <- as.numeric(str)*0.001
  } else { prop <- as.numeric(str)*0.01}
 
  prop
}

# Trace la distribution des indicateurs pour les pluies decroissantes, et l'ellipse des pluies dans le plan des indicateurs
summary.distrib<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean,rev=FALSE){
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE,rean)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE,rean)
  
  precip<-get.precip(nbdays,start,end)
  
  if (descriptor1=="cel" | descriptor2=="cel") {
    descr1<-descr1[-1]
    descr2<-descr2[-1]
    precip<-precip[-1]
  }
  idsort<-sort(precip,index.return=TRUE,decreasing=ifelse(rev,FALSE,TRUE))$ix
  print(range(precip[idsort[10000]]))
  N<-length(precip)
  print(N)
  
  descr<-cbind(descr1,descr2)
  
  meanx.all<-mean(descr[,1],na.rm=TRUE)
  meany.all<-mean(descr[,2],na.rm=TRUE)
  sdx.all<-sd(descr[,1],na.rm=TRUE)
  sdy.all<-sd(descr[,2],na.rm=TRUE)
  
  #print(summary(param))
  #nvec<-c(30,100,500,1000,5000,10000)
  if(!rev){
  nvec<-c(30,100,500,seq(1000,N,by=500))
  } else{nvec<-c(1849,seq(2000,N,by=500))}
  #nvec<-seq(30,N,by=100)
  
  meanx.vec<-NULL
  meany.vec<-NULL
  sdx.vec<-NULL
  sdy.vec<-NULL
  
  for (n in nvec){
    idx1<-idsort[1:n]
    
    meanx.vec<-c(meanx.vec,mean(descr[idx1,1],na.rm=TRUE))
    meany.vec<-c(meany.vec,mean(descr[idx1,2],na.rm=TRUE))
    sdx.vec<-c(sdx.vec,sd(descr[idx1,1],na.rm=TRUE))
    sdy.vec<-c(sdy.vec,sd(descr[idx1,2],na.rm=TRUE))
  }
  graphics.off()
  #quartz(width=12,heigh=4.5);par(mar=c(5.1,5.1,2.1,2.1))
  png(file=paste0(get.dirstr(k,rean),"summary.distrib/summaryplot_",ifelse(rev,"rev_",""),descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),".png"),width=22,heigh=9,units="in",res=300)
  par(mfrow=c(1,3))
  
  plot(descr[,1],descr[,2],col="gray",xlab=descriptor1,ylab=descriptor2)
  lines(meanx.vec,meany.vec)
  points(meanx.vec,meany.vec)
  points(meanx.all,meany.all,col="red")
  #segments(x0=meanx.vec-sdx.vec,y0=meany.vec,x1=meanx.vec+sdx.vec,y1=meany.vec)
  #segments(x0=meanx.vec,y0=meany.vec-sdy.vec,x1=meanx.vec,y1=meany.vec+sdy.vec)
  #segments(x0=meanx.all-sdx.all,y0=meany.all,x1=meanx.all+sdx.all,y1=meany.all,col="red")
  #segments(x0=meanx.all,y0=meany.all-sdy.all,x1=meanx.all,y1=meany.all+sdy.all,col="red")
  
  require(plotrix)
  draw.ellipse(meanx.all,meany.all,sdx.all,sdy.all,border="red")
  for (i in 1:length(meanx.vec)) draw.ellipse(meanx.vec[i],meany.vec[i],sdx.vec[i],sdy.vec[i],border="black")
  
  
  ylim<-range(meanx.all+sdx.all*c(-1,1),meanx.vec+sdx.vec,meanx.vec-sdx.vec)
  plot(nvec,meanx.vec,xlim=c(30,N),ylim=ylim,ylab=descriptors[1])
  segments(x0=nvec,y0=meanx.vec-sdx.vec,x1=nvec,y1=meanx.vec+sdx.vec)
  points(N,meanx.all,col="red")
  segments(x0=N,y0=meanx.all-sdx.all,x1=N,y1=meanx.all+sdx.all,col="red")
  
  ylim<-range(meany.all+sdy.all*c(-1,1),meany.vec+sdy.vec,meany.vec-sdy.vec)
  plot(nvec,meany.vec,xlim=c(30,N),ylim=ylim,ylab=descriptors[2])
  segments(x0=nvec,y0=meany.vec-sdy.vec,x1=nvec,y1=meany.vec+sdy.vec)
  points(N,meany.all,col="red")
  segments(x0=N,y0=meany.all-sdy.all,x1=N,y1=meany.all+sdy.all,col="red")
  
  dev.off()
  #dev2bitmap(file=paste0(get.dirstr(k),"summary.distrib/summaryplot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),".pdf"),width=12,heigh=4.5,type="pdfwrite")
  
}

# Applique un test statistique surla distribution des pluies a differents endroits du plan des indicateurs
test.distrib<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rand=FALSE){
  library(sm)
  library(ks)
  
  if (rand) randstr<-"_rand"
  else randstr<-""
  
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  print(paste0(get.dirstr(k),"fit.empir",get.CVstr(TRUE),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_fix2.Rdata"))
  load(file=paste0(get.dirstr(k),"fit.empir",get.CVstr(TRUE),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_fix2.Rdata"))
  param<-param[,"nbnei"]
  #print(summary(param))
  #print(param)
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE)
  
  precip<-get.precip(nbdays,start,end)
  
  if (descriptor1=="cel" | descriptor2=="cel") {
    descr1<-descr1[-1]
    descr2<-descr2[-1]
    precip<-precip[-1]
    param<-param[-1]
  }
  idsort<-sort(precip,index.return=TRUE,decreasing=TRUE)$ix
  print(range(precip[idsort[10000]]))
  N<-length(precip)
  
  descr<-cbind(descr1,descr2)
  
  meanx.all<-mean(descr[,1],na.rm=TRUE)
  meany.all<-mean(descr[,2],na.rm=TRUE)
  sdx.all<-sd(descr[,1],na.rm=TRUE)
  sdy.all<-sd(descr[,2],na.rm=TRUE)
  
  #print(summary(param))
  nvec<-c(30,100,500,1000,5000,10000)
  pval.mat<-NULL
  graphics.off()
  #quartz(width=18,heigh=8)
  png(file=paste0(get.dirstr(k),"test.distrib/densityplot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),randstr,".png"),width=22,heigh=10,units="in",res=72)
  par(mfrow=c(2,6))
  pval.mat<-NULL
  pval.mat2<-NULL
  
  for (n in nvec){
    if (rand) idx1<-sample(1:N,n)
    else idx1<-idsort[1:n]
    
    if (TRUE){
      pval<-NULL
      pval2<-NULL
      for (i in 1:100){
        #print(i)
        set.seed(i)
        idx2<-sample(1:N,n)
        y<-param[c(idx1,idx2)]
        g<-c(rep(1,length(idx1)),rep(2,length(idx2)))
        tmp<-sm.density.compare(y, g, model="equal",display="none")
        pval<-c(pval,tmp$p)
        pval2<-c(pval2,kde.test(x1=descr[idx1,],x2=descr[idx2,])$pvalue)
      }
      print(pval)
      print(pval2)
      pval.mat<-cbind(pval.mat,pval)
      pval.mat2<-cbind(pval.mat2,pval2)
    }
    plot(descr1,descr2,col=getcol(param),xlab=descriptor1,ylab=descriptor2)
    addscale(param)
    points(descr1[idx1],descr2[idx1],pch=21,bg="black")
    title(n)
    
    plot(density(param[idx1]),main="Density")
    lines(density(param),col="red")
    if (n==nvec[1]){
      legend("topright",legend=c("sample","all"),lty=1,col=c("black","red"))
    }
  }
  dev.off()
  pval.list<-list(pval.sm=pval.mat,pval.ks=pval.mat2)
  save(pval.list,file=paste0(get.dirstr(k),"test.distrib/pval_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),randstr,".Rdata"))
  #dev2bitmap(file=paste0(get.dirstr(k),"test.distrib/densityplot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),randstr,".pdf"),type="pdfwrite",width=18,heigh=8)
  
  colnames(pval.mat)<-colnames(pval.mat2)<-nvec
  quartz(width=10,heigh=10)
  par(mfrow=c(2,2))
  #sm
  boxplot(pval.mat,ylim=c(0,1),ylab="pvalue")
  proppval<-round(apply(pval.mat,2,function(v) sum(v<0.05))/nrow(pval.mat)*100)
  title("sm test")
  plot(nvec,proppval,ylab="% pvalue<0.05",xlab="",ylim=c(0,100))
  lines(nvec,proppval)
  abline(h=5,lty=2)
  
  #ks
  boxplot(pval.mat2,ylim=c(0,1),ylab="pvalue")
  proppval<-round(apply(pval.mat2,2,function(v) sum(v<0.05))/nrow(pval.mat2)*100)
  title("ks test")
  plot(nvec,proppval,ylab="% pvalue<0.05",xlab="",ylim=c(0,100))
  lines(nvec,proppval)
  abline(h=5,lty=2)
  
  dev2bitmap(file=paste0(get.dirstr(k),"test.distrib/pvalplot_",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),randstr,".pdf"),type="pdfwrite",width=10,heigh=10)
}

# Repartition sur la pariode et dans l'annee des plus fortes precipitations
time.extr <- function(start="1950-01-01",end="2011-12-31"){
  
  precip <- get.precip(3,start,end)
  dates <- as.Date(getdates(start,end))
  
  # 62 plus fortes sequences
  dates.extr <- dates[sort(precip,decreasing = T,index.return=T)$ix[1:62]]
  ann <- as.numeric(format(dates.extr,"%Y"))
  mois <- as.numeric(format(dates.extr,"%m"))
  
  pdf(file=paste0("2_Travail/Rresults/time.extr/repartition_extr_",start,"_",end,".pdf"),width = 7,height = 4)
  hist(x = ann,breaks = seq(1950.5,2011.5,1),col="cornflowerblue",border = "royalblue",
       xlab = "Year",main="62 max 3-days precipitations")
  hist(mois,breaks=0:12,col="cornflowerblue",border = "royalblue",
       xlab = "Month",main="62 max 3-days precipitations",xaxt="n")
  axis(side = 1,at = 0.5:11.5,month.abb)
  
  # max annuel
  ann <- format(dates,"%Y")[1:22643]
  pos <- aggregate(precip,by=list(ann),which.max)
  mon <- as.integer(pos[,2]/365.25*12+1)
  hist(mon,breaks = 0:12,col="cornflowerblue",border = "royalblue",
       xlab = "Month",main="Annual max 3-days precipitations",xaxt="n")
  axis(side = 1,at = 0.5:11.5,month.abb)
  graphics.off()
}
