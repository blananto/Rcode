library(ncdf4) # NetCDF
library(fBasics) #timPalette
library(ismev) #gp.fit
library(splancs) #inpip
library(fields) #Tps
library(KernSmooth) #bkde2D
library(plotrix) #draw.circle
library(mev) #egp2.fit
require(akima) #interp: bilinear interpolation
library(zoo) #rollapply
library(RANN) #nn2
library(scoringRules) #crps
library(evd) #fpot

# Calcule et trace les CRPSS pour analogie indicateurs pour les differents couples d'indicateurs et differents rayons 
compare.crps<-function(which="",k=NULL,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="",standardize=TRUE,CV=TRUE,rean){
  descr<-list(
    c("sing","sing"),
    c("cel","sing"),
    c("snei","sing"),
    c("rsing","rsing"),
    c("cel","rsing"),
    c("snei","rsing"),
    c("sing","rsing"),
    c("cel","sing","rsing"),
    c("snei","sing","rsing"),
    c("A","A"))
  
  ndesc<-length(descr)
  
  precip<-get.precip(nbdays,start,end)
  soso<-sort(precip,index.return=TRUE,decreasing=TRUE) # Classement des pluies de la plus forte a la plus faible avec indices
  
  colnam<-c("all ecdf","pos ecdf","p0 binom")
  
  if (which=="rad"){ # Graphique pour les differents rayons
    rads<-paste0("nrn",c("02","05","1","2"))
    radstr.save<-""
    kvec<-k
    kstr.save<-paste0("_k",k)
    distance <-dist
    diststr.save<-dist
    legs<-c("02","05","1","2")
    dirstr.save<-get.dirstr(k,rean)
  }
  if (which=="") { # Graphique pour un rayon donne
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<-k
    kstr.save<-paste0("_k",k)
    distance <-dist
    diststr.save<-dist
    legs<-"05"
    #if (k==1) colnam<-c(colnam,"pos gamma","pos egp")
    dirstr.save<-get.dirstr(k,rean)
  }
  if (which=="dist") { # Graphique pour les differentes distances (pour un rayon donne)
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<-k
    kstr.save<-paste0("_k",k)
    distance <-c("TWS","RMSE")
    diststr.save<-dist # pour stocker dist et l'utiliser pour le chemin du graphique, meme apres la boucle
    legs<-c("TWS","RMSE")
    dirstr.save<-get.dirstr(k,rean)
  }
  if (which=="k") { # dans ce cas, on compare pour un meme rayon les CRPSS des trois k (500,1000 et 500+1000) et enregistrement en amont dans l'arborescence
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<- c(1,2,4)
    kstr.save<-""
    distance <-dist
    diststr.save<-dist
    legs<-kvec
    dirstr.save<-get.dirstr(rean = rean)
  }
  if (which=="k_dist") { # dans ce cas, on compare pour un meme rayon les CRPSS des trois k (500,1000 et 500+1000) et enregistrement en amont dans l'arborescence
    rads<-radtype
    radstr.save<-paste0("_",radtype)
    kvec<- c(1,2)
    kstr.save<-""
    distance <-c("TWS","RMSE")
    diststr.save<-dist
    legs<-c("500hPa_TWS","1000hPa_TWS","500hPa_RMSE","1000hPa_RMSE")
    dirstr.save<-get.dirstr(rean = rean)
  }
  
  
  for (nam in colnam){ # pour chacune des lois de probabilite ajustee
    crps.mat<-NULL
    coln<-NULL
    for(dist in distance){ # pour les distances selectionnees
      for (k in kvec){ # pour les k selectionnes
        for (rad in rads){ # pour les rayons selectionnes
          print(c(nam,str))
          descriptors<-descr
          descriptors<-lapply(descr,paste.descr,"05")
          print(descriptors)
          
          for (i in 1:length(descriptors)){ # Import et recuperation des scores selon les differents couples d'indicateurs pour une distribution
            print(descriptors[[i]])
            print(paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".Rdata"))
            if (substr(descriptors[[i]][1],1,1)=="A") load(file=paste0(get.dirstr(k,rean),"compute_crps-CV.A/05_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
            else load(file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors[[i]],collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",rad,".Rdata"))
            crps.mat<-cbind(crps.mat,crps[,nam])
            coln<-c(coln,paste0(descr[[i]],collapse="-"))
          }
        }
      }
    }
    print(crps.mat[1:10,])
    colnames(crps.mat)<-coln
    
    idx<-which(apply(crps.mat,1,function(v) all(!is.na(v)))) # indice des lignes completes sans NA
    print(length(idx))
    #  idx<-intersect(idx,which(!is.na(fit$loglik.i)))
    #  print(length(idx))
    
    meancrps<-apply(crps.mat[idx,],2,mean) # moyenne du CRPS par couple d'indicateur sur toutes les dates (sans NA sur la ligne)
    
    pos <- soso$ix[soso$x > 0.01] # on ne conserve que les indices des pluies positives
    len <- length(pos)
    idx.0<-intersect(idx,pos[(len-0.05*len):len]) # moyenne du CRPS par couple d'indicateurs sur les 5% des pluies les plus faibles
    meancrps.0<-apply(crps.mat[idx.0,],2,mean)
    
    idx.1<-intersect(idx,soso$ix[1:(62*12)]) # moyenne du CRPS par couple d'indicateurs sur l'equivalent du maximum par mois de pluie (les 12*62 plus fortes pluies)
    meancrps.1<-apply(crps.mat[idx.1,],2,mean)
    
    idx.2<-intersect(idx,soso$ix[1:62]) # moyenne du CRPS par couple d'indicateurs sur l'equivalent du maximum par an de pluie (les 62 plus fortes pluies)
    meancrps.2<-apply(crps.mat[idx.2,],2,mean)
    
    load(file=paste0("2_Travail/",rean,"/Rresults/compute_score_climato/score_mean",nbdays,"day_",start,"_",end,".Rdata")) # Import du score climatologique
    
    #normalize<-mean(score$crps[idx,nam])
    #normalize.1<-mean(score$crps[idx.1,nam])
    #normalize.2<-mean(score$crps[idx.2,nam])
    
    if (substr(nam,1,3)=="pos") { # si "pos" dans le nom du modele
      normalize<-mean(score$crps[idx,"pos ecdf"]) # Moyenne du score climato sur toutes dates non NA
      normalize.0<-mean(score$crps[idx.0,"pos ecdf"]) # sur 5% les plus faibles
      normalize.1<-mean(score$crps[idx.1,"pos ecdf"]) # sur 12*62
      normalize.2<-mean(score$crps[idx.2,"pos ecdf"]) # sur 62
    }
    else{ # sinon, on fait la meme chose?
      normalize<-mean(score$crps[idx,nam])
      normalize.0<-mean(score$crps[idx.0,nam])
      normalize.1<-mean(score$crps[idx.1,nam])
      normalize.2<-mean(score$crps[idx.2,nam])
    }
    
    
    if (seasonal) seastr<-"seasonal"
    else seastr<-"overall"
    
    filename<-paste0(dirstr.save,"compare.crps/",which,match(nam,colnam),"_",diststr.save,"_member",member,kstr.save,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),radstr.save,".png")
    if (substr(nam,1,3)=="pos"){ # si "pos", trois fenetres graphiques
      png(file=filename,width=18,height=8,units="in",res=72)
      par(mar=c(9,4,4,2),mfrow=c(1,4))
    }
    else{ # sinon, une seule fenetre graphique
      png(file=filename,width=7,height=8,units="in",res=72)
      par(mar=c(5.5,4,4,2))
    }
    plot(c(1,ndesc),c(0,max(1-meancrps/normalize)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45)) # Definition des proprietes du graphique
    for (j in 1:length(legs)) { # pour chaque legs, on ajoute une courbe
      points(1:ndesc,1-meancrps[(j-1)*ndesc+(1:ndesc)]/normalize,col=j,pch=19)
      lines(1:ndesc,1-meancrps[(j-1)*ndesc+(1:ndesc)]/normalize,col=j)
    }
    #legend(legend=1:3,col=1:3)
    axis(2)
    axis(1,labels=coln,at=1:length(meancrps),las=3)
    box()
    
    title(paste0(nam," ",seastr))
    grid()
    legend("topleft",legend=legs,col=1:length(legs),pch=19)
    abline(v=1:ndesc,lty=2,col=gray(0.5))
    
    if (substr(nam,1,3)=="pos"){ # si "pos" dans nam, on ajoute deux autres graphiques au png avec les 62*12 puis 12 plus fortes pluies
      plot(c(1,ndesc),c(0,max(1-meancrps.0/normalize.0)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.70)) # graphique 62
      for (j in 1:length(legs)) {
        points(1:ndesc,1-meancrps.0[(j-1)*ndesc+(1:ndesc)]/normalize.0,col=j,pch=1)
        lines(1:ndesc,1-meancrps.0[(j-1)*ndesc+(1:ndesc)]/normalize.0,col=j)
      }
      axis(2)
      axis(1,labels=coln,at=1:length(meancrps.0),las=3)
      box()
      title(paste0(nam," min 5% > 0.01"))
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      grid()
      
      plot(c(1,ndesc),c(0,max(1-meancrps.1/normalize.1)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45)) # graphique 62*12
      for (j in 1:length(legs)) {
        points(1:ndesc,1-meancrps.1[(j-1)*ndesc+(1:ndesc)]/normalize.1,col=j,pch=19)
        lines(1:ndesc,1-meancrps.1[(j-1)*ndesc+(1:ndesc)]/normalize.1,col=j)
      }
      axis(2)
      axis(1,labels=coln,at=1:length(meancrps.1),las=3)
      box()
      title(paste0(nam," monthly max (62*12 largest)"))
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      grid()
      
      plot(c(1,ndesc),c(0,max(1-meancrps.2/normalize.2)),axes=FALSE,xlab="",ylab="",type="n",ylim=c(0,0.45)) # graphique 62
      for (j in 1:length(legs)) {
        points(1:ndesc,1-meancrps.2[(j-1)*ndesc+(1:ndesc)]/normalize.2,col=j,pch=19)
        lines(1:ndesc,1-meancrps.2[(j-1)*ndesc+(1:ndesc)]/normalize.2,col=j)
      }
      axis(2)
      axis(1,labels=coln,at=1:length(meancrps.2),las=3)
      box()
      title(paste0(nam," yearly max (62 largest)"))
      abline(v=1:ndesc,lty=2,col=gray(0.5))
      grid()
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

# Calcul des indicateurs
compute_criteria<-function(k,dist,start="1950-01-01",end="2011-12-31",update=FALSE,rean){
  gc()
  
  print(paste0(get.dirstr(k,rean),"compute_criteria/criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
  
  if (update) {
    load(file=paste0(file=paste0(get.dirstr(k,rean),"compute_criteria/criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
    coln<-colnames(criteria)
  }
  
  dist.vec<-getdist(k,dist,start,end,rean)
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  dates<-getdates(start,end)
  N<-length(dates)
  gc()
  
  U<-c(0,(N-1):1); # U = 0, 22644, 22643, 22642, ...
  sU<-sapply(1:(N-1),function(x) sum(U[1:x])) # on fait la somme de U[1], U[1:2], etc pour obtenir la position de la derniere distance qui separe chaque date
  gc()
  
  if (!update) coln.new<-c("cel","mind","sing05","sing1","sing2","sing5","lsing05","lsing1","lsing2","lsing5","pers05","pers1","pers2","pers5","q05","q1","q2","q5","pcel","pnei05","pnei1","pnei2","pnei5","snei05","snei1","snei2","snei5")
  if (update) {
    coln.new<-c("sing10","lsing10","pers10","q10","pnei10","snei10")
  }
  
  criteria.new<-matrix(NA,ncol=length(coln.new),nrow=N)
  colnames(criteria.new)<-coln.new
  
  for (i in 1:N){
    ddi<-getdist4i(i,dist.vec,N,sU)
    if (seasonal){
      j<-selind_season(i,len=30,dates)
      di<-ddi[j]
    }
    else di<-ddi
    n<-length(di)
    
    gc()
    
    soso<-sort(di,index.return=TRUE) # classement par plus petit score, et donne les positions
    qi05<-di[soso$ix[(0.005*n)]]     # quantile 0.5%
    qi1<-di[soso$ix[(0.01*n)]]       # 1%
    qi2<-di[soso$ix[(0.02*n)]]       # 2%
    qi5<-di[soso$ix[(0.05*n)]]       # 5%
    qi10<-di[soso$ix[(0.1*n)]]       # 10%
    # Donne le rayon du cercle (les 0.5% les plus proches, les 1% les plus proches...)
    
    gc()
    
    idi05<-soso$ix[2:(0.005*n)] # recupere la position des 0.5% les plus proches
    idi1<-soso$ix[2:(0.01*n)]   # des 1% les plus proches
    idi2<-soso$ix[2:(0.02*n)]   # des 2% les plus proches
    idi5<-soso$ix[2:(0.05*n)]   # des 5% les plus proches
    idi10<-soso$ix[2:(0.1*n)]   # des 10% les plus proches
    
    gc()
    
    ri05<-rle(as.integer(di<=qi05)) # calcule le temps que reste le score a l'interieur ou l'exterieur du quantile 0.5%
    ri1<-rle(as.integer(di<=qi1))   # du quantile 1%
    ri2<-rle(as.integer(di<=qi2))   # du quantile 2%
    ri5<-rle(as.integer(di<=qi5))   # du quantile 5%
    ri10<-rle(as.integer(di<=qi10)) # du quantile 10%
    
    tmp<-NULL
    for (cc in coln.new){
      
      # Celerite
      if (cc=="cel") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,ddi[i-1])} # score de la journee avec la veille
      
      # Minimum distance
      if (cc=="mind") tmp<-c(tmp,di[soso$ix[2]]) # score minimum obtenu avec la meilleure journee analogue
      
      # Singularite
      if (cc=="sing05") tmp<-c(tmp,mean(di[idi05])) # moyenne des distances (scores) des 0.5% les plus proches
      if (cc=="sing1") tmp<-c(tmp,mean(di[idi1]))   # des 1% les plus proches
      if (cc=="sing2") tmp<-c(tmp,mean(di[idi2]))   # des 2% les plus proches
      if (cc=="sing5") tmp<-c(tmp,mean(di[idi5]))   # des 5% les plus proches
      if (cc=="sing10") tmp<-c(tmp,mean(di[idi10])) # des 10% les plus proches
      
      # Log Singularite
      if (cc=="lsing05") tmp<-c(tmp,mean(log(di[idi05]))) # moyenne du logarithme des distances des 0.5% les plus proches
      if (cc=="lsing1") tmp<-c(tmp,mean(log(di[idi1])))   # des 1% les plus proches
      if (cc=="lsing2") tmp<-c(tmp,mean(log(di[idi2])))   # des 2% les plus proches
      if (cc=="lsing5") tmp<-c(tmp,mean(log(di[idi5])))   # des 5% les plus proches
      if (cc=="lsing10") tmp<-c(tmp,mean(log(di[idi10]))) # des 10% les plus proches
      
      # Persistence
      if (cc=="pers05") tmp<-c(tmp,mean(ri05$lengths[ri05$values==1])) # moyenne du temps passe a l'interieur du quantile 0.5%
      if (cc=="pers1") tmp<-c(tmp,mean(ri1$lengths[ri1$values==1]))    # du quantile 1%
      if (cc=="pers2") tmp<-c(tmp,mean(ri2$lengths[ri2$values==1]))    # du quantile 2%
      if (cc=="pers5") tmp<-c(tmp,mean(ri5$lengths[ri5$values==1]))    # du quantile 5%
      if (cc=="pers10") tmp<-c(tmp,mean(ri10$lengths[ri10$values==1])) # du quantile 10%
      
      # Quantiles
      if (cc=="q05") tmp<-c(tmp,qi05) # Quantile 0.5%
      if (cc=="q1") tmp<-c(tmp,qi1)   # Quantile 1%
      if (cc=="q2") tmp<-c(tmp,qi2)   # Quantile 2%
      if (cc=="q5") tmp<-c(tmp,qi5)   # Quantile 5%
      if (cc=="q10") tmp<-c(tmp,qi10) # Quantile 10%
      
      # Probabilite celerite: probabilite d'avoir un score en dessous de celui de la veille
      if (cc=="pcel"){if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,ecdf(dim1)(ddi[i-1]))} # calcul de la fonction de repartition empirique (F=i/n, fonction creneau) et probabilite au non depassement de la veille (donne une indication sur son rang dans les analogues)
      
      # Persistence neighbour: probabilite qu'un jour du voisinage soit dans le voisinage de la veille le jour precedent (qu'il suive la meme trajectoire)
      if (cc=="pnei05") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi05-1) %in% idi05m1))} # moyenne de la condition: les veilles des journees analogues sont dans les 0.5% les plus proches de la veille de la journee cible
      if (cc=="pnei1") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi1-1) %in% idi1m1))}    # dans les 1% les plus proches
      if (cc=="pnei2") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi2-1) %in% idi2m1))}    # dans les 2% les plus proches
      if (cc=="pnei5") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi5-1) %in% idi5m1))}    # dans les 5% les plus proches
      if (cc=="pnei10") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi10-1) %in% idi10m1))} # dans les 10% les plus proches
      
      # Persistence neighbour: probabilite qu'un jour du voisinage soit deja dans le voisinage le jour precedent
      if (cc=="snei05") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi05-1) %in% idi05))} # moyenne de la condition: les veilles des journees analogues sont deja dans les 0.5% les plus proches de la journee cible
      if (cc=="snei1") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi1-1) %in% idi1))}    # dans les 1% les plus proches
      if (cc=="snei2") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi2-1) %in% idi2))}    # dans les 2% les plus proches
      if (cc=="snei5") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi5-1) %in% idi5))}    # dans les 5% les plus proches
      if (cc=="snei10") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,mean((idi10-1) %in% idi10))} # dans les 10% les plus proches
      
      # Celerite neighbour: on moyenne la celerite des plus proches voisins
      if (cc=="celnei05") tmp<-c(tmp,mean(criteria[idi05,"cel"],na.rm=TRUE)) # moyenne des celerites des 0.5% les plus proches
      if (cc=="celnei1") tmp<-c(tmp,mean(criteria[idi1,"cel"],na.rm=TRUE))   # des 1% les plus proches
      if (cc=="celnei2") tmp<-c(tmp,mean(criteria[idi2,"cel"],na.rm=TRUE))   # des 2% les plus proches
      if (cc=="celnei5") tmp<-c(tmp,mean(criteria[idi5,"cel"],na.rm=TRUE))   # des 5% les plus proches
      if (cc=="celnei10") tmp<-c(tmp,mean(criteria[idi10,"cel"],na.rm=TRUE)) # des 10% les plus proches
      
      # La singularite relative rsing (ou local dimension) est calculee dans la fonction get.descriptor
      
    }
    
    dim1<-di
    idi05m1<-idi05
    idi1m1<-idi1
    idi2m1<-idi2
    idi5m1<-idi5
    idi10m1<-idi10
    
    criteria.new[i,]<-tmp
    
    if (i %% 50==0) {
      print(i)
      print(criteria.new[(i-49):i,])
    }
    
    dim1<-di
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
  
  save(criteria,file=paste0(file=paste0(get.dirstr(k,rean),"compute_criteria/criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
}

# Calcul des CRPS apres analogie au sens des indicateurs
compute_crps<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean){
  
  print(paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
  
  precip<-get.precip(nbdays,start,end)
  N<-length(precip)
  
  descr <- matrix(NA,N,length(descriptors))
  for (i in 1:length(descriptors)) descr[,i]<-get.descriptor(descriptors[i],k,dist,nbdays,start,end,standardize,rean) # on importe les descripteurs
  
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
  
  load(file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata")) # on importe la loi empirique
  
  crps<-matrix(NA,ncol=5,nrow=N)
  colnames(crps)<-c("all ecdf","pos ecdf","p0 binom","pos gamma","pos egp") # memes colonnes que pour score climato
  for (i in 1:N){
    if (i %%50==0) print(i)
    if (!is.na(descr[i,1]) & !is.na(descr[i,2])){ # si les indicateurs du jour ne sont pas nuls
      tmp<-get.closest(i,descr,precip,CV,nbdays,radtype) # on prend les voisins selon un certain rayon
      pi<-precip[tmp$idx]
      crps[i,"all ecdf"]<-crps_sample(precip[i],pi) # crps all ecdf
      crps[i,"p0 binom"]<-crps_binom(precip[i]==0,1,param[i,"p0"]) # crps p0 binom
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
  save(crps,file=paste0(get.dirstr(k,rean),"compute_crps",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
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

# Calcul des distances de maniere generique
compute_dist.gen<-function(k,dist,start="1950-01-01",end="2011-12-31",rean){
  if (k %in% 1:2) {
    if (dist %in% c("TWS","RMSE")){
      if (dist=="TWS") dist.list<-compute_TWS(k,start,end,rean)
      if (dist=="RMSE") dist.list<-compute_RMSE(k,start,end,rean)
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

# Trouve le jour de NetCDF associe a la date
date_to_number<-function(nc,day,rean){
  which(getdays(nc,rean)==day)
}

# Calcul des parametres de la loi empirique selon le voisinage au sens des indicateurs
fit.empir<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,radtype="nrn05",CV=TRUE,rean){
  
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
      parami<-c(ni,1-ni1/ni,mean(pi[pi>0]),sd(pi[pi>0]),skewness(pi[pi>0]),kurtosis(pi[pi>0]),tmp$radius) # remplissage de la ligne: nbre de voisins, proba jour sec, moyenne et ecart-type des pluies non nulles, coefficient d'assymetrie et d'applatissement, rayon du cercle
    }
    else parami<-rep(NA,7) # si un des deux indicateurs est nul, ligne du tableau parametre vide
    param[i,]<-parami
  }
  colnames(param)<-c("nbnei","p0","mean","sd","skewness","kurtosis","radius")
  
  save(param,file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",paste0(descriptors,collapse="_"),"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(standardize),"_",radtype,".Rdata"))
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

# fit une cdf empirique avec deux niveaux d'analogie: analogie classique pus indicateurs
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
getdata=function(k,day0,day1=day0,rean){
  num0<-date_to_number(nc[[k]],day0,rean)
  num1<-date_to_number(nc[[k]],day1,rean)
  N<-length(num0:num1)
  infowind<-getinfo_window(k,rean)
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
  if(rean == "ERA20C") orig <- "1900"
  substr(as.POSIXct(nc$dim$time$vals*3600,origin=paste0(orig,'-01-01 00:00'),tz="GMT"),1,10) #format "YYYY-MM-DD"
}

# Importe la liste contenant la distance souhaitee
getdist<-function(k,dist,start="1950-01-01",end="2011-12-31",rean){
  load(file=paste0("2_Travail/",rean,"/Rresults/compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
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
getinfo_window=function(k,rean){ 
  # lat et lon correspondent aux degres de longitude et latitude pour lesquels on a les donnees des geopot 500 et 1000
  lon<-nc[[1]]$dim$lon$vals # de -30 a 50 deg E tranches de 2 deg -> 41 valeurs
  lat<-nc[[1]]$dim$lat$vals #de 24 a 72 deg N tranches de 2 deg -> 25 valeurs
  # time=data500$dim$time$vals # 58804 valeurs par tranches de 24h (du 01-01-1851 au 31-12-2011), valeurs a 9h chaque jour
  # parametre fenetre d'analogie
  
  c_lon<-6 # centre de la fenetre longitude
  c_lat<-44 # centre de la fenetre latitude
  
  if(k==1) {# 500hPA
    d_lon<-32 # on prend 32 en longitude pour 500hPa
    d_lat<-16 # on prend 16 en latitude pour 500hPa
  }
  if (k==2){ # 1000 hPA
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
    #print(str2prop(radstr)*length(precip))
    n<-nn2(descr[idef,],t(descr[i,]),str2prop(radstr)*length(precip)) # Nearest Neighbour Search
    radius<-max(n$nn.dists) # le rayon est la distance maximum selectionnee dans n
    idx<-sort(idef[n$nn.idx]) # idx = indice des jours dans le rayon du cercle
  }
  if (CV) idx<-setdiff(idx,(i-nbdays+1):(i+nbdays-1)) # si Cross-Validation, on retire les indices des sequences de nbdays jours dans lesquelles le jour i intervient
  pi<-precip[idx]
  ni0<-sum(pi>0) # nombre de jours pluvieux dans le voisinage selectionne
  if (ni0<20) {
    #we increase the radius to embrace >=20 positive rainfall
    ipos<-which(precip>0 & !is.na(descr[,1]) & !is.na(descr[,2])) # selection des pluies positives avec indicateurs non nuls
    if (CV) ipos<-setdiff(ipos,(i-nbdays+1):(i+nbdays-1)) # si cross-validation, on enleve les sequences incluent le jour i
    n<-nn2(descr[ipos,],t(descr[i,]),20) # recherche de 20 voisins parmi les jours pluvieux
    radius<-max(n$nn.dists)
    #idx<-which(d<=radius^2)
    idx<-sort(ipos[n$nn.idx]) # nouveau voisinage
  }
  lclose<-list(idx=idx,radius=radius)
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
  
}

# Import et mise en forme des indicateurs: calculs supplementaires, moyenne glissante sur trois jours
get.descriptor<-function(descriptor,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",standardize=TRUE,rean){
  
  load(file=paste0(file=paste0(get.dirstr(k,rean),"compute_criteria/criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
  
  if (descriptor=="celcelnei") descr<-criteria[,"cel"]/criteria[,"celnei"] # division de la celerite du jour par la celerite moyenne des voisins
  else if (substr(descriptor,1,5)=="rsing") {
    nb<-substr(descriptor,6,nchar(descriptor))
    descr<-criteria[,paste0("sing",nb)]/criteria[,paste0("q",nb)] # calcul de rsing (local dimension), en divisant la singularite par le quantile souhaite
  }
  else if (substr(descriptor,1,6)=="rlsing") {
    nb<-substr(descriptor,7,nchar(descriptor))
    descr<-criteria[,paste0("lsing",nb)]-log(criteria[,paste0("q",nb)]) # calcul de rlsing (log local dimension)
  }
  else descr<-criteria[,descriptor] # sinon, on selectionne tout simplement la colonne specifiee dans descriptor
  
  if (nbdays>1) descr<-rollapply(descr,width=nbdays,FUN=mean) # moyenne glissante de l'indicateur sur nbdays jours
  
  if (standardize) return(descr/sd(descr,na.rm=TRUE)) # si standardise, on divise l'indicateur par l'ecart type de sa serie
  else return(descr)
}

# Modifie le repertoire d'ecriture des resultats selon k et overall/seasonal 
get.dirstr<-function(k=NULL,rean){
  if (seasonal) {dirstr<-paste0("2_Travail/",rean,"/Rresults/seasonal/")
  }else dirstr<-paste0("2_Travail/",rean,"/Rresults/overall/")
  if (!is.null(k)) dirstr<-paste0(dirstr,"k",k,"/")
  dirstr
}

# Renvoie les indices des nbre events extremes de pluie a partir d'une certaine reference
get.ind.extr <- function(nbre, ref = "1950-01-01", nbdays, start="1950-01-01", end="2011-12-31"){
  
  # Classement
  precip <- get.precip(nbdays, start, end)
  tri <- sort(precip, decreasing = TRUE, index.return = TRUE)
  ind <- tri$ix[1:nbre]
  
  # Reference
  if(ref != "1950-01-01"){
    dt <- get.delta(ref, start)
    ind <- ind + dt
  }
  return(ind)
}

# Importation des donnees de precipitation
get.precip<-function(nbdays,start="1950-01-01",end="2011-12-31"){
  
  precip <- read.csv(file=paste0("2_Travail/Precip/Isere@Grenoble_cum",nbdays,"day_1950-01-01_2011-12-31.csv"))
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

# Calcule les gradients pour 500 (k=1) ou 1000 (k=2) pour tous les jours
grad<-function(k,day0,day1,rean){
  x=getdata(k,day0,day1,rean)  #data500 ou data1000, 3 dims
  gradlon=(makegrad(x,1)) #selon lon
  gradlat=(makegrad(x,2)) #selon lat
  return(list(gradlon=gradlon,gradlat=gradlat))
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
  
  nc<<-list(nc500=nc500,nc1000=nc1000)
}

# Sort la matrice des gradients pour input= data500 ou data1000 sur notre fenetre d'analogie, l=1 longitude, l=2 latitude
makegrad<-function(mat,l){
  dims<-dim(mat)
  if (l==1) gradmat<-mat[2:dims[[1]],,]-mat[1:(dims[[1]]-1),,] # l=1 1ere ligne= difference entre 2eme et 1ere ligne, 2eme ligne=difference entre la 3eme et le 2eme ligne --> une ligne de moins que input
  if (l==2) gradmat<-mat[,2:dims[[2]],]-mat[,1:(dims[[2]]-1),] # l=2 1ere colonne= difference entre 2eme et 1ere colonne, 2eme colonne=difference entre la 3eme et la 2eme colonne  --> une colonne de moins que input
  return(gradmat)
}

# Trace la repartition des parametres empiriques dans le plan des indicateurs (un couple d'indicateur par png)
plot.empir<-function(descriptors,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,rean){
  
  descriptor1<-descriptors[1]
  descriptor2<-descriptors[2]
  
  print(paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata"))
  load(file=paste0(get.dirstr(k,rean),"fit.empir",get.CVstr(CV),"/",descriptor1,"_",descriptor2,"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype,".Rdata")) # importation des parametres de la loi empirique
  param<-param[,colnames(param)!="kurtosis"]
  
  descr1<-get.descriptor(descriptor1,k,dist,nbdays,start,end,standardize=FALSE,rean)
  descr2<-get.descriptor(descriptor2,k,dist,nbdays,start,end,standardize=FALSE,rean)
  
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
    ind.extr <- get.ind.extr(nbre = 3,nbdays = 3)
    ind.extr <- as.vector(ind.extr)
    points(descr1[ind.extr],descr2[ind.extr],pch=19)
    title(colnames(param)[i])
  }
  dev.off()
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

# Ajout d'une chaine de caractere a l'indicateur (05,1,2,5,10 a tous les indicateurs sauf cel)
paste.descr<-function(descr,str){
  descr[descr!="cel"]<- paste0(descr[descr!="cel"],str)
  descr
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
    fit.empir(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, radtype = radtype, start = start, end = end, rean = rean)
    #fit.loglik.gamma(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, radtype = radtype, start = start, end = end, rean = rean)
    #fit.loglik.egp(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, radtype = radtype, start = start, end = end, rean = rean)
    #fit.loglik.p0(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, radtype = radtype, start = start, end = end, rean = rean)
    compute_crps(descriptors = descr[[i]], k = k, dist = dist, nbdays = nbdays, radtype = radtype, start = start, end = end, rean = rean)
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

# Renvoie un numerique a partir d'une chaine de caracteres (0.005, 0.01)
str2prop<-function(str){
  
  if(substr(str,1,1) == "0") {
    prop <- as.numeric(str)*0.001
  } else { prop <- as.numeric(str)*0.01}
 
  prop
}