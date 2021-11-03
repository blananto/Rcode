source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Combine les scores de 500hPa et 1000hPa
combine.dist.k3 <- function(wk1=0.5,wk2=0.5,dist,rean,start,end){
  
  # Imports
  print("Imports")
  dist.k1 <- getdist(k = 1,dist = dist,start = start,end = end,rean = rean,threeday = F,period = "past")
  gc()
  dist.k2 <- getdist(k = 2,dist = dist,start = start,end = end,rean = rean,threeday = F,period = "past")
  gc()
  
  # Traitement
  print("Traitement")
  dist.list <- mapply(function(x,y){as.integer(wk1*x+wk2*y)},x=dist.k1,y=dist.k2)
  
  # Export
  print("Export")
  save(dist.list,file=paste0("2_Travail/1_Past/",rean,"/compute_dist/",dist,"_",rean,"_k",3,"_",start,"_",end,".Rdata"))
  gc()
}

# Comparaison des precip reconstituees avec et sans schaake
compare.schaake <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,season=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  load(paste0(get.dirstr(k,rean,period),"reordering.precip/precip_final_",type,"_",ifelse(type=="gamma",paste0("q",q.threshold,"_"),""),ifelse(season,"season_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana_",nbana,"_",rean,"_k",k,"_",dist,"_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  load(paste0(get.dirstr(k,rean,period),"generate.precip/precip_",type,"_",ifelse(type=="gamma",paste0("q",q.threshold,"_"),""),ifelse(season,"season_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana_",nbana,"_",rean,"_k",k,"_",dist,"_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  precip.ana.sample <- do.call(rbind,precip.ana.sample)
  
  pos <- match(dates.ana,dates.rean) # on ne garde que periode commune avec precip pour comparer à l'observé
  precip.ana.final <- precip.ana.final[pos,]
  precip.ana.sample <- precip.ana.sample[pos,]
  
  # Graphiques autocorrelation
  # Autocorrelation saisonniere 1 jour
  print("Grahiques autocorrelation")
  
  png(filename = paste0(get.dirstr(k,rean,period),"compare.schaake/compare_autocor_",type,"_",ifelse(type=="gamma",paste0("q",q.threshold,"_"),""),ifelse(season,"season_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana_",nbana,"_",rean,"_k",k,"_",dist,"_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".png"),width = 9,height = 4,units = "in",res = 600)
  par(mfrow=c(1,2))
  autocor.month.obs <- aggregate(precip,by=list(substr(dates.ana,6,7)),function(v){acf(v,lag.max=1,plot=F)$acf[2]})
  autocor.month.ana <- apply(precip.ana.sample,2,function(v){aggregate(v,by=list(substr(dates.ana,6,7)),function(w){acf(w,lag.max=1,plot=F)$acf[2]})})
  autocor.month.ana <- do.call(cbind,lapply(autocor.month.ana,function(v){v[,2]}))
  tmp.ana <- apply(autocor.month.ana,1,range)
  autocor.month.schaake <- apply(precip.ana.final,2,function(v){aggregate(v,by=list(substr(dates.ana,6,7)),function(w){acf(w,lag.max=1,plot=F)$acf[2]})})
  autocor.month.schaake <- do.call(cbind,lapply(autocor.month.schaake,function(v){v[,2]}))
  tmp.schaake <- apply(autocor.month.schaake,1,range)
  
  ylim <- c(0,0.5)
  plot(autocor.month.obs,type="l",col="red",lwd=2,ylim=ylim,xaxt="n",xlab="Month",ylab="Autocorrelation (1-day lag)",main=nam2str(bv))
  grid()
  axis(side = 1,at = 1:12,labels = month.abb)
  polygon(c(1:12,12:1), c(tmp.ana[2,],rev(tmp.ana[1,])),
          col = adjustcolor("black",alpha.f=0.2) ,border = F)
  polygon(c(1:12,12:1), c(tmp.schaake[2,],rev(tmp.schaake[1,])),
          col = adjustcolor("blue",alpha.f=0.2) ,border = F)
  lines(autocor.month.obs,col="red",lwd=2)
  
  # Autorrelation plusieurs jours
  autocor.obs <- acf(precip,lag.max=5,plot=F)$acf[-1]
  autocor.ana <- apply(precip.ana.sample,2,function(v){acf(v,lag.max=5,plot=F)$acf[-1]})
  tmp.ana <- apply(autocor.ana,1,range)
  autocor.schaake <- apply(precip.ana.final,2,function(v){acf(v,lag.max=5,plot=F)$acf[-1]})
  tmp.schaake <- apply(autocor.schaake,1,range)
  
  ylim=c(0,0.4)
  plot(autocor.obs,type="l",col="red",lwd=2,ylim=ylim,xlab="Day",ylab="Autocorrelation",main=nam2str(bv))
  grid()
  polygon(c(1:5,5:1), c(tmp.ana[2,],rev(tmp.ana[1,])),
          col = adjustcolor("black",alpha.f=0.2) ,border = F)
  polygon(c(1:5,5:1), c(tmp.schaake[2,],rev(tmp.schaake[1,])),
          col = adjustcolor("blue",alpha.f=0.2) ,border = F)
  lines(autocor.obs,col="red",lwd=2)
  legend("topright",inset=.02,c("Obs","Random","Schaake Shuffle"),col=c("red",adjustcolor("black",alpha.f=0.2),adjustcolor("blue",alpha.f=0.2)),lty=1,lwd=2,bty="n",cex=0.8)
  graphics.off()
  
  # Graphiques cumuls 3 et 5 jours
  print("Grahiques cumuls plusieurs pas de temps")
  png(filename = paste0(get.dirstr(k,rean,period),"compare.schaake/compare_cumuls_",type,"_",ifelse(type=="gamma",paste0("q",q.threshold,"_"),""),ifelse(season,"season_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana_",nbana,"_",rean,"_k",k,"_",dist,"_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".png"),width = 12,height = 8,units = "in",res = 600)
  par(mfcol=c(2,4))
  ti <- c(2,3,4,5)
  nyear <- length(unique(substr(dates.ana,1,4)))
  rp.axis <- c(1,2,5,10,25,60)
  
  for(i in 1:length(ti)){
    
    print(i)
    
    # Donnees lissees
    precip.cum <- rollapply(precip,ti[i],sum)
    precip.ana.cum <- rollapply(precip.ana.sample,ti[i],sum)
    precip.schaake.cum <- rollapply(precip.ana.final,ti[i],sum)
    
    # Mise en forme boxplot
    tmp <- precip.cum
    tmp.ana <- c(precip.ana.cum)
    tmp.schaake <- c(precip.schaake.cum)
    length(tmp) <- length(tmp.ana)
    tab <- cbind(tmp,tmp.ana,tmp.schaake)
    colnames(tab) <- c("Obs","Random","Schaake \nShuffle")
    
    # Boxplot
    boxplot(tab,outline=F,col=c("red",adjustcolor("black",alpha.f=0.2),adjustcolor("blue",alpha.f=0.2)),ylab=paste0(ti[i],"-day precipitation (mm)"),main=nam2str(bv))
    grid();par(new=T)
    boxplot(tab,outline=F,col=c("red",adjustcolor("black",alpha.f=0.2),adjustcolor("blue",alpha.f=0.2)),ylab=paste0(ti[i],"-day precipitation (mm)"),main=nam2str(bv))
    
    # Mise en forme distribution
    rp.sort <- sort(1/(1-ecdf(precip.cum)(precip.cum))/365.25) # periode de retour empirique
    rp.sort[is.infinite(rp.sort)] <- nyear
    precip.sort <- sort(precip.cum)
    
    precip.sort.ana <- apply(precip.ana.cum,2,sort)
    mini.ana <- apply(precip.sort.ana,1,min)
    maxi.ana <- apply(precip.sort.ana,1,max)
    
    precip.sort.schaake <- apply(precip.schaake.cum,2,sort)
    mini.schaake <- apply(precip.sort.schaake,1,min)
    maxi.schaake <- apply(precip.sort.schaake,1,max)
    
    # Distribution
    ylim <- c(0,range(maxi.ana,maxi.schaake)[2])
    plot(rp.sort,precip.sort,log="x",xaxt="n",type="n",ylim=ylim,xlab="Return Period (year)",ylab=paste0(ti[i],"-day precipitation (mm/day)"),main=paste0(nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
    grid(nx = NA,ny = NULL);abline(v=rp.axis,lty=3,col="grey")
    axis(side = 1,at = rp.axis,labels = rp.axis)
    polygon(c(rp.sort,rev(rp.sort)), c(maxi.ana,rev(mini.ana)),col = adjustcolor("black",alpha.f=0.2),border = F)
    polygon(c(rp.sort,rev(rp.sort)), c(maxi.schaake,rev(mini.schaake)),col = adjustcolor("blue",alpha.f=0.2),border = F)
    lines(rp.sort,precip.sort,col="red",lwd=3)
    if(i==1) legend("topleft",inset=.02,c("Obs","Random","Schaake Shuffle"),col=c("red",adjustcolor("black",alpha.f=0.2),adjustcolor("blue",alpha.f=0.2)),lty=1,lwd=2,bty="n")
  }
  
  graphics.off()
  
  # Graphiques frequence et longueur des sequences seches
  print("Grahiques sequences seches")
  png(filename = paste0(get.dirstr(k,rean,period),"compare.schaake/compare_dry_sequences_",type,"_",ifelse(type=="gamma",paste0("q",q.threshold,"_"),""),ifelse(season,"season_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana_",nbana,"_",rean,"_k",k,"_",dist,"_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".png"),width = 5,height = 4,units = "in",res = 600)
  seuil <- 0
  tmp <- rle(precip<=seuil)
  precip.dry <- tmp$lengths[tmp$values]
  precip.ana.dry <- unlist(apply(precip.ana.sample,2,function(v){tmp<-rle(v<=seuil);tmp$lengths[tmp$values]}))
  precip.schaake.dry <- unlist(apply(precip.ana.final,2,function(v){tmp<-rle(v<=seuil);tmp$lengths[tmp$values]}))
  length(precip.dry) <- length(precip.ana.dry) <- length(precip.schaake.dry) <- max(length(precip.ana.dry),length(precip.schaake.dry))
  tab <- cbind(precip.dry,precip.ana.dry,precip.schaake.dry)
  colnames(tab) <- c("Obs","Random","Schaake Shuffle")
  
  boxplot(tab,outline=F,col=c("red",adjustcolor("black",alpha.f=0.2),adjustcolor("blue",alpha.f=0.2)),ylab="Dry sequences length (days)",main=nam2str(bv))
  grid();par(new=T)
  boxplot(tab,outline=F,col=c("red",adjustcolor("black",alpha.f=0.2),adjustcolor("blue",alpha.f=0.2)),ylab="Dry sequences length (days)",main=nam2str(bv))
  graphics.off()
}

# Calcule les RL20 et mean max selon la GEV pour les scenarios analogues
compute.gev.ana <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,start,end,var2=NULL,nbana2=NULL){

  # Dates utiles
  dates <- getdates(start,end)
  year <- as.numeric(unique(substr(dates,1,4)))
  
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import GEV
  load(file = get.path(fit_gev_ana = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,var2 = var2,nbana2 = nbana2,start = start,end = end))
  
  # Estimations RL20/mean et export
  meanGEV.fct<-function(mu,sig,xi){mu+sig/xi*(gamma(1-xi)-1)}
  sais <- c("winter","spring","summer","autumn")
  val.gev <- vector(mode="list",length=length(fit))
  
  for(i in 1:length(fit)){
    print(paste0(i,"/",length(fit)))
    
    for(j in 1:length(sais)){
      tab <- matrix(data = NA,nrow = 2,ncol = 4)
      colnames(tab) <- c("Year","RL20","Mean","pval")
      tab[,1] <- c(year[1],year[length(year)])
      
      ny <- nrow(fit[[i]][[j]]$fitbest$vals)
      tab[1,2] <- qgev(1-1/20,fit[[i]][[j]]$fitbest$vals[1,1],fit[[i]][[j]]$fitbest$vals[1,2],fit[[i]][[j]]$fitbest$vals[1,3])
      tab[2,2] <- qgev(1-1/20,fit[[i]][[j]]$fitbest$vals[ny,1],fit[[i]][[j]]$fitbest$vals[ny,2],fit[[i]][[j]]$fitbest$vals[ny,3])
      tab[1,3] <- meanGEV.fct(fit[[i]][[j]]$fitbest$vals[1,1],fit[[i]][[j]]$fitbest$vals[1,2],fit[[i]][[j]]$fitbest$vals[1,3])
      tab[2,3] <- meanGEV.fct(fit[[i]][[j]]$fitbest$vals[ny,1],fit[[i]][[j]]$fitbest$vals[ny,2],fit[[i]][[j]]$fitbest$vals[ny,3])
      tab[1,4] <- fit[[i]][[j]]$pval
      val.gev[[i]][[j]] <- tab
    }
    names(val.gev[[i]]) <- sais
  }
  
  save(val.gev,file = get.path(compute_gev_ana = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,var2 = var2,nbana2 = nbana2,start = start,end = end))
}

# Calcule les RL20 et mean max selon la GEV pour l'observe
compute.gev.obs <- function(bv="Isere",nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  year <- as.numeric(unique(substr(dates,1,4)))
  
  # Import GEV
  load(file = paste0("2_Travail/1_Past/Rresults/fit.gev.obs/fit_gev_",bv,"_mean",nbdays,"day_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
  
  # Estimations RL20/mean et export
  meanGEV.fct<-function(mu,sig,xi){mu+sig/xi*(gamma(1-xi)-1)}
  sais <- c("winter","spring","summer","autumn")
  val.gev <- vector(mode="list",length=length(sais))
  names(val.gev) <- sais
  
  for(i in 1:length(sais)){
    tab <- matrix(data = NA,nrow = 2,ncol = 4)
    colnames(tab) <- c("Year","RL20","Mean","pval")
    tab[,1] <- c(year[1],year[length(year)])
    
    ny <- nrow(fit[[i]]$fitbest$vals)
    tab[1,2] <- qgev(1-1/20,fit[[i]]$fitbest$vals[1,1],fit[[i]]$fitbest$vals[1,2],fit[[i]]$fitbest$vals[1,3])
    tab[2,2] <- qgev(1-1/20,fit[[i]]$fitbest$vals[ny,1],fit[[i]]$fitbest$vals[ny,2],fit[[i]]$fitbest$vals[ny,3])
    tab[1,3] <- meanGEV.fct(fit[[i]]$fitbest$vals[1,1],fit[[i]]$fitbest$vals[1,2],fit[[i]]$fitbest$vals[1,3])
    tab[2,3] <- meanGEV.fct(fit[[i]]$fitbest$vals[ny,1],fit[[i]]$fitbest$vals[ny,2],fit[[i]]$fitbest$vals[ny,3])
    tab[1,4] <- fit[[i]]$pval
    val.gev[[i]] <- tab
  }
  
  save(val.gev,file = paste0("2_Travail/1_Past/Rresults/compute.gev.obs/compute_gev_",bv,"_mean",nbdays,"day_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
}

# Calcule les RL20 et mean max selon la GEV pour l'observe, pour tous les BVs
compute.gev.obs.all <- function(nbdays=1,start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  year <- as.numeric(unique(substr(dates,1,4)))
  
  # Import GEV
  load(file = paste0("2_Travail/1_Past/Rresults/fit.gev.obs/fit_gev_all_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  # Estimations RL20/mean et export
  meanGEV.fct<-function(mu,sig,xi){mu+sig/xi*(gamma(1-xi)-1)}
  sais <- c("winter","spring","summer","autumn")
  val.gev <- vector(mode="list",length=length(fit))
  
  for(i in 1:length(fit)){
    print(paste0(i,"/",length(fit)))
    for(j in 1:length(sais)){
      print(sais[j])
      tab <- matrix(data = NA,nrow = 2,ncol = 4)
      colnames(tab) <- c("Year","RL20","Mean","pval")
      tab[,1] <- c(year[1],year[length(year)])
      
      ny <- nrow(fit[[i]][[j]]$fitbest$vals)
      tab[1,2] <- qgev(1-1/20,fit[[i]][[j]]$fitbest$vals[1,1],fit[[i]][[j]]$fitbest$vals[1,2],fit[[i]][[j]]$fitbest$vals[1,3])
      tab[2,2] <- qgev(1-1/20,fit[[i]][[j]]$fitbest$vals[ny,1],fit[[i]][[j]]$fitbest$vals[ny,2],fit[[i]][[j]]$fitbest$vals[ny,3])
      tab[1,3] <- meanGEV.fct(fit[[i]][[j]]$fitbest$vals[1,1],fit[[i]][[j]]$fitbest$vals[1,2],fit[[i]][[j]]$fitbest$vals[1,3])
      tab[2,3] <- meanGEV.fct(fit[[i]][[j]]$fitbest$vals[ny,1],fit[[i]][[j]]$fitbest$vals[ny,2],fit[[i]][[j]]$fitbest$vals[ny,3])
      tab[1,4] <- fit[[i]][[j]]$pval
      val.gev[[i]][[j]] <- tab
    }
    names(val.gev[[i]]) <- sais
  }
  save(val.gev,file = paste0("2_Travail/1_Past/Rresults/compute.gev.obs/compute_gev_all_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Convertit les TWS en integer 10^9 pour reduire la memoire utilisee
convert_dist_past <- function(k,dist,rean){
  
  # Import
  print("Import")
  dates.dist <- get.start.end.rean(rean,"past","dist")
  load(paste0("2_Travail/1_Past/",rean,"/compute_dist/",dist,"_",rean,"_k",k,"_",dates.dist[1],"_",dates.dist[2],"_old.Rdata"))
  
  # Conversion
  if (!is.integer(dist.list[[1]])) {
    print("Conversion")
    dist.list <- lapply(dist.list,function(v){as.integer(v*(10^9))})
  }
  gc()
  
  # Sauvegarde
  print("Sauvegarde")
  save(dist.list,file=paste0("2_Travail/1_Past/",rean,"/compute_dist/",dist,"_",rean,"_k",k,"_",dates.dist[1],"_",dates.dist[2],".Rdata"))
  rm(dist.list)
  gc()
}

# Calage loi exponentielle sur les distributions de precip reconstruite par analogie classique
fit.exp <- function(k,dist,nbdays,nbana=0.2,nbmini=10,seuil=0,bv,spazm=T,rean,period="past",season=F,dP=F,moments=T){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2010-12-31")
  
  # Import des precip
  load(get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana))
  
  # Fonction pour maximum de vraisemblance: revient au meme que methode des moments pour loi exponentielle
  if(!moments) nloglik.fct<-function(param,x){
    -sum(dexp(x, rate = param, log = TRUE)) # - car l'optimisation cherche a minimiser cette valeur
  }
  
  # Traitement
  fit <- matrix(NA,ncol=2,nrow=N)
  colnames(fit) <- c("rate","p0")
  
  for (i in 1:N){
    
    if (i %%50==0) print(i)
    
    # Precipitations et stats utiles
    precip.i <- na.omit(precip.ana[[i]]) # on retire les NA
    precip.pos <- precip.i[precip.i>seuil]
    if(length(precip.pos)>nbmini){ # on ne fit pas la loi avec moins de nbmini valeurs positives
      p0 <- 1-length(precip.pos)/length(precip.i)
      mean.pos <- mean(precip.pos)
      rate <- 1/mean.pos
      
      if(!moments){
        opt <- optim(rate,nloglik.fct,method="L-BFGS-B",x=precip.pos,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) 
        rate <- opt$par
      }
      fit[i,]<- c(rate,p0)
    }
  }
  
  save(fit,file = get.path(fit_exp = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Calage loi gamma sur les distributions de precip reconstruite par analogie classique
fit.gamma <- function(k,dist,nbdays,nbana=0.2,nbmini=10,seuil=0,bv,spazm=T,rean,period="past",season=F,dP=F,moments=F,var2=NULL,nbana2=NULL){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2017-12-31")
  
  # Import des precip
  load(get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana,var2 = var2,nbana2 = nbana2))
  
  # Fonction pour maximum de vraisemblance
  if(!moments) nloglik.fct<-function(param,x){
    -sum(dgamma(x, shape=param[2], scale = param[1], log = TRUE)) # - car l'optimisation cherche a minimiser cette valeur
  }
  
  # Traitement
  fit <- matrix(NA,ncol=3,nrow=N)
  colnames(fit) <- c("scale","shape","p0")
  
  for (i in 1:N){
    
    if (i %%50==0) print(i)
    
    # Precipitations et stats utiles
    precip.i <- na.omit(precip.ana[[i]]) # on retire les NA
    precip.pos <- precip.i[precip.i>seuil]
    if(length(precip.pos)>nbmini){ # on ne fit pas la loi avec moins de nbmini valeurs positives
      p0 <- 1-length(precip.pos)/length(precip.i)
      mean.pos <- mean(precip.pos)
      var.pos <-var(precip.pos)
      scale <- var.pos/mean.pos
      shape <- mean.pos^2/var.pos
      
      if(!moments){
        opt <- optim(c(scale,shape),nloglik.fct,method="L-BFGS-B",x=precip.pos,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) 
        scale <- opt$par[1];shape <- opt$par[2]
      }
      fit[i,]<- c(scale,shape,p0)
    }
  }
  
  save(fit,file = get.path(fit_gamma = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,var2 = var2,nbana2 = nbana2,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Calage loi gamma sur les distributions de precip reconstruite par analogie classique, pour tous les BVs
fit.gamma.all <- function(k,dist,nbdays,nbana=0.2,nbmini=10,seuil=0,rean,period="past",season=F,dP=F,moments=F,var2=NULL,nbana2=NULL){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2017-12-31")
  
  # Import des precip
  load(get.path(save_precip_ana_all = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
  
  # Fonction pour maximum de vraisemblance
  if(!moments) nloglik.fct<-function(param,x){
    -sum(dgamma(x, shape=param[2], scale = param[1], log = TRUE)) # - car l'optimisation cherche a minimiser cette valeur
  }
  
  # Traitement
  fit <- vector(mode = "list",length=length(precip.ana))
    
  for(i in 1:length(precip.ana)){
    print(paste0(i,"/",length(precip.ana)))
    fit[[i]] <- matrix(NA,ncol=3,nrow=N)
    colnames(fit[[i]]) <- c("scale","shape","p0")
    
    for (j in 1:N){
      if (j %%50==0) print(paste0(j,"/",N))
      
      # Precipitations et stats utiles
      precip.j <- na.omit(precip.ana[[i]][[j]]) # on retire les NA
      precip.pos <- precip.j[precip.j>seuil]
      if(length(precip.pos)>nbmini){ # on ne fit pas la loi avec moins de nbmini valeurs positives
        p0 <- 1-length(precip.pos)/length(precip.j)
        mean.pos <- mean(precip.pos)
        var.pos <-var(precip.pos)
        scale <- var.pos/mean.pos
        shape <- mean.pos^2/var.pos
        
        if(!moments){
          opt <- optim(c(scale,shape),nloglik.fct,method="L-BFGS-B",x=precip.pos,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) 
          scale <- opt$par[1];shape <- opt$par[2]
        }
        fit[[i]][j,]<- c(scale,shape,p0)
      }
    }
  }
  save(fit,file = get.path(fit_gamma_all = T,k = k,rean = rean,period = period,moments = moments,season = season,dP = dP,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,var2 = var2,nbana2 = nbana2,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Calage loi gamma sur les distributions de precip reconstruite par tirage aleatoire (tirage saisonnier possible)
fit.gamma.random <- function(replic=1000,nbdays,nbana=0.2,nbmini=10,seuil=0,bv,spazm=T,rean,period="past",moments=F,season=T,ncores=3){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k=1)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precip
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  # Fonction pour maximum de vraisemblance
  if(!moments) nloglik.fct<-function(param,x){
    -sum(dgamma(x, shape=param[2], scale = param[1], log = TRUE)) # - car l'optimisation cherche a minimiser cette valeur
  }
  
  # Traitement
  nbnei <- round(nbana*0.01*n,0)
  fit <- matrix(NA,ncol=3,nrow=N)
  colnames(fit) <- c("scale","shape","p0")
  
  if(season){
  ind.sea <- list() # indices de recherche des analogues dans la meme saison
  dat <- getdates("2020-01-01","2020-12-31")
  for(i in 1:length(dat)){
  ind.sea[[i]] <- get.ind.season.ana(dat[i],dates.ana)
  }
  pos.dat <- match(substr(dates.rean,6,10),substr(dat,6,10))
  }
  
  print(paste0("Parallelisation sur ",ncores, " coeurs"))
  outfile <- paste0(get.dirstr(k,rean,period),"fit.gamma.random/calcul.txt")
  print(paste0("Logfile for // loop : ",outfile))
  cl <- makeCluster(ncores, outfile=outfile) 
  registerDoParallel(cl)
  
  fit.random <- foreach (i=1:replic) %dopar%{
    print(paste0(i,"/",replic))
    
    for (j in 1:N){
      # Precipitations et stats utiles
      if(season){
        sampl <- ind.sea[[pos.dat[i]]]
      }else{sampl <- 1:N}
      precip.j <- precip[sample(x = sampl,size = nbnei,replace = F)]
      precip.pos <- precip.j[precip.j>seuil]
      if(length(precip.pos)>nbmini){ # on ne fit pas la loi avec moins de nbmini valeurs positives
        p0 <- 1-length(precip.pos)/length(precip.j)
        mean.pos <- mean(precip.pos)
        var.pos <-var(precip.pos)
        scale <- var.pos/mean.pos
        shape <- mean.pos^2/var.pos
        
        if(!moments){
          opt <- optim(c(scale,shape),nloglik.fct,method="L-BFGS-B",x=precip.pos,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) 
          scale <- opt$par[1];shape <- opt$par[2]
        }
        fit[j,]<- c(scale,shape,p0)
      }
    }
    fit
  }
  
  stopCluster(cl)
  fit <- fit.random
  save(fit,file = get.path(fit_gamma_random = T,rean = rean,period = period,moments = moments,bv = bv,season = season,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Calage loi gev sur une serie de precipitation, pour une saison
fit.gev <- function(precip,sais,start,end){
  
  dates <- getdates(start,end)
  year <- as.numeric(unique(substr(dates,1,4)))
  
  # Extraction du max saisonnier
  pos <- get.ind.season(sais = sais,start = start,end = end)
  precip <- precip[pos$pos.season]
  precip[pos$pos.NA] <- NA
  dim(precip) <- c(pos$l.season,pos$n.season)
  max.vec <- apply(precip,2,function(v){max(v,na.rm=T)})
  
  # Covariable
  if(sais=="winter") year <- year[-1]
  covar<-(year-year[1])/(year[length(year)]-year[1])
  
  # Calage GEV
  fit0<-gev.fit(max.vec,show=F) # modele stationnaire
  fit1<-gev.fit(max.vec,as.matrix(covar),mul=1,show=F) # modele non stationnaire, location varie
  fit2<-gev.fit(max.vec,as.matrix(covar),sigl=1,show=F) # modele non stationnaire, scale varie
  fit3<-gev.fit(max.vec,as.matrix(covar),mul=1,sigl=1,show=F) # modele non stationnaire, location et scale varient
  
  covar<-(year-1985)/(year[length(year)]-year[1]);covar[covar<0]<-0 # a partir d'une certaine date
  fit4<-gev.fit(max.vec,as.matrix(covar),mul=1,show=F) # modele non stationnaire, location varie
  fit5<-gev.fit(max.vec,as.matrix(covar),sigl=1,show=F) # modele non stationnaire, scale varie
  fit6<-gev.fit(max.vec,as.matrix(covar),mul=1,sigl=1,show=F) # modele non stationnaire, location et scale varient
  
  pval1<-1-pchisq(2*(fit0$nllh-fit1$nllh),df=1)
  pval2<-1-pchisq(2*(fit0$nllh-fit2$nllh),df=1)
  pval3<-1-pchisq(2*(fit0$nllh-fit3$nllh),df=2)
  pval4<-1-pchisq(2*(fit0$nllh-fit4$nllh),df=1)
  pval5<-1-pchisq(2*(fit0$nllh-fit5$nllh),df=1)
  pval6<-1-pchisq(2*(fit0$nllh-fit6$nllh),df=2)
  
  idx.min<-which.min(c(pval1,pval2,pval3,pval4,pval5,pval6))
  if (idx.min==1) {fitbest <- fit1; pval <- pval1}
  if (idx.min==2) {fitbest <- fit2; pval <- pval2}
  if (idx.min==3) {fitbest <- fit3; pval <- pval3}
  if (idx.min==4) {fitbest <- fit4; pval <- pval4}
  if (idx.min==5) {fitbest <- fit5; pval <- pval5}
  if (idx.min==6) {fitbest <- fit6; pval <- pval6}
  
  list(fit0=fit0,fitbest=fitbest,pval=pval) # modele stationnaire, meilleur modele non stationnaire, pvalue
}

# Calage loi gev sur les precip analogues
fit.gev.ana <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,start,end){
  
  # Dates utiles
  dates <- getdates(start,end)
  
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  if(schaake){
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
  }else{
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  
  pos <- match(dates,dates.rean)
  precip.ana.final <- precip.ana.final[pos,]
  
  # Fit
  sais <- c("winter","spring","summer","autumn")
  fit <- vector(mode="list",length=ncol(precip.ana.final))
  
  for(i in 1:ncol(precip.ana.final)){
    print(paste0(i,"/",ncol(precip.ana.final)))
    for(j in 1:length(sais)){
      fit[[i]][[j]] <- fit.gev(precip = precip.ana.final[,i],sais = sais[j],start = start,end = end)
    }
    names(fit[[i]]) <- sais
  }
  save(fit,file = get.path(fit_gev_ana = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,start = start,end = end))
}

# Calage loi gev sur les precip analogues, en parallele
fit.gev.ana.par <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,start,end,var2=NULL,nbana2=NULL){
  
  # Dates utiles
  dates <- getdates(start,end)
  
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  if(schaake){
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
  }else{
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  
  pos <- match(dates,dates.rean)
  precip.ana.final <- precip.ana.final[pos,]
  
  # Fit
  sais <- c("winter","spring","summer","autumn")
  ncores <- 6
  
  print(paste0("Parallelisation sur ",ncores, " coeurs"))
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  
  fit <- foreach (i=1:ncol(precip.ana.final),.packages = c("ismev","stats")) %dopar%{
    
    source(paste0("2_Travail/1_Past/",rean,"/fit.gev.ana/utils_fit_gev_ana.R"), encoding = 'UTF-8')
    
    fit.i <- vector(mode="list",length = length(sais))
    for(j in 1:length(sais)){
      fit.i[[j]] <- fit.gev(precip = precip.ana.final[,i],sais = sais[j],start = start,end = end)
    }
    names(fit.i) <- sais
    fit.i
  }
  
  stopCluster(cl)
  save(fit,file = get.path(fit_gev_ana = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,var2 = var2,nbana2 = nbana2,start = start,end = end))
}

# Calage loi gev sur les precip observees
fit.gev.obs <- function(bv="Isere",nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31"){
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  
  # Fit
  sais <- c("winter","spring","summer","autumn")
  fit <- vector(mode="list",length=length(sais))
  names(fit) <- sais
  
  for(i in 1:length(sais)){
    print(sais[i])
    fit[[i]] <- fit.gev(precip = precip,sais = sais[i],start = start,end = end)
  }
  save(fit,file = paste0("2_Travail/1_Past/Rresults/fit.gev.obs/fit_gev_",bv,"_mean",nbdays,"day_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
}

# Calage loi gev sur les precip observees, pour tous les BVs
fit.gev.obs.all <- function(nbdays=1,start="1950-01-01",end="2017-12-31"){
  
  # Import des precipitations
  load("2_Travail/Data/Precip/SPAZM/dailyspazmSecteurHydro_1950-2017.Rdata")
  precip <- matdaily;rm(matdaily)
  
  # Fit
  sais <- c("winter","spring","summer","autumn")
  fit <- vector(mode="list",length=ncol(precip))
  
  for(i in 1:ncol(precip)){
    print(paste0(i,"/",ncol(precip)))
    for(j in 1:length(sais)){
      print(sais[j])
      fit[[i]][[j]] <- fit.gev(precip = precip[,i],sais = sais[j],start = start,end = end)
    }
    names(fit[[i]]) <- sais
  }
  save(fit,file = paste0("2_Travail/1_Past/Rresults/fit.gev.obs/fit_gev_all_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Calage loi gaussienne inverse sur les distributions de precip reconstruite par analogie classique
fit.invgauss <- function(k,dist,nbdays,nbana=0.2,nbmini=10,seuil=0,bv,spazm=T,rean,period="past",season=F,dP=F,moments=T){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2010-12-31")
  
  # Import des precip
  load(get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana))
  
  # Fonction pour maximum de vraisemblance: revient au meme que methode des moments pour loi exponentielle
  if(!moments) nloglik.fct<-function(param,x){
    -sum(dinvGauss(x = x,nu = param[1],lambda = param[2],log = T)) # - car l'optimisation cherche a minimiser cette valeur
  }
  
  # Traitement
  fit <- matrix(NA,ncol=3,nrow=N)
  colnames(fit) <- c("mean","shape","p0")
  
  for (i in 1:N){
    
    if (i %%50==0) print(i)
    print(i)
    # Precipitations et stats utiles
    precip.i <- na.omit(precip.ana[[i]]) # on retire les NA
    precip.pos <- precip.i[precip.i>seuil]
    if(length(precip.pos)>nbmini){ # on ne fit pas la loi avec moins de nbmini valeurs positives
      p0 <- 1-length(precip.pos)/length(precip.i)
      mean.pos <- mean(precip.pos)
      var.pos <- var(precip.pos)
      shape <- mean.pos^3/var.pos
      
      if(!moments){
        opt <- optim(c(mean.pos,shape),nloglik.fct,method="L-BFGS-B",x=precip.pos,lower=c(0.000001,0.000001),upper=c(Inf,Inf)) 
        mean.pos <- opt$par[1];shape <- opt$par[2]
      }
      fit[i,]<- c(mean.pos,shape,p0)
    }
  }
  
  save(fit,file = get.path(fit_invgauss = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Genere les precip par analogie par tirage aleatoire dans les analogues ou dans les loi
generate.precip <- function(n=45,type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,season=F,dP=F,rgamma=F,moments=F,var2=NULL,nbana2=NULL){
  
  # On se limite a 45 series generees car 45 analogues uniquement pour nbana=0.2
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2017-12-31")
  
  # Import des precip analogues
  load(get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana,var2 = var2,nbana2 = nbana2))
  if(length(precip.ana[[1]])>n){
    precip.ana <- lapply(precip.ana,function(v){v[1:n]})
  }
  
  # Initialisation: cas empirique
  precip.ana.sample <- precip.ana
  
  # Generation precip loi gamma
  if(type=="gamma"){
    
    # Import des parametres de la loi gamma
    load(get.path(fit_gamma = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,var2 = var2,nbana2 = nbana2,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    
    # Traitement
    for(i in 1:N){
      if(i%%50==0) print(i)
      if(!is.na(fit[i,"shape"])){
        shape <- fit[i,"shape"]
        scale <- fit[i,"scale"] 
        p0 <- fit[i,"p0"]
        
        if(!rgamma){
          p <- runif(n = n,min = 0,max = q.threshold)
          q <- rep(0,length(p)) # on initialise le vecteur avec precip = 0
          ind <- which(p>p0)
          ppos <- ((p-p0)/(1-p0))[ind]
          pq <- qgamma(p=ppos,shape = shape,scale = scale)
          q[ind] <- pq
        }else{
          pq <- rgamma(n = n*(1-p0),shape = shape,scale = scale)
          q <- c(pq,rep(0,n-length(pq)))
        }
        
        precip.ana.sample[[i]] <- q
      }
    }
  }
  
  # Generation precip loi exponentielle
  if(type=="exp"){
    
    # Import des parametres de la loi gamma
    load(get.path(fit_exp = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    
    # Traitement
    for(i in 1:N){
      if(i%%50==0) print(i)
      if(!is.na(fit[i,"rate"])){
        rate <- fit[i,"rate"]
        p0 <- fit[i,"p0"]
        pq <- rexp(n = n*(1-p0),rate = rate)
        q <- c(pq,rep(0,n-length(pq)))
        }
        precip.ana.sample[[i]] <- q
      }
  }
  
  # Generation precip loi exponentielle
  if(type=="invgauss"){
    
    # Import des parametres de la loi gamma
    load(get.path(fit_invgauss = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    
    # Traitement
    for(i in 1:N){
      if(i%%50==0) print(i)
      if(!is.na(fit[i,"mean"])){
        nu <- fit[i,"mean"]
        shape <- fit[i,"shape"]
        p0 <- fit[i,"p0"]
        pq <- rinvGauss(n = n*(1-p0),nu = nu,lambda = shape)
        q <- c(pq,rep(0,n-length(pq)))
      }
      precip.ana.sample[[i]] <- q
    }
  }
  
  # Sortie
  save(precip.ana.sample,file = get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
}

# Genere les precip par analogie par tirage aleatoire dans les analogues ou dans les loi, pour tous les BVs
generate.precip.all <- function(n=45,type="empir",k,dist,nbdays,nbana=0.2,rean,period="past",q.threshold=0.99,seuil=0,season=F,dP=F,rgamma=F,moments=F,var2=NULL,nbana2=NULL){
  
  # On se limite a 45 series generees car 45 analogues uniquement pour nbana=0.2
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2017-12-31")
  
  # Import des precip analogues
  load(get.path(save_precip_ana_all = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
  if(length(precip.ana[[1]][[1]])>n){
    precip.ana <- lapply(precip.ana,function(v){lapply(v,function(w){w[1:n]})})
  }
  
  # Initialisation: cas empirique
  precip.ana.sample <- precip.ana
  
  # Generation precip loi gamma
  if(type=="gamma"){
    
    # Import des parametres de la loi gamma
    load(get.path(fit_gamma_all = T,k = k,rean = rean,period = period,moments = moments,season = season,dP = dP,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,var2 = var2,nbana2 = nbana2,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    
    # Traitement
    for(i in 1:length(fit)){
      print(paste0(i,"/",length(fit)))
      fit.i <- fit[[i]]
      for(j in 1:N){
        if(j%%50==0) print(paste0(j,"/",N))
        if(!is.na(fit.i[j,"shape"])){
          shape <- fit.i[j,"shape"]
          scale <- fit.i[j,"scale"] 
          p0 <- fit.i[j,"p0"]
          
          if(!rgamma){
            p <- runif(n = n,min = 0,max = q.threshold)
            q <- rep(0,length(p)) # on initialise le vecteur avec precip = 0
            ind <- which(p>p0)
            ppos <- ((p-p0)/(1-p0))[ind]
            pq <- qgamma(p=ppos,shape = shape,scale = scale)
            q[ind] <- pq
          }else{
            pq <- rgamma(n = n*(1-p0),shape = shape,scale = scale)
            q <- c(pq,rep(0,n-length(pq)))
          }
          precip.ana.sample[[i]][[j]] <- q
        }
      }
    }
  }
  
  # Sortie
  save(precip.ana.sample,file = get.path(generate_precip_all = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
}

# Genere les precip par analogie par tirage aleatoire dans les analogues ou dans les loi
generate.precip.random <- function(n=45,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=1,seuil=0,rgamma=F,moments=F,season=T,ncores=3){
  
  # On se limite a 45 series generees car 45 analogues uniquement pour nbana=0.2
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k=1)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2010-12-31")
  
    
  # Import des parametres de la loi gamma
  print(get.path(fit_gamma_random = T,rean = rean,period = period,moments = moments,bv = bv,season = season,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
  load(get.path(fit_gamma_random = T,rean = rean,period = period,moments = moments,bv = bv,season = season,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
  replic <- length(fit) 
  
    # Traitement
  precip.ana.sample <- matrix(NA,nrow = N,ncol = n)
  
  print(paste0("Parallelisation sur ",ncores, " coeurs"))
  outfile <- paste0(get.dirstr(k,rean,period),"generate.precip.random/calcul.txt")
  print(paste0("Logfile for // loop : ",outfile))
  cl <- makeCluster(ncores, outfile=outfile) 
  registerDoParallel(cl)
  
  precip.random <- foreach (i=1:replic) %dopar%{
    
    print(paste0(i,"/",replic))
    fit.i <- fit[[i]]
    
    for(j in 1:N){
      
      shape <- fit.i[j,"shape"]
      scale <- fit.i[j,"scale"] 
      p0 <- fit.i[j,"p0"]
      
      if(!rgamma){
        p <- runif(n = n,min = 0,max = q.threshold)
        q <- rep(0,length(p)) # on initialise le vecteur avec precip = 0
        ind <- which(p>p0)
        ppos <- ((p-p0)/(1-p0))[ind]
        pq <- qgamma(p=ppos,shape = shape,scale = scale)
        q[ind] <- pq
      }else{
        pq <- rgamma(n = n*(1-p0),shape = shape,scale = scale)
        q <- c(pq,rep(0,n-length(pq)))
      }
      precip.ana.sample[j,] <- q
    }
    precip.ana.sample
  }
  
  stopCluster(cl)
  
  # Sortie
  precip.ana.sample <- precip.random
  save(precip.ana.sample,file = get.path(generate_precip_random = T,rean = rean,period = period,moments = moments,bv = bv,season = season,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Renvoie les indices associes a une saison pour l'analogie saisoniere
get.ind.season.ana <- function(dat,dates){
  
  if(substr(dat,6,10)=="02-29") dat <- paste0(substr(dat,1,5),"02-28") # si 29 fevrier, on centre sur la veille
  pos <- which(substr(dates,6,10) == substr(dat,6,10)) # tous les indices de ce jour calendaire
  pos <- as.list(pos)
  pos <- unlist(lapply(pos,function(v){(v-30):(v+30)}))
  pos <- pos[pos>0 & pos <= length(dates)]
  pos
}

# Renvoie les limites en lon et lat pour un bv pour une variable thermodynamique
get.lon.lat.var.bv <- function(bv){
  
  if(bv=="Isere-seul"){
    lim.lon <- c(4,7)# c(6.25,6.5) #c(4.5,7)
    lim.lat <- c(43,46.5) #c(45.25,45.5) #c(44.5,46.5)
  }
  
  if(bv=="Drac-seul"){
    lim.lon <- c(4,7)#c(6,6.25) #c(4.5,7)
    lim.lat <- c(43,46.5)#lim.lat <- c(44.75,45) #c(44.5,46.5)
  }
  
  list(lim.lon=lim.lon,lim.lat=lim.lat)
}

# Renvoie le chemin vers un fichier
get.path <- function(k=NULL,rean=NULL,period=NULL,season=NULL,dP=NULL,dist=NULL,nbdays=NULL,start.rean=NULL,end.rean=NULL,start.ana=NULL,end.ana=NULL,bv=NULL,spazm=NULL,nbana=NULL,seuil=NULL,moments=NULL,type=NULL,rgamma=NULL,q.threshold=NULL,schaake=NULL,random=NULL,short=NULL,gamma=F,exp=F,invgauss,start=NULL,end=NULL,var2=NULL,nbana2=NULL,
                     save_ana=F,save_precip_ana=F,save_precip_ana_all=F,fit_gamma=F,fit_gamma_all=F,fit_exp=F,fit_gamma_random=F,fit_invgauss=F,generate_precip=F,generate_precip_all=F,generate_precip_random=F,reordering_precip=F,plot_fit_precip=F,plot_param_ana=F,plot_precip_ana=F,plot_precip_ana_day=F,plot_cum_ana_sais=F,plot_max_ana_sais=F,plot_distrib_ana_sais=F,fit_gev_ana=F,compute_gev_ana=F,plot_gev=F){
  
  if(save_ana){path <- paste0(get.dirstr(k,rean,period),"save.ana/ana_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),dist,"_",rean,"_k",k,ifelse(!is.null(var2),paste0("_",var2),""),"_mean",nbdays,"day_",start.rean,"_",end.rean,"_ana_",ifelse(short,start.ana,start.rean),"_",ifelse(short,end.ana,end.rean),".Rdata")}
  if(save_precip_ana){path <- paste0(get.dirstr(k,rean,period),"save.precip.ana/precip_ana_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(save_precip_ana_all){path <- paste0(get.dirstr(k,rean,period),"save.precip.ana/precip_ana_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),"allbv_","mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(fit_gamma){path <- paste0(get.dirstr(k,rean,period),"fit.gamma/fit_gamma_",ifelse(moments,"moments_",""),bv,"_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_seuil",seuil,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(fit_gamma_all){path <- paste0(get.dirstr(k,rean,period),"fit.gamma/fit_gamma_",ifelse(moments,"moments_",""),"allbv_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),"mean",nbdays,"day_nbana",nbana,"_seuil",seuil,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(fit_exp){path <- paste0(get.dirstr(k,rean,period),"fit.exp/fit_exp_",ifelse(moments,"moments_",""),bv,"_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_seuil",seuil,"_",rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(fit_gamma_random){path <- paste0(get.dirstr(k,rean,period),"fit.gamma.random/fit_gamma_",ifelse(moments,"moments_",""),bv,ifelse(season,"_season",""),"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_seuil",seuil,"_",rean,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(fit_invgauss){path <- paste0(get.dirstr(k,rean,period),"fit.invgauss/fit_invgauss_",ifelse(moments,"moments_",""),bv,"_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_seuil",seuil,"_",rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(generate_precip){path <- paste0(get.dirstr(k,rean,period),"generate.precip/precip_",type,ifelse(type!="empir" & moments,"_moments",""),"_",ifelse(type!="empir" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type!="empir",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(generate_precip_all){path <- paste0(get.dirstr(k,rean,period),"generate.precip/precip_",type,ifelse(type!="empir" & moments,"_moments",""),"_",ifelse(type!="empir" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type!="empir",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),"allbv_","mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(generate_precip_random){path <- paste0(get.dirstr(k,rean,period),"generate.precip.random/precip_",ifelse(moments,"moments_",""),bv,ifelse(season,"_season",""),"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_seuil",seuil,"_",rean,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(reordering_precip){path <- paste0(get.dirstr(k,rean,period),"reordering.precip/precip_final_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".Rdata")}
  if(plot_fit_precip){path <- paste0(get.dirstr(k,rean,period),"plot.fit.precip/plot_",bv,"_",ifelse(gamma,"gamma_",""),ifelse(exp,"exp_",""),ifelse(invgauss,"invgauss_",""),ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",ifelse(gamma|exp,paste0("seuil",seuil,"_"),""),ifelse((gamma|exp) & moments,"moments_",""),rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".png")}
  if(plot_param_ana){path <- paste0(get.dirstr(k,rean,period),"plot.param.ana/plot_",bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(type=="gamma" & moments,"moments_",""),rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,".png")}
  if(plot_precip_ana){path <- paste0(get.dirstr(k,rean,period),"plot.precip.ana/precip_ana_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".png")}
  if(plot_precip_ana_day){path <- paste0(get.dirstr(k,rean,period),"plot.precip.ana.day/precip_ana_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma",paste0("q",q.threshold,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".png")}
  if(plot_cum_ana_sais){path <- paste0(get.dirstr(k,rean,period),"plot.cum.ana.sais/plot_cum_ana_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".png")}
  if(plot_max_ana_sais){path <- paste0(get.dirstr(k,rean,period),"plot.max.ana.sais/plot_max_ana_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".png")}
  if(plot_distrib_ana_sais){path <- paste0(get.dirstr(k,rean,period),"plot.distrib.ana.sais/plot_distrib_ana_",ifelse(random,"random_",""),type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,"_",start.rean,"_",end.rean,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".png")}
  if(fit_gev_ana){path <- paste0(get.dirstr(k,rean,period),"fit.gev.ana/fit_gev_ana_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start,"_",end,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".Rdata")}
  if(compute_gev_ana){path <- paste0(get.dirstr(k,rean,period),"compute.gev.ana/compute_gev_ana_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start,"_",end,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".Rdata")}
  if(plot_gev){path <- paste0(get.dirstr(k,rean,period),"plot.gev/plot_gev_",type,ifelse(type=="gamma" & moments,"_moments",""),"_",ifelse(type=="gamma" & !rgamma,paste0("q",q.threshold,"_"),""),ifelse(rgamma,"rgamma_",""),ifelse(type=="gamma",paste0("seuil",seuil,"_"),""),ifelse(season,"season_",""),ifelse(dP,"dP_",""),bv,"_",ifelse(spazm,"spazm_",""),"mean",nbdays,"day_nbana",nbana,"_",rean,"_k",k,"_",dist,ifelse(!is.null(var2),paste0("_nbana",nbana2,"_",rean,"_",var2,"_RMSE"),""),"_",start,"_",end,"_ana_",start.ana,"_",end.ana,ifelse(!schaake,"_without_schaake",""),".png")}
  path
}

# A voir si on fait aussi l'etude pour les cumuls saisonniers/cumuls saisonniers par WP
map.trend.cum.obs <- function(){
  
}

# Carte des tendances gev sur tous les bvs
map.trend.gev.obs <- function(type="rl20",nbdays=1,start="1950-01-01",end="2017-12-31",save=F){
  
  # Import des tendances gev
  load(paste0("2_Travail/1_Past/Rresults/compute.gev.obs/compute_gev_all_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  # Traitement pour vecteur trend RL20 ou trend mean
  sais <- c("winter","spring","summer","autumn")
  res <- matrix(data = NA,nrow = length(val.gev),ncol = length(sais))
  colnames(res) <- sais
  colu <- ifelse(type=="rl20",2,3)
  
  for(i in 1:length(sais)){
    res[,i] <- unname(unlist(lapply(val.gev,function(v){x <- v[[i]];(x[2,colu]-x[1,colu])/x[1,colu]*100})))
  }
  
  # Cartes
  if(save){
    png(filename = paste0("2_Travail/1_Past/Rresults/map.trend.gev.obs/map_trend_gev_",type,"_mean",nbdays,"day_",start,"_",end,".png"),width = 6,height = 7,units = "in",res = 600)
    par(mfrow=c(2,2),mar=c(1,1,2,1))
  }
  for(i in 1:length(sais)){
    myimage.plot.bv(vec = res[,i],type = "trend",main=nam2str(sais[i]))
  }
  if(save) graphics.off()
}

# Carte des BV moyens du sud-est de la France, colories par un vecteur (moyenne de precip, tendance d'extreme...)
myimage.plot.bv<-function(vec,type="trend",...){
  
  # Import des donnees carto
  load(file="2_Travail/Data/Carto/SecteurHydro/SecteurHydro_LII.Rdata")#bv.list
  load(file="2_Travail/Data/Carto/SecteurHydro/idBV_alps_SecteurHydro.Rdata") #idx.list
  id2plot<-which(unlist(lapply(idx.list,length))>0)
  
  # Parametres graphiques
  if(type=="trend"){
    leg <- seq(-100,100,20)
    colo<-colorRampPalette(c("darkred","white","darkblue"))(10)
    #colo <- brewer.pal(n = 10, name = "RdBu")
  }
  if(type=="mean"){
    ran <- range(vec)
    leg <- seq(ran[1],ran[2],length.out = 11)
    colo<-colorRampPalette(c("white", "darkblue"))(11)[-1]
  }
  
  breaks <- seq(leg[1],leg[length(leg)],length=length(colo)+1)
  thiscol.fct<-function(val){
    if (is.na(val)) tc<-NA
    else if (val>=breaks[length(breaks)]) tc<-colo[length(colo)]
    else tc<-colo[which(val>=breaks[-length(breaks)] & val<breaks[-1])]
    tc
  }
  
  thiscol<-sapply(vec,thiscol.fct)

  # Graphique
  plot(c(770.5,1070.5),c(1790.5,2165.5),type="n",xlab="",ylab="",axes=FALSE,...)
  for (i in id2plot) {
    for(j in 1:length(bv.list[[i]])){
      polygon(bv.list[[i]][[j]],col=thiscol[match(i,id2plot)])
      #text(x = mean(bv.list[[i]][[j]][,1]),y = mean(bv.list[[i]][[j]][,2]),labels=which(id2plot==i))
    }
  }
  box()
  
  # Legende
  colorlegend(colbar = colo,labels = paste0("     ",round(leg,0)),at = seq(0,1,length.out = length(leg)),
              vertical = T,xlim = c(1040,1060),ylim = c(1790,2165))
}

# Trace les cumuls saisonniers analogie/obs en fonction du temps
plot.cum.ana.sais <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,start,end){
  
  # Dates utiles
  dates <- getdates(start,end)
  
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  
  if(schaake){
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
  }else{
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  
  pos <- match(dates,dates.rean)
  precip.ana.final <- precip.ana.final[pos,] # on ne traite qu'avec la serie reduite
  
  # Traitement et graphique par saison
  sais <- c("winter","spring","summer","autumn","year")
  
  png(filename = get.path(plot_cum_ana_sais = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,start = start,end = end),width = 9,height = 9,units = "in",res = 600)
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  
  for(i in 1:length(sais)){
    
    print(sais[i])
    year <- unique(substr(dates,1,4))
    if(sais[i] == "winter") year <- year[-1] # car pas de premiere saison hiver
    
    # Traitement precip obs
    pos <- get.ind.season(sais = sais[i],start = start,end = end)
    precip.i <- precip[pos$pos.season]
    precip.i[pos$pos.NA] <- NA
    dim(precip.i) <- c(pos$l.season,pos$n.season)
    precip.i <- apply(precip.i,2,sum,na.rm=T)
    
    # Traitement precip analogues
    if(sais[i]!="year"){precip.ana.i <- precip.ana.final[pos$pos.season,]
    }else{precip.ana.i <- apply(precip.ana.final,2,function(v){v[pos$pos.season]})}
    precip.ana.i[pos$pos.NA,] <- NA
    dim(precip.ana.i) <- c(pos$l.season,pos$n.season*ncol(precip.ana.i))
    precip.ana.i <- apply(precip.ana.i,2,sum,na.rm=T)
    dim(precip.ana.i) <- c(length(year),ncol(precip.ana.final))
    
    mini <- apply(precip.ana.i,1,min)
    maxi <- apply(precip.ana.i,1,max)
    q10 <- apply(precip.ana.i,1,quantile,probs=0.1)
    q90 <- apply(precip.ana.i,1,quantile,probs=0.9)
    q25 <- apply(precip.ana.i,1,quantile,probs=0.25)
    q75 <- apply(precip.ana.i,1,quantile,probs=0.75)
    
    # Graphique
    ylim <- c(0,ifelse(sais[i]!="year",1000,2500))
    
    plot(year,precip.i,type="n",ylim=ylim,xlab="Year",ylab="Precipitation (mm)",main=paste0(nam2str(sais[i])," - ",nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
    grid()
    polygon(c(year,rev(year)), c(maxi, rev(mini)),
            col = "gray80",border = F)
    polygon(c(year,rev(year)), c(q10, rev(q90)),
            col = "gray50",border = F)
    polygon(c(year,rev(year)), c(q75, rev(q25)),
            col = "gray30",border = F)
    lines(year,precip.i,col="red",lwd=3)
    #legend("topleft",inset=.02,bty="n",c("Obs","Analog Scenarios (range,IC80%,IC50%)"),col=c("red","gray50"),lwd=2,cex=0.8,horiz = T)
    
  }
  graphics.off()
}

# Trace les distributions par saison des precip analogie/obs
plot.distrib.ana.sais <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,random=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  nyear <- length(unique(substr(dates.ana,1,4)))
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  gc()
  if(schaake & !random){
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
  }
  if(!schaake & !random){
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  if(random){
    print(get.path(generate_precip_random = T,rean = rean,period = period,moments = moments,bv = bv,season = season,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    load(get.path(generate_precip_random = T,rean = rean,period = period,moments = moments,bv = bv,season = season,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    precip.ana.final <- do.call(cbind,precip.ana.sample)
    #precip.ana.final <- precip.ana.sample[[1]]
  }
  gc()
  
  pos <- match(dates.ana,dates.rean)
  precip.ana.final <- precip.ana.final[pos,] # on ne traite qu'avec la serie reduite a la periode des obs
  
  # Traitement et graphique par saison
  sais <- c("winter","spring","summer","autumn","year")
  
  png(filename = get.path(plot_distrib_ana_sais = T,k = k,rean = rean,period = period,random = random,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake),width = 7,height = 9,units = "in",res = 600)
  par(mfrow=c(3,2),mar=c(4,4,2,2),pty="s")
  
  for(i in 1:length(sais)){
    
    print(sais[i])
    pos <- get.ind.season.past(sais = sais[i],start = start.end.ana[1],end = start.end.ana[2])
    rp.sort <- 1/(1-ppoints(length(pos)))
    rp.sort <- rp.sort*nyear/length(pos) # passage en journalier: on divise par le nombre de jour par an
    
    # Traitement precip obs
    precip.i <- sort(precip[pos])
    
    # Traitement precip analogues
    precip.ana.i <- apply(precip.ana.final[pos,],2,sort)
    mini <- apply(precip.ana.i,1,min)
    maxi <- apply(precip.ana.i,1,max)
    q10 <- apply(precip.ana.i,1,quantile,probs=0.1)
    q90 <- apply(precip.ana.i,1,quantile,probs=0.9)
    q25 <- apply(precip.ana.i,1,quantile,probs=0.25)
    q75 <- apply(precip.ana.i,1,quantile,probs=0.75)
    
    # Graphique
    rp.axis <- c(1,5,25,100)
    #ylim <- c(0,max(maxi))
    ylim <- c(0,200)
    plot(rp.sort,precip.i,log="x",xaxt="n",type="n",ylim=ylim,xlab="Return Period (year)",ylab=paste0(nbdays,"-day precipitation (mm/day)"),main=paste0(nam2str(sais[i])," - ",nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
    grid(nx = NA,ny = NULL);abline(v=rp.axis,lty=3,col="grey")
    axis(side = 1,at = rp.axis,labels = rp.axis)
    polygon(c(rp.sort,rev(rp.sort)), c(maxi,rev(mini)),col = "gray80",border = F)
    polygon(c(rp.sort,rev(rp.sort)), c(q90,rev(q10)),col = "gray50",border = F)
    polygon(c(rp.sort,rev(rp.sort)), c(q75,rev(q25)),col = "gray30",border = F)
    lines(rp.sort,precip.i,col="red",lwd=3)
    legend("topleft",bty="n",c("Obs","Analog Scenarios \n(range,IC80%,IC50%)"),col=c("red","gray50"),lwd=2,cex=0.7)
  }
  graphics.off()
}

# Trace les cdf empiriques, gamma, et egpd pour 9 jours aleatoires
plot.fit.precip <- function(k,dist,bv,spazm,nbdays,nbana,rean,period="past",gamma=F,exp=F,invgauss=F,season=T,dP=F,seuil=0,moments=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],start.end.rean[2])
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2010-12-31")
  
  # Import des precip et des lois
  load(get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana))
  
  if(gamma){
  load(get.path(fit_gamma = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
  fit_gamma <- fit
  }
  
  if(exp){
    load(get.path(fit_exp = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    fit_exp <- fit
  }
  
  if(invgauss){
    load(get.path(fit_invgauss = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
    fit_invgauss <- fit
  }
  
  # Graphique
  ind <- which(!apply(fit_gamma,1,function(v){all(is.na(v))}))
  print(paste0("Jours non fittés: ",round((N-length(ind))/N*100,1),"%"))
  
  #ind <- sample(ind,9);print(ind)
  ind <- c(43700,43637,11325,40762,39223,7796,53852,50003,43533)
  
  if(nbdays==1){
    main <- dates.rean[ind]
  }else{
    main <- paste0(dates.rean[ind]," - ",dates.rean[ind+nbdays-1])
  }

  png(filename = get.path(plot_fit_precip = T,k = k,dist = dist,rean = rean,period = period,bv = bv,gamma = gamma,exp = exp,invgauss = invgauss,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,moments = moments,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]),
      width = 6,height = 7,units = "in",res = 600)
  par(mfrow=c(3,3),mar=c(4,4,2,2))
  
  for(i in 1:length(ind)){
    
    # Empirique
    precip <- precip.ana[[ind[i]]]
    xlim <- c(0-diff(range(precip))*0.2,max(precip)+diff(range(precip))*0.4)
    plot(ecdf(precip),main=main[i],xlim=xlim,xlab=paste0(nbdays,"-day ",nam2str(bv)," Precipitation (mm/day)"))
    abline(v = 0,lty=3)
    
    if(gamma){
      # Gamma
      shape <- fit_gamma[ind[i],"shape"]
      scale <- fit_gamma[ind[i],"scale"] 
      p0 <- fit_gamma[ind[i],"p0"]
      p <- seq(round(p0,4)+0.0001,1,0.0001) # on veut tracer la loi de y=p0 à y=1, tous les 0.0001
      ppos <- (p-p0)/(1-p0)
      pq = qgamma(p=ppos,shape = shape,scale = scale) # on extrait les quantiles des probas positives (car gamma loi calee sur precip/quantiles positifs)
      lines(pq,p,col="red")
      segments(0,0,0,p0,col="red")
      abline(h = p0,lty=3)
    }
    
    if(exp){
      # Exponential
      rate <- fit_exp[ind[i],"rate"]
      p0 <- fit_exp[ind[i],"p0"]
      p <- seq(round(p0,4)+0.0001,1,0.0001) # on veut tracer la loi de y=p0 à y=1, tous les 0.0001
      ppos <- (p-p0)/(1-p0)
      pq = qexp(p=ppos,rate=rate) # on extrait les quantiles des probas positives
      lines(pq,p,col="royalblue")
      segments(0,0,0,p0,col="royalblue")
      abline(h = p0,lty=3)
    }
    
    if(invgauss){
      # Inverse Gaussian
      nu <- fit_invgauss[ind[i],"mean"]
      shape <- fit_invgauss[ind[i],"shape"]
      p0 <- fit_invgauss[ind[i],"p0"]
      p <- seq(round(p0,3)+0.001,0.999,0.001) # bizarre: erreur quand proba de 1 demandee
      ppos <- (p-p0)/(1-p0)
      pq = qinvGauss(p=ppos,nu=nu,lambda = shape) # on extrait les quantiles des probas positives
      lines(pq,p,col="darkgreen")
      segments(0,0,0,p0,col="darkgreen")
      abline(h = p0,lty=3)
    }
    
      # Valeurs tirees
      #load("I:/ongoing_Documents/2_Travail/1_Past/20CR-m1/generate.precip/precip_gamma_Isere-seul_spazm_mean1day_nbana_0.2_20CR-m1_k1_TWS_1851-01-01_2014-12-31_ana_1950-01-01_2010-12-31.Rdata")
      #precip.sample <- precip.ana.sample[[ind[i]]]
      #points(precip.sample,ecdf(precip.sample)(precip.sample),col="blue")
      #mean(precip.sample)
      #mean(precip)
  }
  graphics.off()
}

# Trace les RL20 et mean max selon la GEV pour les scenarios analogues et l'observe
plot.gev <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,start,end,var2=NULL,nbana2=NULL){
  
  # Import des GEV observees
  load(file = paste0("2_Travail/1_Past/Rresults/compute.gev.obs/compute_gev_",bv,"_mean",nbdays,"day_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
  gev.obs <- val.gev
  
  # Import des GEV analogues
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  load(file = get.path(compute_gev_ana = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,var2 = var2,nbana2 = nbana2,start = start,end = end))
  gev.ana <- val.gev
  rm(val.gev)
  
  # Graphiques
  sais <- c("winter","spring","summer","autumn")
  png(filename = get.path(plot_gev = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,var2 = var2,nbana2 = nbana2,start = start,end = end),width = 10,height = 6,units = "in",res = 600)
  par(mfcol=c(2,4),mar=c(2,4,2,1))
  
  if(bv=="Isere") {ylim1 <- c(0,120);ylim2 <- c(0,60)}
  if(bv=="Isere-seul") {ylim1 <- c(0,120);ylim2 <- c(0,60)}
  if(bv=="Drac-seul") {ylim1 <- c(0,160);ylim2 <- c(0,80)}
  
  for(i in 1:length(sais)){
    
    # RL 20ans
    plot(gev.obs[[i]][,1],gev.obs[[i]][,2],ylim=ylim1,type="n",xaxt="n",ylab="Precipitation (mm)",main=paste0("RL20 - ",nam2str(sais[i])))
    axis(side = 1,at = gev.obs[[i]][,1],labels = gev.obs[[i]][,1])
    grid(nx=NA,ny=NULL)
    abline(v = gev.obs[[i]][,1],col="grey",lty=3)
    for(j in 1:length(gev.ana)){
      lines(gev.ana[[j]][[i]][,1],gev.ana[[j]][[i]][,2],lty=ifelse(gev.ana[[j]][[i]][1,4]<0.05,1,2))
    }
    lines(gev.obs[[i]][,1],gev.obs[[i]][,2],col="red",lwd=2,lty=ifelse(gev.obs[[i]][1,4]<0.05,1,2))
    
    # Mean
    plot(gev.obs[[i]][,1],gev.obs[[i]][,3],ylim=ylim2,type="n",xaxt="n",ylab="Precipitation (mm)",main=paste0("Mean - ",nam2str(sais[i])))
    axis(side = 1,at = gev.obs[[i]][,1],labels = gev.obs[[i]][,1])
    grid(nx=NA,ny=NULL)
    abline(v = gev.obs[[i]][,1],col="grey",lty=3)
    for(j in 1:length(gev.ana)){
      lines(gev.ana[[j]][[i]][,1],gev.ana[[j]][[i]][,3],lty=ifelse(gev.ana[[j]][[i]][1,4]<0.05,1,2))
    }
    lines(gev.obs[[i]][,1],gev.obs[[i]][,3],col="red",lwd=2,lty=ifelse(gev.obs[[i]][1,4]<0.05,1,2))
  }
  graphics.off()
}

# Trace les cumuls saisonniers analogie/obs en fonction du temps
plot.max.ana.sais <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F,start,end){
  
  # Dates utiles
  dates <- getdates(start,end)
  
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  
  if(schaake){
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
  }else{
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  
  pos <- match(dates,dates.rean)
  precip.ana.final <- precip.ana.final[pos,] # on ne traite qu'avec la serie reduite
  
  # Traitement et graphique par saison
  sais <- c("winter","spring","summer","autumn","year")
  
  png(filename = get.path(plot_max_ana_sais = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake,start = start,end = end),width = 9,height = 9,units = "in",res = 600)
  par(mfrow=c(3,2),mar=c(4,4,2,2))
  
  for(i in 1:length(sais)){
    
    print(sais[i])
    year <- unique(substr(dates,1,4))
    if(sais[i] == "winter") year <- year[-1] # car pas de premiere saison hiver
    
    # Traitement precip obs
    pos <- get.ind.season(sais = sais[i],start = start,end = end)
    precip.i <- precip[pos$pos.season]
    precip.i[pos$pos.NA] <- NA
    dim(precip.i) <- c(pos$l.season,pos$n.season)
    precip.i <- apply(precip.i,2,max,na.rm=T)
    
    # Traitement precip analogues
    if(sais[i]!="year"){precip.ana.i <- precip.ana.final[pos$pos.season,]
    }else{precip.ana.i <- apply(precip.ana.final,2,function(v){v[pos$pos.season]})}
    precip.ana.i[pos$pos.NA,] <- NA
    dim(precip.ana.i) <- c(pos$l.season,pos$n.season*ncol(precip.ana.i))
    precip.ana.i <- apply(precip.ana.i,2,max,na.rm=T)
    dim(precip.ana.i) <- c(length(year),ncol(precip.ana.final))
    
    mini <- apply(precip.ana.i,1,min)
    maxi <- apply(precip.ana.i,1,max)
    q10 <- apply(precip.ana.i,1,quantile,probs=0.1)
    q90 <- apply(precip.ana.i,1,quantile,probs=0.9)
    q25 <- apply(precip.ana.i,1,quantile,probs=0.25)
    q75 <- apply(precip.ana.i,1,quantile,probs=0.75)
    
    # Graphique
    ylim <- c(0,max(maxi))
    
    plot(year,precip.i,type="n",ylim=ylim,xlab="Year",ylab="Precipitation (mm)",main=paste0(nam2str(sais[i])," - ",nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
    grid()
    polygon(c(year,rev(year)), c(maxi, rev(mini)),
            col = "gray80",border = F)
    polygon(c(year,rev(year)), c(q10, rev(q90)),
            col = "gray50",border = F)
    polygon(c(year,rev(year)), c(q75, rev(q25)),
            col = "gray30",border = F)
    lines(year,precip.i,col="red",lwd=3)
    #legend("topleft",inset=.02,bty="n",c("Obs","Analog Scenarios (range,IC80%,IC50%)"),col=c("red","gray50"),lwd=2,cex=0.8,horiz = T)
    
  }
  graphics.off()
}

# Trace les densites des parametres de la gamma par bv pour chaque saison
plot.param.ana <- function(k,dist,bv,spazm,nbdays,nbana,rean,period="past",type="empir",season=T,dP=F,seuil=0,moments=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],start.end.rean[2])
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2010-12-31")
  
  # Import des lois
  load(get.path(fit_gamma = T,k = k,rean = rean,period = period,moments = moments,bv = bv,season = season,dP = dP,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,dist = dist,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
  
  # Graphiques gamma classique
  png(filename = get.path(plot_param_ana = T,k = k,dist = dist,rean = rean,period = period,bv = bv,type = type,spazm = spazm,nbdays = nbdays,nbana = nbana,seuil = seuil,moments = moments,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]),width = 7,height = 6,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(4,4,3,0))
  
  # Shape
  plot(density(rgamma(n = 1000000,shape = 1,scale = 1)),type="l",col=1,main="Gamma",xlab="Precipitation (mm)")
  for(i in 2:5){
    lines(density(rgamma(n = 1000000,shape = i,scale = 1)),type="l",col=i)
  }
  legend("topright",inset=.02,bty="n",legend = paste0("shape=",1:5,", scale=1"),col=1:5,lty=1)
  
  # Scale
  plot(density(rgamma(n = 1000000,shape = 1,scale = 1)),type="l",col=1,main="Gamma",xlab="Precipitation (mm)")
  for(i in 2:5){
    lines(density(rgamma(n = 1000000,shape = 1,scale = i)),type="l",col=i)
  }
  legend("topright",inset=.02,bty="n",legend = paste0("shape=1, scale=",1:5),col=1:5,lty=1)
  
  
  # Graphiques gamma precip bv
  season <- c("winter","spring","summer","autumn")
  
  # Shape
  ind <- get.ind.season.past(sais = season[3],start = start.end.rean[1],end = start.end.rean[2]) # on commence par ete pour ylim
  plot(density(na.omit(fit[ind,"shape"])),ylim=c(0,1.2*max(density(na.omit(fit[ind,"shape"]))$y)),col=1,main=paste0(nam2str(bv)," - Shape"),xlab="Shape")
  for(i in 1:length(season)){
    ind <- get.ind.season.past(sais = season[i],start = start.end.rean[1],end = start.end.rean[2])
    lines(density(na.omit(fit[ind,"shape"])),col=i)
  }
  legend("topright",inset=.02,bty="n",col=1:4,legend = nam2str(season),lty=1)
  
  # Scale
  ind <- get.ind.season.past(sais = season[3],start = start.end.rean[1],end = start.end.rean[2])
  plot(density(na.omit(fit[ind,"scale"])),ylim=c(0,1.2*max(density(na.omit(fit[ind,"scale"]))$y)),col=1,main=paste0(nam2str(bv)," - Scale"),xlab="Scale")
  for(i in 1:length(season)){
    ind <- get.ind.season.past(sais = season[i],start = start.end.rean[1],end = start.end.rean[2])
    lines(density(na.omit(fit[ind,"scale"])),col=i)
  }
  legend("topright",inset=.02,bty="n",col=1:4,legend = nam2str(season),lty=1)
    
  graphics.off()
}

# Trace differents graphiques des precip analogues VS obs
plot.precip.ana <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,schaake=T,season=F,dP=F,rgamma=T,moments=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  if(schaake){
    print(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
  }else{
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  
  pos <- match(dates.ana,dates.rean)
  precip.ana.final <- precip.ana.final[pos,] # on ne traite qu'avec la serie reduite
  
  # Traitement cumul annuel
  precip.ann <- aggregate(precip,by=list(substr(dates.ana,1,4)),sum,na.rm=T)
  precip.ana.ann <- aggregate(precip.ana.final,by=list(substr(dates.ana,1,4)),sum,na.rm=T)
  mini <- apply(precip.ana.ann[,-1],1,min)
  maxi <- apply(precip.ana.ann[,-1],1,max)
  q10 <- apply(precip.ana.ann[,-1],1,quantile,probs=0.1)
  q90 <- apply(precip.ana.ann[,-1],1,quantile,probs=0.9)
  q25 <- apply(precip.ana.ann[,-1],1,quantile,probs=0.25)
  q75 <- apply(precip.ana.ann[,-1],1,quantile,probs=0.75)
  
  # Graphique evolution cumul annuel
  png(filename = get.path(plot_precip_ana = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake),width = 9,height = 6,units = "in",res = 600)
  layout(matrix(1:4,2,2,byrow = T),widths = c(2,1),heights = c(1,1))
  par(mar=c(4,4,2,2))
  ylim <- c(0,2500)
  
  plot(precip.ana.ann[,1],mini,type="n",ylim=ylim,xlab="Year",ylab="Precipitation (mm)",main=paste0(nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
  grid()
  polygon(c(precip.ana.ann[,1], rev(precip.ana.ann[,1])), c(maxi, rev(mini)),
          col = "gray80",border = F)
  polygon(c(precip.ana.ann[,1], rev(precip.ana.ann[,1])), c(q10, rev(q90)),
          col = "gray50",border = F)
  polygon(c(precip.ana.ann[,1], rev(precip.ana.ann[,1])), c(q75, rev(q25)),
          col = "gray30",border = F)
  lines(precip.ann,col="red",lwd=3)
  legend("bottomleft",inset=.02,bty="n",c("Obs","Analog Scenarios \n(range,IC80%,IC50%)"),col=c("red","gray50"),lwd=2,cex=0.8)
  
  # Traitement cycle interannuel (periode commune 1950-2010)
  nyear <- length(unique(substr(dates.ana,1,4)))
  precip.inter <- aggregate(precip,by=list(substr(dates.ana,6,7)),sum)
  precip.inter[,2] <- precip.inter[,2]/nyear # on a le cumul moyen par mois
  
  precip.ana.inter <- aggregate(precip.ana.final,by=list(substr(dates.ana,6,7)),sum)
  precip.ana.inter[,-1] <- precip.ana.inter[,-1]/nyear # on a le cumul moyen par mois
  mini <- apply(precip.ana.inter[,-1],1,min)
  maxi <- apply(precip.ana.inter[,-1],1,max)
  q10 <- apply(precip.ana.inter[,-1],1,quantile,probs=0.1)
  q90 <- apply(precip.ana.inter[,-1],1,quantile,probs=0.9)
  q25 <- apply(precip.ana.inter[,-1],1,quantile,probs=0.25)
  q75 <- apply(precip.ana.inter[,-1],1,quantile,probs=0.75)
  
  # Graphique cycle interannuel
  ylim <- c(0,200)
  plot(precip.inter[,1],precip.inter[,2],type="n",ylim=ylim,xaxt="n",xlab="Month",ylab="Precipitation (mm)",main=paste0(nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
  axis(side = 1,at = 1:12,labels = month.abb)
  grid(nx = NA,ny = NULL);abline(v=1:12,lty=3,col="grey")
  polygon(c(precip.ana.inter[,1], rev(precip.ana.inter[,1])), c(maxi, rev(mini)),
          col = "gray80",border = F)
  polygon(c(precip.ana.inter[,1], rev(precip.ana.inter[,1])), c(q90, rev(q10)),
          col = "gray50",border = F)
  polygon(c(precip.ana.inter[,1], rev(precip.ana.inter[,1])), c(q75, rev(q25)),
          col = "gray30",border = F)
  lines(precip.inter[,1],precip.inter[,2],col="red",lwd=3)
  legend("bottomleft",inset=.02,bty="n",c("Obs","Analog Scenarios \n(range,IC80%,IC50%)"),col=c("red","gray50"),lwd=2,cex=0.8)
  
  # Traitement boxplots saisonniers (en base R car pas possible de combiner ggplot et base R)
  sea <- c("winter","spring","summer","autumn")
  precip.sais <- list()
  
  for(i in 1:length(sea)){
    
    ind <- get.ind.season(sais = sea[i],start =  start.end.ana[1],end =  start.end.ana[2])
    p <- precip[ind$pos.season]
    dim(p) <- c(ind$l.season,ind$n.season)
    p <- apply(p,2,sum)
    
    p.ana <- precip.ana.final[ind$pos.season,]
    dim(p.ana) <- c(ind$l.season,ind$n.season*ncol(p.ana))
    p.ana <- apply(p.ana,2,sum)
    
    precip.sais[[i*2-1]] <- p
    precip.sais[[i*2]] <- p.ana
  }
  
  l.max <- max(unlist(lapply(precip.sais,length)))
  precip.sais <- lapply(precip.sais,function(v){c(v,rep(NA,l.max-length(v)))})
  precip.sais <- do.call(cbind,precip.sais)
  
  # Graphique boxplots saisonniers
  colo <- rep(c("red","grey"),length(sea))
  nam <- c(t(cbind(paste0(nam2str(sea),"\n Obs"),paste0(nam2str(sea),"\n Ana"))))
  boxplot(precip.sais,col=colo,xaxt="n",ylab="Precipitation (mm)",main=paste0(nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
  axis(side = 1,at = 1:8,labels = nam,padj=0.5,cex.axis=0.8)
  
  # Traitement distribution des extremes
  #rp.sort <- sort(1/(1-ecdf(precip)(precip))/365.25) # periode de retour empirique
  #rp.sort[is.infinite(rp.sort)] <- nyear
  rp.sort <- 1/(1-ppoints(length(precip)))/365.25
  precip.sort <- sort(precip)
  
  precip.sort.ana <- apply(precip.ana.final,2,sort)
  mini <- apply(precip.sort.ana,1,min)
  maxi <- apply(precip.sort.ana,1,max)
  q10 <- apply(precip.sort.ana,1,quantile,probs=0.1)
  q90 <- apply(precip.sort.ana,1,quantile,probs=0.9)
  q25 <- apply(precip.sort.ana,1,quantile,probs=0.25)
  q75 <- apply(precip.sort.ana,1,quantile,probs=0.75)
  
  
  # Graphique distribution des extremes
  rp.axis <- c(1,5,25,100)
  ylim <- c(0,120)
  plot(rp.sort,precip.sort,log="x",xaxt="n",type="n",ylim=ylim,xlab="Return Period (year)",ylab=paste0(nbdays,"-day precipitation (mm/day)"),main=paste0(nam2str(bv)," - ",ncol(precip.ana.final)," scenarios - ",rean))
  grid(nx = NA,ny = NULL);abline(v=rp.axis,lty=3,col="grey")
  axis(side = 1,at = rp.axis,labels = rp.axis)
  polygon(c(rp.sort,rev(rp.sort)), c(maxi,rev(mini)),col = "gray80",border = F)
  polygon(c(rp.sort,rev(rp.sort)), c(q90,rev(q10)),col = "gray50",border = F)
  polygon(c(rp.sort,rev(rp.sort)), c(q75,rev(q25)),col = "gray30",border = F)
  lines(rp.sort,precip.sort,col="red",lwd=3)
  legend("topleft",bty="n",c("Obs","Analog Scenarios \n(range,IC80%,IC50%)"),col=c("red","gray50"),lwd=2,cex=0.7)
  graphics.off()
  
}

# Trace differents graphiques des precip analogues VS obs au pas de temps de construction (1j ou 3j)
plot.precip.ana.day <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,schaake=T,season=F,dP=F,moments=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des precipitations
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  if(schaake){
    load(get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
  }else{
    load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
    precip.ana.final <- do.call(rbind,precip.ana.sample)
  }
  
  pos <- match(dates.ana,dates.rean)
  precip.ana.final <- precip.ana.final[pos,]
  
  # Traitement obs
  precip.mean <- aggregate(precip,by=list(substr(dates.ana,6,10)),mean)[,2]
  precip.med <- aggregate(precip,by=list(substr(dates.ana,6,10)),median)[,2]
  precip.q90 <- aggregate(precip,by=list(substr(dates.ana,6,10)),quantile,probs=0.9)[,2]
  precip.max <- aggregate(precip,by=list(substr(dates.ana,6,10)),max)[,2]
  
  # Traitement ana: on regroupe tous les membres pour etre independant du scenario
  calday <- substr(getdates("2020-01-01","2020-12-31"),6,10)
  precip.ana <- list()
  for(i in 1:length(calday)){
    print(i)
    pos <- which(substr(dates.ana,6,10)==calday[i])
    precip.ana[[i]] <- c(precip.ana.final[pos,])
  }
  
  precip.ana.mean <- unlist(lapply(precip.ana,mean))
  precip.ana.med <- unlist(lapply(precip.ana,median))
  precip.ana.q90 <- unlist(lapply(precip.ana,quantile,probs=0.9))
  precip.ana.max <- unlist(lapply(precip.ana,max))
  
  # Graphiques
  png(filename = get.path(plot_precip_ana_day = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,schaake = schaake),width = 12,height = 6,units = "in",res = 600)
  par(mfrow=c(2,2))
  
  plot(precip.mean,type="l",col="red",lwd=3,xlab="Days",ylab="Precipitation (mm)",main="Mean",ylim=range(c(precip.mean,precip.ana.mean)))
  grid()
  lines(precip.mean,col="red",lwd=3)
  lines(precip.ana.mean,col="grey",lwd=2)
  
  plot(precip.med,type="l",col="red",lwd=3,xlab="Days",ylab="Precipitation (mm)",main="Median",ylim=range(c(precip.med,precip.ana.med)))
  grid()
  lines(precip.med,col="red",lwd=3)
  lines(precip.ana.med,col="grey",lwd=2)
  
  plot(precip.q90,type="l",col="red",lwd=3,xlab="Days",ylab="Precipitation (mm)",main="q90",ylim=range(c(precip.q90,precip.ana.q90)))
  grid()
  lines(precip.q90,col="red",lwd=3)
  lines(precip.ana.q90,col="grey",lwd=2)
  
  plot(precip.max,type="l",col="red",lwd=3,xlab="Days",ylab="Precipitation (mm)",main="Max",ylim=range(c(precip.max,precip.ana.max)))
  grid()
  lines(precip.max,col="red",lwd=3)
  lines(precip.ana.max,col="grey",lwd=2)
  
  graphics.off()
}

# Trace la relation entre la singularite et le nbre de fois ou le jour est choisi comme analogue
plot.sing.ana <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,schaake=T,season=F,dP=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist")
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Imports
  pos <- match(dates.ana,dates.rean)
  sing05 <- get.descriptor(descriptor = "sing05",k = k,dist = dist,nbdays = nbdays,start = start.end.rean[1],end = start.end.rean[2],standardize = F,rean = rean,threeday = F,
                           desais = F,period = "past",start.ana = start.end.ana[1],end.ana = start.end.ana[2])
  
  sing05 <- sing05[pos]
  load(paste0(get.dirstr(k,rean,period),"save.ana/ana_",ifelse(season,"season_",""),ifelse(dP,"dP_",""),dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  precip <- get.precip(nbdays,start.end.ana[1],start.end.ana[2],bv,spazm)
  
  # Traitement
  nbnei <- round(nbana*0.01*n,0)
  nei <- lapply(nei,function(v){v[1:nbnei]})
  nei <- unlist(nei)
  co <- rep(0,n)
  count <- table(nei)
  nam <- as.numeric(names(count))
  co[nam] <- count
  
  plot(sing05,co)
  plot(density(co))
}

# Application du Schaake shuffle temporel
reordering.precip <- function(type="empir",k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",q.threshold=0.99,seuil=0,season=F,dP=F,rgamma=F,moments=F,var2=NULL,nbana2=NULL){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k,var = var2)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  start.end.ana <- c("1950-01-01","2017-12-31")
  
  # Import precip
  print("Imports")
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  # Import precip generees
  load(get.path(generate_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
  N <- length(precip.ana.sample)
  
  # Import dates analogues
  load(get.path(save_ana = T,short = T,k = k,rean = rean,period = period,season = season,dP = F,dist = dist,nbdays = nbdays,var2 = var2,nbana2 = nbana2,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
  nbnei <- length(precip.ana.sample[[1]]) # on veut le meme nbre d'analogues que de precip a reordonner
  nei <- lapply(nei,function(v){v[1:nbnei]})
  
  # Schaake shuffle
  # rang.im1p1 <- lapply(nei,function(v){rank(precip[v+1],ties.method = "random")}) # attribution d'un rang aleatoire pour les 0
  # rang.brut <- lapply(precip.ana.sample,function(v){rank(v,ties.method = "random")})
  # new.order <- mapply(function(x,y){match(x,y)},x=rang.im1p1,y=rang.brut,SIMPLIFY = F)
  # precip.ana.final <- t(mapply(function(x,y){x[y]},x=precip.ana.final,y=new.order)) # directement les series de precip qu'on veut
  # precip.ana.final <- unname(rbind(precip.ana.sample.1,precip.ana.final))
  
  # Pas possible de faire sans boucle! Il faut reclasser selon l'ordre reclassé de i-1!
  print("Schaake shuffle")
  precip.ana.final <- precip.ana.sample
  nei.new <- nei
  
  for(i in 2:N){ # on ne reordonne pas le premier pas de temps
    
    if(i%%50==0) print(paste0(i,"/",N))
    
    # rang precip lendemain des analogues de la veille
    rang.im1p1 <- rank(precip[nei.new[[i-1]]+1],ties.method = "random")
    
    # rang precip du jour
    rang.brut <- rank(precip.ana.final[[i]],ties.method = "random")
    
    # nouvel ordre
    pos <- match(rang.im1p1,rang.brut)
    
    # reordonnancement et mise à jour nei
    precip.ana.final[[i]] <- precip.ana.final[[i]][pos]
    nei.new[[i]] <- nei.new[[i]][pos]
  }
  
  precip.ana.final <- do.call(rbind,precip.ana.final)
  
  # Export
  save(precip.ana.final,file = get.path(reordering_precip = T,k = k,rean = rean,period = period,type = type,rgamma=rgamma,q.threshold = q.threshold,seuil = seuil,moments = moments,season = season,dP = dP,bv = bv,spazm = spazm,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana,var2 = var2,nbana2 = nbana2))
}

# Run la methode des analogues
run.analogs <- function(k = 1,dist = "TWS",start = "1851-01-01",end = "2011-12-31",
                        rean = "20CR-m1",period = "past",nbdays = 1,nbana = 0.2,
                        bv = "Isere-seul",spazm = T,season = T,seuil = 0,moments = T,
                        type = "empir",q.threshold = 0.99,schaake = T,var2 = "vv700",nbana2 = 0.2){
  
  # Parametres
  compute_dist=F
  save_ana=F
  save_precip=F
  fit_law=F
  generate_precip=T
  gev=F
  graphics=F
  
  # Calcul des distances
  if(compute_dist){
    print("Calcul des distances")
    compute_dist_gen(k = k,dist = dist,start = start,end = end,rean = rean,period = period)
    if(k==3){combine.dist.k3(wk1 = 0.5,wk2 = 0.5,dist = dist,rean = rean,start = start,end = end)}
  }
  
  # Enregistrement des analogues et des precipitations associees
  if(save_ana){
    print("Enregistrement des dates analogues")
    #save.ana(k = k,dist = dist,nbdays = nbdays,rean = rean,period = period)
    save.ana.par(k = k,dist = dist,nbdays = nbdays,rean = rean,period = period,season = season,ncores = 6)
    #save.ana.sea(k = k,dist = dist,nbdays = nbdays,rean = rean,period = period)
  }
  if(save_precip){
    print("Enregistrement des precipitations analogues")
    #save.precip.ana(k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,season = season,dP = F)
    save.precip.ana.all(k = k,dist = dist,nbdays = nbdays,nbana = nbana,rean = rean,period = period,season = season,dP = F)
    #save.precip.ana.twostage(k = k,dist1 = dist,var = var2,nbdays = nbdays,nbana1 = nbana,nbana2 = nbana2,bv = bv,spazm = spazm,rean = rean,period = period,season = season,dP = F)
    }
  
  # Calage d'une loi parametrique
  if(fit_law){
    print("Calage d'une loi parametrique")
    #fit.gamma(k = k,dist = dist,nbdays = nbdays,nbana = nbana,nbmini = 10,seuil = seuil,bv = bv,spazm = spazm,rean = rean,period = period,season = season,dP = F,moments = moments,var2 = var2,nbana2 = nbana2)
    fit.gamma.all(k = k,dist = dist,nbdays = nbdays,nbana = nbana,nbmini = 10,seuil = seuil,rean = rean,period = period,season = season,dP = F,moments = moments,var2 = var2,nbana2 = nbana2)
    #fit.gamma.random(replic = 1000,nbdays = nbdays,nbana = nbana,nbmini = 10,seuil = seuil,bv = bv,spazm = spazm,rean = rean,period = period,moments = moments,ncores = 10)
    #fit.exp(k = k,dist = dist,nbdays = nbdays,nbana = nbana,nbmini = 10,seuil = seuil,bv = bv,spazm = spazm,rean = rean,period = period,season = season,dP = F,moments = T)
    #fit.invgauss(k = k,dist = dist,nbdays = nbdays,nbana = nbana,nbmini = 10,seuil = seuil,bv = bv,spazm = spazm,rean = rean,period = period,season = season,dP = F,moments = T)
    }
  
  # Generation des series de precipitations avec reodonnancement
  if(generate_precip){
    print("Generation des series de precipitations")
    #generate.precip(n = 50,type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,season = season,dP = F,rgamma = T,moments = moments,var2 = var2,nbana2 = nbana2)
    generate.precip.all(n = 50,type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,season = season,dP = F,rgamma = T,moments = moments,var2 = var2,nbana2 = nbana2)
    #generate.precip.random(n = 45,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = "past",q.threshold = q.threshold,seuil = seuil,rgamma = F,moments = moments,ncores = 10)
    #reordering.precip(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,season = season,dP = F,rgamma = T,moments = moments,var2 = var2,nbana2 = nbana2)
  }
  
  # Calculs GEV
  if(gev){
    #fit.gev.obs(bv = bv,nbdays = nbdays,spazm = spazm,start = start,end = end)
    fit.gev.ana.par(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = T,moments = moments,start = start,end = end,var2 = var2,nbana2 = nbana2)
    #compute.gev.obs(bv = bv,nbdays = nbdays,spazm = spazm,start = start,end = end)
    compute.gev.ana(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = T,moments = moments,start = start,end = end,var2 = var2,nbana2 = nbana2)
    }
  
  # Graphiques de sortie
  if(graphics){
    print("Graphiques")
    #plot.fit.precip(k = k,dist = dist,bv = bv,spazm = spazm,nbdays = nbdays,nbana = nbana,rean = rean,period = period,gamma = T,exp = T,invgauss = T,season = season,dP = F,seuil = seuil,moments = moments)
    #plot.param.ana(k = k,dist = dist,bv = bv,spazm = spazm,nbdays = nbdays,nbana = nbana,rean = rean,period = period,type = type,season = season,dP = F,seuil = seuil,moments = moments)
    #plot.precip.ana(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = F,moments = moments)
    #plot.precip.ana.day(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,schaake = schaake,season = season,dP = F,moments = moments)
    #plot.cum.ana.sais(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = T,moments = moments,start = start,end = end)
    #plot.max.ana.sais(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = T,moments = moments,start = start,end = end)
    #plot.distrib.ana.sais(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = F,moments = moments,random=F)
    #compare.schaake(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,season = season,moments = moments)
    plot.gev(type = type,k = k,dist = dist,nbdays = nbdays,nbana = nbana,bv = bv,spazm = spazm,rean = rean,period = period,q.threshold = q.threshold,seuil = seuil,schaake = schaake,season = season,dP = F,rgamma = T,moments = moments,start = start,end = end,var2 = var2,nbana2 = nbana2)
    }
}

# Run la methode des analogues en boucle
run.analogs.loop <- function(){
  
  # Attention a bien boucler uniquement sur les parametres utiles
  # Attention a bien selectionner uniquement les etapes utiles dans run.analogs
  # Attention a ne pas utiliser k dans les indices de boucle
  
  # Parametres
  k <- c(1,2,3)
  dist <- "TWS"
  start <- "1950-01-01"
  end <- "2017-12-31"
  rean <- c("ERA5")#,"20CR-m0","20CR-m2","ERA20C","ERA5")
  period <- "past"
  nbdays <- c(1)#,3)
  nbana <- c(0.2)
  nbana2 <- NULL
  var2 <- NULL
  bv <- "Isere-seul"#c("Drac-seul","Isere-seul","Isere")
  spazm <- c(T)#,F)
  season <- c(T)#,F)
  seuil <- c(0)#,0.1)
  moments <- c(T)#,F)
  type <- "gamma"#c("empir","invgauss")
  q.threshold <- c(1)#,0.999,1)
  schaake <- c(T)#,F)
  
  for(i in 1:length(k)){
    print(paste0("k=",k[i]))
    for(j in 1:length(rean)){
      print(paste0("rean=",rean[j]))
      for(l in 1:length(nbdays)){
        print(paste0("nbdays=",nbdays[l]))
        for(m in 1:length(nbana)){
          print(paste0("nbana=",nbana[m]))
          for(n in 1:length(bv)){
            print(paste0("bv=",bv[n]))
            for(o in 1:length(spazm)){
              print(paste0("spazm=",spazm[o]))
              for(p in 1:length(season)){
                print(paste0("season=",season[p]))
                for(q in 1:length(seuil)){
                  print(paste0("seuil=",seuil[q]))
                  for(r in 1:length(moments)){
                    print(paste0("moments=",moments[r]))
                    for(s in 1:length(type)){
                      print(paste0("type=",type[s]))
                      for(t in 1:length(q.threshold)){
                        print(paste0("q.threshold=",q.threshold[t]))
                        for(u in 1:length(schaake)){
                          print(paste0("schaake=",schaake[u]))
                          run.analogs(k = k[i],dist = dist,start = start,end = end,
                                      rean = rean[j],period = period,nbdays = nbdays[l],nbana = nbana[m],
                                      bv = bv[n],spazm = spazm[o],season = season[p],seuil = seuil[q],moments = moments[r],
                                      type = type[s],q.threshold = q.threshold[t],schaake = schaake[u],var2 = var2,nbana2 = nbana2)
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

# Renvoie les indices des 5% des journees/sequences les plus analogues (pour les sequences, somme des scores de chaque ieme journee des deux sequences)
# Plus efficace pour nbdays=3 (iterations dependantes)
save.ana <- function(k,dist,nbdays,rean,period="past"){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2021-07-17")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  ind <- match(dates.ana,dates.rean) # indice de l'archive de recherche des analogues
  
  # Import des distances 
  print("Import des distances")
  gc()
  dist.vec<-getdist(k = k,dist = dist,start = start.end.rean[1],end = start.end.rean[2],rean = rean,threeday = F,period = "past")
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  
  U<-c(0,(N-1):1)
  sU<-sapply(1:(N-1),function(x) sum(U[1:x]))
  gc()
  
  # Dates analogues de chaque journee/sequence
  print("Selection des dates analogues")
  nei.long <- vector("list",length=N) # Analogues sur toute periode
  nei.short <- vector("list",length=N) # Analogues sur periode reduite
  di.mat <- NULL
  
  for(i in 1:N){
    gc()
    if (i %% 5==0) {print(paste0(i,"/",N-nbdays+1))}
    
    if(nbdays==1){ # si nbdays = 1
      di <- getdist4i(i,dist.vec,N,sU)
      out <- i # on ne retirera des analogues que le jour j (on accepte veille et lendemain)
    }
    
    if(nbdays>1){# si nbdays > 1
      out <- (i-nbdays+1):(i+nbdays-1) # on retirera les sequences analogues qui chevauchent la sequence i (qui ont au moins 1 jour dans la sequence cible)
      
      if (i==1){
        for (j in 1:nbdays) di.mat<-rbind(di.mat,getdist4i(j,dist.vec,N,sU)) # distances journalieres pour chaque journee de la sequence
      }else {
        di.mat<-di.mat[-1,]
        di.mat<-rbind(di.mat,getdist4i(i+nbdays-1,dist.vec,N,sU)) # distances journalieres pour chaque journee de la sequence
      }
      tmp.mat<-matrix(NA,ncol=N-nbdays+1,nrow=nbdays)
      for (j in (1:nbdays)) tmp.mat[j,]<-di.mat[j,j:(N-nbdays+j)] # decalage des distances pour calcul analogie sur nbdays jours
      di <- apply(tmp.mat,2,sum)
    }
    
    # Indices analogues
    gc()
    soso<-sort(di,index.return=TRUE)
    soso$ix <- soso$ix[!(soso$ix %in% out)] # on retire les journees/sequences trop proches de la cible
    #soso$ix[1:(0.2*N)] # on garde les 20% dans un permier temps
    nei.long[[i]] <- soso$ix[1:(0.05*N)] # on garde les 5% analogues pour etre large
    soso$ix <- soso$ix[soso$ix %in% ind]
    nei.short[[i]] <- soso$ix[1:(0.05*n)] # on garde les 5% analogues pour etre large
    nei.short[[i]] <- nei.short[[i]] - ind[1] + 1 # on remet les indices dans le referentiel de l'archive ou l'on cherche les analogues
  }
  
  # Enregistrement
  print("Enregistrement")
  nei <- nei.long
  save(nei,file=paste0(get.dirstr(k,rean,period),"save.ana/ana_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.rean[1],"_",start.end.rean[2],".Rdata"))
  
  nei <- nei.short
  save(nei,file=paste0(get.dirstr(k,rean,period),"save.ana/ana_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
}

# Classe les analogues de save.ana par proximite en termes de MPD
save.ana.dP <- function(k,dist,nbdays,rean,period="past",nb=111){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist")
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Imports
  load(paste0(get.dirstr(k,rean,period),"save.ana/ana_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  dP <- get.descriptor(descriptor = "dP",k = k,dist = dist,nbdays = nbdays,start = start.end.rean[1],end = start.end.rean[2],standardize = F,rean = rean,threeday = F,
                       desais = F,period = "past",start.ana = start.end.ana[1],end.ana = start.end.ana[2])
  
  # Traitement
  delta <- match(dates.ana,dates.rean)[1]
  nei.tmp <- lapply(nei,function(v){v <- v+delta-1;v[1:nb]})
  #tmp <- c(mapply(function(v,x){dP[v]-dP[x]},v=nei.tmp,x=1:length(dP)))
  #plot(density(na.omit(tmp)))
  #abline(v=0,col="red")
  #quantile(na.omit(tmp),0.39)
  
  fun <- function(ind.i,ind.ana,dP){
    sort(abs(dP[ind.i] - dP[ind.ana]),index.return=T)$ix
  }
  
  
  for(i in 1:N){
    if(i%%50==0) print(i)
    order.i <- fun(i,nei.tmp[[i]],dP)
    nei[[i]] <- nei.tmp[[i]][order.i]
  }
  
  # On repasse en bons indices
  nei <- lapply(nei,function(v){v-delta+1})
  
  # Export
  save(nei,file = paste0(get.dirstr(k,rean,period),"save.ana/ana_dP_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  
}

# Extrait de save.ana les analogues etant dans une fenetre de + ou - 1 mois
save.ana.sea <- function(k,dist,nbdays,rean,period="past"){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import
  load(paste0(get.dirstr(k,rean,period),"save.ana/ana_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  
  # Traitement
  neibon <- mapply(function(x,y){ind <- get.ind.season.ana(y,dates.ana);x[x %in% ind]},x=nei,y=dates.rean)
  
  len <- unlist(lapply(neibon,length)) # si nbre ana dans la saison < 45
  pos <- which(len<45)
  print(paste0("Nbre de jours avec moins de 45 analogues dans + ou - 1 mois: ",length(pos)))
  for(i in pos){
    nei.i <- nei[[i]]
    ind.i <- get.ind.season.ana(dates.rean[i],dates.ana)
    nei.i.out <- nei.i[!(nei.i %in% ind.i)]
    neibon[[i]] <- c(neibon[[i]],nei.i.out[1:(45-len[i])])
  }
  
  # Export
  nei <- neibon
  save(nei,file = paste0(get.dirstr(k,rean,period),"save.ana/ana_season_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  
}

# On retire les indices 22278 et 22279 qui ne doivent pas apparaitre dans les analogues 3j sur periode de precip 1950-2010, sans tout recalculer
save.ana.correct <- function(k,dist,nbdays,rean,period="past"){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist")
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2010-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import
  load(paste0(get.dirstr(k,rean,period),"save.ana/ana_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  # Traitement
  nei <- lapply(nei,function(v){v[v<=n]}) # on retire les indices 22278 et 22279
  # Export
  save(nei,file=paste0(get.dirstr(k,rean,period),"save.ana/ana_",dist,"_",rean,"_k",k,"_mean",nbdays,"day_",start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".Rdata"))
  
}

# Renvoie les indices des 5% des journees/sequences les plus analogues (pour les sequences, somme des scores de chaque ieme journee des deux sequences)
# Plus efficace pour nbdays=1 (iterations independantes)
save.ana.par <- function(k,dist,nbdays,rean,period="past",season=F,ncores=6){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Indices utiles
  ind <- match(dates.ana,dates.rean) # indice de l'archive de recherche des analogues
  
  if(season){
    ind.sea <- list() # indices de recherche des analogues dans la meme saison
    dat <- getdates("2020-01-01","2020-12-31")
    for(i in 1:length(dat)){
      ind.sea[[i]] <- get.ind.season.ana(dat[i],dates.rean)
    }
    pos.dat <- match(substr(dates.rean,6,10),substr(dat,6,10))
  }
  
  # Import des distances 
  print("Import des distances")
  gc()
  dist.vec<-getdist(k = k,dist = dist,start = start.end.rean[1],end = start.end.rean[2],rean = rean,threeday = F,period = "past")
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  
  U<-c(0,(N-1):1)
  sU<-sapply(1:(N-1),function(x) sum(U[1:x]))
  gc()
  
  # Dates analogues de chaque journee/sequence
  print("Selection des dates analogues")
  #nei.long <- vector("list",length=N) # Analogues sur toute periode
  #nei.short <- vector("list",length=N) # Analogues sur periode reduite
  di.mat <- NULL
  
  print(paste0("Parallelisation sur ",ncores, " coeurs"))
  outfile <- paste0(get.dirstr(k,rean,period),"save.ana/calcul.txt")
  print(paste0("Logfile for // loop : ",outfile))
  cl <- makeCluster(ncores, outfile=outfile) 
  registerDoParallel(cl)
  
  # Pour aller voir calcul.txt
  # cmd
  # cd dossier (tapper uniquement I: pour aller sur le DD externe)
  # powershell Get-Content calcul.txt -Wait
  # CTRL+C deux fois pour fermer le .txt
  
  nei <- foreach (i=1:(N-nbdays+1)) %dopar%{
    
    source("2_Travail/1_Past/getdist4i.R", encoding = 'UTF-8')
    
    gc()
    if (i %% 50==0) {print(paste0(i,"/",N-nbdays+1))}
    
    if(nbdays==1){ # si nbdays = 1
      di <- getdist4i(i,dist.vec,N,sU)
      out <- i # on ne retirera des analogues que le jour j (on accepte veille et lendemain)
    }else{# si nbdays > 1
      out <- (i-nbdays+1):(i+nbdays-1) # on retirera les sequences analogues qui chevauchent la sequence i (qui ont au moins 1 jour dans la sequence cible)
      for (j in 1:nbdays) {di.mat<-rbind(di.mat,getdist4i(i+j-1,dist.vec,N,sU)[j:(N-nbdays+j)])} # en parallelisation: oblige d'etre independant de l'iteration i-1 donc getdist4i doit etre fait nbdays fois a chaque iteration
      di <- apply(di.mat,2,sum)
    }
    
    # Indices analogues
    gc()
    soso<-sort(di,index.return=TRUE)$ix
    soso <- soso[!(soso %in% out)] # on retire les journees/sequences trop proches de la cible
    if(season) soso <- soso[soso %in% ind.sea[[pos.dat[i]]]] # on retire les journees qui ne sont pas dans + ou - 30 jours
    if(!season) soso[1:(0.2*N)] # on garde les 20% dans un permier temps
    soso
  }
  
  stopCluster(cl)
  print(paste0("Fin calcul indicateurs a : ",Sys.time()))
  
  # Traitement final
  print("Traitement final")
  nei.short <- lapply(nei,function(v){x <- v[v %in% ind];x <- x-ind[1]+1;x <- x[1:(0.05*n)];x}) # on conserve les 5% les plus proches
  nei.long <- lapply(nei,function(v){v[1:(0.05*N)]})
  gc()
  
  # Enregistrement
  print("Enregistrement")
  nei <- nei.long
  short <- F
  save(nei,file=get.path(save_ana = T,short = short,k = k,rean = rean,period = period,season = season,dP = F,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
  
  nei <- nei.short
  short <- T
  save(nei,file=get.path(save_ana = T,short = short,k = k,rean = rean,period = period,season = season,dP = F,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2]))
}

# Sauvegarde les precipitations analogues pour analogie classique
save.precip.ana <- function(k,dist,nbdays,nbana=0.2,bv,spazm=T,rean,period="past",season=F,dP=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des dates analogues et des precip
  load(get.path(save_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],short=T))
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  # Selection des precip et enregistrement
  nbnei <- round(nbana*0.01*n,0)
  print(paste0("Nombre d'analogues selectionnes: ",nbnei))
  
  precip.ana <- vector("list",length=N)
  
  for(i in 1:N){
    if(i%%50 == 0) print(paste0(i,"/",N))
    precip.ana[[i]] <- precip[nei[[i]][1:nbnei]]
  }
  save(precip.ana,file = get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana))
}

# Sauvegarde les precipitations analogues pour analogie classique, pour tous les BVs
save.precip.ana.all <- function(k,dist,nbdays,nbana=0.2,rean,period="past",season=F,dP=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  # Import des dates analogues
  load(get.path(save_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],short=T))
  
  # Import des precip
  load("2_Travail/Data/Precip/SPAZM/dailyspazmSecteurHydro_1950-2017.Rdata")
  precip <- matdaily;rm(matdaily)
  
  # Selection des precip et enregistrement
  nbnei <- round(nbana*0.01*n,0)
  print(paste0("Nombre d'analogues selectionnes: ",nbnei))
  
  precip.ana <- vector("list",length=ncol(precip))
  
  for(i in 1:ncol(precip)){
    print(paste0(i,"/",ncol(precip)))
    for(j in 1:N){
      if(j%%50 == 0) print(paste0(j,"/",N))
      precip.ana[[i]][[j]] <- unname(precip[nei[[j]][1:nbnei],i])
    }
  }
  save(precip.ana,file = get.path(save_precip_ana_all = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],nbana = nbana))
}

# Sauvegarde les indices et les precipitations analogues pour analogie en 2 etapes
save.precip.ana.twostage <- function(k,dist1="TWS",var="vv700",nbdays,nbana1=0.5,nbana2=0.2,bv,spazm=T,rean,period="past",season=F,dP=F){
  
  # Dates utiles
  start.end.rean <- get.start.end.rean(rean,period,"dist",k)
  dates.rean <- getdates(start.end.rean[1],as.character(as.Date(start.end.rean[2])-nbdays+1))
  N <- length(dates.rean)
  
  end.var <- get.start.end.rean(rean,period,"dist",k,var = var)[2]
  print(end.var)
  dates.var <- getdates(start.end.rean[1],as.character(as.Date(end.var)-nbdays+1))
  M <- length(dates.var)
  
  start.end.ana <- c("1950-01-01","2017-12-31")
  dates.ana <- getdates(start.end.ana[1],as.character(as.Date(start.end.ana[2])-nbdays+1))
  n <- length(dates.ana)
  
  shif <- match(dates.ana,dates.rean)
  
  # Import des dates analogues et des precip
  load(get.path(save_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist1,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = start.end.rean[2],start.ana = start.end.ana[1],end.ana = start.end.ana[2],short=T))
  length(nei) <- length(dates.var)
  precip <- get.precip(nbdays = nbdays,start = start.end.ana[1],end = start.end.ana[2],bv = bv,spazm = spazm)
  
  # Selection premiere etape
  nbnei1 <- round(nbana1*0.01*n,0)
  print(paste0("Nombre d'analogues etape 1: ",nbnei1))
  nei <- lapply(nei,function(v){v[1:nbnei1]})
  
  # Selection deuxieme etape
  nbnei2 <- round(nbana2*0.01*n,0)
  print(paste0("Nombre d'analogues etape 2: ",nbnei2))
  lim <- get.lon.lat.var.bv(bv = bv)
  data <- getdata(k = 1,day0 = start.end.rean[1],day1 = end.var,rean = rean,lim.lon = lim$lim.lon,lim.lat = lim$lim.lat,var = var)
  precip.ana <- vector("list",length=length(nei))
  
  if(shif[1]!=1) nei <- lapply(nei,function(v){v+shif[1]-1}) # si decalage de dates entre rean et archive analogues, on se remet dans le referentiel de la reanalyse pour calcul analogie
  fct<-function(v) sum(abs(v))
  
  for(i in 1:length(nei)){
    if(i%%50==0) print(paste0(i,"/",length(nei)))
    data.i <- array(data = data[,,i],dim = c(dim(data[,,i]),nbnei1))
    data.j <- data[,,nei[[i]]]
    dist <- apply(data.i-data.j,3,function(v) sqrt(mean(v^2)))
    
    # TWS
    #gradloni <- makegrad(data.i,1);gradlati <- makegrad(data.i,2)
    #gradlonj <- makegrad(data.j,1);gradlatj <- makegrad(data.j,2)
    #dif_grad<-apply(gradloni-gradlonj,3,fct)+apply(gradlati-gradlatj,3,fct)
    #max_grad<-apply(pmax(abs(gradloni),abs(gradlonj)),3,sum)+apply(pmax(abs(gradlati),abs(gradlatj)),3,sum)
    #dist <- dif_grad/max_grad/2
    
    ind <- sort(dist,index.return=T)$ix[1:nbnei2]
    if(shif[1]!=1) ind <- ind-shif[1]+1 # on se remet dans le referentiel de l'archive analogue pour precipitations
    nei[[i]] <- nei[[i]][ind]
    precip.ana[[i]] <- precip[nei[[i]]]
  }
  
  # Export
  save(nei,file = get.path(save_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist1,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = end.var,start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana,var2 = var,short = T))
  save(precip.ana,file = get.path(save_precip_ana = T,k = k,rean = rean,period = period,season = season,dP = dP,dist = dist1,nbdays = nbdays,start.rean = start.end.rean[1],end.rean = end.var,start.ana = start.end.ana[1],end.ana = start.end.ana[2],bv = bv,spazm = spazm,nbana = nbana1,var2 = var,nbana2 = nbana2))
}

