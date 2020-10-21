source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Calcul de la correlation avec tous les lissages entre un indicateur et les cumuls de precipitation aux pas de temps journalier et annuels, pour une saison
compute.cor.descr.precip <- function(k,dist,bv="Isere-seul",sais="all",start="1950-01-01",end="2010-12-31",rean){
  
  descr <- c("celnei","dP","mean","singnei","rsingnei")
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  
  # Import indicateur
  if(rean=="20CR") start.file <- "1851-01-01"
  if(rean=="ERA20C") start.file <- "1900-01-01"
  des <- matrix(NA,length(getdates(start.file,end)),length(descr))
  for(i in 1:length(descr)){
  des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = 1,start = start.file,end = end,
                                      standardize = F,rean = rean,threeday = F,period = "past",start.ana = "1979-01-01",
                                      end.ana = "2002-08-31")
  }
  
  if(start.file != start){
    delta <- length(getdates(start.file,start))
    des <- apply(des,2,function(v) v[delta:length(v)])
  }
  
  colnames(des) <- descr
  
  # Import type de temps
  wp.num <- c(1,2,5,8)
  wp.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  wp <- get.wp(nbdays = 1,start = start,end=end,risk = F,bv = "Isere",agreg = T)
  tt <- matrix(NA,length(dates),length(wp.num))
  for(i in 1:length(wp.num)){
    tt[,i] <- wp
    tt[tt[,i] != wp.num[i],i] <- 0
  }
  
  colnames(tt) <- wp.name
  
  # Import NAO
  NAO <- get.nao(start=head(year,1),end = tail(year,1),sais = "all",daily = T)
  NAO <- NAO[,4]
  
  # Import precip
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = F)
  
  # Calcul correlation au pas de temps journalier
  explain <- cbind(des,tt,NAO) # tous les indicateurs et les TT (variables explicatives)
  
  if(sais=="winter"){pos <- substr(dates,6,7) %in% c("12","01","02");explain[!pos,] <- NA; precip[!pos] <- NA}
  if(sais=="spring"){pos <- substr(dates,6,7) %in% c("03","04","05");explain[!pos,] <- NA; precip[!pos] <- NA}
  if(sais=="summer"){pos <- substr(dates,6,7) %in% c("06","07","08");explain[!pos,] <- NA; precip[!pos] <- NA}
  if(sais=="autumn"){pos <- substr(dates,6,7) %in% c("09","10","11");explain[!pos,] <- NA; precip[!pos] <- NA}
  explain.daily <- na.omit(explain)
  precip.daily <- na.omit(precip)
  
  corr <- matrix(NA,365*10,ncol(explain.daily))
  len.sais <- 365
  if(sais!="all") len.sais <- 91
  
  for(i in 1:(len.sais*10)){
    if(i %% 50 == 0) print(i)
    precip.tmp <- rollapply(precip.daily,i,mean,partial=F,fill=NA)
    explain.tmp <- apply(explain.daily,2,function(v) rollapply(v,i,mean,partial=F,fill=NA))
    corr[i,] <- round(cor(precip.tmp,explain.tmp,use="pairwise.complete.obs"),3)  # cor fait directement la correlation par colonne entre un vecteur et une matrice
  }
  
  colnames(corr) <- colnames(explain.daily)
  save(corr,file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_daily_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  
  # Calcul correlation au pas de temps annuel/saisonnier
  precip.yearly <- aggregate(precip,by=list(year),mean,na.rm=T)[,2]
  explain.yearly <- apply(explain,2,function(v) aggregate(v,by=list(year),mean,na.rm=T)[,2])
  
  corr <- matrix(NA,10,ncol(explain.yearly))
  
  for(i in 1:10){
    precip.tmp <- rollapply(precip.yearly,i,mean,partial=F,fill=NA)
    explain.tmp <- apply(explain.yearly,2,function(v) rollapply(v,i,mean,partial=F,fill=NA))
    corr[i,] <- round(cor(precip.tmp,explain.tmp,use="pairwise.complete.obs"),3) # cor fait directement la correlation par colonne entre un vecteur et une matrice
  }
  
  colnames(corr) <- colnames(explain.yearly)
  save(corr,file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_yearly_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  
  # Calcul correlation sur les max annuels/saisonniers 1 jour
  maxi.1 <- aggregate(precip,by=list(year),which.max)
  add <- match(paste0(maxi.1[,1],"-01-01"),dates)
  maxi.1 <- maxi.1[,2]+add-1
  
  precip.max <- precip[maxi.1]
  explain.max <- explain[maxi.1,]
  
  corr <- matrix(NA,10,ncol(explain.max))
  
  for(i in 1:10){
    precip.tmp <- rollapply(precip.max,i,mean,partial=F,fill=NA)
    explain.tmp <- apply(explain.max,2,function(v) rollapply(v,i,mean,partial=F,fill=NA))
    corr[i,] <- round(cor(precip.tmp,explain.tmp,use="pairwise.complete.obs"),3) # cor fait directement la correlation par colonne entre un vecteur et une matrice
  }
  
  colnames(corr) <- colnames(explain.max)
  save(corr,file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_max_1day_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  
  # Calcul correlation sur les max annuels/saisonniers 3 jours
  precip.3 <- rollapply(precip,3,mean)
  explain.3 <- apply(explain,2,function(v) rollapply(v,3,mean))
  maxi.3 <- aggregate(precip.3,by=list(year[1:length(precip.3)]),which.max)
  add <- match(paste0(maxi.3[,1],"-01-01"),dates)
  maxi.3 <- maxi.3[,2]+add-1
  
  precip.max <- precip[maxi.3]
  explain.max <- explain[maxi.3,]
  
  corr <- matrix(NA,10,ncol(explain.max))
  
  for(i in 1:10){
    precip.tmp <- rollapply(precip.max,i,mean,partial=F,fill=NA)
    explain.tmp <- apply(explain.max,2,function(v) rollapply(v,i,mean,partial=F,fill=NA))
    corr[i,] <- round(cor(precip.tmp,explain.tmp,use="pairwise.complete.obs"),3) # cor fait directement la correlation par colonne entre un vecteur et une matrice
  }
  
  colnames(corr) <- colnames(explain.max)
  save(corr,file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_max_3day_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  
}

# Calcul des distances de maniere generique
compute_dist_gen_past<-function(k,dist,start="1950-01-01",end="2011-12-31",rean,period="past"){
  if (k %in% 1:2) {
    if (dist %in% c("TWS","RMSE","RMSE_I","RMSE_II","Mahalanobis")){
      if (dist=="TWS") dist.list<-compute_TWS(k,start,end,rean)
      if (dist=="RMSE") dist.list<-compute_RMSE(k,start,end,rean)
      if (dist=="RMSE_I") dist.list<-compute_RMSE_I(k,start,end,rean)
      if (dist=="RMSE_II") dist.list<-compute_RMSE_II(k,start,end,rean)
      if (dist=="Mahalanobis") dist.list<-compute_mahalanobis(k,start,end,rean)
      save(dist.list,file=paste0(get.dirstr(k,rean,period),"compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
    }
    else{ # si on demande un nTWS, sTWS, nRMSE ou sRMSE, les distances sont normalisees par la moyenne ou l'ecart type des distances
      dist0<-substr(dist,2,nchar(dist))
      load(file=paste0(get.dirstr(k,rean,period),"compute_dist/",dist0,"_member",member,"_Z",Z,"_",start,"_",end,".Rdata"))
      type<-substr(dist,1,1)
      print(type)
      if (type=="n") norm<-mean(unlist(dist.list))
      if (type=="s") norm<-sd(unlist(dist.list))
      print(norm)
      for (i in 1:length(dist.list)) dist.list[[i]]<-dist.list[[i]]/norm
      
      save(dist.list,file=paste0(get.dirstr(k,rean,period),"compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
    }
  }
}

# Calcul des indicateurs
compute_criteria_past<-function(k,dist,start="1950-01-01",end="2011-12-31",start.ana="1950-01-01",end.ana="2011-12-31",update=FALSE,rean,threeday=FALSE,rev=FALSE,period="past"){
  gc()
  
  print(paste0(get.dirstr(k,rean,period),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".Rdata"))
  
  if (update) {
    load(file=paste0(get.dirstr(k,rean,period),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".Rdata"))
    coln<-colnames(criteria)
  }
  
  dist.vec<-getdist(k,ifelse(rev,"TWS",dist),start,end,rean,threeday,period) # si rev, on prend les nei en TWS
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  dates<-getdates(start,ifelse(threeday,as.character(as.Date(end)-3+1),end))
  N<-length(dates)
  gc()
  
  U<-c(0,(N-1):1); # U = 0, 22644, 22643, 22642, ...
  sU<-sapply(1:(N-1),function(x) sum(U[1:x])) # somme cumulee de U: on fait la somme de U[1], U[1:2], etc pour obtenir la position de la derniere distance qui separe chaque date
  gc()
  
  seq.all <- getdates(start,end) # fenetre entiere
  seq.ana <- getdates(start.ana,end.ana) # fenetre de recherche des analogues
  ind <- match(seq.ana,seq.all)
  
  if (!update) {
    coln.new<-c("cel","sing05","q05")
  }
  if (update) {
    coln.new<-c("celnei","singnei","rsingnei")
  }
  
  criteria.new<-matrix(NA,ncol=length(coln.new),nrow=N)
  colnames(criteria.new)<-coln.new
  
  for (i in 1:N){
    
    ddi<-getdist4i(i,dist.vec,N,sU)
    di<-ddi
    n<-length(ind) # longueur de la fenetre de recherche des analogues
    
    gc()
    
    soso<-sort(di,index.return=TRUE) # classement par plus petit score, et donne les positions
    soso$ix <- soso$ix[soso$ix %in% ind] # on ne garde que les plus proches faisant partie de la fenetre de recherche des analogues
    qi05<-di[soso$ix[(0.005*n)]]     # quantile 0.5%
    idi05<-soso$ix[1:(0.005*n)] # recupere la position des 0.5% les plus proches
    idi05 <- idi05[idi05!=i] # le jour cible ne peut etre un analogue de lui-meme
    
    tmp<-NULL
    
    for (cc in coln.new){
      # Celerite
      if (cc=="cel") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,ddi[i-1])}
      if (cc=="celnei") tmp<-c(tmp,mean(criteria[idi05,"cel"],na.rm=TRUE))
      
      ## Singularite
      if (cc=="sing05") tmp<-c(tmp,mean(di[idi05]))
      if (cc=="singnei") tmp <- c(tmp,mean(criteria[idi05,"sing05"],na.rm=TRUE))
      
      #Singularite relative
      if (cc=="q05") tmp<-c(tmp,qi05)
      if (cc=="rsingnei") tmp <- c(tmp,mean(criteria[idi05,"rsing05"],na.rm=TRUE))
      
      # dP
      #if (cc=="dPnei") tmp <- c(tmp,mean(criteria[idi05,"dP"],na.rm=TRUE))
    }
    
    idi05m1<-idi05
    criteria.new[i,]<-tmp
    
    if (i %% 50==0) {print(i)}
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
  
  save(criteria,file=paste0(get.dirstr(k,rean,period),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".Rdata"))
}

# Calcul de la latitude journaliere du jet
compute.lat.jet <- function(gamme=c(5450,5550),k,start="1851-01-01",end="2010-12-31",rean){
  
  # Recuperation longitude et latitude
  nc <- load.nc(rean)$nc500
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  
  info <- getinfo_window(k,rean=rean) # lon et lat de notre fenetre d'analogie
  lon <- lon[info[1,1]:(info[1,1]+info[1,2]-1)]
  lat <- lat[info[2,1]:(info[2,1]+info[2,2]-1)]
  
  # Import
  geo <- getdata(k,start,end,rean)
  
  # Latitude moyenne du jet
  lati <- apply(geo,3,function(mat){nb <- apply(mat,2,function(v) {sum(v > gamme[1] & v < gamme[2])});weighted.mean(lat,nb)})
  save(lati,file=paste0("2_Travail/1_Past/",rean,"/compute.lat.jet/weighted_mean_lat_jet_",start,"_",end,"_btw_",gamme[1],"_and_",gamme[2],".Rdata"))
  
}
  
# Reconstitution des WP de 1850 a 1950 (20CR) puis de 1900 a 1950 (ERA20C)
compute_wp_past <- function(k,dist,start="1851-01-01",end="1947-12-31",start.ana="1948-01-01",end.ana="2010-12-31",rean="20CR"){
  
  dates.useful <- getdates(start,end)
  dates.ana <- getdates(start.ana,end.ana)
  dates.all <- getdates(start,end.ana)
  
  # Import des scores d'analogie
  dist.vec<-getdist(k,dist,start,end.ana,rean,threeday=F,period="past")
  length(dist.vec) <- length(dates.useful) # on reduit aux dates utiles pour limiter memoire
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  N<-length(dates.all) # on veut quand meme les analogues sur toute la periode
  
  U<-c(0,(N-1):1); # U = 0, 22644, 22643, 22642, ...
  sU<-sapply(1:(N-1),function(x) sum(U[1:x])) # somme cumulee de U: on fait la somme de U[1], U[1:2], etc pour obtenir la position de la derniere distance qui separe chaque date
  gc()
  
  ind.ana <- match(dates.ana,dates.all)
  
  # Import des WP aggreges
  wp <- get.wp(nbdays = 1,start = start.ana,end = end.ana,risk = F,agreg = T)
  tt <- vector("list",length=length(dates.useful))
  
  #cl <- makeCluster(nb_cores,outfile="") # create a cluster with n cores
  #registerDoParallel(cl) # register the cluster
  
  for(i in 1:length(dates.useful)){
    
    if (i %% 50==0) {print(i)}
    
    # Distances
    di<-getdist4i(i,dist.vec,N,sU)
    n<-length(ind.ana) # longueur de la fenetre de recherche des analogues
    gc()
    
    soso<-sort(di,index.return=TRUE) # classement par plus petit score, et donne les positions
    soso$ix <- soso$ix[soso$ix %in% ind.ana] # on ne garde que les plus proches faisant partie de la fenetre de recherche des analogues
    idi05<-soso$ix[1:(0.005*n)] # recupere la position des 0.5% les plus proches
    idi05 <- idi05[idi05!=i] 
    
    # Type de temps
    tt[[i]] <- cbind(idi05,wp[idi05-ind.ana[1]+1]) # on va chercher le wp de la bonne journee
  }
  
  #stopCluster(cl) # shut down the cluster
  gc()
  
  save(tt,file=paste0("2_Travail/1_Past/",rean,"/compute_wp_past/wp_",start,"_",end,"_start.ana_",start.ana,"_end.ana_",end.ana,".Rdata"))
}

# Calcul de dP
get.dP <- function(k,nbdays,start="1950-01-01",end="2011-12-31",rean){
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean) 
  des <- apply(geo,3,function(x) max(x)-min(x))
  des <- rollapply(des,nbdays,mean)
}

# Import et mise en forme de la BD RTM de JD
get.event <- function(){
  
  # Import
  data <- read.csv(file = "2_Travail/Data/Event/Synthèse 118 évènements_20200618.csv",header = T,sep = ";")
  as.character(data$RTM.T)
}

# Import de la serie mensuelle de NAO
get.nao <- function(start="1950",end="2019",sais="all",daily=F){
  if(!daily){
  if(sais=="all"){
    nao <- read.csv("2_Travail/Data/NAO/annual_NAO_1865-2019.txt",skip =3,header = T,sep = "")
  }
  
  if(sais!="all"){
    nao <- read.csv("2_Travail/Data/NAO/seasonal_NAO_1865-2019.txt",skip =3,header = T,sep = "")
  if(sais=="winter") nao <- nao[,c("YEAR","DJF")]
  if(sais=="spring") nao <- nao[,c("YEAR","MAM")]
  if(sais=="summer") nao <- nao[,c("YEAR","JJA")]
  if(sais=="autumn") nao <- nao[,c("YEAR","SON")]
    }
  }else{
    nao <- read.csv(file = "2_Travail/Data/NAO/daily_NAO_NOAA.txt",sep="",skip= 4)
  }
  
  year <- seq(as.numeric(start),as.numeric(end))
  nao <- nao[nao[,"YEAR"] %in% year,]
  nao
}

# Import des dates d'evenements 1-fluviaux seuls / 2-torrentiels seuls / 3-concomittants seuls
get.dates.event <- function(type=1,start="1950-01-01",end="2018-12-31"){
  
  
  
  if(type==1){
    dates <- data[data$RTM.T>=1,1]
  }
  
  
  
  
  
  }

# Calcul de la moyenne de l'altitude du geopotentiel
get.mean <- function(k,nbdays,start="1950-01-01",end="2011-12-31",rean){
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean) 
  des <- apply(geo,3,function(x) mean(x))
  des <- rollapply(des,nbdays,mean)
}

# Trace la correlation avec tous les lissages entre un indicateur et les cumuls de precipitation aux pas de temps journalier et annuels, pour une saison
plot.cor.descr.precip <- function(k,dist,bv="Isere-seul",sais="all",start="1950-01-01",end="2010-12-31",rean){
  
  # Parametres
  ind <- c("celnei",
           "singnei",
           "rsingnei",
           "dP",
           "Atlantic",
           "Mediterranean",
           "Northeast",
           "Anticyclonic",
           "NAO")
  
  ind.name <- c("Celerity",
           "Singularity",
           "Relative singularity",
           "MPD",
           "Oceanic",
           "Mediterranean",
           "Northeast",
           "Anticyclonic",
           "NAO")
  
  colo <- c("black",
            "olivedrab3",
            "darkorange",
            "royalblue",
            #"deeppink3",
            "grey",
            "darkgoldenrod3",
            "darkorchid2",
            "red",
            "aquamarine2")
  
  xlim <- c(0,ifelse(sais=="all",3650,900))
  
  # Graphique journalier global
  load(file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_daily_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  corr <- corr[,ind]
  corr.daily <- corr
  sign <- unname(apply(corr,2,function(v) quantile(v,probs=0.5,na.rm=T)/abs(quantile(v,probs=0.5,na.rm=T)))) # on passe en positif les correlations negatives pour mieux comparer
  print(sign)
  lty <- sign; lty[is.na(lty)] <- 1; lty[lty==-1] <- 5 # pointilles pour les correlations negatives
  
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.cor.descr.precip/plot_daily_all_",bv,"_",sais,"_",start,"_",end,".png"),width = 17,height = 12,units = "cm",res=300)
  plot(corr[,1],xlim=xlim,ylim=c(-0.3,1),type="n",xlab="Smoothing length (days)",ylab="Correlation")
  grid(nx=NA,ny=NULL)
  if(sais=="all"){
  abline(v=seq(0,3650,365),lty=3,col="grey")
  } else{abline(v=seq(0,900,90),lty=3,col="grey")}
  
  for(i in 1:ncol(corr)){
    tmp <- corr[,i]*sign[i]
    lines(tmp,col=colo[i],lwd=2,lty=lty[i])
  }
  abline(h=c(0,1),lty=2)
  legend("bottomleft",col=colo,ind.name,lty=1,lwd=2,bty="n",ncol=5,cex=0.6)
  graphics.off()
  
  # Graphique journalier zoom
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.cor.descr.precip/plot_daily_zoom_",bv,"_",sais,"_",start,"_",end,".png"),width = 17,height = 12,units = "cm",res=300)
  plot(corr[,1],xlim=c(0,100),ylim=c(-0.3,1),type="n",xlab="Smoothing length (days)",ylab="Correlation")
  grid()
  abline(v=seq(0,3650,365),lty=3,col="grey")
  
  for(i in 1:ncol(corr)){
    tmp <- corr[,i]*sign[i]
    lines(tmp,col=colo[i],lwd=2,lty=lty[i])
  }
  abline(h=c(0,1),lty=2)
  legend("bottomleft",col=colo,ind.name,lty=1,lwd=2,bty="n",ncol=5,cex=0.6)
  graphics.off()
  
  # Graphique pas de temps annuel
  load(file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_yearly_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  corr <- corr[,ind]
  sign <- unname(apply(corr,2,function(v) quantile(v,probs=0.5,na.rm=T)/abs(quantile(v,probs=0.5,na.rm=T)))) # on passe en positif les correlations negatives pour mieux comparer
  print(sign)
  lty <- sign; lty[is.na(lty)] <- 1; lty[lty==-1] <- 5 # pointilles pour les correlations negatives
  
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.cor.descr.precip/plot_yearly_",bv,"_",sais,"_",start,"_",end,".png"),width = 17,height = 12,units = "cm",res=300)
  plot(corr[,1],ylim=c(-0.3,1),type="n",xlab="Smoothing length (years)",ylab="Correlation")
  grid()
  abline(v=0:10,lty=3,col="grey")
  
  for(i in 1:ncol(corr)){
    tmp <- corr[,i]*sign[i]
    lines(tmp,col=colo[i],lwd=2,lty=lty[i])
  }
  abline(h=c(0,1),lty=2)
  legend("bottomleft",col=colo,ind.name,lty=1,lwd=2,bty="n",ncol=5,cex=0.6)
  graphics.off()
  
  # Graphique mix pas de temps journalier et annuel si on est en saisonnier
  #if(sais != "all"){
  #  corr <- rbind(corr.daily[])
    
  # Je passe d'un lissage journalier a un lissage annuel donc pas la meme chose 
  # Faire un lissage journalier qui ne deborde pas sur la saison? Ce ne serait pas vraiment
  # un lissage: pas de temps journalier. Puis lissage 2 jours mais sans deborder sur la saison suivante, etc
  # jusqu'à avoir une valeur par saison: la moyenne des 90 jours. Axe des abscisses de 1 à 90 donc.
    
  #}
  
  # Graphique max annuel/saisonnier 1 jour
  load(file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_max_1day_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  corr <- corr[,ind]
  sign <- unname(apply(corr,2,function(v) quantile(v,probs=0.5,na.rm=T)/abs(quantile(v,probs=0.5,na.rm=T)))) # on passe en positif les correlations negatives pour mieux comparer
  print(sign)
  lty <- sign; lty[is.na(lty)] <- 1; lty[lty==-1] <- 5 # pointilles pour les correlations negatives
  
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.cor.descr.precip/plot_max_1day_",bv,"_",sais,"_",start,"_",end,".png"),width = 17,height = 12,units = "cm",res=300)
  plot(corr[,1],ylim=c(-0.3,1),type="n",xlab="Smoothing length (years)",ylab="Correlation")
  grid()
  abline(v=0:10,lty=3,col="grey")
  
  for(i in 1:ncol(corr)){
    tmp <- corr[,i]*sign[i]
    lines(tmp,col=colo[i],lwd=2,lty=lty[i])
  }
  abline(h=c(0,1),lty=2)
  legend("bottomleft",col=colo,ind.name,lty=1,lwd=2,bty="n",ncol=5,cex=0.6)
  graphics.off()
  
  # Graphique max annuel/saisonnier 3 jours
  load(file=paste0("2_Travail/1_Past/",rean,"/compute.cor.descr.precip/corr_max_3day_",bv,"_",sais,"_",start,"_",end,".Rdata"))
  corr <- corr[,ind]
  sign <- unname(apply(corr,2,function(v) quantile(v,probs=0,5)/abs(quantile(v,probs=0,5)))) # on passe en positif les correlations negatives pour mieux comparer
  print(sign)
  lty <- sign; lty[is.na(lty)] <- 1; lty[lty==-1] <- 5 # pointilles pour les correlations negatives
  
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.cor.descr.precip/plot_max_3day_",bv,"_",sais,"_",start,"_",end,".png"),width = 17,height = 12,units = "cm",res=300)
  plot(corr[,1],ylim=c(-0.3,1),type="n",xlab="Smoothing length (years)",ylab="Correlation")
  grid()
  abline(v=0:10,lty=3,col="grey")
  
  for(i in 1:ncol(corr)){
    tmp <- corr[,i]*sign[i]
    lines(tmp,col=colo[i],lwd=2,lty=lty[i])
  }
  abline(h=c(0,1),lty=2)
  legend("bottomleft",col=colo,ind.name,lty=1,lwd=2,bty="n",ncol=5,cex=0.6)
  graphics.off()
}

# Trace l'evolution dans le temps d'un indicateur, pour plusieurs reanalyses
plot.descr <- function(descr,k,dist,liss=5,ana.comm=F,align=F,nao=F){
  
  # Reanalyses
  rean <- c(
    "20CR",
    "ERA20C",
    "NCEP",
    "JRA55",
    "ERA40",
    "JRA55C",
    "ERA5")
  
  # Dates sur lesquelles les indicateurs sont calcules
  dates <- list(
    c("1851-01-01","2010-12-31"),
    c("1900-01-01","2010-12-31"),
    c("1950-01-01","2010-12-29"),
    c("1958-01-01","2010-12-31"),
    c("1957-09-01","2002-08-31"),
    c("1972-11-01","2012-12-30"),
    c("1979-01-01","2010-12-31")
  )
  
  # Dates sur lesquelles les analogues sont pioches
  dates.ana <- dates
  if(ana.comm) {
    for(i in 1:length(dates)) {dates.ana[[i]] <- c("1979-01-01","2002-08-31")}
    nam <- paste0("_start.ana_",dates.ana[[1]][1],"_end.ana_",dates.ana[[1]][2])
    }
  
  # Import
  dat <- list()
  des <- list()
  for(i in 1:length(rean)){
    dat[[i]] <- getdates(dates[[i]][1],dates[[i]][2])
    des[[i]] <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = dates[[i]][1],end = dates[[i]][2],
                               standardize = F,rean = rean[i],threeday = F,period = "past",start.ana = dates.ana[[i]][1],
                               end.ana = dates.ana[[i]][2])
    
    if(descr=="cel"){dat[[i]] <- dat[[i]][-1];des[[i]] <- des[[i]][-1]}
    des[[i]] <- rollapply(des[[i]],liss*365,mean,partial=T)
    if(align) des[[i]] <- des[[i]] - des[[i]][dat[[i]]=="2000-01-01"]
  }
  
  if(nao){
    naoi <- get.nao(start = "1950",end = substr(dat[[1]][length(dat[[1]])],1,4))
    annmonth <- seq(as.Date("1950-01-01"),as.Date(dat[[1]][length(dat[[1]])]),by="month") # vecteur annee mois
    pos <- match(substr(annmonth,1,7),substr(dat[[1]],1,7)) # position dans le vecteur date journalier
    indnao <- rep(NA,length(dat[[1]]))
    for(i in 2:length(pos)){indnao[pos[i-1]:(pos[i]-1)] <- naoi[i-1]} # nao repete en journalier pour avoir meme x-axis
    indnao <- rollapply(indnao,liss*365,mean,partial=T) # lissage
  }
  
  # Graphique evolution
  start.xaxis <- trunc(as.numeric(substr(dat[[1]],1,4)[1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- trunc(as.numeric(substr(dat[[1]],1,4)[length(dat[[1]])])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  pos.axis <- match(xaxis,substr(dat[[1]],1,4))
  if(is.na(pos.axis[1])) pos.axis[1] <- -length(getdates(paste0(xaxis[1],"-01-01"),dat[[1]][1]))+1
  
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.descr/plot_",descr,"_evolution_liss=",liss,ifelse(ana.comm,nam,""),ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 15,height = 9,units = "cm",res=300)
  par(mar=c(4,4,0,3))
  plot(des[[1]],type="n",xaxt="n",ylim=range(des),xlab="Year",ylab=nam2str(descr,whole=T))
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = pos.axis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
  abline(v = pos.axis,lty=3,col="grey")
  lines(des[[1]],lwd=2)
  
  for(i in 2:length(rean)){
    lines(match(dat[[i]],dat[[1]]),des[[i]],col=i,lwd=2)
  }
  
  if(nao){
    par(new=T)
    plot(indnao,pch=19,col="grey",type="l",lwd=2,xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  
  legend("topleft",inset=.02,rean,col=1:length(rean),lty=1,lwd=2,bty="n",cex=0.7)
  graphics.off()
  
  # Graphique comparaison
  comb <- t(combn(x = 1:length(rean),m=2)) # toutes les combinaisons possibles
  
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.descr/plot_",descr,"_comparaison",ifelse(ana.comm,nam,""),".png"),width = 30,height = 20,units = "cm",res=300)
  par(mfrow=c(3,7),mar=c(4,4,2,2),pty="s")
  
  for(i in 1:nrow(comb)){
    print(i)
    
    rean.1 <- comb[i,1]
    rean.2 <- comb[i,2]
    
    comm <- intersect(as.character(dat[[rean.1]]),as.character(dat[[rean.2]]))
    pos.1 <- match(as.Date(comm),dat[[rean.1]])
    pos.2 <- match(as.Date(comm),dat[[rean.2]])
    
    plot(des[[rean.1]][pos.1],des[[rean.2]][pos.2],pch=19,cex=0.1,xlim=range(des),ylim=range(des),
         xlab=paste0(nam2str(descr)," ",rean[rean.1]),ylab=paste0(nam2str(descr)," ",rean[rean.2]))
    title(main = paste0(round(length(comm)/365.25,0)," years"))
    abline(0,1,col="red")
    text(range(des)[1],range(des)[2],paste0("R² = ",round(cor(des[[rean.1]][pos.1],des[[rean.2]][pos.2])^2,2)),adj=c(0,1),cex=0.9)
  }
  graphics.off()
}

# Trace l'evolution dans le temps d'un indicateur, pour plusieurs reanalyses et par saison, avec NAO
plot.trend.descr <- function(descr,k,dist,sais="all",liss=5,ana.comm=F,align=F,nao=F){
  
  # Reanalyses
  rean <- c(
    "20CR",
    "ERA20C",
    "NCEP",
    "JRA55",
    "ERA40",
    "JRA55C",
    "ERA5")
  
  # Dates sur lesquelles les indicateurs sont calcules
  dates <- list(
    c("1851-01-01","2010-12-31"),
    c("1900-01-01","2010-12-31"),
    c("1950-01-01","2010-12-29"),
    c("1958-01-01","2010-12-31"),
    c("1957-09-01","2002-08-31"),
    c("1972-11-01","2012-12-30"),
    c("1979-01-01","2010-12-31")
  )
  
  # Dates sur lesquelles les analogues sont pioches
  dates.ana <- dates
  if(ana.comm) {
    for(i in 1:length(dates)) {dates.ana[[i]] <- c("1979-01-01","2002-08-31")}
    nam <- paste0("_start.ana_",dates.ana[[1]][1],"_end.ana_",dates.ana[[1]][2])
  }
  
  # Import
  dat <- list()
  des <- list()
  for(i in 1:length(rean)){
    dat[[i]] <- getdates(dates[[i]][1],dates[[i]][2])
    des[[i]] <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = dates[[i]][1],end = dates[[i]][2],
                               standardize = F,rean = rean[i],threeday = F,period = "past",start.ana = dates.ana[[i]][1],
                               end.ana = dates.ana[[i]][2])
    
    if(descr=="cel"){dat[[i]] <- dat[[i]][-1];des[[i]] <- des[[i]][-1]}
    
    if(sais=="winter"){
      pos <- substr(dat[[i]],6,7) %in% c("12","01","02")
      des[[i]][!pos] <- NA
    }
    if(sais=="spring"){
      pos <- substr(dat[[i]],6,7) %in% c("03","04","05")
      des[[i]][!pos] <- NA
    }
    if(sais=="summer"){
      pos <- substr(dat[[i]],6,7) %in% c("06","07","08")
      des[[i]][!pos] <- NA
    }
    if(sais=="autumn"){
      pos <- substr(dat[[i]],6,7) %in% c("09","10","11")
      des[[i]][!pos] <- NA
    }
    
    dat[[i]] <- substr(dat[[i]],1,4)
    agg <- aggregate(des[[i]],by=list(dat[[i]]),mean,na.rm=T)
    agg <- agg[!is.na(agg[,2]),]
    des[[i]] <- agg[,2]
    dat[[i]] <- agg[,1]
    des[[i]] <- rollapply(des[[i]],liss,mean,partial=F,fill=NA)
    if(align) des[[i]] <- des[[i]] - des[[i]][dat[[i]]=="2000"]
  }
  
  if(nao){
    naoi <- get.nao(start = "1865",end = "2010",sais = sais)
    nao.year <- naoi[,1]
    nao.ind  <- naoi[,2]
    nao.ind <- rollapply(nao.ind,liss,mean,partial=F,fill=NA)
    nao.ind <- c(rep(NA,nao.year[1] - as.numeric(dat[[1]][1])-1),nao.ind)
  }
  
  # Graphique evolution
  start.xaxis <- trunc(as.numeric(dat[[1]][1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- trunc(as.numeric(dat[[1]][length(dat[[1]])])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  pos.axis <- match(xaxis,dat[[1]])
  if(is.na(pos.axis[1])) pos.axis[1] <- -length(xaxis[1]-as.numeric(dat[[1]][1]))+1
  
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.descr/plot_",descr,"_evolution_liss=",liss,ifelse(ana.comm,nam,""),"_",sais,ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 15,height = 9,units = "cm",res=300)
  par(mar=c(4,4,0,3))
  plot(des[[1]],type="n",xaxt="n",ylim=range(des,na.rm = T),xlab="Year",ylab=nam2str(descr,whole=T))
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = pos.axis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
  abline(v = pos.axis,lty=3,col="grey")
  lines(des[[1]],lwd=2)
  
  for(i in 2:length(rean)){
    lines(match(dat[[i]],dat[[1]]),des[[i]],col=i,lwd=2)
  }
  
  if(nao){
    par(new=T)
    plot(nao.ind,pch=19,col="grey",type="l",lwd=2,xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  
  legend("topleft",inset=.02,rean,col=1:length(rean),lty=1,lwd=2,bty="n",cex=0.7)
  graphics.off()
  
}

# Trace l'evolution dans le temps de la latitude du jet saison et type de temps
plot.trend.lat.jet <- function(gamme=c(5450,5550),wp="all",start="1900-01-01",end="2010-12-31",rean,liss=1){
  
  # Import
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  load(paste0("2_Travail/1_Past/",rean,"/compute.lat.jet/weighted_mean_lat_jet_",start,"_",end,"_btw_",gamme[1],"_and_",gamme[2],".Rdata"))
  
  # Traitement
  
  if(wp!="all"){
    load(paste0("2_Travail/1_Past/",rean,"/compute_wp_past/wp_final_",start,"_1947-12-31_start.ana_1948-01-01_end.ana_2010-12-31_n.ana=1.Rdata"))
    lati[tt_final!=wp] <- NA
  }
  
  # Annuel
  lat.ann <- aggregate(lati,by=list(year),mean,na.rm=T)
  lat.ann[,1] <- as.numeric(lat.ann[,1])
  
  # Saisonnier
  sais <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  sais.name <- c("winter","spring","summer","autumn")
  lat.sais <- matrix(NA,length(unique(year)),5)
  lat.sais[,1] <- as.numeric(unique(year))
  for(i in 1:length(sais)){
    tmp <- lati
    tmp[!(substr(dates,6,7) %in% sais[[i]])] <- NA
    lat.sais[,i+1] <- aggregate(tmp,by=list(year),mean,na.rm=T)[,2]
  }
  
  # Lissage
  if(liss!=1){
    lat.ann[,2] <- rollapply(lat.ann[,2],liss,mean,partial=F,fill=NA)
    lat.sais[,-1] <- apply(lat.sais[,-1],2,function(v) rollapply(v,liss,mean,partial=F,fill=NA))
  }
  
  # Graphiques
  # Annuel
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.trend.lat.jet/lat_jet_ann",start,"_",end,"_liss=",liss,"_wp=",wp,".png"),width = 12,height = 9,units = "cm",res=300)
  par(mar=c(4,4,1,3))
  plot(lat.ann,type="n",xlab="Year",ylab="Jet Latitude (°)")
  grid()
  lines(lat.ann,lwd=2)
  abline(lm(lat.ann[,2]~lat.ann[,1]))
  graphics.off()
  
  # Saisonnier
  for(i in 1:4){
    png(filename = paste0("2_Travail/1_Past/",rean,"/plot.trend.lat.jet/lat_jet_",sais.name[i],"_",start,"_",end,"_liss=",liss,"_wp=",wp,".png"),width = 12,height = 9,units = "cm",res=300)
    par(mar=c(4,4,1,3))
    plot(lat.sais[,c(1,i+1)],type="n",xlab="Year",ylab="Jet Latitude (°)")
    grid()
    lines(lat.sais[,c(1,i+1)],lwd=2)
    abline(lm(lat.sais[,i+1]~lat.sais[,1]))
    graphics.off()
  }
  
}

# Trace l'évolution des precip et des max par bv pour la periode 1950-2011
plot.trend.precip <- function(bv="Isere-seul",nbdays,start="1950-01-01",end="2011-12-31",nao=F,liss=5){
  
  # Import
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  ann <- as.numeric(substr(dates,1,4))
  precip <- get.precip(nbdays,start,end,bv,spazm=F)
  
  # Traitement
  cum_an <- aggregate(precip,by=list(ann),sum)
  max_an <- aggregate(precip,by=list(ann),max)
  
  sais <- substr(dates,6,7)
  mois <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  for(i in 1:4){sais[sais %in% mois[[i]]] <- i}
  saison <- c("winter","spring","summer","autumn")
  
  cum_sais <- aggregate(precip,by=list(ann,sais),sum)
  colnames(cum_sais) <- c("year","season","value")
  cum_sais <- as.data.frame(pivot_wider(cum_sais,names_from = season,values_from = value))
  
  max_sais <- aggregate(precip,by=list(ann,sais),max)
  colnames(max_sais) <- c("year","season","value")
  max_sais <- as.data.frame(pivot_wider(max_sais,names_from = season,values_from = value))
  
  if(liss!=1){
    cum_an[,2] <- rollapply(cum_an[,2],liss,mean,partial=T)
    max_an[,2] <- rollapply(max_an[,2],liss,mean,partial=T)
    cum_sais[,2:5] <- apply(cum_sais[,2:5],2,function(v) rollapply(v,liss,mean,partial=T))
    max_sais[,2:5] <- apply(max_sais[,2:5],2,function(v) rollapply(v,liss,mean,partial=T))
  }
  
  if(nao){
    nao_an <- get.nao(start = ann[1],end = ann[length(ann)],sais = "all")
    nao_sais <- vector("list",length = 4)
    for(i in 1:4) nao_sais[[i]] <- get.nao(start = ann[1],end = ann[length(ann)],sais = saison[i])
    if(liss!=1){
      nao_an[,2] <- rollapply(nao_an[,2],liss,mean,partial=T)
      nao_sais <- lapply(nao_sais,function(v) {v[,2] <- rollapply(v[,2],liss,mean,partial=T);return(v)})
    }
  }
  
  # Figures
  
  # Cumul annuel
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_",nbdays,"days_cum_an_liss=",liss,ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
  par(mar=c(4,4,1,3))
  plot(cum_an[,c(1,2)],type="n",xlab="Year",ylab="Precipitation (mm)")
  grid()
  lines(cum_an[,c(1,2)],col="cornflowerblue",lwd=2)
  abline(lm(cum_an[,2]~cum_an[,1]))
  if(nao){
    par(new=T)
    plot(nao_an,col="grey",type="l",xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  graphics.off()
  
  # Max annuel
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_",nbdays,"days_max_an_liss=",liss,ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
  par(mar=c(4,4,1,3))
  plot(max_an[,c(1,2)],type="n",xlab="Year",ylab="Precipitation (mm)")
  grid()
  lines(max_an[,c(1,2)],col="cornflowerblue",lwd=2)
  abline(lm(max_an[,2]~max_an[,1]))
  if(nao){
    par(new=T)
    plot(nao_an,col="grey",type="l",xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  graphics.off()
  
  # Cumuls saisonniers
  for(i in 1:4){
    png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_",nbdays,"days_cum_",saison[i],"_liss=",liss,ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
    par(mar=c(4,4,1,3))
    plot(cum_sais[,c(1,i+1)],type="n",xlab="Year",ylab="Precipitation (mm)")
    grid()
    lines(cum_sais[,c(1,i+1)],col="cornflowerblue",lwd=2)
    abline(lm(cum_sais[,i+1]~cum_sais[,1]))
    if(nao){
      par(new=T)
      plot(nao_sais[[i]],col="grey",type="l",xlab="",ylab="",xaxt="n",yaxt="n")
      abline(h=0,col="grey")
      axis(side = 4)
      mtext("NAOI",side=4,line=2)
    }
    graphics.off()
  }
  
  # Max saisonniers
  for(i in 1:4){
    png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_",nbdays,"days_max_",saison[i],"_liss=",liss,ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
    par(mar=c(4,4,1,3))
    plot(max_sais[,c(1,i+1)],type="n",xlab="Year",ylab="Precipitation (mm)")
    grid()
    lines(max_sais[,c(1,i+1)],col="cornflowerblue",lwd=2)
    abline(lm(max_sais[,i+1]~max_sais[,1]))
    if(nao){
      par(new=T)
      plot(nao_sais[[i]],col="grey",type="l",xlab="",ylab="",xaxt="n",yaxt="n")
      abline(h=0,col="grey")
      axis(side = 4)
      mtext("NAOI",side=4,line=2)
    }
    graphics.off()
  }
  
}

# Trace l'évolution de l'occurrence des WP par saison
plot.trend.wp <- function(start="1900-01-01",end="2010-12-31",rean,liss=1,n.ana=115){
  
  # Import
  load(paste0("2_Travail/1_Past/",rean,"/compute_wp_past/wp_final_",start,"_1947-12-31_start.ana_1948-01-01_end.ana_2010-12-31_n.ana=",n.ana,".Rdata"))
  
  # Traitement
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  
  tt <- c(1,2,5,8)
  tt.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  wp_ann <- matrix(NA,length(unique(year)),5)
  wp_ann[,1] <- as.numeric(unique(year))
  for(i in 1:4){
    wp_ann[,i+1] <- aggregate(tt_final,by=list(year),function(v){sum(v==tt[i])/length(v)*100})[,2]
  }
  
  sais <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  sais.name <- c("winter","spring","summer","autumn")
  wp_sais <- vector("list",length=4)
  for(i in 1:length(sais)){
    tmp <- tt_final
    tmp[!(substr(dates,6,7) %in% sais[[i]])] <- NA
    
    wp_sais[[i]] <- matrix(NA,length(unique(year)),5)
    wp_sais[[i]][,1] <- as.numeric(unique(year))
    for(j in 1:4){
      wp_sais[[i]][,j+1] <- aggregate(tmp,by=list(year),function(v){v <- na.omit(v);sum(v==tt[j])/length(v)*100})[,2]
    }
  }
  
  if(liss!=1){
    wp_ann <- cbind(wp_ann[,1],apply(wp_ann[,-1],2,function(v) rollapply(v,liss,mean,partial=F,fill=NA)))
    wp_sais <- lapply(wp_sais,function(v){v <- cbind(v[,1],apply(v[,-1],2,function(w) rollapply(w,liss,mean,partial=F,fill=NA)))})
  }
  
  # Graphiques
  colo <- c("blue","red","darkgreen","darkgrey")
  start.xaxis <- trunc(as.numeric(year[1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- trunc(as.numeric(year[length(year)])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  
  # Annuel
  png(filename = paste0("2_Travail/1_Past/",rean,"/plot.trend.wp/plot_ann_evolution_",start,"_",end,"_liss=",liss,"_n.ana=",n.ana,".png"),width = 15,height = 9,units = "cm",res=300)
  par(mar=c(4,4,0,3))
  plot(wp_ann[,c(1,2)],type="n",xaxt="n",ylim=c(0,round(max(wp_ann[,2],na.rm=T)*0.1,0)*10+10),xlab="Year",ylab="WP occurence (%)")
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
  abline(v = xaxis,lty=3,col="grey")
  
  for(i in 1:4){
    lines(wp_ann[,c(1,i+1)],col=colo[i],lwd=2)
    abline(lm(wp_ann[,i+1]~wp_ann[,1]),col=colo[i])
  }
  
  legend("topleft",inset=.02,tt.name,col=colo,lty=1,lwd=2,bty="n",cex=0.6,horiz = T)
  graphics.off()

  # Saisons
  for(i in 1:4){
    png(filename = paste0("2_Travail/1_Past/",rean,"/plot.trend.wp/plot_",sais.name[i],"_evolution_",start,"_",end,"_liss=",liss,"_n.ana=",n.ana,".png"),width = 15,height = 9,units = "cm",res=300)
    par(mar=c(4,4,0,3))
    plot(wp_sais[[i]][,c(1,2)],type="n",xaxt="n",ylim=c(0,round(max(wp_sais[[i]][,2],na.rm=T)*0.1,0)*10+10),xlab="Year",ylab="WP occurence (%)")
    grid(ny=NULL,nx=NA)
    axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
    abline(v = xaxis,lty=3,col="grey")
    
    for(j in 1:4){
      lines(wp_sais[[i]][,c(1,j+1)],col=colo[j],lwd=2)
      abline(lm(wp_sais[[i]][,j+1]~wp_sais[[i]][,1]),col=colo[j])
    }
    
    legend("topleft",inset=.02,tt.name,col=colo,lty=1,lwd=2,bty="n",cex=0.6,horiz = T)
    graphics.off()
  }
  
}

# Run des fonctions
run <- function(type=1){
  
  # plot.trend.descr()
  if(type==1){
    descr <- c("cel","celnei","dP","mean","sing05","singnei","rsing05","rsingnei")
    saison <- c("all","winter","spring","summer","autumn")
    nao <- c(T,F)
    
    for(i in 1:length(descr)){
      print(i)
      for(j in 1:length(saison)){
          plot.trend.descr(descr = descr[i],k = 1,dist = "TWS",sais = saison[j],
                           liss = 20,ana.comm = T,align = F,nao = F)
      }
    }
  }
  
  # scatterplot.descr.precip
  if(type==2){
    bv <- c("Isere-seul","Drac-seul")
    typ <- c("cum","max")
    descr <- c("cel","celnei","dP","mean","sing05","singnei","rsing05","rsingnei")
    nbdays <- c(1,3)
    rean <- c("20CR","ERA20C")
    liss <- c(1,5)

    
    for(i in 1:length(bv)){
      for(j in 1:length(typ)){
        for(k in 1:length(descr)){
          for(l in 1:length(nbdays)){
            for(m in 1:length(rean)){
              for(n in 1:length(liss)){
                scatterplot.descr.precip(bv = bv[i],type = typ[j],descr = descr[k],k = 1,dist = "TWS",
                                         nbdays = nbdays[l],start = "1950-01-01",end = "2010-12-31",
                                         rean = rean[m],liss = liss[n])
              }
            }
          }
        }
      }
    }
  }
  
  # scatterplot.descr.nao
  if(type==3){
    descr <- c("cel","celnei","dP","mean","sing05","singnei","rsing05","rsingnei")
    rean <- c("20CR","ERA20C")
    start <- c("1865-01-01","1900-01-01")
    liss <- c(1,5)
    
    for(i in 1:length(descr)){
      for(j in 1:length(rean)){
        for(k in 1:length(liss)){
          scatterplot.descr.nao(descr = descr[i],k = 1,dist = "TWS",start = start[j],end="2010-12-31",
                                rean = rean[j],liss = liss[k])
        }
      }
    }
  }
  
  # compute.cor.descr.precip et plot.cor.descr.precip
  if(type==4){
    rean <- c("20CR","ERA20C")
    bv <- c("Isere-seul","Drac-seul")
    sais <- c("all","winter","spring","summer","autumn")
    
    for(i in rean){
      for(j in bv){
        for(k in sais){
          #compute.cor.descr.precip(k = 1,dist = "TWS",bv = j,sais = k,start = "1950-01-01",
                          #         end = "2010-12-31",rean = i)
          plot.cor.descr.precip(k = 1,dist = "TWS",bv = j,sais = k,start = "1950-01-01",
                                   end = "2010-12-31",rean = i)
        }
      }
    }
  }
  
  if(type==5){
    descr <- c("cel","celnei","dP","mean","sing05","singnei","rsing05","rsingnei")
    wp <- c(1,2,5,8)
    rean <- c("20CR","ERA20C")
    
    for(i in descr){
      for(j in wp){
        for(k in rean){
          scatterplot.descr.wp(descr = i,wp = j,k = 1,dist = "TWS",start = "1950-01-01",end = "2010-12-31",
                               rean = k,liss = 1)
        }
      }
    }
  }
}

# Scatterplot d'un descripteur et de NAO par saison et par reanalyse avec lissage possible
scatterplot.descr.nao <- function(descr,k,dist,start="1865-01-01",end="2010-12-31",rean="20CR",liss=1){
  
  # Import descripteur
  # Reanalyses
  reanalyses <- c(
    "20CR",
    "ERA20C",
    "NCEP",
    "JRA55",
    "ERA40",
    "JRA55C",
    "ERA5")
  
  # Dates sur lesquelles les indicateurs sont calcules
  dates <- list(
    c("1851-01-01","2010-12-31"),
    c("1900-01-01","2010-12-31"),
    c("1950-01-01","2010-12-29"),
    c("1958-01-01","2010-12-31"),
    c("1957-09-01","2002-08-31"),
    c("1972-11-01","2012-12-30"),
    c("1979-01-01","2010-12-31")
  )
  
  pos <- which(rean==reanalyses)
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = dates[[pos]][1],end = dates[[pos]][2],
                             standardize = F,rean = rean,threeday = F,period = "past",start.ana = "1979-01-01",
                             end.ana = "2002-08-31")
  if(start!= dates[[pos]][1]){
    delta <- length(getdates(dates[[pos]][1],start))
    des <- des[delta:length(des)]
  }
  if(end!= dates[[pos]][2]){
    delta <- length(getdates(end,dates[[pos]][2]))
    des <- des[1:(length(des)-delta+1)]
  }
  
  dates <- getdates(start,end)
  
  # Import NAO
  sais.name <- c("all","winter","spring","summer","autumn")
  nao <- vector("list",length=5)
  for(i in 1:5){
    nao[[i]] <- get.nao(start=substr(start,1,4),end=substr(end,1,4),sais = sais.name[i])
  }
  
  # Traitement descripteur
  sais <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  des.final <- vector("list",length=5)
  des.final[[1]] <- aggregate(des,by=list(substr(dates,1,4)),mean)
  for(i in 1:4){
    des.tmp <- des
    pos <- substr(dates,6,7) %in% sais[[i]]
    des.tmp[!pos] <- NA
    des.final[[i+1]] <- aggregate(des.tmp,by=list(substr(dates,1,4)),mean,na.rm=T)
  }
  
  # lissage
  if(liss!=1){
    des.final <- lapply(des.final,function(v){v[,2] <- rollapply(v[,2],liss,mean,partial=F,fill=NA);return(v)})
    nao <- lapply(nao,function(v){v[,2] <- rollapply(v[,2],liss,mean,partial=F,fill=NA);return(v)})
  }
  
  # Scatterplots
  png(filename = paste0("2_Travail/1_Past/",rean,"/scatterplot.descr.nao/plot_",descr,"_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
  par(mfrow=c(2,3),pty="s",mar=c(4,4,4,2))
  for(i in 1:5){
    titre <- paste0(sais.name[i]," (R = ",round(cor(des.final[[i]][,2],nao[[i]][,2],use="pairwise.complete.obs"),2),")")
    plot(des.final[[i]][,2],nao[[i]][,2],pch=19,xlab=descr,ylab="NAOI",main=titre)
  }
  
  graphics.off()
}

# Scatterplot d'un descripteur et des precipitations par saison et par reanalyse avec lissage possible
scatterplot.descr.precip <- function(bv,type="cum",descr,k,dist,nbdays,start="1950-01-01",end="2010-12-31",rean="20CR",liss=1){
  
  # Import precip
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  ann <- as.numeric(substr(dates,1,4))
  precip <- get.precip(nbdays,start,end,bv,spazm=F)
  
  # Traitement precip
  sais <- substr(dates,6,7)
  mois <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  for(i in 1:4){sais[sais %in% mois[[i]]] <- i}
  sais.name <- c("winter","spring","summer","autumn")
  
  precip_ann <- aggregate(precip,by=list(ann),ifelse(type=="cum",sum,max))
  
  precip_sais <- aggregate(precip,by=list(ann,sais),ifelse(type=="cum",sum,max))
  colnames(precip_sais) <- c("year","season","value")
  precip_sais <- as.data.frame(pivot_wider(precip_sais,names_from = season,values_from = value))
  
  if(liss!=1){
    precip_ann[,2] <- rollapply(precip_ann[,2],liss,mean,partial=T)
    precip_sais[,2:5] <- apply(precip_sais[,2:5],2,function(v) rollapply(v,liss,mean,partial=T))
  }
  
  # Import descripteur
  # Reanalyses
  reanalyses <- c(
    "20CR",
    "ERA20C",
    "NCEP",
    "JRA55",
    "ERA40",
    "JRA55C",
    "ERA5")
  
  # Dates sur lesquelles les indicateurs sont calcules
  dat <- list(
    c("1851-01-01","2010-12-31"),
    c("1900-01-01","2010-12-31"),
    c("1950-01-01","2010-12-29"),
    c("1958-01-01","2010-12-31"),
    c("1957-09-01","2002-08-31"),
    c("1972-11-01","2012-12-30"),
    c("1979-01-01","2010-12-31")
  )
  
  pos <- which(rean==reanalyses)
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = dat[[pos]][1],end = dat[[pos]][2],
                        standardize = F,rean = rean,threeday = F,period = "past",start.ana = "1979-01-01",
                        end.ana = "2002-08-31")
  if(start!= dat[[pos]][1]){
    delta <- length(getdates(dat[[pos]][1],start))
    des <- des[delta:length(des)]
  }
  if(end!= dat[[pos]][2]){
    delta <- length(getdates(end,dat[[pos]][2]))
    des <- des[1:(length(des)-delta+1)]
  }
  
  # Traitement descripteur
  des_ann <- aggregate(des,by=list(ann),mean)
  
  des_sais <- aggregate(des,by=list(ann,sais),mean)
  colnames(des_sais) <- c("year","season","value")
  des_sais <- as.data.frame(pivot_wider(des_sais,names_from = season,values_from = value))
  
  if(liss!=1){
    des_ann[,2] <- rollapply(des_ann[,2],liss,mean,partial=T)
    des_sais[,2:5] <- apply(des_sais[,2:5],2,function(v) rollapply(v,liss,mean,partial=T))
  }
  
  # Scatterplots
  png(filename = paste0("2_Travail/1_Past/",rean,"/scatterplot.descr.precip/plot_",bv,"_",descr,"_",type,"_",nbdays,"days_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
  par(mfrow=c(2,3),pty="s",mar=c(4,4,4,2))
  ylab <- paste0("Precip ",ifelse(type=="cum","Accumulation","Maxima")," (mm)")
  
  # Annuel
  titre <- paste0("all"," (R = ",round(cor(des_ann[,2],precip_ann[,2],use="pairwise.complete.obs"),2),")")
  plot(des_ann[,2],precip_ann[,2],pch=19,xlab=descr,ylab=ylab,main=titre)
  
  for(i in 1:4){
    titre <- paste0(sais.name[i]," (R = ",round(cor(des_sais[,i+1],precip_sais[,i+1],use="pairwise.complete.obs"),2),")")
    plot(des_sais[,i+1],precip_sais[,i+1],pch=19,xlab=descr,ylab=ylab,main=titre)
  }
  
  graphics.off()
}

# Scatterplot d'un descripteur et de l'occurrence des WP par saison et par reanalyse avec lissage possible
scatterplot.descr.wp <- function(descr,wp=1,k,dist,start="1950-01-01",end="2010-12-31",rean="20CR",liss=1){
  
  # Import descripteur
  # Reanalyses
  reanalyses <- c(
    "20CR",
    "ERA20C",
    "NCEP",
    "JRA55",
    "ERA40",
    "JRA55C",
    "ERA5")
  
  # Dates sur lesquelles les indicateurs sont calcules
  dates <- list(
    c("1851-01-01","2010-12-31"),
    c("1900-01-01","2010-12-31"),
    c("1950-01-01","2010-12-29"),
    c("1958-01-01","2010-12-31"),
    c("1957-09-01","2002-08-31"),
    c("1972-11-01","2012-12-30"),
    c("1979-01-01","2010-12-31")
  )
  
  pos <- which(rean==reanalyses)
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = dates[[pos]][1],end = dates[[pos]][2],
                        standardize = F,rean = rean,threeday = F,period = "past",start.ana = "1979-01-01",
                        end.ana = "2002-08-31")
  if(start!= dates[[pos]][1]){
    delta <- length(getdates(dates[[pos]][1],start))
    des <- des[delta:length(des)]
  }
  if(end!= dates[[pos]][2]){
    delta <- length(getdates(end,dates[[pos]][2]))
    des <- des[1:(length(des)-delta+1)]
  }
  
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  
  # Traitement descripteur
  sais <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  sais.name <- c("all","winter","spring","summer","autumn")
  des.final <- vector("list",length=5)
  des.final[[1]] <- aggregate(des,by=list(year),mean)
  for(i in 1:4){
    des.tmp <- des
    pos <- substr(dates,6,7) %in% sais[[i]]
    des.tmp[!pos] <- NA
    des.final[[i+1]] <- aggregate(des.tmp,by=list(substr(dates,1,4)),mean,na.rm=T)
  }
  
  # Import WP et traitement
  tt <- get.wp(nbdays = 1,start,end,agreg=T)
  tt_ann <- aggregate(tt,by=list(year),function(v){sum(v==wp)/length(v)*100})
  
  tt_sais <- matrix(NA,length(unique(year)),4)
  
  for(i in 1:4){
    tmp <- tt
    tmp[!(substr(dates,6,7) %in% sais[[i]])] <- NA
    tt_sais[,i] <- aggregate(tmp,by=list(year),function(v){v <- na.omit(v);sum(v==wp)/length(v)*100})[,2]
  }
  
  tt <- cbind(tt_ann,tt_sais)
  
  # lissage
  if(liss!=1){
    des.final <- lapply(des.final,function(v){v[,2] <- rollapply(v[,2],liss,mean,partial=F,fill=NA);return(v)})
    tt[,-1] <- apply(tt[,-1],2,function(v) rollapply(v,liss,mean,partial=T))
    }
  
  # Scatterplots
  png(filename = paste0("2_Travail/1_Past/",rean,"/scatterplot.descr.wp/plot_",descr,"_wp=",wp,"_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
  par(mfrow=c(2,3),pty="s",mar=c(4,4,4,2))
  for(i in 1:5){
    titre <- paste0(sais.name[i]," (R = ",round(cor(des.final[[i]][,2],tt[,i+1],use="pairwise.complete.obs"),2),")")
    plot(des.final[[i]][,2],tt[,i+1],pch=19,xlab=descr,ylab="WP Occurrence (%)",main=titre)
  }
  
  graphics.off()
}

# Scatterplot de NAO et des precipitations par saison avec lissage possible
scatterplot.nao.precip <- function(bv,start="1950-01-01",end="2010-12-31",liss=1){
  
  # Import precip
  dates <- getdates(start,end)
  ann <- as.numeric(substr(dates,1,4))
  precip <- get.precip(nbdays,start,end,bv,spazm=F)
  
  # Traitement precip
  sais <- substr(dates,6,7)
  mois <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  for(i in 1:4){sais[sais %in% mois[[i]]] <- i}
  sais.name <- c("winter","spring","summer","autumn")
  
  precip_ann <- aggregate(precip,by=list(ann),sum)
  
  precip_sais <- aggregate(precip,by=list(ann,sais),sum)
  colnames(precip_sais) <- c("year","season","value")
  precip_sais <- as.data.frame(pivot_wider(precip_sais,names_from = season,values_from = value))
  
  if(liss!=1){
    precip_ann[,2] <- rollapply(precip_ann[,2],liss,mean,partial=T)
    precip_sais[,2:5] <- apply(precip_sais[,2:5],2,function(v) rollapply(v,liss,mean,partial=T))
  }
  
  # Import NAO
  nao <- matrix(NA,length(unique(ann)),6)
  nao[,1] <- unique(ann)
  sais.name <- c("all",sais.name)
  for(i in 1:5){
    nao[,i+1] <- get.nao(start=substr(start,1,4),end=substr(end,1,4),sais = sais.name[i])[,2]
  }
  
  if(liss!=1){
    nao[,-1] <- apply(nao[,-1],2,function(v) rollapply(v,liss,mean,partial=T))
  }
 
  
  # Scatterplots
  png(filename = paste0("2_Travail/1_Past/Rresults/scatterplot.nao.precip/plot_",bv,"_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
  par(mfrow=c(2,3),pty="s",mar=c(4,4,4,2))
  ylab <- "Precip Accumulation (mm)"
  
  # Annuel
  titre <- paste0("all"," (R = ",round(cor(nao[,2],precip_ann[,2],use="pairwise.complete.obs"),2),")")
  plot(nao[,2],precip_ann[,2],pch=19,xlab="NAOI",ylab=ylab,main=titre)
  sais.name <- sais.name[-1]
  for(i in 1:4){
    titre <- paste0(sais.name[i]," (R = ",round(cor(nao[,i+2],precip_sais[,i+1],use="pairwise.complete.obs"),2),")")
    plot(nao[,i+2],precip_sais[,i+1],pch=19,xlab="NAOI",ylab=ylab,main=titre)
  }
  
  graphics.off()
}

# Scatterplot de l'occurrence des WP et des precipitations par saison avec lissage possible
scatterplot.wp.precip <- function(wp=1,bv,start="1950-01-01",end="2010-12-31",liss=1){
  
  # Import precip
  dates <- getdates(start,end)
  ann <- as.numeric(substr(dates,1,4))
  precip <- get.precip(nbdays,start,end,bv,spazm=F)
  
  # Traitement precip
  sais <- substr(dates,6,7)
  mois <- list(c("12","01","02"),c("03","04","05"),c("06","07","08"),c("09","10","11"))
  for(i in 1:4){sais[sais %in% mois[[i]]] <- i}
  sais.name <- c("winter","spring","summer","autumn")
  
  precip_ann <- aggregate(precip,by=list(ann),sum)
  
  precip_sais <- aggregate(precip,by=list(ann,sais),sum)
  colnames(precip_sais) <- c("year","season","value")
  precip_sais <- as.data.frame(pivot_wider(precip_sais,names_from = season,values_from = value))
  
  if(liss!=1){
    precip_ann[,2] <- rollapply(precip_ann[,2],liss,mean,partial=T)
    precip_sais[,2:5] <- apply(precip_sais[,2:5],2,function(v) rollapply(v,liss,mean,partial=T))
  }
  
  # Import WP
  tt <- get.wp(nbdays = 1,start,end,agreg=T)
  
  tt_ann <- aggregate(tt,by=list(ann),function(v){sum(v==wp)/length(v)*100})

  tt_sais <- matrix(NA,length(unique(ann)),5)
  tt_sais[,1] <- unique(ann)
  
  for(i in 1:4){
  tmp <- tt
  tmp[!(substr(dates,6,7) %in% mois[[i]])] <- NA
  tt_sais[,i+1] <- aggregate(tmp,by=list(ann),function(v){v <- na.omit(v);sum(v==wp)/length(v)*100})[,2]
  }
  
  if(liss!=1){
    tt_sais[,-1] <- apply(tt_sais[,-1],2,function(v) rollapply(v,liss,mean,partial=T))
  }
  
  
  # Scatterplots
  png(filename = paste0("2_Travail/1_Past/Rresults/scatterplot.wp.precip/plot_",bv,"_wp=",wp,"_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
  par(mfrow=c(2,3),pty="s",mar=c(4,4,4,2))
  ylab <- "Precip Accumulation (mm)"
  
  # Annuel
  titre <- paste0("all"," (R = ",round(cor(tt_ann[,2],precip_ann[,2],use="pairwise.complete.obs"),2),")")
  plot(tt_ann[,2],precip_ann[,2],pch=19,xlab="WP occurrence (%)",ylab=ylab,main=titre)

  for(i in 1:4){
    titre <- paste0(sais.name[i]," (R = ",round(cor(tt_sais[,i+1],precip_sais[,i+1],use="pairwise.complete.obs"),2),")")
    plot(tt_sais[,i+1],precip_sais[,i+1],pch=19,xlab="WP occurrence (%)",ylab=ylab,main=titre)
  }
  
  graphics.off()
}

# Selection du WP journalier passe en fonction du nombre d'analogues choisis
select.wp.past <- function(rean,n.ana=115){
  
  # Import
  if(rean=="ERA20C") load("2_Travail/1_Past/ERA20C/compute_wp_past/wp_1900-01-01_1947-12-31_start.ana_1948-01-01_end.ana_2010-12-31.Rdata")
  
  # Traitement
  if(n.ana==1) tt_final <- unlist(lapply(tt,function(v)v[1,2]))
  if(n.ana>2) tt_final <- unlist(lapply(tt,function(v) as.numeric(names(which.max(table(v[1:n.ana,2]))))))
  
  # Ajout des WP connus
  wp <- get.wp(nbdays = 1,start = "1948-01-01",end = "2010-12-31",risk = F,agreg = T)
  tt_final <- c(tt_final,wp)
  
  if(rean=="ERA20C") save(tt_final,file=paste0("2_Travail/1_Past/ERA20C/compute_wp_past/wp_final_1900-01-01_1947-12-31_start.ana_1948-01-01_end.ana_2010-12-31_n.ana=",n.ana,".Rdata"))
  
}
  
# Verification de la representativite des types de temps choisis dans la reconstitution passee
verif.rep.past.wp <- function(rean,n.ana=115){
  
  # Import
  if(rean=="ERA20C") load("2_Travail/1_Past/ERA20C/compute_wp_past/wp_1900-01-01_1947-12-31_start.ana_1948-01-01_end.ana_2010-12-31.Rdata")
  
  # Traitement
  tmp <- lapply(tt,function(v) unname(sort(table(v[1:n.ana,2]),decreasing=T)[1])) # recuperation du nombre d'analogues du TT majoritaire
  pourcentage <- mean(unlist(tmp))/n.ana # pourcentage du TT majoritaire
  pourcentage
}

# Travail passé:

# - courbe de 1850 à 2010 de l'évolution des trois indicateurs pour les trois rénanalyses:
#   dans cette configuration, les indicateurs sont calculés en allant chercher les analogues sur
#   toute la periode disponible.

# - nuages de points des indicateurs VS dP coloriés par densité pour les trois reanalyses et 
#   pour des fenêtres de 50 ans (1851-1900; 1901-1950; 1951-2010) et de 25 ans
#   (1851-1875; 1876-1900; 1901-1925; 1926-1950; 1951-1975; 1976-2010).

# - construction cdf precip a partir des analogues au sens des indicateurs, 
#   un a un pour voir leur impact respectif, puis couple singnei/rsingnei: 
#   comment évolue la proba de precip annuelle (RP1an)? On va chercher les plus proches en celerite,
#   singularite et singularite relative dans la periode d'obs de 1950 a 2010 en gardant le fait d'avoir
#   calcule les analogues sur toute la periode de la reanalyse. Non: pour avoir une distribution de 
#   precip pour les analogues classiques, il faut aussi calculer les analogues uniquement sur la periode obs.
