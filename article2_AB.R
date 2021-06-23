source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Combinaison de plusieurs plot.trend.descr
combine.plot.trend.descr <- function(k,dist,liss=5,ana.comm=F,align=F,nao=F){
  
  # Indicateurs
  #ind <- c("cel","sing05","rsing05","dP")
  ind <- "dP"
  namind <- ifelse(length(ind)==1,ind,"alldescr")
  
  # Saisons
  #sais <- "year"
  sais <- c("winter","spring","summer","autumn")
  namsais <- ifelse(length(sais)==1,sais,"allsais")
  
  # Graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.descr/plot_combine_trend_",namind,"_",namsais,"_liss=",liss,ifelse(ana.comm,paste0("_ana_",dates.ana[[1]][1],"_",dates.ana[[1]][2]),""),ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 9,height = 6,units = "in",res=600)
  par(mfrow=c(2,2))
  for(i in 1:length(sais)){
    print(paste0("plot ",i,"/",length(sais)))
    leg <- ifelse(i==1,T,F)
    plot.trend.descr(descr = ind,k = 1,dist = dist,sais = sais[i],liss = liss,
                     ana.comm = ana.comm,align = align,nao = nao,leg = leg,save = F,type="season")
  }
  graphics.off()
}

# Comparaison des indicateurs calcules sur differentes periodes d'analogie
compare.descr.ana <- function(descr,k,dist,rean){
  
  start.end.rean <- get.start.end.rean(rean,"past","criteria")
  dates.rean <- getdates(start.end.rean[1],start.end.rean[2])
  start.end.ana <- c("1950-01-01","2010-12-31")
  
  # Import des 2 versions de l'indicateur
  ind.long <- get.descriptor(descriptor = descr,k = k,dist = dist,
                 nbdays = 1,start = start.end.rean[1],end = start.end.rean[2],standardize=F,
                 rean = rean,period = "past",start.ana = start.end.rean[1],end.ana = start.end.rean[2])
  
  ind.short <- get.descriptor(descriptor = descr,k = k,dist = dist,
                             nbdays = 1,start = start.end.rean[1],end = start.end.rean[2],standardize=F,
                             rean = rean,period = "past",start.ana = start.end.ana[1],end.ana = start.end.ana[2])
  
  # Rappel du nombre d'analogues selectionnes pour chaque periode
  n.ana.long <- round(0.005*length(getdates(start.end.rean[1],start.end.rean[2])),0)
  n.ana.short <- round(0.005*length(getdates(start.end.ana[1],start.end.ana[2])),0)
  print(paste0("Periode ",start.end.rean[1],"-",start.end.rean[2],": ",n.ana.long," analogues selectionnes"))
  print(paste0("Periode ",start.end.ana[1],"-",start.end.ana[2],": ",n.ana.short," analogues selectionnes"))
  
  # Scatterplot et correlation
  nam.long <- paste0(descr," - Analogs in ",substr(start.end.rean[1],1,4),"-",substr(start.end.rean[2],1,4))
  nam.short <- paste0(descr," - Analogs in ",substr(start.end.ana[1],1,4),"-",substr(start.end.ana[2],1,4))
  corr <- round(cor(ind.long,ind.short,use = "pairwise.complete.obs"),3)
  
  png(filename = paste0(get.dirstr(k,rean,"past"),"compare.descr.ana/scatterplot_",descr,"_ana_",
                        start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".png"),
      width = 6,height = 6,units = "in",res=600)
  par(pty="s")
  lim <- range(na.omit(ind.long),na.omit(ind.short))
  plot(ind.long,ind.short,xlim=lim,ylim=lim,main=paste0(nam2str(descr,whole=T)," (R=",corr,")"),xlab=nam.long,ylab=nam.short)
  grid();par(new=T)
  plot(ind.long,ind.short,xlim=lim,ylim=lim,main=paste0(nam2str(descr,whole=T)," (R=",corr,")"),xlab=nam.long,ylab=nam.short)
  abline(0,1,col="red")
  graphics.off()
  
  # Evolution (en annuel)
  ind.long.ann <- aggregate(ind.long,by=list(substr(dates.rean,1,4)),mean,na.rm=T)
  ind.short.ann <- aggregate(ind.short,by=list(substr(dates.rean,1,4)),mean,na.rm=T)
  lim <- range(ind.long.ann[,2],ind.short.ann[,2])
  pos.leg <- ifelse(rean=="ERA20C","topright","bottomright")
  
  png(filename = paste0(get.dirstr(k,rean,"past"),"compare.descr.ana/trend_",descr,"_ana",
                        start.end.rean[1],"_",start.end.rean[2],"_ana_",start.end.ana[1],"_",start.end.ana[2],".png"),
      width = 7,height = 4,units = "in",res=600)
  par(mar=c(4,4,3,0))
  plot(ind.long.ann,type="l",ylim=lim,main=nam2str(descr,whole=T),xlab="Year",ylab=descr)
  grid();par(new=T)
  plot(ind.long.ann,type="l",ylim=lim,main=nam2str(descr,whole=T),xlab="Year",ylab=descr)
  lines(ind.short.ann,col="red")
  legend(pos.leg,lty=1,col=c("black","red"),legend = c(nam.long,nam.short),bty="n",cex=0.8)
  graphics.off()
}

# Comparaison des tendances des indicateurs avec et sans nei (cel, sing, rsing)
compare.descr.nei <- function(descr=c("cel","celnei"),k,dist,rean){
  
  start.end.rean <- get.start.end.rean(rean,"past","criteria")
  dates.rean <- getdates(start.end.rean[1],start.end.rean[2])
  #start.end.ana <- c("1950-01-01","2010-12-31")
  start.end.ana <- start.end.rean
  
  # Import des indicateurs
  ind1 <- get.descriptor(descriptor = descr[1],k = k,dist = dist,
                              nbdays = 1,start = start.end.rean[1],end = start.end.rean[2],standardize=F,
                              rean = rean,period = "past",start.ana = start.end.ana[1],end.ana = start.end.ana[2])
  
  ind2 <- get.descriptor(descriptor = descr[2],k = k,dist = dist,
                         nbdays = 1,start = start.end.rean[1],end = start.end.rean[2],standardize=F,
                         rean = rean,period = "past",start.ana = start.end.ana[1],end.ana = start.end.ana[2])
  
  # Evolution (en annuel)
  ind1.ann <- aggregate(ind1,by=list(substr(dates.rean,1,4)),mean,na.rm=T)
  ind2.ann <- aggregate(ind2,by=list(substr(dates.rean,1,4)),mean,na.rm=T)
  lim <- range(ind1.ann[,2],ind2.ann[,2])
  pos.leg <- ifelse(rean=="ERA20C","topright","bottomright")
  
  png(filename = paste0(get.dirstr(k,rean,"past"),"compare.descr.nei/trend_",descr[1],"_",descr[2],"_ana_",
                        start.end.ana[1],"_",start.end.ana[2],".png"),width = 7,height = 4,units = "in",res=600)
  par(mar=c(4,4,3,0))
  plot(ind1.ann,type="l",ylim=lim,main=nam2str(descr[1],whole=T),xlab="Year",ylab=nam2str(descr[1],whole=T))
  grid();par(new=T)
  plot(ind1.ann,type="l",ylim=lim,main=nam2str(descr[1],whole=T),xlab="Year",ylab=nam2str(descr[1],whole=T))
  abline(lm(ind1.ann[,2]~as.numeric(ind1.ann[,1])))
  lines(ind2.ann,col="red")
  abline(lm(ind2.ann[,2]~as.numeric(ind2.ann[,1])),col="red")
  legend(pos.leg,lty=1,col=c("black","red"),legend = c(descr[1],descr[2]),bty="n",cex=0.8)
  graphics.off()
}

# Calcul des indicateurs
compute_criteria_past_par <-function(k,dist,rean,start="1851-01-01",end="2010-12-31",start.ana="1950-01-01",end.ana="2010-12-31",period="past",update=FALSE,ncores){
  
  # Dates utiles
  dates <- getdates(start,end) # toutes les dates
  dates.ana <- getdates(start.ana,end.ana) # dates dans lesquelles on va cherches les analogues
  start.end.dist <- get.start.end.rean(rean,"past","dist") # debut et fin du fichier distance (dist.list)
  N<-length(dates)
  
  # Import des distances sur la periode souhaitee
  print(paste0("Import des distances ",dist))
  dist.vec<-getdist(k,dist,start,end,rean,threeday=F,period) # si rev, on prend les nei en TWS
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  
  U<-c(0,(N-1):1); # U = 0, 22644, 22643, 22642, ...
  sU<-sapply(1:(N-1),function(x) sum(U[1:x])) # somme cumulee de U: on fait la somme de U[1], U[1:2], etc pour obtenir la position de la derniere distance qui separe chaque date
  gc()
  
  # Import du fichier existant si update, et definition des indicateurs a calculer
  if (update) {
    load(file=paste0(get.dirstr(k,rean,period),"compute_criteria/criteria_",dist,"_",rean,"_k",k,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".Rdata"))
    coln.new<-c("celnei","singnei","rsingnei")
  }
  if (!update) {
    coln.new<-c("cel","sing05","q05")
  }
  
  # Fenetre de recherche des analogues
  ind <- match(dates.ana,dates)
  n<-length(ind)
  
  # Parallelisation & Calcul des indicateurs
  print("Calcul des indicateurs")
  print(paste0("Parallelisation sur ",ncores, " coeurs"))
  
  outfile <- paste0(get.dirstr(k,rean,period),"compute_criteria/calcul.txt")
  print(paste0("Logfile for // loop : ",outfile))
  cl <- makeCluster(ncores, outfile=outfile) 
  registerDoParallel(cl)
  
  # Pour aller voir calcul.txt
  # cmd
  # cd dossier (tapper uniquement I: pour aller sur le DD externe); dir au lieu de ls
  # powershell Get-Content calcul.txt -Wait // sous linux: tail -f calcul.txt
  # CTRL+C deux fois pour fermer le .txt
  
  criteria.new <- foreach (i=1:N,.combine = rbind) %dopar%{
    
    source("2_Travail/1_Past/getdist4i.R", encoding = 'UTF-8')
    if (i %% 500==0) {print(i)}
    di<-getdist4i(i,dist.vec,N,sU)
    gc()
    
    soso<-sort(di,index.return=TRUE) # classement par plus petit score, et donne les positions
    soso$ix <- soso$ix[soso$ix %in% ind] # on ne garde que les plus proches faisant partie de la fenetre de recherche des analogues
    soso$ix <- soso$ix[soso$ix!=i] # on retire la journee concernee de ses propores analogues
    qi05<-di[soso$ix[(0.005*n)]] # quantile 0.5%
    idi05<-soso$ix[1:(0.005*n)] # recupere la position des 0.5% les plus proches
    
    tmp<-NULL
    
    for (cc in coln.new){
      # Celerite
      if (cc=="cel") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,di[i-1])}
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
    
    gc()
    tmp
  }
  
  stopCluster(cl)
  print(paste0("Fin calcul indicateurs a : ",Sys.time()))
  
  # Mise en forme criteria.new
  colnames(criteria.new) <- coln.new
  if(!update) {criteria.new <- criteria.new * (10^-9)}
  
  if (update) {
    criteria<-cbind(criteria,criteria.new)
  } else { criteria<-criteria.new}
  
  save(criteria,file=paste0(get.dirstr(k,rean,period),"compute_criteria/criteria_",dist,"_",rean,"_k",k,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".Rdata"))
}

# Calcul densite de points dans un plan
compute.density <- function(rean,k,descriptors,dist,sais="all",nbdays=3,start="1851-01-01",end="2010-12-31",period="past",start.ana="1851-01-01",end.ana="2010-12-31",quant=F,wp=F,agreg=F,spazm=F){
  
  # Import des descripteurs
  N <- length(getdates(start,end))-nbdays+1
  descr <- matrix(NA,N,length(descriptors))
  
  for(i in 1:length(descriptors)){
    descr[,i] <- get.descriptor(descriptor = descriptors[i],k = k,dist = dist[i],
                                nbdays = nbdays,start = start,end = end,standardize=T,
                                rean = rean,period = period,start.ana = start.ana,end.ana = end.ana)
  }
  
  # Traitement si cas paticulier (saison,quant,wp)
  if(sais!="all"){
    pos <- get.ind.season.past(sais = sais,start = start,end = end)
    pos <- pos[1:(length(pos)-nbdays+1)]
    descr <- descr[pos,]
  }
  
  if(quant){
    descr <- apply(descr,2,function(v) ecdf(v)(v)*100)
  }
  
  rad <- 0.5*sd(descr[,1]) # on prend un demi ecart type de rayon. Si WP, on calcule la densite selon le meme rayon, pour rester dans le referenciel de densite de tout le nuage
  
  if(wp!=F){
    tt <- get.wp(nbdays,start,end,risk=F,bv="Isere",agreg=agreg,spazm = spazm)
    descr <- descr[tt==wp,]
  }
  
  # Calcul du voisinage
  nb <- NULL
  for(i in 1:nrow(descr)){
    if (i %%50==0) print(i)
    count <- nn2(data = descr,query = t(descr[i,]),searchtype = "radius",radius = rad,k = nrow(descr)) # nombre de voisins dans le rayon
    nb[i] <- rowSums(count$nn.idx>0)-1
  }
  print(paste0("min density: ",min(nb)))
  print(paste0("max density: ",max(nb)))
  res <- list(nb = nb,descr1 = descr[,1],descr2 = descr[,2])
  save(res,file = paste0("2_Travail/1_Past/",rean,"/compute.density/nbnei_",
                         descriptors[1],"_",descriptors[2],ifelse(length(descriptors)==3,paste0("_",descriptors[3]),""),
                         ifelse(quant,"_quant",""),ifelse(wp!=F,paste0("_wp",wp),""),
                         ifelse(agreg,"_agreg",""),ifelse(spazm,"_spazm",""),"_mean",nbdays,"day_",start,"_",end,
                         ifelse(sais!="all",paste0("_",sais),""),"_ana_",start.ana,"_",end.ana,".Rdata"))
}

# Calcul de la latitude journaliere du jet
compute.lat.jet <- function(gamme=c(5450,5550),k,start="1851-01-01",end="2010-12-31",rean){
  
  # Recuperation longitude et latitude
  nc <- load.nc(rean)$nc500
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  
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
  
  gc()
  
  save(tt,file=paste0("2_Travail/1_Past/",rean,"/compute_wp_past/wp_",start,"_",end,"_start.ana_",start.ana,"_end.ana_",end.ana,".Rdata"))
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

# Sous-periodes de 30 ans
get.subperiod <- function(){
  
  sub.period <- list(
    c("1861-01-01","1890-12-31"),
    c("1891-01-01","1920-12-31"),
    c("1921-01-01","1950-12-31"),
    c("1951-01-01","1980-12-31"),
    c("1981-01-01","2010-12-31")
  )
  
  sub.period
}

# Difference altitude geopotentiel 1900-1930 et 1970-2000
map.diff.geo <- function(k,rean,sais){
  
  # Import des donnees
  dates.deb <- c("1900-01-01","1929-12-31")
  data.deb <- getdata(k = k,day0 = dates.deb[1],day1 = dates.deb[2],rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
  pos.deb <- get.ind.season.past(sais = sais,start = dates.deb[1],end = dates.deb[2])
  data.deb <- data.deb[,,pos.deb]
  
  dates.fin <- c("1970-01-01","1999-12-31")
  data.fin <- getdata(k = k,day0 = dates.fin[1],day1 = dates.fin[2],rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
  pos.fin <- get.ind.season.past(sais = sais,start = dates.fin[1],end = dates.fin[2])
  data.fin <- data.fin[,,pos.fin]
  
  # Moyennes et difference
  mean.deb <- apply(data.deb,1:2,mean)
  mean.fin <- apply(data.fin,1:2,mean)
  diff.geo <- mean.fin - mean.deb
  
  # Cartes
  nc <- load.nc(rean,var="hgt")
  nc <- nc[[k]]
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  fen <- getinfo_window(k = k,rean=rean,var = "hgt")
  nc_close(nc)
  
  png(filename = paste0("2_Travail/1_Past/Rresults/map.diff.geo/map_diff_k",k,"_",rean,"_",sais,".png"),width = 12,height = 6,units = "in",res=600)
  par(pty="s",mfrow=c(1,3),mar=c(4,6,4,6))
  
  #range(mean.deb,mean.fin)
  breaks <- seq(4850,6100,length.out = 12)
  N <- 11
  lab <- seq(4900,6100,200)
  
  # debut
  image.plot(lon,lat,mean.deb,xlim=c(-15,25),ylim=c(25,65),main="1900-1930",
             col=rev(brewer.pal(n = N, name = "RdBu")),
             xlab="Longitude (°)",ylab="Latitude (°)",
             legend.line=-2.3,breaks = breaks,
             axis.args = list(at=lab,labels=as.character(lab),cex.axis=1.3))
  
  data(wrld_simpl)
  plot(wrld_simpl, add = TRUE)
  rect(xleft = lon[fen[1,1]]-1,ybottom = lat[fen[2,1]]-1,xright = lon[fen[1,1]+fen[1,2]-1]+1,ytop = lat[fen[2,1]+fen[2,2]-1]+1,lwd=2)
  
  # fin
  image.plot(lon,lat,mean.fin,xlim=c(-15,25),ylim=c(25,65),main="1970-2000",
             col=rev(brewer.pal(n = N, name = "RdBu")),
             xlab="Longitude (°)",ylab="Latitude (°)",
             legend.line=-2.3,breaks = breaks,
             axis.args = list(at=lab,labels=as.character(lab),cex.axis=1.3))
  
  plot(wrld_simpl, add = TRUE)
  rect(xleft = lon[fen[1,1]]-1,ybottom = lat[fen[2,1]]-1,xright = lon[fen[1,1]+fen[1,2]-1]+1,ytop = lat[fen[2,1]+fen[2,2]-1]+1,lwd=2)
  
  # diff
  breaks <- seq(-50,50,length.out = 12)
  N <- 11
  lab <- seq(-50,50,10)
  
  image.plot(lon,lat,diff.geo,xlim=c(-15,25),ylim=c(25,65),main="1970-2000 - 1900-1930",
             col=rev(brewer.pal(n = N, name = "RdBu")),
             xlab="Longitude (°)",ylab="Latitude (°)",
             legend.line=-2.3,breaks = breaks,
             axis.args = list(at=lab,labels=as.character(lab),cex.axis=1.3))
  
  plot(wrld_simpl, add = TRUE)
  rect(xleft = lon[fen[1,1]]-1,ybottom = lat[fen[2,1]]-1,xright = lon[fen[1,1]+fen[1,2]-1]+1,ytop = lat[fen[2,1]+fen[2,2]-1]+1,lwd=2)
  
  graphics.off()
}

# Difference altitude geopotentiel 1900-1930 et 1970-2000 pour 2 reanalyses et les 4 saisons
map.diff.geo.combine <- function(k,rean=c("20CR-m1","ERA20C")){
  
  # Import des donnees
  dates.deb <- c("1900-01-01","1929-12-31")
  dates.fin <- c("1970-01-01","1999-12-31")
  season <- c("winter","spring","summer","autumn")
  
  # Figure
  png(filename = paste0("2_Travail/1_Past/Rresults/map.diff.geo/map_diff_combine_allsais_k",k,"_",rean[1],"_",rean[2],".png"),width = 8,height = 5,units = "in",res=600)
  layout(matrix(c(1:8,rep(9,4)),nrow = 3,ncol = 4,byrow = T),widths = c(1,1,1,1),heights = c(1,1,0.3))
  par(mar=c(1,1,2,1),pty="s")
  
  for(i in 1:length(rean)){
    print(rean[i])
    data.deb <- getdata(k = k,day0 = dates.deb[1],day1 = dates.deb[2],rean = rean[i],climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
    data.fin <- getdata(k = k,day0 = dates.fin[1],day1 = dates.fin[2],rean = rean[i],climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
    
    # Pour ajout du rectangle
    nc <- load.nc(rean[i],var="hgt")
    nc <- nc[[k]]
    lon <- nc$dim$lon$vals
    lat <- nc$dim$lat$vals
    fen <- getinfo_window(k = k,rean=rean[i],var = "hgt")
    nc_close(nc)
    
    for(j in 1:length(season)){
      print(season[j])
      pos.deb <- get.ind.season.past(sais = season[j],start = dates.deb[1],end = dates.deb[2])
      pos.fin <- get.ind.season.past(sais = season[j],start = dates.fin[1],end = dates.fin[2])
      
      data.deb.sea <- data.deb[,,pos.deb]
      data.fin.sea <- data.fin[,,pos.fin]
      
      # Moyennes et difference
      mean.deb <- apply(data.deb.sea,1:2,mean)
      mean.fin <- apply(data.fin.sea,1:2,mean)
      diff.geo <- mean.fin - mean.deb
      
      # Cartes
      breaks <- seq(-50,50,length.out = 12)
      N <- 11
      lab <- seq(-50,50,10)
      
      image(lon,lat,diff.geo,xlim=c(-15,25),ylim=c(25,65),main=paste0(rean[i]," - ",nam2str(season[j])),
                 col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
      
      plot(wrld_simpl, add = TRUE)
      rect(xleft = lon[fen[1,1]]-1,ybottom = lat[fen[2,1]]-1,xright = lon[fen[1,1]+fen[1,2]-1]+1,ytop = lat[fen[2,1]+fen[2,2]-1]+1,lwd=2)
    }
  }
  
  # Legende
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = lab,at =  seq(0,1,length.out = N),
              vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
  text(x = 1,y = 0.9,"1970-2000 minus 1900-1930 Geopotential Height (m)",cex=1.2,font=2)
  
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

# Trace l'evolution d'un indicateur, pour plusieurs reanalyses et par saison, avec NAO
plot.trend.descr <- function(descr,k,dist,sais="year",liss=5,ana.comm=F,align=F,nao=F,leg=T,save=F,type="descr"){
  
  # Reanalyses
  rean <- c(
    "20CR-m0",
    "20CR-m1",
    "20CR-m2",
    "ERA20C",
    "ERA5"
    #"NCEP",
    #"JRA55",
    #"ERA40",
    #"JRA55C"
    )
  
  colo <- c(
    "black",
    "lightgray",
    "darkgrey",
    "red",
    "royalblue"
    #"NCEP",
    #"JRA55",
    #"ERA40",
    #"JRA55C"
  )
  
  # Dates utiles
  dates <- list(length=length(rean))
  dat <- dates.ana <- dates
  
  for(i in 1:length(rean)){
    dates[[i]] <- get.start.end.rean(rean = rean[i],period="past",type = "criteria",k = k) # date debut fin rean
    dat[[i]] <- getdates(dates[[i]][1],dates[[i]][2]) # serie de dates
    dates.ana[[i]] <- dates[[i]] # dates ana
    if(ana.comm) dates.ana[[i]] <- c("1950-01-01","2010-12-31") # dates ana
  }
  
  # Import et traitement
  des <- list(length=length(rean))
  if(align) delta <- vector(length=length(rean))
  
  for(i in 1:length(rean)){
    
    # Import de l'indicateur
    ind <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,
                               start = dates[[i]][1],end = dates[[i]][2],
                               standardize = F,rean = rean[i],threeday = F,
                               period = "past",start.ana = dates.ana[[i]][1],
                               end.ana = dates.ana[[i]][2])
    
    # Traitement saisonnier
    pos <- get.ind.season(sais = sais,start = dates[[i]][1],end = dates[[i]][2])
    ind <- ind[pos$pos.season]
    ind[pos$pos.NA] <- NA
    dim(ind) <- c(pos$l.season,pos$n.season)
    rm(pos)
    
    ind <- apply(ind,2,mean,na.rm=T)
    year <- unique(substr(dat[[i]],1,4))
    if(sais=="winter") year <- year[-1] # car pas de premiere saison hiver
    des[[i]] <- data.frame(Year=year,Ind=ind)
    des[[i]][,1] <- as.character(des[[i]][,1])
    
    # Alignement et lissage
    if(align){
      delta[i] <- des[[i]][nrow(des[[i]]),2]
      des[[i]][,2] <- des[[i]][,2] - des[[i]][nrow(des[[i]]),2]
    }
    if(liss!=1){
      des[[i]][,2] <- rollapply(des[[i]][,2],liss,mean,partial=F,fill=NA)
    }
  }
  
  # Si NAO (non utilise)
  if(nao){
    naoi <- get.nao(start = "1865",end = "2010",sais = sais)
    nao.year <- naoi[,1]
    nao.ind  <- naoi[,2]
    nao.ind <- rollapply(nao.ind,liss,mean,partial=F,fill=NA)
    nao.ind <- c(rep(NA,nao.year[1] - as.numeric(dat[[1]][1])-1),nao.ind)
  }
  
  # Graphique
  # Parametres
  start.xaxis <- trunc(as.numeric(des[[1]][1,1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- trunc(as.numeric(des[[1]][nrow(des[[1]]),1])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  pos.axis <- match(xaxis,des[[1]][,1])
  if(is.na(pos.axis[1])) pos.axis[1] <- pos.axis[2]-10
  ylim <- range(unlist(lapply(des,function(v){range(v[,2],na.rm=T)})))
  if(align) {delta <- delta - delta[which(rean=="ERA5")]; ylim <- range(ylim,delta)} # ERA5 en reference
  main <- ifelse(type=="season",nam2str(sais),nam2str(descr,whole=T))
  
  # Initialisation
  if(save) png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.descr/plot_trend_",descr,"_",sais,"_liss=",liss,ifelse(ana.comm,paste0("_ana_",dates.ana[[1]][1],"_",dates.ana[[1]][2]),""),ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 7,height = 5,units = "in",res=600)
  par(mar=c(4,4,2,1))
  plot(des[[1]][,2],type="n",xaxt="n",ylim=ylim,xlab="Year",ylab=nam2str(descr,whole=T,unit = T),main=main)
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = pos.axis,labels = xaxis)
  abline(v = pos.axis,lty=3,col="grey")
  
  # Ajout courbes
  for(i in 1:length(rean)){
    lines(match(des[[i]][,1],des[[1]][,1]),des[[i]][,2],col=colo[i],lwd=2)
  }
  
  # Points delta si align
  if(align) points(rep(tail(pos.axis,1),length(delta)),delta,col=colo,pch=18,cex=1.5)
  
  # NAO
  if(nao){
    par(new=T)
    plot(nao.ind,pch=19,col="grey",type="l",lwd=2,xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  
  # Legende
  if(leg) legend("topleft",inset=.01,nam2str(rean),col=colo,lty=1,lwd=2,bty="n",cex=0.7)
  if(save) graphics.off()
}

# Trace l'evolution d'un indicateur conditionnellement au type de temps et a la saison
plot.trend.descr.cond <- function(descr,k,dist,rean,sais="all",type="mean",start="1950-01-01",end="2010-12-31",start.ana="1950-01-01",end.ana="2010-12-31",pvalue=T){
  
  dates <- getdates(start,end)
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = start,end = end,
                        standardize = F,rean = rean,threeday = F,period = "past",start.ana = start.ana,
                        end.ana = end.ana)
  
  # Import des types de temps
  wp <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Mise en forme
  pos <- get.ind.season(sais = sais,start = start,end = end)$pos.season
  non.pos <- which(!((1:length(dates)) %in% pos))
  
  tab <- data.frame(sais = dates,wp = wp,descr = des)
  tab$sais <- as.character(tab$sais)
  tab$sais[non.pos] <- "0"
  
  if(sais=="winter"){
    tmp <- which(substr(tab$sais,6,7)=="12")
    tab$sais[tmp] <- paste0(as.character(as.numeric(substr(tab$sais[tmp],1,4))+1),substr(tab$sais[tmp],5,10))
    rm(tmp)
  }
  
  tab$sais[pos] <- substr(tab$sais[pos],1,4)
  tab$sais <- as.integer(tab$sais)
  
  # Traitement
  if(type=="mean"){
  res <- aggregate(tab$descr,by=list(tab$sais,tab$wp),mean,na.rm=T)
  }
  
  if(type=="q10"){
  res <- aggregate(tab$descr,by=list(tab$sais,tab$wp),quantile,probs=0.1,na.rm=T)
  }
  
  if(type=="q90"){
  res <- aggregate(tab$descr,by=list(tab$sais,tab$wp),quantile,probs=0.9,na.rm=T)
  }
  
  res <- res[-which(res[,1]==0),]
  
  # Graphique
  wp.unique <- sort(unique(tab$wp))
  namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  png(filename = paste0(get.dirstr(k,rean,"past"),"plot.trend.descr.cond/plot_",descr,"_",sais,"_",type,ifelse(!pvalue,"_nopvalue",""),".png"),width = 8,height = 6,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(4,4,3,0.5))
  
  for(i in 1:length(wp.unique)){
    
    pos <- which(res[,2]==wp.unique[i])
    reg <- lm(res[pos,3]~res[pos,1])
    pval <- round(unname(summary(reg)$coefficients[,4][2]),3)
    main <- ifelse(pvalue,paste0(namflow[i]," (p.value = ",pval,")"),namflow[i])
    
    plot(res[pos,1],res[pos,3],type="l",xlab="Year",ylab=nam2str(descr,whole=T,unit=T),main=main)
    grid();par(new=T)
    plot(res[pos,1],res[pos,3],type="l",xlab="Year",ylab=nam2str(descr,whole=T,unit=T),main=main)
    abline(reg$coefficients[1],reg$coefficients[2],col="red")
  }
  
  graphics.off()
}

# Trace l'evolution d'un indicateur conditionnellement au type de temps et a la saison
plot.trend.descr.cond.regquant <- function(descr,k,dist,rean,sais="all",start="1950-01-01",end="2010-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  dates <- getdates(start,end)
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = start,end = end,
                        standardize = F,rean = rean,threeday = F,period = "past",start.ana = start.ana,
                        end.ana = end.ana)
  
  # Import des types de temps
  wp <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  wp.unique <- sort(unique(wp))
  namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  
  # Mise en forme
  tab <- data.frame(sais = dates,wp = wp,descr = des)
  pos <- get.ind.season(sais = sais,start = start,end = end)$pos.season
  tab <- tab[pos,]
  tab$sais <- as.character(tab$sais)
  
  if(sais=="winter"){
    tmp <- which(substr(tab$sais,6,7)=="12")
    tab$sais[tmp] <- paste0(as.character(as.numeric(substr(tab$sais[tmp],1,4))+1),substr(tab$sais[tmp],5,10))
    rm(tmp)
  }
  
  tab$sais <- substr(tab$sais,1,4)
  tab$sais <- as.integer(tab$sais)
  tab <- as.data.frame(tab)
  
  # Traitement
  res <- vector("list",length=length(wp.unique))
  
  for(i in 1:length(wp.unique)){
    res[[i]] <- rq(formula = descr~sais,tau = c(0.1,0.5,0.9),data = tab,subset = which(wp==wp.unique[i]))
  }
  
  # Graphique
  ylim <- range(tab$descr)
  
  png(filename = paste0(get.dirstr(k,rean,"past"),"plot.trend.descr.cond.regquant/plot_",descr,"_",sais,".png"),width = 8,height = 6,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(4,4,3,0.5))
  
  for(i in 1:length(wp.unique)){
    
    # Utiles
    pos <- tab$wp == wp.unique[i]
    coef <- res[[i]]$coefficients
    fitted <- res[[i]]$fitted.values
    all.year <- unique(res[[i]]$x[,2])
    summa <- summary(res[[i]])
    
    # Les points
    plot(tab$sais[pos],tab$descr[pos],ylim = ylim,pch=19,cex=0.3,xlab="Year",ylab=nam2str(descr,whole=T,unit=T),main=namflow[i])
    grid();par(new=T)
    plot(tab$sais[pos],tab$descr[pos],ylim = ylim,pch=19,cex=0.3,xlab="Year",ylab=nam2str(descr,whole=T,unit=T),main=namflow[i])
    
    # Les regressions
    abline(coef[,1],col="red",lwd=2)
    abline(coef[,2],col="blue",lwd=2)
    abline(coef[,3],col="red",lwd=2)
    
    # Les pvalue
    if(ncol(summa[[1]]$coefficients) == 4){
      shadowtext(x = round(quantile(all.year,0.9),0),y = tail(fitted[,1],1)-0.05*diff(ylim),labels = paste0("pval = ",round(summa[[1]]$coefficients[2,4],4)),font = 2,cex = 0.8,col = "red",bg = "white",r = 0.2)
      shadowtext(x = round(quantile(all.year,0.9),0),y = tail(fitted[,2],1)-0.05*diff(ylim),labels = paste0("pval = ",round(summa[[2]]$coefficients[2,4],4)),font = 2,cex = 0.8,col = "blue",bg = "white",r = 0.2)
      shadowtext(x = round(quantile(all.year,0.9),0),y = tail(fitted[,3],1)-0.05*diff(ylim),labels = paste0("pval = ",round(summa[[3]]$coefficients[2,4],4)),font = 2,cex = 0.8,col = "red",bg = "white",r = 0.2)
    }
  }
  
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

# Trace la relation entre la singularite
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
  
# Run des fonctions
run.past <- function(type=1){
  
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
  
  # scatterplot.descr.wp
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
  
  # compute.density
  if(type==6){
    descr <- list(
      c("cel","dP"),
      c("sing05","dP"),
      c("rsing05","dP")
    )
    
    sub.period <- get.subperiod()
    sais <- c("all","winter","autumn","spring","summer")
    
    for(i in 1:length(descr)){
      print(descr[[i]])
      for(j in 1:length(sub.period)){
        print(sub.period[[j]])
        for(k in 1:length(sais)){
        compute.density(rean = "20CR",k = 1,descriptors = descr[[i]],dist = c("TWS","TWS"),sais = sais[k],
                        nbdays = 3,start = sub.period[[j]][1],end = sub.period[[j]][2],
                        period = "past",start.ana = "1851-01-01",end.ana = "2010-12-31")
        }
      }
    }
  }
  
  # scatterplot.descr.wp
  if(type==7){
    descr <- list(
      c("cel","dP"),
      c("sing05","dP"),
      c("rsing05","dP")
    )

    sais <- c("winter","autumn","spring","summer")
    
    for(i in 1:length(descr)){
      print(descr[[i]])
      for(j in 1:length(sais)){
        print(sais[j])
        scatterplot.descr.subperiod(rean = "20CR",k = 1,descriptors = descr[[i]],dist = c("TWS","TWS"),
                                    sais = sais[j],nbdays = 3,start = "1851-01-01",end = "2010-12-31",
                                    period = "past",start.ana = "1851-01-01",end.ana = "2010-12-31")
      }
    }
  }
  
  # plot.trend.descr.cond
  if(type==8){
    descr <- c("cel","dP","sing05","rsing05")
    sais <- c("winter","autumn","spring","summer")
    type <- c("mean","q10","q90")
    
    for(i in 1:length(descr)){
      print(descr[i])
      for(j in 1:length(sais)){
        print(sais[j])
        for(k in 1:length(type)){
          print(type[k])
          plot.trend.descr.cond(descr = descr[i],k = 1,dist = "TWS",rean = "ERA5",sais = sais[j],type = type[k])
        }
      }
    }
  }
  
  # plot.trend.descr.cond.regquant
  if(type==9){
    descr <- c("cel","dP","sing05","rsing05")
    sais <- c("winter","autumn","spring","summer")
    
    for(i in 1:length(descr)){
      print(descr[i])
      for(j in 1:length(sais)){
        print(sais[j])
         plot.trend.descr.cond.regquant(descr = descr[i],k = 1,dist = "TWS",rean = "ERA5",sais = sais[j])
      }
    }
  }
  
  # compare.descr.ana
  if(type==10){
    descr <- c("cel")
    rean <- c("20CR-m0","20CR-m1","20CR-m2","ERA20C")
    
    for(i in 1:length(descr)){
      print(descr[i])
      for(j in 1:length(rean)){
        print(rean[j])
        compare.descr.ana(descr = descr[i],k = 1,dist = "TWS",rean = rean[j])
      }
    }
  }
  
  # compare.descr.nei
  if(type==11){
    descr <- list(c("cel","celnei"),c("sing05","singnei"),c("rsing05","rsingnei"))
    rean <- c("20CR-m0","20CR-m1","20CR-m2","ERA20C","ERA5")
    
    for(i in 1:length(descr)){
      print(descr[[i]])
      for(j in 1:length(rean)){
        print(rean[j])
        compare.descr.nei(descr = descr[[i]],k = 1,dist = "TWS",rean = rean[j])
      }
    }
  }
  
  # map.diff.geo
  if(type==12){
    rean <- c("20CR-m1","ERA20C")
    sais <- c("year","winter","spring","summer","autumn")
    
    for(i in 1:length(rean)){
      print(rean[i])
      for(j in 1:length(sais)){
        print(sais[j])
        map.diff.geo(k = 1,rean = rean[i],sais = sais[j])
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

# Scatterplot de deux descripteurs pour une sous période de 30 ans, colorie par densite de points
scatterplot.descr.subperiod <- function(rean,k,descriptors,dist,sais="all",nbdays=3,start="1851-01-01",end="2010-12-31",period="past",start.ana="1851-01-01",end.ana="2010-12-31"){
  
  dates <- getdates(start,end)
  
  # Import des densites et indicateurs
  sub.period <- get.subperiod()
  res.sub.period <- vector("list",length = length(sub.period))
  
  for(i in 1:length(sub.period)){
    load(paste0("2_Travail/1_Past/",rean,"/compute.density/nbnei_",
                descriptors[1],"_",descriptors[2],ifelse(length(descriptors)==3,paste0("_",descriptors[3]),""),
                "_mean",nbdays,"day_",sub.period[[i]][1],"_",sub.period[[i]][2],ifelse(sais!="all",paste0("_",sais),""),"_ana_",start.ana,"_",end.ana,".Rdata"))
    res.sub.period[[i]] <- res
  }
  
  # Scatterplots
  descr1.range <- range(lapply(res.sub.period,function(v){range(v$descr1)}))
  descr2.range <- range(lapply(res.sub.period,function(v){range(v$descr2)}))
  nb.range <- range(lapply(res.sub.period,function(v){range(v$nb)}))
  par(pty="s")
  
  for(i in 1:length(sub.period)){
    png(filename = paste0("2_Travail/1_Past/20CR/scatterplot.descr.subperiod/scatterplot_",descriptors[1],"_",descriptors[2],ifelse(sais!="all",paste0("_",sais),""),"_",sub.period[[i]][1],"_",sub.period[[i]][2],"_ana_",start.ana,"_",end.ana,".png"),
        width = 6,height = 6,units = "in",res = 600)
    plot(res.sub.period[[i]]$descr1,
         res.sub.period[[i]]$descr2,
         col=getcol(res.sub.period[[i]]$nb,range = nb.range),
         xlab="",
         ylab="",
         xlim=c((mean(descr1.range)-(descr1.range[2]-descr1.range[1])/2),(mean(descr1.range)+(descr1.range[2]-descr1.range[1])*1.4/2)),
         ylim=c((mean(descr2.range)-(descr2.range[2]-descr2.range[1])/2),(mean(descr2.range)+(descr2.range[2]-descr2.range[1])*1.4/2)),
         main=paste0(substr(sub.period[[i]][1],1,4)," - ",substr(sub.period[[i]][2],1,4)),cex.axis=1.2,cex.main=1.5,xaxt="n",yaxt="n")
    grid();par(new=T)
    plot(res.sub.period[[i]]$descr1,
         res.sub.period[[i]]$descr2,
         col=getcol(res.sub.period[[i]]$nb,range = nb.range),
         xlab="",
         ylab="",
         xlim=c((mean(descr1.range)-(descr1.range[2]-descr1.range[1])/2),(mean(descr1.range)+(descr1.range[2]-descr1.range[1])*1.4/2)),
         ylim=c((mean(descr2.range)-(descr2.range[2]-descr2.range[1])/2),(mean(descr2.range)+(descr2.range[2]-descr2.range[1])*1.4/2)),
         main=paste0(substr(sub.period[[i]][1],1,4)," - ",substr(sub.period[[i]][2],1,4)),cex.axis=1.2,cex.main=1.5,xaxt="n",yaxt="n")
    
    title(xlab=nam2str(descriptors[1],unit=T), line=2.5, cex.lab=1.2)
    title(ylab=nam2str(descriptors[2],unit=T), line=2.5, cex.lab=1.2)
    axis(1);axis(2)
    addscale(vec = c(res.sub.period[[i]]$nb,nb.range),r=0)
    text(x=mean(descr1.range)+(descr1.range[2]-descr1.range[1])*1.1/2,
         y=mean(descr2.range)+(descr2.range[2]-descr2.range[1])*1.3/2,
         paste0(round(min(res.sub.period[[i]]$nb,na.rm=T),2),"-",round(max(res.sub.period[[i]]$nb,na.rm=T),2)))
    graphics.off()
  }
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
