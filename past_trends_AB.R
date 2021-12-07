source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

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

# Carte de l'altitude des mailles ERA5
map.topo.ERA5 <- function(reg=c("small","medium","large")){
  
  start <- end <- "1950-01-01"
  
  # Import
  print("Import")
  topo <- getdata(k = 1,day0 = start,day1 = end,rean = "ERA5",var = "topo")
  fen <- getinfo_window(k = 1,rean = "ERA5",var = "topo")
  nc <- load.nc(rean = "ERA5",var = "topo")
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  
  lon <- lon[fen[1,1]:(fen[1,1]+fen[1,2]-1)]
  lat <- lat[fen[2,1]:(fen[2,1]+fen[2,2]-1)]
  gc()
  
  # Carte
  print("Carte")
  pal <- colorRampPalette(c("darkolivegreen2","chocolate4","white"))
  
  if(reg=="small"){xlim <- c(3,8);ylim <- c(43,47)}
  if(reg=="medium"){xlim <- c(-5,10);ylim <- c(42,52)}
  if(reg=="large"){xlim <- c(-10,18);ylim <- c(36,52)}
  
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.topo.ERA5/map_topo_ERA5_",reg,".png"),width = 9,height = 6,units = "in",res=600)
  par(pty="s",mar=c(4,4,2,1))
  image.plot(lon,lat,topo,xlim=xlim,ylim=ylim,
             col=pal(100),
             xlab="Longitude (°)",ylab="Latitude (°)",main="Altitude ERA5",legend.lab = "Altitude (m)",legend.line = 3.5)
  data("wrld_simpl")
  plot(wrld_simpl, add = TRUE)
  points(x = 5.73,y = 45.18,col = "red",pch=19)
  graphics.off()
}

# Carte des tendances des humidite specifiques
map.trend.sph <- function(z="850",wp=NULL,extr=F,rean="ERA5",start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  var <- paste0("sph",z)
  
  # Import
  print("Import")
  tmp <- getdata(k = 1,day0 = start,day1 = end,rean = rean,var = var,return.lonlat = T)
  sph <- tmp$data;lon <- tmp$lon;lat <- tmp$lat;rm(tmp)
  gc()
  
  # WP
  if(!is.null(wp)){
    tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    sph[,,tt!=wp] <- NA
    
    if(extr){
      ind.extr <- get.ind.max.flow(flow = wp,agreg = T,nbdays = 1,start = start,end = end,spazm = T,supseuil = T,nei = T)
      sph[,,!((1:length(dates)) %in% ind.extr)] <- NA
    }
  }
  
  # Calcul tendances
  ann <- substr(dates,1,4)
  an <- as.numeric(unique(ann))
  
  # Annuel
  print("Calcul Annuel")
  sph.moy <- apply(sph,1:2,mean,na.rm=T)
  sph.ann <- apply(sph,1:2,function(v){aggregate(v,by=list(ann),mean,na.rm=T)[,2]})
  trend.ann <- apply(sph.ann,2:3,function(v){lm(v~an)$coefficients[2]*10}) # unit/10 ans
  trend.ann.rel <- trend.ann/sph.moy*100 # passage en %/10 ans
  
  # Saisonnier
  if(!extr){
    print("Calcul Saisonnier")
    sais <- c("winter","spring","summer","autumn")
    sph.moy.sais <- trend.sais <- trend.sais.rel <- list()
    
    for(i in 1:length(sais)){
      print(paste0(i,"/",length(sais)))
      vec.sais <- rep(NA,length(dates))
      sea <- get.ind.season(sais = sais[i],start = start,end = end)
      vec.sais[sea$pos.season] <- i
      sph.moy.sais[[i]] <- apply(sph[,,sea$pos.season],1:2,mean,na.rm=T)
      sph.sais <- apply(sph,1:2,function(v){aggregate(v,by=list(ann,vec.sais),mean,na.rm=T)[,3]})
      trend.sais[[i]] <- apply(sph.sais,2:3,function(v){lm(v~an)$coefficients[2]*10})
      trend.sais.rel[[i]] <- trend.sais[[i]]/sph.moy.sais[[i]]*100
    }
  }
  
  # Cartes
  print("Cartes")
  N <- 9
  par(pty="s",mar=c(4,4,2,2))
  
  # Annuel
  # Champs moyen
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.sph/",var,"/mean_sph_annual",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
  #breaks <- seq(0,30,length.out = N+1)
  image.plot(lon,lat,sph.moy*1000,xlim=c(-10,22),ylim=c(36,52),
             col=brewer.pal(n = N, name = "BuGn"),#breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly mean",legend.lab = "sph (g.kg-1)",legend.line = 4.5,legend.mar = 8)
  data("wrld_simpl")
  plot(wrld_simpl, add = TRUE)
  graphics.off()
  
  # Tendance absolue
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.sph/",var,"/trend_sph_annual",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
  #breaks <- seq(0,0.45,length.out = N+1)
  image.plot(lon,lat,trend.ann*1000,xlim=c(-10,22),ylim=c(36,52),
             col=brewer.pal(n = N, name = "RdBu"),#breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly absolute trend",legend.lab = "Trend (g.kg-1/10y)",legend.line = 4.5,legend.mar = 8)
  plot(wrld_simpl, add = TRUE)
  graphics.off()
  
  # Saisonnier
  if(!extr){
    for(i in 1:length(sais)){
      
      print(paste0(i,"/",length(sais)))
      
      # Champs moyen
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.sph/",var,"/mean_sph_",sais[i],ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
      #breaks <- seq(0,30,length.out = N+1)
      image.plot(lon,lat,sph.moy.sais[[i]]*1000,xlim=c(-10,22),ylim=c(36,52),
                 col=brewer.pal(n = N, name = "BuGn"),#breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," mean"),legend.lab = "sph (g.kg-1)",legend.line = 4.5,legend.mar = 8)
      data("wrld_simpl")
      plot(wrld_simpl, add = TRUE)
      graphics.off()
      
      # Tendance absolue
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.sph/",var,"/trend_sph_",sais[i],ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
      #breaks <- seq(0,0.45,length.out = N+1)
      image.plot(lon,lat,trend.sais[[i]]*1000,xlim=c(-10,22),ylim=c(36,52),
                 col=brewer.pal(n = N, name = "RdBu"),#breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," absolute trend"),legend.lab = "Trend (g.kg-1/10y)",legend.line = 4.5,legend.mar = 8)
      plot(wrld_simpl, add = TRUE)
      graphics.off()
    }
  }
}

# Carte des tendances de l'humidite atmospherique
map.trend.tcw <- function(wp=NULL,extr=F,rean="ERA5",start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  
  # Import
  print("Import")
  tcw <- getdata(k = 1,day0 = start,day1 = end,rean = rean,var = "tcw")
  fen <- getinfo_window(k = 1,rean = rean,var = "tcw")
  nc <- load.nc(rean = rean,var = "tcw")
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  
  lon <- lon[fen[1,1]:(fen[1,1]+fen[1,2]-1)]
  lat <- lat[fen[2,1]:(fen[2,1]+fen[2,2]-1)]
  gc()
  
  # WP
  if(!is.null(wp)){
    tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    tcw[,,tt!=wp] <- NA
    
    if(extr){
      ind.extr <- get.ind.max.flow(flow = wp,agreg = T,nbdays = 1,start = start,end = end,spazm = T,supseuil = T,nei = T)
      tcw[,,!((1:length(dates)) %in% ind.extr)] <- NA
    }
  }
  
  # Calcul tendances
  ann <- substr(dates,1,4)
  an <- as.numeric(unique(ann))
  
  # Annuel
  print("Calcul Annuel")
  tcw.moy <- apply(tcw,1:2,mean,na.rm=T)
  tcw.ann <- apply(tcw,1:2,function(v){aggregate(v,by=list(ann),mean,na.rm=T)[,2]})
  trend.ann <- apply(tcw.ann,2:3,function(v){lm(v~an)$coefficients[2]*10}) # unit/10 ans
  trend.ann.rel <- trend.ann/tcw.moy*100 # passage en %/10 ans
  
  # Saisonnier
  if(!extr){
    print("Calcul Saisonnier")
    sais <- c("winter","spring","summer","autumn")
    tcw.moy.sais <- trend.sais <- trend.sais.rel <- list()
    
    for(i in 1:length(sais)){
      print(paste0(i,"/",length(sais)))
      vec.sais <- rep(NA,length(dates))
      sea <- get.ind.season(sais = sais[i],start = start,end = end)
      vec.sais[sea$pos.season] <- i
      tcw.moy.sais[[i]] <- apply(tcw[,,sea$pos.season],1:2,mean,na.rm=T)
      tcw.sais <- apply(tcw,1:2,function(v){aggregate(v,by=list(ann,vec.sais),mean,na.rm=T)[,3]})
      trend.sais[[i]] <- apply(tcw.sais,2:3,function(v){lm(v~an)$coefficients[2]*10})
      trend.sais.rel[[i]] <- trend.sais[[i]]/tcw.moy.sais[[i]]*100
    }
  }
  
  # Cartes
  print("Cartes")
  N <- 9
  par(pty="s",mar=c(4,4,2,0.5))
  
  # Annuel
  # Champs moyen
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.tcw/mean_tcw_annual",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
  breaks <- seq(0,30,length.out = N+1)
  image.plot(lon,lat,tcw.moy,xlim=c(-10,22),ylim=c(36,52),
             col=brewer.pal(n = N, name = "BuGn"),breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly mean",legend.lab = "TCW (mm)",legend.line = 3)
  data("wrld_simpl")
  plot(wrld_simpl, add = TRUE)
  graphics.off()
  
  # Tendance absolue
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.tcw/trend_tcw_annual",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
  breaks <- seq(0,0.45,length.out = N+1)
  image.plot(lon,lat,trend.ann,xlim=c(-10,22),ylim=c(36,52),
             col=brewer.pal(n = N, name = "BuGn"),breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly absolute trend",legend.lab = "Trend (mm/10y)",legend.line = 3)
  plot(wrld_simpl, add = TRUE)
  graphics.off()
  
  # Tendance relative
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.tcw/trend_tcw_annual_rel",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
  breaks <- seq(0,3,length.out = N+1)
  image.plot(lon,lat,trend.ann.rel,xlim=c(-10,22),ylim=c(36,52),
             col=brewer.pal(n = N, name = "BuGn"),breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly relative trend",legend.lab = "Trend (%/10y)",legend.line = 3)
  plot(wrld_simpl, add = TRUE)
  graphics.off()
  
  # Saisonnier
  if(!extr){
    for(i in 1:length(sais)){
      
      print(paste0(i,"/",length(sais)))
      
      # Champs moyen
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.tcw/mean_tcw_",sais[i],ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
      breaks <- seq(0,30,length.out = N+1)
      image.plot(lon,lat,tcw.moy.sais[[i]],xlim=c(-10,22),ylim=c(36,52),
                 col=brewer.pal(n = N, name = "BuGn"),breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," mean"),legend.lab = "TCW (mm)",legend.line = 3)
      data("wrld_simpl")
      plot(wrld_simpl, add = TRUE)
      graphics.off()
      
      # Tendance absolue
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.tcw/trend_tcw_",sais[i],ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
      breaks <- seq(0,0.45,length.out = N+1)
      image.plot(lon,lat,trend.sais[[i]],xlim=c(-10,22),ylim=c(36,52),
                 col=brewer.pal(n = N, name = "BuGn"),breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," absolute trend"),legend.lab = "Trend (mm/10y)",legend.line = 3)
      plot(wrld_simpl, add = TRUE)
      graphics.off()
      
      # Tendance relative
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.tcw/trend_tcw_",sais[i],"_rel",ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=600)
      breaks <- seq(0,3,length.out = N+1)
      image.plot(lon,lat,trend.sais.rel[[i]],xlim=c(-10,22),ylim=c(36,52),
                 col=brewer.pal(n = N, name = "BuGn"),breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," relative trend"),legend.lab = "Trend (%/10y)",legend.line = 3)
      plot(wrld_simpl, add = TRUE)
      graphics.off()
    }
  }
}

# Carte des tendances d'une variable choisie aux dates des precipitations extremes sur Isere/Drac
map.trend.var.extr <- function(var="vv700",sais="winter",bv="Isere",wp=1,rean="ERA5",start="1950-01-01",end="2017-12-31",reg="small"){
  
  dates <- getdates(start,end)
  ann <- substr(dates,1,4)
  
  # Region
  tmp <- get.region(reg = reg)
  xlim <- tmp$xlim; ylim <- tmp$ylim; rm(tmp)
  if(reg=="small") load("2_Travail/Data/Carto/SecteurHydro/SecteurHydro_lonlat.Rdata")
  
  # Import
  tmp <- getdata(k = 1,day0 = start,day1 = end,rean = rean,lim.lon = xlim,lim.lat = ylim,var = var,return.lonlat = T)
  data <- tmp$data;lon <- tmp$lon;lat <- tmp$lat;rm(tmp)
  if(substr(var,1,3)=="sph") data <- data*1000 # passage en g/kg
  gc()
  
  # Dates precip extremes
  extr <- get.ind.max.sais(sais = sais,wp = wp,nbdays = 1,start = start,end = end,bv = bv,spazm = T)
  data[,,!((1:length(dates)) %in% extr)] <- NA
  
  # Calcul moyenne et tendances
  data.moy <- apply(data,1:2,mean,na.rm=T)
  data.trend <- apply(data,1:2,function(v){lm(v~as.numeric(ann))$coefficients[2]*10}) # unit/10 ans
  
  # Carte moyenne et tendance
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.var.extr/map_extr_",var,"_",sais,"_wp",wp,"_",bv,"_",start,"_",end,"_",reg,".png"),width = ifelse(reg!="small",14,12),height = 6,units = "in",res=600)
  par(mfrow=c(1,2),mar=c(4,4,2,4))
  
  param <- get.param.map(field = data.moy,var = var,type = "Mean")
  image.plot(lon,lat,data.moy,xlim=xlim,ylim=ylim,breaks=param$breaks,col=param$col,
             xlab="Longitude (°)",ylab="Latitude (°)",main=param$main,legend.lab = param$leg,legend.line = 3,legend.mar = 12)
  points(x = 5.73,y = 45.18,pch=19,cex=1.5)
  if(reg!="small"){
    data("wrld_simpl")
    plot(wrld_simpl, add = TRUE)
  }else{
    for (i in 1:length(bv.list.lonlat)) {
      for(j in 1:length(bv.list.lonlat[[i]])){
        par(new=F)
        polygon(bv.list.lonlat[[i]][[j]])
      }
    }
  }
  
  param <- get.param.map(field = data.trend,var = var,type = "Trend")
  image.plot(lon,lat,data.trend,xlim=xlim,ylim=ylim,breaks=param$breaks,col=param$col,
             xlab="Longitude (°)",ylab="Latitude (°)",main=param$main,legend.lab = param$leg,legend.line = 3,legend.mar = 12)
  points(x = 5.73,y = 45.18,pch=19,cex=1.5)
  if(reg!="small"){
    data("wrld_simpl")
    plot(wrld_simpl, add = TRUE)
  }else{
    for (i in 1:length(bv.list.lonlat)) {
      for(j in 1:length(bv.list.lonlat[[i]])){
        par(new=F)
        polygon(bv.list.lonlat[[i]][[j]])
      }
    }
  }
  graphics.off()
}

# Carte des tendances des vitesses verticales
map.trend.vv <- function(z="850",wp=NULL,extr=F,rean="ERA5",start="1950-01-01",end="2017-12-31",reg="large"){
  
  dates <- getdates(start,end)
  var <- paste0("vv",z)
  
  # Import
  print("Import")
  vv <- getdata(k = 1,day0 = start,day1 = end,rean = rean,var = var)
  fen <- getinfo_window(k = 1,rean = rean,var = var)
  nc <- load.nc(rean = rean,var = var)
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  
  lon <- lon[fen[1,1]:(fen[1,1]+fen[1,2]-1)]
  lat <- lat[fen[2,1]:(fen[2,1]+fen[2,2]-1)]
  gc()
  
  # WP
  if(!is.null(wp)){
    tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    vv[,,tt!=wp] <- NA
    
    if(extr){
      ind.extr <- get.ind.max.flow(flow = wp,agreg = T,nbdays = 1,start = start,end = end,spazm = T,supseuil = T,nei = T)
      vv[,,!((1:length(dates)) %in% ind.extr)] <- NA
    }
  }
  
  # Calcul tendances
  ann <- substr(dates,1,4)
  an <- as.numeric(unique(ann))
  
  # Annuel
  print("Calcul Annuel")
  vv.moy <- apply(vv,1:2,mean,na.rm=T)
  vv.ann <- apply(vv,1:2,function(v){aggregate(v,by=list(ann),mean,na.rm=T)[,2]})
  trend.ann <- apply(vv.ann,2:3,function(v){lm(v~an)$coefficients[2]*10}) # unit/10 ans
  trend.ann.rel <- trend.ann/vv.moy*100 # passage en %/10 ans
  
  # Saisonnier
  if(!extr){
    print("Calcul Saisonnier")
    sais <- c("winter","spring","summer","autumn")
    vv.moy.sais <- trend.sais <- trend.sais.rel <- list()
    
    for(i in 1:length(sais)){
      print(paste0(i,"/",length(sais)))
      vec.sais <- rep(NA,length(dates))
      sea <- get.ind.season(sais = sais[i],start = start,end = end)
      vec.sais[sea$pos.season] <- i
      vv.moy.sais[[i]] <- apply(vv[,,sea$pos.season],1:2,mean,na.rm=T)
      vv.sais <- apply(vv,1:2,function(v){aggregate(v,by=list(ann,vec.sais),mean,na.rm=T)[,3]})
      trend.sais[[i]] <- apply(vv.sais,2:3,function(v){lm(v~an)$coefficients[2]*10})
      trend.sais.rel[[i]] <- trend.sais[[i]]/vv.moy.sais[[i]]*100
    }
  }
  
  # Cartes
  print("Cartes")
  N <- 9
  if(reg=="small"){xlim <- c(3,8);ylim <- c(43,47)}
  if(reg=="medium"){xlim <- c(-5,10);ylim <- c(42,52)}
  if(reg=="large"){xlim <- c(-10,18);ylim <- c(36,52)}
  par(pty="s",mar=c(4,4,2,0.5))
  
  # Annuel
  # Champs moyen
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.vv/",var,"/mean_vv_annual",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,ifelse(reg!="large",paste0("_",reg),""),".png"),width = 8,height = 6,units = "in",res=600)
  #breaks <- seq(0,30,length.out = N+1)
  image.plot(lon,lat,vv.moy,xlim=xlim,ylim=c(36,52),
             col=brewer.pal(n = N, name = "RdBu"),#breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly mean",legend.lab = "vv (Pa.s-1)",legend.line = 3)
  data("wrld_simpl")
  plot(wrld_simpl, add = TRUE)
  if(type!="large") points(x = 5.73,y = 45.18,col = "red",pch=19)
  graphics.off()
  
  # Tendance absolue
  png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.vv/",var,"/trend_vv_annual",ifelse(!is.null(wp),paste0("_wp=",wp),""),ifelse(extr,"_extr",""),"_",start,"_",end,ifelse(reg!="large",paste0("_",reg),""),".png"),width = 8,height = 6,units = "in",res=600)
  #breaks <- seq(0,0.45,length.out = N+1)
  image.plot(lon,lat,trend.ann,xlim=xlim,ylim=c(36,52),
             col=brewer.pal(n = N, name = "RdBu"),#breaks = breaks,
             xlab="Longitude (°)",ylab="Latitude (°)",main="Yearly absolute trend",legend.lab = "Trend (Pa.s-1/10y)",legend.line = 3)
  plot(wrld_simpl, add = TRUE)
  if(type!="large") points(x = 5.73,y = 45.18,col = "red",pch=19)
  graphics.off()
  
  # Saisonnier
  if(!extr){
    for(i in 1:length(sais)){
      
      print(paste0(i,"/",length(sais)))
      
      # Champs moyen
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.vv/",var,"/mean_vv_",sais[i],ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,ifelse(reg!="large",paste0("_",reg),""),".png"),width = 8,height = 6,units = "in",res=600)
      #breaks <- seq(0,30,length.out = N+1)
      image.plot(lon,lat,vv.moy.sais[[i]],xlim=xlim,ylim=c(36,52),
                 col=brewer.pal(n = N, name = "RdBu"),#breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," mean"),legend.lab = "vv (Pa.s-1)",legend.line = 3)
      data("wrld_simpl")
      plot(wrld_simpl, add = TRUE)
      if(type!="large") points(x = 5.73,y = 45.18,col = "red",pch=19)
      graphics.off()
      
      # Tendance absolue
      png(filename = paste0("2_Travail/1_Past/",rean,"/map.trend.vv/",var,"/trend_vv_",sais[i],ifelse(!is.null(wp),paste0("_wp=",wp),""),"_",start,"_",end,ifelse(reg!="large",paste0("_",reg),""),".png"),width = 8,height = 6,units = "in",res=600)
      #breaks <- seq(0,0.45,length.out = N+1)
      image.plot(lon,lat,trend.sais[[i]],xlim=xlim,ylim=c(36,52),
                 col=brewer.pal(n = N, name = "RdBu"),#breaks = breaks,
                 xlab="Longitude (°)",ylab="Latitude (°)",main=paste0(nam2str(sais[i])," absolute trend"),legend.lab = "Trend (Pa.s-1/10y)",legend.line = 3)
      plot(wrld_simpl, add = TRUE)
      if(type!="large") points(x = 5.73,y = 45.18,col = "red",pch=19)
      graphics.off()
    }
  }
}

# Trace l'evolution dans le temps de la latitude du jet saison et type de temps
plot.trend.lat.jet <- function(gamme=c(5450,5550),wp="all",start="1900-01-01",end="2010-12-31",rean,liss=1){
  
  # Import
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  load(paste0("2_Travail/1_Past/",rean,"/compute.lat.jet/weighted_mean_lat_jet_",start,"_",end,"_btw_",gamme[1],"_and_",gamme[2],".Rdata"))
  
  # Traitement
  if(wp!="all" & rean!="ERA5"){
    load(paste0("2_Travail/1_Past/",rean,"/compute_wp_past/wp_final_",start,"_1947-12-31_start.ana_1948-01-01_end.ana_2010-12-31_n.ana=1.Rdata"))
    lati[tt_final!=wp] <- NA
  }
  if(wp!="all"){
    tt_final <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
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
plot.trend.precip <- function(bv="Isere-seul",nbdays,spazm=F,start="1950-01-01",end="2019-12-31",nao=F,liss=1){
  
  # Import
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  ann <- as.numeric(substr(dates,1,4))
  precip <- get.precip(nbdays,start,end,bv,spazm=spazm)
  
  # Traitement
  cum_an <- aggregate(precip,by=list(ann),sum)
  max_an <- aggregate(precip,by=list(ann),max)
  
  sais <- c("winter","spring","summer","autumn")
  vec.sais <- rep(NA,length(dates))
  
  for(i in 1:4){
    sea <- get.ind.season(sais = sais[i],start = start,end = end)
    vec.sais[sea$pos.season] <- i
  }
  
  cum_sais <- aggregate(precip,by=list(ann,vec.sais),sum)
  colnames(cum_sais) <- c("year","season","value")
  cum_sais <- as.data.frame(pivot_wider(cum_sais,names_from = season,values_from = value))
  
  max_sais <- aggregate(precip,by=list(ann,vec.sais),max)
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
    for(i in 1:4) nao_sais[[i]] <- get.nao(start = ann[1],end = ann[length(ann)],sais = sais[i])
    if(liss!=1){
      nao_an[,2] <- rollapply(nao_an[,2],liss,mean,partial=T)
      nao_sais <- lapply(nao_sais,function(v) {v[,2] <- rollapply(v[,2],liss,mean,partial=T);return(v)})
    }
  }
  
  # Figures
  
  if(nbdays==1){
    # Cumul annuel
    png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_cum_an_",head(ann,1),"_",tail(ann,1),"_liss=",liss,ifelse(spazm,"_spazm",""),ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
    par(mar=c(4,4,2,3))
    plot(cum_an[,c(1,2)],type="n",xlab="Year",ylab="Precipitation (mm)",main="Year")
    grid()
    lines(cum_an[,c(1,2)],col="cornflowerblue",lwd=2)
    reg <- lm(cum_an[,2]~cum_an[,1])
    abline(reg,col="red",lwd=2)
    text(quantile(cum_an[,1],0.8),quantile(cum_an[,2],0.8),paste0("pvalue=",round(unname(summary(reg)$coefficients[,4][2]),3)),col="red",font=2)
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
      png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_cum_",sais[i],"_",head(ann,1),"_",tail(ann,1),"_liss=",liss,ifelse(spazm,"_spazm",""),ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
      par(mar=c(4,4,2,3))
      plot(cum_sais[,c(1,i+1)],type="n",xlab="Year",ylab="Precipitation (mm)",main=nam2str(sais[i]))
      grid()
      lines(cum_sais[,c(1,i+1)],col="cornflowerblue",lwd=2)
      reg <- lm(cum_sais[,i+1]~cum_sais[,1])
      abline(reg,col="red",lwd=2)
      text(quantile(cum_sais[,1],0.8),quantile(cum_sais[,i+1],0.8),paste0("pvalue=",round(unname(summary(reg)$coefficients[,4][2]),3)),col="red",font=2)
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
  
  # Max annuel
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_max_an_",nbdays,"days_",head(ann,1),"_",tail(ann,1),"_liss=",liss,ifelse(spazm,"_spazm",""),ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
  par(mar=c(4,4,2,3))
  plot(max_an[,c(1,2)],type="n",xlab="Year",ylab="Precipitation (mm)",main="Year")
  grid()
  lines(max_an[,c(1,2)],col="cornflowerblue",lwd=2)
  reg <- lm(max_an[,2]~max_an[,1])
  abline(reg,col="red",lwd=2)
  text(quantile(max_an[,1],0.8),quantile(max_an[,2],0.8),paste0("pvalue=",round(unname(summary(reg)$coefficients[,4][2]),3)),col="red",font=2)
  if(nao){
    par(new=T)
    plot(nao_an,col="grey",type="l",xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  graphics.off()
  
  # Max saisonniers
  for(i in 1:4){
    png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip/",bv,"_max_",sais[i],"_",nbdays,"days_",head(ann,1),"_",tail(ann,1),"_liss=",liss,ifelse(spazm,"_spazm",""),ifelse(nao,"_nao",""),".png"),width = 12,height = 9,units = "cm",res=300)
    par(mar=c(4,4,2,3))
    plot(max_sais[,c(1,i+1)],type="n",xlab="Year",ylab="Precipitation (mm)",main=nam2str(sais[i]))
    grid()
    lines(max_sais[,c(1,i+1)],col="cornflowerblue",lwd=2)
    reg <- lm(max_sais[,i+1]~max_sais[,1])
    abline(reg,col="red",lwd=2)
    text(quantile(max_sais[,1],0.8),quantile(max_sais[,i+1],0.8),paste0("pvalue=",round(unname(summary(reg)$coefficients[,4][2]),3)),col="red",font=2)
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

# Trace l'evolution du SPH a differentes altitudes en un point donne
plot.trend.sph.point <- function(pt.lon=6,pt.lat=45,rean="ERA5",start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  ann <- substr(dates,1,4)
  lev <- c("925","850","700","500")
  
  # Import des donnees de sph au point
  sph <- matrix(data = NA,nrow = length(dates),ncol = length(lev))
  for(i in 1:length(lev)){
    sph[,i] <- getdata(k = 1,day0 = start,day1 = end,rean = rean,pt.lon = pt.lon,pt.lat = pt.lat,var = paste0("sph",lev[i]))
  }
  
  # Traitement
  # Annuel
  print("Calcul Annuel")
  sph.ann <- aggregate(sph,by=list(ann),mean,na.rm=T)
  
  # Max annuel de precip Atlantique
  extr <- get.ind.max.flow(flow = 1,agreg = T,nbdays = 1,start = start,end = end,spazm = T,supseuil = F,nei = F)
  plot(as.Date(dates[extr]),sph[extr,1]*1000,type="l",ylim=c(0,12))
  grid()
  for(i in 1:length(lev)){
    lines(as.Date(dates[extr]),sph[extr,i]*1000,col=colo[i])
    abline(lm(sph[extr,i]*1000~as.Date(dates[extr])),col=colo[i])
  }
  
  # Saisonnier
  #if(!extr){
  #  print("Calcul Saisonnier")
  #  sais <- c("winter","spring","summer","autumn")
  #  sph.moy.sais <- trend.sais <- trend.sais.rel <- list()
  #  
  #  for(i in 1:length(sais)){
  #    print(paste0(i,"/",length(sais)))
  #    vec.sais <- rep(NA,length(dates))
  #    sea <- get.ind.season(sais = sais[i],start = start,end = end)
  #    vec.sais[sea$pos.season] <- i
  #    sph.moy.sais[[i]] <- apply(sph[,,sea$pos.season],1:2,mean,na.rm=T)
  #    sph.sais <- apply(sph,1:2,function(v){aggregate(v,by=list(ann,vec.sais),mean,na.rm=T)[,3]})
  #    trend.sais[[i]] <- apply(sph.sais,2:3,function(v){lm(v~an)$coefficients[2]*10})
  #    trend.sais.rel[[i]] <- trend.sais[[i]]/sph.moy.sais[[i]]*100
  #  }
  #}
  
  
  # Graphiques
  colo <- c("red","darkorange","darkgreen","darkblue")
  
  plot(as.numeric(sph.ann[,1]),sph.ann[,2]*1000,type="l",ylim=c(0.5,6.5))
  grid()
  for(i in 1:length(lev)){
    lines(as.numeric(sph.ann[,1]),sph.ann[,i+1]*1000,col=colo[i])
    reg <- lm(sph.ann[,i+1]*1000~as.numeric(sph.ann[,1]))
    abline(reg,col=colo[i])
    text(2010,tail(sph.ann[,i+1]*1000,1)*1.05,paste0(round(reg$coefficients[2]*10/mean(sph.ann[,i+1]*1000,na.rm=T)*100,4),"g/kg/10y"))
    print(paste0("pvalue=",summary(reg)$coefficients[,4][2]))
  }
  
  sph.ann <- t(sph.ann[,-1])
  plot(sph.ann[,1]*1000,c(925,850,700,500),type="l",xlim=c(0.5,6.5),ylim = rev(range(c(925,850,700,500))),log="y")
  for(i in 1:ncol(sph.ann)){
    lines(sph.ann[,i]*1000,c(925,850,700,500))
  }
}

# Trace l'évolution de l'occurrence des WP reconstitues dans le passe par saison
plot.trend.wp.past <- function(start="1900-01-01",end="2010-12-31",rean,liss=1,n.ana=115){
  
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
run.past.trends <- function(type=1){
  
  # plot.trend.precip
  if(type==1){
    bv <- c("Isere-seul","Drac-seul","Isere")
    nbdays=c(1,3)
    spazm=F
    start <- "1950-01-01"
    end <- "2019-12-31"
    
    for(i in 1:length(bv)){
      print(bv[i])
      for(j in 1:length(nbdays)){
        print(nbdays[j])
        plot.trend.precip(bv = bv[i],nbdays = nbdays[j],spazm = spazm,start = start,end = end,nao = F,liss = 1)
      }
    }
  }
  
  # map.trend.var.extr
  if(type==2){
    bv <- c("Isere","Isere-seul","Drac-seul")
    wp <- c(1,2)
    sais <- c("spring","autumn","winter")
    var <- c("rh700","rh850","vv500","vv700","vv850","vv925","sph500","sph700","sph850","sph925","tcw","t700","t850")
    reg <- c("small")#,"large")
    
    for(i in 1:length(bv)){
      print(bv[i])
      for(j in 1:length(wp)){
        print(wp[j])
        for(k in 1:length(sais)){
          print(sais[k])
          for(l in 1:length(var)){
            print(var[l])
            for(m in 1:length(reg)){
              print(reg[m])
              map.trend.var.extr(var = var[l],sais = sais[k],bv = bv[i],wp = wp[j],rean = "ERA5",start = "1950-01-01",end = "2017-12-31",reg = reg[m])
              gc()
            }
          }
        }
      }
    }
  }
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