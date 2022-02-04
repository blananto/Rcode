# Installation et chargement des packages utiles
#install.packages("ncdf4") # NetCDF

library("fields") # image.plot()
library("maptools") # wrld_simpl
library("ncdf4") # NetCDF
library("RColorBrewer") # brewer.pal

# Definition du repertoire de travail
#setwd("2_Travail/")

# Repere l'indice d'un netcdf associe a une date donne
date_to_number<-function(nc,day,rean,climat=NULL){
  which(getdays(nc,rean,climat)==day)
}

# Extrait les donnees d'un netcdf (k=1 pour 500hPa, k=2 pour 1000hPa) pour notre fenetre d'analogie
getdata<-function(k,day0,day1=day0,rean=NULL,climat=NULL,run=1,large_win=F,small_win=F,all=F,ssp=NULL,var="hgt"){
  
  # Import du netcdf
  nc <- load.nc(rean,var,climat,run,ssp)
  if(var=="hgt"){ nc <- nc[[k]]}
  
  # Definition du premier et dernier pas de temps a extraire
  num0<-date_to_number(nc,day0,rean,climat)
  num1<-date_to_number(nc,day1,rean,climat)
  N<-length(num0:num1)
  
  # Definition de la fenetre spatiale a extraire
  infowind<-getinfo_window(k,large_win,small_win,all,rean,var,climat,run,ssp)
  
  # Extraction
  data <- ncvar_get(nc = nc,varid=var,start=c(infowind[1,1],infowind[2,1],num0),count=c(infowind[1,2],infowind[2,2],N))
  nc_close(nc)
  data
}

# Defini une sequence de dates journalieres
getdates<-function(start="1950-01-01",end="2011-12-31"){
  as.character(seq(from=as.Date(start),to=as.Date(end),by="days"))
}

# Date associee a chaque pas de temps dans le fichier NetCDF
getdays<-function(nc,rean,climat=NULL){
  
  # Les reanalyses: origines differentes et pas de temps horaire
  if(is.null(climat)){
    
    # Origine
    if(substr(rean,1,4) == "20CR") orig <- "1800"
    if(substr(rean,1,4) == "NCEP") orig <- "1950"
    if(substr(rean,1,3) == "JRA") orig <- "1800"
    if(substr(rean,1,3) == "ERA") orig <- "1900"
    
    # Vecteur temps
    ti <- nc$dim$time$vals*3600
    substr(as.POSIXct(ti,origin=paste0(orig,'-01-01 00:00'),tz="GMT"),1,10) #format "YYYY-MM-DD"
    
    # Les modeles de climat: origine en 1850 et pas de temps journalier, aussi bien pour historical que pour ScenarioMIP
  }else{
    # Origine
    orig <- "1850"
    
    # Vecteur temps
    ti <- nc$dim$time$vals
    as.character(as.Date(ti,origin=as.Date(paste0(orig,"-01-01")))) #format "YYYY-MM-DD"
  }
}

# Repere les indices d'un netcdf pour definir la fenetre spatiale d'analogie
getinfo_window<-function(k,large_win = F,small_win=F,all=F,rean=NULL,var="hgt",climat=NULL,run=1,ssp=NULL){ 
  
  # large_win = T pour 1000hPa sur grande fenetre
  # small_win = T pour PWAT sur une fenetre reduite
  # all = T pour avoir tout le domaine du fichier netcdf
  
  # Import
  nc <- load.nc(rean,var,climat,run,ssp)
  if(var == "hgt"){ nc <- nc[[k]]}
  
  # Coord
  lon<-nc$dim$lon$vals
  lat<-nc$dim$lat$vals
  nc_close(nc)
  
  # Fenetre spatiale
  if(all){# Complete
    infolon <- c(1,length(lon)) # indice du 1er point de grille et nbre de points de grille
    infolat <- c(1,length(lat))
  }else{ # Incomplete
    
    c_lon<-6 # centre de la fenetre longitude
    c_lat<-44 # centre de la fenetre latitude
    
    if(k==1 | large_win) {# 500hPA
      d_lon<-32 # on prend 32 en longitude pour 500hPa
      d_lat<-16 # on prend 16 en latitude pour 500hPa
    }
    
    if (k==2 & !large_win){ # 1000 hPA
      d_lon<-16 # on prend 16 en longitude pour 1000hPa
      d_lat<-8  # on prend 8 en latitude pour 1000hPa
    }
    
    if(!is.null(rean)){
      if(k==1 & small_win & rean=="20CR") {# PWAT
        d_lon<-1
        c_lat<-46
        d_lat<-1
      }
      
      if(k==1 & small_win & rean=="ERA5") {# TCW
        c_lon<-6.25
        c_lat<-45.25
        d_lon<-1.5
        d_lat<-1.5
      }
    }
    
    deb_lon <- which.min(abs(lon - (c_lon-d_lon/2)))
    fin_lon <- which.min(abs(lon - (c_lon+d_lon/2)))
    len_lon <- fin_lon - deb_lon + 1
    
    deb_lat <- which.min(abs(lat - (c_lat-d_lat/2)))
    fin_lat <- which.min(abs(lat - (c_lat+d_lat/2)))
    len_lat <- fin_lat - deb_lat + 1
    
    infolon <- c(deb_lon,len_lon) # indice du 1er point de grille et nbre de points de grille
    infolat <- c(deb_lat,len_lat)
  }
  
  return(rbind(infolon,infolat))# matrice 2x2 
  
}

# Charge les fichiers NetCDF 500 hPa (et 1000 hPa pour certains) d'une reanayse ou d'un run d'un modele de climat
load.nc<-function(rean = NULL,var="hgt",climat=NULL,run=1,ssp=NULL){
  
  # Import de la reanalyse souhaitee
  if(!is.null(rean)){
    # Reanalyses brutes, pour differentes variables
    if(rean == "20CR" & var == "hgt"){
      member <- 1
      nc500<-nc_open(paste0("2_Travail/Data/Reanalysis/20CR/Membre_",member,"/20Crv2c_Membre_",member,"_HGT500_1851-2011_daily.nc"))
      nc1000<-nc_open(paste0("2_Travail/Data/Reanalysis/20CR/Membre_",member,"/20Crv2c_Membre_",member,"_HGT1000_1851-2011_daily.nc"))
    }
    
    if(rean == "20CR" & var == "pwat"){
      nc<-nc_open("2_Travail/Data/Reanalysis/20CR/PWAT/20Crv2_EnsembleMean_PWAT_1851-2014_daily.nc")
    }
    
    if(rean == "20CR" & var == "uwind"){
      nc <-nc_open("2_Travail/Data/Reanalysis/20CR/WIND/20Crv2c_EnsembleMean_UWIND_500_1851-2014_daily.nc")
    }
    
    if(rean == "20CR" & var == "vwind"){
      nc <-nc_open("2_Travail/Data/Reanalysis/20CR/WIND/20Crv2c_EnsembleMean_VWIND_500_1851-2014_daily.nc")
    }
    
    if(rean == "ERA20C"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/ERA20C/ERA20C_HGT500_1900_2010_daily.nc")
      nc1000<-nc_open("2_Travail/Data/Reanalysis/ERA20C/ERA20C_HGT1000_1900_2010_daily.nc")
    }
    
    if(rean == "ERA20C_18"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/ERA20C/ERA20C_HGT500_1900_2010_18h.nc")
      nc1000<-nc_open("2_Travail/Data/Reanalysis/ERA20C/ERA20C_HGT1000_1900_2010_18h.nc")
    }
    
    if(rean == "ERA20C_4daily"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/ERA20C/ERA20C_HGT500_1900_2010_4daily.nc")
      nc1000<-nc_open("2_Travail/Data/Reanalysis/ERA20C/ERA20C_HGT1000_1900_2010_4daily.nc")
    }
    
    if(rean == "NCEP"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/NCEP/NCEP_HGT500_1950_2010_daily.nc")
      nc1000<-nc_open("2_Travail/Data/Reanalysis/NCEP/NCEP_HGT1000_1950_2010_daily.nc")
    }
    
    if(rean == "ERA5" & var == "hgt"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/ERA5/HGT500/ERA5_HGT500_1950-01-01_2020-04-30_daily.nc")
      nc1000<-nc500
    }
    
    if(rean == "ERA5" & var == "tcw"){
      nc<-nc_open("2_Travail/Data/Reanalysis/ERA5/TCW/ERA5_TCW_1950_2019_daily.nc")
      names(nc$dim) <- c("lon","lat","time","bnds")
      names(nc$var)[2] <- "tcw"
      nc$var$tcw$name <- "tcw"
    }
    
    if(rean == "ERA5" & var == "uwind"){
      nc <-nc_open("2_Travail/Data/Reanalysis/ERA5/WIND/ERA5_UWIND_1950_2019_daily.nc")
      names(nc$dim) <- c("lon","lat","time","bnds")
      names(nc$var)[2] <- "uwind"
      nc$var$uwind$name <- "uwind"
    }
    
    if(rean == "ERA5" & var == "vwind"){
      nc <-nc_open("2_Travail/Data/Reanalysis/ERA5/WIND/ERA5_VWIND_1950_2019_daily.nc")
      names(nc$dim) <- c("lon","lat","time","bnds")
      names(nc$var)[2] <- "vwind"
      nc$var$vwind$name <- "vwind"
    }
    
    if(rean == "ERA40"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/ERA40/ERA40_HGT500_1957-09-01_2002-08-31_daily.nc")
      nc1000<-nc_open("2_Travail/Data/Reanalysis/ERA40/ERA40_HGT1000_1957-09-01_2002-08-31_daily.nc")
    }
    
    if(rean == "JRA55"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/JRA55/JRA55_HGT500_1958_2019_daily.nc")
      nc1000<-nc500
    }
    
    if(rean == "JRA55C"){
      nc500<-nc_open("2_Travail/Data/Reanalysis/JRA55C/JRA55C_HGT500_1972-11-01_2012-12-31_daily.nc")
      nc1000<-nc500
    }
    
    # Reanalyses interpolees sur les grilles des modeles de climat, pour la Z500
    if(rean == "20CR-cnrm"){
      nc500 <-nc_open("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Reanalysis/20CR/20Crv2c_Membre_1_HGT500_1851-2011_daily_remap_CNRM-CM6.nc")
      nc1000 <- nc500
    }
    
    if(rean == "20CR-ipsl"){
      nc500 <-nc_open("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Reanalysis/20CR/20Crv2c_Membre_1_HGT500_1851-2011_daily_remap_IPSL-CM6A.nc")
      nc1000 <- nc500
    }
    
    if(rean == "ERA20C-cnrm"){
      nc500 <-nc_open("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Reanalysis/ERA20C/ERA20C_HGT500_1900_2010_daily_remap_CNRM-CM6.nc")
      nc1000 <- nc500
    }
    
    if(rean == "ERA20C-ipsl"){
      nc500 <-nc_open("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Reanalysis/ERA20C/ERA20C_HGT500_1900_2010_daily_remap_IPSL-CM6A.nc")
      nc1000 <- nc500
    }
    
    if(rean == "ERA5-cnrm"){
      nc500 <-nc_open("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Reanalysis/ERA5/ERA5_HGT500_1950-01-01_2020-04-30_daily_remap_CNRM-CM6.nc")
      nc1000 <- nc500
    }
    
    if(rean == "ERA5-ipsl"){
      nc500 <-nc_open("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Reanalysis/ERA5/ERA5_HGT500_1950-01-01_2020-04-30_daily_remap_IPSL-CM6A.nc")
      nc1000 <- nc500
    }
    
    # Import du modele de climat souhaite
  }else{
    
    if(climat == "cnrm" & is.null(ssp)){
      nc500<-nc_open(paste0("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Climate_models/CNRM-CM6-1/Historical/zg500_AERday_CNRM-CM6-1_historical_r",as.character(run),"i1p1f2_gr_18500101-20141231_Europe.nc"))
      names(nc500$var)[3] <- "hgt"
      nc500$var$hgt$name <- "hgt"
      nc1000<-nc500
    }
    
    if(climat == "cnrm" & !is.null(ssp)){
      nc500<-nc_open(paste0("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Climate_models/CNRM-CM6-1/ScenarioMIP/zg500_AERday_CNRM-CM6-1_ssp",as.character(ssp),"_r",as.character(run),"i1p1f2_gr_20150101-21001231_Europe.nc"))
      names(nc500$var)[3] <- "hgt"
      nc500$var$hgt$name <- "hgt"
      nc1000<-nc500
    }
    
    if(climat == "ipsl" & is.null(ssp)){
      nc500<-nc_open(paste0("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Climate_models/IPSL-CM6A-LR/Historical/zg500_AERday_IPSL-CM6A-LR_historical_r",as.character(run),"i1p1f1_gr_18500101-20141231_Europe.nc"))
      names(nc500$var)[3] <- "hgt"
      nc500$var$hgt$name <- "hgt"
      nc1000<-nc500
    }
    
    if(climat == "ipsl" & !is.null(ssp)){
      nc500<-nc_open(paste0("9_Encadrement/Stage_Jules_Boulard_M2_2021/Travail/Data/Climate_models/IPSL-CM6A-LR/ScenarioMIP/zg500_AERday_IPSL-CM6A-LR_ssp",as.character(ssp),"_r",as.character(run),"i1p1f1_gr_20150101-21001231_Europe.nc"))
      names(nc500$var)[3] <- "hgt"
      nc500$var$hgt$name <- "hgt"
      nc1000<-nc500
    }
  }
  
  if(var=="hgt") nc<-list(nc500=nc500,nc1000=nc1000)
  nc
}

# Carte de geopotentiel d'un jour donne
map.geo <- function(date,rean="20CR",climat=NULL,run=1,k,nbdays=1,save=F,win=F,let=F,leg=T,iso=F,condens=F,ssp=NULL){
  
  # Dates concernees
  start <- date
  end <- as.character(as.Date(date)+nbdays-1)
  
  # Import des donnees
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean,climat = climat,run = run,large_win = F,small_win = F,all = T,var = "hgt",ssp = ssp)
  
  # Import lon/lat et fenetre
  nc <- load.nc(rean,climat,run,var="hgt",ssp)
  nc <- nc[[k]]
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  
  fen <- getinfo_window(k = k,rean=rean,climat=climat,run=run,var = "hgt")
  gc()
  
  # Parametres graphiques
  if(k==1){ 
    breaks <- seq(4850,6100,length.out = 12)
    N <- 11
    lab <- seq(4900,6100,200)
    lev <- seq(4900,6100,100)
  }else{
    breaks <- seq(-300,400,length.out = 8)
    N <- 7
    lab <- seq(-300,400,100)
    lev <- lab
  }
  
  # Carte
  if(save) {
    png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.geo/",date,"_k",k,"_",nbdays,"day.png"),
        width = ifelse(nbdays==3,800,600),height = 250,units = "px")
    layout(matrix(1:nbdays,1,nbdays))
  }
  
  par(pty="s")
  if(nbdays==3){
    #par(mar=c(6,4,6,4))
    par(mar=c(0,0,0,7))
  } 
  cex <- ifelse(nbdays==3,1.8,1.5)
  
  for(i in 1:nbdays){
    
    if(nbdays==1){
      z <- geo
    }else{z <- geo[,,i]}
    
    # Carte geopotentiel
    if(leg){
      if(!condens){
        image.plot(lon,lat,z,xlim=c(-15,25),ylim=c(25,65),
                   col=rev(brewer.pal(n = N, name = "RdBu")),
                   xlab="Longitude (°)",ylab="Latitude (°)",main=as.Date(date)+i-1,
                   legend.line=-2.3, cex.axis=cex, cex.lab=cex, cex.main=cex,
                   breaks = breaks,axis.args = list(at=lab,labels=as.character(lab),cex.axis=1.3))
      }else{
        image.plot(lon,lat,z,xlim=c(-15,25),ylim=c(25,65),
                   col=rev(brewer.pal(n = N, name = "RdBu")),
                   xlab="",ylab="",main="",xaxt="n",yaxt="n",
                   breaks = breaks,axis.args = list(at=lab,labels=as.character(lab)))
      }
    }else{
      if(!condens){
        image(lon,lat,z,xlim=c(-15,25),ylim=c(25,65),
              col=rev(brewer.pal(n = N, name = "RdBu")),
              xlab="Longitude (°)",ylab="Latitude (°)",main=as.Date(date)+i-1,
              cex.axis=cex, cex.lab=cex, cex.main=cex,
              breaks = breaks)
      }else{
        image(lon,lat,z,xlim=c(-15,25),ylim=c(25,65),
              col=rev(brewer.pal(n = N, name = "RdBu")),
              xlab="",ylab="",main="",xaxt="n",yaxt="n",
              breaks = breaks)
      }
    }
    
    # Isohypses
    if(iso) contour(x=lon,y=lat,z=z, levels=lev, drawlabels=F, lty=1, lwd=1, add=TRUE, col="black")
    
    # Carte monde
    data(wrld_simpl)
    plot(wrld_simpl, add = TRUE)
    
    # Region d'etude
    points(6,45,col="red",pch=19)
    
    # Fenetre d'analogie
    delta <- abs(lon[1]-lon[2])/2
    if(win) rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    
    # Lettre
    if(i==1 & let!=F) mtext(let, side=3, at=-30,line = 2,cex=1.5)
    
    # Date dans la carte
    if(condens) shadowtext(5,61,as.Date(date)+i-1,font=2,cex=1.7,col="black",bg="white",r=0.3)
    
    box()
  }
  
  if(save) graphics.off()
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