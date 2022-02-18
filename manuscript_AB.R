library(raster) # database France
library(mapproj) # database Europe

# setwd("../../")
# Site utile: https://dimension.usherbrooke.ca/pages/32

# Combinaison de plot.descr.flood
combine.plot.descr.flood <- function(nbdays=1,rean="20CR-m1",start = "1851-01-01",end = "2010-12-31"){
  
  png(filename = paste0(get.dirstr(k = 1,rean = rean,period = "past"),"plot.descr.flood/combine_plot_descr_flood_",nbdays,"day_",start,"_",end,".png"),width = 8,height = 5,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(2,4,2,1))
  season <- c("winter","spring","summer","autumn")
  for(i in 1:length(season)){
    plot.descr.flood(sais = season[i],nbdays = nbdays,rean = rean,start = start,end = end)
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.precip.gev.wp
combine.plot.trend.precip.gev.wp <- function(bv=c("Isere-seul","Drac-seul"),nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31",rain=F){
  
  # graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip.gev.wp/plot_trend_precip_gev_wp_",ifelse(rain,"rain_",""),bv[1],"_",bv[2],"_",nbdays,"day_",ifelse(spazm,"spazm_",""),start,"_",end,".png"),width = 10,height = 6,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,5),2:6,rep(7,5),8:12),nrow = 4,ncol = 5,byrow = T),widths = c(0.1,1,1,1,1),heights = c(0.1,1,0.1,1))
  
  for(i in 1:length(bv)){
    
    # Titre
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
    text(1,1,nam2str(bv[i]),font=2,cex=2.5)
    
    # Axe
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,ifelse(rain,"Rain (mm/day)","Precipitation (mm/day)"),srt=90,cex=1.3)
    
    # Plot
    par(mar=c(3,2,2,0))
    plot.trend.precip.gev.wp(bv = bv[i],nbdays = nbdays,spazm = spazm,start = start,end = end,leg = ifelse(i==1,T,F),rain = rain)
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.temp
combine.plot.trend.temp <- function(bv=c("Isere-seul","Drac-seul"),nbdays=1,wp="all",start="1950-01-01",end="2017-12-31",rain=F){
  
  # graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.temp/plot_trend_temp_",ifelse(rain,"rain_",""),bv[1],"_",bv[2],"_wp",wp,"_",nbdays,"day_",start,"_",end,".png"),width = 10,height = 6,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,5),2:6,rep(7,5),8:12),nrow = 4,ncol = 5,byrow = T),widths = c(0.1,1,1,1,1),heights = c(0.1,1,0.1,1))
  
  for(i in 1:length(bv)){
    
    # Titre
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
    text(1,1,nam2str(bv[i]),font=2,cex=2.5)
    
    # Axe
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,"Temperature (°C)",srt=90,cex=1.3)
    
    # Plot
    par(mar=c(3,2,2,0))
    plot.trend.temp(bv = bv[i],nbdays = nbdays,wp = wp,start = start,end = end,leg = ifelse(i==1,T,F),rain = rain)
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.temp.wp
combine.plot.trend.temp.wp <- function(wp=c(1,2),bv=c("Isere-seul","Drac-seul"),nbdays=1,start="1950-01-01",end="2017-12-31"){
  
  # graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.temp.wp/plot_trend_temp_wp",bv[1],"_",bv[2],"_wp",paste(wp,collapse="_wp"),"_",nbdays,"day_",start,"_",end,".png"),width = 10,height = 6,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,5),2:6,rep(7,5),8:12),nrow = 4,ncol = 5,byrow = T),widths = c(0.1,1,1,1,1),heights = c(0.1,1,0.1,1))
  
  for(i in 1:length(bv)){
    
    # Titre
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
    text(1,1,nam2str(bv[i]),font=2,cex=2.5)
    
    # Axe
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,"Temperature (°C)",srt=90,cex=1.3)
    
    # Plot
    par(mar=c(3,2,2,0))
    plot.trend.temp.wp(bv = bv[i],nbdays = nbdays,wp = wp,start = start,end = end,leg = ifelse(i==1,T,F))
  }
  graphics.off()
}
  
# Comparaison de MPD a la distance des centres d'action
compare.dP.dist <-  function(k,rean,nbdays,start="1950-01-01",end="2017-12-31"){
  
  # Import des indicateurs
  dP <- get.descriptor(descriptor = "dP",k = k,dist = "TWS",nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,threeday = F,desais = F,period = "past",start.ana = "1950-01-01",end.ana = "2010-12-31")
  dP_grad <- get.dP.grad(k = k,nbdays = nbdays,start = start,end = end,rean = rean)
  dP.grad <- dP_grad[,1]
  di <- dP_grad[,2]
  lon.min <- dP_grad[,3]
  lon.max <- dP_grad[,4]
  lat.min <- dP_grad[,5]
  lat.max <- dP_grad[,6]
  rm(dP_grad)
  
  # Traitement: position des centres d'action au bord de la fenetre
  pos.min <- unique(sort(c(which(lon.min %in% c(-10,22)),which(lat.min %in% c(36,52)))))
  pos.max <- unique(sort(c(which(lon.max %in% c(-10,22)),which(lat.max %in% c(36,52)))))
  pos <- intersect(pos.min,pos.max)
  
  classes <- seq(0,800,100)
  per <- rep(NA,length(classes)-1)
  nb.class <- per
  
  for(i in 1:(length(classes)-1)){
    pos.class <- which(dP>classes[i] & dP<classes[i+1])
    per[i] <- length(intersect(pos.class,pos))/length(pos.class)*100
    nb.class[i] <- round(length(pos.class)/length(dP)*100,1)
    
  }

  png(filename = paste0(get.dirstr(k,rean,period="present"),"compare.dP.dist/compare_dP_dist_k",k,"_mean",nbdays,"day_",rean,"_",start,"_",end,".png"),width = 6,height = 4,units = "in",res = 600)
  par(mar=c(4,5,1,0))
  barplot(height = per,space=0,ylim=c(0,100),xlab="MPD (m)",ylab="Percentage of Z500 min/max \n on the contour of the domain (%)",col="cornflowerblue",font.lab=2)
  grid();par(new=T)
  barplot(height = per,space=0,ylim=c(0,100),col="cornflowerblue")
  axis(side = 1,at = 0:8,labels = classes)
  text(x = 0.5:7.5,y=5,paste0(nb.class,"%"),font=2)
  box()
  graphics.off()
}

# Comparaison de MPD a MPD version gradient
compare.dP.grad <- function(k,rean,nbdays,start="1950-01-01",end="2017-12-31",explain=F,qua=F,example=F){
  
  # Import des indicateurs
  dP <- get.descriptor(descriptor = "dP",k = k,dist = "TWS",nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,threeday = F,desais = F,period = "past",start.ana = "1950-01-01",end.ana = "2010-12-31")
  dP_grad <- get.dP.grad(k = k,nbdays = nbdays,start = start,end = end,rean = rean)
  dP.grad <- dP_grad[,1]
  di <- dP_grad[,2]
  lon.min <- dP_grad[,3]
  lon.max <- dP_grad[,4]
  rm(dP_grad)
  
  # Traitement et Graphique Scatterplot
  corr <- cor(dP,dP.grad)
  if(explain){
    diag <- sqrt(32^2+16^2)
    larg <- 16
    pos.diag <- which(di==diag)
    pos.larg <- which(di==larg)
  }
  
  png(filename = paste0(get.dirstr(k,rean,period="present"),"compare.dP.grad/compare_dP_dPgrad_k",k,"_mean",nbdays,"day_",rean,"_",start,"_",end,ifelse(explain,"_explain",""),ifelse(qua,"_qua",""),ifelse(example,"_example",""),".png"),width = 5,height = 5,units = "in",res = 600)
  par(pty="s",mar=c(4,4,0,0))
  plot(dP,dP.grad,xlab="MPD (m)",ylab="MPD_norm (m/°)",pch=19,cex=0.7)
  grid();par(new=T)
  plot(dP,dP.grad,xlab="MPD (m)",ylab="MPD_norm (m/°)",pch=19,cex=0.7)
  text(200,40,bquote("r"^2*~"="~.(round(corr^2,2))),cex=1.2,font=2)
  
  if(explain){
    points(dP[pos.diag],dP.grad[pos.diag],col="red",pch=19,cex=0.7)
    points(dP[pos.larg],dP.grad[pos.larg],col="blue",pch=19,cex=0.7)
  }
  if(qua){
    abline(h=quantile(dP.grad,probs=0.75),col="red",lwd=2)
    abline(v=quantile(dP,probs=0.75),col="red",lwd=2)
    text(450,quantile(dP.grad,probs=0.8),"q75",col="red",font=2)
  }
  if(example){
    points(dP[10574],dP.grad[10574],col="darkgreen",pch=19)
    text(dP[10574]+30,dP.grad[10574]+1,"1",col="darkgreen",font=2,cex=1.5)
    points(dP[11695],dP.grad[11695],col="darkgreen",pch=19)
    text(dP[11695]+30,dP.grad[11695]+1,"2",col="darkgreen",font=2,cex=1.5)
  }
  graphics.off()
  
  # Traitement et Graphique densité des longitudes des pressions min pour les flux/MPD forts
  pos.dP <- which(dP>quantile(dP,probs=0.8))
  pos.dPgrad <- which(dP.grad>quantile(dP.grad,probs=0.8))
  
  # Min
  png(filename = paste0(get.dirstr(k,rean,period="present"),"compare.dP.grad/compare_lon_pressure_min_q80_k",k,"_mean",nbdays,"day_",rean,"_",start,"_",end,".png"),width = 7,height = 5,units = "in",res = 600)
  plot(density(lon.min[pos.dP]),ylim=c(0,0.1),xlab="Longitude (°)",main="Longitude du minimum de pression - MPD>q80 et MPD_grad>q80",lwd=2)
  grid()
  lines(density(lon.min[pos.dPgrad]),col="red",lwd=2)
  legend("topright",inset=.02,bty="n",col=c("black","red"),c("MPD","MPD_grad"),lty=1,lwd=2)
  graphics.off()
  
  # Max
  png(filename = paste0(get.dirstr(k,rean,period="present"),"compare.dP.grad/compare_lon_pressure_max_q80_k",k,"_mean",nbdays,"day_",rean,"_",start,"_",end,".png"),width = 7,height = 5,units = "in",res = 600)
  plot(density(lon.max[pos.dP]),xlab="Longitude (°)",main="Longitude du maximum de pression - MPD>q80 et MPD_grad>q80",lwd=2)
  grid()
  lines(density(lon.max[pos.dPgrad]),col="red",lwd=2)
  legend("topright",inset=.02,bty="n",col=c("black","red"),c("MPD","MPD_grad"),lty=1,lwd=2)
  graphics.off()
  
  # Traitement et Graphique flux et MPD des precip extremes Atlantiques en fonction de la longitude
  pos.extr <- get.ind.max.flow(flow = 1,agreg = T,nbdays = nbdays,start = start,end = end,spazm = T,supseuil = T,nei = T)
  
  # MPD
  png(filename = paste0(get.dirstr(k,rean,period="present"),"compare.dP.grad/dP_lon_pressure_min_precip_extr_k",k,"_mean",nbdays,"day_",rean,"_",start,"_",end,".png"),width = 7,height = 5,units = "in",res = 600)
  plot(lon.min,dP,pch=19,cex=0.5,xlab="Longitude (°)",ylab="MPD (m)",main=paste0("Precip Max Atlantique - ",nbdays," jours"))
  grid()
  points(lon.min,dP,pch=19,cex=0.5)
  abline(h=quantile(dP,probs=0.8),col="red",lwd=2)
  text(22.7,quantile(dP,probs=0.9),"q80",col="red",font=2,cex=0.7)
  points(lon.min[pos.extr],dP[pos.extr],pch=19,col="blue")
  graphics.off()
  
  # MPD_grad
  png(filename = paste0(get.dirstr(k,rean,period="present"),"compare.dP.grad/dPgrad_lon_pressure_min_precip_extr_k",k,"_mean",nbdays,"day_",rean,"_",start,"_",end,".png"),width = 7,height = 5,units = "in",res = 600)
  plot(lon.min,dP.grad,pch=19,cex=0.5,xlab="Longitude (°)",ylab="MPD_grad (m/°)",main=paste0("Precip Max Atlantique - ",nbdays," jours"))
  grid()
  points(lon.min,dP.grad,pch=19,cex=0.5)
  abline(h=quantile(dP.grad,probs=0.8),col="red",lwd=2)
  text(22.7,quantile(dP.grad,probs=0.9),"q80",col="red",font=2,cex=0.7)
  points(lon.min[pos.extr],dP.grad[pos.extr],pch=19,col="blue")
  graphics.off()
}

# Calcul de la latitude du max de vent zonal
compute.lat.maxwind <- function(var="uwind",rean="ERA5",start="1950-01-01",end="2019-12-31"){
  
  # Import
  data <- getdata(k = 1,day0 = start,day1 = end,rean = rean,climat = NULL,return.lonlat = T,var = var)
  lon <- data$lon;lat <- data$lat;data <- data$data
  
  # Traitement
  zonal.mean <- apply(data,3,colMeans)
  lat.max <- apply(zonal.mean,2,function(v){lat[which.max(v)]})
  save(lat.max,file = paste0(get.dirstr(k = 1,rean = rean,period = "past"),"compute.lat.maxwind/lat_max_zonal_wind",start,"_",end,".Rdata"))
}

# Renvoie un tableau de l'occurence des WP agreges par saison et annee
get.occ.wp.agreg <- function(start="1950-01-01",end="2019-12-31"){
  
  # Import
  wp <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Traitement
  season <- c("winter","spring","summer","autumn")
  tab <- matrix(data = NA,nrow = 5,ncol = 4)
  tab[1,] <- round(table(wp)/length(wp)*100,0) # year
  
  for(i in 1:length(season)){
    ind <- get.ind.season.past(sais = season[i],start = start,end = end,nbdays = 1)
    tab[i+1,] <- round(table(wp[ind])/length(ind)*100,0)
  }
  
  # Export
  write.csv(tab,file="2_Travail/0_Present/Rresults/get.occ.wp.agreg/tab_occ_wp_agreg.csv",row.names = F)
}

# Trace la carte de l'Europe avec frontières Auer et Grenoble
image.europe.auer <- function(){
  
  # Import
  ligne=read.csv("2_Travail/Data/Carto/Shape_CRSM.csv")
  
  # Mise en forme
  ligne.indiv <- unique(ligne[,1])
  pos.rect <- match(ligne.indiv,ligne[,1])
  rectang <- ligne[pos.rect,2:3]
  rectang <- apply(rectang,2,range)
  ligne <- ligne[-pos.rect,]
  
  # Carte
  png("2_Travail/0_Present/Rresults/image.europe.auer/map_europe_auer.png",width = 6,height = 4,units = "in",res = 600)
  par(mar=c(0,0,0,0))
  map(database= "world", ylim=c(39,52), xlim=c(-5,21), col=alpha("darkgreen",alpha = 0.3),fill=T)
  rect(xleft = rectang[1,1],ybottom = rectang[1,2],xright = rectang[2,1],ytop = rectang[2,2],lwd=2,border="black",lty=2)
  for(i in 1:length(ligne.indiv)){
    lines(ligne[ligne[,1]==ligne.indiv[i],2],ligne[ligne[,1]==ligne.indiv[i],3],lwd=3,col="purple")
  }
  points(x = 5.73,y = 45.18,pch=19,cex=1.5,col="red")
  box()
  graphics.off()
}

# Trace la carte de la France avec Grenoble
image.france <- function(){
  
  # Import
  FranceFormes <- getData(name="GADM", country="FRA", level=0)
  
  # Carte
  png("2_Travail/0_Present/Rresults/image.france/map_france_grenoble.png",width = 5,height = 5,units = "in",res = 600)
  par(mar=c(0,0,0,0))
  plot(FranceFormes)
  points(x = 5.73,y = 45.18,pch=19,cex=1.5,col="red"); text(x = 5.73-0.4,y = 45.18-0.4,"Grenoble",col="red",font=4) # Grenoble
  points(x = 2.35,y = 48.86,pch=19); text(x = 2.35+0.4,y = 48.86+0.4,"Paris",font=3) # Paris
  points(x = 4.83,y = 45.76,pch=19); text(x = 4.83+0.4,y = 45.76+0.4,"Lyon",font=3) # Lyon
  points(x = 5.38,y = 43.30,pch=19); text(x = 5.38+0.4,y = 43.30+0.4,"Marseille",font=3) # Marseille
  graphics.off()
}

# Champs composite d'une variable sur deux periodes et difference, pour un wp et une saison
map.abs.diff.var.wp <- function(var="t850",wp=1,sais="winter",rean="ERA5",start="1950-01-01",end="2017-12-31",extr=F,bv=NULL,extr.wp=T,wet=F){
  
  k <- 1
  dates <- getdates(start,end)
  
  # Import des donnees
  print("Import")
  data <- getdata(k = k,day0 = start,day1 = end,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = var,return.lonlat = T)
  lon <- data$lon;lat <- data$lat;data <- data$data
  if(substr(var,1,3)=="sph") data <- data*1000 # passage en g/kg
  if(var=="slp") data <- data/100 # passage en hPa
  gc()
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4))
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  per <- rep(1,length(dates))
  per[as.Date(dates)>=sep] <- 2
  
  # WP et saison
  tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  ind <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = 1)
  
  # Traitement
  print("Traitement")
  if(!extr & !wet){
    pos.per1 <- intersect(ind,which(per==1 & tt==wp))
    pos.per2 <- intersect(ind,which(per==2 & tt==wp))
  }
  if(extr){
    ind.max <- get.ind.max.sais(sais = sais,wp = ifelse(extr.wp,wp,"all"),nbdays = 1,start = start,end = end,bv = bv,spazm = T)
    pos.per1 <- intersect(ind.max,which(per==1))
    pos.per2 <- intersect(ind.max,which(per==2))
  }
  if(wet){
    ind.wet <- get.ind.precip.sais(sais = sais,wp = wp,nbdays = 1,start = start,end = end,bv = bv,spazm = T,rain = F)
    pos.per1 <- intersect(ind.wet,which(per==1))
    pos.per2 <- intersect(ind.wet,which(per==2))
  }
  
  data.per1 <- data[,,pos.per1]
  data.per2 <- data[,,pos.per2]
  mean.per1 <- apply(data.per1,1:2,mean)
  mean.per2 <- apply(data.per2,1:2,mean)
  diff.per <- mean.per2 - mean.per1
  
  # Cartes
  print("Cartes")
  data(wrld_simpl)
  
  png(filename = paste0(get.dirstr(k,rean,"past"),"map.abs.diff.var.wp/map_abs_diff_",var,ifelse(k==2,"1000",""),"_",sais,"_wp",wp,ifelse(extr,paste0("_extr_",bv,ifelse(!extr.wp,"_wpall","")),""),ifelse(wet,"_wet",""),".png"),width = 18,height = 5.5,units = "in",res=600)
  par(mfrow=c(1,3),mar=c(0,0,3,10),pty="s")
  
  # Periode 1
  param <- get.param.map(ran = range(mean.per1,mean.per2),var = var,type = "Mean")
  image.plot(lon,lat,mean.per1,xlim=c(-15,25),ylim=c(25,65),main=period1,
        col=param$col,xaxt="n",yaxt="n",xlab="",ylab="",breaks = param$breaks,
        legend.lab = param$leg,legend.line =-3,legend.mar = 10,legend.cex=1.5,cex.main=3,legend.width=2,axis.args=list(cex.axis=2))
  plot(wrld_simpl, add = TRUE)
  box()
  
  # Periode 2
  image.plot(lon,lat,mean.per2,xlim=c(-15,25),ylim=c(25,65),main=period2,
             col=param$col,xaxt="n",yaxt="n",xlab="",ylab="",breaks = param$breaks,
             legend.lab = param$leg,legend.line =-3,legend.mar = 10,legend.cex=1.5,cex.main=3,legend.width=2,axis.args=list(cex.axis=2))
  plot(wrld_simpl, add = TRUE)
  box()
  
  # Difference
  param <- get.param.map(ran = range(diff.per),var = var,type = "Diff")
  image.plot(lon,lat,diff.per,xlim=c(-15,25),ylim=c(25,65),main=paste0(period2," minus ",period1),
             col=param$col,xaxt="n",yaxt="n",xlab="",ylab="",breaks = param$breaks,
             legend.lab = param$leg,legend.line =-3,legend.mar = 10,legend.cex=1.5,cex.main=3,legend.width=2,axis.args=list(cex.axis=2))
  plot(wrld_simpl, add = TRUE)
  box()
  
  graphics.off()
}

# Carte de geopotentiel 500hPa pour une journée avec min et max pour explication MPD
map.geo.dP <- function(date="1963-03-09",rean="ERA5"){
  
  # Import
  data <- getdata(k = 1,day0 = date,day1 = date,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = F,ssp = NULL,return.lonlat = T,var = "hgt")
  
  # Traitement: pt min et max
  pos.min <- which(data$data == min(data$data), arr.ind = TRUE)
  pos.max <- which(data$data == max(data$data), arr.ind = TRUE)
  print(min(data$data))
  print(max(data$data))
  
  # Carte
  png(filename = paste0("2_Travail/0_Present/Rresults/map.geo.dP/map_k1_dP_",date,"_",rean,".png"),width = 6,height = 5,units = "in",res = 600)
  par(pty="s",mar=c(4.5,4,3,1))
  map.geo(date = date,rean = rean,climat = NULL,run = 1,k = 1,nbdays = 1,save = F,win = T,let = F,leg = T,iso = T,wind = F,condens = F,ssp = NULL)
  points(data$lon[pos.min[1]],data$lat[pos.min[2]],pch=21,bg="blue",cex=1.5)
  points(data$lon[pos.max[1]],data$lat[pos.max[2]],pch=21,bg="red",cex=1.5)
  graphics.off()
}

# Carte de geopotentiel 500hPa pour une journée et plusieurs reanalyses
map.geo.rean <- function(date="1993-10-08",rean=c("20CR","ERA20C","ERA5")){
  
  png(filename = paste0("2_Travail/0_Present/Rresults/map.geo.rean/map_k1_",date,"_",paste(rean,collapse="_"),".png"),width = 8,height = 3,units = "in",res = 600)
  layout(mat = matrix(1:(length(rean)+1),nrow = 1,ncol = length(rean)+1),widths = c(rep(1,length(rean)),0.25),heights = rep(1,length(rean)+1))
  par(mar=c(0,0,3,1))
  
  # Cartes
  for(i in 1:length(rean)){
    map.geo(date = date,rean = rean[i],climat = NULL,run = 1,k = 1,nbdays = 1,save = F,win = F,let = F,leg = F,iso = T,wind = F,condens = T,ssp = NULL)
    title(rean[i],cex.main=2)
  }
  # Legende
  N <- 11
  leg <- seq(4900,6100,200)
  ran <- seq(4850,6100)
  
  par(pty="m",mar=c(1,0,4,0))
  plot(1,1,xaxt="n",yaxt="n",bty="n",type="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = paste0("        ",leg),at = ecdf(ran)(leg),
              vertical = T,xlim = c(0.1,0.5),ylim = c(0,1),cex=1.4)
  graphics.off()
}

# Carte de geopotentiel 500hPa pour plusieurs journées avec score TWS
map.geo.tws <- function(dat=c("1961-07-02","1986-08-10","2013-01-26"),nbdays=2,rean="ERA5"){
  
  dates <- getdates(start = "1950-01-01",end = "2017-12-31")
  
  # Import TWS jours consecutifs (celerite)
  cel <- get.descriptor(descriptor = "cel",k = 1,dist = "TWS",nbdays = 1,start = "1950-01-01",end = "2017-12-31",standardize = F,rean = rean,
                        threeday = F,desais = F,period = "present")
  
  # Cartes
  nb.dates <- length(dat)
  
  png(filename = paste0("2_Travail/0_Present/Rresults/map.geo.tws/map_k1_",paste(dat,collapse="_"),"_",nbdays,"day_",rean,".png"),width = 6,height = 7,units = "in",res = 600)
  layout(mat = matrix(data = 1:((nbdays+2)*nb.dates),nrow = nb.dates,ncol = nbdays+2,byrow = T),widths = c(rep(1,nbdays),0.3,0.5),heights = rep(1,nb.dates))
  
  for(i in 1:nb.dates){
    
    par(mar=c(1,0.1,0,0.5))
    
    # Cartes
    for(j in 1:nbdays){
      map.geo(date = as.character(as.Date(dat[i])+j-1),rean = rean,climat = NULL,run = 1,k = 1,nbdays = 1,save = F,win = T,let = F,leg = F,iso = T,wind = F,condens = T,ssp = NULL)
    }
    
    # Legende
    N <- 11
    leg <- seq(4900,6100,200)
    ran <- seq(4850,6100)
    
    par(pty="m",mar=c(2,0,1,0))
    plot(1,1,xaxt="n",yaxt="n",bty="n",type="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
    colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
                labels = paste0("        ",leg),at = ecdf(ran)(leg),
                vertical = T,xlim = c(0.1,0.5),ylim = c(0,1),cex=1.4)
    
    # TWS et percentile de celerite
    per <- round(ecdf(cel)(cel[which(dates==as.character(as.Date(dat[i])+nbdays-1))]),2)
    plot(1,1,xaxt="n",yaxt="n",bty="n",type="n",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1))
    text(0.5,0.5,paste0("TWS = ",round(cel[which(dates==as.character(as.Date(dat[i])+nbdays-1))],2),"\n (q=",per,")"),font=2,cex=1.7)
    print(per)
  }
  graphics.off()
}

# Percentile des indicateurs aux dates de crues de l'Isère ou du Drac
plot.descr.flood <- function(sais="winter",nbdays=1,rean="20CR-m1",start = "1851-01-01",end = "2010-12-31"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  colnames(des) <- nam2str(descr)
  
  for(i in 1:length(descr)){
    des[,i] <- get.descriptor(descriptor = descr[i],k = 1,dist = "TWS",nbdays = nbdays,start = start,end = end,
                              standardize = F,rean = rean,threeday = F,desais = F,period = "past",start.ana = start,
                              end.ana = end)
    des[,i] <- ecdf( des[,i])( des[,i])*100
  }
  
  # Import des dates de crues
  bd <- read_bd(path = "2_Travail/Data/Event/BD-RTM-IGE_20201208_finale.csv")
  bd <- bd[bd$type=="R"|bd$type=="TR",] # que riviere ou torrent+riviere
  rownames(bd) <- 1:nrow(bd)
  bd <- bd[-c(27,41),] # car romanche seule
  flood <- bd$date
  
  if(nbdays>1){
    flood <- flood-nbdays+1
    pos <- which(!is.na(bd$fin_hydro))
    flood[pos] <- bd$fin_hydro[pos]-nbdays+1 # on prend la sequence s'arretant a fin_hydro
  }
  
  flood <- as.character(flood)
  flood.sais <- flood[which(flood %in% dates[get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)])]
  
  # Graphique
  season <- c("winter","spring","summer","autumn")
  colo <- c("cornflowerblue","darkgreen","red","burlywood1")
  main <- paste0(nam2str(sais)," (",length(flood.sais)," floods)")
  
  des <- des[match(flood.sais,dates),]
  boxplot(des,ylim=c(0,100),col=colo[sais==season],main=main,ylab="Percentile of descriptor value (%)")
  grid();par(new=T)
  boxplot(des,ylim=c(0,100),col=colo[sais==season],main=main,ylab="Percentile of descriptor value (%)")
  
  # Ajout crues
  if(sais=="autumn"){
    pos <- which(substr(flood.sais,1,6)=="1859-1")
    text(1:4,des[pos,],"2",font=2,col="red")
    pos <- which(substr(flood.sais,1,6)=="1928-1")
    text(1:4,des[pos,],"3",font=2,col="red")
    #pos <- which(substr(flood.sais,1,7)=="1968-09")
    #text(1:4,des[pos,],"3",font=2,col="red")
  }
  if(sais=="spring"){
    pos <- which(substr(flood.sais,1,7)=="1856-05")
    text(1:4,des[pos,],"1",font=2,col="red")
    pos <- which(substr(flood.sais,1,7)=="2001-03")
    text(1:4,des[pos,],"4",font=2,col="darkblue")
    pos <- which(substr(flood.sais,1,7)=="2008-05")
    text(1:4,des[pos,],"5",font=2,col="darkblue")
    pos <- which(substr(flood.sais,1,7)=="2010-05")
    text(1:4,des[pos,],"6",font=2,col="darkblue")
    
  }
  #if(sais=="summer"){
  #  pos <- which(substr(flood.sais,1,9)=="1914-07-2")
  #  text(1:4,des[pos,],"8",font=2)
  #  pos <- which(substr(flood.sais,1,7)=="1955-06")
  #  text(1:4,des[pos,],"9",font=2)
  #}
}

# Trace differentes caracteristiques des crues de riviere de la BD RTM-IGE
plot.event.river <- function(){
  
  # Montrer que la saison
  # Separer en Isere seul, Drac seul, Isere+Drac: colorier la barre de chaque mois!
  
  # Import
  bd <- read_bd(path = "2_Travail/Data/Event/BD-RTM-IGE_20201208_finale.csv")
  
  # Traitement
  bd <- bd[bd$type=="R"|bd$type=="TR",] # que riviere ou torrent+riviere
  bd <- bd[-c(27,41),] # car romanche seule
  rownames(bd) <- 1:nrow(bd)
  #bd <- bd[1:39,]
  mon <- as.numeric(substr(bd$date,6,7))
  
  # Isere seule, Drac seul, ou les deux
  bv <- rep(NA,length(mon))
  bv[!is.na(bd$Categorie_Coeur_Is) & !is.na(bd$Categorie_Coeur_Dr)] <- "Isere+Drac"
  bv[!is.na(bd$Categorie_Coeur_Is) & is.na(bd$Categorie_Coeur_Dr)] <- "Isere"
  bv[is.na(bd$Categorie_Coeur_Is) & !is.na(bd$Categorie_Coeur_Dr)] <- "Drac"
  bv[is.na(bv)] <- bd$autre_site_cite[is.na(bv)]
  bv[bv=="Isere,Drac"] <- "Isere+Drac"
  bv[c(10,18,32,39)] <- "Isere+Drac"
  bv[c(27,28,40)] <- "Isere"
  
  # Mise en forme
  tmp <- data.frame(Month=mon,BV=bv,Count=rep(1,length(mon)))
  res <- aggregate(tmp$Count,by=list(tmp$Month,tmp$BV),sum)
  colnames(res) <- colnames(tmp)
  res$BV <- factor(res$BV,levels = c("Isere+Drac","Drac","Isere"))
  res$Month <- factor(res$Month)
  
  # Graphique saison et BV
  colo <- c("grey","burlywood1","cornflowerblue")
  print(paste0("Nombre de crues: ",length(mon)))
  print(paste0("Nombre de crues Isere+Drac: ",sum(bv=="Isere+Drac")))
  print(paste0("Nombre de crues Isere: ",sum(bv=="Isere")))
  print(paste0("Nombre de crues Drac: ",sum(bv=="Drac")))
  
  png(filename = "2_Travail/0_Present/Rresults/plot.event.river/plot_event_sais_bv.png",width = 10,height = 6,units = "in",res = 600)
  par(mar=c(3,4,0.5,0.5))
  ggplot(res, aes(fill=BV, y=Count, x=Month)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 15,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 15,face = "bold"),
          axis.text.x = element_text(size=12),axis.text.y = element_text(size=12),
          plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=18),
          legend.key.si = unit(1,"cm"),legend.text = element_text(size=15),
          legend.title = element_text(hjust=0.5,vjust=1,size = 18,face = "bold"),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    geom_bar(position="stack",stat = "identity",width = 1,colour="black",size=1)+
    scale_fill_manual(values=colo)+
    scale_x_discrete(breaks=1:12, labels=month.abb)+
    grids(axis = "xy",color = "grey",linetype = "dashed")
  graphics.off()
  
  #hist(mon,breaks=0:12,xaxt="n",col="cornflowerblue",border = "royalblue",xlab = "Month",main="") # attention à mai sur-représenté avec les décennales
  #grid(ny=NULL,nx=NA)
  #hist(mon,breaks=0:12,xaxt="n",col="cornflowerblue",border = "royalblue",xlab = "Month",main="",add=T) # attention à mai sur-représenté avec les décennales
  #axis(side = 1,at = 0.5:11.5,month.abb)
  #
  ## Graphique intensite
  #int.is <- aggregate(bd$Categorie_Coeur_Is,by=list(mon),mean,na.rm=T)[,2]
  #int.dr <- aggregate(bd$Categorie_Coeur_Dr,by=list(mon),mean,na.rm=T)[,2]
  #barplot(int.is,xaxt="n",col="cornflowerblue",border = "royalblue",xlab = "Month",main="")
  #barplot(int.dr,xaxt="n",col="cornflowerblue",border = "royalblue",xlab = "Month",main="")
  #
  #q.is <- aggregate(bd$Q_Is_bast,by=list(mon),mean,na.rm=T)[,2]
  #q.dr <- aggregate(bd$Q_Drac,by=list(mon),mean,na.rm=T)[,2]
  #barplot(q.is,xaxt="n",col="cornflowerblue",border = "royalblue",xlab = "Month",main="")
  #barplot(q.dr,xaxt="n",col="cornflowerblue",border = "royalblue",xlab = "Month",main="")
  
  # Graphique evenements cumules
  ann <- seq(1850,2019)
  
  bd.hist <- bd[bd$Coeur_Lang==1,]
  ann.event <- as.numeric(substr(bd.hist$date,1,4))
  cum.event <- rep(0,length(ann))
  cum.event[match(ann.event,ann)] <- 1
  cum.event <- cumsum(cum.event)
  
  bd.rtm <- bd[bd$RTM_I==1,]
  ann.event <- as.numeric(substr(bd.rtm$date,1,4))
  cum.event <- rep(0,length(ann))
  cum.event[match(ann.event,ann)] <- 1
  cum.event <- cumsum(cum.event)
  
  plot(ann,cum.event,type="l")
  lines(ann,cum.event,col="blue")
}

# Trace differents graphiques de latitude du max du vent zonal, par wp
plot.lat.maxwind <- function(wp=1,rean="ERA5",start="1950-01-01",end="2019-12-31"){
  
  dates <- getdates(start,end)
  ann <- unique(substr(dates,1,4))
  
  # Imports
  load(paste0(get.dirstr(k = 1,rean = rean,period = "past"),"compute.lat.maxwind/lat_max_zonal_wind",start,"_",end,".Rdata"))
  tt <- get.wp(nbdays = 1,start = start,end = end,agreg = T,spazm = T)
  pos.tt <- which(tt==wp)
  if(wp=="all"){pos.tt <- 1:length(dates)}
  
  # Graphique densité
  season <- c("winter","spring","summer","autumn")
  png(filename = paste0(get.dirstr(k = 1,rean = rean,period = "past"),"plot.lat.maxwind/plot_density_max_zonal_wind_wp",wp,"_",start,"_",end,".png"),width = 9,height = 7,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(4,4,2,0))
  for(i in 1:length(season)){
    pos.sea <- get.ind.season.past(sais = season[i],start = start,end = end,nbdays = 1)
    pos.i <- intersect(pos.tt,pos.sea)
    plot(density(lat.max[pos.i]),xlab="Latitude (°)",main=nam2str(season[i]),cex.main=2)
    grid();par(new=T)
    plot(density(lat.max[pos.i]),xlab="Latitude (°)",main=nam2str(season[i]),cex.main=2)
  }
  graphics.off()
  
  # Graphique tendance
  if(wp!="all"){lat.max[tt!=wp] <- NA}
  png(filename = paste0(get.dirstr(k = 1,rean = rean,period = "past"),"plot.lat.maxwind/plot_trend_max_zonal_wind_wp",wp,"_",start,"_",end,".png"),width = 9,height = 7,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(2,4,2,4))
  
  for(i in 1:length(season)){
    # Calcul
  pos.i <- get.ind.season(sais = season[i],start = start,end = end)
  tmp <- lat.max[pos.i$pos.season]
  tmp[pos.i$pos.NA] <- NA
  dim(tmp) <- c(pos.i$l.season,pos.i$n.season)
  vec.i <- apply(tmp,2,mean,na.rm=T)
  if(season[i]=="winter"){vec.i <- c(NA,vec.i)}
  
  # Plot
  plot(ann,vec.i,type="l",xlab="",ylab="Latitude (°)",main=nam2str(season[i]),cex.main=2)
  grid();par(new=T)
  plot(ann,vec.i,type="l",xlab="",ylab="Latitude (°)",main=nam2str(season[i]),cex.main=2)
  abline(lm(vec.i~as.numeric(ann),na.action = "na.exclude"))
  
  # NAO
  shif <- round(mean(range(vec.i,na.rm=T),na.rm=T),0) # decalage a appliquer à NAO pour être dans le graphique
  nao <- get.nao(start = substr(start,1,4),end = substr(end,1,4),sais = season[i],daily = F)
  lines(nao[,1],nao[,2]+shif,col="darkblue")
  abline(lm((nao[,2]+shif)~as.numeric(nao[,1])),col="darkblue")
  axis(side = 4,at = pretty(range(vec.i,na.rm=T)),labels = pretty(range(vec.i,na.rm=T))-shif,col.axis="darkblue")
  mtext(text = "NAO",side = 4,srt=90,line = 2.5,col="darkblue")
  text(quantile(as.numeric(ann),0.1),min(vec.i,na.rm=T)+0.5,paste0("R=",round(cor(vec.i,nao[,2],use = "pairwise.complete.obs"),2)))
  }
  graphics.off()
}
  
# Trace des percentiles de precip associes aux max annuels de debit sur l'Isere a St Gervais
plot.quant.rain <- function(nbdays=3,spazm=F,start = "1969-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  
  # Imports
  precip <- get.precip(bv = "Isere",spazm=spazm,nbdays = nbdays,start = start,end = end)
  crues <- read.table("2_Travail/0_Present/Rresults/image.cumul/crues_1969_2017_Isere_StGervais.txt",header = F)
  crues <- as.data.frame(crues)
  
  # Traitement
  event <- as.character(as.Date(crues[,1],format="%d/%m/%Y")-nbdays+1) # on prend les pluies de j-3,j-2,j-1 de la crue (car 8hj a 7hj+1)
  event.ann <- as.numeric(substr(event,1,4))
  event.month <- as.numeric(substr(event,6,7))
  ind   <- match(event,dates)
  
  distrib <- ecdf(precip)(precip)*100
  dis.ind <- distrib[ind]
  
  debit <- crues[,2]
  
  # Graphique
  colo <- colorRampPalette(c("darkblue","lightgreen","purple"))(12)
  
  png(file = paste0("2_Travail/0_Present/Rresults/image.cumul/Quant_rain_flood_last_",nbdays,"days",ifelse(spazm,"_spazm",""),".png"),width = 6,height = 4,units = "in",res = 600)
  par(mar=c(2,4,0.5,0))
  plot(event.ann,dis.ind,pch=21,ylim=c(0,100),ylab=paste0(nbdays,"-day precipitation percentile (%)"),cex=debit/mean(debit)*2)
  grid()
  abline(h=90,col="red",lty=2)
  points(event.ann,dis.ind,pch=21,bg=colo[event.month],cex=debit/mean(debit)*2)
  colorlegend(colbar = colo,labels = 1:12,at = seq(0,1,length.out = 12),xlim = c(1992,2017),ylim = c(0,10),vertical = F)
  text(2004.5,15,"Month")
  legend("bottomleft",inset=.02,pch=21,bty="n",
         pt.cex=c(max(debit/mean(debit)*2),median(debit/mean(debit)*2),min(debit/mean(debit)*2)),
         paste0(c(max(debit,na.rm=T),median(debit,na.rm=T),min(debit,na.rm=T))," m3/s"))
  graphics.off()
  
  # Stats
  print(paste0("Nombre de max annuels: ",length(dis.ind)))
  print(paste0("Nombre superieur a q90: ",sum(dis.ind>90)))
  print(paste0("Nombre superieur a q95: ",sum(dis.ind>95)))
  print(paste0("Nombre superieur a q98: ",sum(dis.ind>98)))
  print(paste0("Nombre superieur a q99: ",sum(dis.ind>99)))
}

# Trace la saisonalité de thetawp avec max saisonniers de precip ou pluie par dessus
#plot.sais.thetawp <- function(bv="Isere-seul",rean="ERA5",start="1950-01-01",end="2017-12-31"){
#  
#  dates <- getdates(start,end)
#  dates.pos <- unique(substr(dates,6,10))
#  
#  # Import thetawp à Grenoble (Ecoutoux à ces coordonnees)
#  thetawp <- getdata(k = 1,day0 = start,day1 = end,rean = rean,pt.lon = 5.75,pt.lat = 45.25,var = "thetawp")
#  
#  # Traitement
#  pos.flood <- match(c("2001-03-22","2008-05-30","2010-05-31","2015-05-02"),dates)
#  
#  plot.sais.quantile(dates = dates,vec = thetawp,y.label = "Thetawp (°C)",main = "",colo = "red")
#  season <- c("winter","spring","summer","autumn")
#  for(i in 1:length(season)){
#    pos.i <- get.ind.max.sais(sais = season[i],wp = "all",nbdays = 1,start = start,end = end,bv = bv,spazm = T,rain = T)
#    x.i <- match(substr(dates[pos.i],6,10),dates.pos)
#    points(x.i,thetawp[pos.i],pch=19,col="red",cex=0.7)
#    if((pos.flood-1) %in% pos.i){
#      pos.ii <- which(pos.i %in% (pos.flood-1))
#      points(x.i[pos.ii],thetawp[pos.i[pos.ii]],pch=19,col="darkgreen")
#      }
#  }
#  
#  
#  
#}

# Graphique d'evolution des max saisonniers de precip avec droites GEV (sans enregistrement)
plot.trend.precip.gev.wp <- function(bv="Isere-seul",nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31",leg=T,rain=F){
  
  dates <- getdates(start,end)
  ann <- as.numeric(unique(substr(dates,1,4)))
  
  # Imports
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm,rain = rain)
  load(paste0("2_Travail/1_Past/Rresults/fit.gev.obs/fit_gev_",ifelse(rain,"rain_",""),bv,"_mean",nbdays,"day",ifelse(spazm,"_spazm",""),"_",start,"_",end,".Rdata"))
  wp <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = spazm)
  
  # Traitement et Graphique
  sais <- c("winter","spring","summer","autumn")
  meanGEV.fct<-function(mu,sig,xi){mu+sig/xi*(gamma(1-xi)-1)}
  
  for(i in 1:length(sais)){
    
    # Traitement max saisonnier
    pos <- get.ind.max.sais(sais = sais[i],wp = "all",nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm,rain = rain)
    max.vec <- precip[pos]
    if(sais[i]=="winter"){pos <- c(NA,pos);max.vec <- c(NA,max.vec)}
    
    # Graphique max saisonnier
    plot(ann,max.vec,type="l",xlab="",ylab="",ylim=c(0,100),main=nam2str(sais[i]))
    grid();par(new=T)
    plot(ann,max.vec,type="l",xlab="",ylab="",ylim=c(0,100),main=nam2str(sais[i]))
    
    # Points WP
    colo <- rep("black",length(pos))
    colo[wp[pos]==1] <- "cornflowerblue"
    colo[wp[pos]==2] <- "burlywood1"
    per <- round(table(wp[pos])/length(pos)*100,0)
    
    points(ann,max.vec,pch=21,bg=colo)
    legend("topleft",c(paste0("Atlantic (",per[1],"%)"),paste0("Mediterranean (",per[2],"%)")),pt.bg=c("cornflowerblue","burlywood1"),pch=21,bty="n")
    
    # Traitement GEV
    fit.i <- fit[[i]]
    rl20 <- mea <- NA
    
    for(j in 1:nrow(fit.i$fitbest$vals)){
      rl20[j] <- qgev(1-1/20,fit.i$fitbest$vals[j,1],fit.i$fitbest$vals[j,2],fit.i$fitbest$vals[j,3])
      mea[j] <- meanGEV.fct(fit.i$fitbest$vals[j,1],fit.i$fitbest$vals[j,2],fit.i$fitbest$vals[j,3])
    }
    if(sais[i]=="winter"){rl20 <- c(NA,rl20);mea <- c(NA,mea)}
    
    lty.i <- ifelse(fit.i$pval<0.05,1,2)
    
    # Graphique GEV
    lines(ann,mea,col="blue",lty=lty.i,lwd=2)
    lines(ann,rl20,col="red",lty=lty.i,lwd=2)
    if(i==1 & leg){legend("topright",c("RL20","Mean"),lty=c(1,1),lwd=2,col=c("red","blue"),bty="n")}
    #if(i==1 & leg){legend("topright",c("Mean"),lty=1,lwd=2,col="blue",bty="n")}
  }
}

# Graphique d'evolution de la temperature saisonniere et aux dates des max, par wp
plot.trend.temp <- function(bv="Isere-seul",nbdays=1,wp="all",start="1950-01-01",end="2017-12-31",leg=T,rain=F){
  
  dates <- getdates(start,end)
  ann <- as.numeric(unique(substr(dates,1,4)))
  
  # Imports
  temp <- get.temp(nbdays = nbdays,start = start,end = end,bv = bv)
  temp.qua <- ecdf(temp)(temp)
  tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  if(wp!="all"){temp[tt!=wp] <- NA}
  
  # Traitement et graphique
  season <- c("winter","spring","summer","autumn")
  
  for(i in 1:length(season)){
    
    # Moyennes saisonnieres
    pos.i <- get.ind.season(sais = season[i],start = start,end = end)
    tmp <- temp[pos.i$pos.season]
    tmp[pos.i$pos.NA] <- NA
    dim(tmp) <- c(pos.i$l.season,pos.i$n.season)
    vec.i <- apply(tmp,2,mean,na.rm=T)
    if(season[i]=="winter"){vec.i <- c(NA,vec.i)}
    
    # Dates des max de precip
    pos.max.i <- get.ind.max.sais(sais = season[i],wp = wp,nbdays = nbdays,start = start,end = end,bv = bv,spazm = T,rain = rain)
    if(season[i]=="winter"){pos.max.i <- c(NA,pos.max.i)}
    #print(temp.qua[pos.max.i])
    #print(temp[pos.max.i])
    
    # Plot
    plot(ann,vec.i,type="l",xlab="",ylab="Temperature (°C)",main=nam2str(season[i]),ylim=range(c(temp[pos.max.i],vec.i),na.rm = T))
    grid();par(new=T)
    plot(ann,vec.i,type="l",xlab="",ylab="Temperature (°C)",main=nam2str(season[i]),ylim=range(c(temp[pos.max.i],vec.i),na.rm = T))
    abline(lm(vec.i~as.numeric(ann)))
    
    lines(ann,temp[pos.max.i],col="red")
    abline(lm(temp[pos.max.i]~as.numeric(ann)),col="red")
    if(i==1 & leg){legend("bottomright",c("Mean",ifelse(rain,"Max Rain","Max Precip")),lty=c(1,1),lwd=2,col=c("black","red"),bty="n")}
  }
}

# Graphique d'evolution de la temperature saisonniere pour differents wp
plot.trend.temp.wp <- function(wp=c(1,2),bv="Isere-seul",nbdays=1,start="1950-01-01",end="2017-12-31",leg=T){
  
  dates <- getdates(start,end)
  ann <- as.numeric(unique(substr(dates,1,4)))
  
  # Imports
  temp <- get.temp(nbdays = nbdays,start = start,end = end,bv = bv)
  tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Traitement et graphique
  season <- c("winter","spring","summer","autumn")
  vec.sais <- vector(mode="list",length = length(season))
  wp.all <- c(1,2,5,8)
  wp.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  colo <- c("blue","red","darkgreen","dimgrey")
  
  for(i in 1:length(season)){ # 2 boucles j (1 traitement, 1 graphique) pour avoir le range des valeurs pour ylim du plot
    
    # Moyennes saisonnieres
    pos.i <- get.ind.season(sais = season[i],start = start,end = end)
    tmp <- temp[pos.i$pos.season]
    tmp[pos.i$pos.NA] <- NA
    dim(tmp) <- c(pos.i$l.season,pos.i$n.season)
    vec.i <- apply(tmp,2,mean,na.rm=T)
    if(season[i]=="winter"){vec.i <- c(NA,vec.i)}
    
    vec.sais[[i]][[1]] <- vec.i
    
    # Moyennes saisonnieres par WP
    for(j in 1:length(wp)){
      temp.j <- temp
      temp.j[tt!=wp[j]] <- NA
      
      tmp <- temp.j[pos.i$pos.season]
      tmp[pos.i$pos.NA] <- NA
      dim(tmp) <- c(pos.i$l.season,pos.i$n.season)
      vec.i <- apply(tmp,2,mean,na.rm=T)
      if(season[i]=="winter"){vec.i <- c(NA,vec.i)}
      vec.sais[[i]][[j+1]] <- vec.i
    }
    
    # Graphique
    plot(ann,vec.sais[[i]][[1]],type="l",xlab="",ylab="Temperature (°C)",main=nam2str(season[i]),ylim=range(unlist(vec.sais[[i]]),na.rm = T))
    grid();par(new=T)
    plot(ann,vec.sais[[i]][[1]],type="l",xlab="",ylab="Temperature (°C)",main=nam2str(season[i]),ylim=range(unlist(vec.sais[[i]]),na.rm = T))
    abline(lm(vec.sais[[i]][[1]]~as.numeric(ann)),lwd=2)
    
    # Boucle autres WP
    for(j in 1:length(wp)){
      lines(ann,vec.sais[[i]][[j+1]],col=colo[wp.all==wp[j]])
      abline(lm(vec.sais[[i]][[j+1]]~as.numeric(ann)),col=colo[wp.all==wp[j]],lwd=2)
    }
    if(i==1 & leg){legend("bottom",c("All",wp.name[match(wp,wp.all)]),lty=1,lwd=2,col=c("black",colo[match(wp,wp.all)]),bty="n")}
  }
}