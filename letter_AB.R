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

# Calcul correlation entre descripteurs (et NAO et WP) et precip au pas de temps journalier et annuel, par saison (de 1j à 90j).
# Les lissages journaliers ne se chevauchent pas entre differentes saisons. Une saison d'hiver est a cheval sur deux annees (on perd un hiver).
compute.cor.descr.precip.clean <- function(k,dist,bv="Isere-seul",sais="winter",spazm=T,start="1950-01-01",end="2017-12-31",rean){
  
  dates <- getdates(start,end)
  nyear <- length(unique(substr(dates,1,4)))
  descr <- c("celnei","singnei","rsingnei","dP")
  
  # Import indicateur
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  colnames(des) <- nam2str(descr,whole=T)
  for(i in 1:length(descr)){
  des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = 1,start = start,end = end,
                              standardize = F,rean = rean,threeday = F,period = "present")
  }
  
  # Import NAO
  NAO <- get.nao(start = substr(start,1,4),end = substr(end,1,4),sais = "all",daily = T)[,4]
  des <- cbind(des,NAO)
  
  # Import WP
  wp <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = spazm)
  wp <- cbind(wp,wp,wp,wp)
  colnames(wp) <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  wp[wp[,1]!=1,1] <- 0;wp[wp[,1]==1,1] <- 1
  wp[wp[,2]!=2,2] <- 0;wp[wp[,2]==2,2] <- 1
  wp[wp[,3]!=5,3] <- 0;wp[wp[,3]==5,3] <- 1
  wp[wp[,4]!=8,4] <- 0;wp[wp[,4]==8,4] <- 1
  des <- cbind(des,wp)
  
  # Import precip
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = spazm)
  
  # Dates de la saison d'interet
  if(sais=="winter") {sta <- "12-01"; N <- 90}
  if(sais=="spring") {sta <- "03-01"; N <- 92}
  if(sais=="summer") {sta <- "06-01"; N <- 92}
  if(sais=="autumn") {sta <- "09-01"; N <- 91}
  
  pos.sta <- which(substr(dates,6,10)== sta)
  if(sais=="winter"){pos.sta <- pos.sta[-length(pos.sta)]; nyear <- nyear-1} # on retire le dernier hiver incomplet

  pos.all <- pos.sta
  for(i in 1:(N-1)){ # saisons de longueur 90 jours
    pos.all <- sort(c(pos.all,pos.sta+i))
  }
  
  # Tableaux propres pour calcul correlation
  precip.sais <- precip[pos.all]
  dim(precip.sais) <- c(N,nyear)
  
  corr.daily <- matrix(data = NA,nrow = N,ncol = ncol(des))
  colnames(corr.daily) <- colnames(des)
  
  corr.yearly <- matrix(data = NA,nrow = 10,ncol = ncol(des))
  colnames(corr.yearly) <- colnames(des)
  
  for(i in 1:ncol(des)){
    print(i)
    # Indicateur
    des.sais <- des[pos.all,i]
    dim(des.sais) <- c(N,nyear)
    
    # Calcul correlation journaliere
    for(j in 1:N){
      des.tmp <- des.sais
      precip.tmp <- precip.sais
      
      des.j <- as.vector(apply(des.tmp,2,function(v){rollapply(v,j,mean)}))
      precip.j <- as.vector(apply(precip.tmp,2,function(v){rollapply(v,j,mean)}))
      corr.daily[j,i] <- cor(des.j,precip.j)
    }
    
    # Calcul correlation annuelle
    des.tmp <- apply(des.sais,2,mean)
    precip.tmp <- apply(precip.sais,2,mean)
    
    for(j in 1:10){
      des.j <- rollapply(des.tmp,j,mean)
      precip.j <- rollapply(precip.tmp,j,mean)
      corr.yearly[j,i] <- cor(des.j,precip.j)
    }
  }
  
  corr <- list(corr.daily=corr.daily,corr.yearly=corr.yearly)
  save(corr,file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.cor.descr.precip.clean/",bv,"_",sais,".Rdata"))
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

# Regroupement des plot de correlation en fonction du lissage pour les deux BV
plot.cor.all <- function(bv1="Isere-seul",bv2="Drac-seul",rean){
  
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/plot.cor.descr.precip.clean/plot_combine_",bv1,"_",bv2,".png"),width = 8,height = 9,units = "in",res=600)
  layout(matrix(c(rep(1,2),rep(2,2),rep(3,2),rep(4,2),5:8,rep(9,2),rep(10,2),11:14,rep(15,2),rep(16,2),17:20,rep(21,2),rep(22,2),23:30,rep(31,4)),nrow = 11,ncol = 4,byrow = T),widths = c(0.34,0.66,0.34,0.66),heights = c(0.3,0.13,1,0.13,1,0.13,1,0.13,1,0.15,0.4))
  
  # Noms des BV
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.5,0.5,nam2str(bv1),cex=3,bg=2)
  plot.new()
  text(0.5,0.5,nam2str(bv2),cex=3,bg=2)
  
  # Graphiques
  sais <- c("winter","spring","summer","autumn")
  bv <- c(bv1,bv2)
  for(i in 1:length(sais)){
    # Nom saison x 2
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,sais[i],cex=2,bg=2)
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,sais[i],cex=2,bg=2)
    
    for(j in 1:length(bv)){
      plot.cor.descr.precip.clean(k = 1,bv = bv[j],sais = sais[i],rean = rean,save = F)
    }
  }
  
  # xlab
  for(i in 1:2){
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.7,0.5,"Smoothing length\n(days)")
  
  par(mar=c(0,0,0,0))
  plot.new()
  text(0.5,0.5,"Smoothing length\n(years)")
  }
  
  # Legende
  colo <- c("black",
            "olivedrab3",
            "darkorange",
            "royalblue",
            "grey",
            "darkgoldenrod3",
            "darkorchid2",
            "red",
            "aquamarine2")
  
  ind.name <- c("Celerity",
                "Singularity",
                "Relative singularity",
                "MPD",
                "NAO",
                "Oceanic",
                "Mediterranean",
                "Northeast",
                "Anticyclonic")
  
  par(mar=c(0,3,0,0))
  plot.new()
  legend("center",legend = ind.name,col = colo,bty="n",lty=1,lwd=2,ncol=5,bg=2,cex=1.2)
  
  graphics.off()
}

# Graphique des correlations interannuelles par saison pour les deux BVs
plot.cor.interannual <- function(bv1="Isere-seul",bv2="Drac-seul",rean){
  
  # Import des correlations
  sais <- c("winter","spring","summer","autumn")
  bv <- c(bv1,bv2)
  
  # BV1
   corr.bv1 <- matrix(NA,4,9)
    for(i in 1:length(sais)){
      load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/compute.cor.descr.precip.clean/",bv1,"_",sais[i],".Rdata"))
      corr.bv1[i,] <- corr$corr.yearly[1,]
    }
   colnames(corr.bv1) <- colnames(corr$corr.yearly)
   sign.bv1 <- apply(corr.bv1,2,sign)
   pch.bv1 <- sign.bv1; pch.bv1 <- apply(pch.bv1,2,function(v){v[v==-1] <- 19;v[v==1] <- 15;return(v)})
   
   # BV2
   corr.bv2 <- matrix(NA,4,9)
   for(i in 1:length(sais)){
     load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/compute.cor.descr.precip.clean/",bv2,"_",sais[i],".Rdata"))
     corr.bv2[i,] <- corr$corr.yearly[1,]
   }
   colnames(corr.bv2) <- colnames(corr$corr.yearly)
   sign.bv2 <- apply(corr.bv2,2,sign)
   pch.bv2 <- sign.bv2; pch.bv2 <- apply(pch.bv2,2,function(v){v[v==-1] <- 19;v[v==1] <- 15;return(v)})
   
   # Graphique
   colo <- c("black",
             "olivedrab3",
             "darkorange",
             "royalblue",
             "grey",
             "darkgoldenrod3",
             "darkorchid2",
             "red",
             "aquamarine2")
   
   ind.name <- c("Celerity",
                 "Singularity",
                 "Relative singularity",
                 "MPD",
                 "NAO",
                 "Oceanic",
                 "Mediterranean",
                 "Northeast",
                 "Anticyclonic")
   
   png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/plot.cor.interannual/plot_",bv1,"_",bv2,".png"),width = 7,height = 3,units = "in",res=600)
   layout(matrix(1:3,1,3),widths = c(1,1,0.5),heights = c(1,1,1,1))
   
   par(mar=c(4,4,2,1))
   plot(1:4,corr.bv1[,1]*sign.bv1[,1],xlab="",xaxt="n",ylim=c(0.5,0.8),ylab="Correlation",pch=pch.bv1[,1],col=colo[i],main=nam2str(bv1))
   grid(nx=NA,ny=NULL)
   abline(v=1:4,lty=3,col="grey")
   axis(side = 1,at = 1:4,labels = sais)
   for(i in 1:ncol(corr.bv1)){
     points(1:4,corr.bv1[,i]*sign.bv1[,i],pch=pch.bv1[,i],col=colo[i],cex=2)
   }
   
   plot(1:4,corr.bv2[,1]*sign.bv2[,1],xlab="",xaxt="n",ylim=c(0.5,0.8),ylab="Correlation",pch=pch.bv2[,1],col=colo[i],main=nam2str(bv2))
   grid(nx=NA,ny=NULL)
   abline(v=1:4,lty=3,col="grey")
   axis(side = 1,at = 1:4,labels = sais)
   for(i in 1:ncol(corr.bv1)){
     points(1:4,corr.bv2[,i]*sign.bv2[,i],pch=pch.bv2[,i],col=colo[i],cex=2)
   }
   
   par(mar=c(0,1,0,0))
   plot.new()
   legend("center",legend = ind.name,col = colo,bty="n",pch=19,bg=2,cex=1.2)
   graphics.off()
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

# Trace les correlations propres à la fois en journalier et en annuel
plot.cor.descr.precip.clean <- function(k,bv="Isere-seul",sais="winter",rean,save=F){
  
  # Import des correlations
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.cor.descr.precip.clean/",bv,"_",sais,".Rdata"))
  corr.daily <- corr$corr.daily
  corr.yearly <- corr$corr.yearly
  
  # Signe des correlations
  sign <- unname(apply(corr.daily,2,function(v) quantile(v,probs=0.5,na.rm=T)/abs(quantile(v,probs=0.5,na.rm=T)))) # on passe en positif les correlations negatives pour mieux comparer
  ltyp <- sign; ltyp[ltyp==-1] <- 2
  
  # Couleurs
  colo <- c("black",
            "olivedrab3",
            "darkorange",
            "royalblue",
            "grey",
            "darkgoldenrod3",
            "darkorchid2",
            "red",
            "aquamarine2")
  
  # Graphique
  if(save){
    png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.cor.descr.precip.clean/plot_",bv,"_",sais,".png"),width = 600,height = 300,units = "px")
    layout(matrix(c(rep(1,2),2:3),nrow = 2,ncol = 2,byrow = T),widths = c(0.3,0.7),heights = c(0.1,0.9))
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,sais,cex=2,bg=2)
    }
  
  par(mar=c(ifelse(save,4,2),5,0,0))
  plot(1:nrow(corr.daily),corr.daily[,1]*sign[1],log="x",type="l",lty=ltyp[1],
       ylim=c(-0.2,1),xaxt="n",xlab="Smoothing length (days)",ylab="Correlation")
  grid(nx=NA,ny=NULL)
  axis(side = 1,at = c(1,10,100),labels = c(1,10,100))
  abline(v=c(1,10,100),lty=3,col="grey")
  for(i in 1:ncol(corr.daily)){
    lines(1:nrow(corr.daily),corr.daily[,i]*sign[i],lty=ltyp[i],col=colo[i],lwd=2)
  }
  
  par(mar=c(ifelse(save,4,2),1,0,0))
  plot(1:nrow(corr.yearly),corr.yearly[,1]*sign[1],type="l",lty=ltyp[1],
       ylim=c(-0.2,1),xaxt="n",yaxt="n",xlab="Smoothing length (years)",ylab="")
  grid(nx=NA,ny=NULL)
  axis(side = 1,at = seq(1,10,1),labels = seq(1,10,1))
  abline(v=seq(1,10,1),lty=3,col="grey")
  for(i in 1:ncol(corr.yearly)){
    lines(1:nrow(corr.yearly),corr.yearly[,i]*sign[i],lty=ltyp[i],col=colo[i],lwd=2)
  }
  rect(xleft = 0.8,ybottom = -0.2,xright = 1.2,ytop = 1,lty = 2)
  
  if(save) graphics.off()
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

# Scatterplot de 2 indicateurs VS precip par saison au pas de temps interannuel
scatterplot.descr.precip.interannual <- function(descr=c("dP","singnei"),bv="Isere-seul",k,dist,start="1950-01-01",end="2017-12-31",rean="ERA5",spazm=T){
  
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  dec <- which(substr(dates,6,7)=="12")
  year[dec] <- as.character(as.numeric(year[dec])+1)
  
  # Import indicateur
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  colnames(des) <- nam2str(descr,whole=T)
  for(i in 1:length(descr)){
    des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = 1,start = start,end = end,
                              standardize = F,rean = rean,threeday = F,period = "present")
  }
  
  # Import precip
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = spazm)
  
  # Moyennes interannuelles
  sais <- c("winter","spring","summer","autumn")
  sta <- c("12-01","03-01","06-01","09-01")
  N <- c(90,92,92,91)
  ind.sais <- vector(length=length(dates))
  
  for(i in 1:length(sais)){
      pos.sta <- which(substr(dates,6,10)== sta[i])
      if(sais[i]=="winter"){pos.sta <- pos.sta[-length(pos.sta)]}
      pos.all <- pos.sta
      for(j in 1:(N[i]-1)){ # saisons de longueur 90 jours
        pos.all <- sort(c(pos.all,pos.sta+j))
      }
    ind.sais[pos.all] <- i
  }
  
  precip.inter <- aggregate(precip,by=list(as.numeric(year),ind.sais),sum)
  colnames(precip.inter) <- c("year","season","val")
  precip.inter <- pivot_wider(precip.inter,names_from = "season",values_from = "val")
  precip.inter <- precip.inter[,-2] # on enleve saison 0
  precip.inter <- as.matrix(precip.inter[sort(as.matrix(precip.inter[,1]),index.return=T)$ix,-1]) # bon ordre
  
  ind.inter <- apply(des,2,function(v){aggregate(v,by=list(as.numeric(year),ind.sais),mean)})
  colnames(ind.inter[[1]]) <- colnames(ind.inter[[2]]) <- c("year","season","val")
  ind.inter <- lapply(ind.inter,function(mat){pivot_wider(mat,names_from = "season",values_from = "val")})
  
  des1.inter <- ind.inter[[1]][,-2]
  des1.inter <- as.matrix(des1.inter[sort(as.matrix(des1.inter[,1]),index.return=T)$ix,-1])
  des2.inter <- ind.inter[[2]][,-2]
  des2.inter <- as.matrix(des2.inter[sort(as.matrix(des2.inter[,1]),index.return=T)$ix,-1])

  # Graphique
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/scatterplot.descr.precip.interannual/scatterplot_",bv,"_",descr[1],"_",descr[2],ifelse(spazm,"_spazm",""),".png"),width = 8,height = 5,units = "in",res=400)
  par(pty="s",mfcol=c(2,4))
  for(i in 1:4){
  plot(des1.inter[,i],precip.inter[,i],pch=19,xlab=nam2str(descr[1],whole=T), ylab="Precipitation (mm)",
       main=paste0(sais[i]," (R = ",round(cor(na.omit(des1.inter[,i]),na.omit(precip.inter[,i])),2),")"))
  plot(des2.inter[,i],precip.inter[,i],pch=19,xlab=nam2str(descr[2],whole=T), ylab="Precipitation (mm)",
        main=paste0(sais[i]," (R = ",round(cor(na.omit(des2.inter[,i]),na.omit(precip.inter[,i])),2),")"))
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
  png(filename = paste0("2_Travail/1_Present/Rresults/scatterplot.nao.precip/plot_",bv,"_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
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
  png(filename = paste0("2_Travail/1_Present/Rresults/scatterplot.wp.precip/plot_",bv,"_wp=",wp,"_",start,"_",end,"_liss=",liss,".png"),width = 15,height = 10,units = "cm",res=300)
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