source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Calcul des champs composites pour une saison, pour une saison/année, et anomalie de cette année/année par rapport à saison
compute.composite.sais <- function(sais="winter",year="1995",k=1,var="hgt",start="1950-01-01",end="2017-12-31",rean){
  
  dates <- getdates(start,end)
  
  # Dates de la saison d'interet
  print("Dates")
  if(sais=="winter") {sta <- "12-01"; N <- 90}
  if(sais=="spring") {sta <- "03-01"; N <- 92}
  if(sais=="summer") {sta <- "06-01"; N <- 92}
  if(sais=="autumn") {sta <- "09-01"; N <- 91}
  
  pos.sta <- which(substr(dates,6,10)== sta)
  yr <- year
  
  if(sais=="winter"){
    pos.sta <- pos.sta[-length(pos.sta)] # on retire le dernier hiver incomplet
    yr <- as.character(as.numeric(yr)-1) # pour deb.year: on prend le premier decembre de l'annee d'avant
    } 
  
  pos.all <- pos.sta
  for(i in 1:(N-1)){ # saisons de longueur 90 jours
    pos.all <- sort(c(pos.all,pos.sta+i))
  }
  
  # Dates de la saison et de l'annee d'interet
  deb.year <- which(dates==paste0(yr,"-",sta))
  pos.year <- deb.year:(deb.year+N-1)
  
  # Import
  print("Import reanalyse")
  field <- getdata(k = k,day0 = start,day1 = end,rean = rean,large_win = F,small_win = F,all = T,var = var)
  gc()
  
  # Traitement
  print("Calcul composite")
  comp <- apply(field[,,pos.all],1:2,mean)
  comp.year <- apply(field[,,pos.year],1:2,mean)
  anom.year <- comp.year - comp
  gc()
  
  # Export
  res <- list(comp = comp, comp.year = comp.year, anom.year = anom.year)
  save(res, file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.composite.sais/composite_",sais,"_",year,"_k",k,"_",start,"_",end,".Rdata"))
  
}

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
  descr <- c("celnei","singnei","rsingnei","dP","dP_grad","cel","sing05","rsing05")
  
  # Import indicateur
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  colnames(des) <- descr
  for(i in 1:length(descr)){
  des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = 1,start = start,end = end,
                              standardize = F,rean = rean,threeday = F,period = "present")
  }
  
  # Import NAO
  NAO <- get.nao(start = substr(start,1,4),end = substr(end,1,4),sais = "all",daily = T)[,4]
  NAO_norm <- get.nao(start = substr(start,1,4),end = substr(end,1,4),sais = "all",daily = T,normalized = T)[,4]
  des <- cbind(des,NAO,NAO_norm)
  
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
  tmp <- get.ind.season(sais = sais,start = start,end = end)
  pos.season <- tmp$pos.season
  n.season <- tmp$n.season
  l.season <- tmp$l.season
  rm(tmp)
  
  # Tableaux propres pour calcul correlation
  precip.sais <- precip[pos.season]
  dim(precip.sais) <- c(l.season,n.season)
  
  corr.daily <- matrix(data = NA,nrow = l.season,ncol = ncol(des))
  colnames(corr.daily) <- colnames(des)
  
  corr.yearly <- matrix(data = NA,nrow = 10,ncol = ncol(des))
  colnames(corr.yearly) <- colnames(des)
  
  for(i in 1:ncol(des)){
    print(i)
    # Indicateur
    des.sais <- des[pos.season,i]
    dim(des.sais) <- c(l.season,n.season)
    
    # Calcul correlation journaliere
    for(j in 1:l.season){
      des.tmp <- des.sais
      precip.tmp <- precip.sais
      
      des.j <- as.vector(apply(des.tmp,2,function(v){rollapply(v,j,mean,na.rm=T)})) # na.rm = T car 2 valeurs de NA dans NAO_norm
      precip.j <- as.vector(apply(precip.tmp,2,function(v){rollapply(v,j,mean,na.rm=T)}))
      corr.daily[j,i] <- cor(des.j,precip.j,use = "pairwise.complete.obs")
    }
    
    # Calcul correlation annuelle
    des.tmp <- apply(des.sais,2,mean,na.rm=T)
    precip.tmp <- apply(precip.sais,2,mean,na.rm=T)
    
    for(j in 1:10){
      des.j <- rollapply(des.tmp,j,mean,na.rm=T)
      precip.j <- rollapply(precip.tmp,j,mean,na.rm=T)
      corr.yearly[j,i] <- cor(des.j,precip.j,use = "pairwise.complete.obs")
    }
  }
  
  corr <- list(corr.daily=corr.daily,corr.yearly=corr.yearly)
  save(corr,file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.cor.descr.precip.clean/",bv,"_",sais,ifelse(spazm,"_spazm",""),".Rdata"))
}

# Calcul des valeurs d'indicateur (et precip) moyennees par saison
compute.descr.interannual <- function(k,dist,bv="Isere",spazm=T,start="1950-01-01",end="2017-12-31",rean){
  
  dates <- getdates(start,end)
  year <- substr(dates,1,4)
  pos.dec <- which(substr(dates,6,7)=="12")
  year[pos.dec] <- as.character(as.numeric(year[pos.dec])+1) # on attribue les mois de decembre a l'annee suivante
  year[which(year==as.numeric(substr(end,1,4))+1)] <- substr(end,1,4)
  nyear <- length(unique(year))
  descr <- c("celnei","singnei","rsingnei","dP","cel","sing05","rsing05")
  
  # Import indicateur
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  colnames(des) <- descr
  for(i in 1:length(descr)){
    des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = 1,start = start,end = end,
                              standardize = F,rean = rean,threeday = F,period = "present")
  }
  
  # Import NAO
  NAO <- get.nao(start = substr(start,1,4),end = substr(end,1,4),sais = "all",daily = T)[,4]
  NAO_norm <- get.nao(start = substr(start,1,4),end = substr(end,1,4),sais = "all",daily = T,normalized = T)[,4]
  des <- cbind(des,NAO,NAO_norm)
  
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
  
  # Moyennes interannuelles saisonnieres ou non
  sais <- c("winter","spring","summer","autumn")
  sta <- c("12-01","03-01","06-01","09-01")
  N <- c(90,92,92,91)
  ind.sais <- vector(length=length(dates))
  
  for(i in 1:length(sais)){
    pos.sta <- which(substr(dates,6,10)== sta[i])
    if(sais[i]=="winter"){pos.sta <- pos.sta[-length(pos.sta)]}
    pos.all <- pos.sta
    for(j in 1:(N[i]-1)){ # longueur des saisons
      pos.all <- sort(c(pos.all,pos.sta+j))
    }
    ind.sais[pos.all] <- i
  }
  
  precip.inter <- aggregate(precip,by=list(as.numeric(year),ind.sais),sum)
  colnames(precip.inter) <- c("year","season","val")
  precip.inter <- pivot_wider(precip.inter,names_from = "season",values_from = "val")
  precip.inter <- precip.inter[,-2] # on enleve saison 0
  precip.inter <- as.matrix(precip.inter[sort(as.matrix(precip.inter[,1]),index.return=T)$ix,-1]) # bon ordre
  
  precip.inter.ann <- aggregate(precip,by=list(substr(dates,1,4)),sum)
  precip.inter <- cbind(precip.inter,precip.inter.ann[,2])
  
  ind.inter <- apply(des,2,function(v){aggregate(v,by=list(as.numeric(year),ind.sais),mean,na.rm=T)})
  for(i in 1:length(ind.inter)) colnames(ind.inter[[i]]) <- c("year","season","val")
  ind.inter <- lapply(ind.inter,function(mat){pivot_wider(mat,names_from = "season",values_from = "val")})
  for(i in 1:length(ind.inter)){
    ind.inter[[i]] <- as.matrix(ind.inter[[i]][sort(as.matrix(ind.inter[[i]][,1]),index.return=T)$ix,-c(1,2)])
  }
  
  ind.inter.ann <- apply(des,2,function(v){aggregate(v,by=list(substr(dates,1,4)),mean,na.rm=T)})
  for(i in 1:length(ind.inter)) ind.inter[[i]] <- cbind(ind.inter[[i]],ind.inter.ann[[i]][,2])
  
  # Export
  save(precip.inter,file = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/precip_interannual_",bv,"_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
  save(ind.inter,file = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  
}

# Calcul de la longueur moyenne (persistence) des types de temps
compute.length.wp <- function(agreg=T,start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  
  # Import des wp
  wp <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Calcul des longueurs moyennes
  wp.unique <- sort(unique(wp))
  len.wp <- list(length=length(wp.unique))
  
  for(i in 1:length(wp.unique)){
    
    # Longueur des sequences du type de temps
    len <- rle(wp==wp.unique[i])
    len <- len$lengths[len$values]
    
    # Dates associees au debut de la sequence
    dates.wp <- as.numeric(as.Date(dates[wp==wp.unique[i]]))
    tmp <- dates.wp[2:length(dates.wp)] - dates.wp[1:(length(dates.wp)-1)] # difference des dates: si 1, journees consecutives
    tmp <- c(2,tmp) # on selectionne forcement la premiere date de dates.wp, donc on lui attribue un chiffre different de 1
    dates.wp.unique <- as.character(as.Date(dates.wp[tmp!=1]))
    
    # Annee accociee au debut de la sequence
    ann.wp.unique <- substr(dates.wp.unique,1,4)
    
    # Saison associee au debut de la sequence
    season.wp.unique <- substr(dates.wp.unique,6,7)
    season.wp.unique[season.wp.unique %in% c("03","04","05")] <- "spring"
    season.wp.unique[season.wp.unique %in% c("06","07","08")] <- "summer"
    season.wp.unique[season.wp.unique %in% c("09","10","11")] <- "autumn"
    season.wp.unique[season.wp.unique %in% c("12","01","02")] <- "winter"
    
    # Stockage
    len.wp[[i]] <- list(len=len,
                        dates=dates.wp.unique,
                        year=ann.wp.unique,
                        season= season.wp.unique)
    print(mean(len))
  }
  
  # plot des tendances
  for(i in 1:length(wp.unique)){
    tmp <- aggregate(len.wp[[i]]$len,by=list(len.wp[[i]]$year,len.wp[[i]]$season),mean)
    plot(as.numeric(tmp[tmp[,2]=="autumn",1]),tmp[tmp[,2]=="autumn",3],pch=19,main=paste0("WP = ",wp.unique[i]," autumn"))
    grid()
    abline(lm(tmp[tmp[,2]=="autumn",3]~as.numeric(tmp[tmp[,2]=="autumn",1])))
    
    plot(as.numeric(tmp[tmp[,2]=="winter",1]),tmp[tmp[,2]=="winter",3],pch=19,main=paste0("WP = ",wp.unique[i]," winter"))
    grid()
    abline(lm(tmp[tmp[,2]=="winter",3]~as.numeric(tmp[tmp[,2]=="winter",1])))
    
    plot(as.numeric(tmp[tmp[,2]=="spring",1]),tmp[tmp[,2]=="spring",3],pch=19,main=paste0("WP = ",wp.unique[i]," spring"))
    grid()
    abline(lm(tmp[tmp[,2]=="spring",3]~as.numeric(tmp[tmp[,2]=="spring",1])))
    
    plot(as.numeric(tmp[tmp[,2]=="summer",1]),tmp[tmp[,2]=="summer",3],pch=19,main=paste0("WP = ",wp.unique[i]," summer"))
    grid()
    abline(lm(tmp[tmp[,2]=="summer",3]~as.numeric(tmp[tmp[,2]=="summer",1])))
    }
  
  # Interprétation des graphes de tendance
  # Pas de tendance sur la longueur des sequences oceaniques, 
  # si ce n'est une légère reduction de la duree des sequences oceaniques en automne.
  # Reduction de la duree du Med en hiver! (saison ou il se fait bouffer par l'oceanique, donc ca va dans le bon sens).
  # Surtout à partir des années 2000! Et aussi legerement au printems.
  # Rien sur le Nord Est.
  # Sequences anticycloniques legerement plus longues au printemps et plus long en été!
  
  # le plus interessant: Med qui raccourci en hiver (marqué recemment), et anticyclonique qui rallonge en ete.
}

# Corrplot entre les differents indicateurs pour une saison ou l'année entière
corrplot.descr <- function(k,dist,sais="year",start="1950-01-01",end="2017-12-31",rean,save=T){
  
  # Import des indicateurs interannuels
  ind <- c("sing05","dP","NAO_norm","Atlantic","Mediterranean","Northeast","Anticyclonic")
  season <- c("winter","spring","summer","autumn","year")
  pos <- which(season==sais)
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  
  # Traitement et tracé de la figure
  if(save){
    png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/corrplot.descr/corrplot_descr_",sais,"_",start,"_",end,".png"),width = 10,height = 6,units = "in",res = 600)
  }

    tab <- matrix(NA,nrow(ind.inter[[1]]),length(ind.inter))
    for(j in 1:length(ind.inter)) tab[,j] <- ind.inter[[j]][,pos]
    colnames(tab) <- names(ind.inter)
    tab <- tab[,ind] # selection des indicateurs
    tab.cor <- cor(tab,use="pairwise.complete.obs")
    colnames(tab.cor) <- rownames(tab.cor) <- nam2str(colnames(tab.cor),whole=T)
    
    # Corrplot
    col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061")) # definition de la palette par defeut pour pouvoir l'inverser ensuite
    corrplot(corr = tab.cor,method="circle",type="upper",diag=F,addCoef.col="black",col=rev(col2(200)),
             title=nam2str(sais),addgrid.col = "black",tl.col="black",tl.srt=45,mar=c(0,0,1,0),cex.main=2,oma=c(0,0,0,0),
             cl.align.text="l")

  if(save) graphics.off()
}

# Corrplot combiné: deux saisons
corrplot.descr.combine <- function(k,dist,sais=c("year","winter"),start="1950-01-01",end="2017-12-31",rean){

  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/corrplot.descr/corrplot_descr_",sais[1],"_",sais[2],"_",start,"_",end,".png"),width = 10,height = 4,units = "in",res = 600)
  par(mfrow=c(1,2))
  corrplot.descr(k = k,dist = dist,sais = sais[1],start = start,end = end,rean = rean,save = F)
  corrplot.descr(k = k,dist = dist,sais = sais[2],start = start,end = end,rean = rean,save = F)
  graphics.off()
}

# Corrplot entre les differents indicateurs pour les 4 saisons
corrplot.descr.season <- function(k,dist,start="1950-01-01",end="2017-12-31",rean){
  
  # Import des indicateurs interannuels
  ind <- c("sing05","dP","NAO","Atlantic","Mediterranean","Northeast","Anticyclonic")
  sais <- c("winter","spring","summer","autumn")
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  
  # Traitement et tracé de la figure
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/corrplot.descr/corrplot_descr_season_",start,"_",end,".png"),width = 10,height = 10,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(0,0,3,0))
  
  for(i in 1:length(sais)){
    
    # Extraction de la liste et calcul correlation
    tab <- matrix(NA,nrow(ind.inter[[1]]),length(ind.inter))
    for(j in 1:length(ind.inter)) tab[,j] <- ind.inter[[j]][,i]
    colnames(tab) <- names(ind.inter)
    tab <- tab[,ind] # selection des indicateurs
    tab.cor <- cor(tab,use="pairwise.complete.obs")
    colnames(tab.cor) <- rownames(tab.cor) <- nam2str(colnames(tab.cor),whole=T)
    
    # Corrplot
    corrplot(corr = tab.cor,method="circle",type="lower",diag=F,addCoef.col="black",
             title=sais[i],addgrid.col = "black",tl.col="black",tl.srt=45,mar=c(0,0,1,0),cex.main=2)
  }
  
  graphics.off()
  
  
}

# Import de la serie mensuelle de NAO
get.nao <- function(start="1950",end="2019",sais="all",daily=F,normalized=F){
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
    if(!normalized){
    nao <- read.csv(file = "2_Travail/Data/NAO/daily_NAO_NOAA.txt",sep="",skip= 4)
    }else{
      nao <- read.csv(file = "2_Travail/Data/NAO/daily_normalized_NAO_NOAA.txt",sep="",skip= 2)
    }
  }
  
  year <- seq(as.numeric(start),as.numeric(end))
  nao <- nao[nao[,"YEAR"] %in% year,]
  nao
}

# Carte composite des anomalies par wp ou groupement de wp, pour des composites deja calcules car trop lourds (ERA5)
map.composite.sais <- function(sais="winter",descr="sing05",years=c(NA,NA),k,start,end,rean,leg=T,win=T,iso=T){
  
  ann.deb <- as.numeric(substr(start,1,4))
  
  # Import des indicateurs interannuels (pour annee des min et max)
  des <- c("celnei","singnei","rsingnei","dP","cel","sing05","rsing05","NAO","Atlantic","Mediterranean","Northeast","Anticyclonic")
  season <- c("winter","spring","summer","autumn")
  
  if(is.na(years[1])){
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  des.inter <- ind.inter[[which(des==descr)]][,which(sais == season)]
  ann.min <- ann.deb + which.min(des.inter) - 1
  ann.max <- ann.deb + which.max(des.inter) - 1
  }else{
    ann.min <- years[1]
    ann.max <- years[2]
  }
  
  # Import des composites
  load(file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.composite.sais/composite_",sais,"_",ann.min,"_k",k,"_",start,"_",end,".Rdata"))
  res.min <- res
  
  load(file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.composite.sais/composite_",sais,"_",ann.max,"_k",k,"_",start,"_",end,".Rdata"))
  res.max <- res
  
  # Parametres graphiques
  if(rean=="ERA5"){
    lon <- seq(-20,30,0.25)
    lat <- seq(25,70,0.25)
    lon.wind <- seq(-30,30,0.25)
    lat.wind <- seq(27,70,0.25)
    delta=0.25
  }
  
  # Pour couleur et legende
  breaks <- seq(4850,6100,length.out = 12)
  N <- 11
  lab <- seq(4900,6100,200)
  pos.lab <- ecdf(seq(4850,6100))(lab)
  
  fen <- getinfo_window(k = k,large_win = F,small_win = F,all = F,rean = rean,var = "hgt")
  
  # Carte
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.composite.sais/composite_",sais,"_",descr,"_k",k,"_",start,"_",end,".png"),width=8,height=5,units = "in",res=300)
 
  layout(matrix(c(1:3,1,4:5),2,3,byrow=T))
  par(pty="s",mar=c(0,0,2,6))
  cex <- 1.5
  ann <- as.character(c(ann.min,ann.max))
  res <- list(res.min,res.max)
  
  # Composite saison toute periode
  image(lon,lat,res[[1]][[1]],xlim=c(-15,25),ylim=c(25,65),
        col=rev(brewer.pal(n = N, name = "RdBu")),
        xlab="",ylab="",main="",xaxt="n",yaxt="n",
        cex.axis=cex, cex.lab=cex, cex.main=cex,
        breaks = breaks)
  title(main=nam2str(sais),cex=3)
  data(wrld_simpl)
  plot(wrld_simpl, add = TRUE)
  points(6,45,col="red",pch=19)
  if(win) rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
  if(iso) contour(x=lon,y=lat,z=res[[1]][[1]], levels=breaks, drawlabels=F, lty=1, lwd=1, add=TRUE, col="black")
  if(leg){
    par(xpd=T)
    colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = lab,at =  pos.lab,
              vertical = T,cex=1.4,xlim=c(26,29),ylim=c(25,65),align = "l")
    par(xpd=F)
  }
  
  for(i in 1:2){
    for(j in c(3:2)){
      
      if(j==3){ 
        breaks <- seq(-200,200,length.out = 12)
        N <- 11
        lab <- seq(-200,200,50)
        pos.lab <- ecdf(seq(-200,200))(lab)
      }else{
        breaks <- seq(4850,6100,length.out = 12)
        N <- 11
        lab <- seq(4900,6100,200)
        pos.lab <- ecdf(seq(4850,6100))(lab)
      }

        image(lon,lat,res[[i]][[j]],xlim=c(-15,25),ylim=c(25,65),
              col=rev(brewer.pal(n = N, name = "RdBu")),
              xlab="",ylab="",main="",xaxt="n",yaxt="n",
              cex.axis=cex, cex.lab=cex, cex.main=cex,
              breaks = breaks)
        title(main=paste0(nam2str(sais)," ",ann[i],ifelse(j==3," (anomalies)","")),cex=3)
        
        data(wrld_simpl)
        plot(wrld_simpl, add = TRUE)
        points(6,45,col="red",pch=19)
        if(win) rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
        if(iso) contour(x=lon,y=lat,z=res[[i]][[j]], levels=breaks, drawlabels=F, lty=1, lwd=1, add=TRUE, col="black")
        if(leg){
          par(xpd=T)
          colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
                      labels = lab,at =  pos.lab,
                      vertical = T,cex=1.4,xlim=c(26,29),ylim=c(25,65),align = "l")
          par(xpd=F)
        }
    }
  }
  graphics.off()
}

# Cartes composite tous wp
map.wp <- function(wp=c(1,2),rean="ERA5",anomalies=T,start="1950-01-01",end="2019-12-31"){
  
  # Source fct map.composite.wp.light dans article1_AB.R
  source('I:/ongoing_Documents/2_Travail/Rcode/article1_AB.R', encoding = 'UTF-8')
  
  l <- length(wp)
  
  # Figure
  png(filename = paste0(get.dirstr(1,rean,"present"),"map.wp/map_comp_wp_",paste(wp,collapse = "_"),"_",rean,ifelse(anomalies,"_anomalies",""),"_",start,"_",end,".png"),width = l*2.5,height = 4,units = "in",res = 600)
  layout(matrix(1:(l+1),ncol=l+1,nrow=1),heights = 1,widths = c(rep(1,l),0.3))
  
  # Cartes composite wp
  par(mar=c(1,1,1,1))
  for(i in 1:length(wp)){
    map.composite.wp.light(wp = wp[i],k = 1,start = start,end = end,rean = rean,
                           leg = F,win = T,let = F,iso = F,agreg = T,anomalies = anomalies,wind = T)
  }
  
  # Legende
  N <- 11
  if(anomalies){
    leg <- seq(-200,200,50)
  }else{
    leg <- seq(4900,6100,200)
  }
  
  par(pty="m",mar=c(1,0,1,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(0,1),ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = paste0("        ",leg),at =  seq(0, 1,length.out = length(leg)),
              vertical = T,xlim = c(0.1,0.4),ylim = c(0.25,0.75),cex=1.4)
  
  graphics.off()
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
      plot.cor.descr.precip.clean(k = 1,bv = bv[j],sais = sais[i],rean = rean,save = F,spazm = T)
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

# Regroupement des plot de correlation en fonction du lissage pour un seul BV
plot.cor.one <- function(bv="Isere",spazm=T,rean){
  
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/plot.cor.descr.precip.clean/plot_combine_",bv,ifelse(spazm,"_spazm",""),".png"),width = 8,height = 6,units = "in",res=600)
  layout(matrix(c(rep(1,2),rep(2,2),3:6,rep(7,2),rep(8,2),9:16,rep(17,4)),nrow = 6,ncol = 4,byrow = T),widths = c(0.34,0.66,0.34,0.66),heights = c(0.15,1,0.15,1,0.15,0.2))
  
  # Graphiques
  sais <- c("winter","spring","summer","autumn")
  ind <- c("NAO_norm","Atlantic","Mediterranean","Northeast","Anticyclonic","dP","sing05")
  colo <- c("black","blue","red","darkgreen","dimgrey","purple","darkorange")
  
  for(i in c(1,3)){
    # Nom saison x 2
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,nam2str(sais[i]),cex=2,font=2)
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,nam2str(sais[i+1]),cex=2,font=2)
    
    plot.cor.descr.precip.clean(ind = ind,colo = colo,k = 1,bv = bv,sais = sais[i],spazm = spazm,rean = rean,save = F)
    plot.cor.descr.precip.clean(ind = ind,colo = colo,k = 1,bv = bv,sais = sais[i+1],spazm = spazm,rean = rean,save = F)
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
  colo <- c("purple","darkorange","black","blue","red","darkgreen","dimgrey")
  ind.name <- c("MPD","Singularity","NAO","Atlantic","Mediterranean","Northeast","Anticyclonic")
  lwidth <- c(4,4,2,2,2,2,2)
  
  par(mar=c(0,3,0,0))
  plot.new()
  legend(0.2,1,legend = ind.name,col = colo,bty="n",lty=1,lwd=lwidth,ncol=5,bg=2,cex=1.2)
  
  graphics.off()
}

# Graphique des correlations interannuelles par saison pour les deux BVs
plot.cor.interannual <- function(bv1="Isere-seul",bv2="Drac-seul",rean){
  
  # Import des correlations
  sais <- c("winter","spring","summer","autumn")
  bv <- c(bv1,bv2)
  
  # BV1
   corr.bv1 <- matrix(NA,4,12)
    for(i in 1:length(sais)){
      load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/compute.cor.descr.precip.clean/",bv1,"_",sais[i],".Rdata"))
      corr.bv1[i,] <- corr$corr.yearly[1,]
    }
   colnames(corr.bv1) <- colnames(corr$corr.yearly)
   sign.bv1 <- apply(corr.bv1,2,sign)
   pch.bv1 <- sign.bv1; pch.bv1 <- apply(pch.bv1,2,function(v){v[v==-1] <- 19;v[v==1] <- 15;return(v)})
   
   # BV2
   corr.bv2 <- matrix(NA,4,12)
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
             "orange",
             "yellow",
             "purple",
             "grey",
             "darkgoldenrod3",
             "darkorchid2",
             "red",
             "aquamarine2")
   
   ind.name <- c("Celerity",
                 "Singularity",
                 "Relative singularity",
                 "MPD",
                 "cel",
                 "sing",
                 "rqing",
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

# Graphique des correlations interannuelles par saison pour un seul BV
plot.cor.interannual.one <- function(bv="Isere",spazm=T,rean){
  
  # Indicateurs qu'on garde
  ind <- c("NAO_norm","Atlantic","Mediterranean","Northeast","Anticyclonic","dP","dP_grad","sing05")
  sign <- c(-1,1,-1,-1,-1,1,1,-1) # pour etre coherent avec figure lissage
  
  # Import des correlations
  sais <- c("winter","spring","summer","autumn")
  
  corr.bv <- matrix(NA,4,14)
  for(i in 1:length(sais)){
    load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/compute.cor.descr.precip.clean/",bv,"_",sais[i],ifelse(spazm,"_spazm",""),".Rdata"))
    corr.bv[i,] <- corr$corr.yearly[1,]
  }
  colnames(corr.bv) <- colnames(corr$corr.yearly)
  corr.bv <- corr.bv[,ind]
  corr.bv <- t(apply(corr.bv,1,function(v) v*sign))

  lty.bv <- sign; lty.bv[lty.bv==-1] <- 2
  
  # Graphique
  colo <- c("black","blue","red","darkgreen","dimgrey","purple","mediumpurple","darkorange")
  colo.leg <- c("purple","mediumpurple","darkorange","black","blue","red","darkgreen","dimgrey")
  
  ind.name <- c("NAO","Atlantic","Mediterranean","Northeast","Anticyclonic","MPD","MPD_grad","Singularity") # MPD et sing à la fin pour être au dessus des autres lignes
  ind.name.leg <- c("MPD","MPD_grad","Singularity","NAO","Atlantic","Mediterranean","Northeast","Anticyclonic") # Au debut dans la legende
  
  lwidth <- c(2,2,2,2,2,4,4,4)
  lwidth.leg <- c(4,4,4,2,2,2,2,2)
  
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/plot.cor.interannual/plot_",bv,ifelse(spazm,"_spazm",""),".png"),width = 6,height = 4,units = "in",res=600)
  layout(matrix(1:2,1,2),widths = c(0.7,0.3))
  
  par(mar=c(3,4,0,0))
  plot(1:4,corr.bv[,1],type="n",col=colo[1],lty=lty.bv[1],ylim=c(-0.2,1),xlab="",xaxt="n",ylab="Correlation")
  grid(nx=NA,ny=NULL)
  abline(v=1:4,lty=3,col="grey")
  axis(side = 1,at = 1:4,labels = sais)
  for(i in 1:ncol(corr.bv)){
    lines(1:4,corr.bv[,i],lty=lty.bv[i],lwd=lwidth[i],col=colo[i])
  }
  
  par(mar=c(0,1,0,0))
  plot.new()
  legend("center",legend = ind.name.leg,col = colo.leg,bty="n",lwd=lwidth.leg,bg=2)
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
plot.cor.descr.precip.clean <- function(ind,colo,k,bv="Isere-seul",sais="winter",spazm=T,rean,save=F){
  
  # Import des correlations
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.cor.descr.precip.clean/",bv,"_",sais,ifelse(spazm,"_spazm",""),".Rdata"))
  corr.daily <- corr$corr.daily
  corr.daily <- corr.daily[,ind]
  corr.yearly <- corr$corr.yearly
  corr.yearly <- corr.yearly[,ind]
  
  # Signe des correlations
  #sign <- sign(corr.yearly[1,]) # on garde le signe de la correlation interannuelle
  sign <- unname(apply(corr.yearly,2,function(v) quantile(v,probs=0.5,na.rm=T)/abs(quantile(v,probs=0.5,na.rm=T)))) # on passe en positif les correlations negatives pour mieux comparer
  ltyp <- sign; ltyp[ltyp==-1] <- 2
  
  # Graphique
  if(save){
    png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.cor.descr.precip.clean/plot_",bv,"_",sais,ifelse(spazm,"_spazm",""),".png"),width = 600,height = 300,units = "px")
    layout(matrix(c(rep(1,2),2:3),nrow = 2,ncol = 2,byrow = T),widths = c(0.3,0.7),heights = c(0.1,0.9))
    par(mar=c(0,0,0,0))
    plot.new()
    text(0.5,0.5,sais,cex=2,bg=2)
    }
  
  par(mar=c(ifelse(save,4,2),5,0,0))
  plot(1:nrow(corr.daily),corr.daily[,1]*sign[1],log="x",type="n",lty=ltyp[1],
       ylim=c(-0.2,1),xaxt="n",xlab="Smoothing length (days)",ylab="Correlation")
  grid(nx=NA,ny=NULL)
  axis(side = 1,at = c(1,10,100),labels = c(1,10,100))
  abline(v=c(1,10,100),lty=3,col="grey")
  for(i in 1:ncol(corr.daily)){
    lines(1:nrow(corr.daily),corr.daily[,i]*sign[i],lty=ltyp[i],col=colo[i],lwd=ifelse(ind[i] %in% c("dP","sing05"),4,2))
  }
  
  par(mar=c(ifelse(save,4,2),1,0,0))
  plot(1:nrow(corr.yearly),corr.yearly[,1]*sign[1],type="n",lty=ltyp[1],
       ylim=c(-0.2,1),xaxt="n",yaxt="n",xlab="Smoothing length (years)",ylab="")
  grid(nx=NA,ny=NULL)
  axis(side = 1,at = seq(1,10,1),labels = seq(1,10,1))
  abline(v=seq(1,10,1),lty=3,col="grey")
  for(i in 1:ncol(corr.yearly)){
    lines(1:nrow(corr.yearly),corr.yearly[,i]*sign[i],lty=ltyp[i],col=colo[i],lwd=ifelse(ind[i] %in% c("dP","sing05"),4,2))
  }
  #rect(xleft = 0.8,ybottom = -0.2,xright = 1.2,ytop = 1,lty = 2)
  
  if(save) graphics.off()
}

# Plot correlation cumul saisonnier VS fort plusieurs quantiles pour une saison
plot.cor.precip <- function(bv,sais,nbdays=1,start="1950-01-01",end="2017-12-31",spazm){
  
  dates <- getdates(start,end)
  
  # Import
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = spazm)
  
  # Mise en forme
  tmp <- get.ind.season(sais = sais,start = start,end = end)
  pos.season <- tmp$pos.season
  n.season <- tmp$n.season
  l.season <- tmp$l.season
  rm(tmp)
  
  precip <- precip[pos.season]
  dim(precip) <- c(l.season,n.season)
  
  # Si plusieurs jours
  if(nbdays!=1){
    precip <- apply(precip,2,function(v){rollapply(v,nbdays,mean)})
    l.season <- nrow(precip)
  }
  
  # Traitement
  cumul.precip <- apply(precip,2,sum)
  
  quant <- seq(0.05,0.95,0.05)
  quant.precip <- matrix(nrow=length(quant),ncol=n.season)
  cor.precip <- vector(length=length(quant))
  for(i in 1:length(quant)){
    tmp <- apply(precip,2,function(v){quantile(v,probs=quant[i])})
    quant.precip[i,] <- tmp
    cor.precip[i] <- cor(tmp,cumul.precip)
  }
  
  # Graphique
  png(filename = paste0("2_Travail/0_Present/Rresults/plot.cor.precip/plot_cor_precip_",bv,"_",sais,"_",nbdays,"days_",start,"_",end,ifelse(spazm,"_spazm",""),".png"),width = 8,height = 6,units = "in",res=600)
  plot(quant,cor.precip,type="l",ylim=c(0,1),lwd=2,xlab=paste0("Percentile of ",nbdays,"-day precipitation"),
       ylab=paste0("Correlation with ",nam2str(sais)," Precipitation Accumulation"),main=nam2str(sais))
  grid()
  graphics.off()
}

# Code temporaire pour tracer des densites d'indicateurs sur une periode en journalier
plot.density.descr <- function(descr="sing05",start="1994-12-01",end="1995-02-28",k,dist,rean){
  
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = "1950-01-01",end = "2017-12-31",
                        standardize = F,rean = rean,threeday = F,period = "present")
  
  dates <- getdates(start = "1950-01-01",end = "2017-12-31")
  des.winter <- des[which()]
  des.per <- des[which(dates==start):which(dates==end)]
  plot(density(des.per))
  des.per <- des[which(dates=="2013-12-01"):which(dates=="2014-01-28")]
  lines(density(des.per),col="red")
  min(des.per)
  quantile(des.per,probs=0.1)
  
  
}

# Run
run.letter <- function(type=1){
  
  if(type==1){# compute.cor.descr.precip.clean
    bv <- "Isere"# c("Isere","Isere-seul","Drac-seul")
    sais <- c("winter","spring","summer","autumn")
    
    for(i in 1:length(bv)){
      for(j in 1:length(sais)){
        compute.cor.descr.precip.clean(k = 1,dist = "TWS",bv = bv[i],sais = sais[j],spazm = F,start = start,end=end,rean = rean)
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
  
  # compute.cor.descr.precip et plot.cor.descr.precip
  if(type==3){
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
  
  # scatterplot.descr.nao
  if(type==4){
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
  
  # scatterplot.descr.wp
  if(type==6){
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
  
}

# Scatterplot pour lettre: 2 ind VS precip et ind VS ind
scatterplot.all <- function(k,dist,bv,descr=c("dP","sing05"),start="1950-01-01",end="2017-12-31",spazm=T,rean){
  
  ind <- c("sing05","dP","NAO","Atlantic","Anticyclonic")
  sais <- c("winter","spring","summer","autumn")
  
  # Import
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/precip_interannual_",bv,"_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  
  # Scatterplot ind/precip par saison
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/scatterplot.all/scatterplot_all_",bv,"_",descr[1],"_",descr[2],"_",start,"_",end,ifelse(spazm,"_spazm",""),".png"),width = 8,height = 6,units = "in",res=600)
  par(pty="s",mfrow=c(4,4),mar=c(3,2,2,0))
  
  des1 <- ind.inter[[which(names(ind.inter)==descr[1])]]
  for(i in 1:4){
    plot(des1[,i],precip.inter[,i],pch=19,xlab="", ylab="Precipitation (mm)",
         main=paste0(sais[i]," (R = ",round(cor(des1[,i],precip.inter[,i],use="pairwise.complete.obs"),2),")"))
    title(xlab=paste0(nam2str(descr[1],whole=T),ifelse(descr[1]=="MPD"," (m)","")),line=2)
    if(i==1) mtext(text = "a)",side = 3,line=1,adj = -1)
  }
  
  des2 <- ind.inter[[which(names(ind.inter)==descr[2])]]
  for(i in 1:4){
    plot(des2[,i],precip.inter[,i],pch=19,xlab="", ylab="Precipitation (mm)",
         main=paste0(sais[i]," (R = ",round(cor(des2[,i],precip.inter[,i],use="pairwise.complete.obs"),2),")"))
    title(xlab=paste0(nam2str(descr[2],whole=T),ifelse(descr[2]=="MPD"," (m)","")),line=2)
    if(i==1) mtext(text = "b)",side = 3,line=1,adj = -1)
    }
  
  # Scatterplot des differentes combinaisons d'indicateurs
  comb <- list(
    c("dP","NAO"),
    c("dP","sing05"),
    c("dP","Atlantic"),
    c("dP","Anticyclonic"),
    c("sing05","NAO"),
    c("sing05","Atlantic"),
    c("sing05","Anticyclonic")
  )
  
  pal <- colorRampPalette(c("white","blue"))
  pt1 <- which.min(des1[,1])
  colo <- rep("black",nrow(des1));colo[pt1]="red"
  
  for(i in 1:length(comb)){
    print(i)
    des1 <- ind.inter[[which(names(ind.inter)==comb[[i]][1])]][,5] # on ne prend que l'hiver
    des2 <- ind.inter[[which(names(ind.inter)==comb[[i]][2])]][,5]
    plot(des1,des2,pch=21,xlab="",ylab=paste0(nam2str(comb[[i]][2],whole=T),ifelse(comb[[i]][2]=="MPD"," (m)","")),
         main=paste0("winter"," (R = ",round(cor(des1,des2,use="pairwise.complete.obs"),2),")"),
         bg=pal(n = length(des1))[rank(precip.inter[,5])],col=colo)
    title(xlab=paste0(nam2str(comb[[i]][1],whole=T),ifelse(comb[[i]][1]=="MPD"," (m)","")),line=2)
    if(i==1) mtext(text = "c)",side = 3,line=1,adj = -1)
  }
  
  # Legende
  plot.new()
  title(main = "Precipitations")
  mini <- min(precip.inter[,1],na.rm=T);maxi <- max(precip.inter[,1],na.rm=T);mil <- (maxi-mini)/2+mini
  leg <- paste0(as.character(round(c(mini,mil,maxi),0))," mm")
  color.legend(0,0.1,0.2,0.9,legend = leg,gradient="y",rect.col = pal(length(des1)),align="rb",cex=0.8,font=2,pos=4)
  
  graphics.off()
}

# Scatterplot entre descripteurs
scatterplot.descr <- function(k,dist,start="1950-01-01",end="2017-12-31",rean){
  
  # Import indicateurs
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
       
  # Scatterplot des differentes combinaisons d'indicateurs
  comb <- list(
    c("sing05","dP"),
    c("sing05","Atlantic"),
    c("dP","Atlantic"),
    c("dP","NAO_norm")
  )
  
  year <- seq(as.numeric(substr(start,1,4)),as.numeric(substr(end,1,4)))
  colo <- rep("black",length(year))
  pos.1 <- which(year==1989);pos.2 <- which(year==1995)
  colo[pos.1] <- "red";colo[pos.2] <- "deepskyblue3"
  
  letter <- c("a)","b)","c)","d)")
  
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/scatterplot.descr/scatterplot_descr_",start,"_",end,".png"),width = 7,height = 2,units = "in",res=600)
  par(mfrow=c(1,4),mar=c(3,4,1.5,0),pty="s")
  for(i in 1:length(comb)){
    print(i)
    des1 <- ind.inter[[which(names(ind.inter)==comb[[i]][1])]][,1] # on ne prend que l'hiver
    des2 <- ind.inter[[which(names(ind.inter)==comb[[i]][2])]][,1]
    plot(des1,des2,pch=19,xlab="",ylab="",col=colo,main=paste0("R = ",round(cor(des1,des2,use="pairwise.complete.obs"),2)))
    grid()
    points(des1,des2,pch=19,col=colo)
    title(xlab=paste0(nam2str(comb[[i]][1],whole=T),ifelse(comb[[i]][1]=="dP"," (m)","")),line=2)
    title(ylab=paste0(nam2str(comb[[i]][2],whole=T),ifelse(comb[[i]][2]=="dP"," (m)","")),line=2)
    mtext(text = letter[i],side = 3,line=1,adj = -0.2)
  }
  
  graphics.off()
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

# Scatterplot d'un max/quantile de descripteur et des max de precipitations par saison
scatterplot.descr.precip.extr <- function(bv="Isere",descr=c("dP","sing05"),k,dist,nbdays=1,rean="ERA5",start="1950-01-01",end="2017-12-31",qua=F,prob=0.95,dat=F){

  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Imports
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = T)
  des <- matrix(NA,length(dates),length(descr))
  for(i in 1:length(descr)){
    des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,
                              standardize = F,rean = rean,threeday = F,period = "present")
  }
  
  # Traitement et graphique
  season <- c("winter","spring","summer","autumn")
  
  png(filename = paste0(get.dirstr(k = k,rean = rean,period = "present"),"scatterplot.descr.precip.extr/scatterplot_precip_extr_",bv,"_",ifelse(qua,paste0("q",prob*100),"max"),"_",descr[1],"_",descr[2],"_",nbdays,"day_",ifelse(dat,"dat_",""),start,"_",end,".png"),width = 8,height = 4,units = "in",res = 600)
  par(mfrow=c(2,4),pty="s",mar=c(4,4,2,0))
  
  for(i in 1:length(descr)){
    for(j in 1:length(season)){
      print(season[j])
      
      if(dat){
        ind <- get.ind.max.sais(sais = season[j],wp = "all",nbdays = nbdays,start = start,end = end,bv = bv,spazm = T)
        precip.j <- precip[ind]
        des.j <- des[ind,i]
      }else{
      tmp <- get.ind.season(sais = season[j],start = start,end = end)
      
      # Traitement precip
      precip.j <- precip[tmp$pos.season]
      precip.j[tmp$pos.NA] <- NA
      dim(precip.j) <- c(tmp$l.season,tmp$n.season)
      precip.j <- apply(precip.j,2,max,na.rm=T)
      
      # Traitement descr
      des.j <- des[tmp$pos.season,i]
      des.j[tmp$pos.NA] <- NA
      dim(des.j) <- c(tmp$l.season,tmp$n.season)
      if(!qua){
        if(descr[i]=="dP"){
          des.j <- apply(des.j,2,max,na.rm=T)
        }else{des.j <- apply(des.j,2,min,na.rm=T)}
      }else{des.j <- apply(des.j,2,function(v){quantile(v,probs=ifelse(descr[i]=="dP",prob,1-prob),na.rm=T)})}
      }
      # Graphique
      corre <- round(cor(des.j,precip.j,use = "pairwise.complete.obs"),2)
      plot(des.j,precip.j,pch=19,xlab=paste0(ifelse(!qua,ifelse(descr[i]=="dP","Max","Min"),paste0("q",ifelse(descr[i]=="dP",prob,1-prob)*100))," ",nbdays,"-day ",nam2str(descr[i],whole=T)),ylab=paste0("Max ",nbdays,"-day Precipitation (mm)"),main=paste0(nam2str(season[j])," (R=",corre,")"),ylim=range(precip))
      grid();par(new=T)
      plot(des.j,precip.j,pch=19,xlab=paste0(ifelse(!qua,ifelse(descr[i]=="dP","Max","Min"),paste0("q",ifelse(descr[i]=="dP",prob,1-prob)*100))," ",nbdays,"-day ",nam2str(descr[i],whole=T)),ylab=paste0("Max ",nbdays,"-day Precipitation (mm)"),main=paste0(nam2str(season[j])," (R=",corre,")"),ylim=range(precip))
    }
  }
  graphics.off()
}

# Scatterplot de 2 indicateurs VS precip par saison au pas de temps interannuel
scatterplot.descr.precip.interannual <- function(descr=c("dP","sing05"),bv="Isere",k,dist,start="1950-01-01",end="2017-12-31",rean="ERA5",spazm=T){
  
  des <- c("celnei","singnei","rsingnei","dP","cel","sing05","rsing05","NAO","Atlantic","Mediterranean","Northeast","Anticyclonic")
  sais <- c("winter","spring","summer","autumn")
  
  # Import
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/precip_interannual_",bv,"_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  des1.inter <- ind.inter[[which(des==descr[1])]]
  des2.inter <- ind.inter[[which(des==descr[2])]]
  
  # Xlim et Ylim
  xlim1 <- range(des1.inter,na.rm=T)
  xlim2 <- range(des2.inter,na.rm=T)
  ylim <- range(precip.inter[,1:4],na.rm=T)
  
  # Affichage des annees min et max de precip, des1 et des2
  ann.start <- as.numeric(substr(start,1,4))
  print("Precipitation")
  print(apply(precip.inter,2,function(v){paste0("Annee min = ",which.min(v)+ann.start-1," et Annee max = ",which.max(v)+ann.start-1)}))
  
  print(nam2str(descr[1],whole=T))
  print(apply(des1.inter,2,function(v){paste0("Annee min = ",which.min(v)+ann.start-1," et Annee max = ",which.max(v)+ann.start-1)}))
  
  print(nam2str(descr[2],whole=T))
  print(apply(des2.inter,2,function(v){paste0("Annee min = ",which.min(v)+ann.start-1," et Annee max = ",which.max(v)+ann.start-1)}))
  
  # Graphique
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/scatterplot.descr.precip.interannual/scatterplot_",bv,"_",descr[1],"_",descr[2],ifelse(spazm,"_spazm",""),".png"),width = 8,height = 4,units = "in",res=400)
  par(pty="s",mfcol=c(2,4),mar=c(3,4,2,0))
  for(i in 1:4){
    if(i==1){
      year <- seq(as.numeric(substr(start,1,4)),as.numeric(substr(end,1,4)))
      colo <- rep("black",length(year))
      pos.1 <- which(year==1955) #pos.1 <- which(year==1989);pos.2 <- which(year==1995)
      colo[pos.1] <- "deepskyblue3" #colo[pos.1] <- "red";colo[pos.2] <- "blue"
    }else{colo <- "black"}
    
    plot(des1.inter[,i],precip.inter[,i],pch=19,main=paste0(nam2str(sais[i])," (R = ",round(cor(na.omit(des1.inter[,i]),na.omit(precip.inter[,i])),2),")"),
         col=colo,xlab="",ylab="",xlim=xlim1,ylim=ylim)
    grid()
    points(des1.inter[,i],precip.inter[,i],pch=19,col=colo)
    title(xlab = paste0(nam2str(descr[1],whole=T),ifelse(descr[1]=="dP"," (m)","")),line=2)
    title(ylab = "Precipitation (mm)",line=2)
    
    plot(des2.inter[,i],precip.inter[,i],pch=19,main=paste0(nam2str(sais[i])," (R = ",round(cor(na.omit(des2.inter[,i]),na.omit(precip.inter[,i])),2),")"),
         col=colo,xlab="",ylab="",xlim=xlim2,ylim=ylim)
    grid()
    points(des2.inter[,i],precip.inter[,i],pch=19,col=colo)
    title(xlab = paste0(nam2str(descr[2],whole=T),ifelse(descr[2]=="dP"," (m)","")),line=2)
    title(ylab = "Precipitation (mm)",line=2)
  }
  graphics.off()
}

# Scatterplot de 2 indicateurs VS precip par saison au pas de temps interannuel
scatterplot.descr.precip.interannual.color <- function(descr=c("dP","singnei"),colo=c("NAO","Atlantic"),bv="Isere-seul",k,dist,start="1950-01-01",end="2017-12-31",rean="ERA5",spazm=T){
  
  des <- c("celnei","singnei","rsingnei","dP","cel","sing05","rsing05","NAO","Atlantic","Mediterranean","Northeast","Anticyclonic")
  sais <- c("winter","spring","summer","autumn")
  
  # Import
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/precip_interannual_",bv,"_",ifelse(spazm,"spazm_",""),start,"_",end,".Rdata"))
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  
  des1.inter <- ind.inter[[which(des==descr[1])]]
  des2.inter <- ind.inter[[which(des==descr[2])]]
  des1.colo <- ind.inter[[which(des==colo[1])]]
  des2.colo <- ind.inter[[which(des==colo[2])]]
  
  # Affichage des annees min et max de precip, des1 et des2
  ann.start <- as.numeric(substr(start,1,4))
  print("Precipitation")
  print(apply(precip.inter,2,function(v){paste0("Annee min = ",which.min(v)+ann.start-1," et Annee max = ",which.max(v)+ann.start-1)}))
  
  print(nam2str(descr[1],whole=T))
  print(apply(des1.inter,2,function(v){paste0("Annee min = ",which.min(v)+ann.start-1," et Annee max = ",which.max(v)+ann.start-1)}))
  
  print(nam2str(descr[2],whole=T))
  print(apply(des2.inter,2,function(v){paste0("Annee min = ",which.min(v)+ann.start-1," et Annee max = ",which.max(v)+ann.start-1)}))
  # Graphique
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/scatterplot.descr.precip.interannual/scatterplot_color_",bv,"_",descr[1],"_",descr[2],ifelse(spazm,"_spazm",""),".png"),width = 8,height = 5,units = "in",res=400)
  par(pty="s",mfcol=c(2,4))
  for(i in 1:4){
    plot(des1.inter[,i],precip.inter[,i],pch=19,xlab=nam2str(descr[1],whole=T), ylab="Precipitation (mm)",
         main=paste0(sais[i]," (R = ",round(cor(na.omit(des1.inter[,i]),na.omit(precip.inter[,i])),2),")"),
         col=getcol(des1.colo[,i]))
    addscale(des1.colo[,i])
    plot(des2.inter[,i],precip.inter[,i],pch=19,xlab=nam2str(descr[2],whole=T), ylab="Precipitation (mm)",
         main=paste0(sais[i]," (R = ",round(cor(na.omit(des2.inter[,i]),na.omit(precip.inter[,i])),2),")"),
         col=getcol(des2.colo[,i]))
    addscale(des2.colo[,i])
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

# Scatterplot de 2 indicateurs au pas de temps interannuel, pour l'annee et pour l'hiver (retour reviewers article 1)
scatterplot.two.descr <- function(descr=c("dP","NAO_norm"),k,dist,spazm=T,start="1950-01-01",end="2017-12-31",rean){
  
  # Import
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.descr.interannual/descr_interannual_k",k,"_",dist,"_",start,"_",end,".Rdata"))
  nam <- names(ind.inter)
  
  # Traitement et graphique
  des1.ann <- ind.inter[[descr[1]]][,5]
  des1.win <- ind.inter[[descr[1]]][,1]
  cor.ann <- round(cor(des1.ann,des2.ann,use = "pairwise.complete.obs"),2)
  
  des2.ann <- ind.inter[[descr[2]]][,5]
  des2.win <- ind.inter[[descr[2]]][,1]
  cor.win <- round(cor(des1.win,des2.win,use = "pairwise.complete.obs"),2)
  
  
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/scatterplot.two.descr/scatterplot_",descr[1],"_",descr[2],"_k",k,"_",dist,"_",start,"_",end,".png"),width=6,height=3,unit="in",res=600)
  par(mfrow=c(1,2),pty="s",mar=c(4,4,2,1))
  plot(des1.ann,des2.ann,pch=19,xlab=nam2str(descr[1],unit=T),ylab=nam2str(descr[2],unit=T),main=paste0("Year (R=",cor.ann,")"))
  grid();par(new=T)
  plot(des1.ann,des2.ann,pch=19,xlab=nam2str(descr[1],unit=T),ylab=nam2str(descr[2],unit=T),main=paste0("Year (R=",cor.ann,")"))
  
  plot(des1.win,des2.win,pch=19,xlab=nam2str(descr[1],unit=T),ylab=nam2str(descr[2],unit=T),main=paste0("Winter (R=",cor.win,")"))
  grid();par(new=T)
  plot(des1.win,des2.win,pch=19,xlab=nam2str(descr[1],unit=T),ylab=nam2str(descr[2],unit=T),main=paste0("Winter (R=",cor.win,")"))
  
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

# Comprehension du fort cumul de precip sur l'hiver 1954/1955
understant.winter.1955 <- function(bv="Isere",start="1950-01-01",end="2017-12-31",spazm=T){
  
  dates <- getdates(start,end)
  yr <- unique(substr(dates,1,4))
  nyear <- length(yr)
  
  # Import precip
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = spazm)
  
  # Import indicateurs: singularite et MPD
  dP <- get.descriptor(descriptor = "dP",k = 1,dist = "TWS",nbdays = 1,start = start,end = end,standardize = F,rean = "ERA5")
  sing05 <- get.descriptor(descriptor = "sing05",k = 1,dist = "TWS",nbdays = 1,start = start,end = end,standardize = F,rean = "ERA5")
  cel <- get.descriptor(descriptor = "cel",k = 1,dist = "TWS",nbdays = 1,start = start,end = end,standardize = F,rean = "ERA5")
  
  # Import types de temps
  wp <- get.wp(nbdays = 1,start = start,end = end,agreg=T,spazm=spazm)
  
  # Mise en forme hiver
  sta <- "12-01"
  N <- 90
  
  pos.sta <- which(substr(dates,6,10)== sta)
  pos.sta <- pos.sta[-length(pos.sta)]
  yr <- yr[-1]
  nyear <- nyear-1
  
  pos.all <- pos.sta
  for(i in 1:(N-1)){ # saisons de longueur 90 jours
    pos.all <- sort(c(pos.all,pos.sta+i))
  }
  
  precip.winter <- precip[pos.all]
  dim(precip.winter) <- c(N,nyear)
  colnames(precip.winter) <- yr
  
  dP.winter <- dP[pos.all]
  dim(dP.winter) <- c(N,nyear)
  colnames(dP.winter) <- yr
  
  sing.winter <- sing05[pos.all]
  dim(sing.winter) <- c(N,nyear)
  colnames(sing.winter) <- yr
  
  wp.winter <-wp[pos.all]
  dim(wp.winter) <- c(N,nyear)
  colnames(wp.winter) <- yr
  
  # Distributions
  boxplot(precip.winter,col="royalblue")
  boxplot(dP.winter,col="red")
  boxplot(sing.winter,col="darkorange")
  
  table(wp.winter[,5])/nrow(wp.winter)
  plot(yr,apply(wp.winter,2,function(v){sum(v==1)}),col="blue",pch=19);grid()
  plot(yr,apply(wp.winter,2,function(v){sum(v==2)}),col="red",pch=19);grid()
  plot(yr,apply(wp.winter,2,function(v){sum(v==5)}),col="grey",pch=19);grid()
  plot(yr,apply(wp.winter,2,function(v){sum(v==8)}),col="black",pch=19);grid()
  
  # Nombres de jours superieurs a 20mm
  nb.20 <- apply(precip.winter,2,function(v){sum(v>20)})
  plot(yr,nb.20,pch=19,xlab="year",ylab="Number of days > 20mm");grid()
  
  nb.40 <- apply(precip.winter,2,function(v){sum(v>40)})
  plot(nb.40)
  
  nb.60 <- apply(precip.winter,2,function(v){sum(v>60)})
  plot(nb.60)
  
  # Cumul de precip en dessous d'un seuil de MPD
  med <- median(dP.winter)
  pos <- apply(dP.winter,2,function(v){which(v<med)})
  res <- NULL
  for(i in 1:ncol(dP.winter)){
    res[i] <- sum(precip.winter[pos[[i]],i])
  }
  plot(yr,res,pch=19)
  
  # Precip par mois
  precip.month <- aggregate(precip,by=list(substr(dates,6,7),substr(dates,1,4)),sum)
  dP.month <- aggregate(dP,by=list(substr(dates,6,7),substr(dates,1,4)),mean)
  sing.month <- aggregate(sing05,by=list(substr(dates,6,7),substr(dates,1,4)),mean)
  cel.month <- aggregate(cel,by=list(substr(dates,6,7),substr(dates,1,4)),mean)
  
  par(pty="s")
  png(filename = "2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/dP_precip_december.png",width = 4,height = 4,units = "in",res=300)
  plot(dP.month[dP.month[,1]=="12",3],precip.month[precip.month[,1]=="12",3],pch=19,xlab="MPD (m)",ylab="Precipitation (mm)",main="December")
  points(dP.month[dP.month[,1]=="12",3][6],precip.month[precip.month[,1]=="12",3][6],pch=19,col="red")
  graphics.off()
  
  png(filename = "2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/dP_precip_january.png",width = 4,height = 4,units = "in",res=300)
  plot(dP.month[dP.month[,1]=="01",3],precip.month[precip.month[,1]=="01",3],pch=19,xlab="MPD (m)",ylab="Precipitation (mm)",main="January")
  points(dP.month[dP.month[,1]=="01",3][6],precip.month[precip.month[,1]=="01",3][6],pch=19,col="red")
  #points(dP.month[dP.month[,1]=="01",3][46],precip.month[precip.month[,1]=="01",3][46],pch=19,col="blue")
  graphics.off()
  
  png(filename = "2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/dP_precip_february.png",width = 4,height = 4,units = "in",res=300)
  plot(dP.month[dP.month[,1]=="02",3],precip.month[precip.month[,1]=="02",3],pch=19,xlab="MPD (m)",ylab="Precipitation (mm)",main="February")
  points(dP.month[dP.month[,1]=="02",3][6],precip.month[precip.month[,1]=="02",3][6],pch=19,col="red")
  graphics.off()
  
  png(filename = "2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/sing_precip_december.png",width = 4,height = 4,units = "in",res=300)
  plot(sing.month[sing.month[,1]=="12",3],precip.month[precip.month[,1]=="12",3],pch=19,xlab="Singularity",ylab="Precipitation (mm)",main="December")
  points(sing.month[sing.month[,1]=="12",3][6],precip.month[precip.month[,1]=="12",3][6],pch=19,col="red")
  graphics.off()
  
  png(filename = "2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/sing_precip_january.png",width = 4,height = 4,units = "in",res=300)
  plot(sing.month[sing.month[,1]=="01",3],precip.month[precip.month[,1]=="01",3],pch=19,xlab="Singularity",ylab="Precipitation (mm)",main="January")
  points(sing.month[sing.month[,1]=="01",3][6],precip.month[precip.month[,1]=="01",3][6],pch=19,col="red")
  graphics.off()
  
  png(filename = "2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/sing_precip_february.png",width = 4,height = 4,units = "in",res=300)
  plot(sing.month[sing.month[,1]=="02",3],precip.month[precip.month[,1]=="02",3],pch=19,xlab="Singularity",ylab="Precipitation (mm)",main="February")
  points(sing.month[sing.month[,1]=="02",3][6],precip.month[precip.month[,1]=="02",3][6],pch=19,col="red")
  graphics.off()
  
  plot(cel.month[cel.month[,1]=="01",3],precip.month[precip.month[,1]=="01",3],pch=19,xlab="Celerity",ylab="Precipitation (mm)",main="January")
  points(cel.month[cel.month[,1]=="01",3][6],precip.month[precip.month[,1]=="01",3][6],pch=19,col="red")

  
  # Boxplot du mois de janvier
  # Marais barométrique: faible MPD et singulier!
  # Fortes précipitations concentrées du 10 au 18 avec fort MPD. Le reste est un marais!
  
  # Mise en forme janvier
  dP.jan <- dP[substr(dates,6,7)=="01"]
  dim(dP.jan) <- c(31,nyear+1)
  colnames(dP.jan) <- c("1950",yr)
  
  sing.jan <- sing05[substr(dates,6,7)=="01"]
  dim(sing.jan) <- c(31,nyear+1)
  colnames(sing.jan) <- c("1950",yr)
  
  cel.jan <- cel[substr(dates,6,7)=="01"]
  dim(cel.jan) <- c(31,nyear+1)
  colnames(cel.jan) <- c("1950",yr)
  
  precip.jan <- precip[substr(dates,6,7)=="01"]
  dim(precip.jan) <- c(31,nyear+1)
  colnames(precip.jan) <- c("1950",yr)
  
  plot(seq(as.Date("1955-01-01"),as.Date("1955-01-31"),"days"),precip.jan[,6],type="l",xlab="Date",ylab="Precipitation (mm)")
  grid()
  
  boxplot(dP.jan,col="red")
  
  max=apply(dP.jan,2,function(v){max(v)-min(v)})
  plot(max)
  
  min=apply(dP.jan,2,min)
  plot(min)
  
  plot(density(dP.jan[,6]),ylim=c(0,0.007))
  for(i in 1:ncol(dP.jan)) lines(density(dP.jan[,i]))
  lines(density(dP.jan[,6]),col="red",lwd=3)
  
  q90 <- apply(dP.jan,2,function(v){quantile(v,probs=0.99)})
  plot(q90)
  
  q10 <- apply(dP.jan,1,quantile,probs=0.1)
  med <- apply(dP.jan,1,quantile,probs=0.5)
  q90 <- apply(dP.jan,1,quantile,probs=0.9)
  
  q10 <- apply(sing.jan,1,quantile,probs=0.1)
  med <- apply(sing.jan,1,quantile,probs=0.5)
  q90 <- apply(sing.jan,1,quantile,probs=0.9)
  
  plot(seq(as.Date("1955-01-01"),as.Date("1955-01-31"),"days"),sing.jan[,6],type="l",xlab="Date",ylab="Singularity")
  lines(seq(as.Date("1955-01-01"),as.Date("1955-01-31"),"days"),q10,lty=3)
  lines(seq(as.Date("1955-01-01"),as.Date("1955-01-31"),"days"),q90,lty=3)
  lines(seq(as.Date("1955-01-01"),as.Date("1955-01-31"),"days"),med,lty=3)
  lines(seq(as.Date("1955-01-01"),as.Date("1955-01-31"),"days"),sing.jan[,8],col="blue")
  
  plot(apply(dP.jan,2,function(v){sum(v>400 & v<500)}),pch=19)
  
  dP.jan.liss <- apply(dP.jan,2,rollapply,width=10,mean)
  precip.jan.liss <- apply(precip.jan,2,rollapply,width=10,sum)
  plot(dP.jan.liss[,6],precip.jan.liss[,6])
  
  # Autocorrelation
  autocorr <- apply(dP.jan,2,function(v){min(acf(v,plot=F)$acf)})
  autocorr.pos <- apply(dP.jan,2,function(v){which.min(acf(v,plot=F)$acf)}) 
  
  # Le fait que les jours de forts MPD soient à la suite
  
  png("2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/dP_january55.png",width = 6,height = 5,units = "in",res=300)
  par(mar=c(5,5,5,5))
  plot(1:31,dP.jan[,6],type="l",xlab="Days",ylab="MPD (m)",ylim=c(100,1000),main="MPD")
  points(11:14,dP.jan[11:14,6],pch=19,col="red")
  par(new=T)
  plot(1:31,precip.jan[,6],ylim=rev(c(0,100)),type="h",col="blue",lwd=5,yaxt="n",xlab="",ylab="")
  axis(side = 4,at=seq(0,60,20),labels =seq(0,60,20))
  mtext(text = "Precipitation (mm)",side = 4,line = 3,at = 30)
  graphics.off()
  
  png("2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/sing05_january55.png",width = 6,height = 5,units = "in",res=300)
  par(mar=c(5,5,5,5))
  plot(1:31,sing.jan[,6],type="l",xlab="Days",ylab="Singularity",ylim=c(0.1,0.4),main="Singularity")
  points(11:14,sing.jan[11:14,6],pch=19,col="red")
  par(new=T)
  plot(1:31,precip.jan[,6],ylim=rev(c(0,100)),type="h",col="blue",lwd=5,yaxt="n",xlab="",ylab="")
  axis(side = 4,at=seq(0,60,20),labels =seq(0,60,20))
  mtext(text = "Precipitation (mm)",side = 4,line = 3,at = 30)
  graphics.off()
  
  png("2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/cel_january55.png",width = 6,height = 5,units = "in",res=300)
  par(mar=c(5,5,5,5))
  plot(1:31,cel.jan[,6],type="l",xlab="Days",ylab="Celerity",ylim=c(0.1,0.6),main="Celerity")
  points(11:14,cel.jan[11:14,6],pch=19,col="red")
  par(new=T)
  plot(1:31,precip.jan[,6],ylim=rev(c(0,100)),type="h",col="blue",lwd=5,yaxt="n",xlab="",ylab="")
  axis(side = 4,at=seq(0,60,20),labels =seq(0,60,20))
  mtext(text = "Precipitation (mm)",side = 4,line = 3,at = 30)
  graphics.off()
  
  logic <- dP.jan[,6]>400
  max(rle(logic)$lengths[which(rle(logic)$values)])
  
  tmp=apply(dP.jan,2,function(v){max(rle(v>450)$length[which(rle(v>450)$values)])})
  plot(1:68,tmp)
  
  tmp=apply(sing.jan,2,function(v){max(rle(v<0.2)$length[which(rle(v<0.2)$values)])})
  plot(1:68,tmp)
  
  tmp=apply(dP.jan,2,function(v){v>500})
  tmp1=apply(sing.jan,2,function(v){v<0.21})
  tmp2=apply(cel.jan,2,function(v){v<0.15})
  
  tmp3=tmp
  for(i in 1:ncol(tmp3)){
    tmp3[,i]=tmp[,i] & tmp1[,i] & tmp2[,i]
  }
  
  tmp4=apply(tmp3,2,function(v){max(rle(v)$length[which(rle(v)$values)])})
  plot(1:68,tmp4)
  
  # Celerite
  png("2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/cel_min_january.png",width = 6,height = 5,units = "in",res=300)
  cele <- apply(cel.jan,2,min,na.rm=T)
  plot(1950:2017,cele,pch=19,xlab="Years",ylab="Min january celerity",main="Celerity")
  grid()
  points(1950:2017,cele,pch=19)
  graphics.off()
  
  png("2_Travail/0_Present/ERA5/Rresults/overall/k1/understand.winter.1955/cel_q10_january.png",width = 6,height = 5,units = "in",res=300)
  cele <- apply(cel.jan,2,function(v){quantile(v,probs=0.1,na.rm=T)})
  plot(1950:2017,cele,pch=19,xlab="Years",ylab="Q10 january celerity",main="Celerity")
  grid()
  points(1950:2017,cele,pch=19)
  graphics.off()
  
  # Lissage sur 4 jours
  cele.liss <- apply(cel.jan,2,function(v){rollapply(v,4,mean)})
  cele <- apply(cele.liss,2,min,na.rm=T)
  plot(1950:2017,cele,pch=19,xlab="Years",ylab="Min january celerity")
  
  
  # Importance d'une faible celerite et position du flux exacte.
  
  
  
}
