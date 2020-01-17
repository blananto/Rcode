source('~/2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Carte des geopotentiels aux 4 coins du plan celnei TWS - celnei RMSE
map.corner <- function(k, rean){
  
  # Dates situees aux coins
  dates <- c("1981-12-12","1992-04-16","1975-07-09","1953-10-19")
  lettre <- c("a)","b)","c)","d)")
  
  # Cartes
  pdf(file = paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/map.corner/map_corner_celnei.pdf"),width = 10,height = 10)
  par(mfcol=c(4,3),mar=c(3,2,3,2))
  for(i in 1:4){
    map.geo(date = dates[i],rean = rean,k = k,nbdays = 3,save = F,win = T)
  }
  graphics.off()
  
}

# Saisonnalite de dP et d'un indicateur TWS et RMSE
plot.sais.all <- function(descr,k,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  # Calcul dP
  load.nc(rean)
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean)
  dates <- getdates(start,end)
  dP <- apply(geo,3,function(x) max(x)-min(x))
  dP <- rollapply(dP,nbdays,mean)
  
  # Import descr
  descr.tws <- get.descriptor(descr,k,"TWS",nbdays,start,end,F,rean)
  descr.rmse <- get.descriptor(descr,k,"RMSE",nbdays,start,end,F,rean)
  
  # dates
  dates.all <- dates[-((length(dates)-nbdays+2):length(dates))]
  dates.chro <- getdates(start = "1950-01-01",end = "1956-01-01")
  ind.chro <- match(dates.chro,dates)
  
  # Figure
  pdf(file = paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/plot.sais.all/plot_sais_all_",descr,".pdf"),width = 8,height = 6)
  par(mfrow=c(2,3),mar=c(2.5,2,2.5,2))
  plot.sais.quantile(dates = dates.all,vec = dP,label = "Pressure Gradient (m)")
  title("Pressure gradient")
  mtext("a)",side=3,adj=0,line=0.8,font=1)
  plot.sais.quantile(dates = dates.all,vec = descr.tws,label = paste0(descr," TWS"))
  title(paste0(descr," TWS"))
  mtext("b)",side=3,adj=0,line=0.8,font=1)
  plot.sais.quantile(dates = dates.all,vec = descr.rmse,label = paste0(descr," RMSE"))
  title(paste0(descr," RMSE"))
  mtext("c)",side=3,adj=0,line=0.8,font=1)
  plot.chronique(dates = dates.chro,vec = dP[ind.chro],label = "Pressure Gradient (m)",liss = T,extr = T)
  mtext("d)",side=3,adj=0,line=0.8,font=1)
  plot.chronique(dates = dates.chro,vec = descr.tws[ind.chro],label = paste0(descr," TWS"),liss = T,extr = T)
  mtext("e)",side=3,adj=0,line=0.8,font=1)
  plot.chronique(dates = dates.chro,vec = descr.rmse[ind.chro],label = paste0(descr," RMSE"),liss = T,extr = T)
  mtext("f)",side=3,adj=0,line=0.8,font=1)
  graphics.off()
  
}

# Plan des indicateurs colorie par densite, gradient de pression, et saison
plot.empir.dP.sais <- function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F),coin=F){
  
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
  
  # Points aleatoires dans le plan
  ale <- sample(1:length(descr1),62)
  
  # Points au coin du plan
  pos.coin <- which(getdates(start,end) %in% c("1981-12-12","1975-07-09","1992-04-16","1953-10-19"))
  
  # Points des extremes de precipitations
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  
  # Creation du plot
  pdf(file=paste0(path1,"plot.empir.dP.sais",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,".pdf"),width=8,height=3.5) # manip substr pour enlever le std
  par(mfrow=c(1,3))
  
  gamme <- c(0,4445)
  param <- param[,c("nbnei")]
  
  # Densite
  plot(descr1,descr2,
       col=getcol(param,range = gamme),
       xlab=paste0(namdescr[1]," ",dist[1]),
       ylab=paste0(namdescr[2]," ",dist[2]),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
  addscale(vec = c(param,gamme))
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.8,
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
       paste0(round(min(param,na.rm=T),2),"-",round(max(param,na.rm=T),2)))
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white") 
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  mtext("a)",side=3,adj=0,line=0.8,font=1)
  
  # Gradient de pression
  geo <- getdata(k = k[1],day0 = start,day1 = end,rean = rean[1]) 
  deltaP <- apply(geo,3,function(x) max(x)-min(x))
  deltaP <- rollapply(deltaP,nbdays,mean)
  
  plot(descr1,descr2,
       col=getcol(deltaP),
       xlab=paste0(namdescr[1]," ",dist[1]),
       ylab=paste0(namdescr[2]," ",dist[2]),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
  addscale(vec = round(deltaP,0))
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.8,
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2*0.95,
       paste0(round(min(deltaP,na.rm=T),0),"-",round(max(deltaP,na.rm=T),0)))
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white")
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  mtext("b)",side=3,adj=0,line=0.8,font=1)
  
  # Saison
  dates <- getdates(start,end)
  sea <- dates
  sea[which(substr(dates,6,7) %in% c("06","07","08"))] <- 1
  sea[which(substr(dates,6,7) %in% c("09","10","11"))] <- 2
  sea[which(substr(dates,6,7) %in% c("12","01","02"))] <- 3
  sea[which(substr(dates,6,7) %in% c("03","04","05"))] <- 4
  sea <- as.numeric(sea)
  colo <- c("red","darkorange","blue","olivedrab3")
    
  plot(descr1,descr2,
         col=colo[sea],
         xlab=paste0(namdescr[1]," ",dist[1]),
         ylab=paste0(namdescr[2]," ",dist[2]),
         ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2))
  legend("topleft",c("summer","autumn","winter","spring"),col=colo,pch=1,bty="n",ncol=2)
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white")
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  mtext("c)",side=3,adj=0,line=0.8,font=1)
  
  graphics.off()
  
}

# Desaisonalisation de dP et d'un indicateur
plot.desais.dP.descr <- function(descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean){
  
  dates <- getdates(start,end)
  if(nbdays!=1) dates <- dates[-((length(dates)-nbdays+2):(length(dates)))]
  
  # Calcul de dP
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean) 
  deltaP <- apply(geo,3,function(x) max(x)-min(x))
  deltaP <- rollapply(deltaP,nbdays,mean)
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean)
  
  # Traitement: desaisonalisation et quantiles
  tab <- cbind(deltaP,des)
  qua <- as.data.frame(matrix(NA,nrow(tab),4))
  
  for(i in 1:2){
  sais <- aggregate(tab[,i],by=list(substr(dates,6,10)),mean)
  pos <- match(substr(dates,6,10),sais[,1])
  sais.chro <- sais[pos,2]
  desais <- tab[,i]-sais.chro
  qua[,i] <- ecdf(tab[,i])(tab[,i])*100
  qua[,i+2] <- ecdf(desais)(desais)*100
  }
  qua <- qua[,c(1,3,2,4)]
  colnames(qua) <- c("deltaP","desais deltaP",descr,paste0("desais ",descr))
  
  # Mise en forme
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  tab.final <- pivot_longer(data = qua[ind.extr,],1:4,names_to = "name",values_to = "quantile")
  tab.final$name <- factor(tab.final$name,colnames(qua))
  
  # Boxplot
  ggplot(tab.final, aes(x=name, y=quantile, fill=name)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0,0,0,1),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "none")+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_x_discrete(labels=c(expression(paste(Delta,"p")),expression(paste("Sea. adj. ",Delta,"p")),descr,paste("Sea. adj. ",descr)))+
    scale_fill_manual(values=c("purple","purple","cornflowerblue","cornflowerblue"))+
    geom_vline(xintercept=2.5, linetype="dashed")+
    xlab("")+
    ylab("Percentile (%)")
  
  ggsave(filename = paste0("2_Travail/20CR/Rresults/overall/k",k,"/plot.desais.dP.descr/plot_desais_dP",descr,"_",nbdays,"day_",start,"_",end,"_",bv,".png"),width = 12,height = 6,units="cm",dpi = 200)
  graphics.off()
}


