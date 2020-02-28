source('~/2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')


# classification manuelle des circulations selon leur distance TWS ou RMSE
classif.perso <- function(k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  # Import des distances
  dist <- score.to.mat(k,dist,nbdays,start,end,rean)
  
  # Sequences d'interet et centres predefinis
  des <- get.dP(k,nbdays,start,end,rean)
  ind <- sort(des,decreasing=T,index.return=T)$ix[1:62]
  centers <- c(43,47,58) # position predefinie des sequences centres (NW, SW, N)
  
  # Graphique des distances
  for(i in 1:length(centers)){
    di <- dist[ind[centers[i]],ind[-centers[i]]]
    plot(sort(di),rep(1,length(di)))
    tail(sort(di))
  }
  
}

# Combine deux fonctions en une
combine.functions <- function(fun,descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean){
  
  # plot.desais.descr
  if(fun==1){
    png(filename = paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/plot.desais.descr/plot_desais_",descr[1],"_",descr[2],"_",
                          dist,"_member",member,"_k",k,"mean",nbdays,"day_",start,"_",end,".png"),width = 550,height = 250,units = "px")
    plot.desais.1 <- plot.desais.descr(descr[1],k,dist,nbdays,start,end,rean)
    plot.desais.2 <- plot.desais.descr(descr[2],k,dist,nbdays,start,end,rean)
    grid.arrange(plot.desais.1,plot.desais.2,nrow=1,ncol=2)
    graphics.off()
  }
  
  # plot.empir.dp.sais
  if(fun==2){
    png(file=paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/plot.empir.dP.sais-CV/plot_combine.png"),width=8,height=9,units = "in",res=1200) # manip substr pour enlever le std
    layout(matrix(c(rep(1,3),2:4,rep(5,3),6:8,rep(9,3),10:12),nrow = 6,ncol = 3,byrow = T),widths = rep(1,3),heights = c(0.1,1.2,0.1,1.2,0.1,1.2))
    par(pty="s")
    nam <- c("Celerity","Singularity","Relative singularity")
    let <- list(c("a)","b)","c)"),c("d)","e)","f)"),c("g)","h)","i)"))
    
    for(i in 1:length(descr)){
      # Nom de l'indicateur
      par(mar=c(0,0,0,0),pty="m")
      plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
      text(1,0.5,nam[i],cex=2,font=2)
      
      # Graphiques
      par(mar=c(3,4,0,3),pty="s")
      plot.empir.dP.sais(c(rean,rean),c(k,k),c(descr[i],descr[i]),dist,nbdays,start,end,coin=ifelse(descr[i]=="celnei",T,F),save=F,let=let[[i]])
    }
    graphics.off()
  }
  
  # plot.quant.descr
  if(fun==3){
    a <- plot.quant.descr(descr[1],nbdays,start,end,rean,save=F)
    b <- plot.quant.descr(descr[2],nbdays,start,end,rean,save=F)
    c <- plot.quant.descr(descr[3],nbdays,start,end,rean,save=F)
    d <- get_legend(a)
    final <- ggarrange(a, b, c, as_ggplot(d), nrow=2,ncol=2, legend = "none",
                       labels = c("Celerity","Singularity","Relative singularity"),hjust = c(-3.3,-2.2,-1.1))
    ggsave(filename = paste0("2_Travail/",rean,"/Rresults/overall/plot.quant.descr/plot_combine.png"),plot = final,width = 10,height = 6)
    graphics.off()
  }
  
}

# Calcul de dP
get.dP <- function(k,nbdays,start="1950-01-01",end="2011-12-31",rean){
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean) 
  des <- apply(geo,3,function(x) max(x)-min(x))
  des <- rollapply(des,nbdays,mean)
}

# Carte des geopotentiels aux 4 coins du plan celnei TWS - celnei RMSE
map.corner <- function(k, rean){
  
  # Dates situees aux coins
  dates <- c("1981-12-12","1992-04-16","1975-07-09","1953-10-19")
  lettre <- c("a)","b)","c)","d)")
  
  # Parametres de la legende
  if(k==1){
    leg <- as.character(seq(4900,6100,200))
    N=11
  } else{
    leg <- as.character(seq(-300,400,100))
    N=7
  }
  
  # Cartes
  png(file = paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/map.corner/map_corner_celnei.png"),width = 8,height = 14,units = "in",res = 1200)
  layout(matrix(c(1:12,rep(13,3)),nrow = 5,ncol = 3,byrow = T),widths = rep(1,3),heights = c(rep(1,4),0.5))

  for(i in 1:4){
    for(j in 1:3){
      date <- as.character(as.Date(dates[i])+j-1)
      map.geo(date = date,rean = rean,k = k,nbdays = 1,save = F,win = T,let = ifelse(j==1,lettre[i],F),leg=F)
    }
  }
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = leg,at =  seq(0, 1,length.out = length(leg)),
              vertical = F,xlim = c(0.7,1.3),ylim = c(0.3,0.6),cex=1.4)
  text(x = 1,y = 0.75,"Geopotential height (m)",cex=1.6)
  graphics.off()
  
}

# Carte des geopotentiels de l'extreme et d'une sequence seche avec leur 100e analogue
map.extr.dry <- function(k, rean){
  
  # Recherche des dates
  dates <- NULL
  
  dates[1] <- getdates()[get.ind.extr(nbre = 2,ref = "1950-01-01",nbdays = 1)[2]]
  dates[2] <- getdates()[get.ana(date = dates[1],rank = 100,ref = "1950-01-01",k = k,dist = "TWS",nbdays = 1,rean=rean)$ind]
  tmp <- which(get.precip(1)==0)
  dates[3] <- getdates()[sample(tmp,1)]
  dates[4] <- getdates()[get.ana(date = dates[3],rank = 100,ref = "1950-01-01",k = k,dist = "TWS",nbdays = 1,rean=rean)$ind]
  
  lettre <- c("a)","b)","c)","d)")
  
  # Parametres de la legende
  if(k==1){
    leg <- as.character(seq(4900,6100,200))
    N=11
  } else{
    leg <- as.character(seq(-300,400,100))
    N=7
  }
  
  # Cartes
  png(file = paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/map.extr.dry/map_extr_dry.png"),width = 6,height = 7,units = "in",res = 1200)
  layout(matrix(c(1:4,rep(5,2)),nrow = 3,ncol = 2,byrow = T),widths = rep(1,2),heights = c(rep(1,2),0.5))
  
  for(i in 1:4){
      map.geo(date = dates[i],rean = rean,k = k,nbdays = 1,save = F,win = T,let = lettre[i],leg=F)
  }
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = leg,at =  seq(0, 1,length.out = length(leg)),
              vertical = F,xlim = c(0.7,1.3),ylim = c(0.3,0.6),cex=1.4)
  text(x = 1,y = 0.75,"Geopotential height (m)",cex=1.6)
  graphics.off()
  
}

# Desaisonnalisation d'un indicateur (dP possible)
plot.desais.descr <- function(descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean){
  
  dates <- getdates(start,end)
  if(nbdays!=1) dates <- dates[-((length(dates)-nbdays+2):(length(dates)))]
  
  if(descr=="dP"){
    des <- get.dP(k,nbdays,start,end,rean)
  }else{
    # Import de l'indicateur
    des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean)
  }
  
  # Traitement: desaisonalisation et quantiles
  desais <- compute.desais(des,dates)
  
  qua <- as.data.frame(matrix(NA,length(des),2))
  qua[,1] <- ecdf(des)(des)*100
  qua[,2] <- ecdf(desais)(desais)*100
  colnames(qua) <- c("Original","Seasonal adjusted")
  
  # Mise en forme
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  tab.final <- pivot_longer(data = qua[ind.extr,],1:2,names_to = "name",values_to = "quantile")
  tab.final$name <- factor(tab.final$name,colnames(qua))
  
  # Boxplot
  colo <- ifelse(rep(descr=="dP",2),rep("grey",2),rep("cornflowerblue",2))
  main <- ifelse(descr=="dP","Pressure gradient",paste0(descr," ",dist))
  
  plot.desais.descr <- ggplot(tab.final, aes(x=name, y=quantile, fill=name)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1,1,0,1),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=10),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "none")+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=colo)+
    xlab("")+
    ylab("Percentile (%)")+
    ggtitle(main)
  
  plot.desais.descr
}

# plot la correlation entre le gradient de pression et un indicateur
plot.dP.descr <- function(descriptor,k,dist,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  # Import
  load.nc(rean)
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean)
  descr <- get.descriptor(descriptor = descriptor,k = k,dist = dist,nbdays = nbdays,
                          start = start,end = end,standardize = F,rean = rean)
  dates <- getdates(start,end)
  length(dates) <- length(descr)
  sais.descr <- aggregate(descr,by=list(substr(dates,6,10)),mean)[,2]
  desais.descr <- compute.desais(descr,dates)
  
  # dP
  dP <- get.dP(k,nbdays,start,end,rean)
  sais.dP <- aggregate(dP,by=list(substr(dates,6,10)),mean)[,2]
  desais.dP <- compute.desais(dP,dates)
  
  # Correlation
  cor <- round(cor(dP,descr)^2,2)
  lm  <- lm(descr~dP)
  
  cor.sais <- round(cor(sais.dP,sais.descr)^2,2)
  lm.sais  <- lm(sais.descr~sais.dP)
  
  cor.desais <- round(cor(desais.dP,desais.descr)^2,2)
  lm.desais  <- lm(desais.descr~desais.dP)
  
  # Graphiques
  png(file = paste0(get.dirstr(k,rean),"plot.dP.descr/plot_dP_",descriptor,"_",dist,"_member",member,
                    "_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width = 8,height = 3.5,units = "in",res = 1200)
  par(pty="s",mfrow=c(1,3))
  
  # Brut
  sea <- dates
  sea[which(substr(sea,6,7) %in% c("06","07","08"))] <- 1
  sea[which(substr(sea,6,7) %in% c("09","10","11"))] <- 2
  sea[which(substr(sea,6,7) %in% c("12","01","02"))] <- 3
  sea[which(substr(sea,6,7) %in% c("03","04","05"))] <- 4
  sea <- as.numeric(sea)
  colo <- c("red","darkorange","blue","olivedrab3")
  
  plot(dP,descr,col=colo[sea],xlab="Pressure gradient (m)",ylab=paste0(descriptor," ",dist),main="Original")
  abline(lm$coefficients[1],lm$coefficients[2],lwd=2)
  text(quantile(dP,0.1),quantile(descr,0.001),paste0("R² = ",cor),font=2)
  legend("topright",c("summer","autumn","winter","spring"),col=colo,pch=1,bty="n",ncol=2,cex=0.9)
  
  
  # Sais
  sea <- sort(unique(substr(dates,6,10)))
  sea[which(substr(sea,1,2) %in% c("06","07","08"))] <- 1
  sea[which(substr(sea,1,2) %in% c("09","10","11"))] <- 2
  sea[which(substr(sea,1,2) %in% c("12","01","02"))] <- 3
  sea[which(substr(sea,1,2) %in% c("03","04","05"))] <- 4
  sea <- as.numeric(sea)
  colo <- c("red","darkorange","blue","olivedrab3")
  
  plot(sais.dP,sais.descr,col=colo[sea],xlab="Pressure gradient (m)",ylab=paste0(descriptor," ",dist),main="Interannual")
  abline(lm.sais$coefficients[1],lm.sais$coefficients[2],lwd=2)
  text(quantile(sais.dP,0.3),quantile(sais.descr,0.05),paste0("R² = ",cor.sais),font=2)
  legend("topright",c("summer","autumn","winter","spring"),col=colo,pch=1,bty="n",ncol=2,cex=0.9)
  
  # Desais
  plot(desais.dP,desais.descr,col="darkgray",xlab="Pressure gradient (m)",ylab=paste0(descriptor," ",dist),main="Seasonal adjusted")
  abline(lm.desais$coefficients[1],lm.desais$coefficients[2],lwd=2)
  text(quantile(desais.dP,0.01),quantile(desais.descr,0.001),paste0("R² = ",cor.desais),font=2)
  
  graphics.off()
  
  # noise
  #noise <- apply(geo,3,function(x) sd(as.vector(x)))
  #noise <- rollapply(noise,nbdays,mean)
  #noise <- noise/dP
  #cor.no <- round(cor(noise,descr)^2,2)
  #lm.no <- lm(descr~noise)
  #sais.noise <- aggregate(noise,by=list(substr(dates[-c(length(dates)-1,length(dates))],6,10)),mean)[,2]
  
  # Graphique noise
  #png(filename = paste0(get.dirstr(k,rean),"plot.dP.noise/plot_noise_",descriptor,"_",dist,"_member",member,
  #                      "_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width = 800,height = 400,units = "px")
  #par(pty="s")
  #par(mfrow=c(1,2))
  #plot(noise,descr,type="n",xlab="Noise (-)",ylab=paste0(descriptor," ",dist))
  #points(noise[summ],descr[summ],pch=19,col="red",cex=0.2)
  #points(noise[aut],descr[aut],pch=19,col="darkorange",cex=0.2)
  #points(noise[win],descr[win],pch=19,col="blue",cex=0.2)
  #points(noise[spr],descr[spr],pch=19,col="olivedrab3",cex=0.2)
  #abline(lm.no$coefficients[1],lm.no$coefficients[2],lwd=1.5)
  #title(paste0("R² = ",cor.no))
  #legend(ifelse(dist=="RMSE","bottomright","topright"),c("summer","autumn","winter","spring"),
  #       col=c("red","darkorange","blue","olivedrab3"),pch=19,bty="n")
  #
  #plot(sais.noise,sais.descr,pch=19,cex=0.5,
  #     xlab="Interannual Noise (-)",ylab=paste0("Interannual ",descriptor," ",dist))
  #abline(lm(sais.descr~sais.noise),lwd=1.5)
  #title(paste0("R² = ",round(cor(sais.noise,sais.descr)^2,2)))
  #
  #graphics.off()
  
  # Graphique dP - noise
  #cor.b <- round(cor(dP,noise)^2,2)
  #lm.b <- lm(noise~dP)
  #
  #png(filename = paste0(get.dirstr(k,rean),"plot.dP.noise/plot_dP_noise_member",member,
  #                      "_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width = 400,height = 400,units = "px")
  #par(pty="s")
  #plot(dP,noise,type="n",xlab="Max pressure gradient (m)",ylab="Noise (-)")
  #points(dP[summ],noise[summ],pch=19,col="red",cex=0.2)
  #points(dP[aut],noise[aut],pch=19,col="darkorange",cex=0.2)
  #points(dP[win],noise[win],pch=19,col="blue",cex=0.2)
  #points(dP[spr],noise[spr],pch=19,col="olivedrab3",cex=0.2)
  #abline(lm.b$coefficients[1],lm.b$coefficients[2],lwd=1.5)
  #title(paste0("R² = ",cor.b))
  #graphics.off()
}

# Plan des indicateurs colorie par densite, gradient de pression, et saison
plot.empir.dP.sais <- function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F),coin=F,save=T,let=F){
  
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
  if(save){
    pdf(file=paste0(path1,"plot.empir.dP.sais",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,".pdf"),width=8,height=3.5) # manip substr pour enlever le std
    par(mfrow=c(1,3))
  }
  
  gamme <- c(0,4445)
  param <- param[,c("nbnei")]
  
  # Densite
  plot(descr1,descr2,
       col=getcol(param,range = gamme),
       xlab=paste0(namdescr[1]," ",dist[1]),
       ylab=paste0(namdescr[2]," ",dist[2]),
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3))
  addscale(vec = c(param,gamme))
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T)),
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.15,
       paste0(round(min(param,na.rm=T),2),"-",round(max(param,na.rm=T),2)))
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white") 
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  if(length(let)!=1) mtext(let[1],side=3,adj=0,line=0.8,font=1)
  
  # Gradient de pression
  deltaP <- get.dP(k[1],nbdays,start,end,rean[1])
  
  plot(descr1,descr2,
       col=getcol(deltaP),
       xlab=paste0(namdescr[1]," ",dist[1]),
       ylab=paste0(namdescr[2]," ",dist[2]),
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3))
  addscale(vec = round(deltaP,0))
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white")
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  if(length(let)!=1) mtext(let[2],side=3,adj=0,line=0.8,font=1)
  
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
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),  
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3))
  legend("topleft",c("summer","autumn","winter","spring"),col=colo,pch=1,bty="n",ncol=2)
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white")
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  if(length(let)!=1) mtext(let[3],side=3,adj=0,line=0.8,font=1)
  
  if(save) graphics.off()
  
}

# Saisonnalite de dP et d'un indicateur TWS et RMSE
plot.sais.all <- function(descr,k,nbdays=3,start="1950-01-01",end="2011-12-31",rean){
  
  dates <- getdates(start,end)
  
  # Calcul dP
  dP <- get.dP(k,nbdays,start,end,rean)
  
  # Import descr
  descr.tws <- get.descriptor(descr,k,"TWS",nbdays,start,end,F,rean)
  descr.rmse <- get.descriptor(descr,k,"RMSE",nbdays,start,end,F,rean)
  
  # dates
  dates.all <- dates[-((length(dates)-nbdays+2):length(dates))]
  dates.chro <- getdates(start = "1950-01-01",end = "1955-01-01")
  ind.chro <- match(dates.chro,dates)
  
  # Figure
  png(file = paste0("2_Travail/",rean,"/Rresults/overall/k",k,"/plot.sais.all/plot_sais_all_",descr,".png"),width = 8,height = 6,units = "in",res = 1200)
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
  plot.chronique(dates = dates.chro,vec = dP[ind.chro],label = "Pressure Gradient (m)",liss = F,interan = T,extr = T)
  mtext("d)",side=3,adj=0,line=0.8,font=1)
  plot.chronique(dates = dates.chro,vec = descr.tws[ind.chro],label = paste0(descr," TWS"),liss = F,interan = T,extr = T)
  mtext("e)",side=3,adj=0,line=0.8,font=1)
  plot.chronique(dates = dates.chro,vec = descr.rmse[ind.chro],label = paste0(descr," RMSE"),liss = F,interan = T,extr = T)
  mtext("f)",side=3,adj=0,line=0.8,font=1)
  graphics.off()
  
}



