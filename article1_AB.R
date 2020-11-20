source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Fonction utile
#smoothScatter(descr1,descr2,colramp = colorRampPalette(c("white","blue","darkgreen","gold","darkorange","red")),nbin = 100,bandwidth = 0.5)

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
combine.functions <- function(fun,descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean,quant=F){
  
  # plot.desais.descr
  if(fun==1){
    png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.desais.descr/plot_desais_",descr[1],"_",descr[2],"_",
                          dist,"_member",member,"_k",k,"mean",nbdays,"day_",start,"_",end,".png"),width = 550,height = 250,units = "px")
    plot.desais.1 <- plot.desais.descr(descr[1],k,dist,nbdays,start,end,rean)
    plot.desais.2 <- plot.desais.descr(descr[2],k,dist,nbdays,start,end,rean)
    grid.arrange(plot.desais.1,plot.desais.2,nrow=1,ncol=2)
    graphics.off()
  }
  
  # plot.empir.dp.sais
  if(fun==2){
      png(file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.empir.dP.sais-CV/plot_combine",ifelse(quant,"_quant",""),".png"),width=8,height=9,units = "in",res=1200)
      layout(matrix(c(rep(1,3),2:4,rep(5,3),6:8,rep(9,3),10:12),nrow = 6,ncol = 3,byrow = T),widths = rep(1,3),heights = c(0.1,1.2,0.1,1.2,0.1,1.2))
      let <- list(c("a)","b)","c)"),c("d)","e)","f)"),c("g)","h)","i)"))

    par(pty="s")
    nam <- c("Celerity","Singularity","Relative singularity")
    
    for(i in 1:length(descr)){
      # Nom de l'indicateur
      par(mar=c(0,0,0,0),pty="m")
      plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
      text(1,0.5,nam[i],cex=2,font=2)
      
      # Graphiques
      par(mar=c(3,4,0,3),pty="s")
      plot.empir.dP.sais(c(rean,rean),c(k,k),c(descr[i],descr[i]),dist,nbdays,start,end,coin=ifelse(descr[i]=="celnei",F,F),save=F,let=let[[i]],quant = quant)
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
    ggsave(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/plot.quant.descr/plot_combine.png"),plot = final,width = 10,height = 6)
    graphics.off()
  }
}

# Combinaison des boxplots par bv pour article
combine.quant.bv <- function(){
  
  # Parametres
  bv1 <- "Isere-seul"
  bv2 <- "Drac-seul"
  k <- 1
  dist <- "TWS"
  nbdays <- 3
  start <- "1950-01-01"
  end <- "2017-12-31"
  rean <- "ERA5"
  comm <- F
  spazm <- T
  supseuil <- T
  
  # Graphiques
  p <- compare.descr.bv(bv1 = bv1,bv2 = bv2,k = k,dist = dist,start = start,end = end,rean = rean,comm = F,nbdays = nbdays,spazm = spazm,save=F,legend=F,supseuil = supseuil)
  q <- compare.density.bv.all(bv = c(bv1,bv2),k = k,dist = dist,nbdays = nbdays,start = start,end = end,rean = rean,quant = F,save = F,legend = T,spazm = spazm,supseuil = supseuil)
  ggarrange(p,q,labels = c("a)","b)"),ncol = 2,nrow = 1)
  ggsave(filename = paste0(get.dirstr(k,rean),"combine.quant.bv/plot_",bv1,"_",bv2,"_",nbdays,"day_",start,"_",end,ifelse(comm,"_commun",""),ifelse(spazm,"_spazm",""),ifelse(supseuil,"_supseuil",""),".png"),width = 30,height = 11,units="cm",dpi = 200)
  graphics.off()
}

# Combinaison des boxplots par wp pour article
combine.quant.wp <- function(){
  
  # Parametres
  wp1 <- 1
  wp2 <- 2
  agreg <- T
  k <- 1
  dist <- "TWS"
  nbdays <- 3
  start <- "1950-01-01"
  end <- "2017-12-31"
  rean <- "ERA5"
  spazm <- T
  supseuil <- T
  
  # Graphiques
  p <- compare.descr.flow(flow = c(wp1,wp2),agreg = agreg,k = k,dist = dist,nbdays = nbdays,start = start,end = end,rean = rean,spazm = spazm,save=F,legend=F,supseuil=supseuil)
  q <- compare.density.wp.all(wp = c(wp1,wp2),agreg = agreg,k = k,dist = dist,nbdays = nbdays,start = start,end = end,rean = rean,quant = F,save = F,legend = T,spazm=spazm,supseuil=supseuil)
  ggarrange(p,q,labels = c("a)","b)"),ncol = 2,nrow = 1)
  ggsave(filename = paste0(get.dirstr(k,rean),"combine.quant.wp/plot_wp",wp1,"_wp",wp2,"_",nbdays,"day_",start,"_",end,ifelse(spazm,"_spazm",""),ifelse(supseuil,"_supseuil",""),".png"),width = 32,height = 11,units="cm",dpi = 200)
  graphics.off()
}

# Comparaison de la densite des max annuels de precip par rapport a la densite de tous les points, pour deux BVs
compare.density.bv <- function(descr=c("celnei","dP"),bv=c("Isere-seul","Drac-seul"),k,dist,nbdays,start,end,rean,quant=F,save=F){
  
  # Import des densites de pts pour un couple d'indicateur
  load(paste0(get.dirstr(k,rean),"compute.density/nbnei_",ifelse(quant,"quant_",""),descr[1],"_",descr[2],"_",
                ifelse(length(descr)==3,paste0(descr[3],"_"),""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  
  # Max annuels de precip
  ind.bv1 <- get.ind.max(type = "year",nbdays,start,end,bv = bv[1],spazm = F)
  ind.bv2 <- get.ind.max(type = "year",nbdays,start,end,bv = bv[2],spazm = F)
  
  # Mise en forme
  data <- data.frame(Ind = c(rep("all",length(nb)),rep(nam2str(bv[1]),length(ind.bv1)),rep(nam2str(bv[2]),length(ind.bv2))),
                     Nbnei = c(nb,nb[ind.bv1],nb[ind.bv2]))
  data$Ind <- factor(data$Ind,levels = c(nam2str(bv[1]),nam2str(bv[2]),"all"))
  
  # boxplot
  ggplot(data, aes(x=Ind, y=Nbnei, fill=Ind)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "none",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=1,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1","orange"))+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")+
    xlab("BV")+
    ylab("nbnei")
  
  if(save) ggsave(filename = paste0(get.dirstr(k,rean),"compare.density.bv/compare_density_",descr[1],"_",descr[2],"_",
                           ifelse(length(descr)==3,paste0(descr[3],"_"),""),ifelse(quant,"quant_",""),bv[1],"_",bv[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png"))
  
}

# Comparaison de la densite de pts dans le plan pour les max annuels, pour les 3 combaisons de descripteurs en meme temps (donc nbnei en percentile)
compare.density.bv.all <- function(bv=c("Isere-seul","Drac-seul"),k,dist,nbdays,start,end,rean,quant=F,save=F,legend=T,spazm=F,supseuil=F){
  
  descr <- list(
    c("celnei","dP"),
    c("singnei","dP"),
    c("rsingnei","dP")
  )
  
  # Import des densites de pts pour un couple d'indicateur
  nbnei <- as.data.frame(matrix(NA,length(getdates(start,end))-nbdays+1,length(descr)))
  for(i in 1:length(descr)){
    load(paste0(get.dirstr(k,rean),"compute.density/nbnei_",ifelse(quant,"quant_",""),descr[[i]][1],"_",descr[[i]][2],"_",
                ifelse(length(descr[[i]])==3,paste0(descr[[i]][3],"_"),""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    nbnei[,i] <- ecdf(nb)(nb)*100
    colnames(nbnei)[i] <- paste0(descr[[i]][1],"-MPD")
  }
  
  # Max annuels de precip et mise en forme
  if(!supseuil){
    ind.bv1 <- get.ind.max(type = "year",nbdays,start,end,bv = bv[1],spazm = spazm)
  }else{ind.bv1 <- get.ind.extr(bv[1],nbdays,start,end,nei=T,spazm,seuil="qua") }
  nbnei.1 <- cbind(rep(nam2str(bv[1]),length(ind.bv1)),nbnei[ind.bv1,])
  colnames(nbnei.1)[1] <- "Catchment"
  
  if(!supseuil){
    ind.bv2 <- get.ind.max(type = "year",nbdays,start,end,bv = bv[2],spazm = spazm)
  }else{ind.bv2 <- get.ind.extr(bv[2],nbdays,start,end,nei=T,spazm,seuil="qua") }
  nbnei.2 <- cbind(rep(nam2str(bv[2]),length(ind.bv2)),nbnei[ind.bv2,])
  colnames(nbnei.2)[1] <- "Catchment"
  
  data <- rbind(nbnei.1,nbnei.2)
  data <- pivot_longer(data,2:4,names_to = "descr",values_to = "percentile")
  data$descr <- factor(data$descr,levels=colnames(nbnei))
  
  # Boxplot
  p <- ggplot(data, aes(x=descr, y=percentile, fill=Catchment)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=1,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = ifelse(legend,"right","none"),legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=1,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=c(1.5,2.5), linetype="dashed")+
    xlab("Descriptor plane")+
    ylab("Percentile of number of neighbors (%)")
  
  if(save) {
    ggsave(filename = paste0(get.dirstr(k,rean),"compare.density.bv.all/compare_density_",ifelse(quant,"quant_",""),bv[1],"_",bv[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,ifelse(spazm,"_spazm",""),ifelse(supseuil,"_supseuil",""),".png"),plot=p,width = 8,height = 4)
    graphics.off()
  }else{p}
}

# Comparaison de la densite des max annuels de precip par rapport a la densite des points du meme WP
compare.density.wp <- function(descr=c("celnei","dP"),wp=c(1,2),agreg=T,k,dist,nbdays,start,end,rean,quant=F){
  
  # Import des densites de pts pour un couple d'indicateur
  nbnei <- list()
  for(i in 1:length(wp)){
  load(paste0(get.dirstr(k,rean),"compute.density/nbnei_",ifelse(quant,"quant_",""),"wp",wp[i],"_",descr[1],"_",descr[2],"_",
              ifelse(length(descr)==3,paste0(descr[3],"_"),""),ifelse(agreg,"agreg_",""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nbnei[[i]] <- nb
  }
  
  # Max annuels de precip
  ind.wp1 <- get.ind.max.flow(flow = wp[1],agreg,nbdays,start,end)
  ind.wp2 <- get.ind.max.flow(flow = wp[2],agreg,nbdays,start,end)
  tt <- get.wp(nbdays,start,end,agreg=agreg)
  
  # Mise en forme
  data <- data.frame(Ind = c(rep(paste0("All wp",wp[1]),length(nbnei[[1]])),rep(paste0("wp",wp[1]),length(ind.wp1)),
                             rep(paste0("All wp",wp[2]),length(nbnei[[2]])),rep(paste0("wp",wp[2]),length(ind.wp2))),
                     Nbnei = c(nbnei[[1]],nbnei[[1]][match(ind.wp1,which(tt==wp[1]))],
                               nbnei[[2]],nbnei[[2]][match(ind.wp2,which(tt==wp[2]))]))
  data$Ind <- factor(data$Ind,levels = c(paste0("wp",wp[1]),paste0("All wp",wp[1]),paste0("wp",wp[2]),paste0("All wp",wp[2])))
  
  # boxplot
  ggplot(data, aes(x=Ind, y=Nbnei, fill=Ind)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "none",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=1,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c(rep("cornflowerblue",2),rep("burlywood1",2)))+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")+
    xlab("BV")+
    ylab("nbnei")
  
  ggsave(filename = paste0(get.dirstr(k,rean),"compare.density.wp/compare_density_",descr[1],"_",descr[2],"_",
                           ifelse(length(descr)==3,paste0(descr[3],"_"),""),ifelse(quant,"quant_",""),"wp",wp[1],"_wp",wp[2],"_",ifelse(agreg,"agreg_",""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png"))
  
}

# Comparaison des densités pour les max annuels en fonction de leur flow, mais pour tous 3 couples de descripteurs en meme temps
compare.density.wp.all <- function(wp=c(1,2),agreg=T,k,dist,nbdays,start,end,rean,quant=F,save=T,legend=T,spazm=F,supseuil=F){
  
  descr <- list(
    c("celnei","dP"),
    c("singnei","dP"),
    c("rsingnei","dP")
  )
  
  if(wp[1]==1 & wp[2]==2 & agreg){
    namwp <- c("Zonal","Meridional")
  }else{namwp <- wp}
  
  # Imports
  tt <- get.wp(nbdays,start,end,agreg=agreg,bv="Isere",spazm=spazm)
  tt1 <- which(tt==wp[1])
  tt2 <- which(tt==wp[2])
  
  ind.wp1 <- get.ind.max.flow(flow = wp[1],agreg,nbdays,start,end,spazm,supseuil)
  ind.wp2 <- get.ind.max.flow(flow = wp[2],agreg,nbdays,start,end,spazm,supseuil)
  
  # Traitement separe des tt car longueur differente
  nbnei.1 <- as.data.frame(matrix(NA,length(tt1),length(descr)))
  for(i in 1:length(descr)){
    load(paste0(get.dirstr(k,rean),"compute.density/nbnei_",ifelse(quant,"quant_",""),"wp",wp[1],"_",descr[[i]][1],"_",descr[[i]][2],"_",
                ifelse(length(descr[[i]])==3,paste0(descr[[i]][3],"_"),""),ifelse(agreg,"agreg_",""),dist,ifelse(spazm,"_spazm",""),"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    nbnei.1[,i] <- ecdf(nb)(nb)*100
    colnames(nbnei.1)[i] <- paste0(descr[[i]][1],"-MPD")
  }
  
  nbnei.1 <- cbind(rep(namwp[1],length(ind.wp1)),nbnei.1[match(ind.wp1,tt1),])
  colnames(nbnei.1)[1] <- "Flow"
  
  nbnei.2 <- as.data.frame(matrix(NA,length(tt2),length(descr)))
  for(i in 1:length(descr)){
    load(paste0(get.dirstr(k,rean),"compute.density/nbnei_",ifelse(quant,"quant_",""),"wp",wp[2],"_",descr[[i]][1],"_",descr[[i]][2],"_",
                ifelse(length(descr[[i]])==3,paste0(descr[[i]][3],"_"),""),ifelse(agreg,"agreg_",""),dist,ifelse(spazm,"_spazm",""),"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    nbnei.2[,i] <- ecdf(nb)(nb)*100
    colnames(nbnei.2)[i] <- paste0(descr[[i]][1],"-MPD")
  }
  
  nbnei.2 <- cbind(rep(namwp[2],length(ind.wp2)),nbnei.2[match(ind.wp2,tt2),])
  colnames(nbnei.2)[1] <- "Flow"
  
  data <- rbind(nbnei.1,nbnei.2)
  data <- pivot_longer(data,2:4,names_to = "descr",values_to = "percentile")
  data$descr <- factor(data$descr,levels=colnames(nbnei.1)[-1])
  
  # boxplot
  p <- ggplot(data, aes(x=descr, y=percentile, fill=Flow)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=1,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = ifelse(legend,"right","none"),legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=0.5,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=c(1.5,2.5), linetype="dashed")+
    xlab("Descriptor plane")+
    ylab("Percentile of number of neighbors (%)")
  
  if(save) {
    ggsave(filename = paste0(get.dirstr(k,rean),"compare.density.wp.all/compare_density_",ifelse(quant,"quant_",""),wp[1],"_",wp[2],"_",dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,ifelse(spazm,"_spazm",""),ifelse(supseuil,"_supseuil",""),".png"),plot=p)
    graphics.off()
  }else{p}
}

# Comparaison des percentiles d'indicateurs pour deux bassins versants
compare.descr.bv <- function(bv1,bv2,descr=c("celnei","singnei","rsingnei","dP"),k,dist,nbdays,start,end,rean,comm=F,flow=F,save=T,legend=T,spazm=F,supseuil=F){
  
  # Import indicateurs
  mat <- NULL
  nam.descr <- NULL
  
  for(i in 1:length(descr)){
    mat <- cbind(mat,get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,rean=rean))
    nam.descr[i] <- ifelse(descr[i]=="dP","MPD",descr[i])
    colnames(mat)[ncol(mat)] <- nam.descr[i]
  }
  
  if(flow!=F){
    namflow <- ifelse(flow==1,"zonal","meridional")
    wp <- get.wp(nbdays,start,end,risk = F,bv = "Isere",agreg = T)
    mat <- apply(mat,2,function(v) {v[wp!=flow]=NA;return(v)})
  }
  
  #pwat <- get.pwat(k,nbdays,start,end,rean)
  #mat <- cbind(mat,pwat)
  mat <- as.data.frame(apply(mat,2,function(v) ecdf(v)(v)*100))
  
  # Extremes
  if(!supseuil){
  ind.extr1 <- get.ind.max(type="year",nbdays,start,end,bv1,spazm=spazm)
  ind.extr2 <- get.ind.max(type="year",nbdays,start,end,bv2,spazm=spazm)
  }else{
    ind.extr1 <- get.ind.extr(bv = bv1,nbdays,start,end,nei = T,spazm,seuil = "qua")
    ind.extr2 <- get.ind.extr(bv = bv2,nbdays,start,end,nei = T,spazm,seuil = "qua")
  }
  
  if(comm){
    commun <- which(ind.extr1==ind.extr2)
    print(paste0("Max annuels en commun: ",length(commun)))
    ind.extr1 <- ind.extr1[-commun]
    ind.extr2 <- ind.extr2[-commun]
  }
  
  bv <- rep(nam2str(bv1),length(ind.extr1))
  mat1 <- as.data.frame(cbind(bv,mat[ind.extr1,]))
  bv <- rep(nam2str(bv2),length(ind.extr2))
  mat2 <- as.data.frame(cbind(bv,mat[ind.extr2,]))
  mat <- rbind(mat1,mat2)
  
  mat <- pivot_longer(as.data.frame(mat),2:5,names_to = "descr",values_to = "percentile")
  mat$descr <- factor(mat$descr,levels = nam.descr)#c(descr,"pwat"))
  
  # Graphique
  p <- ggplot(mat, aes(x=descr, y=percentile, fill=bv)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=1,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = ifelse(legend,"right","none"),legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=1,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")+
    xlab("Atmospheric descriptor")+
    ylab("Percentile of descriptor value (%)")+
    labs(fill="Catchment")#,title = namdescr
  
  if(save) {
    ggsave(filename = paste0(get.dirstr(k,rean),"compare.descr.bv/plot_",bv1,"_",bv2,"_",nbdays,"day_",start,"_",end,ifelse(comm,"_commun",""),ifelse(flow!=F,paste0("_",namflow),""),ifelse(spazm,"_spazm",""),ifelse(supseuil,"_supseuil",""),".png"),plot=p,width = 20,height = 11,units="cm",dpi = 200)
    graphics.off()
  }else{p}
}

# Comparaison des percentiles d'indicateurs pour deux flux
compare.descr.flow <- function(flow=c(1,2),agreg=T,descr=c("celnei","singnei","rsingnei","dP"),k,dist,nbdays,start,end,rean,spazm=F,save=T,legend=T,supseuil=F){
  
  # Import indicateurs
  mat <- NULL
  nam.descr <- NULL
  
  for(i in 1:length(descr)){
    mat <- cbind(mat,get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,rean=rean))
    nam.descr[i] <- ifelse(descr[i]=="dP","MPD",nam2str(descr[i]))
    colnames(mat)[ncol(mat)] <- nam.descr[i]
  }
  
  # Traitement
  namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  namflow <- namflow[flow]
  ind.1 <- get.ind.max.flow(flow[1],agreg=agreg,nbdays,start,end,spazm,supseuil)
  ind.2 <- get.ind.max.flow(flow[2],agreg=agreg,nbdays,start,end,spazm,supseuil)
  
  wp <- get.wp(nbdays,start,end,agreg=T)
  mat.1 <- mat
  mat.1[wp!=flow[1],] <- NA
  mat.1 <- apply(mat.1,2,function(v) {ecdf(v)(v)*100})
  mat.1 <- mat.1[ind.1,]
  
  mat.2 <- mat
  mat.2[wp!=flow[2],] <- NA
  mat.2 <- apply(mat.2,2,function(v) {ecdf(v)(v)*100})
  mat.2 <- mat.2[ind.2,]
  
  mat <- as.data.frame(rbind(mat.1,mat.2))
  mat <- cbind(c(rep(namflow[1],length(ind.1)),rep(namflow[2],length(ind.2))),mat)
  colnames(mat)[1] <- "Flow"
  mat <- pivot_longer(mat,2:5,names_to = "descr",values_to = "percentile")
  mat$Flow <- factor(mat$Flow,levels = namflow[flow])
  mat$descr <- factor(mat$descr,levels = nam.descr)
  
  # Graphique
  p <- ggplot(mat, aes(x=descr, y=percentile, fill=Flow)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=1,size = 12,face = "bold"),axis.text.x = element_text(size=12),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = ifelse(legend,"right","none"),legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=0.4,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")+
    xlab("Atmospheric descriptor")+
    ylab("Percentile of descriptor value (%)")+
    labs(fill="Flow")#,title = namdescr
  
  if(save){
    ggsave(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compare.descr.flow/plot_",namflow[1],"_",namflow[2],"_",nbdays,"day_",start,"_",end,ifelse(spazm,"_spazm",""),ifelse(supseuil,"_supseuil",""),".png"),plot=p,width = 20,height = 11,units="cm",dpi = 200)
    graphics.off()
  }else{p}
}

# Calcul densite de points dans un plan
compute.density <- function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",quant=F,wp=F,agreg=F,spazm=F){
  
  # Import des descripteurs
  N <- length(getdates(start,end))-nbdays+1
  descr <- matrix(NA,N,length(descriptors))
  
  for(i in 1:length(descriptors)){
    descr[,i] <- get.descriptor(descriptors[i],k,dist[i],nbdays,start,end,standardize=T,rean)
  }
  
  if(quant){
    descr <- apply(descr,2,function(v) ecdf(v)(v)*100)
  }
  
  rad <- 0.5*sd(descr[,1]) # on prend un demi ecart type de rayon. Si WP, on calcule la densite selon le meme rayon, pour rester dans le referenciel de densite de tout le nuage
  
  if(wp!=F) {tt <- get.wp(nbdays,start,end,risk=F,bv="Isere",agreg=agreg,spazm = spazm); descr <- descr[tt==wp,]}
  
  # Calcul du voisinage
  nb <- NULL
  for(i in 1:nrow(descr)){
    if (i %%50==0) print(i)
    count <- nn2(data = descr,query = t(descr[i,]),searchtype = "radius",radius = rad,k = nrow(descr)) # nombre de voisins dans le rayon (tmp$idef: indices ou les deux descr sont non NA)
    nb[i] <- rowSums(count$nn.idx>0)-1
  }
  print(paste0("min density: ",min(nb)))
  print(paste0("max density: ",max(nb)))
  save(nb,file = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/compute.density/nbnei_",ifelse(quant,"quant_",""),ifelse(wp!=F,paste0("wp",wp,"_"),""),descriptors[1],"_",descriptors[2],ifelse(length(descriptors)==3,paste0("_",descriptors[3]),""),"_",ifelse(agreg,"agreg_",""),ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),ifelse(spazm,"_spazm",""),"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
}

# Comparaison des indicateurs calcules en moyenne journaliere ou sur les champs instantanes a 18h (ERA20C)
compare.ERA20C.ERA20C18 <- function(descr=c("celnei","singnei","rsingnei","dP"),k,dist,nbdays,start="1950-01-01",end="2010-12-31"){
  
  # Parametres graphiques
  png(filename = paste0("2_Travail/0_Present/Rresults/compare.ERA20C.ERA20C18/scatterplot_k",k,"_",dist,"_mean",nbdays,"days_",start,"_",end,".png"),width = 7,height = 7,units = "in",res = 300)
  par(mfrow=c(2,2),pty="s",mar=c(4,2,3,0))
  
  # Boucle sur les indicateurs
  for(i in 1:length(descr)){
    
    # Import
    ind.era20c <- get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,
                                 rean = "ERA20C",period = "present")
    ind.era20c18 <- get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,
                                 rean = "ERA20C_18",period = "present")
    ran <- range(ind.era20c,ind.era20c18)
    
    # Scatterplot
    plot(ind.era20c,ind.era20c18,xlim=ran,ylim=ran,xlab="ERA20C daily",ylab="ERA20C 1800 UTC",main=nam2str(descr[i],whole = T),pch=19)
    abline(0,1,col="red",lwd=2)
    text(ran[2]-1/3*(ran[2]-ran[1]),ran[2]-2/3*(ran[2]-ran[1]),paste0("R2 = ",round(cor(ind.era20c,ind.era20c18),2)),font=2)
  }
  graphics.off()
}

# Compte le % de sequences d'un meme flux en dessous d'un certain quantile de deux descripteurs
count.flow.descr <- function(flow,agreg,qua,descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean){
  
  # Import et traitement
  des1 <- get.descriptor(descr[1],k,dist,nbdays,start,end,standardize = F,rean)
  des1 <- ecdf(des1)(des1)
  
  des2 <- get.descriptor(descr[2],k,dist,nbdays,start,end,standardize = F,rean)
  des2 <- ecdf(des2)(des2)
  
  wp <- get.wp(nbdays,start,end,agreg=agreg)
  des1 <- des1[wp==flow]
  des2 <- des2[wp==flow]
  
  # Sortie
  print(nam2str(descr[1],whole=T))
  print(paste0("En dessous: ",round(table(des1<qua)[2]/length(des1)*100,1),"%"))
  print(paste0("Au dessus: ",round(table(des1<qua)[1]/length(des1)*100,1),"%"))
  
  print(nam2str(descr[2],whole=T))
  print(paste0("En dessous: ",round(table(des2<qua)[2]/length(des2)*100,1),"%"))
  print(paste0("Au dessus: ",round(table(des2<qua)[1]/length(des2)*100,1),"%"))
  
  plot(ecdf(des1))
  lines(ecdf(des2),col="red")
}

# Genere le tableau d'altitude de la region par maille de 75m (BD Alti France)
get.bdalti.region <- function(){
  
  # Dalles necessaires
  dalles <- c("0900_6450","0900_6525","0900_6600")
  X <- NULL
  
  Y <- NULL
  
  Z <- NULL
  
  for(i in 1:length(dalles)){
  # Import des mailles d'intérêt et extraction du pixel/altitude
  dalle <- read.asciigrid(paste0("2_Travail/Data/Carto/BDALTIV2_2-0_75M_ASC_LAMB93-IGN69_FRANCE_2020-04-28.7z/BDALTIV2_2-0_75M_ASC_LAMB93-IGN69_FRANCE_2020-04-28/BDALTIV2/1_DONNEES_LIVRAISON_2020-05-00201/BDALTIV2_MNT_75M_ASC_LAMB93_IGN69_FRANCE/BDALTIV2_75M_FXX_",dalles[i],"_MNT_LAMB93_IGN69.asc"))
  proj4string(dalle) <- CRS("+init=epsg:2154")
  dalle <- spTransform(dalle,CRS("+init=epsg:27572"))
  xi <- unique(as.matrix(coordinates(dalle))[,1])
  yi <- unique(rev(as.matrix(coordinates(dalle))[,2]))
  alt <- dalle@data[,1]
  dim(alt) <- c(length(xi),length(yi))
  alt <- alt[,ncol(alt):1]
  
  X=xi
  Y=c(Y,yi)
  Z=cbind(Z,alt)
  }
  
  # Dalles necessaires
  dalles2 <- c("0975_6450","0975_6525","0975_6600")
  X2 <- NULL
  
  Y2 <- NULL
  
  Z2 <- NULL
  
  for(i in 1:length(dalles2)){
    # Import des mailles d'intérêt et extraction du pixel/altitude
    dalle <- read.asciigrid(paste0("2_Travail/Data/Carto/BDALTIV2_2-0_75M_ASC_LAMB93-IGN69_FRANCE_2020-04-28.7z/BDALTIV2_2-0_75M_ASC_LAMB93-IGN69_FRANCE_2020-04-28/BDALTIV2/1_DONNEES_LIVRAISON_2020-05-00201/BDALTIV2_MNT_75M_ASC_LAMB93_IGN69_FRANCE/BDALTIV2_75M_FXX_",dalles2[i],"_MNT_LAMB93_IGN69.asc"))
    xi <- unique(as.matrix(coordinates(dalle))[,1])
    yi <- unique(rev(as.matrix(coordinates(dalle))[,2]))
    alt <- dalle@data[,1]
    dim(alt) <- c(length(xi),length(yi))
    alt <- alt[,ncol(alt):1]
    
    X2=xi
    Y2=c(Y2,yi)
    Z2=cbind(Z2,alt)
  }
  
  Xfin <- c(X,X2)
  Yfin <- Y
  Zfin <- rbind(Z,Z2)
  
  bdalti <- list(Fx=Xfin,Fy=Yfin,Fz=Zfin)
  
  # conversion
  XY <- cbind(c(Xfin,rep(1049925,1000)),Yfin)
  
  crs <- CRS("+init=epsg:2154")
  p <- SpatialPoints(XY, proj4string=crs)
  g <- spTransform(p, CRS("+init=epsg:27572"))
  fin=coordinates(g)
  
  bdalti <- list(Fx=fin[1:2000,1],Fy=fin[,2],Fz=Zfin)
  
}

# Correlation en valeur et en rang des indicateurs
get.cor.descr <- function(k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean,period="present"){
  
  #Indicateurs
  descr <- list(
    c("celnei","dP"),
    c("singnei","dP"),
    c("rsingnei","dP")
  )
  
  for(i in 1:length(descr)){
    des1 <- get.descriptor(descr[[i]][1],k,dist,nbdays,start,end,standardize = F,rean,period=period)
    des2 <- get.descriptor(descr[[i]][2],k,dist,nbdays,start,end,standardize = F,rean,period=period)
    print(descr[[i]])
    print(paste0("Classique: ",round(cor(des1,des2,method="pearson"),2)))
    print(paste0("Rang: ",round(cor(des1,des2,method="spearman"),2)))
  }
  
}

# Calcul de dP
get.dP <- function(k,nbdays,start="1950-01-01",end="2011-12-31",rean){
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean) 
  des <- apply(geo,3,function(x) max(x)-min(x))
  des <- rollapply(des,nbdays,mean)
}

# Max annuels de precip issus de zonal ou meridional (Drac et Isere)
get.ind.max.flow <- function(flow,agreg,nbdays,start,end,spazm=F,supseuil=F){
  
  # Import des max annuels des deux BVs
  if(!supseuil){
  ind1 <- get.ind.max(type = "year",nbdays,start,end,bv="Isere-seul",spazm)
  ind2 <- get.ind.max(type = "year",nbdays,start,end,bv="Drac-seul",spazm)
  }else{
    ind1 <- get.ind.extr(bv = "Isere-seul",nbdays,start,end,nei = T,spazm,seuil = "qua")
    ind2 <- get.ind.extr(bv = "Drac-seul",nbdays,start,end,nei = T,spazm,seuil = "qua")
  }
  
  # WPs
  wp <- get.wp(nbdays,start,end,risk=F,bv="Isere",agreg=agreg,spazm = spazm)
  ind <- sort(c(ind1[wp[ind1]==flow],ind2[wp[ind2]==flow]))
  unique(ind) # pour ne compter qu'une seule fois les dates doublons
}

# Calcul de pwat moyen journalier
get.pwat <- function(k=1,nbdays,start="1950-01-01",end="2011-12-31",rean,small_win){
  
  # Import
  if(rean=="20CR") vari <- "pwat"
  if(rean=="ERA5") vari <- "tcw"
  pwat <- getdata(k = k,day0 = start,day1 = end,rean = rean,var = vari,small_win=small_win) 
  
  # Moyenne sur la fenetre
  if(small_win & rean=="20CR"){# si small_win et 20CR, tableau et non array car seulement 2 pts de grille
  des <- apply(pwat,2,function(x) mean(x))
  }else{des <- apply(pwat,3,function(x) mean(x))}
  
  # Moyenne glissante sur nbdays
  des <- rollapply(des,nbdays,mean) 
}

# Carte composite des anomalies par wp ou groupement de wp
map.composite.wp <- function(wp,k,start,end,rean,leg=T,win=T,let=F,iso=T,agreg=F,anomalies=T,wind=F){
  
  # Import wp
  tt <- get.wp(nbdays = 1,start,end,risk=F,bv = "Isere",agreg=agreg)
  ind <- which(tt == wp)
  
  # Import des donnees
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean,large_win = F,small_win = F,all = T,var = "hgt")
  
  # Import lon/lat et fenetre
  nc <- load.nc(rean,var="hgt")
  nc <- nc[[k]]
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)
  
  fen <- getinfo_window(k = k,rean=rean,var = "hgt")
  
  # Vent
  if(wind){
    # Import des donnees
    uwind <- getdata(k = k,day0 = start,day1 = end,rean = rean,large_win = F,small_win = F,all = T,var = "uwind")
    vwind <- getdata(k = k,day0 = start,day1 = end,rean = rean,large_win = F,small_win = F,all = T,var = "vwind")
    
    # Import lon/lat
    nc.wind <- load.nc(rean,var="uwind")
    lon.wind <- nc.wind$dim$lon$vals
    lat.wind <- nc.wind$dim$lat$vals
    
    # Selection d'1 pt sur N uniquement
    N <- 3
    if(rean=="ERA5") N <- 24
    
    uwind <- uwind[seq(1,length(lon.wind),N),seq(1,length(lat.wind),N),]
    vwind <- vwind[seq(1,length(lon.wind),N),seq(1,length(lat.wind),N),]
    lon.wind <- lon.wind[seq(1,length(lon.wind),N)]
    lat.wind <- lat.wind[seq(1,length(lat.wind),N)]
    
    # Max de vent pour reglage de la taille des fleches
    if(rean=="20CR"){
    maxi.all <- get.min.max.wind(k,start="1950-01-01",end="2011-12-31",rean)[2]
    }
    if(rean=="ERA5"){
      maxi.all <-100
    }
  }
  
  # Traitement: composites
  geo.wp <- geo[,,ind]
  comp <- apply(geo.wp,1:2,mean)
  
  if(wind){
  uwind.wp <- uwind[,,ind]
  vwind.wp <- vwind[,,ind]
  comp.uwind <- apply(uwind.wp,1:2,mean)
  comp.vwind <- apply(vwind.wp,1:2,mean)
  maxi.fen <- max(sqrt(comp.uwind^2+comp.vwind^2))
  }
  
  # Traitement: anomalies
  if(anomalies){
    comp_all <- apply(geo,1:2,mean)
    comp <- comp - comp_all
    
    if(wind){
    comp.uwind.all <- apply(uwind,1:2,mean)
    comp.vwind.all <- apply(vwind,1:2,mean)
    comp.uwind <- comp.uwind - comp.uwind.all
    comp.vwind <- comp.vwind - comp.vwind.all
    maxi.fen <- max(sqrt(comp.uwind^2+comp.vwind^2))
    }
    }
  
  # Parametres graphiques
  if(anomalies){
    breaks <- seq(-200,200,length.out = 12)
    N <- 11
    lab <- seq(-200,200,50)
    lev <- seq(-200,200,100)
  }else{
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
  }
  
  # Noms des circulations
  if(agreg){namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  }else{namflow <- c("Atlantic\nWave","Steady\nOceanic","Southwest\nCirculation","South\nCirculation",
                     "Northeast\nCirculation","East\nReturn","Central\nDepression","Anticyclonic")}
  namflow <- namflow[wp]
  
  # Carte
  par(pty="s")
  cex <- 1.5
    
    if(leg){
      image.plot(lon,lat,comp,xlim=c(-15,25),ylim=c(25,65),
                 col=rev(brewer.pal(n = N, name = "RdBu")),
                 xlab="",ylab="",main="",xaxt="n",yaxt="n",
                 legend.line=-2.3, cex.axis=cex, cex.lab=cex, cex.main=cex,
                 breaks = breaks,axis.args = list(at=lab,labels=as.character(lab),cex.axis=1.3))
    }else{
      image(lon,lat,comp,xlim=c(-15,25),ylim=c(25,65),
            col=rev(brewer.pal(n = N, name = "RdBu")),
            xlab="",ylab="",main="",xaxt="n",yaxt="n",
            cex.axis=cex, cex.lab=cex, cex.main=cex,
            breaks = breaks)
    }
    
    data(wrld_simpl)
    plot(wrld_simpl, add = TRUE)
    points(6,45,col="red",pch=19)
    if(win) rect(xleft = lon[fen[1,1]]-1,ybottom = lat[fen[2,1]]-1,xright = lon[fen[1,1]+fen[1,2]-1]+1,ytop = lat[fen[2,1]+fen[2,2]-1]+1,lwd=2)
    if(let!=F) mtext(let, side=3, at=-20,line = 0,cex=1.4)
    if(iso) contour(x=lon,y=lat,z=comp, levels=lev, drawlabels=F, lty=1, lwd=1, add=TRUE, col="black")
    
    # Vent 500hPa
    if(wind){
      coord <- as.matrix(expand.grid(lon.wind,lat.wind))
      wi <- as.matrix(cbind(as.vector(comp.uwind),as.vector(comp.vwind)))
      
      arrow.plot(a1 = coord,a2 = wi,length=0.05,xpd=F,arrow.ex = 0.4*maxi.fen/maxi.all,lwd=2)
    }
    
    shadowtext(5,61,namflow,font=2,cex=2,col="black",bg="white",r=0.3)
    
    box()
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
  png(file = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.corner/map_corner_celnei.png"),width = 8,height = 14,units = "in",res = 1200)
  layout(matrix(c(1:12,rep(13,3)),nrow = 5,ncol = 3,byrow = T),widths = rep(1,3),heights = c(rep(1,4),0.5))

  for(i in 1:4){
    for(j in 1:3){
      date <- as.character(as.Date(dates[i])+j-1)
      map.geo(date = date,rean = rean,k = k,nbdays = 1,save = F,win = T,let = ifelse(j==1,lettre[i],F),leg=F,iso = T)
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

# Cartes d'une gamme croissante d'un descripteur pour gif
map.descr.gif <- function(type="all",descr,k,dist,nbdays,start="1950-01-01",end="2011-12-31",rean,period="present"){
  
  # Import
  des <- get.descriptor(descr,k,dist,nbdays,start,end,standardize = F,rean,period=period)
  dates <- getdates(start,end)
  
  # Traitement
  tri <- sort(des,decreasing=ifelse(type=="max",T,F),index.return=T)$ix
  if(type=="all"){
  ran <- tri[seq(1,length(tri),length.out = 60)]
  }else{
    ran <- tri[1:60]
  }
  
  # Cartes
  for(i in 1:length(ran)){
    print(i)
    png(filename = paste0(get.dirstr(k,rean,period),"map.descr.gif/",descr,"_",type,"/",descr,"_",i,"_",type,"_k",k,"_",nbdays,"day.png"),
                       width = ifelse(nbdays==3,1050,350),height = 350,units = "px")
    layout(matrix(1:nbdays,1,nbdays))
    map.geo(date = dates[ran[i]],rean = rean,k = k,nbdays = nbdays,save=F,win = T,let = F,leg = F,iso = T)
    graphics.off()
  }
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
  png(file = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.extr.dry/map_extr_dry.png"),width = 6,height = 7,units = "in",res = 1200)
  layout(matrix(c(1:4,rep(5,2)),nrow = 3,ncol = 2,byrow = T),widths = rep(1,2),heights = c(rep(1,2),0.5))
  
  for(i in 1:4){
      map.geo(date = dates[i],rean = rean,k = k,nbdays = 1,save = F,win = T,let = lettre[i],leg=F,iso = T)
  }
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = leg,at =  seq(0, 1,length.out = length(leg)),
              vertical = F,xlim = c(0.7,1.3),ylim = c(0.3,0.6),cex=1.4)
  text(x = 1,y = 0.75,"Geopotential height (m)",cex=1.6)
  graphics.off()
  
}

# Carte des geopotentiels des sequences min et max d'un descripteur
map.min.max <- function(descr,k,dist,nbdays,start,end,rean,poster=F){
  
  # min et max
  des <- get.descriptor(descr,k,dist,nbdays,start,end,T,rean,F)
  ind.min <- which.min(des)
  ind.max <- which.max(des)
  
  # Parametres de la legende
  if(k==1){
    ran <- seq(4850,6100)
    leg <- as.character(seq(4900,6100,200))
    N=11
  } else{
    ran <- seq(-300,400)
    leg <- as.character(seq(-300,400,100))
    N=7
  }
  
  # Cartes
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.min.max/map_min_max_",descr,ifelse(!poster,"","_poster"),".png"),width = 6,height = 6,units = "in",res = 1200)
  layout(matrix(c(1:6,rep(7,3)),nrow = 3,ncol = 3,byrow = T),widths = rep(1,3),heights = c(rep(1,2),0.5))
  if(poster) par(mar=c(0.5,0.5,0.5,0.5))

  for(i in 1:2){
    ind <- ifelse(i==1,ind.min,ind.max)
    for(j in 1:3){
      date <- getdates()[ind+j-1]
      if(!poster){
        map.geo(date = date,rean = rean,k = k,nbdays = 1,save = F,win = T,let = F,leg=F,iso = T)
      }else{
        map.geo.condens(date = date,rean = rean,k = k,nbdays = 1,save = F,win = T,let = F,leg=F,iso = T)
      }
    }
  }
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = leg,at =  ecdf(ran)(leg),
              vertical = F,xlim = c(0.7,1.3),ylim = c(0.55,0.85),cex=1.4)
  text(x = 1,y = 0.95,"Geopotential height (m)",cex=1.6)
  graphics.off()
  
}

# Carte des geopotentiels des sequences min et max de plusieurs descripteur
map.min.max.all <- function(descr=c("celnei","singnei","rsingnei"),k,dist,nbdays,start,end,rean,wind=F){
  
  # Parametres graphiques
  nam <- nam2str(descr,whole = T)
  let <- c("a)","b)","c)","d)","e)","f)")
  
  if(k==1){
    ran <- seq(4850,6100)
    leg <- as.character(seq(4900,6100,200))
    N=11
  } else{
    ran <- seq(-300,400)
    leg <- as.character(seq(-300,400,100))
    N=7
  }
  
  # Graphiques
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.min.max.all/map_min_max_",paste(descr,collapse = "_"),".png"),width = 5.5,height = 12,units = "in",res = 200)
  #layout(matrix(c(rep(1,7),2:8,rep(9,7),10:16,rep(17,7),18:24,rep(25,7)),nrow = 7,ncol = 7,byrow = T),widths = c(rep(1,3),0.4,rep(1,3)),heights = c(0.2,1.2,0.2,1.2,0.2,1.2,0.5))
  layout(matrix(c(rep(1,3),2:7,rep(8,3),9:14,rep(15,3),16:21,rep(22,3)),nrow = 10,ncol = 3,byrow = T),widths = rep(1,3),heights = c(0.3,1.2,1.2,0.3,1.2,1.2,0.3,1.2,1.2,0.7))
  
  for(i in 1:length(descr)){
    
    # Nom de l'indicateur
    par(mar=c(0,0,0,0),pty="m")
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,nam[i],cex=2.5,font=2)
    
    # min et max
    des <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,
                          start = start,end = end,standardize = F,rean = rean,threeday = F)
    ind.min <- which.min(des)
    ind.max <- which.max(des)
    
    # Cartes
    par(mar=c(0.3,2,0.3,0))
    
    # min
    ind <- ind.min
      for(j in 1:3){
        date <- getdates(start,end)[ind+j-1]
        map.geo(date = date,rean = rean,k = k,nbdays = 1,save = F,win = T,let = F,leg=F,iso = F,wind = wind,condens = T)
        if(j==1) mtext(let[(i-1)*2+1], side=3, at=-25,line = -3,cex=1.2)
        }
    #plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    
    # max
    ind <- ind.max
    for(j in 1:3){
      date <- getdates(start,end)[ind+j-1]
      map.geo(date = date,rean = rean,k = k,nbdays = 1,save = F,win = T,let = F,leg=F,iso = F,wind = wind,condens = T)
      if(j==1) mtext(let[(i-1)*2+2], side=3, at=-25,line = -3,cex=1.2)
      }
  }
  
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = leg,at =  ecdf(ran)(leg),
              vertical = F,xlim = c(0.7,1.3),ylim = c(0.25,0.65),cex=1.4)
  text(x = 1,y = 0.85,"Geopotential height (m)",cex=1.6)
  
  graphics.off()
}

# Carte de l'eau precipitable d'un jour donne
map.pwat <- function(date,rean,k,nbdays=1,save=F,win=F,let=F,leg=T){
  
  # Import des donnees
  nc <- load.nc(rean,var = "pwat")
  fen <- getinfo_window(k,var = "pwat",rean = rean)
  
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  num0<-date_to_number(nc,date,rean)
  nc_close(nc)
  
  pwat <- ncvar_get(nc,varid="pwat",start=c(1,1,num0),count=c(length(lon),length(lat),nbdays))
  dim(pwat) <- c(length(lon),length(lat),nbdays)
  
  # Parametres graphiques
    breaks <- seq(0,57,length.out = 10)
    N <- 9
    lab <- seq(0,57,10)
  
  # Carte
  if(save) {
    png(filename = paste0("2_Travail/0_Present/20CR/Rresults/overall/k",k,"/map.pwat/",date,"_k",k,"_",nbdays,"day.png"),
        width = ifelse(nbdays==3,1050,350),height = 350,units = "px")
    layout(matrix(1:nbdays,1,nbdays))
  }
  
  par(pty="s")
  if(nbdays==3) par(mar=c(6,4,6,4))
  cex <- ifelse(nbdays==3,1.8,1.5)
  
  for(i in 1:nbdays){
    
    if(leg){
      image.plot(lon,lat,pwat[,,i],xlim=c(-20,25),ylim=c(25,70),
                 col=brewer.pal(n = N, name = "BuGn"),
                 xlab="Longitude (°)",ylab="Latitude (°)",main=as.Date(date)+i-1,
                 legend.line=-2.3, cex.axis=cex, cex.lab=cex, cex.main=cex,
                 breaks = breaks,axis.args = list(at=lab,labels=as.character(lab),cex.axis=1.3))
    }else{
      image(lon,lat,pwat[,,i],xlim=c(-20,25),ylim=c(25,70),
            col=brewer.pal(n = N, name = "BuGn"),
            xlab="Longitude (°)",ylab="Latitude (°)",main=as.Date(date)+i-1,
            cex.axis=cex, cex.lab=cex, cex.main=cex,
            breaks = breaks)
    }
    
    data(wrld_simpl)
    plot(wrld_simpl, add = TRUE)
    points(6,45,col="red",pch=19)
    if(win) rect(xleft = lon[fen[1,1]]-1,ybottom = lat[fen[2,1]]-1,xright = lon[fen[1,1]+fen[1,2]-1]+1,ytop = lat[fen[2,1]+fen[2,2]-1]+1,lwd=2)
    if(i==1 & let!=F) mtext(let, side=3, at=-30,line = 2,cex=1.5)
    box()
  }
  
  if(save) graphics.off()
  
}

# Carte composite de geopotentiel et d'anomalies pour deux WP
map.wp.flow <- function(flow=c(1,2),agreg=T,k,rean,start="1950-01-01",end="2011-12-31",wind=F){
  
  # Noms des circulations
  if(agreg){namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  }else{namflow <- c("Atlantic\nWave","Steady\nOceanic","Southwest\nCirculation","South\nCirculation",
             "Northeast\nCirculation","East\nReturn","Central\nDepression","Anticyclonic")}
  namflow <- namflow[flow]
  
  # Graphique
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/map.wp.flow/map_wp_flow_",namflow[1],"_",namflow[2],"_wp","_member",member,"_k",k,"_",start,"_",end,".png"),width = 6,height = 5,units = "in",res=200)
  layout(matrix(1:6,2,3,byrow=F),width=c(1.3,1.3,0.4))
  par(pty="s",mar=c(0,3,1,1))
  let <- list(c("a)","c)"),c("b)","d)"))
  let <- list(c(F,F),c(F,F))
  
  # Cartes composites
  for(i in 1:length(flow)){
    map.composite.wp(flow[i],k,start,end,rean,leg=F,let=let[[i]][1],agreg=agreg,anomalies = F,wind=wind,iso = T)
    map.composite.wp(flow[i],k,start,end,rean,leg=F,let=let[[i]][2],agreg=agreg,anomalies = T,wind=wind,iso = T)
  }
  
  # Legende cartes
  par(pty="m",mar=c(1,0,2,0))
  
  if(k==1){ 
    N <- 11
    leg <- seq(4900,6100,200)
    ran <- seq(4850,6100)
  }else{
    N <- 7
    leg <- seq(-300,400,100)
    ran <- seq(-300,400)
  }
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(0,1),ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = paste0("        ",leg),at = ecdf(ran)(leg),
              vertical = T,xlim = c(0.1,0.4),ylim = c(0,1),cex=1.4)
  
  N <- 11
  leg <- seq(-200,200,50)
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(0,1),ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = paste0("        ",leg),at =  seq(0, 1,length.out = length(leg)),
              vertical = T,xlim = c(0.1,0.4),ylim = c(0,1),cex=1.4)
  
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
  precip <- get.precip(nbdays,start,end)
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end,bv="Isere-seul")
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
    ylim(c(0,100))+
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
  
  # Extremes
  precip <- get.precip(nbdays,start,end)
  ind.extr <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end)
  
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
  points(dP[ind.extr],descr[ind.extr],pch=21,bg="white")
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
  points(desais.dP[ind.extr],desais.descr[ind.extr],pch=21,bg="white")
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

# Plan des indicateurs colorie par densite pour deux bv differents
plot.empir.bv <- function(bv1,bv2,rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F),spazm=F,save=T,let=F,quant=F,sea=F,pwat=F,comm=F,max=T,supseuil=F){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/0_Present/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/0_Present/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import descripteurs
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=F,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=F,rean[2],threeday[2])

  if(quant){
    descr1 <- ecdf(descr1)(descr1)*100
    descr2 <- ecdf(descr2)(descr2)*100
  }
  
  # Import nbnei
  if(TRUE %in% threeday) {comp <- paste0("_threeday",which(threeday==T))
  } else {comp <- ""}
  
  if(!quant){
    load(file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,comp,".Rdata"))
    nb<-param[,"nbnei"]
  }else{
    load(file = paste0("2_Travail/0_Present/",rean[1],"/Rresults/overall/k",k[1],"/compute.density/nbnei_",ifelse(quant,"quant_",""),descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  }
  
  # Nom propre des indicateurs
  cond <- substr(descriptors,nchar(descriptors)-1,nchar(descriptors)) == substr(radtype,nchar(radtype)-1,nchar(radtype))
  descriptors[cond] <- substr(descriptors[cond],1,nchar(descriptors[cond])-2)
  namdescr <- nam2str(descriptors, cloud = TRUE,whole = T)
  
  # Points aleatoires dans le plan
  ale <- sample(1:length(descr1),62)
  
  # Extremes
  if(!supseuil){
  ind.extr1 <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end,bv=bv1,spazm = spazm)
  ind.extr2 <- get.ind.max(type = "year",nbdays = nbdays,start = start,end = end,bv=bv2,spazm= spazm)
  }else{
    ind.extr1 <- get.ind.extr(bv = bv1,nbdays = nbdays,start = start,end = end,nei = T,spazm = spazm,seuil = "qua")
    ind.extr2 <- get.ind.extr(bv = bv2,nbdays = nbdays,start = start,end = end,nei = T,spazm = spazm,seuil = "qua")
  }
  
  if(comm){
  commun <- which(ind.extr1==ind.extr2)
  print(paste0("Annual max in common: ",length(commun)))
  ind.extr1 <- ind.extr1[-commun]
  ind.extr2 <- ind.extr2[-commun]
  }
  
  # Import pwat
  if(pwat){
    pw <- get.pwat(k = 1,nbdays,start,end,rean[1])
    pw <- ecdf(pw)(pw)*100
  }
  
  # Calcul saison
  if(sea){
    dates <- getdates(start,end)
    se <- dates
    se[which(substr(dates,6,7) %in% c("06","07","08"))] <- 1
    se[which(substr(dates,6,7) %in% c("09","10","11"))] <- 2
    se[which(substr(dates,6,7) %in% c("12","01","02"))] <- 3
    se[which(substr(dates,6,7) %in% c("03","04","05"))] <- 4
    se <- as.numeric(se)
    colo <- c("red","darkorange","blue","olivedrab3")
  }
  
  # Creation du plot
  if(save){
    png(filename=paste0(path1,"plot.empir.bv",get.CVstr(CV),"/plot_",bv1,"_",bv2,"_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,ifelse(quant,"_quant",""),ifelse(sea,"_sea",""),ifelse(pwat,"_pwat",""),ifelse(comm,"_comm",""),ifelse(max,"","_without_max"),ifelse(supseuil,"_supseuil",""),".png"),width=8,height=5,units = "in",res=1200) # manip substr pour enlever le std
    par(mfrow=c(1,2),pty="s")
  }
  
  # Parametres graphiques
  if(quant){gamme <- c(0,3756)
  }else{gamme <- c(0,4319)}
  
  if(sea){
    bg1 <- colo[se[ind.extr1]]
    bg2 <- colo[se[ind.extr2]]
  }else{bg1 <- bg2 <- "white"}
  
  if(pwat){
    cex1 <- pw[ind.extr1]/50
    cex2 <- pw[ind.extr2]/50
  }else{cex1 <- cex2 <- 1.2}
  
  # bv1
  plot(descr1,descr2,
       col=getcol(nb,range = gamme),
       xlab="",
       ylab="",
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3),
       main=nam2str(bv1),cex.axis=1.2,cex.main=1.5,xaxt="n",yaxt="n")
  title(xlab=paste0(namdescr[1],ifelse(quant," percentile ","")), line=2.5, cex.lab=1.2)
  title(ylab=paste0(namdescr[2],ifelse(quant," percentile ","")), line=2.5, cex.lab=1.2)
  if(quant){
  axis(1,seq(0,100,20),seq(0,100,20))
  axis(2,seq(0,100,20),seq(0,100,20))
  }else{axis(1);axis(2)}
  addscale(vec = c(nb,gamme),r=0)
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.9,
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2,
       paste0(round(min(nb,na.rm=T),2),"-",round(max(nb,na.rm=T),2)))
  #points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  if(max) points(descr1[ind.extr1],descr2[ind.extr1],pch=22,bg=bg1,cex=cex1) 
  if(length(let)!=1) mtext(let[1],side=3,adj=-0.2,line=-0.5,font=1,cex=1.2)
  
  # bv2
  plot(descr1,descr2,
       col=getcol(nb,range = gamme),
       xlab="",
       ylab="",
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3),
       main=nam2str(bv2),cex.axis=1.2,cex.main=1.5,xaxt="n",yaxt="n")
  title(xlab=paste0(namdescr[1],ifelse(quant," percentile ","")), line=2.5, cex.lab=1.2)
  title(ylab=paste0(namdescr[2],ifelse(quant," percentile ","")), line=2.5, cex.lab=1.2)
  if(quant){
    axis(1,seq(0,100,20),seq(0,100,20))
    axis(2,seq(0,100,20),seq(0,100,20))
  }else{axis(1);axis(2)}
  addscale(vec = c(nb,gamme),r=0)
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T))*0.9,
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.2,
       paste0(round(min(nb,na.rm=T),2),"-",round(max(nb,na.rm=T),2)))
  #points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  if(max) points(descr1[ind.extr2],descr2[ind.extr2],pch=22,bg=bg2,cex=cex2) 
  if(length(let)!=1) mtext(let[2],side=3,adj=-0.2,line=-0.5,font=1,cex=1.2)
  
  if(save) graphics.off()
  
}

# plot.empir.bv combine
plot.empir.bv.combine <- function(){
  
  # Parametres
  bv1 <- "Isere-seul"
  bv2 <- "Drac-seul"
  rean <- c("ERA5","ERA5")
  k <- c(1,1)
  descr <- list(
    c("celnei","dP"),
    c("singnei","dP"),
    c("rsingnei","dP")
  )
  nam <- c("Celerity vs. MPD","Singularity vs. MPD","Relative singularity vs. MPD")
  dist <- c("TWS","TWS")
  nbdays=3
  start="1950-01-01"
  end="2017-12-31"
  let <- list(c("a)","b)"),c("c)","d)"),c("e)","f)"))
  radtype="nrn05"
  CV=TRUE;threeday=c(F,F);spazm=T;save=F;quant=T;sea=F;pwat=F;comm=F;max=T;supseuil=T
  
  # Graphique
  png(file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.empir.bv-CV/plot_combine_",bv1,"_",bv2,ifelse(spazm,"_spazm",""),ifelse(quant,"_quant",""),ifelse(supseuil,"_supseuil",""),".png"),width=6,height=9,units = "in",res=200)
  layout(matrix(c(rep(1,2),2:3,rep(4,2),5:6,rep(7,2),8:9),nrow = 6,ncol = 2,byrow = T),widths = rep(1,2),heights = c(0.1,1.2,0.1,1.2,0.1,1.2))

  for(i in 1:length(descr)){
    # Nom de l'indicateur
    par(mar=c(0,0,0,0),pty="m")
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,nam[i],cex=2,font=2)
    
    # Graphiques
    par(mar=c(4,1.5,1.5,0),pty="s")
    plot.empir.bv(bv1 = bv1,bv2 = bv2,rean = rean,k = k,descriptors = descr[[i]],dist = dist,
                  nbdays = nbdays,start = start,end = end,radtype = radtype,
                  CV = CV,threeday = threeday,spazm = spazm,save = save,let = let[[i]],
                  quant = quant,sea = sea,pwat = pwat,comm = comm,max = max,supseuil = supseuil)
  }
  graphics.off()
}

# Plan des indicateurs colorie par densite, gradient de pression, et saison
plot.empir.dP.sais <- function(rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F),coin=F,save=T,let=F,quant=F){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/0_Present/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/0_Present/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import precip
  precip <- get.precip(nbdays,start,end)
  
  # Import descripteurs
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=F,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=F,rean[2],threeday[2])
  
  if(quant){
    descr1 <- ecdf(descr1)(descr1)*100
    descr2 <- ecdf(descr2)(descr2)*100
  }
  
  # Import nbnei
  if(TRUE %in% threeday) {comp <- paste0("_threeday",which(threeday==T))
  } else {comp <- ""}
  
  if(!quant){
  load(file=paste0(path1,"fit.empir",get.CVstr(CV),"/",path2,comp,".Rdata"))
  nb<-param[,"nbnei"]
  }else{
    load(file = paste0("2_Travail/0_Present/",rean[1],"/Rresults/overall/k",k[1],"/compute.density/nbnei_",ifelse(quant,"quant_",""),descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  }
  
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
  pdf(file=paste0(path1,"plot.empir.dP.sais",get.CVstr(CV),"/plot_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),comp,ifelse(quant,"_quant",""),".pdf"),width=8,height=3.5) # manip substr pour enlever le std
  par(mfrow=c(1,3))
  }
  
  gamme <- c(0,4094)
  
  # Densite
  plot(descr1,descr2,
       col=getcol(nb,range = gamme),
       xlab=paste0(ifelse(quant,"Percentile ",""),namdescr[1]," ",ifelse(descriptors[1]!="dP",dist[1],"")),
       ylab=paste0(ifelse(quant,"Percentile ",""),namdescr[2]," ",ifelse(descriptors[2]!="dP",dist[2],"")),
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3))
  addscale(vec = c(nb,gamme))
  text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T)),
       y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.15,
       paste0(round(min(nb,na.rm=T),2),"-",round(max(nb,na.rm=T),2)))
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white") 
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  if(length(let)!=1) mtext(let[1],side=3,adj=0,line=0.8,font=1)
  
  # Gradient de pression
  deltaP <- get.dP(k[1],nbdays,start,end,rean[1])
  
  plot(descr1,descr2,
       col=getcol(deltaP),
       xlab=paste0(ifelse(quant,"Percentile ",""),namdescr[1]," ",ifelse(descriptors[1]!="dP",dist[1],"")),
       ylab=paste0(ifelse(quant,"Percentile ",""),namdescr[2]," ",ifelse(descriptors[2]!="dP",dist[2],"")),
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3))
  addscale(vec = round(deltaP,0))
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
         xlab=paste0(ifelse(quant,"Percentile ",""),namdescr[1]," ",ifelse(descriptors[1]!="dP",dist[1],"")),
         ylab=paste0(ifelse(quant,"Percentile ",""),namdescr[2]," ",ifelse(descriptors[2]!="dP",dist[2],"")),
       xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),  
       ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3))
  legend("topleft",c("summer","autumn","winter","spring"),col=colo,pch=1,bty="n",ncol=2)
  points(descr1[ale],descr2[ale],pch=19,cex=0.9)
  points(descr1[ind.extr],descr2[ind.extr],pch=21,bg="white")
  if(coin) points(descr1[pos.coin],descr2[pos.coin],pch=19,cex=0.8,col="red")
  if(length(let)!=1) mtext(let[3],side=3,adj=0,line=0.8,font=1)
  
  if(save) graphics.off()
  
}

# Plan des indicateurs colorie par WP (si save, on met les 8 wp)
plot.empir.wp <- function(wp=1:8,rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F),spazm=F,save=T,let=F,quant=F,agreg=F,title="",supseuil=F){
  
  # Definition du repertoire de travail (lecture et ecriture)
  if(rean[1] != rean[2]){ 
    path1 <- paste0("2_Travail/0_Present/Comparaison_reanalyses/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",rean[1],"_",rean[2],"_",dist[1],"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  if(k[1] != k[2]){
    path1 <- paste0("2_Travail/0_Present/",rean[1],"/Rresults/overall/")
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_k",k[1],"_k",k[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  } 
  if(rean[1] == rean[2] & k[1] == k[2]){
    path1 <- get.dirstr(k[1],rean[1])
    path2 <- paste0(descriptors[1],"_",descriptors[2],"_",ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,get.stdstr(TRUE),"_",radtype)
  }
  
  # Import descripteurs
  descr1<-get.descriptor(descriptors[1],k[1],dist[1],nbdays,start,end,standardize=F,rean[1],threeday[1])
  descr2<-get.descriptor(descriptors[2],k[2],dist[2],nbdays,start,end,standardize=F,rean[2],threeday[2])
  
  if(quant){
    descr1 <- ecdf(descr1)(descr1)*100
    descr2 <- ecdf(descr2)(descr2)*100
  }
  
  # Import wp
  tt <- get.wp(nbdays,start,end,risk = F,bv = "Isere",agreg=agreg,spazm = spazm)
  
  # Nom propre des indicateurs
  cond <- substr(descriptors,nchar(descriptors)-1,nchar(descriptors)) == substr(radtype,nchar(radtype)-1,nchar(radtype))
  descriptors[cond] <- substr(descriptors[cond],1,nchar(descriptors[cond])-2)
  namdescr <- nam2str(descriptors, cloud = TRUE,whole = T) # whole permet d'avoir le nom complet de l'indicateur
  
  # Creation du plot
  if(save){
    png(file=paste0(path1,"plot.empir.wp",get.CVstr(CV),"/plot_wp_",substr(path2,1,nchar(path2)-10),substr(path2,nchar(path2)-5,nchar(path2)),ifelse(quant,"_quant",""),ifelse(supseuil,"_supseuil",""),".png"),width=6,height=12,units = "in",res=1200)
    par(mfrow=c(4,2),pty="s")
  }
  
  if(agreg){
    #namflow <- c("Zonal","Meridional","North-East","Anticyclonic")
    namflow <- c("Atlantic","Mediterranean","Continental","Anticyclonic")
  }else{namflow <- c("Atlantic\nWave","Steady\nOceanic","Southwest\nCirculation","South\nCirculation",
                     "Northeast\nCirculation","East\nReturn","Central\nDepression","Anticyclonic")}
  namflow <- namflow[wp]
  
  gamme <- list(rep(c(0,2965),length(wp)))
  if(agreg & quant) gamme <- list(c(0,3142),c(0,1505))
  if(agreg & !quant) gamme <- list(c(0,2386),c(0,1405))
  
  for(i in wp){
    # Import densite de pts
    load(file = paste0("2_Travail/0_Present/",rean[1],"/Rresults/overall/k",k[1],"/compute.density/nbnei_",ifelse(quant,"quant_",""),"wp",i,"_",descriptors[1],"_",descriptors[2],"_",ifelse(agreg,"agreg_",""),ifelse(dist[1]!=dist[2],paste0(dist[1],"_",dist[2]),dist[1]),ifelse(spazm,"_spazm",""),"_member",member,"_k",k[1],"_mean",nbdays,"day_",start,"_",end,".Rdata"))
    occ <- round(sum(tt==i)/length(tt)*100,0)
    cex <- 1.5
    
    # Import des max
    ind.max <- get.ind.max.flow(flow=i,agreg,nbdays,start,end,supseuil=supseuil,spazm = spazm)
    
    # plot
    plot(descr1,descr2,
         col="grey",
         xlab="",
         ylab="",
         xlim=c((min(descr1,na.rm=T)+max(descr1,na.rm=T))/2-((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2),(min(descr1,na.rm=T)+max(descr1,na.rm=T))/2+((max(descr1,na.rm=T)-min(descr1,na.rm=T))*1.3/2)),
         ylim=c(min(descr2,na.rm=T),min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.3),
         main=namflow[i], cex.main=cex,xaxt="n",yaxt="n")
    title(xlab=paste0(namdescr[1],ifelse(quant," percentile ","")), line=2.5, cex.lab=1.2)
    title(ylab=paste0(namdescr[2],ifelse(quant," percentile ","")), line=2.5, cex.lab=1.2)
    if(quant){
    axis(1,seq(0,100,20),seq(0,100,20))
    axis(2,seq(0,100,20),seq(0,100,20))
    }else{axis(1);axis(2)}
    points(descr1[tt==i],descr2[tt==i],col=getcol(c(nb,gamme[[i]])))
    points(descr1[ind.max],descr2[ind.max],pch=22,bg="white",cex=1.2)
    addscale(vec = c(nb,gamme[[i]]))
    text(x=min(descr1,na.rm=T)+(max(descr1,na.rm=T)-min(descr1,na.rm=T)),
         y=min(descr2,na.rm=T)+(max(descr2,na.rm=T)-min(descr2,na.rm=T))*1.15,
         paste0(round(min(nb,na.rm=T),2),"-",round(max(nb,na.rm=T),2)))
    if(length(let)!=1) mtext(let[i],side=3,adj=-0.2,line=-0.5,font=1,cex=1.2)
  }
  if(save) graphics.off()
  
}

# plot.empir.wp combine
plot.empir.wp.combine <- function(){
  
  # Parametres
  wp <- c(1,2)
  namflow <- c("Atlantic","Mediterranean")
  rean <- c("20CR","20CR")
  k <- c(1,1)
  descr <- list(
    c("celnei","dP"),
    c("singnei","dP"),
    c("rsingnei","dP")
  )
  nam <- c("Celerity vs. MPD","Singularity vs. MPD","Relative singularity vs. MPD")
  
  dist <- c("TWS","TWS")
  nbdays=3
  start="1950-01-01"
  end="2011-12-31"
  let <- list(c("a)","b)"),c("c)","d)"),c("e)","f)"))
  radtype="nrn05"
  agreg=T
  CV=TRUE;threeday=c(F,F);spazm=T;save=F;quant=T;sea=F;pwat=F;comm=F;supseuil=T
  
  # Graphique
  png(file=paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.empir.wp-CV/plot_combine_",namflow[1],"_",namflow[2],ifelse(spazm,"_spazm",""),ifelse(quant,"_quant",""),ifelse(supseuil,"_supseuil",""),".png"),width=6,height=9,units = "in",res=200)
  layout(matrix(c(rep(1,2),2:3,rep(4,2),5:6,rep(7,2),8:9),nrow = 6,ncol = 2,byrow = T),widths = rep(1,2),heights = c(0.1,1.2,0.1,1.2,0.1,1.2))
  
  for(i in 1:length(descr)){
    # Nom de l'indicateur
    par(mar=c(0,0,0,0),pty="m")
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,nam[i],cex=2,font=2)
    
    # Graphiques
    par(mar=c(4,1.5,1.5,0),pty="s")
    plot.empir.wp(wp,rean,k,descr[[i]],dist,nbdays,start,end,radtype,
                  CV,threeday,spazm,save,let[[i]],quant,agreg,supseuil=supseuil)
  }
  graphics.off()
}

# Combinaison des deux graphiques plot.period.ana() et plot.month.ana()
plot.fig.ana <- function(Q="05",k,nbdays,dist,start="1950-01-01",end="2011-12-31",rean){
  
  plot1 <- plot.month.ana(Q,k,nbdays,dist,start,end,rean)
  plot2 <- plot.period.ana(Q,k,nbdays,dist,start,end,rean)
  plot.final <- ggarrange(plot1, NULL, plot2, nrow=1,ncol=3,widths = c(1,0.05,1.5),labels = c("a)","","b)"))
  ggsave(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k1/plot.fig.ana/plot_fig_ana_mean",nbdays,"days.png"),plot = plot.final,width = 10,height = 4.5)
}
  
# Repartition des mois des analogues par rapport au mois cible
plot.month.ana <- function(Q="05",k,nbdays,dist,start="1950-01-01",end="2011-12-31",rean){
  
  dates <- getdates(start,end)
  
  # Import des analogues
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/save.nei.A/nei_",
              dist,"_member1_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei <- nei[[which(names(nei)==Q)]]
  
  # Traitement
  fonc <- function(v){
    tmp <- as.numeric(substr(dates[v],6,7)) # on recupere les mois
    tmp1 <- tmp[-1]-tmp[1] # on soustrait au mois du jour
    tmp1[tmp1>6] <- tmp1[tmp1>6]-12 # si >6, on enleve 12 (donc mois vont aller de -5 a +6)
    tmp1[tmp1< (-6)] <- tmp1[tmp1< (-6)]+12 # si <-6, on ajoute 12 (donc mois vont aller de -6 a +5)
    tmp1
  }
  
  nei.mon <- unlist(lapply(nei,fonc)) # longueur de 22645*112 analogues
  
  # Histogramme
  count <- hist(x = nei.mon,breaks=seq(-6.5,6.5,1),plot=F)$counts
  print(sum(count[6:8])/sum(count)*100) # analogues compris entre m-1 et m+1
  print(sum(count[5:9])/sum(count)*100) # analogues compris entre m-2 et m+2
  
  # Classique
  #hist(x=nei.mon,axes=F,breaks=seq(-6.5,6.5,1),col="azure3",xlab="Month(Target Day) - Month(Analog Days)",freq=F,main="")
  #grid(nx=NA,ny=NULL)
  #abline(v = c(-6.5:6.5), col = "grey", lty = "dotted",lwd=1)
  #hist(x=nei.mon,axes=F,breaks=seq(-6.5,6.5,1),col="azure3",xlab="Month(Target Day) - Month(Analog Days)",freq=F,main="",add=T)
  #lines(c(-6.6,6.5),c(0,0))
  #axis(2)
  #axis(1,at =-6:6,tick = F,line = -1)
  
  # ggplot2
  nei.mon <- as.data.frame(nei.mon)
  ggplot(nei.mon,aes(nei.mon)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12,colour="black",vjust=0),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "right",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_text(hjust=0.5,vjust=2,size = 12,face = "bold"))+
    geom_histogram(aes(y = stat(count) / sum(count)),breaks=seq(-6.5,6.5,1),col="darkblue",fill="grey")+
    xlab("Month(Analog Day) - Month(Target Day)")+
    ylab("Frequency")+
    scale_x_continuous(breaks=c(-6:6),labels=c(-6:6),minor_breaks = F)

    
}

# Boxplot de la repartition des analogues qu'on va chercher pour periode coupee en deux (homogene ou pas?)
plot.period.ana <- function(Q="05",k,nbdays,dist,start="1950-01-01",end="2011-12-31",rean){
  
  dates <- getdates(start,end)
  
  # Import des analogues
  load(paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/save.nei.A/nei_",
              dist,"_member1_k",k,"_mean",nbdays,"day_",start,"_",end,".Rdata"))
  nei <- nei[[which(names(nei)==Q)]]
  
  # Traitement
  yr <- as.numeric(unique(substr(dates,1,4))) # annees
  mid <- which(dates==paste0(as.character(yr[1]+round(length(yr)/2,0)-1),"-12-31")) # position du milieu de la periode
  
  nei.inf <- unlist(lapply(nei,function(v){sum(v[-1]<mid)/length(v[-1])*100})) # % des analogues de chaque jour dans la premiere periode
  nei.sup <- unlist(lapply(nei,function(v){sum(v[-1]>mid)/length(v[-1])*100})) # % des analogues de chaque jour dans la premiere periode
  
  # Mise en forme
  period1 <- paste0(start," -\n",dates[mid])
  period2 <- paste0(as.character(as.Date(dates[mid])+1)," -\n ",end)
  
  target_period <- c(rep(period1,mid),rep(period2,length(nei)-mid))
  mat <- as.data.frame(cbind(target_period,nei.inf,nei.sup))
  colnames(mat) <- c("target_period",period1,period2)
  mat <- pivot_longer(mat,2:3,names_to = "ana_period",values_to = "ana")
  mat$target_period <- factor(mat$target_period)
  mat$ana_period <- factor(mat$ana_period)
  mat$ana <- as.numeric(as.character(mat$ana))
  
  # Boxplot
  ggplot(mat, aes(x=target_period, y=ana, fill=ana_period)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12,colour="black",vjust=0),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "right",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_blank())+#element_text(hjust=0.5,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=1.5, linetype="dashed")+
    ylim(0,100)+
    xlab("Target period")+
    ylab("Analog days (%)")+
    labs(fill="Analog period")
}

# Graphique temporel de la repartition des max
plot.period.max <- function(bv1="Isere-seul",bv2="Drac-seul",nbdays,start="1950-01-01",end="2011-12-31",spazm=F){
  
  dates <- as.Date(getdates(start,end))
  
  # Max
  max.bv1 <- sort(get.ind.extr(bv1,nbdays,start,end,nei = T,spazm,seuil = "qua"))
  max.bv2 <- sort(get.ind.extr(bv2,nbdays,start,end,nei = T,spazm,seuil = "qua"))
  
  vect1 <- rep(0,length(dates));vect2 <- vect1
  vect1[max.bv1] <- 1; vect2[max.bv2] <- 1
  
  ann.1 <- aggregate(vect1,by=list(substr(dates,1,4)),sum)
  ann.2 <- aggregate(vect2,by=list(substr(dates,1,4)),sum)
  
  
  # Graphique
  png(filename = paste0("2_Travail/0_Present/Rresults/plot.period.max/plot.period.max_",bv1,"_",bv2,ifelse(spazm,"_spazm",""),".png"),
      width = 5,height = 2.5,units = "in",res = 300)
  
  # Nouvelle version
  layout(matrix(1:3,nrow = 3,ncol = 1,byrow = T),heights=c(1,1,0.1))
  par(mar=c(2,4,2,0))
  plot(ann.1,type="l",ylim=c(0,6),col="cornflowerblue",pty="n", main=nam2str(bv1),xlab="",ylab="Number of event")
  grid()
  lines(ann.1,col="cornflowerblue")
  
  plot(ann.2,type="l",ylim=c(0,6),col="burlywood1",pty="n", main=nam2str(bv2),xlab="Year",ylab="Number of event")
  grid()
  lines(ann.2,col="burlywood1")
  
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",bty="n")
  text(1.035,1,"Year")
  graphics.off()
  
  # Ancienne version
  #plot(ann.1[,1],rep(1,nrow(ann.1)),ylim=c(0,3),type="n",yaxt="n",xlab="Year",ylab="")
  #grid(ny=NA,nx=NULL)
  #abline(h=c(1,2),col="grey",lty=3)
  #par(new=T)
  #points(ann.1[,1],rep(2,nrow(ann.1)),cex=sqrt(ann.1[,2])*0.8,bg="cornflowerblue",pch=21)
  #points(ann.2[,1],rep(1,nrow(ann.2)),cex=sqrt(ann.2[,2])*0.8,bg="burlywood1",pch=21)
  #legend("topleft",bty="n",c(nam2str(bv1),nam2str(bv2)),pch=19,col=c("cornflowerblue","burlywood1"),cex=1.2,pt.cex=1.2,horiz=T,title = "Catchment")
  #legend("topright",bty="n",as.character(c(1,3,6)),pch=21,pt.cex=c(1*0.8,sqrt(3)*0.8,sqrt(6)*0.8),horiz=T,title="Number of event",cex=1.2)
  #graphics.off()
}
  
# Graphique de la saisonnalite de pwat et du boxplot de comparaison pour deux bvs
plot.pwat <- function(bv1,bv2,start,end,rean,spazm=c(F,F)){
  
  # Imports
  dates <- getdates(start,end)
  pwat <- get.pwat(k = 1,nbdays = 1,start,end,rean)
  ind.1 <- get.ind.max(type = "year",nbdays = 3,start = start,end = end,bv = bv1,spazm = spazm[1])
  ind.2 <- get.ind.max(type = "year",nbdays = 3,start = start,end = end,bv = bv2,spazm = spazm[2])
  
  # Traitement boxplot
  pw <- ecdf(pwat)(pwat)*100
  desais <- compute.desais(pwat,dates)
  desais <- ecdf(desais)(desais)*100
  
  bv <- rep(ifelse(spazm[1]==spazm[2],nam2str(bv1),paste0(nam2str(bv1)," (spazm=", spazm[1],")")),length(ind.1))
  mat1 <- as.data.frame(cbind(bv,pw[ind.1],desais[ind.1]))
  
  bv <- rep(ifelse(spazm[1]==spazm[2],nam2str(bv2),paste0(nam2str(bv2)," (spazm=", spazm[2],")")),length(ind.2))
  mat2 <- as.data.frame(cbind(bv,pw[ind.2],desais[ind.2]))
  mat <- rbind(mat1,mat2)
  colnames(mat) <- c("bv","PW","sea. adj. PW")
  
  mat <- pivot_longer(mat,2:3,names_to = "type",values_to = "percentile")
  mat$percentile <- as.numeric(as.character(mat$percentile))
  mat$type <- factor(mat$type,levels = c("PW","sea. adj. PW"))
  
  # Traitement saisonalite
  comm <- ind.1==ind.2
  ind.1 <- ind.1[!comm]
  ind.2 <- ind.2[!comm]
  ind.comm <- ind.1[comm]
  year <- substr(seq(as.Date("2020-01-01"),as.Date("2020-12-31"),"days"),6,10) # pour ajouter les pts sur la saisonalite
  
  # Graphiques
  png(filename = paste0("2_Travail/0_Present/Rresults/plot.pwat/plot_pwat_",bv1,"_",bv2,"_spazm_",spazm[1],"_",spazm[2]),width = 11,height = 5,units = "in",res = 1200)
  par(mfrow=c(1,2))
  
  plot.sais.quantile(dates = dates,vec = pw,label = "PW percentile (%)")
  #points(match(substr(dates[ind.1],6,10),year),pw[ind.1],pch=22,bg="cornflowerblue",cex=1.3)
  #points(match(substr(dates[ind.2],6,10),year),pw[ind.2],pch=22,bg="burlywood1",cex=1.3)
  #points(match(substr(dates[ind.comm],6,10),year),pw[ind.comm],pch=22,bg="olivedrab3",cex=1.3)
  #legend("topleft",inset=.02,bty="n",fill=c("cornflowerblue","burlywood1","olivedrab3"),
  #        c(nam2str(bv1),nam2str(bv2),"common"))
  
  plot.new() # manip pour ajouter un ggplot au plot classique          
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp1 <-plotViewport(c(1,0,3,0))
  
  p <- ggplot(mat, aes(x=type, y=percentile, fill=bv)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12,colour="black",vjust=0),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "right",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_blank())+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=1.5, linetype="dashed")+
    xlab("")+
    ylab("Percentile (%)")+
    labs(fill="Catchment")
  print(p,vp=vp1)
  graphics.off()
}

# Boxplot de pwat pour les max de precip zonaux et meridionaux
plot.pwat.flow <- function(nbdays,start,end,rean){
  
  # Imports
  dates <- getdates(start,end)
  pwat <- get.pwat(k = 1,nbdays,start,end,rean)
  ind.1 <- get.ind.max.flow(flow=1,agreg=T,nbdays,start,end)
  ind.2 <- get.ind.max.flow(flow=2,agreg=T,nbdays,start,end)
  
  # Traitement
  wp <- get.wp(nbdays,start,end,agreg=T)
  pwat.1 <- pwat
  pwat.1[wp!=1] <- NA
  pwat.1 <- ecdf(pwat.1)(pwat.1)*100
  
  pwat.2 <- pwat
  pwat.2[wp!=2] <- NA
  pwat.2 <- ecdf(pwat.2)(pwat.2)*100
  
  mat <- data.frame(Zonal=pwat.1[ind.1],Meridional=c(pwat.2[ind.2],rep(NA,length(ind.1)-length(ind.2))))
  mat <- pivot_longer(mat,1:2,names_to = "flow",values_to = "percentile")
  mat <- na.omit(mat)
  mat$flow <- factor(mat$flow,levels=c("Zonal","Meridional"))
  
  # Graphique
  ggplot(mat, aes(x=flow, y=percentile)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0.5,0,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12,colour="black",vjust=0),
          axis.text.y = element_text(size=12))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue",fill=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=1.5, linetype="dashed")+
    xlab("")+
    ylab("PW Percentile (%)")
  
  ggsave(filename = "2_Travail/0_Present/Rresults/plot.pwat.flow/plot_pwat_flow.png",width = 5,height = 4)
  graphics.off()
}

# Graphique de la saisonnalite de pwat avec les max annuels des deux influences Z et M
plot.pwat.sais <- function(flow=c(1,2),agreg=T,start="1950-01-01",end="2011-12-31",rean,spazm=F,small_win=F,supseuil=F){
  
  # Imports
  dates <- getdates(start,end)
  year <- substr(seq(as.Date("2020-01-01"),as.Date("2020-12-31"),"days"),6,10) # pour ajouter les pts sur la saisonalite
  
  pw <- get.pwat(k = 1,nbdays = 1,start,end,rean,small_win = small_win)
  
  ind.1 <- get.ind.max.flow(flow = flow[1],agreg = agreg,nbdays = 3,start = start,end = end,spazm = spazm,supseuil = supseuil)
  ind.2 <- get.ind.max.flow(flow = flow[2],agreg = agreg,nbdays = 3,start = start,end = end,spazm = spazm,supseuil = supseuil)
  
  # Graphique
  png(filename = paste0("2_Travail/0_Present/Rresults/plot.pwat.sais/plot_pwat_sais",ifelse(spazm,"_spazm",""),ifelse(small_win,"_small_win",""),ifelse(supseuil,"_supseuil",""),"_",rean,".png"),width = 7,height = 4,units = "in",res = 200)
  par(mar=c(2,4,1,1))
  plot.sais.quantile(dates = dates,vec = pw,y.label = "PW (mm)",colo="red")
  points(match(substr(dates[ind.1],6,10),year),pw[ind.1],pch=21,bg="cornflowerblue",cex=1.2)
  points(match(substr(dates[ind.2],6,10),year),pw[ind.2],pch=21,bg="burlywood1",cex=1.2)
  abline(h=median(pw[ind.1]),lwd=2,lty=2,col="cornflowerblue")
  abline(h=median(pw[ind.2]),lwd=2,lty=2,col="burlywood1")
  legend("topleft",inset=.02,bty="n",pt.bg=c("cornflowerblue","burlywood1"),pch=21,pt.cex = 1.2,
          c("Atlantic","Mediterranean"))
  #legend(5,27,bty="n",col=c("cornflowerblue","burlywood1"),lty=2,lwd=2,
  #       c("Zonal (median)","Meridional (median)"))
  graphics.off()
}

# Boxplot des percentiles des indicateurs avant et apres desaisonalisation
plot.quant.descr.desais <- function(descr=c("celnei","singnei","rsingnei","dP"),k,dist,nbdays,start,end,rean,bv,bvcomm=F){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Indicateurs
  tab <- matrix(NA,length(dates),length(descr)*2)
  
  for(i in 1:length(descr)){
    tmp <- get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,rean=rean)
    tab[,(i*2-1)] <- tmp
    tab[,i*2] <- compute.desais(tmp,dates)
  }
  
  # Humidite
  pwat <- get.pwat(k,nbdays,start,end,rean)
  desais <- compute.desais(pwat,dates)
  pwat <- cbind(pwat,desais)
  
  # Percentile
  tab <- cbind(tab,pwat)
  tab <- apply(tab,2,function(x) ecdf(x)(x)*100)
  colnames(tab)[seq(1,length(descr)*2,by=2)] <- descr
  colnames(tab)[seq(2,length(descr)*2,by=2)] <- "desais"
  
  # Graphique
  ind <- get.ind.max(type = "year",nbdays,start,end,bv)
  if(bvcomm!=F){
    indbis <- get.ind.max(type = "year",nbdays,start,end,bvcomm)
    comm <- which(ind==indbis)
    ind <- ind[-comm]
  }
  
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.quant.descr.desais/plot_",bv,"_",nbdays,"day_",start,"_",end,ifelse(bvcomm!=F,"_comm",""),".png"),width = 700,height = 500,units = "px")
  boxplot(tab[ind,],col="cornflowerblue",main=bv)
  abline(v=seq(2.5,length(descr)*2+0.5,by=2),lty=2,col="grey")
  graphics.off()
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
  png(file = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.sais.all/plot_sais_all_",descr,".png"),width = 8,height = 6,units = "in",res = 1200)
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

# Saisonalite des indicateurs de maniere propre pour article
plot.sais.descr.clean <- function(k,dist,nbdays=1,start="1950-01-01",end="2011-12-31",rean){
  
  descr <- c("celnei","singnei","rsingnei","dP")
  nam <- c("Celerity","Singularity","Relative singularity","MPD")
  dates <- getdates(start,end)
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.sais.descr.clean/plot_sais_",dist,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width = 8,height = 6,units = "in",res=300)
  par(mfrow=c(2,2),mar=c(3,4,2,2))
  
  for(i in 1:length(descr)){
    des <- get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,rean=rean)
    plot.sais.quantile(dates = dates,vec = des,y.label = paste0(nam2str(descr[i]),ifelse(i==4," (m)","")),main=nam[i])
  }
  graphics.off()
}

# Saisonnalite des extremes
plot.sais.extr <- function(nbdays,start="1950-01-01",end="2011-12-31",bv1,bv2,comm=F,spazm=F,supseuil=T){
  
  # Extremes
  dates <- getdates(start,end)
  if(!supseuil){
  month.bv1 <- as.numeric(substr(dates[get.ind.max(type = "year",nbdays,start,end,bv1,spazm)],6,7))
  month.bv2 <- as.numeric(substr(dates[get.ind.max(type = "year",nbdays,start,end,bv2,spazm)],6,7))
  }else{
    month.bv1 <- as.numeric(substr(dates[get.ind.extr(bv1,nbdays,start,end,nei = T,spazm,seuil = "qua")],6,7))
    month.bv2 <- as.numeric(substr(dates[get.ind.extr(bv2,nbdays,start,end,nei = T,spazm,seuil = "qua")],6,7))
  }
  
  if(comm){
    commun <- which(get.ind.max(type = "year",nbdays,start,end,bv1,spazm=T)==get.ind.max(type = "year",nbdays,start,end,bv2,spazm=T))
    month.bv1 <- month.bv1[-commun]
    month.bv2 <- month.bv2[-commun]
  }
  
  # Cumuls mensuels
  ny <- as.numeric(substr(dates[length(dates)],1,4)) - as.numeric(substr(dates[1],1,4)) + 1
  cum.bv1 <- aggregate(get.precip(nbdays=1,start,end,bv1,spazm),by=list(substr(dates,6,7)),function(v) sum(v)/ny)
  cum.bv2 <- aggregate(get.precip(nbdays=1,start,end,bv2,spazm),by=list(substr(dates,6,7)),function(v) sum(v)/ny)
  
  # Graphique
  fac1 <- max(cum.bv1[,2])/(max(table(month.bv1))*0.8) # facteur de division des cumuls mensuels pour rentrer dans l'histogramme
  fac2 <- max(cum.bv2[,2])/(max(table(month.bv2))*0.8)
  fac <- round(mean(c(fac1,fac2)),0)
  
  png(filename = paste0("2_Travail/0_Present/Rresults/plot.sais.extr/plot_ann_max_",bv1,"_",bv2,ifelse(spazm,"_spazm",""),ifelse(comm,"_comm",""),ifelse(supseuil,"_qua",""),".png"),width = 9,height = 3.5,units = "in",res=200)
  par(mfrow=c(1,2),lwd=2,mar=c(2,3,3,6))
  
  hist(month.bv1,axes=F,breaks=0:12,col="azure3",
       xlab="",ylab="",main=nam2str(bv1))
  abline(h = axis(2), col = "grey", lty = "dotted",lwd=1)
  abline(v = c(2,5,8,11), col = "grey", lty = "dotted",lwd=1)
  hist(month.bv1,axes=F,breaks=0:12,col="azure3",
       xlab="Month",ylab="Count",ylim=c(0,ifelse(comm,12,15)),main=nam2str(bv1),add=T)
  lines(c(0,12),c(0,0))
  axis(2)
  title(ylab="Count", line=2)
  text(0.5:11.5,par("usr")[3]-0.3, labels = month.abb, srt = -45, pos = 1, xpd = TRUE)
  points(0.5:11.5,cum.bv1[,2]/fac,pch=19)
  lines(0.5:11.5,cum.bv1[,2]/fac)
  axis(side = 4,at = axis(2),labels = axis(2)*fac)
  text(15,par("usr")[1]++mean(axis(2)),"Precipitation (mm)", srt=90,xpd=T)
  
  hist(month.bv2,axes=F,breaks=0:12,col="azure3",
      xlab="",ylab="",main=nam2str(bv2))
  abline(h = axis(2), col = "grey", lty = "dotted",lwd=1)
  abline(v = c(2,5,8,11), col = "grey", lty = "dotted",lwd=1)
  hist(month.bv2,axes=F,breaks=0:12,col="azure3",
       xlab="Month",ylab="Count",ylim=c(0,ifelse(comm,12,15)),main=nam2str(bv2),add=T)
  lines(c(0,12),c(0,0))
  axis(2)
  title(ylab="Count", line=2)
  text(0.5:11.5,par("usr")[3]-0.3, labels = month.abb, srt = -45, pos = 1, xpd = TRUE)
  points(0.5:11.5,cum.bv2[,2]/fac,pch=19)
  lines(0.5:11.5,cum.bv2[,2]/fac)
  axis(side = 4,at = axis(2),labels = axis(2)*fac)
  text(15,par("usr")[1]+mean(axis(2)),"Precipitation (mm)", srt=90,xpd=T)
  
  graphics.off()
}

# Figure des cartes composites et nuages de pts WP
plot.wp <- function(wp=c(1,2),rean,k,descriptors,dist,nbdays=3,start="1950-01-01",end="2011-12-31",radtype="nrn05",CV=TRUE,threeday=c(F,F),quant=F,agreg=F,title=""){
  
  # Graphique
  png(filename = paste0("2_Travail/0_Present/",rean,"/Rresults/overall/k",k,"/plot.wp/plot_wp_",descriptors[1],"_",descriptors[2],"_wp",wp[1],"_wp",wp[2],"_",ifelse(agreg,"agreg_",""),dist,"_member",member,"_k",k,"_mean",nbdays,"day_",start,"_",end,".png"),width = 8,height = 7,units = "in",res=1200)
  layout(matrix(1:6,2,3,byrow=F),width=c(1.3,1.3,0.4))
  par(pty="s",mar=c(5,7,6,1))
  let <- list(c("a)","c)"),c("b)","d)"))
  
  # Cartes composites
  for(i in 1:length(wp)){
    map.composite.wp(wp[i],k[1],start,end,rean[1],leg=F,let=let[[i]][1],agreg=agreg,title=title[i])
    plot.empir.wp(wp[i],rean,k,descriptors,dist,nbdays,start,end,radtype,CV,threeday,save=F,let=let[[i]][2],quant,agreg=agreg,title=title[i])
  }
  
  # Legende cartes
  leg <- seq(-200,200,50)
  N <- 11

  par(pty="m",mar=c(5,0,6,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",xlim=c(0,1),ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = paste0("        ",leg),at =  seq(0, 1,length.out = length(leg)),
            vertical = T,xlim = c(0.2,0.4),ylim = c(0,1),cex=1.4)
  #text(x = 1,y = 0.75,"Geopotential height (m)",cex=1.6)
  
  graphics.off()
}

# Run les fonctions repetitives
run.article <- function(type=1){
  
  # Type = 1: compute.density()
  if(type==1){
    # Parametres
    descr <- list(
      c("celnei","dP"),
      c("singnei","dP"),
      c("rsingnei","dP")
    )
    wp <- c(1,2)
    #wp <- F
    #wp <- 1:8
    #spazm <- c(F,T)
    spazm <- T
    quant <- c(T,F)
    
    # Fonction
    for(i in 1:length(descr)){
      for(j in 1:length(wp)){
        for(k in 1:length(spazm)){
          for(l in 1:length(quant)){
            compute.density(rean="ERA5",k = 1,descriptors = descr[[i]],dist = c("TWS","TWS"),nbdays = 3,
                            start = "1950-01-01",end="2017-12-31",quant = quant[l],wp = wp[j],agreg = T,spazm = spazm[k])
          }
        }
      }
    }
  }
  
  # Type = 2: illustr.precip.seq()
  if(type==2){
    bv <- c("Isere","Isere-seul","Drac-seul")
    nbdays <- c(2,3,4,5)
    spazm <- c(F,T)
    nei <- c(F,T)
    
    for(i in 1:length(bv)){
      for(j in 1:length(nbdays)){
        for(k in 1:length(spazm)){
          for(l in 1:length(nei)){
            illustr.precip.seq(bv = bv[i],nbdays = nbdays[j],
                               start = "1950-01-01",end="2011-12-31",spazm = spazm[k],nei = nei[l])
          }
        }
      }
    }
  }
}

