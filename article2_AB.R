source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Combninaison de compare.descr.max.flow
combine.compare.descr.max.flow <- function(wp=c(1,2),sais=c("winter","autumn"),bv=c("Isere-seul","Drac-seul"),descr=c("cel","sing05","rsing05","dP"),k,dist,nbdays,start,end,rean,spazm=T,all=F){
  
  # Graphiques
  p <- compare.descr.max.flow(wp = wp,sais = sais[1],bv = bv[1],descr = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,rean = rean,spazm = spazm,all = all)
  q <- compare.descr.max.flow(wp = wp,sais = sais[2],bv = bv[2],descr = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,rean = rean,spazm = spazm,all = all)
  
  # Combinaison
  ggarrange(p,q,ncol = 2, nrow = 1,widths = c(1,1),common.legend = TRUE, legend = "bottom")
  ggsave(filename = paste0(get.dirstr(k,rean,"past"),"compare.descr.max.flow/plot_descr_max_flow_k",k,"_mean",nbdays,"day_",bv[1],"_",sais[1],"_",bv[2],"_",sais[2],"_",start,"_",end,".png"),width = 12,height = 4)
  graphics.off()
}

# Combinaison de map.diff.geo et map.diff.rean
combine.map.diff.geo.rean <- function(k,rean=c("20CR-m1","ERA20C"),signif=T){
  
  # Parametres
  per <- vector(mode = "list",length = 2)
  per[[1]] <- c("1900-01-01","1929-12-31")
  per[[2]] <- c("1970-01-01","1999-12-31")
  per.name <- per
  per.name <- lapply(per.name,function(v){v[2]<-as.numeric(substr(v[2],1,4))+1;paste(substr(v,1,4),collapse = "-")})
  
  rean.maprean <- rean
  rean.maprean[2] <- paste0(rean.maprean[2],"_regrid_20CR")
  
  # Cartes
  png(filename = paste0("2_Travail/1_Past/Rresults/map.diff.geo/map_diff_rean_k",k,"_",rean[1],"_",rean[2],".png"),width = 8,height = 10,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,4),2:5,rep(6,4),7:10,rep(11,4),12:15,rep(16,4),17:20,rep(21,4)),nrow = 9,ncol = 4,byrow = T),widths = c(1,1,1,1),heights = c(0.1,1,0.1,1,0.1,1,0.1,1,0.4))
  
  # map.diff.rean
  for(i in 1:length(per)){
  
    # Titre
    par(pty="m",mar=c(0,0,0,0))
    plot(1,1,type="n",xaxt="n",yaxt="n",bty="n",xlim=c(0,2),ylim=c(0,2))
    text(1,1,paste0(rean[1]," minus ",rean[2]," - ",per.name[[i]]),cex=2,font=2)
    
    # Cartes
    par(mar=c(0.5,0.5,1.5,0.5),pty="s")
    map.diff.rean(k = k,rean = rean.maprean,start = per[[i]][1],end = per[[i]][2],signif = T)
  }
  
  # map.diff.geo
  for(i in 1:length(rean)){
    
    # Titre
    par(pty="m",mar=c(0,0,0,0))
    plot(1,1,type="n",xaxt="n",yaxt="n",bty="n",xlim=c(0,2),ylim=c(0,2))
    text(1,1,paste0(rean[i]," - ",per.name[[2]]," minus ",per.name[[1]]),cex=2,font=2)
    
    # Cartes
    par(mar=c(0.5,0.5,1.5,0.5),pty="s")
    map.diff.geo(k = k,rean = rean[i],signif = signif,save = F)
  }
  
  # Legende
  N <- 11
  lab <- seq(-75,75,25)
  
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = lab,at =  seq(0,1,length.out = length(lab)),
              vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
  text(x = 1,y = 0.9,"Geopotential Height difference (m)",cex=1.2,font=2)
  graphics.off()
}

# Combinaison de map.diff.geo.wp
combine.map.diff.geo.wp <- function(wp=c(8,1,2),k,rean,start="1950-01-01",end="2017-12-31"){
  
  tt <- c(1,2,5,8)
  wp.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  png(filename = paste0(get.dirstr(k,rean,"past"),"map.diff.geo.wp/map_diff_combine_k",k,"_wp",paste(wp,collapse = "_"),".png"),width = 8,height = 5.5,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,4),2:5,rep(6,4),7:10,rep(11,4)),nrow = 5,ncol = 4,byrow = T),widths = c(1,1,1,1),heights = c(0.1,1,0.1,1,0.4))
  
  for(i in 1:length(wp)){
  # Titre
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",bty="n",xlim=c(0,2),ylim=c(0,2))
  text(1,1,wp.name[tt==wp[i]],cex=2,font=2)
  
  # Cartes
  par(mar=c(0.5,0.5,1.5,0.5),pty="s")
  map.diff.geo.wp(wp = wp[i],k = k,rean = rean,start = start,end = end,save = F,signif = T)
  }
  
  # Legende
  dates <- getdates(start,end)
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4))
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  
  N <- 11
  lab <- seq(-60,60,15)
  
  
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = lab,at =  seq(0,1,length.out = length(lab)),
              vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
  text(x = 1,y = 0.9,paste0(period2," minus ",period1," Geopotential Height (m)"),cex=1.2,font=2)
  graphics.off()
  
}

# Combinaison de plusieurs plot.descr.density.subperiod
combine.plot.descr.density.subperiod <- function(wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  descr <- c("cel","sing05","rsing05","dP")
  pl <- vector(mode="list",length=length(descr))
  for(i in 1:length(descr)){
    pl[[i]] <- plot.descr.density.subperiod(descr = descr[i],wp = wp,k = k,dist = dist,nbdays = nbdays,
                                            rean = rean,start = start,end = end,start.ana = start.ana,end.ana = end.ana,add.nb.density = ifelse(i==1,T,F))
  }
  
  ggarrange(plotlist = pl,
            ncol = 2, nrow = 2,
            widths = c(1,1),
            common.legend = TRUE, legend = "bottom",
            align = "hv")
  
  if(is.null(wp)) wp <- "none"
  ggsave(filename = paste0(get.dirstr(k,rean,"past"),"plot.descr.density.subperiod/plot_combine_density_k",k,"_mean",nbdays,"day_wp=",wp,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),width = 10,height = 10)
  graphics.off()
}

# Combinaison de plusieurs plot.descr.rain.norain
combine.plot.descr.rain.norain <- function(bv=c("Isere-seul","Drac-seul"),wp=c(1,2),sais=c("spring","winter"),k,dist,nbdays=1,rean,spazm,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  # Graphiques
  pl <- list()
  for(i in 1:length(bv)){
    pl[[i]] <- plot.descr.rain.norain(bv = bv[i],wp = wp[i],sais = sais[i],k = k,dist = dist,nbdays = nbdays,rean = rean,
                                      spazm = spazm,start = start,end = end,start.ana = start.ana,end.ana = end.ana,save = F)
    pl[[i]] <- pl[[i]]+labs(title=paste(nam2str(bv),wp.name[wp],nam2str(sais),sep = " - ")[i])
  }
  
  # Combinaison
  wp.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  main <- paste(nam2str(bv),wp.name[wp],nam2str(sais),sep = " - ")
  pl.final <- ggarrange(plotlist = pl,ncol = 2,nrow = 1,font.label = list(size=20,face="bold"),legend = "right",common.legend = T)
  ggsave(filename = paste0(get.dirstr(k,rean,period = "past"),"plot.descr.rain.norain/plot_descr_",paste(bv,collapse = "_"),"_",paste(sais,collapse = "_"),"_wp",paste(wp,collapse = "_"),"_",nbdays,"days",ifelse(spazm,"_spazm",""),"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),plot = pl.final,width = 14,height = 5)
}

# Combinaison de plusieurs plot.descr.violin.subperiod
combine.plot.descr.violin.subperiod <- function(k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  # Parametres
  wp <- c(1,2,5,8)
  wp.name <- c("All","Atlantic","Mediterranean","Northeast","Anticyclonic")
  descr <- c("cel","sing05","rsing05","dP")
  
  # Graphiques individuels et par types de temps
  pl <- vector(mode="list",length=(length(wp)+1)*length(descr))
  pl.wp <- vector(mode="list",length=length(wp)+1)
  
  for(i in 1:(length(wp)+1)){
    
    if(i==1){wp.i <- NULL
    }else{wp.i<-wp[i-1]}
    
    for(j in 1:length(descr)){
      pos <- length(descr)*(i-1)+j
      print(pos)
      pl[[pos]] <- plot.descr.violin.subperiod(descr = descr[j],wp = wp.i,k = k,dist = dist,nbdays = nbdays,rean = rean,
                                               start = start,end = end,start.ana = start.ana,end.ana = end.ana,add.nb.box = F)#ifelse(j==1,T,F))
    }
  }
    
  for(i in 1:(length(wp)+1)){
    pos.i <- length(descr)*(i-1)+(1:4)
    if(i==5){
      leg.i <- "bottom";comm.leg.i <- T
    }else{leg.i <- F;comm.leg.i <- F}
    
    pl.wp[[i]] <-  ggarrange(plotlist=pl[pos.i], ncol = 4, nrow = 1,
                             legend=leg.i,common.legend = comm.leg.i,labels=NULL)
    pl.wp[[i]] <- annotate_figure(pl.wp[[i]], top = text_grob(wp.name[i], face = "bold", size = 18))
  }
    
  # Graphique complet
  pl.final <- ggarrange(plotlist=pl.wp, ncol = 1, nrow = 5, legend = F, heights = c(1,1,1,1,1.2))
  pl.final <- annotate_figure(pl.final, left = text_grob("Percentile of descriptor value (%)", face = "bold", size = 16,rot=90))
  
  ggsave(filename = paste0(get.dirstr(k,rean,"past"),"plot.descr.violin.subperiod/plot_combine_violon_k",k,"_mean",nbdays,"day_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),
         plot = pl.final,height = 24,width = 22)
  graphics.off()
}

# Combinaison de plusieurs plot.sais.violin.subperiod
combine.plot.sais.violin.subperiod <- function(wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",combine=F){
  
  # Parametres
  sais <- c("winter","spring","summer","autumn")
  
  # Graphiques
  pl <- vector(mode="list",length=length(sais))
  for(i in 1:length(sais)){
    pl[[i]] <- plot.sais.violin.subperiod(sais = sais[i],wp = wp,k = k,dist = dist,nbdays = nbdays,
                                 rean = rean,start = start,end = end,start.ana = start.ana,end.ana = end.ana,add.nb.box = T)
  }
  
  if(!combine){
    pl.final <- ggarrange(plotlist=pl, ncol = 2, nrow = 2, common.legend = T, legend = "bottom")
    pl.final <- annotate_figure(pl.final, left = text_grob("Percentile of descriptor value (%)", face = "bold", size = 16,rot=90))
    if(is.null(wp)) wp <- "none"
    ggsave(filename = paste0(get.dirstr(k,rean,"past"),"plot.sais.violin.subperiod/plot_combine_violon_k",k,"_mean",nbdays,"day_wp=",wp,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),
           plot = pl.final,height = 8,width = 14)
    graphics.off()
  }
  
  # Si second wp
  if(combine){
    tt <- c(1,2,5,8)
    wp.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
    wp2 <- 2
    
    pl.final <- ggarrange(plotlist=pl, ncol = 2, nrow = 2, legend = F)
    pl.final <- annotate_figure(pl.final, top = text_grob(wp.name[tt==wp], face = "bold", size = 28))
    
    pl2 <- vector(mode="list",length=length(sais))
    for(i in 1:length(sais)){
      pl2[[i]] <- plot.season.violin.subperiod(sais = sais[i],wp = wp2,k = k,dist = dist,nbdays = nbdays,
                                              rean = rean,start = start,end = end,start.ana = start.ana,end.ana = end.ana,add.nb.box = T)
    }
    pl.final2 <- ggarrange(plotlist=pl2, ncol = 2, nrow = 2, common.legend = T, legend = "bottom")
    pl.final2 <- annotate_figure(pl.final2, top = text_grob(wp.name[tt==wp2], face = "bold", size = 28))
    
    pl.combine <- ggarrange(pl.final,pl.final2,ncol = 1,nrow = 2,heights = c(1,1.1))
    pl.combine <- annotate_figure(pl.combine, left = text_grob("Percentile of descriptor value (%)", face = "bold", size = 20,rot=90))
    ggsave(filename = paste0(get.dirstr(k,rean,"past"),"plot.sais.violin.subperiod/plot_combine_violon_k",k,"_mean",nbdays,"day_wp=",wp,"_wp=",wp2,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),
           plot = pl.combine,height = 18,width = 14)
    graphics.off()
  }
}

# Combinaison de plusieurs plot.season.diff.density.subperiod
combine.plot.sais.diff.density.subperiod <- function(wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",ref="all"){
  
  # Parametres
  sais <- c("winter","spring","summer","autumn")
  
  # Graphiques
  pl <- vector(mode="list",length=length(sais))
  for(i in 1:length(sais)){
    pl[[i]] <- plot.sais.diff.density.subperiod(sais = sais[i],wp = wp,k = k,dist = dist,nbdays = nbdays,
                                                  rean = rean,start = start,end = end,start.ana = start.ana,end.ana = end.ana,add.nb.box = T,ref = ref)
  }
  
  pl.final <- ggarrange(plotlist=pl, ncol = 2, nrow = 2, common.legend = T, legend = "bottom")
  pl.final <- annotate_figure(pl.final, left = text_grob("Percentile of descriptor value (%)", face = "bold", size = 16,rot=90))
  if(is.null(wp)) wp <- "none"
  ggsave(filename = paste0(get.dirstr(k,rean,"past"),"plot.sais.diff.density.subperiod/plot_combine_diff_density_k",k,"_mean",nbdays,"day_wp=",wp,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,"_ref=",ref,".png"),
         plot = pl.final,height = 8,width = 14)
  graphics.off()
}

# Combinaison de plusieurs plot.season.diff.nbdays.subperiod
combine.plot.sais.diff.nbdays.subperiod <- function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",ref="all"){
  
  pl <- vector(mode="list",length=length(sais))
  
  for(i in 1:length(sais)){
    pl[[i]] <- plot.sais.diff.nbdays.subperiod(sais = sais[i],wp = wp,k = k,dist = dist,nbdays = nbdays,rean = rean,
                                      start = start,end = end,start.ana = start.ana,end.ana = end.ana,ref = ref)
  }
  
  pl.final <- ggarrange(plotlist=pl, ncol = length(sais), nrow = 1, common.legend = T, legend = "bottom")
  ggsave(filename = paste0(get.dirstr(k,rean,"past"),"plot.sais.diff.nbdays.subperiod/plot_combine_diff_nbdays_k",k,"_",paste(sais,collapse="_"),"_wp=",wp,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,"_ref=",ref,".png"),
         plot = pl.final,height = 5,width = 10)
  graphics.off()
}

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
  dates.ana <- c("1950-01-01","2010-12-31")
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.descr/plot_combine_trend_",namind,"_",namsais,"_liss=",liss,ifelse(ana.comm,paste0("_ana_",dates.ana[1],"_",dates.ana[2]),""),ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 9,height = 6,units = "in",res=600)
  par(mfrow=c(2,2))
  for(i in 1:length(sais)){
    print(paste0("plot ",i,"/",length(sais)))
    leg <- ifelse(i==1,T,F)
    plot.trend.descr(descr = ind,k = 1,dist = dist,sais = sais[i],liss = liss,
                     ana.comm = ana.comm,align = align,nao = nao,leg = leg,save = F,type="season")
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.descr
combine.plot.trend.descr.all <- function(k,dist,liss=1,ana.comm=F,align=F,nao=F){
  
  # Indicateurs
  ind <- c("cel","sing05","rsing05","dP")
  namind <- ifelse(length(ind)==1,ind,"alldescr")
  
  # Saisons
  sais <- c("winter","spring","summer","autumn")
  namsais <- ifelse(length(sais)==1,sais,"allsais")
  
  # Graphique
  dates.ana <- c("1950-01-01","2010-12-31")
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.descr/plot_combine_trend_",namind,"_",namsais,"_liss=",liss,ifelse(ana.comm,paste0("_ana_",dates.ana[1],"_",dates.ana[2]),""),ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 12,height = 12,units = "in",res=600)
  layout(mat = matrix(data = 1:25,nrow = 5,ncol = 5,byrow = 5),widths = c(0.2,1,1,1,1),heights = c(0.2,1,1,1,1))
  
  # Titres
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",bty="n")
  for(i in 1:length(sais)){
    plot(1,1,type="n",xaxt="n",yaxt="n",bty="n")
    text(1,1,nam2str(sais[i]),cex=2,font=2)
  }
  

  # ylab et graphiques
  for(i in 1:length(ind)){
    print(ind[i])
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",xaxt="n",yaxt="n",bty="n")
    text(1,1.1,nam2str(ind[i],whole = T,unit = T),cex=2,font=2,srt=90)
    par(mar=c(3,1.5,0,1))
    for(j in 1:length(sais)){
      print(paste0("plot ",j,"/",length(sais)))
      leg <- F; if(i==1 & j==1) leg <- T
      plot.trend.descr(descr = ind[i],k = 1,dist = dist,sais = sais[j],liss = liss,
                       ana.comm = ana.comm,align = align,nao = nao,leg = leg,save = F,type="season")
    }
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.descr.extr.combine, avec ajout de plot.trend.alldescr.nbr.qua
combine.plot.trend.descr.extr.alldescr <- function(bv=c("Isere-seul","Drac-seul"),sais=c("winter","autumn"),k,dist,rean,nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  png(filename = paste0(get.dirstr(k,rean,period="past"),"plot.trend.descr.extr/plot_",paste(bv,collapse = "_"),"_",paste(sais,collapse = "_"),"_mean",nbdays,"day_",rean,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),width = 9,height = 4,units = "in",res=600)
  par(mfrow=c(1,length(bv)),mar=c(2,4,2,0.5))
  
  # plot.trend.descr.extr.combine
  for(i in 1:length(bv)){
    plot.trend.descr.extr.alldescr(bv = bv[i],sais = sais[i],k = k,dist = dist,rean = rean,nbdays = nbdays,spazm = spazm,start = start,end = end,start.ana = start.ana,end.ana = end.ana,leg=ifelse(i==1,T,F))
  }
  
  # plot.trend.descr.nbr.qua
  #for(i in 1:length(sais)){
  #  plot.trend.alldescr.nbr.qua(qua = 0.25,k = k,dist = dist,nbdays = nbdays,rean = rean,sais = sais[i],wp = "all",start = start,end = end,start.ana = start.ana,end.ana = end.ana,save = F)
  #}
  graphics.off()
}

# Combinaison de plusieurs plot.trend.descr.sais
combine.plot.trend.descr.sais <- function(k,dist,rean,nbdays=1,start="1950-01-01",end="2010-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  # Indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  
  # Graphique
  png(filename = paste0(get.dirstr(k,rean,period="past"),"plot.trend.descr.sais/plot_combine_",rean,"_mean",nbdays,"day_",start,"_",end,".png"),width = 9,height = 6,units = "in",res=600)
  par(mfrow=c(2,2))
  for(i in 1:length(descr)){
    print(paste0("plot ",i,"/",length(descr)))
    leg <- ifelse(i==1,T,F)
    plot.trend.descr.sais(descr = descr[i],k = k,dist = dist,rean = rean,nbdays = nbdays,
                          start = start,end = end,start.ana = start.ana,end.ana = end.ana,leg = leg,save = F)
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.precip.gev
combine.plot.trend.precip.gev <- function(bv=c("Isere-seul","Drac-seul"),nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31"){
  
  # graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip.gev/plot_trend_precip_gev_",bv[1],"_",bv[2],"_",nbdays,"day_",ifelse(spazm,"spazm_",""),start,"_",end,".png"),width = 10,height = 6,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,5),2:6,rep(7,5),8:12),nrow = 4,ncol = 5,byrow = T),widths = c(0.1,1,1,1,1),heights = c(0.1,1,0.1,1))
  
  for(i in 1:length(bv)){
    
    # Titre
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
    text(1,1,nam2str(bv[i]),font=2,cex=2.5)
    
    # Axe
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,"Precipitation (mm)",srt=90,cex=1.3)
    
    # Plot
    par(mar=c(3,2,2,0))
    plot.trend.precip.gev(bv = bv[i],nbdays = nbdays,spazm = spazm,start = start,end = end,leg = ifelse(i==1,T,F))
  }
  graphics.off()
}

# Combinaison de plusieurs plot.trend.precip.wp.combine
combine.plot.trend.precip.wp <- function(bv=c("Isere-seul","Drac-seul"),spazm=T,start="1950-01-01",end="2017-12-31",norm=T){
  
  # graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip.wp.combine/plot_trend_precip_",bv[1],"_",bv[2],"_",start,"_",end,ifelse(spazm,"_spazm",""),ifelse(norm,"_norm",""),".png"),width = 10,height = 6,units = "in",res=600)
  layout(mat = matrix(data = c(rep(1,5),2:6,rep(7,5),8:12),nrow = 4,ncol = 5,byrow = T),widths = c(0.1,1,1,1,1),heights = c(0.1,1,0.1,1))
  
  for(i in 1:length(bv)){
    
    # Titre
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",bty="n",xaxt="n",yaxt="n")
    text(1,1,nam2str(bv[i]),font=2,cex=2.5)
    
    # Axe
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    text(1,0.5,"Precipitation (mm/day)",srt=90,cex=1.3)
    
    # Plot
    par(mar=c(3,2,2,0))
    plot.trend.precip.wp.combine(bv = bv[i],spazm = spazm,start = start,end = end,norm = norm,leg = ifelse(i==1,T,F),save = F)
  }
  graphics.off()
}
  
# Combinaison de plusieurs combine.plot.trend.wp
combine.plot.trend.wp <- function(nbdays=1,start="1950-01-01",end="2017-12-31"){
  
  # Saisons
  sais <- c("winter","spring","summer","autumn")
  
  # Graphique occurence
  print("Occurence")
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.wp/plot_combine_occurence_wp_mean",nbdays,"day_",start,"_",end,".png"),width = 9,height = 6,units = "in",res=600)
  par(mfrow=c(2,2),mar=c(4,4,2,1))
  for(i in 1:length(sais)){
    print(paste0("plot ",i,"/",length(sais)))
    leg <- ifelse(i==1,T,F)
    plot.trend.wp(type = "occurence",sais = sais[i],nbdays = nbdays,start = start,end = end,leg = leg,save = F)
  }
  graphics.off()
  
  # Graphique persistance
  print("Persistance")
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.wp/plot_combine_persistence_wp_mean",nbdays,"day_",start,"_",end,".png"),width = 9,height = 6,units = "in",res=600)
  par(mfrow=c(2,2),mar=c(4,4,2,1))
  for(i in 1:length(sais)){
    print(paste0("plot ",i,"/",length(sais)))
    leg <- ifelse(i==1,T,F)
    plot.trend.wp(type = "persistence",sais = sais[i],nbdays = nbdays,start = start,end = end,leg = leg,save = F)
  }
  graphics.off()
  
  # Graphique occurence et persistance
  print("Occurence et Persistance")
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.wp/plot_combine_occurence_persistence_wp_mean",nbdays,"day_",start,"_",end,".png"),width = 10,height = 6,units = "in",res=600)
  layout(mat = matrix(data = 1:10,nrow = 2,ncol = 5,byrow = F),widths = c(0.1,rep(1,4)),heights = c(1,1))
  
  par(mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  text(1,0.5,"WP Occurrence (%)",srt=90,cex=1.3)
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  text(1,0.5,"WP Persistence (days)",srt=90,cex=1.3)
  
  par(mar=c(2,2,2,0))
  for(i in 1:length(sais)){
    print(paste0("plot ",i,"/",length(sais)))
    leg <- ifelse(i==1,T,F)
    plot.trend.wp(type = "occurence",sais = sais[i],nbdays = nbdays,start = start,end = end,leg = leg,save = F)
    plot.trend.wp(type = "persistence",sais = sais[i],nbdays = nbdays,start = start,end = end,leg = leg,save = F)
  }
  graphics.off()
}

# Combinaison de plusieurs plot.TWS.crossed
combine.plot.TWS.crossed <- function(k,start="1851-01-01",end="2010-12-31"){
  
  # Saisons
  sais <- c("winter","spring","summer","autumn")
  
  # Graphique
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.TWS.crossed/plot_combine_TWS_k",k,".png"),width = 9,height = 6,units = "in",res=600)
  par(mfrow=c(2,2))
  for(i in 1:length(sais)){
    print(paste0("plot ",i,"/",length(sais)))
    leg <- ifelse(i==1,T,F)
    plot.TWS.crossed(k = k,start = start,end = end,sais = sais[i],leg = leg,save = F)
  }
  graphics.off()
}

# Combinaison de scatterplot.density.subperiod
combine.scatterplot.density.subperiod <- function(sais,sais2=NULL,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  # Definition des couples d'indicateurs
  descr <- list(
    c("cel","dP"),
    c("sing05","dP"),
    c("rsing05","dP")
  )
  
  # Graphiques
  pl <- vector(mode="list",length=length(sais))
  
  for(i in 1:length(descr)){
    pl[[i]] <- scatterplot.density.subperiod(descr = descr[[i]],sais = sais,wp = wp,
                                             k = k,dist = dist,nbdays = nbdays,rean = rean,
                                             start = start,end = end,start.ana = start.ana,end.ana = end.ana)
  }
  
  if(!is.null(sais2)){
    for(i in 1:length(descr)){
      pl[[i+length(descr)]] <- scatterplot.density.subperiod(descr = descr[[i]],sais = sais2,wp = wp,
                                                             k = k,dist = dist,nbdays = nbdays,rean = rean,
                                                             start = start,end = end,start.ana = start.ana,end.ana = end.ana)
    }
  }
  
  # Export
  pl.final <- ggarrange(plotlist=pl, ncol = 3, nrow = ifelse(!is.null(sais2),2,1), common.legend = T, legend = "bottom")
  if(!is.null(sais2)){
    pl.final <- pl.final+geom_text(aes(x = 0.5,y=0.965,label=nam2str(sais),fontface=2),size=10)
    pl.final <- pl.final+geom_text(aes(x = 0.5,y=0.51,label=nam2str(sais2),fontface=2),size=10)
  }
  
  ggsave(filename = paste0(get.dirstr(k,rean,"past"),"scatterplot.density.subperiod/scatterplot_combine_density_k",k,"_mean",nbdays,"day_wp=",wp,"_",sais,"_",ifelse(!is.null(sais2),paste0(sais2,"_"),""),start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),
         plot = pl.final,height = ifelse(!is.null(sais2),9,4.5),width = 12)
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

# Comparaison des percentiles d'indicateurs pour deux flux
compare.descr.max.flow <- function(wp=c(1,2),sais="winter",bv="Isere-seul",descr=c("cel","sing05","rsing05","dP"),k,dist,nbdays,start,end,rean,spazm=T,all=F){
  
  # Import indicateurs
  mat <- NULL
  nam.descr <- NULL
  
  for(i in 1:length(descr)){
    mat <- cbind(mat,get.descriptor(descr[i],k,dist,nbdays,start,end,standardize = F,rean=rean,desais=F))
    nam.descr[i] <- nam2str(descr[i],whole=T)
    colnames(mat)[ncol(mat)] <- nam.descr[i]
  }
  
  # Traitement
  tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = spazm)
  
  namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  namflow <- namflow[wp]
  
  ind <- get.ind.max.sais(sais = sais,wp = "all",nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  ind.1 <- intersect(ind,which(tt==wp[1]))
  ind.2 <- intersect(ind,which(tt==wp[2]))
  
  mat.1 <- mat
  if(!all) mat.1[tt!=wp[1],] <- NA
  mat.1 <- apply(mat.1,2,function(v) {ecdf(v)(v)*100})
  mat.1 <- mat.1[ind.1,]
  
  mat.2 <- mat
  if(!all) mat.2[tt!=wp[2],] <- NA
  mat.2 <- apply(mat.2,2,function(v) {ecdf(v)(v)*100})
  mat.2 <- mat.2[ind.2,]
  
  mat <- as.data.frame(rbind(mat.1,mat.2))
  mat <- cbind(c(rep(namflow[1],length(ind.1)),rep(namflow[2],length(ind.2))),mat)
  colnames(mat)[1] <- "Flow"
  mat <- pivot_longer(mat,2:ncol(mat),names_to = "descr",values_to = "percentile")
  mat$Flow <- factor(mat$Flow,levels = namflow[wp])
  mat$descr <- factor(mat$descr,levels = nam.descr)
  
  # Graphique
  p <- ggplot(mat, aes(x=descr, y=percentile, fill=Flow)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0,0,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "right",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=16,colour="black"),
          legend.title = element_blank())+#element_text(hjust=0.4,vjust=2,size = 12,face = "bold"))+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c("cornflowerblue","burlywood1"))+
    geom_vline(xintercept=c(1.5,2.5,3.5), linetype="dashed")+ # ,2.5,3.5
    xlab("")+
    ylab("Percentile of descriptor value (%)")+
    ylim(c(-5,100))+
    labs(fill="Flow",title = paste0(nam2str(bv)," - ",nam2str(sais)))+
    annotate("text",x=c(0.8,1.2),y=-5,label=c(length(ind.1),length(ind.2)),size=4)
  p
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

# Calcule la direction du vent sur un point de grille (ici Prapoutel)
compute.angle.wind <- function(lon=6,lat=45.25,k=1,rean="ERA5",start="1950-01-01",end="2017-12-31"){
  
  # Import
  uwind <- getdata(k = k,day0 = start,day1 = end,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = F,ssp = NULL,pt.lon = lon,pt.lat = lat,return.lonlat = F,var = "uwind")
  vwind <- getdata(k = k,day0 = start,day1 = end,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = F,ssp = NULL,pt.lon = lon,pt.lat = lat,return.lonlat = F,var = "vwind")
  
  # Angle
  phi <- atan(vwind/uwind) # en radians
  phi <- phi/(2*pi)*360 # en degres
  phi[uwind<0] <- 180 + phi[uwind<0]
  phi[uwind>0 & vwind<0] <- 360 + phi[uwind>0 & vwind<0]
  range(phi)
  
  # Export
  save(phi,file = paste0(get.dirstr(k = k,rean = rean,period = "past"),"compute.angle.wind/angle_wind_k",k,"_lon_",lon,"_lat_",lat,"_",start,"_",end,".Rdata"))
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

# Nbre de jour < ou > a un quantile d'indicateurs entre deux sous-périodes, pour une saison et pour un type de temps
compute.sais.diff.density.subperiod <- function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",ref="all"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  #descr <- c("cel","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  
  for(i in 1:length(descr)){
    des.i <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                            threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
    des[,i] <- des.i
    if(ref=="all") des[,i] <- ecdf(des.i)(des.i)*100
  }
  
  tab <- as.data.frame(cbind(dates,des))
  colnames(tab) <- c("Dates",nam2str(descr,whole=T))
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab[pos.NA.tt,-1] <- NA
  }
  
  if(ref=="wp") tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Saison
  pos <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos.NA.sais <- which(!(1:length(dates)) %in% pos)
  tab[pos.NA.sais,-1] <- NA
  
  if(ref=="wpsais") tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-1,1,4))
  period2 <- paste0(substr(sep+1,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Nettoyage et mise en forme du tableau
  pos <- (1:length(dates))[apply(tab[,2:5],1,function(v){!all(is.na(v))})]
  tab <- tab[pos,]
  
  # Nbre de jours
  qua <- seq(10,50,10)
  res <- data.frame(Percentile=as.factor(x = qua),Per1=NA,Per1.rel=NA,Per2=NA,Per2.rel=NA,Change=NA,Change.rel=NA)
  colnames(res)[c(3,5)] <- c(period1,period2)
  
  len.per1 <- sum(tab$Period==period1) # pour normaliser
  len.per2 <- sum(tab$Period==period2)
  
  for(i in 1:length(qua)){
    per1 <- sum(tab[,2]<qua[i] & tab[,3]<qua[i] & tab[,4]<qua[i] & tab[,5]>100-qua[i] & tab$Period==period1)
    per1.rel <- per1/len.per1*100
    per2 <- sum(tab[,2]<qua[i] & tab[,3]<qua[i] & tab[,4]<qua[i] & tab[,5]>100-qua[i] & tab$Period==period2)
    per2.rel <- per2/len.per2*100
    change <- round(((per2/per1)-1)*100,1)
    change.rel <- round(((per2.rel/per1.rel)-1)*100,1)
    res[i,2:7] <- c(per1,per1.rel,per2,per2.rel,change,change.rel)
    print(paste0("Nbre period 1: ",per1.rel))
    print(paste0("Nbre period 2: ",per2.rel))
    print(paste0("% de variation: ",change.rel,"%"))
  }
  res <- pivot_longer(data = res,cols=c(3,5),values_to="Value",names_to="Period")
  res
}

# Calcule la difference de moyenne des indicateurs entre deux sous-périodes, pour une saison et pour un type de temps
compute.sais.diff.mean.subperiod <- function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",ref="all"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  #descr <- c("cel","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  
  for(i in 1:length(descr)){
    des.i <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                            threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
    des[,i] <- des.i
    if(ref=="all") des[,i] <- ecdf(des.i)(des.i)*100
  }
  
  tab <- as.data.frame(cbind(dates,des))
  colnames(tab) <- c("Dates",nam2str(descr,whole=T))
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab[pos.NA.tt,-1] <- NA
  }
  
  if(ref=="wp") tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Saison
  pos <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos.NA.sais <- which(!(1:length(dates)) %in% pos)
  tab[pos.NA.sais,-1] <- NA
  
  if(ref=="wpsais") tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-1,1,4))
  period2 <- paste0(substr(sep,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Nettoyage et mise en forme du tableau
  pos <- (1:length(dates))[apply(tab[,2:5],1,function(v){!all(is.na(v))})]
  tab <- tab[pos,]
  tab <- pivot_longer(data = tab,cols = 2:5,names_to = "Descr",values_to = "Value")
  tab$Value <- as.numeric(as.character(tab$Value))
  tab$Descr <- factor(tab$Descr,levels = nam2str(descr,whole=T))
  pos.ok <- which(!is.na(tab$Value)) # pour retirer celerite 1er jour
  tab <- tab[pos.ok,]
  
  # Calcul
  print(aggregate(tab$Value,by=list(tab$Period,tab$Descr),mean))
}
  
# Calcul des scores TWS jour à jour entre 2 jeux de données
compute_TWS_crossed <- function(k,start="1851-01-01",end="2010-12-31",rean1,rean2){
  
  N <- length(getdates(start,end))
  
  # Import des gradients des 2 jeux de données
  gradlist<-grad(k,start,end,rean1)
  gradlon1 <- gradlist$gradlon
  gradlat1 <- gradlist$gradlat
  
  gradlist<-grad(k,start,end,rean2)
  gradlon2 <- gradlist$gradlon
  gradlat2 <- gradlist$gradlat
  
  rm(gradlist)
  gc()
  
  # Calcul du TWS
  # Numerateur
  fct <- function(v) sum(abs(v))
  dif_grad <- apply(gradlon1-gradlon2,3,fct)+apply(gradlat1-gradlat2,3,fct) # on calcule la somme des differences entre les gradients de chaque point de grille du jour i et le gradient de chaque point de grille de tous les jours j
    
  # Denominateur
  max_grad <- apply(pmax(abs(gradlon1),abs(gradlon2)),3,sum)+apply(pmax(abs(gradlat1),abs(gradlat2)),3,sum) # on calcule la somme des gradients maximums pour chaque point de grille entre ceux du jour i et ceux de tous les jours j 
  gc()
  
  # TWS
  dist.vec <- dif_grad/max_grad/2
  
  # Export
  save(dist.vec,file = paste0("2_Travail/1_Past/Rresults/compute.TWS.crossed/TWS_",rean1,"_",rean2,"_k",k,"_",start,"_",end,".Rdata"))
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

# Calcule le nombre de jours superieurs/inferieurs à un quantile d'indicateur pour chaque saison d'une periode
get.ind.qua <- function(descr,qua=0.2,lower=T,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                        threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
  des <- ecdf(des)(des)
  if(lower){
    pos <- which(des<qua)
  }else{
    pos <- which(des>(1-qua))
  }
  pos
}

# Calcul une moyenne d'indicateur sur toute la periode et a partir de toutes les reanalyses (pour cel en reference dans plot.TWS.crossed)
get.mean.descr.allrean <- function(k,descr,sais,ana.comm=T){
  
  # Reanalyses
  rean <- c(
    "20CR-m0",
    "20CR-m1",
    "20CR-m2",
    "ERA20C",
    "ERA5"
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
  
  for(i in 1:length(rean)){
    
    # Import de l'indicateur
    des[[i]] <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,
                               start = dates[[i]][1],end = dates[[i]][2],
                               standardize = F,rean = rean[i],threeday = F,
                               period = "past",start.ana = dates.ana[[i]][1],
                               end.ana = dates.ana[[i]][2])
    if(i %in% c(4,5)) des[[i]] <- c(rep(NA,length(des[[i-1]])-length(des[[i]])),des[[i]])
  }
  des <- do.call(cbind,des)
  
  # Traitement saisonnier
  ind <- get.ind.season.past(sais = sais,start = dates[[1]][1],end = dates[[1]][2])
  des <- des[ind,]
  moy <- mean(apply(des,1,mean,na.rm=T),na.rm=T)
  moy
}

# Calcule le nombre de jours superieurs/inferieurs à un quantile d'indicateur pour chaque saison d'une periode
get.nbr.qua <- function(descr,qua=0.2,lower=T,wp=NULL,k,dist,nbdays=1,rean,sais,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                        threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.wp <- which(tt!=wp)
    des[pos.NA.wp] <- NA
  }
  
  # Traitement saisonnier
  pos <- get.ind.season(sais = sais,start = start,end = end)
  des <- des[pos$pos.season]
  des[pos$pos.NA] <- NA
  des.qua <- quantile(des,probs=qua,na.rm=T)
  
  # Nombre de jours
  dim(des) <- c(pos$l.season,pos$n.season)
  count <- apply(des,2,function(v){ifelse(lower,sum(v<des.qua,na.rm=T),sum(v>des.qua,na.rm=T))})
  count
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

# Altitude geopotentiel 1900-1930 et 1970-2000 pour 2 reanalyses et les 4 saisons (pour reviwer #2)
map.abs.geo <- function(k,rean=c("20CR-m1","ERA20C")){
  
  # Import des donnees
  dates.deb <- c("1900-01-01","1929-12-31")
  dates.fin <- c("1970-01-01","1999-12-31")
  season <- c("winter","spring","summer","autumn")
  
  for(i in 1:length(season)){
    
    print(season[i])
    pos.deb <- get.ind.season.past(sais = season[i],start = dates.deb[1],end = dates.deb[2])
    pos.fin <- get.ind.season.past(sais = season[i],start = dates.fin[1],end = dates.fin[2])
    
    # Figure
    png(filename = paste0("2_Travail/1_Past/Rresults/map.abs.geo/map_abs_combine_",season[i],"_k",k,"_",rean[1],"_",rean[2],".png"),width = 6,height = 6,units = "in",res=600)
    layout(matrix(c(1:4,rep(5,2)),nrow = 3,ncol = 2,byrow = T),widths = c(1,1),heights = c(1,1,0.3))
    par(mar=c(1,1,2,1),pty="s")
    data(wrld_simpl)
    
    for(j in 1:length(rean)){
      print(rean[j])
      data.deb <- getdata(k = k,day0 = dates.deb[1],day1 = dates.deb[2],rean = rean[j],climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
      data.fin <- getdata(k = k,day0 = dates.fin[1],day1 = dates.fin[2],rean = rean[j],climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
      
      data.deb.sea <- data.deb[,,pos.deb]
      data.fin.sea <- data.fin[,,pos.fin]
      
      # Moyennes
      mean.deb <- apply(data.deb.sea,1:2,mean)
      mean.fin <- apply(data.fin.sea,1:2,mean)
      
      # Pour ajout du rectangle
      nc <- load.nc(rean[j],var="hgt")
      nc <- nc[[k]]
      lon <- nc$dim$lon$vals
      lat <- nc$dim$lat$vals
      fen <- getinfo_window(k = k,rean=rean[j],var = "hgt")
      delta <- 0.25/2 # on force la grille de ERA5 pour tracer le meme rectangle tout le temps
      nc_close(nc)
      
      # Cartes
      breaks <- seq(4850,6100,length.out = 12)
      N <- 11
      lab <- seq(4900,6100,200)
        
      image(lon,lat,mean.deb,xlim=c(-15,25),ylim=c(25,65),main=paste0(nam2str(rean[j])," - 1900-1930"),
              col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
      plot(wrld_simpl, add = TRUE)
      rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
      
      image(lon,lat,mean.fin,xlim=c(-15,25),ylim=c(25,65),main=paste0(nam2str(rean[j])," - 1970-2000"),
            col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
      plot(wrld_simpl, add = TRUE)
      rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    }
    # Legende
    par(pty="m",mar=c(0,0,0,0))
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
                labels = lab,at =  seq(0,1,length.out = length(lab)),
                vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
    text(x = 1,y = 0.9,"Geopotential Height (m)",cex=1.2,font=2)
    graphics.off()
  }
}

# Altitude geopotentiel 1950-1983 et 1984-2017 pour 2 wp et les 4 saisons (pour moi)
map.abs.geo.wp <- function(wp=1,k,rean,start="1950-01-01",end="2017-12-31",extr=F,bv=NULL){
  
  dates <- getdates(start,end)
  
  # Import des donnees
  print("Import")
  data <- getdata(k = k,day0 = start,day1 = end,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
  gc()
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4))
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  per <- rep(1,length(dates))
  per[as.Date(dates)>=sep] <- 2
  
  # WP
  tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Pour lon/lat et ajout du rectangle
  nc <- load.nc(rean,var="hgt")
  nc <- nc[[k]]
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  fen <- getinfo_window(k = k,rean=rean,var = "hgt")
  delta <- abs(lon[1]-lon[2])/2
  nc_close(nc)
  
  # Figure par saison
  season <- c("winter","spring","summer","autumn")
  sais <- rep(NA,length(dates))
  data(wrld_simpl)
  
  for(i in 1:length(season)){
    print(season[i])
    ind <- get.ind.season.past(sais = season[i],start = start,end = end,nbdays = 1)
    sais[ind] <- season[i]
    
    # Donnees
    if(!extr){
    pos.per1 <- which(per==1 & sais==season[i] & tt==wp)
    pos.per2 <- which(per==2 & sais==season[i] & tt==wp)
    }else{
      ind.max <- get.ind.max.sais(sais = season[i],wp = wp,nbdays = 1,start = start,end = end,bv = bv,spazm = T)
      pos.per1 <- intersect(ind.max,which(per==1))
      pos.per2 <- intersect(ind.max,which(per==2))
    }
    
    data.per1 <- data[,,pos.per1]
    data.per2 <- data[,,pos.per2]
    
    # Moyennes
    mean.per1 <- apply(data.per1,1:2,mean)
    mean.per2 <- apply(data.per2,1:2,mean)
    
    # Cartes
    png(filename = paste0(get.dirstr(k,rean,"past"),"map.abs.geo.wp/map_abs_",season[i],"_k",k,"_wp",wp,ifelse(extr,paste0("_extr_",bv),""),".png"),width = 6,height = 4,units = "in",res=600)
    layout(matrix(c(1:2,rep(3,2)),nrow = 2,ncol = 2,byrow = T),widths = c(1,1),heights = c(1,0.3))
    par(mar=c(1,1,2,1),pty="s")
    
    breaks <- seq(4850,6100,length.out = 12)
    N <- 11
    lab <- seq(4900,6100,200)
    
    image(lon,lat,mean.per1,xlim=c(-15,25),ylim=c(25,65),main=period1,
          col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
    plot(wrld_simpl, add = TRUE)
    rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    
    image(lon,lat,mean.per2,xlim=c(-15,25),ylim=c(25,65),main=period2,
          col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
    plot(wrld_simpl, add = TRUE)
    rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    
    # Legende
    par(pty="m",mar=c(0,0,0,0))
    plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
    colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
                labels = lab,at =  seq(0,1,length.out = length(lab)),
                vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
    text(x = 1,y = 0.9,"Geopotential Height (m)",cex=1.2,font=2)
    graphics.off()
  }
}

# Difference altitude geopotentiel 1900-1930 et 1970-2000 pour 2 reanalyses et les 4 saisons
map.diff.geo <- function(k,rean=c("20CR-m1","ERA20C"),signif=F,save=T){
  
  # Import des donnees
  dates.deb <- c("1900-01-01","1929-12-31")
  dates.fin <- c("1970-01-01","1999-12-31")
  season <- c("winter","spring","summer","autumn")
  
  # Figure
  if(save){
  png(filename = paste0("2_Travail/1_Past/Rresults/map.diff.geo/map_diff_combine_allsais_k",k,"_",rean[1],"_",rean[2],".png"),width = 8,height = 5,units = "in",res=600)
  layout(matrix(c(1:8,rep(9,4)),nrow = 3,ncol = 4,byrow = T),widths = c(1,1,1,1),heights = c(1,1,0.3))
  par(mar=c(1,1,2,1),pty="s")
  }
  data(wrld_simpl)
  
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
    delta <- 0.25/2 # on force la grille de ERA5 pour tracer le meme rectangle tout le temps
    nc_close(nc)
    
    for(j in 1:length(season)){
      print(season[j])
      pos.deb <- get.ind.season.past(sais = season[j],start = dates.deb[1],end = dates.deb[2])
      pos.fin <- get.ind.season.past(sais = season[j],start = dates.fin[1],end = dates.fin[2])
      
      data.deb.sea <- data.deb[,,pos.deb]
      data.fin.sea <- data.fin[,,pos.fin]
      
      # Test de significativite
      if(signif){
        selec <- ifelse(substr(rean[i],1,4)=="20CR",1,2) # on ne prend que 1 pt sur 2 pr 20CR et 1 sur 4 pour ERA20C
        ind.lon <- seq(1,length(lon),by=selec) 
        ind.lat <- seq(2,length(lat),by=selec)
        lon.red <- lon[ind.lon]
        lat.red <- lat[ind.lat]
        
        sig <- matrix(data = NA,nrow = length(lon.red),ncol = length(lat.red))
        for(l in 1:length(lon.red)){
          for(m in 1:length(lat.red)){
            sig[l,m] <- as.numeric(t.test(data.deb.sea[ind.lon[l],ind.lat[m],],data.fin.sea[ind.lon[l],ind.lat[m],])$p.value < 0.05)
          }
        }
        sig.res <- expand.grid(lon.red,lat.red)
        colnames(sig.res) <- c("lon","lat")
        dim(sig) <- nrow(sig.res)
        sig.res <- cbind(sig.res,sig)
      }
      
      # Moyennes et difference
      mean.deb <- apply(data.deb.sea,1:2,mean)
      mean.fin <- apply(data.fin.sea,1:2,mean)
      diff.geo <- mean.fin - mean.deb
      
      # Cartes
      breaks <- seq(-75,75,length.out = 12)
      N <- 11
      lab <- seq(-75,75,25)
      
      image(lon,lat,diff.geo,xlim=c(-15,25),ylim=c(25,65),main=nam2str(season[j]),
                 col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
      
      plot(wrld_simpl, add = TRUE)
      if(signif){points(sig.res$lon,sig.res$lat,cex=sig.res$sig/6,pch=20)} # differences de moyennes significatives, en petits points
      rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    }
  }
  
  # Legende
  if(save){
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = lab,at =  seq(0,1,length.out = length(lab)),
              vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
  text(x = 1,y = 0.9,"1970-2000 minus 1900-1930 Geopotential Height (m)",cex=1.2,font=2)
  
  graphics.off()
  }
}

# Difference altitude geopotentiel 1950-1983 et 1984-2017 pour 2 wp et les 4 saisons
map.diff.geo.wp <- function(wp=1,k,rean,start="1950-01-01",end="2017-12-31",save=T,signif=F){
  
  dates <- getdates(start,end)
  
  # Import des donnees
  print("Import")
  data <- getdata(k = k,day0 = start,day1 = end,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
  gc()
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4))
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  per <- rep(1,length(dates))
  per[as.Date(dates)>=sep] <- 2
  
  # Saison
  season <- c("winter","spring","summer","autumn")
  sais <- rep(NA,length(dates))
  for(i in 1:length(season)){
    ind <- get.ind.season.past(sais = season[i],start = start,end = end,nbdays = 1)
    sais[ind] <- season[i]
  }
  
  # WP
  tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Figure
  # Pour lon/lat et ajout du rectangle
  nc <- load.nc(rean,var="hgt")
  nc <- nc[[k]]
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  fen <- getinfo_window(k = k,rean=rean,var = "hgt")
  delta <- abs(lon[1]-lon[2])/2
  nc_close(nc)
  
  if(save) {
    png(filename = paste0(get.dirstr(k,rean,"past"),"map.diff.geo.wp/map_diff_k",k,"_wp",wp,".png"),width = 8,height = 3,units = "in",res=600)
    layout(matrix(c(1:4,rep(5,4)),nrow = 2,ncol = 4,byrow = T),widths = c(1,1,1,1),heights = c(1,0.3))
    par(mar=c(1,1,2,1),pty="s")
  }
  data(wrld_simpl)
  
  for(i in 1:length(season)){
    print(season[i])
    
    # Donnees
    pos.per1 <- which(per==1 & sais==season[i] & tt==wp)
    pos.per2 <- which(per==2 & sais==season[i] & tt==wp)
    
    data.per1 <- data[,,pos.per1]
    data.per2 <- data[,,pos.per2]
    
    # Test de significativite
    if(signif){
      ind.lon <- seq(1,length(lon),by=8) # on ne prend que 1 pt sur 8 pour ERA5
      ind.lat <- seq(4,length(lat),by=8)
      lon.red <- lon[ind.lon]
      lat.red <- lat[ind.lat]
      
      sig <- matrix(data = NA,nrow = length(lon.red),ncol = length(lat.red))
      for(j in 1:length(lon.red)){
        for(k in 1:length(lat.red)){
          sig[j,k] <- as.numeric(t.test(data.per1[ind.lon[j],ind.lat[k],],data.per2[ind.lon[j],ind.lat[k],])$p.value < 0.05)
        }
      }
      sig.res <- expand.grid(lon.red,lat.red)
      colnames(sig.res) <- c("lon","lat")
      dim(sig) <- nrow(sig.res)
      sig.res <- cbind(sig.res,sig)
    }
    
    # Moyennes et difference
    mean.per1 <- apply(data.per1,1:2,mean)
    mean.per2 <- apply(data.per2,1:2,mean)
    diff.per <- mean.per2 - mean.per1
    
    # Cartes
    breaks <- seq(-60,60,length.out = 12)
    N <- 11
    lab <- seq(-60,60,15)
    
    image(lon,lat,diff.per,xlim=c(-15,25),ylim=c(25,65),main=nam2str(season[i]),
          col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
    plot(wrld_simpl, add = TRUE)
    if(signif){points(sig.res$lon,sig.res$lat,cex=sig.res$sig/6,pch=20)} # differences de moyennes significatives, en petits points
    rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    gc()
  }
  
  # Legende
  if(save){
  par(pty="m",mar=c(0,0,0,0))
  plot(1,1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",ylim=c(0,1))
  colorlegend(colbar = rev(brewer.pal(n = N, name = "RdBu")),
              labels = lab,at =  seq(0,1,length.out = length(lab)),
              vertical = F,xlim = c(0.8,1.2),ylim = c(0.25,0.65),cex=1.4)
  text(x = 1,y = 0.9,paste0(period2," minus ",period1," Geopotential Height (m)"),cex=1.2,font=2)
  graphics.off()
  }
}

# Difference altitude geopotentiel 1900-1930 et 1970-2000 pour 2 reanalyses et les 4 saisons
map.diff.rean <- function(k,rean=c("20CR-m1","ERA20C_regrid_20CR"),start="1900-01-01",end="1929-12-31",signif=F){
  
  # Import des donnees
  data.rean1 <- getdata(k = k,day0 = start,day1 = end,rean = rean[1],climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
  data.rean2 <- getdata(k = k,day0 = start,day1 = end,rean = rean[2],climat = NULL,run = 1,large_win = F,small_win = F,all = T,ssp = NULL,var = "hgt")
  
  nc <- load.nc(rean[1],var="hgt")
  nc <- nc[[k]]
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  fen <- getinfo_window(k = k,rean=rean[1],var = "hgt")
  delta <- 0.25/2 # on force la grille de ERA5 pour tracer le meme rectangle tout le temps
  nc_close(nc)
  
  # Cartes
  season <- c("winter","spring","summer","autumn")
  data(wrld_simpl)
  
  for(i in 1:length(season)){
    
    print(season[i])
    pos.sais <- get.ind.season.past(sais = season[i],start = start,end = end)
    data.rean1.sea <- data.rean1[,,pos.sais]
    data.rean2.sea <- data.rean2[,,pos.sais]
      
      # Test de significativite
      if(signif){
        selec <- ifelse(substr(rean[1],1,4)=="20CR",1,2) # on ne prend que 1 pt sur 2 pr 20CR et 1 sur 4 pour ERA20C
        ind.lon <- seq(1,length(lon),by=selec) 
        ind.lat <- seq(2,length(lat),by=selec)
        lon.red <- lon[ind.lon]
        lat.red <- lat[ind.lat]
        
        sig <- matrix(data = NA,nrow = length(lon.red),ncol = length(lat.red))
        for(l in 1:length(lon.red)){
          for(m in 1:length(lat.red)){
            sig[l,m] <- as.numeric(t.test(data.rean1.sea[ind.lon[l],ind.lat[m],],data.rean2.sea[ind.lon[l],ind.lat[m],])$p.value < 0.05)
          }
        }
        sig.res <- expand.grid(lon.red,lat.red)
        colnames(sig.res) <- c("lon","lat")
        dim(sig) <- nrow(sig.res)
        sig.res <- cbind(sig.res,sig)
      }
      
      # Moyennes et difference
      mean.rean1 <- apply(data.rean1.sea,1:2,mean)
      mean.rean2 <- apply(data.rean2.sea,1:2,mean)
      diff.geo <- mean.rean1 - mean.rean2
      
      # Cartes
      breaks <- seq(-75,75,length.out = 12)
      N <- 11
      lab <- seq(-75,75,25)
      
      image(lon,lat,diff.geo,xlim=c(-15,25),ylim=c(25,65),main=nam2str(season[i]),
            col=rev(brewer.pal(n = N, name = "RdBu")),xaxt="n",yaxt="n",xlab="",ylab="",breaks = breaks)
      
      plot(wrld_simpl, add = TRUE)
      if(signif){points(sig.res$lon,sig.res$lat,cex=sig.res$sig/6,pch=20)} # differences de moyennes significatives, en petits points
      rect(xleft = lon[fen[1,1]]-delta,ybottom = lat[fen[2,1]]-delta,xright = lon[fen[1,1]+fen[1,2]-1]+delta,ytop = lat[fen[2,1]+fen[2,2]-1]+delta,lwd=2)
    }
}

# Graphique des directions de vent pour deux sous periodes
plot.angle.wind <- function(wp="all",nbdays=1,lon=6,lat=45.25,k=1,rean="ERA5",start="1950-01-01",end="2017-12-31",extr=F,bv="Isere-seul"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import
  load(paste0(get.dirstr(k = k,rean = rean,period = "past"),"compute.angle.wind/angle_wind_k",k,"_lon_",lon,"_lat_",lat,"_",start,"_",end,".Rdata"))
  
  if(nbdays>1){ # petite manip pour avoir des moyennes coherentes (5° et 355° ne donnent pas 180° mais 0°)
  phi.tab <- matrix(data = NA,nrow = length(dates),ncol = nbdays)
  for(i in 1:nbdays){phi.tab[,i] <- phi[i:(length(phi)-nbdays+i)]}
  phi <- apply(phi.tab,1,mean)
  pos <- which(apply(phi.tab,1,function(v){if(diff(range(v))>180){return(1)}else{return(0)}})==1)
  phi[pos] <- phi[pos]+180 # si ecart de plus de 180°, on veut l'angle oppose
  phi[phi>360] <- phi[phi>360]-360
  }
  
  # Traitement par classe
  classes <- seq(0,330,30)
  phi.cla <- phi
  for(i in 1:length(phi)){phi.cla[i] <- which.min(abs(phi[i]-classes))}
  
  # WP
  ind.wp <- 1:length(phi)
  if(wp!="all"){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    ind.wp <- which(tt == wp)
  }
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  ind.per1 <- which(as.Date(dates)<as.Date(sep))
  ind.per2 <- which(as.Date(dates)>=as.Date(sep))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4)) # - 10 et +10 jours pour etre sur de sauter une annee entre end period1 et start period 2
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  
  # Graphique par saison
  season <- c("winter","spring","summer","autumn")
  colo <- c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3])
  
  png(filename = paste0(get.dirstr(k = k,rean = rean,period = "past"),"plot.angle.wind/angle_wind_k",k,"_lon_",lon,"_lat_",lat,"_wp",wp,"_",nbdays,"day_",ifelse(extr,paste0("max_sais_",bv,"_"),""),start,"_",end,".png"),width = 6,height = 5,units = "in",res = 600)
  par(mfrow=c(2,2),mar=c(0,0,2,1))
  
  for(i in 1:length(season)){
    print(season[i])
    
    ind.sea <- get.ind.season.past(sais = season[i],start = start,end = end,nbdays = nbdays)
    ind.sea.wp <- intersect(ind.sea,ind.wp)
    if(extr){
      ind.extr <- get.ind.max.sais(sais = season[i],wp = wp,nbdays = nbdays,start = start,end = end,bv = bv)
      ind.sea.wp <- intersect(ind.sea.wp,ind.extr)
    }
    
    data <- as.data.frame(matrix(data = 0,nrow = 4,ncol = length(classes)))
    count.per1 <- table(phi.cla[intersect(ind.per1,ind.sea.wp)])
    count.per2 <- table(phi.cla[intersect(ind.per2,ind.sea.wp)])
    data[3,as.numeric(names(count.per1))] <- count.per1
    data[4,as.numeric(names(count.per2))] <- count.per2
    data[1,] <- rep(max(data,ncol(data)))
    data[2,] <- rep(0,ncol(data))
    rownames(data)[3:4] <- c(period1,period2)
    colnames(data) <- classes
    data <- cbind(data[,-(1:3)],data[,1:3]) # manip pour avoir rose des sables dans le bon sens
    
    radarchart(df = data,pcol = colo,pfcol=alpha(colo,0.3),plwd=3,plty=1,
               cglty=1, cglcol="grey", cglwd=0.8,title=nam2str(season[i]))
    legend(x=0.7, y=1.4, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colo, cex=0.8, pt.cex=2)
  }
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

# Densites des indicateur pour deux sous-periodes par saison et pour un type de temps
plot.descr.density.subperiod <- function(descr,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",add.nb.density=F){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                        threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
  tab <- data.frame(Dates=dates,Descr=des,Season=NA,Period=NA)
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab$Descr[pos.NA.tt] <- NA
  }
  
  # Saison
  sais <- c("winter","spring","summer","autumn")
  for(i in 1:length(sais)){
    pos <- get.ind.season.past(sais = sais[i],start = start,end = end,nbdays = nbdays)
    tab$Season[pos] <- nam2str(sais[i])
  }
  tab$Season <- factor(tab$Season,levels = nam2str(sais))
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-1,1,4))
  period2 <- paste0(substr(sep,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Graphiques Densité
  pl <- ggplot(tab, aes(x=Descr, fill=Period)) +
    theme(plot.margin = unit(c(0.5,0.5,1,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12,colour="black",vjust=0),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=4,face="bold",size=14),
          legend.position = "bottom",legend.title = element_text(hjust=-1,vjust=0.6,size = 16,face = "bold"),
          legend.text = element_text(size=12,colour="black"))+
    geom_density(alpha=0.4)+
    facet_wrap(~Season) +
  xlab(nam2str(descr,whole=T)) + 
    ylab("Density")+
    ggtitle(nam2str(descr,whole=T))+
    scale_fill_manual(values = c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]))
    
  # Pour eviter chevauchement x-ticks dans le combine
  if(descr=="cel") pl <- pl+scale_x_continuous(breaks = seq(from = 0.1,to = 0.4,by = 0.1)) 
  if(descr=="sing05") pl <- pl+scale_x_continuous(breaks = seq(from = 0.2,to = 0.4,by = 0.1)) 
  if(descr=="rsing05") pl <- pl+scale_x_continuous(breaks = c(0.9,0.95))
  
  # Extraction et ajout du nombre de points dans les densites
 if(add.nb.density){
  pl.data <- ggplot_build(pl)$data[[1]]
  n <- aggregate(pl.data$n,by=list(pl.data$group,pl.data$PANEL),unique)[,3]
  dim(n) <- c(2,length(sais))
  n[1,] <- paste0("n_per1=",n[1,]);n[2,] <- paste0("n_per2=",n[2,])
  n <- apply(n,2,function(v){paste0(v[1],"\n",v[2])})
  
  dat_text <- data.frame(
    label = n,
    Season   = nam2str(sais),
    Period = "1950-1983"
  )
  pl <- pl + geom_text(
    data    = dat_text,
    mapping = aes(x = Inf, y = Inf, label = label),
    hjust   = 1,
    vjust   = 1,
    size = 3
  )
 }
  
  # Calcul et tracé de la significativite de la difference de distribution: Kolmogorov-Smirnov Test
  test <- vector(length=length(sais))
  
  for(i in 1:length(sais)){
    x <- tab$Descr[tab$Season==nam2str(sais[i]) & tab$Period==period1]
    y <- tab$Descr[tab$Season==nam2str(sais[i]) & tab$Period==period2]
    test[i] <- ks.test(x = x,y = y)$p.value < 0.05
  }
  
  pl+geom_rect(data=subset(tab,Season %in% nam2str(sais)[test]),aes(fill=Season),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,colour="red",fill=NA,alpha=0.9)
}

# Boxplot des indicateurs pour les journees humides (>1mm) ou seches sur un BV, pour une saison et un wp
plot.descr.rain.norain <- function(bv,wp,sais,k,dist,nbdays=1,rean,spazm,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",save=T){
  
  dates <- getdates(start,end)
  descr <- c("cel","sing05","rsing05","dP")
  
  # Imports
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  for(i in 1:length(descr)){
  des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,
                            rean = rean,threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
  }
  colnames(des) <- nam2str(descr,whole = T)
  tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = bv,agreg = T,spazm = spazm)
  
  # Traitement
  pos.wp <- which(tt==wp)
  des[!(1:length(dates)) %in% pos.wp,] <- NA
  des <- apply(des,2,function(v){ecdf(v)(v)*100}) # percentile par rapport aux journees du meme wp
  
  rain.norain <- as.numeric(precip>1)
  rain.norain[rain.norain==1] <- "Wet";rain.norain[rain.norain==0] <- "Dry"
  des <- cbind(as.data.frame(des),rain.norain)
  
  pos.sea <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos <- intersect(pos.sea,pos.wp)
  des <- des[pos,]
  
  des <- pivot_longer(data = as.data.frame(des),1:4,names_to = "descr",values_to = "val")
  des$rain.norain <- factor(des$rain.norain)
  des$descr <- factor(des$descr,levels = nam2str(descr,whole = T))
  
  # Graphique
  pl <- ggplot(des, aes(x=descr, y=val, fill=rain.norain)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,0,0,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "right",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=16,colour="black"),
          legend.title = element_blank())+
    stat_boxplot(geom = "errorbar",col="darkblue",position=position_dodge(width = 0.75),width=0.3) +
    geom_boxplot(outlier.shape = NA,col="darkblue")+
    scale_fill_manual(values=c(brewer.pal(n = 11, name = "RdBu")[7],brewer.pal(n = 11, name = "RdBu")[9]))+
    ylab("Percentile of descriptor value (%)")
  if(save){
  ggsave(filename = paste0(get.dirstr(k,rean,period = "past"),"plot.descr.rain.norain/plot_descr_",bv,"_",sais,"_wp",wp,"_",nbdays,"days",ifelse(spazm,"_spazm",""),"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),plot = pl,width = 8,height = 5)
  graphics.off()
  }else{pl}
}

# Graphique du nombre d'indicateurs en dessous d'un seuil sur une saison pour deux sous-periodes et pour un type de temps (comptage d'extremes de precipitation dans ces seuils)
plot.season.count.subperiod <-  function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",add.nb.box=F){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  
  for(i in 1:length(descr)){
    des[,i] <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                            threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
    #des[,i] <- ecdf(des.i)(des.i)*100
  }
  
  tab <- as.data.frame(cbind(dates,des))
  colnames(tab) <- c("Dates",nam2str(descr,whole=T))
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab[pos.NA.tt,-1] <- NA
  }
  
  pos <- apply(tab[,2:5],1,function(v){!all(is.na(v))})
  tab[,2:5] <- apply(tab[,2:5],2,function(v){v <- na.omit(as.numeric(as.character(v)));dis <- ecdf(v)(v)*100;res <- rep(NA,nrow(tab));res[pos] <- dis;res})
  
  # Saison
  pos <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos.NA.sais <- which(!(1:length(dates)) %in% pos)
  tab[pos.NA.sais,-1] <- NA
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-1,1,4))
  period2 <- paste0(substr(sep,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Nettoyage et mise en forme du tableau
  tab <- pivot_longer(data = tab,cols = 2:5,names_to = "Descr",values_to = "Value")
  tab$Value <- as.numeric(as.character(tab$Value))
  tab$Descr <- factor(tab$Descr,levels = nam2str(descr,whole=T))
  
  # Precipitations extremes
  extr <- get.ind.max.flow(flow = wp,agreg = T,nbdays = nbdays,start = start,end = end,spazm = T,supseuil = T,nei = T)
  dates[extr]
  tmp <- tab[tab$Dates %in% dates[extr],]
  
  pos <- which(tab$Descr=="Celerity" & tab$Value<25)
  pos1 <- which(tab$Descr=="Singularity" & tab$Value<25)
  pos2 <- which(tab$Descr=="Relative Singularity" & tab$Value<25)
  pos3 <- which(tab$Descr=="MPD" & tab$Value>75)
  pos <- sort(c(pos,pos1,pos2,pos3))
  tab <- tab[pos,]
  
  # Graphiques Violon + boxplot
  dodge <- position_dodge(width = 0.9)
  
  pl <- ggplot(tab, aes(x=Descr, y=Value, fill=Period)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0,0,-0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "bottom",legend.key.size = unit(1.5,"cm"),
          legend.title = element_text(hjust=0.5,vjust=0.5,size = 18,face = "bold"),
          legend.text = element_text(size=16,colour="black"))+
    geom_bar(stat="identity")+
    scale_fill_manual(values = alpha(c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]),0.4))+
    #ylim(-5,100)+
    xlab("")+
    ylab("")+
    ggtitle(nam2str(sais,whole=T))
  
  # Extraction et ajout du nombre de points dans les densites
  if(add.nb.box){
    pl.data <- ggplot_build(pl)$data[[1]]
    n <- aggregate(pl.data$n,by=list(pl.data$group,pl.data$PANEL),unique)[,3][1:2]
    
    pl <- pl + annotate("text",x=c(0.8,1.2),y=-5,label=as.character(n),size=5)
  }
  
  # Calcul et tracé de la significativite de la difference de distribution: Kolmogorov-Smirnov Test et Anderson-Darling Test
  test.ks <- test.ad <- vector(length=length(descr))
  
  for(i in 1:length(descr)){
    x <- tab$Value[tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period1]
    y <- tab$Value[tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period2]
    test.ks[i] <- ks.test(x = x,y = y)$p.value < 0.05
    test.ad[i] <- ad.test(x,y)$ad[1,3] < 0.05
  }
  
  res <- test.ks+test.ad
  pos <- 1:4;pos <- pos[res>=1]
  for(i in pos){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="red",fill=NA,size=1,linetype=ifelse(res[i]==2,1,2))}
  #pos.ks <- pos.ad <- 1:4
  #pos.ks <- pos.ks[test.ks];pos.ad <- pos.ad[test.ad]
  #for(i in pos.ks){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="red",fill=NA,size=1)}
  #for(i in pos.ad){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="purple",fill=NA,size=1,linetype=3)}
  
  pl
}

# Singularité en fonction de la latitude du min de Z500, par WP
plot.sing.lat <- function(rean="ERA5",start="1950-01-01",end="2019-12-31"){
  
  # Import des donnees
  sing <- get.descriptor(descriptor = "sing05",k = 1,dist = "TWS",nbdays = 1,start = start,end = end,standardize = F,rean = rean,threeday = F,desais = F,period = "past",start.ana = "1950-01-01",end.ana = "2010-12-31")
  data <- getdata(k = 1,day0 = start,day1 = end,rean = rean,climat = NULL,run = 1,large_win = F,small_win = F,all = F,ssp = NULL,return.lonlat = T,var = "hgt")
  lat <- data$lat
  data <- data$data
  wp <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  wp.unique <- sort(unique(wp))
  
  # Traitement
  pos.min <- apply(data,3,function(mat){lat[which(mat == min(mat), arr.ind = TRUE)[2]]})
  
  # Graphique
  namflow <- c("Atlantic","Mediterranean","North-East","Anticyclonic")
  png(filename = paste0(get.dirstr(k = 1,rean = rean,period = "past"),"plot.sing.lat/plot_sing_lat_wp_",start,"_",end,".png"),width = 6,height = 6,units = "in",res = 600)
  par(pty="s",mfrow=c(2,2),mar=c(4,4,2,0))
  for(i in 1:length(wp.unique)){
    pos.i <- which(wp==wp.unique[i])
    corr <- round(cor(sing[pos.i],pos.min[pos.i],use = "pairwise.complete.obs"),2)
    plot(sing[pos.i],pos.min[pos.i],ylab="min Z500 latitude (°)",xlab="sing",main=paste0(namflow[i]," (R=",corr,")"),pch=19,cex=0.5)
    grid();par(new=T)
    plot(sing[pos.i],pos.min[pos.i],ylab="min Z500 latitude (°)",xlab="sing",main=paste0(namflow[i]," (R=",corr,")"),pch=19,cex=0.5)
    abline(lm(pos.min[pos.i]~sing[pos.i]),col="red",lwd=2)
  }
  graphics.off()
}

# Boxplot d'un indicateur pour deux sous-periodes par saison et pour un type de temps
plot.descr.violin.subperiod <-  function(descr,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",add.nb.box=F){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                        threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
  tab <- data.frame(Dates=dates,Descr=ecdf(des)(des)*100,Season=NA,Period=NA)
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab$Descr[pos.NA.tt] <- NA
  }
  
  # Saison
  sais <- c("winter","spring","summer","autumn")
  for(i in 1:length(sais)){
    pos <- get.ind.season.past(sais = sais[i],start = start,end = end,nbdays = nbdays)
    tab$Season[pos] <- nam2str(sais[i])
  }
  tab$Season <- factor(tab$Season,levels = nam2str(sais))
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-1,1,4))
  period2 <- paste0(substr(sep,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Graphiques Violon + boxplot
  dodge <- position_dodge(width = 0.9)
  
  pl <- ggplot(tab, aes(x=Season, y=Descr, fill=Period)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0,0,-0.5,-0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=12,colour="black",vjust=0),
          axis.text.y = element_text(size=12),plot.title = element_text(hjust = 0.5,vjust=0,face="bold",size=12),
          legend.position = "right",legend.key.size = unit(1.5,"cm"),legend.text = element_text(size=12),
          legend.title = element_blank())+
    geom_violin(position=dodge,col="black")+
    geom_boxplot(position=dodge,width=0.2,alpha=0)+
    scale_fill_manual(values = alpha(c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]),0.4))+
    ylim(-5,100)+
    xlab("")+
    ylab("")+
    ggtitle(nam2str(descr,whole=T))
  
  # Extraction et ajout du nombre de points dans les densites
  if(add.nb.box){
    pl.data <- ggplot_build(pl)$data[[1]]
    n <- aggregate(pl.data$n,by=list(pl.data$group,pl.data$PANEL),unique)[,3]
    
    pl <- pl + annotate("text",x=c(0.8,1.2,1.8,2.2,2.8,3.2,3.8,4.2),y=-5,label=as.character(n),size=4)
  }
  
  # Calcul et tracé de la significativite de la difference de distribution: Kolmogorov-Smirnov Test
  test <- vector(length=length(sais))
  
  for(i in 1:length(sais)){
    x <- tab$Descr[tab$Season==nam2str(sais[i]) & tab$Period==period1]
    y <- tab$Descr[tab$Season==nam2str(sais[i]) & tab$Period==period2]
    test[i] <- ks.test(x = x,y = y)$p.value < 0.05
  }
  
  pos <- 1:4; pos <- pos[test]
  for(i in pos){
  pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="red",fill=NA,size=1)
  }
  
  pl
}

# Boxplot des indicateurs d'une saison pour deux sous-periodes par indicateur et pour un type de temps
plot.sais.violin.subperiod <-  function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",add.nb.box=F){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  
  for(i in 1:length(descr)){
    des.i <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                              threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
    des[,i] <- ecdf(des.i)(des.i)*100
  }
  
  tab <- as.data.frame(cbind(dates,des))
  colnames(tab) <- c("Dates",nam2str(descr,whole=T))
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab[pos.NA.tt,-1] <- NA
  }
  
  # Saison
  pos <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos.NA.sais <- which(!(1:length(dates)) %in% pos)
  tab[pos.NA.sais,-1] <- NA
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-1,1,4))
  period2 <- paste0(substr(sep,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Nettoyage et mise en forme du tableau
  pos <- (1:length(dates))[apply(tab[,2:5],1,function(v){!all(is.na(v))})]
  tab <- tab[pos,]
  tab <- pivot_longer(data = tab,cols = 2:5,names_to = "Descr",values_to = "Value")
  tab$Value <- as.numeric(as.character(tab$Value))
  tab$Descr <- factor(tab$Descr,levels = nam2str(descr,whole=T))
  
  # Graphiques Violon + boxplot
  dodge <- position_dodge(width = 0.9)
  
  pl <- ggplot(tab, aes(x=Descr, y=Value, fill=Period)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0,0,-0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "bottom",legend.key.size = unit(1.5,"cm"),
          legend.title = element_text(hjust=0.5,vjust=0.5,size = 18,face = "bold"),
          legend.text = element_text(size=16,colour="black"))+
    geom_violin(position=dodge,col="black")+
    geom_boxplot(position=dodge,width=0.2,alpha=0)+
    scale_fill_manual(values = c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]))+
    ylim(-5,100)+
    xlab("")+
    ylab("")+
    ggtitle(nam2str(sais,whole=T))
  
  # Extraction et ajout du nombre de points dans les densites
  if(add.nb.box){
    pl.data <- ggplot_build(pl)$data[[1]]
    n <- aggregate(pl.data$n,by=list(pl.data$group,pl.data$PANEL),unique)[,3][1:2]
    
    pl <- pl + annotate("text",x=c(0.8,1.2),y=-5,label=as.character(n),size=5)
  }
  
  # Calcul et tracé de la significativite de la difference de distribution: Kolmogorov-Smirnov Test et Anderson-Darling Test
  test.ks <- test.ad <- vector(length=length(descr))
  
  for(i in 1:length(descr)){
    x <- tab$Value[tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period1]
    y <- tab$Value[tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period2]
    test.ks[i] <- ks.test(x = x,y = y)$p.value < 0.05
    test.ad[i] <- ad.test(x,y)$ad[1,3] < 0.05
  }
  
  res <- test.ks+test.ad
  pos <- 1:4;pos <- pos[res>=1]
  for(i in pos){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="red",fill=NA,size=1,linetype=ifelse(res[i]==2,1,2))}
  #pos.ks <- pos.ad <- 1:4
  #pos.ks <- pos.ks[test.ks];pos.ad <- pos.ad[test.ad]
  #for(i in pos.ks){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="red",fill=NA,size=1)}
  #for(i in pos.ad){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=-1,ymax=101,colour="purple",fill=NA,size=1,linetype=3)}
  
  pl
}

# Differences de distribution des indicateurs entre deux sous-périodes, pour une saison et pour un type de temps
plot.sais.diff.density.subperiod <- function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",add.nb.box=F,ref="all"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  descr <- c("cel","sing05","rsing05","dP")
  #descr <- c("cel","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  
  for(i in 1:length(descr)){
    des.i <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                            threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
    des[,i] <- des.i
    if(ref=="all") des[,i] <- ecdf(des.i)(des.i)*100
  }
  
  tab <- as.data.frame(cbind(dates,des))
  colnames(tab) <- c("Dates",nam2str(descr,whole=T))
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = ifelse(end>as.Date("2017-12-31"),F,T))
    pos.NA.tt <- which(tt!=wp)
    tab[pos.NA.tt,-1] <- NA
  }
  
  if(ref=="wp") tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Saison
  pos <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos.NA.sais <- which(!(1:length(dates)) %in% pos)
  tab[pos.NA.sais,-1] <- NA
  
  if(ref=="wpsais") tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4)) # - 10 et +10 jours pour etre sur de sauter une annee entre end period1 et start period 2
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Nettoyage et mise en forme du tableau
  pos <- (1:length(dates))[apply(tab[,2:5],1,function(v){!all(is.na(v))})]
  #pos <- (1:length(dates))[apply(tab[,2:3],1,function(v){!all(is.na(v))})]
  tab <- tab[pos,]
  tab <- pivot_longer(data = tab,cols = 2:5,names_to = "Descr",values_to = "Value")
  #tab <- pivot_longer(data = tab,cols = 2:3,names_to = "Descr",values_to = "Value")
  tab$Value <- as.numeric(as.character(tab$Value))
  tab$Descr <- factor(tab$Descr,levels = nam2str(descr,whole=T))
  pos.ok <- which(!is.na(tab$Value)) # pour retirer celerite 1er jour
  tab <- tab[pos.ok,]
  
  # Graphiques boxplot
  dodge <- position_dodge(width = 0.9)
  tab.segment <- data.frame(x = c(1,2,3,4), y = rep(0,4), xend = c(1,2,3,4), yend = rep(100,4),Period=NA,Descr=NA,Value=NA)
  #tab.segment <- data.frame(x = c(1,2), y = rep(0,4), xend = c(1,2), yend = rep(100,4),Period=NA,Descr=NA,Value=NA)
  
  
  pl <- ggplot(tab, aes(x=Descr, y=Value, fill=Period)) + 
    theme_bw()+
    theme(plot.margin = unit(c(0.5,0,0,-0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "bottom",legend.key.size = unit(1.5,"cm"),
          legend.title = element_text(hjust=0.5,vjust=0.5,size = 18,face = "bold"),
          legend.text = element_text(size=16,colour="black"))+
    geom_boxplot(position=dodge,width=0.2)+
    scale_fill_manual(values = c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]))+
    ylim(-5,100)+
    xlab("")+
    ylab("")+
    ggtitle(nam2str(sais,whole=T))+
    geom_segment(data=tab.segment,aes(x=x,y=y,xend=xend,yend=yend), linetype="dashed")
  
  # Ajout densites
  value <- as.numeric(as.matrix(tab$Value))
  res <- list()
  
  for(i in 1:length(descr)){
    # Difference de densite
    pos1 <- which(tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period1)
    pos2 <- which(tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period2)
    dens1 <- density(value[pos1],from=0,to=100,n=101)
    dens2 <- density(value[pos2],from=0,to=100,n=101)
    x <- c(0,dens2$y-dens1$y,0)
    y <- c(0,dens1$x,100)
    
    # Mise en forme
    x <- x*35; x <- x+i
    tmp <- data.frame(X=x,Y=y,Descr=rep(nam2str(descr[i],whole=T),length(x)),Value=rep(NA,length(x))) # ajout des colonnes factices Descr et Value obligatoire
    
    # Graphique densites
    pl <- pl +
      geom_path(data = tmp, aes(x = X, y = Y,fill=NULL))
    
    # Stockage des positions pour remplissage
    Xmin <- tmp$X;Xmin[which(tmp$X>i)] <- i
    Xmax <- tmp$X;Xmax[which(tmp$X<i)] <- i
    res[[i]] <- data.frame(Xmin=Xmin,Xmax=Xmax)
    #res[[i]] <- data.frame(Xmin.prec=rep(0.5,length(x)),Xmin=Xmin,Xmax=Xmax,Y=tmp$Y,Descr=tmp$Descr,Value=tmp$Value)
    #if(i>1) res[[i]][,"Xmin.prec"] <- res[[i-1]][,"Xmax"]
  }
  
  # Graphique remplissage: serait plus concis sous forme de boucle, mais mauvais remplissage en boucle avec geom_ribbon...
  
  #for(i in 1:length(descr)){
  #  pl <- pl + geom_ribbon(data=res[[i]],aes(xmin=Xmin.prec,xmax=Xmin,y=Y),fill=NA)
  #  pl <- pl + geom_ribbon(data=res[[i]],aes(xmin=Xmin,xmax=i,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[9])
  #  pl <- pl + geom_ribbon(data=res[[i]],aes(xmin=i,xmax=Xmax,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[3])
  #}
  
  # Mise en forme
  res <- do.call(cbind,res)
  res <- cbind(tmp$Descr,tmp$Value,tmp$Y,res)
  colnames(res)[1:3] <- c("Descr","Value","Y")
  colnames(res)[6:7] <- paste0(colnames(res)[6:7],"2")
  colnames(res)[8:9] <- paste0(colnames(res)[8:9],"3")
  colnames(res)[10:11] <- paste0(colnames(res)[10:11],"4")
  
  # Graphique
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmin,xmax=1,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[9])
  pl <- pl + geom_ribbon(data=res,aes(xmin=1,xmax=Xmax,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[3])
 
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmax,xmax=Xmin2,y=Y),fill=NA)
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmin2,xmax=2,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[9])
  pl <- pl + geom_ribbon(data=res,aes(xmin=2,xmax=Xmax2,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[3])
  
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmax2,xmax=Xmin3,y=Y),fill=NA)
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmin3,xmax=3,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[9])
  pl <- pl + geom_ribbon(data=res,aes(xmin=3,xmax=Xmax3,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[3])
  
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmax3,xmax=Xmin4,y=Y),fill=NA)
  pl <- pl + geom_ribbon(data=res,aes(xmin=Xmin4,xmax=4,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[9])
  pl <- pl + geom_ribbon(data=res,aes(xmin=4,xmax=Xmax4,y=Y),fill=brewer.pal(n = 11, name = "RdBu")[3])

  # Extraction et ajout du nombre de points dans les densites
  if(add.nb.box){
    n <- c(sum(tab$Descr==nam2str(descr[1],whole=T) & tab$Period==period1),
           sum(tab$Descr==nam2str(descr[1],whole=T) & tab$Period==period2))
    
    pl <- pl + annotate("text",x=c(0.8,1.2),y=-5,label=as.character(n),size=5)
  }
  
  # Calcul et tracé de la significativite de la difference de distribution: Kolmogorov-Smirnov Test et Anderson-Darling Test
  test.ks <- test.ad <- vector(length=length(descr))
  
  for(i in 1:length(descr)){
    x <- tab$Value[tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period1]
    y <- tab$Value[tab$Descr==nam2str(descr[i],whole=T) & tab$Period==period2]
    test.ks[i] <- ks.test(x = x,y = y)$p.value < 0.05
    test.ad[i] <- ad.test(x,y)$ad[1,3] < 0.05
  }
  
  res <- test.ks+test.ad
  pos <- 1:4;pos <- pos[res>=1]
  for(i in pos){pl <- pl+geom_rect(xmin=i-0.5,xmax=i+0.5,ymin=0,ymax=100,colour="red",fill=NA,size=1,linetype=ifelse(res[i]==2,1,2))}
  
  pl
}

# Graphique Nbre de jour < ou > a un quantile d'indicateurs entre deux sous-périodes, pour une saison et pour un type de temps
plot.sais.diff.nbdays.subperiod <- function(sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",ref="all"){
  
  # Calcul des nombres de jours
  tab <- compute.sais.diff.density.subperiod(sais=sais,wp=wp,k=k,dist=dist,nbdays=nbdays,rean=rean,start=start,end=end,start.ana=start.ana,end.ana=end.ana,ref=ref)
  
  # Graphique
  #if(nbdays==1){
  #  ylab <- Number of days"
  #  main <- "1 day"
  #}else{
  #  ylab <- paste0("Number of ",nbdays,"-day sequences")
  #  main <- paste0(nbdays," days")
  #}
  
  main <- nam2str(sais)
  ylab <- "% of days"
  ylimit <- c(0,30)
  
  ggplot(data = tab,aes(Percentile,Value, fill=Period))+
    theme_bw()+
    theme(plot.margin = unit(c(1,0.5,0.5,0.5),"cm"),axis.title.x = element_text(vjust=-2,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=3,size = 12,face = "bold"),axis.title.y.right = element_text(vjust=3,size = 12,face = "bold",color = "darkgreen"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),axis.text.y.right = element_text(size=13,color="darkgreen"),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "bottom",legend.key.size = unit(1,"cm"),
          legend.title = element_text(hjust=0.5,vjust=0.5,size = 18,face = "bold"),
          legend.text = element_text(size=16,colour="black"))+
    geom_col(position="dodge")+
    scale_y_continuous(limits=ylimit,sec.axis = sec_axis(~ . *5, name = "% of change"))+
    geom_line(aes(x = c(1,1,2,2,3,3,4,4,5,5),y=Change.rel/5),color="darkgreen",size=1.5)+
    geom_text(aes(label = round(Value,1)),position = position_dodge(0.9),vjust = -0.5,fontface="bold")+
    xlab("Percentile of descriptor value (%)")+
    ylab(ylab)+
    ggtitle(main)+
    scale_fill_manual(values = c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]))
}

# Trace l'evolution du nombre de jours en dessous/dessus d'un quantile pour tous les indicateurs de maniere groupee, pour toutes les ciculations, pour une saison
plot.trend.alldescr.nbr.qua <- function(qua=0.25,k,dist,nbdays=1,rean,sais,wp="all",start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",save=T){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  year <- unique(substr(dates,1,4))
  if(sais=="winter") year <- year[-1]
  
  # Import des indices des jours sup/inf le quantile
  descr <- c("cel","sing05","rsing05","dP")
  inf <- c(T,T,T,F)
  quant <- c(qua,qua,qua,1-qua)
  ind.qua <- vector(mode="list",length=length(descr))
  
  for(i in 1:length(descr)){
  ind.qua[[i]] <- get.ind.qua(descr = descr[i],qua = qua,lower = inf[i],k = k,dist = dist,nbdays = nbdays,rean = rean,start = start,end = end,start.ana = start.ana,end.ana = end.ana)
  }
  ind.qua.all <- Reduce(intersect,ind.qua)
  
  # Cb de precip extreme dans ces caracteristiques?
  extr <- get.ind.max.sais(sais = sais,wp = wp,nbdays = nbdays,start = start,end = end,bv = "Drac-seul",spazm = T)
  match.extr <- intersect(extr,ind.qua.all)
  print(paste0("Nombre extremes dans ces caracteristiques: ",length(match.extr),"/",length(extr)," soit ",round(length(match.extr)/length(extr)*100,2),"%"))
  
  # WP
  if(wp!="all"){
  tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = spazm)
  }
  
  # Traitement saisonnier
  vec.qua <- rep(0,length(dates))
  vec.qua[ind.qua.all] <- 1
  
  if(wp!="all"){vec.qua[tt!=wp] <- 0}
  
  pos <- get.ind.season(sais = sais,start = start,end = end)
  vec.qua <- vec.qua[pos$pos.season]
  vec.qua[pos$pos.NA] <- NA
  dim(vec.qua) <- c(pos$l.season,pos$n.season)
  res <- apply(vec.qua,2,sum,na.rm=T)
  
  # Graphique
  if(save){
  png(filename = paste0(get.dirstr(k = k,rean = rean,period = "past"),"plot.trend.alldescr.nbr.qua/plot_",sais,"_q",qua,"_wp",wp,"_",nbdays,"days_",rean,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),width = 6,height = 4,units = "in",res = 600)
  par(mar=c(2,4,2,0))
  }
  
  plot(year,res,type="n",ylab="Number of days",main=nam2str(sais),cex.main=1.5)
  grid()
  abline(h=0,lty=3)
  par(new=T)
  plot(year,res,type="l",ylab="Number of days",main=nam2str(sais),cex.main=1.5)
  
  reg <- glm(res~as.numeric(year),family = "poisson")
  pval <- summary(reg)$coefficients[2,4]
  lty <- ifelse(pval<0.05,1,2)
  lines(year,reg$fitted.values,col="blue",lwd=2,lty=lty)
  text(2010,max(res),paste0("pval = ",round(pval,2)),font=2,col="blue")
  if(save){graphics.off()}
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
    
    # Lissage et alignement
    if(liss!=1){
      des[[i]][,2] <- rollapply(des[[i]][,2],liss,mean,partial=F,fill=NA)
    }
    
    if(align){
      delta[i] <- des[[i]][max(which(!is.na(des[[i]][,2]))),2]
      des[[i]][,2] <- des[[i]][,2] - delta[i]
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
  #main <- ifelse(type=="season",nam2str(sais),nam2str(descr,whole=T))
  main <- ""
  
  ylim <- range(unlist(lapply(des,function(v){range(v[,2],na.rm=T)})))
  if(align) {delta <- delta - delta[which(rean=="ERA5")]; ylim <- range(ylim,delta)} # ERA5 en reference
  
  if(liss==1 & align){
  if(descr=="cel") ylim <- c(-0.1,0.06)
  if(descr=="sing05") ylim <- c(-0.06,0.08)
  if(descr=="rsing05") ylim <- c(-0.02,0.015)
  if(descr=="dP") ylim <- c(-150,100)
  }
  
  if(liss==5 & align){
  if(descr=="cel") ylim <- c(-0.08,0.02)
  if(descr=="sing05") ylim <- c(-0.03,0.03)
  if(descr=="rsing05") ylim <- c(-0.01,0.015)
  if(descr=="dP") ylim <- c(-80,40)
  }
  
  # Initialisation
  if(save) {
    png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.descr/plot_trend_",descr,"_",sais,"_liss=",liss,ifelse(ana.comm,paste0("_ana_",dates.ana[[1]][1],"_",dates.ana[[1]][2]),""),ifelse(align,"_align",""),ifelse(nao,"_nao",""),".png"),width = 7,height = 5,units = "in",res=600)
    par(mar=c(4,4,2,1))
  }
  plot(des[[1]][,2],type="n",xaxt="n",yaxt="n",ylim=ylim,xlab="Year",ylab=nam2str(descr,whole=T,unit = T),main=main)
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = pos.axis,labels = xaxis,cex.axis=1.5)
  axis(side = 2,cex.axis=1.5)
  abline(v = pos.axis,lty=3,col="grey")
  
  # Ajout courbes
  for(i in 1:length(rean)){
    lines(match(des[[i]][,1],des[[1]][,1]),des[[i]][,2],col=colo[i],lwd=2)
  }
  
  # Points delta si align
  if(align) points(rep(tail(pos.axis,1),length(delta)),delta,col=colo,pch=18,cex=2)
  
  # NAO
  if(nao){
    par(new=T)
    plot(nao.ind,pch=19,col="grey",type="l",lwd=2,xlab="",ylab="",xaxt="n",yaxt="n")
    abline(h=0,col="grey")
    axis(side = 4)
    mtext("NAOI",side=4,line=2)
  }
  
  # Legende
  if(leg) legend("bottomright",inset=.01,nam2str(rean),col=colo,lty=1,lwd=2,bty="n",cex=1.5)
  if(save) graphics.off()
}

# Trace l'evolution d'un indicateur conditionnellement au type de temps et a la saison
plot.trend.descr.cond <- function(descr,k,dist,rean,sais="all",type="mean",start="1950-01-01",end="2010-12-31",start.ana="1950-01-01",end.ana="2010-12-31",pvalue=T){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
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
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
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
  
  png(filename = paste0(get.dirstr(k,rean,"past"),"plot.trend.descr.cond.regquant/plot_",descr,"_",sais,"_",start,"_",end,"_ana_",start.ana,"_",end.ana,".png"),width = 8,height = 6,units = "in",res = 600)
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
    if(ncol(summa[[1]]$coefficients) == 4){
    abline(coef[,1],col="red",lwd=2,lty=ifelse(summa[[1]]$coefficients[2,4]<0.05,1,2))
    abline(coef[,2],col="blue",lwd=2,lty=ifelse(summa[[2]]$coefficients[2,4]<0.05,1,2))
    abline(coef[,3],col="red",lwd=2,lty=ifelse(summa[[3]]$coefficients[2,4]<0.05,1,2))
    }
    
    # Les pvalue
    if(ncol(summa[[1]]$coefficients) == 4){
      shadowtext(x = round(quantile(all.year,0.9),0),y = tail(fitted[,1],1)-0.05*diff(ylim),labels = paste0("pval = ",round(summa[[1]]$coefficients[2,4],4)),font = 2,cex = 0.8,col = "red",bg = "white",r = 0.2)
      shadowtext(x = round(quantile(all.year,0.9),0),y = tail(fitted[,2],1)-0.05*diff(ylim),labels = paste0("pval = ",round(summa[[2]]$coefficients[2,4],4)),font = 2,cex = 0.8,col = "blue",bg = "white",r = 0.2)
      shadowtext(x = round(quantile(all.year,0.9),0),y = tail(fitted[,3],1)-0.05*diff(ylim),labels = paste0("pval = ",round(summa[[3]]$coefficients[2,4],4)),font = 2,cex = 0.8,col = "red",bg = "white",r = 0.2)
    }
  }
  
  graphics.off()
}

# Trace l'evolution des indicateurs aux dates des precipitations extremes pour l'annee et les 4 saisons, par BV ou par WP
plot.trend.descr.extr <- function(bv="Isere-seul",wp=1,descr,k,dist,rean,nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  dat <- as.Date(dates)
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = start,end = end,
                        standardize = F,rean = rean,threeday = F,period = "past",start.ana = start.ana,
                        end.ana = end.ana)
  
  # Graphique
  sais <- c("year","winter","spring","summer","autumn")
  
  png(filename = paste0(get.dirstr(k,rean,period="past"),"plot.trend.descr.extr/plot_",bv,"_wp",wp,"_mean",nbdays,"day_",descr,"_",rean,"_",start,"_",end,".png"),width = 7,height = 7,units = "in",res=600)
  par(mfrow=c(3,2),mar=c(2,4,2,1))
  
  for(i in 1:length(sais)){
    sea <- unique(get.ind.season(sais = sais[i],start = start,end = end)$pos.season)
    extr <- get.ind.max.sais(sais = sais[i],wp = wp,nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
    pos <- intersect(sea,extr)
    plot(dat[pos],des[pos],pch=19,xlab="Year",ylab=nam2str(descr,whole=T),main=nam2str(sais[i]),xlim=range(dat))
    grid();par(new=T)
    plot(dat[pos],des[pos],pch=19,xlab="",ylab="",xlim=range(dat))
    reg <- lm(des[pos]~dat[pos])
    abline(reg,col="red")
    text(x = quantile(as.numeric(dat[pos]),0.5),y=max(des[pos]),paste0("pvalue=",round(summary(reg)$coefficients[,4][2],3)),col="red",font=2)
  }
  graphics.off()
}

# Trace l'evolution des indicateurs aux dates des precipitations extremes pour l'annee pour un bv, une saison, mais tous les inicateurs en percentiles
plot.trend.descr.extr.alldescr <- function(bv="Isere-seul",sais,k,dist,rean,nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",leg=T){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  ann <- as.numeric(unique(substr(dates,1,4)))
  
  # Import de l'indicateur
  descr <- c("cel","sing05","rsing05","dP")
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  for(i in 1:length(descr)){
    des.i <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,
                        standardize = F,rean = rean,threeday = F,period = "past",start.ana = start.ana,
                        end.ana = end.ana)
    des[,i] <- ecdf(des.i)(des.i)*100
  }
  
  # Graphique
  colo <- colo <- c("dimgrey","red","darkgreen","blue")
  extr <- get.ind.max.sais(sais = sais,wp = "all",nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  if(sais=="winter"){extr <- c(NA,extr)}
  
  plot(ann,des[extr,1],col=colo[1],type="l",xlab="Year",ylab="Percentile of descriptor value (%)",main=paste0(nam2str(bv)," - ",nam2str(sais)),cex.main=1.3,ylim=c(0,115),font.lab=2)
  grid()
  abline(h=c(0,100),lty=3)
  for(i in 1:length(descr)){
    lines(ann,des[extr,i],col=colo[i])
    reg <- lm(des[extr,i]~ann)
    lty.i <- ifelse(summary(reg)$coefficients[,4][2]<0.05,1,2)
    abline(reg,col=colo[i],lty=lty.i,lwd=2)
  }
  if(leg){legend("top",nam2str(descr,whole = T),col=colo,lty=rep(1,4),lwd=2,bty="n",ncol=2,cex=0.9)}
}

# Trace l'evolution du nombre de jours en dessous/dessus d'un quantile pour tous les indicateurs, pour toutes les ciculations puis circulations Atl et Med, pour une saison
plot.trend.descr.nbr.qua <- function(qua=0.25,lower=T,wp1=1,wp2=2,k,dist,nbdays=1,rean,sais,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31",leg=T,save=T){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  year <- unique(substr(dates,1,4))
  if(sais=="winter") year <- year[-1]
  
  # Import des nombres de jours sup/inf le quantile
  descr <- c("cel","sing05","rsing05","dP")
  inf <- c(T,T,T,F)
  quant <- c(qua,qua,qua,1-qua)
  nbr <- nbr.wp1 <- nbr.wp2 <- matrix(NA,nrow = length(year),ncol = length(descr))
  
  for(i in 1:length(descr)){
  nbr[,i] <- get.nbr.qua(descr = descr[i],qua = quant[i],lower = inf[i],wp = NULL,k = k,dist = dist,nbdays = nbdays,rean = rean,sais = sais,
                     start = start,end = end,start.ana = start.ana,end.ana = end.ana)
  nbr.wp1[,i] <- get.nbr.qua(descr = descr[i],qua = quant[i],lower = inf[i],wp = wp1,k = k,dist = dist,nbdays = nbdays,rean = rean,sais = sais,
                         start = start,end = end,start.ana = start.ana,end.ana = end.ana)
  nbr.wp2[,i] <- get.nbr.qua(descr = descr[i],qua = quant[i],lower = inf[i],wp = wp2,k = k,dist = dist,nbdays = nbdays,rean = rean,sais = sais,
                         start = start,end = end,start.ana = start.ana,end.ana = end.ana)
  }
  
  # Graphique
  start.xaxis <- trunc(as.numeric(year[1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- round(as.numeric(year[length(year)])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  
  if(save) png(filename = paste0(get.dirstr(k,rean,period="past"),"plot.trend.descr.sais/plot_",descr,"_",rean,"_mean",nbdays,"day_",start,"_",end,".png"),width = 7,height = 5,units = "in",res=600)
  par(mar=c(4,4,2,1))
  
  # All
  ran <- range(nbr);ylim <- c(ran[1],ran[1]+diff(ran)*1.2);rm(ran)
  plot(year,nbr[,1],type="n",xaxt="n",xlab="Year",ylab="Number of day < q0.2",main=paste0(nam2str(sais)," - All"),ylim=ylim)
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
  abline(v = xaxis,lty=3,col="grey")
  
  for(i in 1:4){
    lines(year,nbr[,i],col="gray40",lwd=2,lty=i)
  }
  
  # wp1
  ran <- range(nbr.wp1);ylim <- c(ran[1],ran[1]+diff(ran)*1.2);rm(ran)
  plot(year,nbr.wp1[,1],type="n",xaxt="n",xlab="Year",ylab="Number of day < q0.2",main=paste0(nam2str(sais)," - WP1"),ylim=ylim)
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
  abline(v = xaxis,lty=3,col="grey")
  
  for(i in 1:4){
    lines(year,nbr.wp1[,i],col="lightblue",lwd=2,lty=i)
  }
  
  # wp2
  ran <- range(nbr.wp2);ylim <- c(ran[1],ran[1]+diff(ran)*1.2);rm(ran)
  plot(year,nbr.wp2[,1],type="n",xaxt="n",xlab="Year",ylab="Number of day < q0.2",main=paste0(nam2str(sais)," - WP2"),ylim=ylim)
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
  abline(v = xaxis,lty=3,col="grey")
  
  for(i in 1:4){
    lines(year,nbr.wp2[,i],col="red",lwd=2,lty=i)
    
  }
  
  tmp <- lm(nbr.wp1[,4]~as.numeric(year))
  c(tmp$coefficients,round(unname(summary(tmp)$coefficients[,4][2]),3))
  
}

# Trace l'evolution d'un indicateur conditionnellement a la saison
plot.trend.descr.sais <- function(descr,k,dist,rean,nbdays=1,start="1950-01-01",end="2010-12-31",start.ana="1950-01-01",end.ana="2010-12-31",leg=T,save=T){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  year <- unique(substr(dates,1,4))
  
  # Import de l'indicateur
  des <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                        threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
  
  # Traitement saisonnier
  sais <- c("winter","spring","summer","autumn")
  des.sais <- matrix(NA,nrow = length(year),ncol = length(sais))
  reg <- vector(mode = "list",length = 4)
  
  for(i in 1:length(sais)){
    # Moyenne saisonniere
    pos <- get.ind.season(sais = sais[i],start = start,end = end)
    tmp <- des[pos$pos.season]
    tmp[pos$pos.NA] <- NA
    dim(tmp) <- c(pos$l.season,pos$n.season)
    tmp <- apply(tmp,2,mean,na.rm=T)
    if(sais[i]=="winter") tmp <- c(NA,tmp)
    des.sais[,i] <- tmp
    rm(tmp)
    
    # Regression
    tmp <- lm(des.sais[,i]~as.numeric(year))
    reg[[i]] <- c(tmp$coefficients,round(unname(summary(tmp)$coefficients[,4][2]),3)) # intersept, slope, pvalue
    gc()
  }
  
  # Graphique
  colo <- c("blue","darkgreen","red","darkorange")
  start.xaxis <- trunc(as.numeric(year[1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- round(as.numeric(year[length(year)])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  ran <- range(des.sais,na.rm=T);ylim <- c(ran[1],ran[1]+diff(ran)*1.2);rm(ran)
  
if(save) png(filename = paste0(get.dirstr(k,rean,period="past"),"plot.trend.descr.sais/plot_",descr,"_",rean,"_mean",nbdays,"day_",start,"_",end,".png"),width = 7,height = 5,units = "in",res=600)
    par(mar=c(4,4,2,1))
    plot(year,des.sais[,1],type="n",xaxt="n",xlab="Year",ylab=nam2str(descr,whole=T),main=nam2str(descr,whole=T),ylim=ylim)
    grid(ny=NULL,nx=NA)
    axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
    abline(v = xaxis,lty=3,col="grey")
    
    for(i in 1:4){
      lines(year,des.sais[,i],col=colo[i],lwd=2)
      abline(reg[[i]][1],reg[[i]][2],col=colo[i],lwd=2,lty=ifelse(reg[[i]][3]<0.05,1,3))
    }
    
    if(leg) legend("top",nam2str(sais),col=colo,lty=1,lwd=2,bty="n",cex=0.9,ncol=2)
    if(save) graphics.off()
}

# Graphique d'evolution des max saisonniers de precip avec droites GEV (sans enregistrement)
plot.trend.precip.gev <- function(bv="Isere-seul",nbdays=1,spazm=T,start="1950-01-01",end="2017-12-31",leg=T){
  
  dates <- getdates(start,end)
  ann <- as.numeric(unique(substr(dates,1,4)))
  
  # Imports
  precip <- get.precip(nbdays = nbdays,start = start,end = end,bv = bv,spazm = spazm)
  load(paste0("2_Travail/1_Past/Rresults/fit.gev.obs/fit_gev_",bv,"_mean",nbdays,"day",ifelse(spazm,"_spazm",""),"_",start,"_",end,".Rdata"))
  
  # Traitement et Graphique
  sais <- c("winter","spring","summer","autumn")
  meanGEV.fct<-function(mu,sig,xi){mu+sig/xi*(gamma(1-xi)-1)}
  
  for(i in 1:length(sais)){
    
    # Traitement max saisonnier
    pos <- get.ind.season(sais = sais[i],start = start,end = end)
    precip.i <- precip[pos$pos.season]
    precip.i[pos$pos.NA] <- NA
    dim(precip.i) <- c(pos$l.season,pos$n.season)
    max.vec <- apply(precip.i,2,function(v){max(v,na.rm=T)})
    if(sais[i]=="winter"){max.vec <- c(NA,max.vec)}
    
    # Graphique max saisonnier
    plot(ann,max.vec,type="l",xlab="",ylab="",ylim=c(0,100),main=nam2str(sais[i]))
    grid();par(new=T)
    plot(ann,max.vec,type="l",xlab="",ylab="",ylim=c(0,100),main=nam2str(sais[i]))
    
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
    #lines(ann,rl20,col="red",lty=lty.i,lwd=2)
    #if(i==1 & leg){legend("topright",c("RL20","Mean"),lty=c(1,1),lwd=2,col=c("red","blue"),bty="n")}
    if(i==1 & leg){legend("topright",c("Mean"),lty=1,lwd=2,col="blue",bty="n")}
  }
}

# Trace l'evolution des cumuls saisonniers issus des differentes influences atmospheriques
plot.trend.precip.wp <- function(bv="Isere-seul",wp=1,spazm=T,start="1950-01-01",end="2017-12-31"){
  
  dates <- getdates(start,end)
  year <- unique(substr(dates,1,4))
  
  # Import precip
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = spazm)
  
  # Import wp
  tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = bv,agreg = T,spazm = spazm)
  precip[tt!=wp] <- NA
  
  # Traitement saisonnier
  sais <- c("winter","spring","summer","autumn")
  precip.sais <- matrix(NA,nrow = length(year),ncol = length(sais))
  reg <- vector(mode = "list",length = 4)
  
  for(i in 1:length(sais)){
    # Moyenne saisonniere
    pos <- get.ind.season(sais = sais[i],start = start,end = end)
    tmp <- precip[pos$pos.season]
    tmp[pos$pos.NA] <- NA
    dim(tmp) <- c(pos$l.season,pos$n.season)
    tmp <- apply(tmp,2,sum,na.rm=T)
    if(sais[i]=="winter") tmp <- c(NA,tmp)
    precip.sais[,i] <- tmp
    rm(tmp)
    
    # Regression
    tmp <- lm(precip.sais[,i]~as.numeric(year))
    reg[[i]] <- c(tmp$coefficients,round(unname(summary(tmp)$coefficients[,4][2]),3)) # intersept, slope, pvalue
    gc()
  }
  
  # Graphiques
  wp.num <- c(1,2,5,8)
  wp.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip.wp/plot_trend_precip_",bv,"_wp",wp,"_",start,"_",end,ifelse(spazm,"_spazm",""),".png"),width = 9,height = 5,units = "in",res=600)
  par(mfrow=c(2,2),mar=c(4,4,2,1))
  for(i in 1:4){
  plot(year,precip.sais[,i],type="n",xlab="Year",ylab="Precipitation (mm)",main=paste0(nam2str(sais[i])," - pvalue=",reg[[i]][3]))
  grid()
  lines(year,precip.sais[,i],lwd=2)
  abline(reg[[i]][1],reg[[i]][2],lwd=2,col="red")
  }
  graphics.off()
}

# Trace l'evolution des cumuls saisonniers issus des differentes influences atmospheriques, dans un seul graphique
plot.trend.precip.wp.combine <- function(bv="Isere-seul",spazm=T,start="1950-01-01",end="2017-12-31",norm=T,leg=T,save=T){
  
  dates <- getdates(start,end)
  pos.dec <- which(substr(dates,6,7)=="12")
  dates[pos.dec] <- paste0(as.character(as.numeric(substr(dates[pos.dec],1,4))+1),substr(dates[pos.dec],5,10)) # decembre = annee suivante
  ann <- as.numeric(substr(dates,1,4))
  
  # Import precip
  precip <- get.precip(nbdays = 1,start = start,end = end,bv = bv,spazm = spazm)
  
  # Import wp
  tt.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  tt <- get.wp(nbdays = 1,start = start,end = end,risk = F,bv = bv,agreg = T,spazm = spazm)
  
  # Saison
  sais <- c("winter","spring","summer","autumn")
  vec.sais <- rep(NA,length(dates))
  l.sais <- rep(NA,4)
  
  for(i in 1:length(sais)){
    pos <- get.ind.season(sais = sais[i],start = start,end = end)
    vec.sais[pos$pos.season] <- sais[i]
    l.sais[i] <- pos$l.season
  }
  
  # Aggregation et mise en forme
  res <- aggregate(precip,by=list(vec.sais,ann,tt),sum,na.rm=T)
  colnames(res) <- c("Season","Year","WP","Value")
  if(norm){ # on normalise le cumul de chaque annee par l'occurence de ce wp cette annee là
  occ <- aggregate(tt,by=list(vec.sais,ann,tt),table)
  res$Value <- res$Value/occ$x
  }
  res <- pivot_wider(res,names_from = 3,values_from = 4)
  res[,-c(1,2)] <- apply(res[,-c(1,2)],2,function(v){v[is.na(v)] <- 0;v}) # on met cumul à 0 si combinaison saison/WP manquante
  
  Tot <- aggregate(precip,by=list(vec.sais,ann),sum,na.rm=T)
  l.sais <- l.sais[match(Tot$Group.1,sais)]
  Tot <- Tot$x/l.sais
  res <- cbind(res,Tot)
  
  # Traitement et graphique
  colo <- c("blue","red","darkgreen","dimgrey","black")
  if(save){
    png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.precip.wp.combine/plot_trend_precip_",bv,"_",start,"_",end,ifelse(spazm,"_spazm",""),ifelse(norm,"_norm",""),".png"),width = 9,height = 6,units = "in",res=600)
    par(mfrow=c(2,2),mar=c(2,4,2,1))
  }
  
  for(i in 1:length(sais)){
    pos <- which(res$Season==sais[i])
    res.i <- data.matrix(res[pos,-1])
    if(i==1) res.i <- rbind(c(1950,rep(NA,5)),res.i) # pour avoir meme xaxis que les autres saisons
    res.i.liss <- apply(res.i,2,function(v){unname(rollapply(v,5,mean))})
    
    plot(res.i.liss[,1],res.i.liss[,6],type="n",ylim=c(0,ifelse(norm,10,max(res.i.liss[,-1],na.rm=T))),ylab="Precipitation (mm)",main=nam2str(sais[i]),cex.main=1.5)
    grid()
    
    for(j in c(2,3,6)){
      # Graphiques
      lines(res.i.liss[,1],res.i.liss[,j],col=colo[j-1],lwd=2)
      # Regression sur le non lisse
      reg <- lm(res.i[,j]~res.i[,1])
      abline(reg,col=colo[j-1],lwd=2,lty=ifelse(summary(reg)$coefficients[,4][2]<0.1,1,2))
    }
    if(i==2 & leg) legend("top",legend = c(tt.name[1:2],"Total"),col = colo[c(1,2,5)],lwd=2,lty=1,bty="n",ncol=2,cex=1)
  }
  if(save){graphics.off()}
}

# Trace l'évolution de l'occurrence des WP par saison
plot.trend.wp <- function(type="occurence",sais,nbdays=1,start="1950-01-01",end="2017-12-31",leg=T,save=T){
  
  # Import
  wp <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
  
  # Traitement
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  year <- unique(substr(dates,1,4))
  if(sais=="winter") year <- year[-1]
  
  pos <- get.ind.season(sais = sais,start = start,end = end)
  wp <- wp[pos$pos.season]
  wp[pos$pos.NA] <- NA
  dim(wp) <- c(pos$l.season,pos$n.season)
  
  occ <- pers <- reg.occ <- reg.pers <- vector(mode = "list",length = 4)
  tt <- c(1,2,5,8)
  tt.name <- c("Atlantic","Mediterranean","Northeast","Anticyclonic")
  
  for(i in 1:4){
    # Occurence
    if(type=="occurence"){
      occ[[i]] <- apply(wp,2,function(v){sum(v==tt[i],na.rm=T)/length(v)*100})
      tmp <- lm(occ[[i]]~as.numeric(year))
      reg.occ[[i]] <- c(tmp$coefficients,round(unname(summary(tmp)$coefficients[,4][2]),3)) # intersept, slope, pvalue
    }
    
    # Persistance
    if(type=="persistence"){
      pers[[i]] <- apply(wp,2,function(v){tmp <- rle(v);mean(tmp$length[tmp$values==tt[i]],na.rm=T)})
      tmp <- lm(pers[[i]]~as.numeric(year))
      reg.pers[[i]] <- c(tmp$coefficients,round(unname(summary(tmp)$coefficients[,4][2]),3)) # intersept, slope, pvalue
    }
  }
  
  # Graphiques
  colo <- c("blue","red","darkgreen","dimgrey")
  start.xaxis <- trunc(as.numeric(year[1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- round(as.numeric(year[length(year)])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  
  # Occurence
  if(type=="occurence"){
    if(save){
      png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.wp/plot_occurence_wp_",sais,"_mean",nbdays,"day_",start,"_",end,".png"),width = 7,height = 5,units = "in",res=600)
      par(mar=c(4,4,2,1))
    }
    plot(year,occ[[1]],type="n",xaxt="n",xlab="Year",ylab="WP Occurence (%)",ylim=c(0,100),main=nam2str(sais))
    grid(ny=NULL,nx=NA)
    axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
    abline(v = xaxis,lty=3,col="grey")
    
    for(i in 1:4){
      lines(year,occ[[i]],col=colo[i],lwd=2)
      abline(reg.occ[[i]][1],reg.occ[[i]][2],col=colo[i],lwd=2,lty=ifelse(reg.occ[[i]][3]<0.05,1,3))
    }
    
    if(leg) legend("top",tt.name,col=colo,lty=1,lwd=2,bty="n",cex=1,ncol=2)
    if(save) graphics.off()
  }
  
  # Persistance
  if(type=="persistence"){
    if(save){
      png(filename = paste0("2_Travail/1_Past/Rresults/plot.trend.wp/plot_persistence_wp_",sais,"_mean",nbdays,"day_",start,"_",end,".png"),width = 7,height = 5,units = "in",res=600)
      par(mar=c(4,4,2,1))
    }
    plot(year,pers[[1]],type="n",xaxt="n",xlab="Year",ylab="WP Persistence (days)",ylim=c(0,11),main=nam2str(sais))
    grid(ny=NULL,nx=NA)
    axis(side = 1,at = xaxis,labels = xaxis) # petite manip pour ajouter 1850 a xaxis
    abline(v = xaxis,lty=3,col="grey")
    
    for(i in 1:4){
      lines(year,pers[[i]],col=colo[i],lwd=2)
      abline(reg.pers[[i]][1],reg.pers[[i]][2],col=colo[i],lwd=2,lty=ifelse(reg.pers[[i]][3]<0.05,1,3))
    }
    
    if(leg) legend("top",tt.name,col=colo,lty=1,lwd=2,bty="n",cex=1,ncol=2)
    if(save) graphics.off()
  }
}

# Calcul des scores TWS jour à jour entre 2 jeux de données
plot.TWS.crossed <- function(k,start="1851-01-01",end="2010-12-31",sais,leg=T,save=F){
  
  # Import des TWS
  rean <- list(
    c("20CR-m0","20CR-m1"),
    c("20CR-m0","20CR-m2"),
    c("20CR-m1","20CR-m2"),
    c("20CR-m0","ERA20C_regrid_20CR"),
    c("20CR-m1","ERA20C_regrid_20CR"),
    c("20CR-m2","ERA20C_regrid_20CR"),
    c("20CR-m0","ERA5_regrid_20CR"),
    c("20CR-m1","ERA5_regrid_20CR"),
    c("20CR-m2","ERA5_regrid_20CR"),
    c("ERA20C","ERA5_regrid_ERA20C")
  )
  
  dates <- list(
    c(start,end),
    c(start,end),
    c(start,end),
    c("1900-01-01",end),
    c("1900-01-01",end),
    c("1900-01-01",end),
    c("1950-01-01",end),
    c("1950-01-01",end),
    c("1950-01-01",end),
    c("1950-01-01",end)
  )
  
  colo <- c("gray0","gray30","gray60","red4","red","red3","steelblue4","steelblue1","steelblue","purple")
  
  dist.list <- vector(mode="list",length(length(rean)))
  
  for(i in 1:length(rean)){
    load(paste0("2_Travail/1_Past/Rresults/compute.TWS.crossed/TWS_",rean[[i]][1],"_",rean[[i]][2],"_k",k,"_",dates[[i]][1],"_",dates[[i]][2],".Rdata"))
    dist.list[[i]] <- dist.vec
  }
  
  # Traitement saisonnier
  for(i in 1:length(rean)){
    tmp <- dist.list[[i]]
    pos <- get.ind.season(sais = sais,start = dates[[i]][1],end = dates[[i]][2])
    tmp <- tmp[pos$pos.season]
    tmp[pos$pos.NA] <- NA
    dim(tmp) <- c(pos$l.season,pos$n.season)
    rm(pos)
    
    tmp <- apply(tmp,2,mean,na.rm=T)
    year <- unique(substr(getdates(dates[[i]][1],dates[[i]][2]),1,4))
    if(sais=="winter") year <- year[-1] # car pas de premiere saison hiver
    dist.list[[i]] <- data.frame(Year=year,TWS=tmp)
    dist.list[[i]][,1] <- as.character(dist.list[[i]][,1])
  }
  
  # Graphique
  # Parametres
  start.xaxis <- trunc(as.numeric(dist.list[[1]][1,1])*0.1)*10 # manip pour avoir une annee "ronde" inferieure (si 1958 -> 1950)
  end.xaxis <- trunc(as.numeric(dist.list[[1]][nrow(dist.list[[1]]),1])*0.1)*10 # manip pour avoir une annee "ronde"
  xaxis <- seq(start.xaxis,end.xaxis,10)
  pos.axis <- match(xaxis,dist.list[[1]][,1])
  if(is.na(pos.axis[1])) pos.axis[1] <- pos.axis[2]-10
  #ylim <- c(0,max(unlist(lapply(dist.list,function(v){max(v[,2],na.rm=T)}))))
  ylim <- c(0,0.4)
  main <- nam2str(sais)
  ylab <- "TWS"
  
  # Initialisation
  if(save) png(filename = paste0("2_Travail/1_Past/Rresults/plot.TWS.crossed/plot_TWS_k",k,"_",sais,".png"),width = 7,height = 5,units = "in",res=600)
  par(mar=c(4,4,2,1))
  plot(dist.list[[1]][,2],type="n",xaxt="n",ylim=ylim,xlab="Year",ylab=ylab,main=main)
  grid(ny=NULL,nx=NA)
  axis(side = 1,at = pos.axis,labels = xaxis)
  abline(v = pos.axis,lty=3,col="grey")
  
  abline(h=0,col="black",lwd=2,lty=2)
  #abline(h=get.mean.descr.allrean(k = k,descr = "cel",sais = sais,ana.comm = T),col="black",lwd=2,lty=2)
  abline(h=0.28,col="black",lwd=2,lty=2) # TWS du 1978-12-12 (figure 1 article 2)
  abline(h=0.18,col="black",lwd=2,lty=2) # TWS du 1978-12-12 (figure 1 article 2)
  
  # Ajout courbes
  for(i in 1:length(rean)){
    lines(match(dist.list[[i]][,1],dist.list[[1]][,1]),dist.list[[i]][,2],col=colo[i],lwd=2)
  }
  
  # Legende
  if(leg) legend("topleft",unlist(lapply(rean,function(v){paste(nam2str(v[1]),nam2str(v[2]),sep = " // ")})),col=colo,lty=1,lwd=2,bty="n",cex=0.65,ncol = 2,text.width = 60)
  if(save) graphics.off()
}

# Run des fonctions
run.article2 <- function(type=1){
  
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
  
  # compute.density
  if(type==2){
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
  
  # plot.trend.descr.cond
  if(type==3){
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
  if(type==4){
    descr <- c("cel","dP","sing05","rsing05")
    sais <- c("winter","autumn","spring","summer")
    
    for(i in 1:length(descr)){
      print(descr[i])
      for(j in 1:length(sais)){
        print(sais[j])
        plot.trend.descr.cond.regquant(descr = descr[i],k = 1,dist = "TWS",rean = "ERA5",sais = sais[j],start = "1950-01-01",end = "2010-12-31")
        plot.trend.descr.cond.regquant(descr = descr[i],k = 1,dist = "TWS",rean = "ERA5",sais = sais[j],start = "1950-01-01",end = "2017-12-31")
      }
    }
  }
  
  # compare.descr.ana
  if(type==5){
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
  if(type==6){
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
  if(type==7){
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
  
  # plot.descr.density.subperiod
  if(type==8){
    sais <- c("winter","spring","summer","autumn")
    descr <- c("cel","sing05","rsing05","dP")
    wp <- c(NULL,1,2)
    
    for(i in 1:length(sais)){
      print(sais[i])
      for(j in 1:length(descr)){
        print(descr[j])
        for(k in 1:length(wp)){
          print(wp[k])
          plot.descr.density.subperiod(descr = descr[i],wp = wp[k],k = 1,dist = "TWS",nbdays = 1,rean = "ERA5",sais = sais[j])
        }
      }
    }
  }
  
  # plot.trend.precip.wp
  if(type==9){
    
    bv <- c("Isere-seul","Drac-seul","Isere")
    wp <- c(1,2,5,8)
    spazm <- c(T,F)
    
    for(i in 1:length(bv)){
      print(bv[i])
      for(j in 1:length(wp)){
        print(wp[j])
        for(k in 1:length(spazm)){
          print(spazm[k])
          plot.trend.precip.wp(bv = bv[i],wp = wp[j],spazm = spazm[k],start = "1950-01-01",end = ifelse(spazm[k],"2017-12-31","2019-12-31"))
        }
      }
    }
  }
  
  # plot.trend.descr.extr
  if(type==10){
    descr <- c("cel","sing05","rsing05","dP")
    bv <- c("Isere-seul","Drac-seul","Isere")
    wp <- "all" #c(1,2)
    nbdays <- 1#c(1,3)
    
    for(i in 1:length(descr)){
      print(descr[i])
      for(j in 1:length(bv)){
        print(bv[j])
        for(k in 1:length(wp)){
          print(wp[k])
          for(l in 1:length(nbdays)){
            print(nbdays[l])
            plot.trend.descr.extr(bv = bv[j],wp = wp[k],descr = descr[i],k = 1,dist = "TWS",rean = "ERA5",nbdays = nbdays[l],spazm = T,
                                  start = "1950-01-01",end = "2017-12-31",start.ana = "1950-01-01",end.ana = "2017-12-31")
          }
        }
      }
    }
  }
  
  if(type==11){
    bv <- c("Isere-seul","Drac-seul","Isere")
    sais <- c("winter","spring","summer","autumn")
    wp <- c(1,2)
    
    for(i in 1:length(bv)){
      print(bv[i])
      for(j in 1:length(sais)){
        print(sais[j])
        for(k in 1:length(wp)){
          print(wp[k])
          plot.descr.rain.norain(bv = bv[i],wp = wp[k],sais = sais[j],k = 1,dist = "TWS",nbdays = 1,
                                 rean = "ERA5",spazm = T,start = "1950-01-01",end = "2017-12-31",
                                 start.ana = "1950-01-01",end.ana = "2010-12-31")
        }
      }
    }
  }
}

# Scatterplot des densites des indicateurs dans un plan 2D pour deux sous periodes, pour une influence et une saison
scatterplot.density.subperiod <- function(descr=c("cel","sing05"),sais,wp=NULL,k,dist,nbdays=1,rean,start="1950-01-01",end="2017-12-31",start.ana="1950-01-01",end.ana="2010-12-31"){
  
  dates <- getdates(start,as.character(as.Date(end)-nbdays+1))
  
  # Import des indicateurs
  des <- matrix(data = NA,nrow = length(dates),ncol = length(descr))
  
  for(i in 1:length(descr)){
    des.i <- get.descriptor(descriptor = descr[i],k = k,dist = dist,nbdays = nbdays,start = start,end = end,standardize = F,rean = rean,
                            threeday = F,desais = F,period = "past",start.ana = start.ana,end.ana = end.ana)
    des[,i] <- des.i
  }
  
  tab <- as.data.frame(cbind(dates,des))
  colnames(tab) <- c("Dates","Descr1","Descr2")
  
  # Import des types de temps
  if(!is.null(wp)){
    tt <- get.wp(nbdays = nbdays,start = start,end = end,risk = F,bv = "Isere",agreg = T,spazm = T)
    pos.NA.tt <- which(tt!=wp)
    tab[pos.NA.tt,-1] <- NA
  }
  
  #tab[,-1] <- apply(tab[,-1],2,function(v){ecdf(as.numeric(v))(v)*100})
  
  # Saison
  pos <- get.ind.season.past(sais = sais,start = start,end = end,nbdays = nbdays)
  pos.NA.sais <- which(!(1:length(dates)) %in% pos)
  tab[pos.NA.sais,-1] <- NA
  
  # Periode
  sep <- as.Date(median(as.numeric(as.Date(dates))))
  period1 <- paste0(substr(start,1,4),"-",substr(sep-10,1,4))
  period2 <- paste0(substr(sep+10,1,4),"-",substr(end,1,4))
  tab$Period[as.Date(dates)<sep] <- period1
  tab$Period[as.Date(dates)>=sep] <- period2
  
  # Nettoyage et mise en forme du tableau
  pos <- (1:length(dates))[apply(tab[,2:3],1,function(v){!all(is.na(v))})]
  tab <- tab[pos,]
  tab$Descr1 <- as.numeric(as.character(tab$Descr1))
  tab$Descr2 <- as.numeric(as.character(tab$Descr2))
  tab$Period <- factor(tab$Period,levels=c(period1,period2))
  
  # Graphique
  pl <- ggplot(tab, aes(x=Descr1, y=Descr2)) + 
    theme_bw()+
    theme(plot.margin = unit(c(1.5,0.5,0.5,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
          axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
          axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
          legend.position = "bottom",legend.key.size = unit(1.5,"cm"),
          legend.title = element_text(hjust=0.5,vjust=0.5,size = 18,face = "bold"),
          legend.text = element_text(size=16,colour="black"))+
    geom_point(aes(color=Period),size=0.6)+
    geom_density_2d(aes(colour=Period),size=1,bins=5,contour_var = "density")+
    scale_color_manual(values = c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]))+
    xlab(nam2str(descr[1],whole=T,unit=T))+
    ylab(nam2str(descr[2],whole=T,unit=T))
  
  # Test densite 2D
  tab.per1 <- as.matrix(tab[tab$Period==period1,2:3])
  tab.per2 <- as.matrix(tab[tab$Period==period2,2:3])
  test <- kde.test(na.omit(tab.per1),na.omit(tab.per2))
  
  pl <- pl+
    geom_text(aes(x=quantile(Descr1,0.97,na.rm=T),y=quantile(Descr2,0.99,na.rm=T),label=paste0("pvalue=",round(test$pvalue,2))))
  
  #pl <- ggplot(tab, aes(x=Descr1, y=Descr2)) + 
  #  theme_bw()+
  #  theme(plot.margin = unit(c(0,0,0.5,0.5),"cm"),axis.title.x = element_text(vjust=-4,size = 12,face = "bold"),
  #        axis.title.y = element_text(vjust=4,size = 12,face = "bold"),axis.text.x = element_text(size=13,colour="black",vjust=0),
  #        axis.text.y = element_text(size=13),plot.title = element_text(hjust = 0.5,vjust=3,face="bold",size=18),
  #        legend.position = "bottom",legend.key.size = unit(1.5,"cm"),
  #        legend.title = element_text(hjust=0.5,vjust=0.5,size = 18,face = "bold"),
  #        legend.text = element_text(size=16,colour="black"))+
  #  stat_density_2d(aes(color=Period))+
  #  geom_point(aes(color=Period))+
  #  scale_color_manual(values = c(brewer.pal(n = 11, name = "RdBu")[9],brewer.pal(n = 11, name = "RdBu")[3]))+
  #  xlab(nam2str(descr[1],whole=T,unit=T))+
  #  ylab(nam2str(descr[2],whole=T,unit=T))
  
  pl
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
