source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Calcul des distances de maniere generique
compute_dist_gen_past<-function(k,dist,start="1950-01-01",end="2011-12-31",rean,period="past"){
  if (k %in% 1:2) {
    if (dist %in% c("TWS","RMSE","RMSE_I","RMSE_II","Mahalanobis")){
      if (dist=="TWS") dist.list<-compute_TWS(k,start,end,rean)
      if (dist=="RMSE") dist.list<-compute_RMSE(k,start,end,rean)
      if (dist=="RMSE_I") dist.list<-compute_RMSE_I(k,start,end,rean)
      if (dist=="RMSE_II") dist.list<-compute_RMSE_II(k,start,end,rean)
      if (dist=="Mahalanobis") dist.list<-compute_mahalanobis(k,start,end,rean)
      save(dist.list,file=paste0(get.dirstr(k,rean,period),"compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
    }
    else{ # si on demande un nTWS, sTWS, nRMSE ou sRMSE, les distances sont normalisees par la moyenne ou l'ecart type des distances
      dist0<-substr(dist,2,nchar(dist))
      load(file=paste0(get.dirstr(k,rean,period),"compute_dist/",dist0,"_member",member,"_Z",Z,"_",start,"_",end,".Rdata"))
      type<-substr(dist,1,1)
      print(type)
      if (type=="n") norm<-mean(unlist(dist.list))
      if (type=="s") norm<-sd(unlist(dist.list))
      print(norm)
      for (i in 1:length(dist.list)) dist.list[[i]]<-dist.list[[i]]/norm
      
      save(dist.list,file=paste0(get.dirstr(k,rean,period),"compute_dist/",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
    }
  }
}

# Calcul des indicateurs
compute_criteria_past<-function(k,dist,start="1950-01-01",end="2011-12-31",update=FALSE,rean,threeday=FALSE,rev=FALSE,period="past"){
  gc()
  
  print(paste0(get.dirstr(k,rean,period),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata"))
  
  if (update) {
    load(file=paste0(file=paste0(get.dirstr(k,rean,period),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
    coln<-colnames(criteria)
  }
  
  dist.vec<-getdist(k,ifelse(rev,"TWS",dist),start,end,rean,threeday,period) # si rev, on prend les nei en TWS
  gc()
  dist.vec<-unlist(dist.vec)
  gc()
  dates<-getdates(start,ifelse(threeday,as.character(as.Date(end)-3+1),end))
  N<-length(dates)
  gc()
  
  U<-c(0,(N-1):1); # U = 0, 22644, 22643, 22642, ...
  sU<-sapply(1:(N-1),function(x) sum(U[1:x])) # on fait la somme de U[1], U[1:2], etc pour obtenir la position de la derniere distance qui separe chaque date
  gc()
  
  if (!update) {
    coln.new<-c("cel","sing05","q05")
  }
  if (update) {
    coln.new<-c("celnei","singnei","rsingnei")
  }
  
  criteria.new<-matrix(NA,ncol=length(coln.new),nrow=N)
  colnames(criteria.new)<-coln.new
  
  for (i in 1:N){
    
    ddi<-getdist4i(i,dist.vec,N,sU)
    di<-ddi
    n<-length(di)
    
    gc()
    
    soso<-sort(di,index.return=TRUE) # classement par plus petit score, et donne les positions
    qi05<-di[soso$ix[(0.005*n)]]     # quantile 0.5%
    idi05<-soso$ix[2:(0.005*n)] # recupere la position des 0.5% les plus proches
    
    tmp<-NULL
    
    for (cc in coln.new){
      # Celerite
      #if (cc=="cel") {if (i==1) tmp<-c(tmp,NA) else tmp<-c(tmp,ddi[i-1])}
      if (cc=="celnei") tmp<-c(tmp,mean(criteria[idi05,"cel"],na.rm=TRUE))
      
      ## Singularite
      #if (cc=="sing05") tmp<-c(tmp,mean(di[idi05]))
      if (cc=="singnei") tmp <- c(tmp,mean(criteria[idi05,"sing05"],na.rm=TRUE))
      
      #Singularite relative
      #if (cc=="q05") tmp<-c(tmp,qi05)
      if (cc=="rsingnei") tmp <- c(tmp,mean(criteria[idi05,"rsing05"],na.rm=TRUE))
      
      # dP
      #if (cc=="dPnei") tmp <- c(tmp,mean(criteria[idi05,"dP"],na.rm=TRUE))
    }
    
    idi05m1<-idi05
    criteria.new[i,]<-tmp
    
    if (i %% 50==0) {print(i)}
    gc()
  }
  
  if (update) {
    criteria<-cbind(criteria,criteria.new)
    colnames(criteria)<-c(coln,coln.new)
  }
  else {
    criteria<-criteria.new
    colnames(criteria)<-coln.new
  }
  
  save(criteria,file=paste0(file=paste0(get.dirstr(k,rean,period),"compute_criteria/",ifelse(threeday,"3day_",""),"criteria_",dist,"_member",member,"_k",k,"_",start,"_",end,".Rdata")))
}

# Calcul de dP
get.dP <- function(k,nbdays,start="1950-01-01",end="2011-12-31",rean){
  geo <- getdata(k = k,day0 = start,day1 = end,rean = rean) 
  des <- apply(geo,3,function(x) max(x)-min(x))
  des <- rollapply(des,nbdays,mean)
}

# Trace l'evolution dans le temps d'un indicateur, pour les trois reanalyses
plot.evolution.descr <- function(descr,k,dist,liss=5){
  
  # Variables
  rean <- c("20CR","ERA20C","NCEP")
  dates <- list(
    c("1851-01-01","2010-12-31"),
    c("1900-01-01","2010-12-31"),
    c("1950-01-01","2010-12-29")
  )
  
  # Import
  des <- list()
  for(i in 1:length(rean)){
    des[[i]] <- get.descriptor(descriptor = descr,k = k,dist = dist,nbdays = 1,start = dates[[i]][1],end = dates[[i]][2],
                               standardize = F,rean = rean[i],threeday = F,period = "past")
    des[[i]] <- rollapply(des[[i]],liss*365,mean,partial=T)
  }
  
  # Graphique
  dat <- seq(as.Date(dates[[1]][1]),as.Date(dates[[1]][2]),"days")
  pos <- c(1,match(as.character(seq(1850,2010,10)),substr(as.character(dat),1,4))[-1])
  
  
  png(filename = paste0("2_Travail/1_Past/Rresults/plot.evolution.descr/plot_evolution_",descr,"_liss=",liss,".png"),width = 8,height = 5,units = "in",res = 200)
  plot(dat,des[[1]],type="n",ylim=range(des),
       xlab="Year",ylab=nam2str(descr))
  
  
  
  axis(1,pos,as.character(seq(1850,2010,10)))
  lines(seq(as.Date(dates[[2]][1]),as.Date(dates[[2]][2]),"days"),des[[2]],col="red",lwd=2)
  lines(seq(as.Date(dates[[3]][1]),as.Date(dates[[3]][2]),"days"),des[[3]],col="blue")
  grid()
  graphics.off()
}

# Travail passé:

# - courbe de 1850 à 2010 de l'évolution des trois indicateurs pour les trois rénanalyses:
#   dans cette configuration, les indicateurs sont calculés en allant chercher les analogues sur
#   toute la periode disponible.

# - nuages de points des indicateurs VS dP coloriés par densité pour les trois reanalyses et 
#   pour des fenêtres de 50 ans (1851-1900; 1901-1950; 1951-2010) et de 25 ans
#   (1851-1875; 1876-1900; 1901-1925; 1926-1950; 1951-1975; 1976-2010).

# - construction cdf precip a partir des analogues au sens des indicateurs, 
#   un a un pour voir leur impact respectif, puis couple singnei/rsingnei: 
#   comment évolue la proba de precip annuelle (RP1an)? On va chercher les plus proches en celerite,
#   singularite et singularite relative dans la periode d'obs de 1950 a 2010 en gardant le fait d'avoir
#   calcule les analogues sur toute la periode de la reanalyse.
