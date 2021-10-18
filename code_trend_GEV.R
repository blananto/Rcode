library(ismev)
library(evd)

# Lancer le code ligne par ligne pour l'adapter
# Verifier les indices d'hiver, pour voir si plus coherent de ne lancer le calcul que sur les hivers complets, comme pour lettre et article 2

meanGEV.fct<-function(mu,sig,xi){mu+sig/xi*(gamma(1-xi)-1)}

nbday<-1 #mettre nbday<-3 pour les precip sur 3j 

tab<-NULL
tab.rown<-NULL

yrange<-c(1950,2010)

par(mfrow=c(2,2))


for (region in c("Isere-seul","Drac-seul")){
  for (db in c("tps","spazm")){
    
    load(paste0("/media/blanchej/TOSHIBA EXT/bckp2019/Work/InProgress/Antoine Blanc/Rcode/precip_",nbday,"j/precip_",region,"_",nbday,"day_",yrange[1],"-01-01_",yrange[2],"-12-31_",db,".Rdata"))
    print(unique(SPAZM_BV_big_gre$BV))
      
    date<-seq(as.Date(paste0("01/01/",yrange[1]),"%d/%m/%Y"),as.Date(paste0("31/12/",yrange[2]),"%d/%m/%Y"),by=1) 
    if (nbday==3) date<-date[-c(1,length(date))]
    
    years<-as.numeric(substr(date,1,4))
    months<-as.numeric(substr(date,6,7))
    years[months==12]<-years[months==12]+1
    
    for (season in c("MAM","JJA","SON","DJF")){
      
      print(paste0("season ",season))
      
      if (season=="MAM") mm<-3:5
      if (season=="JJA") mm<-6:8
      if (season=="SON") mm<-9:11
      if (season=="DJF") mm<-c(12,1,2)
      
      
      max.vec<-rep(NA,length(yrange[1]:yrange[2]))
      names(max.vec)<-yrange[1]:yrange[2]
      for (yy in yrange[1]:yrange[2]){
        if (season!="DJF" | (season=="DJF" & yy>yrange[1])){
          wh<-which(years==yy & months %in% mm)
          max.vec[as.character(yy)]<-max(precip[wh],na.rm=TRUE)
        }
      }
      print(yrange)
      print(max.vec)
      plot(yrange[1]:yrange[2],max.vec,xlab="",ylab="max (mm)");title(season)
      
      covar<-(yrange[1]:yrange[2]-yrange[1])/(yrange[2]-yrange[1])

      idx<-which(is.na(max.vec))
      if (length(idx)>0) {max.vec<-max.vec[-idx];covar<-covar[-idx]}
      ny<-length(covar)
      fit0<-gev.fit(max.vec,show=F)
      fit1<-gev.fit(max.vec,as.matrix(covar),mul=1,show=F)
      fit2<-gev.fit(max.vec,as.matrix(covar),sigl=1,show=F)
      fit3<-gev.fit(max.vec,as.matrix(covar),mul=1,sigl=1,show=F)
      
      pval1<-1-pchisq(2*(fit0$nllh-fit1$nllh),df=1)
      pval2<-1-pchisq(2*(fit0$nllh-fit2$nllh),df=1)
      pval3<-1-pchisq(2*(fit0$nllh-fit3$nllh),df=2)
      
      idx.min<-which.min(c(pval1,pval2,pval3))
      print(idx.min)
      pval.min<-min(c(pval1,pval2,pval3))
      print(pval.min)
      print(c(pval1,pval2,pval3))
      
      if (idx.min==1) fitbest<-fit1
      if (idx.min==2) fitbest<-fit2
      if (idx.min==3) fitbest<-fit3
      print(fitbest$mle)
      
      RL20.first<-qgev(1-1/20,fitbest$vals[1,1],fitbest$vals[1,2],fitbest$vals[1,3])
      RL20.last<-qgev(1-1/20,fitbest$vals[ny,1],fitbest$vals[ny,2],fitbest$vals[ny,3])
      mean.first<-meanGEV.fct(fitbest$vals[1,1],fitbest$vals[1,2],fitbest$vals[1,3])
      mean.last<-meanGEV.fct(fitbest$vals[ny,1],fitbest$vals[ny,2],fitbest$vals[ny,3])
      
      tab<-rbind(tab,c(RL20.first,RL20.last,(RL20.last-RL20.first)/RL20.first,mean.first,mean.last,(mean.last-mean.first)/mean.first,pval.min))
    }
    
    tab.rown<-c(tab.rown,paste0(region," ",db," ",c("MAM","JJA","SON","DJF")))
  }
}
rownames(tab)<-tab.rown
colnames(tab)<-c(paste0("RL20 in ",yrange[1]),paste0("RL20 in ",yrange[2]),"% of change RL20",paste0("mean in ",yrange[1]),paste0("mean in ",yrange[2]),"% of change mean","pval")
tab

write.csv(tab,file=paste0("/media/blanchej/TOSHIBA EXT/bckp2019/Work/InProgress/Antoine Blanc/Rcode/precip_",nbday,"j/trends.csv"))
