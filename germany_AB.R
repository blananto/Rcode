source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')

# Extrait les caracteristiques des LSC de la crue
study.lsc.flood <- function(date,nbdays,descr,start="1950-01-01",end="2021-07-17",rean="ERA5"){
  
  dates <- getdates(start,end)
  
  # Imports
  ind <- get.descriptor(descriptor = descr,k = 1,dist = "TWS",nbdays = nbdays,start = start,end = end,
                        standardize = F,rean = rean,threeday = F,desais = F,period = "past",start.ana = start,end.ana = end)
  ind <- ecdf(ind)(ind)
  
  load(paste0("2_Travail/1_Past/",rean,"/save.ana/ana_TWS_",rean,"_k1_mean",nbdays,"day_",start,"_",end,"_ana_",start,"_",end,".Rdata"))
  
  # Traitement et sortie
  pos <- which(dates==date)
  descr.date <- ind[pos]
  descr.ana <- ind[nei[[pos]][1:(length(dates)*0.002)]+1]
  print(dates[nei[[pos]][1:50]])
  print(ind[nei[[pos]][1:50]])
  tmp <- nei[[pos]][1:50]
  
  boxplot(descr.ana,ylim=c(0,1))
  points(descr.date,col="red",pch=19)
  return(tmp)
}

# Trace la carte des analogues de la sequence de crue
plot.ana.flood <- function(date="2021-07-12",nbdays=3,descr,start="1950-01-01",end="2021-07-17",rean="ERA5"){
  
  dates <- getdates(start,end)
  load(paste0("2_Travail/1_Past/",rean,"/save.ana/ana_TWS_",rean,"_k1_mean",nbdays,"day_",start,"_",end,"_ana_",start,"_",end,".Rdata"))
  
  # Cartes
  pos <- which(dates==date)
  dates.ana <- dates[nei[[pos]][1:50]]
  
  for(i in 1:length(dates.ana)){
    print(paste0(i,"/",length(dates.ana),";",dates.ana[i]))
    map.geo(date = dates.ana[i],rean = "ERA5",climat = NULL,run = 1,k = 1,nbdays = 3,save = T,
            win = T,let = F,leg = T,iso = T,wind = F,condens = F,ssp = F)
  }
}
