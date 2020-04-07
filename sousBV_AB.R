source('2_Travail/Rcode/utils_AB.R', encoding = 'UTF-8')
source('2_Travail/Rcode/article1_AB.R', encoding = 'UTF-8')

# Compparaison des max annuels de precip
compare.extr <- function(bv=c("Isere-seul","Drac-seul"),nbdays=3,start="1950-01-01",end="2011-12-31"){
  
  # Import des precip
  precip1 <- get.precip(nbdays,start,end,bv[1])
  precip2 <- get.precip(nbdays,start,end,bv[2])
  
  # Max annuels
  ind1 <- get.ind.max(precip1,type = "year",nbdays,start,end)
  ind2 <- get.ind.max(precip2,type = "year",nbdays,start,end)
  
  # 62 max
  ind1 <- get.ind.extr(62,nei=T,bv=bv[1])
  ind2 <- get.ind.extr(62,nei=T,bv=bv[2])
  
  # Distribution
  plot(ecdf(precip1[ind1]))
  plot(ecdf(precip2[ind2]),col="red",add=T)
  
  
}
  
# Generation du trace propre des deux BVs, ajout du petit BV manquant
create.drac.isere <- function(){
  
  # BV complet Isere@Grenoble
  bord<-read.csv("2_Travail/Data/Carto/borders-lambert93-lambertII-isere.csv",sep=";")
  bord<-cbind(bord$XLII.m,bord$YLII.m)
  bv_all <- as(bord,"gpc.poly")
  plot(bv_all,poly.args = list(col = 1))
  
  # BV incomplet Drac@Grenoble
  bord1<-read.csv("2_Travail/Data/Carto/border_Isere@Grenoble_old.csv",sep=",")
  bord1 <- bord1[bord1[,1] %in% c(670,19,35,124,248,176,596,104),] # Drac
  ui<-unique(bord1[,1])
  
  a1<-as(bord1[bord1[,1]==unique(bord1[,1])[1],2:3],"gpc.poly")
  a2<-as(bord1[bord1[,1]==unique(bord1[,1])[2],2:3],"gpc.poly")
  bv_drac_old<-union(a1,a2)
  
  for (i in ui[-(1:2)]) {
    a1<-as(bord1[bord1[,1]==i,2:3],"gpc.poly")
    bv_drac_old<-union(a1,bv_drac_old)
  }
  
  plot(bv_drac_old,poly.args = list(col = 2),add=T)
  
  # Difference des deux pour BVs seuls
  delta <- setdiff(bv_all,bv_drac_old) # setdiff, intersect, union tres utiles pour manipuler les gpc.poly!
  bv_isere_seul <- as(cbind(delta@pts[[1]]$x,delta@pts[[1]]$y),"gpc.poly")
  plot(bv_isere_seul,poly.args = list(col = 3),add=T)
  
  bv_drac_manquant <- as(cbind(delta@pts[[2]]$x,delta@pts[[2]]$y),"gpc.poly")
  plot(bv_drac_manquant,poly.args = list(col = 4),add=T)
  bv_drac_seul <- union(bv_drac_old,bv_drac_manquant)
  plot(bv_drac_seul,poly.args = list(col = 5),add=T)
  
  # Export des csv
  isere <- cbind(bv_isere_seul@pts[[1]]$x,bv_isere_seul@pts[[1]]$y)
  colnames(isere) <- c("XLII.m","YLII.m")
  write.table(x = isere,file = "2_Travail/Data/Carto/borders-lambertII-isere-seul.csv",sep = ";",
            row.names = T, col.names = T)
  
  drac <- cbind(bv_drac_seul@pts[[1]]$x,bv_drac_seul@pts[[1]]$y)
  colnames(drac) <- c("XLII.m","YLII.m")
  write.table(x = drac,file = "2_Travail/Data/Carto/borders-lambertII-drac-seul.csv",sep = ";",
             row.names = T, col.names = T)
  
  # Ajout du sous BV manquant au csv old
  manquant <- cbind(bv_drac_manquant@pts[[1]]$x,bv_drac_manquant@pts[[1]]$y)
  manquant <- cbind (1,manquant) # on lui donne le numero 1
  colnames(manquant) <- c("Polygone","XLII.m","YLII.m")
  
  bord1 <- read.csv("2_Travail/Data/Carto/border_Isere@Grenoble_small_bv.csv",sep=";")
  colnames(bord1) <- c("Polygone","XLII.m","YLII.m")
  bord1 <- rbind(bord1,manquant)
  write.table(x = bord1,file = "2_Travail/Data/Carto/border_Isere@Grenoble_small_bv.csv",sep = ";",
              row.names = T, col.names = T)
}

# Groupement des small BV en sous BV de tailles equivalentes
create.sousBV <- function(){
  
  png("2_Travail/Rresults/create.sousBV/map_smallBV_sousBV.png",width = 700,height = 700,units = "px")
  # Fond de carte
  load(file=paste("2_Travail/Data/Carto/griddata_1x1_IsereSavoieHautesAlpes.Rdata",sep=""))
  Fx<-griddata$Fx
  Fy<-griddata$Fy
  Fz<-griddata$Fz*1000
  image.plot(Fx,Fy,Fz,col=gray(seq(0.1,0.99,length=100)),xlab="X (km) - Lambert II extended",ylab="Y (km) - Lambert II extended",legend.line=-2.3, cex.axis=1.3, cex.lab=1.3)
  
  # Trace des small BV
  bord <- read.csv("2_Travail/Data/Carto/border_Isere@Grenoble_small_bv.csv",sep=";")
  small_bv <- unique(bord[,1])
  for(i in 1:length(small_bv)){
    data <- bord[bord[,1]==small_bv[i],]
    center.x <- mean(data[,2])
    center.y <- mean(data[,3])
    lines(data[,2:3],col="blue",lwd=3,lty=2)
    text(center.x,center.y,small_bv[i],col="white",font=2)
  }
  
  # Regroupement en sous BV
  sous_bv <- list()
  sous_bv[[1]] <- c(797,577,621,706,87,567)
  sous_bv[[2]] <- 136
  sous_bv[[3]] <- 104
  sous_bv[[4]] <- c(596,176,248,124,35,1)
  sous_bv[[5]] <- c(180,1092,116,792,670,19,1065,934)
  names(sous_bv) <- c("tarentaise","maurienne","romanche","drac","gresivaudan")
  
  res <- list()
  for(i in 1:length(sous_bv)){
    bv <- sous_bv[[i]]
    a<-as(bord[bord[,1]==bv[1],2:3],"gpc.poly")
    if(length(bv)>1){
      for(j in 2:length(bv)){
        a2 <- as(bord[bord[,1]==bv[j],2:3],"gpc.poly")
        a <- union(a,a2)
      }
    }
    res[[i]] <- cbind(a@pts[[1]]$x,a@pts[[1]]$y)
    print(paste0(names(sous_bv)[i],": ",round(geometry::polyarea(a@pts[[1]]$x,a@pts[[1]]$y),0)," km2"))
    lines(res[[i]],col="red",lwd=3)

    # Export csv
    tab <- res[[i]]
    colnames(tab) <- c("XLII.m","YLII.m")
    rownames(tab) <- 1:nrow(tab)
    write.table(x = tab,file = paste0("2_Travail/Data/Carto/borders-lambertII-",names(sous_bv)[i],".csv"),sep = ";",
                row.names = T, col.names = T)
  }
  graphics.off()
}

# Recherche du bv manquant (branche gauche du Y grenoblois) et l'ajoute aux small bv
find.missing.bv <- function(){
  
  # Carte de base
  load(file=paste("2_Travail/Data/Carto/griddata_1x1_IsereSavoieHautesAlpes.Rdata",sep=""))
  Fx<-griddata$Fx
  Fy<-griddata$Fy
  Fz<-griddata$Fz*1000
  image.plot(Fx,Fy,Fz,col=gray(seq(0.1,0.99,length=100)),xlab="X (km) - Lambert II extended",ylab="Y (km) - Lambert II extended",legend.line=-2.3, cex.axis=1.3, cex.lab=1.3)
  
  # Sous bv de Juliette
  pol<-readShapePoly("2_Travail/Data/Carto/SousSecteurHydro/SousSecteurHydro_FXX.shp")
  
  for (i in 1:length(pol@polygons)) {
    tmp<-pol@polygons[[i]]@Polygons[[1]]@coords;
    dfr <- data.frame(x = tmp[,1], y = tmp[,2]);
    coordinates(dfr) <- ~ x + y;
    proj4string(dfr) <- CRS("+proj=lcc +lat_1=49 +lat_2=44 +lat_0=46.5 +lon_0=3 +x_0=700000 +y_0=6600000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs");
    newdfr<-spTransform(dfr,CRS("+proj=lcc +lat_1=46.8 +lat_0=46.8 +lon_0=0 +k_0=0.99987742 +x_0=600000 +y_0=2200000 +a=6378249.2 +b=6356515 +towgs84=-168,-60,320,0,0,0,0 +pm=paris +units=m +no_defs"));
    center.x <- mean(newdfr@coords[,1]/1000)
    center.y <- mean(newdfr@coords[,2]/1000)
    lines(newdfr@coords[,1]/1000,newdfr@coords[,2]/1000,col="white",lwd=2)
    text(center.x,center.y,i,col="white")
    if(i==934) bon <- cbind(newdfr@coords[,1]/1000,newdfr@coords[,2]/1000)
  } # notre bv manquant: numero 1065, 934 aussi?
  
  bord <- read.csv("2_Travail/Data/Carto/border_Isere@Grenoble_small_bv.csv",sep=";")
  bon <- cbind(rep(934,nrow(bon)),bon)
  colnames(bon) <- colnames(bord)
  bord <- rbind(bord,bon)
  write.table(x = bord,file = "2_Travail/Data/Carto/border_Isere@Grenoble_small_bv.csv",sep = ";",
              row.names = T, col.names = T)
}

# Interpolation TPS 1x1km2 puis moyenne
make.precip<-function(start="1950-01-01",end="2011-12-31") {
  
  # Import coordonnees bv et donnees precip
  bord<-read.csv("2_Travail/Data/Carto/borders-lambertII-isere-seul.csv",sep=";")
  load("2_Travail/Data/Precip/data_precip.Rdata") #precip journalieres
  
  # Indice des 1er et derniers jours dans data precip 
  id1<-which(rownames(data$precip)==paste0(substr(start,1,4),substr(start,6,7),substr(start,9,10)))
  id2<-which(rownames(data$precip)==paste0(substr(end,1,4),substr(end,6,7),substr(end,9,10)))
  
  # Quadrillage 1x1km2 rectangulaire
  grids<-merge((seq(min(bord[,1]),max(bord[,1]),by=1)),(seq(min(bord[,2]),max(bord[,2]),by=1))) #plot(grids)
  
  # Points de grille faisant partie du bv
  idx <- inpip(grids,bord) #points(grids[idx,],col="red")
  
  # Interpolation
  ivec<-id1:id2
  precip<-rep(NA,length(ivec))
  for (i in ivec) {
    if (i %%50==0) print(i)
    fit <- Tps(data$info_stations[,c("XLII.km.","YLII.km.")],data$precip[i,]/10,give.warnings=F) # sur tous les pluvios
    tmp <- predict(fit,grids[idx,])  # estimations des points etant dans le bv
    tmp[tmp<0]=0  # si interpolation negetive => 0
    precip[i]<-mean(tmp,na.rm = T) # precip: moyenne des pts de grille sur le bv
  }
  precip
}

