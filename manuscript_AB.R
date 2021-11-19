library(raster) # database France
library(mapproj) # database Europe

# setwd("../../")
# Site utile: https://dimension.usherbrooke.ca/pages/32

# Trace la carte de l'Europe avec frontières Auer et Grenoble
image.europe.auer <- function(){
  
  # Import
  ligne=read.csv("2_Travail/Data/Carto/Shape_CRSM.csv")
  
  # Mise en forme
  ligne.indiv <- unique(ligne[,1])
  pos.rect <- match(ligne.indiv,ligne[,1])
  rectang <- ligne[pos.rect,2:3]
  rectang <- apply(rectang,2,range)
  ligne <- ligne[-pos.rect,]
  
  # Carte
  png("2_Travail/0_Present/Rresults/image.europe.auer/map_europe_auer.png",width = 6,height = 4,units = "in",res = 600)
  par(mar=c(0,0,0,0))
  map(database= "world", ylim=c(39,52), xlim=c(-5,21), col=alpha("darkgreen",alpha = 0.3),fill=T)
  rect(xleft = rectang[1,1],ybottom = rectang[1,2],xright = rectang[2,1],ytop = rectang[2,2],lwd=2,border="black",lty=2)
  for(i in 1:length(ligne.indiv)){
    lines(ligne[ligne[,1]==ligne.indiv[i],2],ligne[ligne[,1]==ligne.indiv[i],3],lwd=3,col="purple")
  }
  points(x = 5.73,y = 45.18,pch=19,cex=1.5,col="red")
  box()
  graphics.off()
}

# Trace la carte de la France avec Grenoble
image.france <- function(){
  
  # Import
  FranceFormes <- getData(name="GADM", country="FRA", level=0)
  
  # Carte
  png("2_Travail/0_Present/Rresults/image.france/map_france_grenoble.png",width = 5,height = 5,units = "in",res = 600)
  par(mar=c(0,0,0,0))
  plot(FranceFormes)
  points(x = 5.73,y = 45.18,pch=19,cex=1.5,col="red"); text(x = 5.73-0.4,y = 45.18-0.4,"Grenoble",col="red",font=4) # Grenoble
  points(x = 2.35,y = 48.86,pch=19); text(x = 2.35+0.4,y = 48.86+0.4,"Paris",font=3) # Paris
  points(x = 4.83,y = 45.76,pch=19); text(x = 4.83+0.4,y = 45.76+0.4,"Lyon",font=3) # Lyon
  points(x = 5.38,y = 43.30,pch=19); text(x = 5.38+0.4,y = 43.30+0.4,"Marseille",font=3) # Marseille
  graphics.off()
}




