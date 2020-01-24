setwd("C:/Users/blananto/Documents")

# Parametres ----
rean     <- "20CR"
member   <- 1
seasonal <- FALSE
k        <- 1
dist     <- "TWS"
nbdays   <- 3
N        <- "02" # %age de voisins selectionnes pour l'analogie classique
Q        <- "05" # %age de voisins selectionnes pour construire les indicateurs 
M        <- "nrn05" # %age de voisins selectionnes pour l'analogie indicateurs
start    <- "1950-01-01"
end      <- "2011-12-31"

# Import des sorties de la reanalyse et calcul des distances ----
load.nc(rean = rean)
compute_dist.gen(k = k, dist = dist, start = start, end = end, rean = rean)

# Calcul du score climatologique ----
compute_score_climato(nbdays = nbdays, start = start, end = end, rean = rean)

# Analogie classique ----
save.nei.A(k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)
fit.empir.A(rad = N, k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)
fit.loglik.p0.A(rad = N, k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)
compute_crps.A(rad = N, k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)
compare.crps.A(k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)

# Analogie indicateurs ----
compute_criteria(k = k, dist = dist, start = start, end = end, rean = rean, update = TRUE)
run(k = k, dist = dist, nbdays = nbdays, str = Q, radtype = M, start = start, end = end, rean = rean) # possibilite d'ajuster plusieurs lois, de modifier les couples d'indicateurs
compare.crps(which = "", k = k, dist = dist, nbdays = nbdays,radtype = M, start = start, end = end, rean = rean)
compare.crps.sais(k = k, dist = dist, nbdays = nbdays, start = start, end = end, radtype = M, rean = rean)

# Visualisation dans le plan des indicateurs ----
descr <- c("sing05","rsing05")
plot.empir(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end, radtype = M, rean = rean, empir=TRUE, obs=FALSE)
plot.empir.mean(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end,
                radtype = M, rean = rean, ref="1900-01-01")

# Version generique
fit.empir(rean = c("20CR","20CR"), k = c(1,1), descriptors = c("rsing05","rsing05"), dist = c("TWS","RMSE"),
          nbdays = 3, start = start, end = end, radtype = M)
plot.empir.clean.obs(rean = c("20CR","20CR"), k = c(1,1), descriptors = c("sing05","sing05"), dist = c("TWS","RMSE"),
           nbdays = 3, start = start, end = end, radtype = M,dP = T,coin = T)
compute_crps(descriptors = c("rsing05","rsing05"), k = c(1,1), dist = c("TWS","RMSE"), nbdays = 3,
             start = start, end = end, radtype = M, rean = c("20CR","20CR"))
compare.crps(which = "", k = k, dist = dist, nbdays = nbdays, start = start, end = end, radtype = M, rean = rean)

# Analogie en deux etapes: selection analogie classique puis sous selection indicateurs ----
descr <- c("sing05","rsing05")

fit.empir.TL(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end,
             radAna = "10", radInd = "05", rean = rean)
compute_crps_TL(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end,
                radAna = "10", radInd = "05", rean = rean)
compare.crps.TL(k = k, dist = dist, nbdays = nbdays, start = start, end = end,
                radAna = "1", radInd = "05", rean = rean)

# Ensuite, les trois fonctions d'un coup avec differents rayons
run.TL(k = k, dist = dist, nbdays = nbdays, str = Q, radAna = "05", radInd = "02", start = start, end = end, rean = rean)


# Analogie double combinee ----
save.score.A(k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)
fit.empir.comb(descriptors = c("cel","sing05"), k = k, dist = dist, nbdays = nbdays, start = start,
               end = end, p = 1, rean = rean)

run.comb(k = k, dist = dist, nbdays = nbdays, str = Q, start = start, end = end, p = 0, rad = "05", rean = rean)
compare.crps.comb(k = k, dist = dist, nbdays = nbdays, start = start, end = end, p = 0, rad = "05", rean = rean)

# Tests distribution
summary.distrib(descriptors = c("persnei","singnei"),k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean, rev = TRUE)
plot.distrib(rean = rean,k = k,descriptors = c("persnei","singnei"),dist = dist,nbdays = nbdays,start = start,end = end)
plot.crps(rean = rean,k = k,descriptors = c("singnei","rsingnei"),dist = dist,nbdays = nbdays,start = start,end = end)
image.region()
image.cumul()
