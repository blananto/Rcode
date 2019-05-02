setwd("C:/Users/blanca/Documents")

# Parametres ----
rean     <- "ERA20C"
member   <- 1
seasonal <- FALSE
k        <- 1
dist     <- "TWS"
nbdays   <- 3
N        <- "02" # %age de voisins selectionnes pour l'analogie classique
Q        <- "05" # %age de voisins selectionnes pour construire les indicateurs 
M        <- "nrn05" # %age de voisins selectionnes pour l'analogie indicateurs
start    <- "1950-01-01"
end      <- "2010-12-31"

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
compute_criteria(k = k, dist = dist, start = start, end = end, rean = rean)
run(k = k, dist = dist, nbdays = nbdays, str = Q, radtype = M, start = start, end = end, rean = rean) # possibilite d'ajuster plusieurs lois, de modifier les couples d'indicateurs
compare.crps(which = "k_dist", dist = "TWS_RMSE", nbdays = nbdays,radtype = M, start = start, end = end, rean = rean)


# Visualisation dans le plan des indicateurs ----
descr <- c("sing05","rsing05")
plot.empir.mean(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end,
                radtype = M, rean = rean, ref="1900-01-01")

# Analogie en deux etapes: selection analogie classique puis sous selection indicateurs ----
descr <- c("sing05","rsing05")

fit.empir.TL(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end,
             radAna = "10", radInd = "05", rean = rean)
compute_crps_TL(descriptors = descr, k = k, dist = dist, nbdays = nbdays, start = start, end = end,
                radAna = "10", radInd = "05", rean = rean)
compare.crps.TL(k = k, dist = dist, nbdays = nbdays, start = start, end = end,
                radAna = "10", radInd = "05", rean = rean)

# Ensuite, les trois fonctions d'un coup avec differents rayons
run.TL(k = k, dist = dist, nbdays = nbdays, str = Q, radAna = "05", radInd = "02", start = start, end = end, rean = rean)


# Analogie double combinee ----
save.score.A(k = k, dist = dist, nbdays = nbdays, start = start, end = end, rean = rean)
fit.empir.comb(descriptors = c("cel","sing05"), k = k, dist = dist, nbdays = nbdays, start = start,
               end = end, p = 1, rean = rean)

run.comb(k = k, dist = dist, nbdays = nbdays, str = Q, start = start, end = end, p = 0, rad = "05", rean = rean)
compare.crps.comb(k = k, dist = dist, nbdays = nbdays, start = start, end = end, p = 0, rad = "05", rean = rean)
