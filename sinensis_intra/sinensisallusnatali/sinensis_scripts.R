####COORDINATES FILTERING

library(spThin)
sinensis<-read.csv(file="sinensis_raw.csv", header=TRUE, sep=",")
sinensist<-thin(
sinensis,
lat.col = "decimalLatitude",
long.col = "decimalLongitude",
spec.col = "scientific_name",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="all_new", #####CHOOSE OUTPUT FOLDER, HERE I PUT ONE I HAD IN MY HARDRIVE
out.base = "thinned_sinensis_raw",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)

###LOADING USEFUL PACKAGES (BATES ET AL. 2020)

library(ade4)
library(sf)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(maptools)
library(dismo)
library(ecospat)
library(geosphere)
library(countrycode)
library(hypervolume)
library(outliers)
library(class)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)
library(ggpubr)
library(grDevices)
library(extrafont)
library(factoextra)

###LOADING POINTS

sinensis<-read.csv(file="/run/media/brcuser/MDVBook/possible_pate/new_riel/all_new/thinned_sinensis_raw_thin1.csv",header=TRUE)
coordinates(sinensis)<-c("longitude","latitude")
sinensis_all<-na.omit(sinensis)

####LOADING RASTERS

lista_variables2 <-list.files(pattern='*.tif')
datafiles <- Sys.glob("*.tif")
a <- stack(lista_variables2)

####CROP NATURAL POINTS

possnat<-extent(101.589483, 144.396512, 17.665114, 20.128202)
sinensissp<-SpatialPoints(sinensis_all)
sinensispossnat<-crop(sinensissp,possnat)

###CROP ALIEN POINTS (US)

usa<-extent(-155.678081,  -60.249979, 25.073150, 48.8996)
sinensisusa<-crop(sinensissp,usa)

###ATTACHING VARIABLES TO POINTS

sinensis_usab<-raster::extract(a,sinensisusa)
sinensis_usac<-cbind(sinensisusa@coords,sinensis_usab)
sinensis_usac<-na.omit(sinensis_usac)
sinensis_natb<-raster::extract(a,sinensispossnat)
sinensis_natc<-cbind(sinensispossnat@coords,sinensis_natb)
sinensis_natc<-na.omit(sinensis_natc)
write.csv(sinensis_natc, "nat_values_sinensis.csv")
write.csv(sinensis_usac, "ex_values_sinensis.csv")

###I THEN COMBINED MANUALLY THE TWO DATASET ON LIBREOFFICE CALC, AFTER ADDING THE "STATUS" WITH THE "N" AND "E" COLUMN

###LOAD THE DATASET

sinensis_values<-read.csv(file="all_values_sinensis.csv",header=TRUE)

###PUT THE FUNCTIONS FROM BATES ET AL. 2020

"overlap.eq.gen" <- function(repi, z1, z2) {
sim.o.D <- NaN
while (is.nan(sim.o.D)) {
# overlap on one axis
##z1$sp is the environmental values for occurences of the species
occ.pool <- c(z1$sp, z2$sp)  # pool of random occurrences
occ.pool <- occ.pool[sample(length(occ.pool))]
row.names(occ.pool)<-c() # remove rownames
rand.row <- base::sample(1:length(occ.pool), length(z1$sp)) # random reallocation of occurrences to dataset
sp1.sim <- occ.pool[rand.row]
sp2.sim <- occ.pool[-rand.row]
z1.sim <- ecospat.grid.clim.dyn(z1$glob, z1$glob1, data.frame(sp1.sim), R = length(z1$x), th.sp=0, th.env=0)  # gridding in this function there is something
z2.sim <- ecospat.grid.clim.dyn(z2$glob, z2$glob1, data.frame(sp2.sim), R = length(z2$x), th.sp=0, th.env=0)
o.i <- ecospat.niche.overlap(z1.sim, z2.sim, cor = F)  # overlap between random and observed niches FAlSE because all environments
sim.o.D <- o.i$D  # storage of overlaps
sim.o.I <- o.i$I
}
return(c(sim.o.D, sim.o.I))
}
"niche.equivilency" <- function (z1, z2, rep, alternative = "greater", ncores = 1)
{
R <- length(z1$x) #100
l <- list()
obs.o <- ecospat.niche.overlap(z1, z2, cor = FALSE) 	#observed niche overlap, FAlSE because all environments
if (ncores == 1) {
sim.o <- as.data.frame(matrix(unlist(lapply(1:rep, overlap.eq.gen, z1=z1, z2=z2)), byrow = TRUE, ncol = 2))
}
colnames(sim.o) <- c("D", "I")
l$sim <- sim.o
l$obs <- obs.o
if (alternative == "greater") {
l$p.D <- (sum(sim.o$D >= obs.o$D) + 1)/(length(sim.o$D) +  1)
l$p.I <- (sum(sim.o$I >= obs.o$I) + 1)/(length(sim.o$I) + 1)
}
if (alternative == "lower") {
l$p.D <- (sum(sim.o$D <= obs.o$D))/(length(sim.o$D)) ###why is it +1
l$p.I <- (sum(sim.o$I <= obs.o$I) + 1)/(length(sim.o$I) + 1)
}
return(l)
}
"map_overlap" <- function(y, status){
y_E <- y[status == "E"]
y_N <- y[status == "N"]
dens_E <- density(y_E)
dens_N <- density(y_N)
xlim <- range(dens_E$x, dens_N$x) #x axis range
ylim <- range(0, dens_E$y, dens_N$y)
ylim[2] <- ylim[2] + 0.1
hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = "T. sinensis", xlab = "Between-Class PCA 1", ylab = "Density", cex.lab=1.5, cex.axis=1.5)
hist(y_N, proba = T, add = T, col = 'white', density = 10, angle = 45)
polygon(dens_E, density = -1, col = col_E)
polygon(dens_N, density = -1, col = col_N)
mat <- cbind(dens_E$x, dens_E$y)
mat <- rbind(mat, mat[1,])
pol_E <- st_polygon(list(mat))
mat <- cbind(dens_N$x, dens_N$y)
mat <- rbind(mat, mat[1,])
pol_N <- st_polygon(list(mat))
pol_inter <- st_intersection(st_buffer(pol_E,0),  st_buffer(pol_N,0))
plot(pol_inter, add = T, col = "grey", ylim=c(0,1))
rug(y_N, side = 1, line = -0.15, col = col_N, tck = 0.03, lwd=2)
rug(y_E, side = 1, line = -0.15, col = col_E, tck = 0.03, lwd=2)
}

col_E="yellow"
col_N="green"

###MAKE THE ANALYSES (BATES ET AL. 2020)

pca.env <- dudi.pca(sinensis_values[,c(8,9,13,14,15,16,18,21,22)], scannf = F, nf = 10) ##so can first make a pca of all, look at eigen values ect.
status <- as.factor(sinensis_values$status)
bet1 <- bca(pca.env, status, scan = FALSE, nf = 9)
#bet1$ratio #between-class ratio
## explore correlaitons to each bio variable to axis
correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis.
correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
correlation_to_bio_variables <- t(correlation_to_bio_variables)
## save interesting values
between_class <- as.data.frame(bet1$ls) #data frame of the first axis of between class analysis
data <-  cbind(status, between_class) #binding the two together
#the scores for the whole environment
scores.globclim <- data$CS1
#scores for the presences
scores.sp.ex <- data$CS1[which(data$status=='E')]
scores.sp.nat <-data$CS1[which(data$status=='N')]
#scores for the whole study region, same as presences as no considering whole study area
scores.clim.ex <- data$CS1[which(data$status=='E')]
scores.clim.nat <-data$CS1[which(data$status=='N')]
bet1 <- bca(pca.env, status, scan = FALSE, nf = 9)
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
#glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only,
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
## Overlap
# D overlap
D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap.
# Niche expansion : delimiting niche categories and quantifying niche dynamics in analogue climates
niche.dyn.whole <- ecospat.niche.dyn.index (grid.clim.nat, grid.clim.inv, intersection = NA) #if intersection= NA means no intersection of analogous envirnoment
niche.dyn.whole$dynamic.index.w
## Equivilency Test
eq.test <- NA
eq.test <- niche.equivilency(grid.clim.nat, grid.clim.inv, rep=1000, alternative = "lower", ncores = 1)
row <- cbind.data.frame(D.overlap, eq.test$p.D, niche.dyn.whole$dynamic.index.w[[1]],niche.dyn.whole$dynamic.index.w[[2]], niche.dyn.whole$dynamic.index.w[[3]])
names(row) <- c('D.overlap', 'eq.test', 'expansion', 'stability', 'Unfilling')
niche_expansion_1<-NULL
niche_expansion_1 <- rbind(niche_expansion_1, row)
View(niche_expansion_1)
write.csv(niche_expansion_1, "ne_sinensis.csv")
pdf(file="sinensis_asia_vs_usa_greenyellow.pdf")
map_overlap(bet1$ls[,1], status)
dev.off()

###NICHE DYNAMICS PLOT

pdf(file="niche_dyn_sinensis.pdf")
ecospat.plot.niche.dyn(grid.clim.nat, grid.clim.inv, quant = 0, title = "", name.axis1 = "Axis 1",
name.axis2 = "Axis 2", interest = 1, col.unf =
"green", col.exp = "red", col.stab = "blue", colZ1 =
"green3", colZ2 = "red3", transparency = 70)
dev.off()

###PCA PLOTS

pdf(file = "sinensis_pca_2ax.pdf")
fviz_pca_var(pca.env,
col.var = "contrib", # Color by contributions to the PC
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE     # Avoid text overlapping
)
dev.off()
pdf(file = "sinensis_corr.pdf")
fviz_contrib(pca.env, choice="var", axes = 1:2, title="T. sinensis (all)")
dev.off()

###MANTEL
pdf(file="mantel_sinensis.pdf")
ecospat.mantel.correlogram(dfvar=sinensis_values,colxy=2:3, n=100, colvar=c(8,9,13,14,15,16,18,21,22), max=50000, nclass=10, nperm=1000)
dev.off()

####MAP

library(tidyverse)
world_coordinates <- map_data("world")
world<-world_coordinates
sinensis_n<-sinensis_values[sinensis_values$status == "N", ]
sinensis_e<-sinensis_values[sinensis_values$status == "E", ]

pdf(file="sinensis_map.pdf")
ggplot() +
geom_map(
data = world, map = world,
aes(long, lat, map_id = region),
color = "grey", fill= "white"
)+
geom_point(
data = sinensis_n,
aes(decimalLongitude, decimalLatitude, color = "green"),
alpha = 1
) +
geom_point(
data = sinensis_e,
aes(decimalLongitude, decimalLatitude, color = "yellow"),
alpha = 1
)
