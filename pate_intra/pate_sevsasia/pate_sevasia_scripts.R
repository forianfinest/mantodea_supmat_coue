###TRIMMING THE POINTS

####notice that i removed some dubious points from this species before trimming

library(spThin)
pate<-read.csv(file="pate_rielaboration.csv", header=TRUE, sep=",")
View(pate)
patet<-thin(
pate,
lat.col = "decimalLatitude",
long.col = "decimalLongitude",
spec.col = "species",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/possible_pate/new_riel", #####CHOOSE OUTPUT FOLDER
out.base = "thinned_pate",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)

####READ COORDINATES

pate<-read.csv(file="thinned_pate_thin1.csv",header=TRUE,sep = ",")
pate<-na.omit(pate)
coordinates(pate) <- c("longitude", "latitude")
patesp<-SpatialPoints(pate)

###CROP ACCORDING TO AREA

###FIRST ASIA

natpate<-extent(92.144139967402, 149.097264967402, -12.096830673867386, 41.068145322583234)
pate_nat<-crop(patesp,natpate)
write.csv(pate_nat@coords, "pate_nat_tocheck.csv")

###THEN NEPAL, INDIA AND PHILIPPINES

se<-extent(80.34533479121058, 94.05627229121056, 19.699873289716056, 29.128531129347728)
pate_boh<-crop(patesp,se)
View(pate_boh)
write.csv(pate_boh@coords, "points_seandphi.csv")


###GET THE BIOCLIM VARIABLES

lista_variables2 <-list.files(pattern='*.tif')
datafiles <- Sys.glob("*.tif")
a <- stack(lista_variables2)
setwd("/run/media/brcuser/MDVBook/possible_pate/new_riel/all_new")

###LINK VARIABLES TO AREAS

areyoupate<-read.csv(file="points_seandphi.csv",header=TRUE,sep = ",")
coordinates(areyoupate)<-c("longitude","latitude")
sesp<-SpatialPoints(areyoupate)
patesephib<-raster::extract(a,sesp)
sephic<-cbind(areyoupate,patesephib)

patesure<-read.csv(file="pate_nat_nophillipines.csv",header=TRUE,sep = ",")
View(patesure)
coordinates(patesure)<-c("longitude","latitude")
patenatsp<-SpatialPoints(patesure)
patenatb<-raster::extract(a,patenatsp)
patenatc<-cbind(patesure,patenatb)

write.csv(patenatc,"patenat_val.csv")
write.csv(sephic,"sephi_val.csv")

###THEN, I EDITED MANUALLY THE CSV AND MERGED THEM TOGETHER

###RUN THE ANALYSES

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
col_N <- "green"
col_E <- "yellow"
"map_overlap" <- function(y, status){
y_E <- y[status == "SE"]
y_N <- y[status == "N"]
dens_E <- density(y_E)
dens_N <- density(y_N)
xlim <- range(dens_E$x, dens_N$x) #x axis range
ylim <- range(0, dens_E$y, dens_N$y)
ylim[2] <- ylim[2] + 0.1
hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = "patellifera (Asia vs SEA)", xlab = "Between-Class PCA 1", ylab = "Density", cex.lab=1.5, cex.axis=1.5)
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
asia_se<-read.csv(file="check_pate.csv",header=TRUE)
pca.env <- dudi.pca(asia_se[,c(8,9,13,14,15,16,18,21,22)], scannf = F, nf = 10) ##so can first make a pca of all, look at eigen values ect.
status <- as.factor(asia_se$status)
status <- as.factor(status)
bet1 <- bca(pca.env, status, scan = FALSE, nf = 9)
correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis.
correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
correlation_to_bio_variables <- t(correlation_to_bio_variables)
correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis.
correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
correlation_to_bio_variables <- t(correlation_to_bio_variables)
between_class <- as.data.frame(bet1$ls) #data frame of the first axis of between class analysis
data <-  cbind(status, between_class) #binding the two together
#the scores for the whole environment
scores.globclim <- data$CS1
#scores for the presences
scores.sp.ex <- data$CS1[which(data$status=='SE')]
scores.sp.nat <-data$CS1[which(data$status=='N')]
#scores for the whole study region, same as presences as no considering whole study area
scores.clim.ex <- data$CS1[which(data$status=='E')]
scores.clim.nat <-data$CS1[which(data$status=='N')]
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
#glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only,
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap.
## Equivilency Test
eq.test <- NA
eq.test <- niche.equivilency(grid.clim.nat, grid.clim.inv, rep=1000, alternative = "lower", ncores = 1)
View(eq.test)
row <- cbind.data.frame(D.overlap, eq.test$p.D)
names(row) <- c('D.overlap', 'eq.test')
write.csv(row,"pate_nvsse.csv")
pdf(file="posspate_overlap.pdf")
map_overlap(bet1$ls[,1], status)
dev.off()

###PCA PLOTS

pdf(file = "pate_check_pca_2ax.pdf")
fviz_pca_var(pca.env,
col.var = "contrib", # Color by contributions to the PC
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE     # Avoid text overlapping
)
dev.off()
pdf(file = "pate_check_corr.pdf")
fviz_contrib(pca.env, choice="var", axes = 1:2, title="H. patellifera (Asia vs SE)")
dev.off()
