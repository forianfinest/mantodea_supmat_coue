####NOTE: ON INATURALIST, THE TRANSCAUCASIAN AND THE GIANT ASIAN MANTISES ARE SPLIT

library(spThin)
tenui<-read.csv(file="tenuidentata_raw.csv", header=TRUE, sep=",")
View(tenui)
tenuit<-thin(
tenui,
lat.col = "decimalLatitude",
long.col = "decimalLongitude",
spec.col = "scientific_name",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/possible_pate/new_riel/all_new", #####CHOOSE OUTPUT FOLDER
out.base = "thinned_tenui_raw",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)

transcaucasica<-read.csv(file="trancaucasica_raw.csv", header=TRUE, sep=",")
View(transcaucasica)
transcaut<-thin(
transcaucasica,
lat.col = "decimalLatitude",
long.col = "decimalLongitude",
spec.col = "scientific_name",
thin.par=10,
reps=100,
locs.thinned.list.return = FALSE,
write.files = TRUE,
max.files = 1,
out.dir="/run/media/brcuser/MDVBook/possible_pate/new_riel/all_new", #####CHOOSE OUTPUT FOLDER
out.base = "thinned_transcaucasica_raw",
write.log.file = TRUE,
log.file = "final",
verbose = TRUE
)

###LOAD SOME LIBRARIES

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

###SET POINTS
transcauc<-read.csv(file="thinned_transcaucasica_raw_thin1.csv",header=TRUE,sep = ",")
transcauc<-na.omit(transcauc)
coordinates(transcauc) <- c("decimalLongitude", "decimalLatitude")
transcaucsp<-SpatialPoints(transcauc)
tenui<-read.csv(file="thinned_tenui_raw_thin1.csv",header=TRUE,sep = ",")
tenui<-na.omit(tenui)
coordinates(tenui) <- c("decimalLongitude", "decimalLatitude")
tenuisp<-SpatialPoints(tenui)

###EXTRACT VARIABLES
lista_variables2 <-list.files(pattern='*.tif')
datafiles <- Sys.glob("*.tif")
a <- stack(lista_variables2)
transcaucb<-raster::extract(a,transcaucsp)
transcauc_c<-cbind(transcauc,transcaucb)
transcauc_c<-na.omit(transcauc_c)
tenuib<-raster::extract(a,tenuisp)
tenuic<-cbind(tenui,tenuib)
tenuic<-na.omit(tenuic)

###GIVEN THAT I AM INTERESTED IN ANALYZING THE NATURAL POPULATIONS, I CROP

transcauc_ext<-extent(29.858549, 93.928067,20.626158, 53.612592)
transcauc_og<-crop(transcaucsp,transcauc_ext)
transcauc_ogvalues<-raster::extract(a, transcauc_og)
transcauc_ogvalues<-na.omit(transcauc_ogvalues)
tenui_ext<-extent(29.858549, 98.050660, 4.126667, 53.612592)
tenui_og<-crop(tenuisp,tenui_ext)
tenui_ogvalues<-raster::extract(a, tenui_og)
tenui_ogvalues<-na.omit(tenui_ogvalues)
write.csv(tenui_ogvalues, "tenui_original.csv")
write.csv(tenui_ogvalues, "transcauc_original.csv")


###I THEN MERGED THE TWO TOGETHER MANUALLY

###I ALSO MADE THE CSV FOR THE TWO SPECIES, TO EDIT FOR INTRASPECIES ANALYSES
write.csv(tenuic, "tenui.csv")
write.csv(transcauc_c, "transcauc.csv")

###NOW, LET'S ADD THE FUNCTIONS

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
col_E <- "green" 
col_N <- "yellow"
"map_overlap" <- function(y, status){
y_E <- y[status == "tenui"]
y_N <- y[status == "trans"]
dens_E <- density(y_E)
dens_N <- density(y_N)
xlim <- range(dens_E$x, dens_N$x) #x axis range
ylim <- range(0, dens_E$y, dens_N$y)
ylim[2] <- ylim[2] + 0.1
hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = "tenuidentata vs transcaucasica", xlab = "Between-Class PCA 1", ylab = "Density", cex.lab=1.5, cex.axis=1.5)
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
check<-read.csv(file="transvstenui.csv",header=TRUE,sep = ",")
status <- as.factor(check$status)
pca.env <- dudi.pca(check[,c(8,9,13,14,15,16,18,21,22)], scannf = F, nf = 10) ##so can first make a pca of all, look at eigen values ect
status <- as.factor(status)
bet1 <- bca(pca.env, status, scan = FALSE, nf = 9)
#bet1$ratio #between-class ratio
## explore correlaitons to each bio variable to axis
correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis.
correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
correlation_to_bio_variables <- t(correlation_to_bio_variables)
between_class <- as.data.frame(bet1$ls) #data frame of the first axis of between class analysis
data <-  cbind(status, between_class) #binding the two together
#the scores for the whole environment
scores.globclim <- data$CS1
#scores for the presences
scores.sp.ex <- data$CS1[which(data$status=='tenui')]
scores.sp.nat <-data$CS1[which(data$status=='trans')]
#scores for the whole study region, same as presences as no considering whole study area
scores.clim.ex <- data$CS1[which(data$status=='tenui')]
scores.clim.nat <-data$CS1[which(data$status=='trans')]
grid.clim.nat <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1= scores.sp.nat, sp=scores.sp.nat, R=100, th.sp=0)
#glob= all background dataframe, glob1= of the native/invasive envinoment range of the species (for each one), envirnomental variable for occurances only,
grid.clim.inv <- ecospat.grid.clim.dyn(glob=scores.globclim, glob1=scores.sp.ex , sp=scores.sp.ex, R=100, th.sp=0)
## Overlap
# D overlap
D.overlap <- ecospat.niche.overlap (grid.clim.nat, grid.clim.inv, cor=F)$D #niche overlap.
eq.test <- NA
eq.test <- niche.equivilency(grid.clim.nat, grid.clim.inv, rep=1000, alternative = "lower", ncores = 1)
View(eq.test)
row <- cbind.data.frame(D.overlap, eq.test$p.D)
names(row) <- c('D.overlap', 'eq.test')
write.csv(row, "transvstenui_metrics.csv")
pdf(file="greentenui_yellowtrans.pdf")
map_overlap(bet1$ls[,1], status)
dev.off()

###PCA PLOTS

pdf(file = "transvstenui_pca_2ax.pdf")
fviz_pca_var(pca.env,
col.var = "contrib", # Color by contributions to the PC
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE     # Avoid text overlapping
)
dev.off()
pdf(file = "transvstenui_corr.pdf")
fviz_contrib(pca.env, choice="var", axes = 1:2, title="trans vs tenui")
dev.off()
