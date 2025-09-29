###I MANUALLY MERGED THE TWO DATASETS AND THEN RUN THE ANALYSES, AS DONE WITH THE INTRA-SPECIFIC ONES AND THE FULL DATASET

sinesis_angui<-read.csv(file="limited_anguisinensis.csv",header=TRUE)


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

"map_overlap" <- function(y, status){
y_E <- y[status == "sinensis"]
y_N <- y[status == "angustipennis"]
dens_E <- density(y_E)
dens_N <- density(y_N)
xlim <- range(dens_E$x, dens_N$x) #x axis range
ylim <- range(0, dens_E$y, dens_N$y)
ylim[2] <- ylim[2] + 0.1
hist(y_E, proba = T,  xlim = xlim, ylim = ylim, col= 'white', density = 10, angle = 135, main = "sinensis vs angustipennis", xlab = "Between-Class PCA 1", ylab = "Density", cex.lab=1.5, cex.axis=1.5)
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
pca.env <- dudi.pca(sinesis_angui[,c(8,9,13,14,15,16,18,21,22)], scannf = F, nf = 10) ##so can first make a pca of all, look at eigen values ect.
status <- as.factor(sinesis_angui$status)
status <- as.factor(status)
bet1 <- bca(pca.env, status, scan = FALSE, nf = 9)
#bet1$ratio #between-class ratio
correlation_to_pca.env <- as.data.frame (bet1$as)#see the correlation with the normal pca for the axis.
correlation_to_pca <- t(correlation_to_pca.env) #this is a data frame of the contribution of each origional pca axis to the between-pca axis
correlation_to_bio_variables <- as.data.frame (bet1$co) #to see the contribution of each worldclim BIO variable
correlation_to_bio_variables <- t(correlation_to_bio_variables)
#ecospat.plot.contrib(bet1$co)
## save interesting values
between_class <- as.data.frame(bet1$ls) #data frame of the first axis of between class analysis
data <-  cbind(status, between_class) #binding the two together
#the scores for the whole environment
scores.globclim <- data$CS1
#scores for the presences
scores.sp.ex <- data$CS1[which(data$status=='sinensis')]
scores.sp.nat <-data$CS1[which(data$status=='angustipennis')]
#scores for the whole study region, same as presences as no considering whole study area
scores.clim.ex <- data$CS1[which(data$status=='sinensis')]
scores.clim.nat <-data$CS1[which(data$status=='angustipennis')]
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
row <- cbind.data.frame(D.overlap, eq.test$p.D)
names(row) <- c('D.overlap', 'eq.test')
write.csv(row, "limited_overlap_sinensis_angustipennis.csv")
pdf(file="limited_overla_sinensis_angui_overlap.pdf")
map_overlap(bet1$ls[,1], status)
dev.off()

###PCA PLOTS

pdf(file = "anguivssinensis_trimm_contr.pdf")
fviz_contrib(pca.env, choice="var", axes = 1:2, title="Tenodera in US (Anderson 2019)")
dev.off()

pdf(file = "anguivssinensis_trimm_pca_2ax.pdf")
fviz_pca_var(pca.env,
col.var = "contrib", # Color by contributions to the PC
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
repel = TRUE     # Avoid text overlapping
)
dev.off()

###IF EVERYTHING SEEMS TOO FAMILIAR, IT IS BECAUSE I JUST CHANGED THE DATASET AND RUN THE SAME ANALYSES AS THE FULL ONE
