###I LOADED THE CSV OF THE SPECIES BEFORE PLOTTING

library(tidyverse)

pate_e<-pate_all[pate_all$status == 'E',]
pate_asia<-pate_asiavsse[pate_asiavse$status == 'N',]
pate_se<-pate_asiavsse[pate_asiavsse$status == 'SE',]

pdf(file="final_map_pate.pdf")
ggplot() +
geom_map(
data = world, map = world,
aes(long, lat, map_id = region),
color = "grey", fill= "white"
)+
geom_point(
data = pate_e,
aes(decimalLongitude, decimalLatitude, color = "green"),
alpha = 1
) +
geom_point(
data = pate_asia,
aes(decimalLongitude, decimalLatitude, color = "yellow"),
alpha = 1
)+
geom_point(
data = pate_se,
aes(decimalLongitude, decimalLatitude, color = "red"),
alpha = 1
)
dev.off()
