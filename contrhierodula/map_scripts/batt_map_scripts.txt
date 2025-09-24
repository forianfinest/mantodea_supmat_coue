###I LOADED THE CSV OF THE SPECIES BEFORE PLOTTING

library(tidyverse)

batt_e<-batt_all[batt_all$status == 'E',]
tenui<-transvstenui[transvstenui$status == 'tenui',]
trans_test<-transvstenui[transvstenui$status == 'trans',]
pdf(file="final_map_batt.pdf")
ggplot() +
geom_map(
data = world, map = world,
aes(long, lat, map_id = region),
color = "grey", fill= "white"
)+
geom_point(
data = tenui,
aes(decimalLongitude, decimalLatitude, color = "green"),
alpha = 1
) +
geom_point(
data = trans_test,
aes(decimalLongitude, decimalLatitude, color = "yellow"),
alpha = 1
)+
geom_point(
data = batt_e,
aes(decimalLongitude, decimalLatitude, color = "red"),
alpha = 1
) +
coord_sf(xlim = c(-7, 98), ylim = c(3, 54), expand = FALSE)
dev.off()
