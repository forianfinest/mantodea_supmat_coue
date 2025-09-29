###BASICALLY, FOR DOING THIS I LOADED ALL THE DATASET, MERGED THEM TOGETHER, CONVERTED THE FINAL BIG TABLE AS A MATRIX AND CALCULATED VIFS

library(usdm)
allagain<-rbind(batt,sinensis,angui,pate)
allagain<-allagain[4:22]
allagain<-as.matrix(allagain)
yah<-vifstep(allagain, th=5)
write.csv(yah@results, "finalvif.csv")
