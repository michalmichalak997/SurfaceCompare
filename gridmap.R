library(dplyr)
library(ggplot2)
library(ggpubr)


tab0<-read.csv("quaternary_oxford.txt", sep=";", dec=".")
tab1<-read.csv("oxford_callovian.txt", sep=";", dec=".")
tab2<-read.csv("czestochowa_koscieliska.txt", sep=";", dec=".")

tab0<-dplyr::filter(tab0, DOC<0.97)
tab1<-dplyr::filter(tab1, DOC<0.97)
tab2<-dplyr::filter(tab2, DOC<0.97)

nrow(tab0)
nrow(tab1)
nrow(tab2)

table(tab0$Dip_ang==1.2)
table(tab1$Dip_ang==1.2)
table(tab2$Dip_ang==1.2)

tab0$classification<-rep(4, nrow(tab0))
head(tab0)
tab0[tab0$Dip_ang<1.2 & tab0$Dip_dir>20 & tab0$Dip_dir<70,]$classification<-1
tab0[tab0$Dip_ang>1.2 & tab0$Dip_dir>20 & tab0$Dip_dir<70,]$classification<-2
tab0[tab0$Dip_dir>=70 & tab0$Dip_dir<=225,]$classification<-3
nrow(tab0)
head(tab0)

tab1$classification<-rep(4, nrow(tab1))
head(tab1)
tab1[tab1$Dip_ang<1.2 & tab1$Dip_dir>20 & tab1$Dip_dir<70,]$classification<-1
tab1[tab1$Dip_ang>1.2 & tab1$Dip_dir>20 & tab1$Dip_dir<70,]$classification<-2
tab1[tab1$Dip_dir>=70 & tab1$Dip_dir<=225,]$classification<-3
nrow(tab1)
head(tab1)


tab2$classification<-rep(4, nrow(tab2))
head(tab2)
tab2[tab2$Dip_ang<1.2 & tab2$Dip_dir>20 & tab2$Dip_dir<70,]$classification<-1
tab2[tab2$Dip_ang>1.2 & tab2$Dip_dir>20 & tab2$Dip_dir<70,]$classification<-2
tab2[tab2$Dip_dir>=70 & tab2$Dip_dir<=225,]$classification<-3
nrow(tab2)
head(tab2)


write.table(x=tab0, file="tab0_classified.txt", col.names = T, row.names = F, sep=";", dec=".")
write.table(x=tab1, file="tab1_classified.txt", col.names = T, row.names = F, sep=";", dec=".")
write.table(x=tab2, file="tab2_classified.txt", col.names = T, row.names = F, sep=";", dec=".")

#Gridding: first surface

gridpliknazwa<-"grid_locate"
nclusterpliknazwa<-"tab1_classified"
gridded <- read.table(paste0(gridpliknazwa, ".txt"), header=TRUE, sep = ";")
head(gridded)
nbottom <- read.table(paste0(nclusterpliknazwa,".txt"), header=TRUE, sep = ";")
table(nbottom$classification)
head(nbottom)
nmerged <- merge(x = gridded, y = nbottom, by = c("IDT1", "IDT2", "IDT3"), all.y = TRUE)
head(nmerged)

nl_pocz<- nrow(nmerged)
nmerged <- filter(nmerged, !is.na(px) )
nl_filtr <- nrow(nmerged)
nl_filtr/nl_pocz*100
colnames(nmerged)
ekspert_pal<- c("green", "black", "blue", "magenta")


max_x<-max(max(nbottom$X1 ), max(nbottom$X2 ), max(nbottom$X3 ))
min_x<-min(min(nbottom$X1 ), min(nbottom$X2 ), min(nbottom$X3 ))
max_y<-max(max(nbottom$Y1 ), max(nbottom$Y2 ), max(nbottom$Y3 ))
min_y<-min(min(nbottom$Y1 ), min(nbottom$Y2 ), min(nbottom$Y3 ))
coeff <- (max_x-min_x)/(max_y-min_y)
coeff

gekspert1<-ggplot(nmerged, aes(x=py, y=px, col=factor(classification) ))+
  geom_point(size=1.0)+
  scale_color_manual("Cluster", values =  ekspert_pal)+ggtitle("Expert-guided partition \n Contact separating Oxfordian sediments from Callovian sediments")


#Gridding: Second surface
gridpliknazwa<-"grid_locate"
nclusterpliknazwa<-"tab2_classified"
gridded <- read.table(paste0(gridpliknazwa, ".txt"), header=TRUE, sep = ";")
head(gridded)
nbottom <- read.table(paste0(nclusterpliknazwa,".txt"), header=TRUE, sep = ";")
table(nbottom$classification)
head(nbottom)
nmerged <- merge(x = gridded, y = nbottom, by = c("IDT1", "IDT2", "IDT3"), all.y = TRUE)
head(nmerged)

nl_pocz<- nrow(nmerged)
nmerged <- filter(nmerged, !is.na(px) )
nl_filtr <- nrow(nmerged)
nl_filtr/nl_pocz*100
colnames(nmerged)
ekspert_pal<- c("green", "black", "blue", "magenta")


max_x<-max(max(nbottom$X1 ), max(nbottom$X2 ), max(nbottom$X3 ))
min_x<-min(min(nbottom$X1 ), min(nbottom$X2 ), min(nbottom$X3 ))
max_y<-max(max(nbottom$Y1 ), max(nbottom$Y2 ), max(nbottom$Y3 ))
min_y<-min(min(nbottom$Y1 ), min(nbottom$Y2 ), min(nbottom$Y3 ))
coeff <- (max_x-min_x)/(max_y-min_y)
coeff

gekspert2<-ggplot(nmerged, aes(x=py, y=px, col=factor(classification) ))+
  geom_point(size=1.0)+
  scale_color_manual("Cluster", values =  ekspert_pal)+ggtitle("Expert-guided partition \n Contact separating Częstochowa Ore-Bearing Clay Fm \n from Kościeliska Beds")

gekspert2

#Gridding common plot

tiff(paste0(nclusterpliknazwa,"_", Sys.Date(), "_","common07072022.tiff"), units="in", width = 10*coeff, height=10, res=300)
ggarrange( gekspert1, gekspert2, ncol=1, nrow=2, labels=c("A", "B"), common.legend = TRUE)
dev.off()

