

###Attempting to make a single figure for Bacterial and Fungal NMDS#####

#set working directory

setwd("~/Desktop/Big_Sur_Project/R_Analysis")


#Reset R's Brain
rm(list=ls())


#Packages to load 
library(tidyverse)
library(ggpubr)
library(EcolUtils)
library(SPECIES)
library(vegan)
library(BiodiversityR)
library(scales)
library("nlme")
library("multcomp")
library("stringr")
library("plyr")

#read in data
bacteriadata <- read.csv( "16S_betadiversity_figure.csv", row.names=1)
fungidata <- read.csv("ITS_betadiversity.csv", row.names=1)

#add in domain column
bacteriadata$Domain <- rep("Bacteria",nrow(bacteriadata)) 
fungidata$Domain <- rep("Fungi",nrow(fungidata )) 

#bindtogether
microbedata <- rbind(bacteriadata,fungidata)



betadiversitydata_Burnplots <- microbedata[which(microbedata$Burn=="Burned"), ]
betadiversitydata_unburnplot <- microbedata[which(microbedata$Burn=="Unburned"), ]

#Change factors for Bacteria and Fungi Datasets
# for Bacteria

str(bacteriadata)

bacteriadata$Burn<- as.factor(bacteriadata$Burn)

levels(bacteriadata$Burn)

bacteriadata$Burn <- relevel(bacteriadata$Burn, ref = "Unburned", "Burned")

bacteriadata$Fire <- as.factor(bacteriadata$Fire)

levels(bacteriadata$Fire)

bacteriadata$Fire <- relevel(bacteriadata$Fire, ref = "Prefire", "Postfire")
#for Fungi

str(fungidata)

fungidata$Burn <- as.factor(fungidata$Burn)

levels(fungidata$Burn)

fungidata$Burn <- relevel(fungidata$Burn, ref = "Unburned", "Burned")

fungidata$Fire <- as.factor(fungidata$Fire)

levels(fungidata$Fire)

fungidata$Fire <- relevel(fungidata$Fire, ref = "Prefire", "Postfire")

#NMDS for bacteria 
NMDS1_bacteria <- ggplot(bacteriadata, aes(x=NMDS1, y=NMDS2, col=Fire, shape=Burn, group=Plot)) +
  scale_shape_manual(values=c(16, 17)) +
   geom_point(size=2) +
  theme_bw() +
  #labs(title = "Burned Plots Bray-Curtis Dissimilarity")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
 labs(title="Bacteria") +
 theme(plot.title = element_text(size =18)) +
theme(legend.position = "none", #put legend under graph,
        legend.title = element_blank(),
       strip.text.x = element_text(size=15), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=10,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Postfire","Prefire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1"))+ #add in manual colors for 
guides(color=FALSE, shape=FALSE)

total_bacteria_fig <- NMDS1_bacteria

pdf("Total_Bacteria_NMDS.pdf", height = 8, width = 8)
total_bacteria_fig
dev.off()

bacteria_figs <- NMDS1_bacteria + facet_wrap(~Burn, scales = "free")

bacteria_figs
#NMDS for fungi
NMDS2_fungi <- ggplot(fungidata, aes(x=NMDS1, y=NMDS2, col=Fire, shape=Burn, group=Plot)) +
  scale_shape_manual(values=c(16, 17)) +
  geom_point(size=2) +
  theme_bw() +
 # labs(title = "Unburned Plots Bray-Curtis Dissimilarity")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  labs(title="Fungi") +
  theme(plot.title = element_text(size =18)) +
   theme(legend.position = "bottom", #put legend under graph
        legend.title = element_blank(),
        strip.text.x = element_text(size=15), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=10,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=15), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Pre-Fire","Post-Fire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) #add in manual colors

total_fungi_fig <- NMDS2_fungi

pdf("Total_Fungi_NMDS.pdf", height = 8, width = 8)
total_fungi_fig
dev.off()

fungi_figs <- NMDS2_fungi + facet_wrap(~Burn, scales = "free")
fungi_figs
#arrange using ggpubr package  (for example to put fungi and fungal plots side by side
NMDSfigures <- ggarrange(bacteria_figs,fungi_figs, ncol = 1, nrow = 2)

NMDSfigures

  pdf("Combined_NMDS_Miseq.pdf", height=10, width=8)
NMDSfigures
dev.off()


library("patchwork")

NMDS_Total <- (total_bacteria_fig | total_fungi_fig) +  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom', legend.text = element_text(size=18))

NMDS_Total

pdf("NMDS_Total.pdf", height=8, width=12)
NMDS_Total
dev.off()


