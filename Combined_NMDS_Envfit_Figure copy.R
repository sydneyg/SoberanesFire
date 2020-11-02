#This script is for combining 4 ggplots into one figure to create a figure of
#the Soberanes unburned NMDS plots and the Envfit/NMDS for the burned plots.

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
library(ggrepel)


###Bacteria Unburned NMDS

#read in data
bacteriadata <- read.csv( "16S_betadiversity_figure.csv", row.names=1)

#Seperate out the Unburned plot info

bacteriadata_unburnplot <- bacteriadata[which(bacteriadata$Burn=="Unburned"), ]

#Fix factors for Bacteria Datasets

str(bacteriadata_unburnplot)

bacteriadata_unburnplot$Burn<- as.factor(bacteriadata_unburnplot$Burn)

levels(bacteriadata_unburnplot$Burn)

levels(bacteriadata_unburnplot$Burn) <- c(levels(bacteriadata_unburnplot$Burn), "Burned")

levels(bacteriadata_unburnplot$Burn)

bacteriadata_unburnplot$Fire <- as.factor(bacteriadata_unburnplot$Fire)

levels(bacteriadata_unburnplot$Fire)

bacteriadata_unburnplot$Fire <- relevel(bacteriadata_unburnplot$Fire, ref = "Prefire", "Postfire")

#NMDS for bacteria 
NMDS1_bacteria <- ggplot(bacteriadata_unburnplot, aes(x=NMDS1, y=NMDS2, col=Fire, shape=Burn, group=Fire)) +
  scale_shape_manual(values=c(16, 17)) +
  geom_point(size=5) +
  theme_bw() +
  #labs(title = "Burned Plots Bray-Curtis Dissimilarity")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  labs(title="(A)   Unburned Bacteria") +
  theme(plot.title = element_text(size =18)) +
  theme(legend.position = "bottom", #put legend under graph,
        legend.title = element_blank(),
        strip.text.x = element_text(size=15), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=10),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=15), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Pre-fire","Post-fire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) +
  guides(color=FALSE)


NMDS1_bacteria

#### Envfit for Bacteria Burned Plots

nmds_bac<- read.csv("16S_envfit_nmds.csv", row.names=1, check.names= FALSE)

nmds_bac$Fire <- as.factor(nmds_bac$Fire)

levels(nmds_bac$Fire)
nmds_bac$Fire <- relevel(nmds_bac$Fire, ref = "Prefire", "Postfire")

phylum_fits_sig_pval_0.001 <-read.csv("phylumfits16S.csv",row.names=1)
phylum_fits_sig_pval_0.001$x<-0
phylum_fits_sig_pval_0.001$y<-0


set.seed(1)

Bac_Envfit <-ggplot()+
  geom_point(data=nmds_bac,aes(x=NMDS1,y=NMDS2,shape=Burn,color=Fire),size=5)+
  scale_shape_manual(values=17) +
  stat_ellipse(data=nmds_bac, aes(x=NMDS1,y=NMDS2, color=Fire)) + 
  geom_text_repel(data=phylum_fits_sig_pval_0.001,aes(x=NMDS1,y=NMDS2,label=phylum),size=5,fontface = "italic", nudge_x = c(1, -0.3, 0.3, 0.5, 1, -0.2, 0.4,0.3, 0.3, 0.3, 0.1), direction = "y")+
  geom_segment(data=phylum_fits_sig_pval_0.001,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1) +
  # geom_path(data=nmds, aes(x=NMDS1,y=NMDS2, color=Fire)) +
  #  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
  #                    values=c("cornflowerblue", "tomato1")) +
  scale_color_manual(labels=c("Pre-fire","Post-fire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) +
  theme_bw()+
  labs(title="(B)   Burned Bacteria") +
  theme(plot.title = element_text(size =18)) +
  theme(legend.position = "bottom", #put legend under graph
        legend.title = element_blank(),
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=15))+
  theme(axis.text.x= element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-2,2)+
  ylim(-2,2)+
  guides(color = guide_legend(order = 2),
         shape = guide_legend(order = 1))

Bac_Envfit

#### NMDS for Fungi Unburned Plot

#read in data

fungidata <- read.csv("ITS_betadiversity.csv", row.names=1)

#Seperate out the Unburned plot info

fungidata_unburnplot <- fungidata[which(fungidata$Burn=="Unburned"), ]


#Fix factors for Fungi Datasets

str(fungidata_unburnplot)

fungidata_unburnplot$Burn <- as.factor(fungidata_unburnplot$Burn)

levels(fungidata_unburnplot$Burn)

fungidata_unburnplot$Burn <- relevel(fungidata_unburnplot$Burn, ref = "Unburned", "Burned")

fungidata_unburnplot$Fire <- as.factor(fungidata_unburnplot$Fire)

levels(fungidata_unburnplot$Fire)

fungidata_unburnplot$Fire <- relevel(fungidata_unburnplot$Fire, ref = "Prefire", "Postfire")

#NMDS for fungi
NMDS2_fungi <- ggplot(fungidata_unburnplot, aes(x=NMDS1, y=NMDS2, col=Fire, shape=Burn, group=Fire)) +
  scale_shape_manual(values=16) +
  geom_point(size=5) +
  theme_bw() +
  # labs(title = "Unburned Plots Bray-Curtis Dissimilarity")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  labs(title="(C)   Unburned Fungi") +
  theme(plot.title = element_text(size =18)) +
  theme(legend.position = "none", #put legend under graph
        legend.title = element_blank(),
        strip.text.x = element_text(size=15), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=10),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=15), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Pre-Fire","Post-Fire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1"))+
  guides(color=FALSE, shape=FALSE)

NMDS2_fungi


####Envfit Figure for Fungi Burned Plots

#Read in the data files

nmds_fun <- read.csv("ITS1_envfit_nmds.csv", row.names=1)

nmds_fun$Fire <- as.factor(nmds_fun$Fire)

levels(nmds_fun$Fire)
nmds_fun$Fire <- relevel(nmds_fun$Fire, ref = "Prefire", "Postfire")


genera_fits_sig_pval_0.001_top_10 <- read.csv("genera_fits_sig_pval_0.001_top10.csv")

set.seed(1)

Fun_Envfit <-ggplot()+
  geom_point(data=nmds_fun,aes(x=NMDS1,y=NMDS2,color=Fire, shape=Burn),size=5)+
  scale_shape_manual(values=17) +
  stat_ellipse(data=nmds_fun, aes(x=NMDS1,y=NMDS2, color=Fire)) + 
  geom_text_repel(data=genera_fits_sig_pval_0.001_top_10,aes(x=NMDS1,y=NMDS2,label=genera),size=5,fontface = "italic", nudge_x = c(0.1, -0.5, 10, 0.105, 0.5, 0.5, 0.5, 0.1, 0.5, -0.2),
                  nudge_y = c(0.3, 0.05, -0.2, -0.05, 0, 0, 0, 0, 0, 0),  direction = "y") + 
  geom_segment(data=genera_fits_sig_pval_0.001_top_10,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1) +
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) +
  theme_bw()+
  labs(title="(D)   Burned Fungi") +
  theme(plot.title = element_text(size =18)) +
  theme(legend.position = "none", #put legend under graph
        legend.title = element_blank(),
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=15))+
  theme(axis.text.x= element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-2,2)+
  ylim(-2,2)+
  guides(color=FALSE, shape=FALSE)

Fun_Envfit


###Knitting them all together using Patchwork######

devtools::install_github("thomasp85/patchwork")

library("patchwork")

#Now to stich them together

Full_Fig <- ((NMDS1_bacteria + Bac_Envfit)/(NMDS2_fungi + Fun_Envfit)) +  plot_layout(guides = 'collect') & 
  theme(legend.position = 'bottom') #&
  #scale_shape_manual(limits = range(c(bacteriadata_unburnplot$Burn, nmds_bac$Burn,
   #                                   fungidata_unburnplot$Burn, nmds_fun$Burn)), values = 16,17)
Full_Fig


pdf("NMDS_Evfit_Combined.pdf", height=10, width=14)
Full_Fig
dev.off()

