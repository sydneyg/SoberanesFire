#May 16, 2019

#Make the figure of bacterial and fungal richness in same panel

#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Desktop/Big_Sur_Project/R_Analysis")

#Packages to load 
library(tidyverse)
library(ggpubr)


#read in data
bacteriadata <- read.csv( "bacteriarichnessdata.csv", row.names=1)
fungidata <- read.csv("fungalrichnessdata.csv", row.names=1)

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

fungidata$Burn <- factor(fungidata$Burn, levels=c("Unburned","Burned"))

levels(fungidata$Burn)

fungidata$Burn <- relevel(fungidata$Burn, ref = "Unburned", "Burned")

fungidata$Fire <- as.factor(fungidata$Fire)

levels(fungidata$Fire)

fungidata$Fire <- relevel(fungidata$Fire, ref = "Prefire", "Postfire")



A<-ggplot(bacteriadata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,0,0),
               vjust=c(1.5,3,1.5,0,5,15),
               label = c("a","a","a","","b","b"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Bacteria") +
  theme_bw() + #make black and white
  ylab("Mean Bacterial Richness") +  #change yaxis label
  theme(legend.position = "NULL", #put legend under graph
        text = element_text(size=20),
        legend.title = element_blank(),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=20, angle=90,hjust=0.5),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled
A

A %>%
  ggexport(filename = "bacterial_richnesspreandpostfire_MiSeq_meanandSE_one.pdf")

B<-ggplot(fungidata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,-0.5,0,0),
               vjust=c(5.5,4.5,8,8.25,7.5,3),
               label = c("a","a","a","a","b","b"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Fungi") +
  theme_bw() + #make black and white
  ylab("Mean Fungal Richness") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=20),
        legend.title = element_blank(),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=20, angle=90,hjust=0.5),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled
B

B %>%
  ggexport(filename = "fungalrichnesspreandpostfire_MiSeq_meanandSE_one.pdf")



#arrange using ggpubr package  (for example to put bacteria and fungal plots side by side
richnessfigures <- ggarrange(A,B, ncol = 1, nrow = 2)

richnessfigures

#other option is to bind together fungi and bacteria together, add a column for bacteria v fungi, then make the plot in ggplot 2 and facet wrap by fungi vs bacteria

pdf("BacterialANDFungalrichness_meanandSE_one_column_2.pdf", height=15, width=7)
richnessfigures
dev.off()
