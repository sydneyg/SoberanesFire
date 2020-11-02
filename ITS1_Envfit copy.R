########Attempting to Use Envfit on the Soberanes Fire ITS1 Sequences based on Sydney's Code#########

#####Source of code is sydney_code.R#########

#set working directory
setwd("~/Desktop/Big_Sur_Project/R_Analysis")

#Reset R's Brain
rm(list=ls())


install.packages("ggrepel")
library(ggrepel)
#make the otu file look like Sydney's
otus<-read.csv("otu_ITS.e12089.csv", row.names=1)

transformed_otus<-t(otus)

write.csv(transformed_otus, "otu_ITS.e12089_rarified_with_tax.csv")

#####Now will edit this transformed table in Excel to add taxonomy info onto it
#####So that it looks like Sydney's "otu_table_pyro_e3950.csv"

#read the fixed file back in under the "otus" name 

#####STARTING WITH BURNED
otus<-read.csv("otu_ITS_rarified_with_tax_burned.csv")
dim(otus)
colnames(otus)
#saving taxonomic info in a different object
otu_taxonomy<-otus[,c(1,47)]

colnames(otu_taxonomy)
# i am using melt to reorganize data, this is more my own preference
install.packages("reshape2")
library(reshape2)
otus_melt<-melt(otus[,-c(47)], id="X")

#now casting it so that rows are samples and columns are otus
otu_cast<-dcast(otus_melt, variable~X, variable.var="value",fun.aggregate=mean,fill=0)

##Looks like for mapping file, the metadata file should work

# reading in mapping file
sample_info<-read.csv("SoberanesFiremetadata_ITS_burned.csv")

#merging sample info and otu counts

merged_dataset<-merge(sample_info, otu_cast,by.x="Sample",by.y="variable")

levels(merged_dataset$Fire)
merged_dataset$Fire <- relevel(merged_dataset$Fire, ref = "Prefire", "Postfire")

write.csv(merged_dataset, "envfit_merged_dataset_burned.csv")

merged_dataset <- read.csv("envfit_merged_dataset_burned.csv", row.names =1)

# ok, now I'm going to do an NMDS with relative abundances...there is some nesting going on here, I have some better explanations in a walkthrough https://github.com/ryanjw/COBS_all_metaG_analysis/tree/master/walkthrough

library(vegan)

nmds_object<-metaMDS(decostand(merged_dataset[,-c(1:45)], "total"), k=3,autotransform=FALSE)


#extracting scores for plotting

nmds<-data.frame(merged_dataset[,1:45],scores(nmds_object))

write.csv(nmds, "ITS1_envfit_nmds.csv")

str(nmds)

nmds$Fire <- as.factor(nmds$Fire)

levels(nmds$Fire)
nmds$Fire <- relevel(nmds$Fire, ref = "Prefire", "Postfire")

###Sydney's original data classified by $Pyrocosms, I am changing this. I think mine should probably
###based on by $Plot

hull_data<-data.frame()
for(i in 1:length(unique(nmds$Plot))){
  new_row<-nmds[nmds$Plot==as.vector(unique(nmds$Plot)[i]),][chull(nmds[nmds$Plot==as.vector(unique(nmds$Plot)[i]),c("NMDS1","NMDS2")]),]
  hull_data<-rbind(hull_data,new_row)
}

levels(hull_data$Fire)
hull_data$Fire <- relevel(hull_data$Fire, ref = "Prefire", "Postfire")

# I am going to plot this in ggplot2 for simplicity (for me), it should be easy to change this code up for whatever you'd like to do
#install.packages("ggplot2")
library(ggplot2)

ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4)+theme_bw()+theme(aspect.ratio=1)

# Ok, now we are getting to what you really want.  Let's summarize your data by the order just to give an example, let me know if you need help changin it up

#first I am doing some messy stuff to parse out the genera from your list
otu_taxonomy$genera<-NA
for(i in 1:dim(otu_taxonomy)[1]){
  otu_taxonomy$genera[i]<-unlist(strsplit(as.vector(otu_taxonomy$taxon[i]),split="; ", fixed=TRUE))[6]
}
otu_taxonomy$genera<-as.factor(otu_taxonomy$genera)

head(otu_taxonomy)

# merging this info into the melted data frame
merged_melted_data<-merge(otus_melt,otu_taxonomy[,c(1,3)],by="X" )

dim(otus_melt)
dim(otu_taxonomy)
library(plyr)
merged_genera<-ddply(merged_melted_data, .(variable,genera),summarise, sums=sum(value))
genera_cast<-dcast(merged_genera, variable~genera, value.var="sums",fun.aggregate=mean,fill=0)

head(genera_cast)
summary(genera_cast)
dim(genera_cast)
colnames(genera_cast)
row.names(genera_cast)
write.csv(genera_cast, "genera_cast.csv")
genera_cast <- read.csv("genera_cast.csv", row.names=1)
dim(nmds)

#now using envfit and pulling out relevant info - genera
fitted<-envfit(nmds_object, genera_cast,permutations=999)
genera_fits<-data.frame(fitted$vectors$arrows,fitted$vectors$r,fitted$vectors$pvals)
genera_fits$genera<-rownames(genera_fits)
names(genera_fits)[3:4]<-c("r","pval")

# adjusting pvalues based on whatever you want
install.packages("fdrtool")
library(fdrtool)
?p.adjust
genera_fits$p.adj<-p.adjust(genera_fits$pval,"BH")
ggplot(genera_fits)+geom_point(aes(x=pval,y=p.adj))+geom_hline(yintercept=0.05)

write.csv(genera_fits, "genera_fits.csv")
#####Following Comment's from Sydney's Code, May NOT be true on Mine
# only one significant at a p.adj < 0.05, easy to plot! I will plot a few more just to show you how I would do it
genera_fits_sig_pval_0.01<-subset(genera_fits, p.adj < 0.06)

genera_fits_sig_pval_0.01<-subset(genera_fits, pval < 0.01)
genera_fits_sig_pval_0.01<-data.frame(genera_fits_sig_pval_0.01)
write.csv(genera_fits_sig_pval_0.01, "generafitsfungi.csv")
genera_fits_sig_pval_0.01 <-read.csv("generafitsfungi.csv",row.names=1)
genera_fits_sig_pval_0.01$x<-0
genera_fits_sig_pval_0.01$y<-0

# you may need to play around with jitter and size for the text to avoid over-plotting
pdf("Soberanes_ITS1_NMDS_genera_pval0.01_burned.pdf")
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4)+
  #scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
   #                 values=c("cornflowerblue", "tomato1")) +
  #scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
   #                  values=c("cornflowerblue", "tomato1")) +
  geom_text(data=genera_fits_sig_pval_0.01,aes(x=NMDS1,y=NMDS2,label=genera, colour=Phylum,),size=3, position=position_jitter(height=.17,width=.25))+
  geom_segment(data=genera_fits_sig_pval_0.01,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)
  
dev.off()

p <- ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,color=Fire),size=4)+
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                   values=c("cornflowerblue", "tomato1")) +
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) +
  theme_bw()
p

p+ geom_text(data=genera_fits_sig_pval_0.01,aes(x=NMDS1,y=NMDS2,label=genera, colour=Phylum),size=3, position=position_jitter(height=.17,width=.25))+
  geom_segment(data=genera_fits_sig_pval_0.01,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)




#I want to try to make one with p<0.001 (DID NOT WORK)
#p<0.01 is too messy to be of use

genera_fits_sig_pval_0.001<-subset(genera_fits, pval < 0.002)
genera_fits_sig_pval_0.001<-data.frame(genera_fits_sig_pval_0.001)
genera_fits_sig_pval_0.001$x<-0
genera_fits_sig_pval_0.001$y<-0

write.csv(genera_fits_sig_pval_0.001, "genera_fits_sig_pval_0.001.csv")

# you may need to play around with jitter and size for the text to avoid over-plotting
pdf("Soberanes_ITS1_NMDS_genera_pval0.005_Burned.pdf")
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4)+
  geom_text(data=genera_fits_sig_pval_0.005,aes(x=NMDS1,y=NMDS2,label=genera),size=3,position=position_jitter(height=.17,width=.25))+
  geom_segment(data=genera_fits_sig_pval_0.005,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)
dev.off()

#p<0.005 is too messy to be of use

genera_fits_sig_pval_0.002<-subset(genera_fits, pval < 0.002)
genera_fits_sig_pval_0.002<-data.frame(genera_fits_sig_pval_0.002)
genera_fits_sig_pval_0.002$x<-0
genera_fits_sig_pval_0.002$y<-0

# you may need to play around with jitter and size for the text to avoid over-plotting

#p<0.0011 is too messy to be of use

genera_fits_sig_pval_0.00101<-subset(genera_fits, pval < 0.00101)
genera_fits_sig_pval_0.00101<-data.frame(genera_fits_sig_pval_0.00101)
genera_fits_sig_pval_0.00101$x<-0
genera_fits_sig_pval_0.00101$y<-0

# you may need to play around with jitter and size for the text to avoid over-plotting
pdf("Soberanes_ITS1_NMDS_genera_pval0.00101_Burned.pdf")
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4)+
  geom_text(data=genera_fits_sig_pval_0.00101,aes(x=NMDS1,y=NMDS2,label=genera),size=3,position=position_dodge(width=0.5))+
  geom_segment(data=genera_fits_sig_pval_0.00101,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)
dev.off()
#p<0.00101 is too messy to be of use

genera_fits_sig_pval_0.001001<-subset(genera_fits, pval < 0.001001)
genera_fits_sig_pval_0.001001<-data.frame(genera_fits_sig_pval_0.001001)
genera_fits_sig_pval_0.001001$x<-0
genera_fits_sig_pval_0.001001$y<-0

# you may need to play around with jitter and size for the text to avoid over-plotting
pdf("Soberanes_ITS1_NMDS_genera_pval0.001001noadjust.pdf")
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Burn,group=Plot),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4)+
  geom_text(data=genera_fits_sig_pval_0.001001,aes(x=NMDS1,y=NMDS2,label=genera),size=3,position=position_jitter(height=.17,width=.25))+
  geom_segment(data=genera_fits_sig_pval_0.001001,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)
dev.off()

genera_fits_sig_pval_0.001001_top_10<-genera_fits_sig_pval_0.001001[1:10, ]

#Figure of top 10

pdf("Soberanes_ITS1_NMDS_genera_pval0.001001_top_10.pdf")
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Burn,group=Plot),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4)+
  geom_text(data=genera_fits_sig_pval_0.001001_top_10,aes(x=NMDS1,y=NMDS2,label=genera),size=3,position=position_jitter(height=.17,width=.25))+
  geom_segment(data=genera_fits_sig_pval_0.001001_top_10,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)
dev.off()

#Figure of top 5

genera_fits_sig_pval_0.001001_top_5<-genera_fits_sig_pval_0.001001[1:5, ]

pdf("Soberanes_ITS1_NMDS_genera_pval0.001001_top_5.pdf")
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Burn,group=Plot),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=4) +
  xlim(-1.5,1.5) +
  geom_text(data=genera_fits_sig_pval_0.001001_top_5,aes(x=NMDS1,y=NMDS2,label=genera),size=3,position=position_jitter(height=.17,width=.25))+
  geom_segment(data=genera_fits_sig_pval_0.001001_top_5,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1)
dev.off()

######Now to make the figure better looking####

write.csv(genera_fits_sig_pval_0.001001_top_10, "ITS1_Envfit_top_10.csv")
genera_fits_sig_pval_0.001001_top_10<- read.csv("ITS1_Envfit_top_10.csv")


levels(hull_data$Fire)
hull_data$Fire <- relevel(hull_data$Fire, ref = "Prefire", "Postfire")


set.seed(1)

pdf("Soberanes_ITS1_Envfit_Figure.pdf", height = 10, width = 12)
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=5)+
  geom_text_repel(data=genera_fits_sig_pval_0.001001_top_10,aes(x=NMDS1,y=NMDS2,label=genera),size=5, nudge_x = c(0.1, -0.2, 0.5, 0.105, 0.5, 0.1, 0.5, 0.1, 0.5, -0.2),
                  nudge_y = c(0.3, 0, -0.2, -0.05, 0, 0, 0, 0, 0, 0),  direction = "y") + 
  geom_segment(data=genera_fits_sig_pval_0.001001_top_10,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1) +
scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                 values=c("cornflowerblue", "tomato1")) +
scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) +
theme(legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
theme(axis.text.x= element_text(size = 18)) +
theme(axis.text.y = element_text(size = 18)) +
theme(axis.title.x = element_text(size = 18)) +
theme(axis.title.y = element_text(size = 18)) +
  xlim(-1.8,1.8)+
  ylim(-2,1.5)
dev.off()

##########Trying to make the figure look like the NMDS figures#####

genera_fits_sig_pval_0.001_top_10 <- read.csv("genera_fits_sig_pval_0.001_top10.csv")

set.seed(1)

str(nmds$Fire)


Envfit_fig <-ggplot()+
   geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,color=Fire),size=5, shape=17)+
  geom_text_repel(data=genera_fits_sig_pval_0.001_top_10,aes(x=NMDS1,y=NMDS2,label=genera),size=5, nudge_x = c(0.1, -0.5, 10, 0.105, 0.5, 0.5, 0.5, 0.1, 0.5, -0.2),
                  nudge_y = c(0.3, 0.05, -0.2, -0.05, 0, 0, 0, 0, 0, 0),  direction = "y") + 
  geom_segment(data=genera_fits_sig_pval_0.001_top_10,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1) +
 # geom_path(data=nmds, aes(x=NMDS1,y=NMDS2, color=Fire)) +
#  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
#                    values=c("cornflowerblue", "tomato1")) +
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) +
  theme_bw()+
  theme(legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position = "none") +
  theme(axis.text.x= element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  stat_ellipse(data=nmds, aes(x=NMDS1,y=NMDS2, color=Fire)) + 
  xlim(-2,2)+
  ylim(-2,2)


pdf("Fungi_Envfit_fig.pdf", height=8, width=8)
Envfit_fig
dev.off()
