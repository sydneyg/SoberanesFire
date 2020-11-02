######Making Envfit of 16S Info Following Same Code from ITS1_Envfit.R############

#set working directory
setwd("~/Desktop/Big_Sur_Project/R_Analysis")

#Reset R's Brain
rm(list=ls())

#make the otu file look like Sydney's
otus<-read.csv("otu_16S.e10367.csv", row.names=1, check.names=FALSE)

transformed_otus<-t(otus)

write.csv(transformed_otus, "otu_16S_rarified_with_tax.csv")

#####Now will edit this transformed table in Excel to add taxonomy info onto it
#####So that it looks like Sydney's "otu_table_pyro_e3950.csv"

#read the fixed file back in under the "otus" name 

#####STARTING WITH BURNED
otus<-read.csv("otu_16S_rarified_with_tax_burned.csv", check.names=FALSE)
dim(otus)
colnames(otus)
#saving taxonomic info in a different object
otu_taxonomy<-otus[,c(1,50)]

colnames(otu_taxonomy)
# i am using melt to reorganize data, this is more my own preference
install.packages("reshape2")
library(reshape2)
otus_melt<-melt(otus[,-c(50)], id="X")

#now casting it so that rows are samples and columns are otus
otu_cast<-dcast(otus_melt, variable~X, variable.var="value",fun.aggregate=mean,fill=0)

##Looks like for mapping file, the metadata file should work

# reading in mapping file
sample_info<-read.csv("SoberanesFiremetadata_16S_burned.csv")

#merging sample info and otu counts

merged_dataset<-merge(sample_info, otu_cast,by.x="Sample",by.y="variable")

levels(merged_dataset$Fire)
merged_dataset$Fire <- relevel(merged_dataset$Fire, ref = "Prefire", "Postfire")

write.csv(merged_dataset, "envfit_merged_dataset_burned_16S.csv")

merged_dataset <- read.csv("envfit_merged_dataset_burned_16S.csv", row.names=1)


# ok, now I'm going to do an NMDS with relative abundances...there is some nesting going on here, I have some better explanations in a walkthrough https://github.com/ryanjw/COBS_all_metaG_analysis/tree/master/walkthrough

library(vegan)

nmds_object<-metaMDS(decostand(merged_dataset[,-c(1:48)], "total"), k=3,autotransform=FALSE)

#extracting scores for plotting

nmds<-data.frame(merged_dataset[,1:48],scores(nmds_object))

write.csv(nmds, "16S_envfit_nmds.csv")

nmds<- read.csv("16S_envfit_nmds.csv", row.names=1, check.names= FALSE)

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

#first I am doing some messy stuff to parse out the phylum from your list
otu_taxonomy$phylum<-NA
for(i in 1:dim(otu_taxonomy)[1]){
  otu_taxonomy$phylum[i]<-unlist(strsplit(as.vector(otu_taxonomy$taxon[i]),split="; ", fixed=TRUE))[2]
}
otu_taxonomy$phylum<-as.factor(otu_taxonomy$phylum)

head(otu_taxonomy)

# merging this info into the melted data frame
merged_melted_data<-merge(otus_melt,otu_taxonomy[,c(1,3)],by="X" )

dim(otus_melt)
dim(otu_taxonomy)
library(plyr)
merged_phylum<-ddply(merged_melted_data, .(variable,phylum),summarise, sums=sum(value))
phylum_cast<-dcast(merged_phylum, variable~phylum, value.var="sums",fun.aggregate=mean,fill=0)

#now using envfit and pulling out relevant info - phylum
fitted<-envfit(nmds_object, phylum_cast,permutations=9999)
phylum_fits<-data.frame(fitted$vectors$arrows,fitted$vectors$r,fitted$vectors$pvals)
phylum_fits$phylum<-rownames(phylum_fits)
names(phylum_fits)[3:4]<-c("r","pval")

# adjusting pvalues based on whatever you want
install.packages("fdrtool")
library(fdrtool)
?p.adjust
phylum_fits$p.adj<-p.adjust(phylum_fits$pval,"BH")
ggplot(phylum_fits)+geom_point(aes(x=pval,y=p.adj))+geom_hline(yintercept=0.05)

#Grabbing everyting with padj>0.06
phylum_fits_sig_pval_0.05<-subset(phylum_fits, p.adj < 0.05)

#phylum_fits_sig_pval_0.01<-subset(phylum_fits, pval < 0.01)
phylum_fits_sig_pval_0.05<-data.frame(phylum_fits_sig_pval_0.05)
write.csv(phylum_fits_sig_pval_0.05, "phylumfitsfungi.csv")
#phylum_fits_sig_pval_0.01 <-read.csv("phylumfitsfungi.csv",row.names=1)
phylum_fits_sig_pval_0.05$x<-0
phylum_fits_sig_pval_0.05$y<-0

#padj of 0.05 not strict enough

phylum_fits_sig_pval_0.001<-subset(phylum_fits, p.adj < 0.001)
phylum_fits_sig_pval_0.001<-data.frame(phylum_fits_sig_pval_0.001)
#write.csv(phylum_fits_sig_pval_0.01, "phylumfitsfungi.csv")
#phylum_fits_sig_pval_0.01 <-read.csv("phylumfitsfungi.csv",row.names=1)
phylum_fits_sig_pval_0.001$x<-0
phylum_fits_sig_pval_0.001$y<-0



#still too many hits, going to constrain the r value

phylum_fits_sig_pval_0.001<-subset(phylum_fits_sig_pval_0.001, r > 0.44)
phylum_fits_sig_pval_0.001<-data.frame(phylum_fits_sig_pval_0.001)
write.csv(phylum_fits_sig_pval_0.001, "phylumfits16S.csv")
phylum_fits_sig_pval_0.001 <-read.csv("phylumfits16S.csv",row.names=1)
phylum_fits_sig_pval_0.001$x<-0
phylum_fits_sig_pval_0.001$y<-0



# you may need to play around with jitter and size for the text to avoid over-plotting

#using ggrepel instead of jitter

library(ggrepel)

set.seed(1)

pdf("Soberanes_16S_Envfit_Figure.pdf", height = 10, width = 12)
ggplot()+
  geom_polygon(data=hull_data,aes(x=NMDS1,y=NMDS2,fill=Fire,group=Fire),alpha=0.3)+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,shape=Plot,colour=Fire),size=5)+
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) +
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) +
  geom_text_repel(data=phylum_fits_sig_pval_0.001,aes(x=NMDS1,y=NMDS2,label=phylum),size=5,nudge_x = c(1, -0.3, 0.3, 0.5, 1, -0.2, 0.1,0.3, 0.3, 0.1, 0.1), direction = "y")+
  geom_segment(data=phylum_fits_sig_pval_0.001,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1) +
theme(legend.title = element_text(size = 20),
      legend.text = element_text(size = 20)) +
  theme(axis.text.x= element_text(size = 18)) +
  theme(axis.text.y = element_text(size = 18)) +
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  xlim(-2,2)+
  ylim(-2,2)
dev.off()


levels(hull_data$Burn )
hull_data$Burn <- relevel(hull_data$Burn, ref = "Unburned", "Burned")

levels(hull_data$Fire)
hull_data$Fire <- relevel(hull_data$Fire, ref = "Prefire", "Postfire")






Envfit_fig <-ggplot()+
  geom_point(data=nmds,aes(x=NMDS1,y=NMDS2,color=Fire),size=5, shape=17)+
  geom_text_repel(data=phylum_fits_sig_pval_0.001,aes(x=NMDS1,y=NMDS2,label=phylum),size=5, nudge_x = c(1, -0.3, 0.3, 0.5, 1, -0.2, 0.1,0.3, 0.3, 0.1, 0.1), direction = "y")+
                 # nudge_y = c(0.3, 0.05, -0.2, -0.05, 0, 0, 0, 0, 0, 0),  direction = "y") + 
  geom_segment(data=phylum_fits_sig_pval_0.001,aes(x=x,y=y,xend=NMDS1,yend=NMDS2))+theme_bw()+theme(aspect.ratio=1) +
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

Envfit_fig

pdf("Bacteria_Envfit_fig.pdf", height=8, width=8)
Envfit_fig
dev.off()

