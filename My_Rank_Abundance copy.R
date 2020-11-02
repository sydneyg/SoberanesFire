#set working directory
setwd("~/Desktop/Big_Sur_Project/R_Analysis")

#Reset R's Brain
rm(list=ls())

#load libraries
library(vegan)
library(BiodiversityR)
library(tidyverse)
library(plyr)
#read in non rarefied otu table from bags
#take the file where I removed the zeroes and also fixed some of the names
#read in dataframes
otu_ITS <- read.csv("otu_ITS_with_tax.csv", row.names=1)
predictor <-read.csv("SoberanesFiremetadata_ITS.csv", row.names=1)

#####################################################################
#######Organize data
#######################################################################

dim(otu_ITS) 

#look at no DNA controls
noDNA <- which(colnames(otu_ITS) %in%c("ITS1NODNA","ITS1PCRNC"))

#make otu table of no DNA controls
OTU_noDNA <-otu_ITS[ ,noDNA]
colSums(OTU_noDNA)

#remove no DNA controls
otu_ITS_2 <-otu_ITS[ ,-noDNA]
otu_ITS_2

dim(otu_ITS_2)
#remove taxonomy info
otu_ITS2 <- otu_ITS_2[,-c(72:81)]
#make this data frame presence absence by converting everything greater than 0 to a 1
head(otu_ITS2)
colnames(otu_ITS2)
otu_ITS2[otu_ITS2 > 0] <- 1 
head(otu_ITS2)

#transform dataframe
otu.trans <- t(otu_ITS2)

## We will create a new data frame for the last two columns as it provides different set of information. We'll call it "taxon" data.frame
taxon <-otu_ITS_2[,c(72:81)]

colnames(taxon)

###Moving on without Doing above step, previously done in other pipeline####
##Discard all samples in dd that don't appear otu table
predictor2 <- predictor[row.names(predictor) %in% row.names(otu.trans ), ]
otu.trans2 <- otu.trans[row.names(otu.trans ) %in% row.names(predictor2), ]
nrow(predictor)
nrow(otu.trans2)
#check that names match
row.names(predictor2) == row.names(otu.trans2)

#####################################################################
#######Make rank abundance curves for bags BASED on FREQUENCY (number of samples present in)
#######################################################################

#USe BioDiversityR function rankabundance to make curve
RankAbun <- rankabundance(otu.trans)
RankAbun
rankabunplot(RankAbun , scale='abundance', addit=FALSE, specnames=c(1,2,3))
rankabunplot(RankAbun , scale='logabun', addit=FALSE, specnames=c(1:30), 
             srt=45, ylim=c(1,100))

names(predictor2)
rankabuncomp(RankAbun, y=predictor2, factor='Burn', scale='proportion', legend=FALSE)

write.csv(RankAbun, "ITS_rankabun.csv")
#get taxonomy information on rank abundance 
class(taxon)
class(taxon$Species)
taxon$Species<- as.character(taxon$Species)
row.names(RankAbun) == taxon$id
taxon$id <- as.character(taxon$id)
#check if they have the same nunber of rows
nrow(RankAbun)
nrow(taxon)
#make them match
match.names3 <- match(row.names(RankAbun), taxon$id)
taxon2 <- taxon[match.names3, ]
taxon2$id ==row.names(RankAbun) 

RankAbun.tax <- cbind(RankAbun,taxon2)
RankAbun.tax

RankAbun.tax_top35 <-RankAbun.tax[1:35, ] 

write.csv(RankAbun.tax_top35, "RankAbund.tax_ITS_top35_frequency.csv")

RankAbun.tax_top35 <- read.csv("RankAbund.tax_ITS_top35_frequency.csv")

rankabundance_fig <- ggplot(data=RankAbun.tax_top35 , aes(x=rank, y=abundance)) + geom_point()+ 
  theme_bw()+
  labs(x=" ", y="Fungal Number of samples") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
  theme(strip.text.x = element_text(size = 6, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        axis.title=element_text(size=16)) #make y axis label larger  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +

rankabundance_fig
#change x axis tick labels
#https://stackoverflow.com/questions/20529252/changing-x-axis-tick-labels-in-r-using-ggplot2

names(RankAbun.tax_top35)
head(RankAbun.tax_top35)
rankabundance_fig1 <- ggplot(data=RankAbun.tax_top35 , aes(x=reorder(Species, rank), y=abundance, fill=Phylum)) + geom_bar(stat = "identity")+
  theme_bw()+
  labs(x=" ", y="Fungal number of samples") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
  theme(strip.text.x = element_text(size = 6, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text.x=element_text(size=12,angle=70, hjust=1),  #change size angle and justification of x axis labels
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        axis.title=element_text(size=16)) #make y axis label larger  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +

rankabundance_fig1


pdf("RankAbundance_fungi_overall_frequency.pdf", height=5, width=8)
rankabundance_fig1
dev.off()

#####################################################################
#######separate by plots 
######################################################################
levels(predictor2$Plot)
Plot601 <- which(predictor2$Plot=="Plot601")
predictor2_601<- predictor2[Plot601, ]

Plot603 <- which(predictor2$Plot=="Plot603")
predictor2_603<- predictor2[Plot603, ]

Plot58 <- which(predictor2$Plot=="Plot58")
predictor2_58<- predictor2[Plot58, ]

#####################################################################
#######separate by burn
######################################################################
levels(predictor2_601$Fire)
Plot601_prefire <- which(predictor2_601$Fire=="Prefire")
predictor2_601_prefire <- predictor2_601[Plot601_prefire, ]

Plot601_postfire <- which(predictor2_601$Fire=="Postfire")
predictor2_601_postfire<- predictor2_601[Plot601_postfire, ]


#separate out otu tables from different seasons
## "match.names" gives you the order of rows in predictor that matches the order in otu. e.g. otu.trans[1, ] = predictor[120, ]. Then we can use this to rearrange the rows, so that they are in the same order and drop the ones that are not in otu.trans

#remove samples from otu table not in spring table
otu.trans_601_prefire <- otu.trans[row.names(otu.trans)  %in% row.names(predictor2_601_prefire), ]
#check that names match
row.names(otu.trans_601_prefire ) == row.names(predictor2_601_prefire)

#remove samples from otu table not in spring table
otu.trans_601_postfire <- otu.trans[row.names(otu.trans)  %in% row.names(predictor2_601_postfire), ]
#check that names match
row.names(otu.trans_601_postfire) == row.names(predictor2_601_postfire)

##############make rank abundance for the seasons
#rank abundance for bags spring1
RankAbun.601.prefire<- rankabundance(otu.trans_601_prefire)
#plot it
rankabunplot(RankAbun.601.prefire, scale='abundance', addit=FALSE, specnames=c(1,2,3))

#rank abundance for bags spring 2
RankAbun.601.postfire <- rankabundance(otu.trans_601_postfire)
#plot it
rankabunplot(RankAbun.601.postfire, scale='abundance', addit=FALSE, specnames=c(1,2,3))

#make function to get taxonomy on it
rankabundtaxonomyfunction <- function(df, taxa, season) {
  #add OTU info to dataframe
  df <- as.data.frame(df)
  df$id <- row.names(df)
  #make them match
  taxa$id <- as.character(taxa$id)
  match.names <- match(df$id, taxa$id)
  taxanew <-  taxa[match.names, ]
  #check that rows match
  taxanew$id ==row.names(df) 
  #bind together into dataframe and add season info
  cbind(df,taxanew, rep(season,nrow(df)))
}

names(taxon)[2] <- "id"

#apply function to four seasons (still need to figure out apply all)
RankAbun.601.prefire.tax <- rankabundtaxonomyfunction(RankAbun.601.prefire,taxon, "601prefire")

RankAbun.601.postfire.tax <- rankabundtaxonomyfunction(RankAbun.601.postfire,taxon, "601postfire")

#bind them together to one dataframe top 25
RankAbund.601.prepostfire_top30 <- rbind(RankAbun.601.prefire.tax[1:30, ],RankAbun.601.postfire.tax[1:30, ])
write.csv(RankAbund.601.prepostfire_top30, "RankAbund.601.prepostfire_top30.csv")

RankAbund.601.prepostfire_top30 <- read.csv("RankAbund.601.prepostfire_top30.csv")

names(RankAbund.601.prepostfire_top30)[19] <- "Genus"
names(RankAbund.601.prepostfire_top30)[18] <- "Family"
names(RankAbund.601.prepostfire_top30)[9] <-"id2"

names(RankAbund.601.prepostfire_top30)[21] <- "Fire"
#re order seasons 

levels(RankAbund.601.prepostfire_top30$Fire )[levels(RankAbund.601.prepostfire_top30$Fire )=="601prefire"] <- "Pre-fire"
levels(RankAbund.601.prepostfire_top30$Fire )[levels(RankAbund.601.prepostfire_top30$Fire )=="601postfire"] <- "Post-fire"



levels(RankAbund.601.prepostfire_top30$Fire )
RankAbund.601.prepostfire_top30$Fire <- relevel(RankAbund.601.prepostfire_top30$Fire, ref = "Pre-fire", "Post-fire")

write.csv(RankAbund.601.prepostfire_top30, "RankAbund_plot601_top30_2.csv")

RankAbund.601.prepostfire_top30 <- read.csv("RankAbund_plot601_top30_2.csv")
#need to finish fixing spcies2
#make facet wrapped figure by season
rankabundance_601_prepostfire <- ggplot(data=RankAbund.601.prepostfire_top30 , aes(x=reorder(Genus,rank), y=abundance, fill=Phylum)) + geom_bar(stat = "identity")+
  theme_bw()+
  labs(x=" ", y="Plot 601 Number of Samples") + #change y axis label to "Temp C" and remove "Date" for x axis and change legend title
  theme(legend.position="bottom", #put legend on bottom
        strip.text.x = element_text(size = 14, colour = "black"), #make season labels bigger
        axis.text.x=element_text(size=18, angle = 50, hjust = 1),  #change size angle and justification of x axis labels
        legend.text = element_text(size=16),
        legend.title = element_text(size=16),
        axis.text.y=element_text(size=14),  #make y axis tick sizes bigger
        axis.title=element_text(size=16)) #make y axis label larger  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1)+  #add in standard deviation error bars +
rankabundance_601_prepostfire

rankabundance_fig2<- rankabundance_601_prepostfire  + facet_wrap (~Fire, nrow=2, ncol=1)

pdf("RankAbundance_plot601_facetwrap_frequency.pdf",height=8, width=12)
rankabundance_fig2
dev.off()




