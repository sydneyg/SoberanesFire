#May 6, 2019

#this is the OTU table of ITS from the Illumina MiSeq

#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Dropbox/StatsandProgramming/SoberanesFire/")

#source in functions
source('/Users/glassmandylan/Desktop/R_Functions/gettaxondd.R', chdir = TRUE)
source('/Users/glassmandylan/Desktop/R_Functions/getrowsums.R', chdir = TRUE)
source('/Users/glassmandylan/Desktop/R_Functions/pvalueadonis.R', chdir = TRUE)
source('/Users/glassmandylan/Desktop/R_Functions/pvaluelegend.R', chdir = TRUE)
source('/Users/glassmandylan/Desktop/R_Functions/avgdist2.R', chdir = TRUE)

#Packages to install
#install.packages("SPECIES")
#install.packages("vegan")
#install.packages("BiodiversityR")
#install.packages("scales")
#install EcoUtils to get the mean of the rarefied OTU table many times
#install.packages("devtools")
#install_github("GuillemSalazar/EcolUtils")
#install.packages (“ggpubr”)
#install.packages(“nlme”
#install.packages(“multcomp”)

#Packages to load 
library(devtools)
library(tidyverse)
library(EcolUtils)
library(SPECIES)
library(vegan)
library(BiodiversityR)
library(scales)
library(ggpubr)
#load library for nonlinear mixed effect models
library(nlme)
library(multcomp)
library("stringr")
library("plyr")

#read in dataframes
otu_ITS <- read.csv("otu_ITS_with_tax.csv", row.names=1)

####Filter Taxonomy
#make last column into a data frame
lastcolumn <- ncol(otu_ITS)
taxon <- otu_ITS[ ,lastcolumn]
id <- row.names(otu_ITS)
dd <- data.frame(id,taxon)

#divide last row into subsets
otu_taxa1 <- ldply(str_split(string = dd$taxon, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
  names(otu_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otu_taxa2 <- as.data.frame(lapply(otu_taxa1, gsub, pattern=" ", replacement=""))
  dd <- cbind(dd[,1:2 ],otu_taxa2)

  rownames(dd) <- id

dim(dd)
dim(otu_ITS)

#attach to the otu table
otu_ITS_2 <- cbind(otu_ITS, dd)


#look at no DNA controls
noDNA <- which(colnames(otu_ITS) %in%c("ITS1NODNA","ITS1PCRNC"))

#make otu table of no DNA controls
OTU_noDNA <-otu_ITS[ ,noDNA]
colSums(OTU_noDNA)

#remove no DNA controls
otu_ITS_2 <-otu_ITS[ ,-noDNA]
otu_ITS_2

#for example, extarctg 11 from OTU4184 if it's 0 or < 11 zero it, otherwise substract 11
#just some minimal spillover of high abundance OTUs, no big deal

########alpha diversity analysis
#transform table first
dim(otu_ITS_2)
head(otu_ITS_2)
row.names(otu_ITS_2)
otu_ITS_filt_notax_trans <- t(otu_ITS_2[ ,1:71])

#employ getrowsums
getrowsums(otu_ITS_filt_notax_trans)
#now to move on, need ot remove taxonomy files
#something seems off with fungal data - 4 samples in prefire 601 sequenced poorly and no DNA and no PCR controls have tons of reads
#contamination with fungi from dna extractions? from pcr? error in illumina miseq step? error in bioinformatics?

#remove samples with not enough sequences
otu_ITS_filt_notax_trans_2 <- otu_ITS_filt_notax_trans[rowSums(otu_ITS_filt_notax_trans)>12088,]

#rarefy OTU table to remove low abundance sequences but keep rest, normalize to 1,000 seq/sample, normalize 1x
otu_ITS.e12089 <- rrarefy(otu_ITS_filt_notax_trans_2, 12089)


otu_ITS.e12089 <- read.csv("otu_ITS.e12089.csv", row.names=1)
#use EcoUtils function to normalize 100 times and get mean
otu_ITS.e12089_mean <- rrarefy.perm(otu_ITS_filt_notax_trans_2, sample = 12089, n = 100, round.out = T)

#test if this matters - 97% correlated, it doesnt matter
mantel(otu_ITS.e12089_mean, otu_ITS.e12089)
plot(otu_ITS.e12089_mean ~ otu_ITS.e12089)


#remove OTUs that got zero reads after normalization
zeroes <- which(colSums(otu_ITS.e12089)==0)
#renove these columns 
otu_ITS.e12089_nozeroes <- otu_ITS.e12089[ ,-zeroes]
dim(otu_ITS.e12089_nozeroes) #2,889 OTUs

write.csv(otu_ITS.e12089_nozeroes, "otu_ITS.e12089_nozeroes.csv")

otu_ITS.e12089_nozeroes <- read.csv("otu_ITS.e12089_nozeroes.csv", row.names=1)

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(otu_ITS.e12089_nozeroes, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,3500))

#figure out how to plot number of OTUs against number of sequences


#now calculate richness with BioDiversityR package
#run function once to see what it does 
OTU_ITS_richness <- estimateR(otu_ITS.e12089_nozeroes)
#run function on transformed table to switch rows and columns
OTU_ITS_richness<- t(estimateR(otu_ITS.e12089_nozeroes))

#make a function to get estimates and make plots
#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_fungalrichness_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao1", col=alpha("red", 0.5),pch=16)
  #perform correlation test
  cortest1 <- cor.test(estimates2[,2],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest1$estimate, cortest1$p.value)
  mtext("C",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  #perform correlation test
  cortest2 <- cor.test(estimates2[,4],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest2$estimate, cortest2$p.value)
  mtext("D",side=3,adj=0)
  dev.off()
  
}


#run function
estimates_plot_function(otu_ITS.e12089_nozeroes,"otu_e12089")




##############################
#other ways to calculate species richness
##############################

# get species richness fo EMF not rarefied
otu.H <- diversity(otu_ITS.e12089_nozeroes) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(otu_ITS.e12089_nozeroes, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
otu_ITS_e12089_richness <- cbind(otu.richness,OTU_ITS_richness)

write.csv(otu_ITS_e12089_richness, "otu_ITS_e12089_richness.csv")

otu_ITS_e12089_richness <- read.csv("otu_ITS_e12089_richness.csv", row.names=1)

#test if data are normal distirbuted, if significant not normal
shapiro.test(otu_ITS_e12089_richness$S.obs) #not significantly different from normal
histogram(otu_ITS_e12089_richness$S.obs)
qqnorm(otu_ITS_e12089_richness$S.obs)

########################################################################################
########Make some figures
##################################################################
#read in dataframe with METADATA - add in info about cardinal directions and meters
metadata <- read.csv("SoberanesFiremetadata_ITS.csv", row.names=1)

#check that names match, if they dont match, use a matching function
row.names(otu_ITS_e12089_richness) == row.names(metadata)


fungidata <- cbind(otu_ITS_e12089_richness,metadata)

write.csv(fungidata, "fungalrichnessdata.csv")

fungidata$Fire <- factor(fungidata$Fire, levels=c("Prefire","Postfire"))
fungidata$Burn <- factor(fungidata$Burn, levels=c("Unburned","Burned"))


fungibyplot <- ggplot(fungidata, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire, alpha=Burn))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of fungal species (OTUs)") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=1),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
 # stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
  #             label = c("a","a","a","b","a","b"), vjust = -0.5, col="black") +
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                      values=c("cornflowerblue", "tomato1"))  #add in manual colors for points/lines
#figure out how to get tukey letters above
 #stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
  #           label = c("a","a","b"), vjust = -0.5, col="black")

#now you have to figure out how to visualize the fact that Plot58 is not burned etc
fungibyplot


fungibyplot%>%
  ggexport(filename = "fungalrichnesspreandpostfire_MiSeq.pdf")

#what's the % reduction from pre to post fire in the 2 plots?
#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.

fungi_plots <- ddply(fungidata, c("Plot","Fire"), summarise,
                T1_mean = mean(S.obs, na.rm=TRUE),
                T1_sd = sd(S.obs, na.rm=TRUE),
                T1_n = sum(!is.na( S.obs)),
                T1_se = T1_sd/sqrt(T1_n)
)

head(fungi_plots)

#this is based on means but try to figure out with standard error
plot58percentreduction <- (100*(267-286)/267) #7% Increase
plot601percentreduction <- (100*(229-143)/229) #37.55% 
plot603percentreduction <- (100*(321-96)/321) #70.09% 


########################################################################################################
#---------------------------------Do stats on richness: lme ------------------------------------------
########################################################################################################

#test if data are normal distirbuted, if significant not normal, it's not significant, normal enough
shapiro.test(fungidata$S.obs)
histogram(fungidata$S.obs)
#given this distribution, what transformation to use to make it normal so we can use lme? or what family to use for glmer?
#do we use just Plot as a random effect or does pre and post fire also need to be a random effect?
#fungidata$S.obs_log <- log(fungidata$S.obs)
#shapiro.test(fungidata$S.obs_log)
#histogram(fungidata$S.obs_log)

#used linear mixed effect model from packaage nlme to run nested anova
#load library for nonlinear mixed effect models
library(nlme)
names(fungidata)

# species richness response variable
lmeY <- lme(S.obs ~ Fire*Burn, random = ~ 1 | Plot, data = fungidata)
model1 <- anova(lmeY)
model1 
summary(model1)


########################################################################################################
#---------------------------------Do stats on richness: Tukey ------------------------------------------
########################################################################################################

#lme on deadtanoak 2010 #sig effects of dead tanoak and plot by dead tanoak interaction
model1<- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = fungidata)
anova(model1)
tuk_T1 <- glht(model1, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T1 <- cld(tuk_T1) 
tuk.cld.T1
#plot(tuk.cld.T1, ylab="S obs fungi", col=c("cornflowerblue", "tomato1"))
plot(tuk.cld.T1, col=c("cornflowerblue", "tomato1"))


model2<- lme(S.obs ~ Burn, random = ~ 1 | Plot, data = fungidata)
anova(model2)
tuk_T2 <- glht(model2, linfct = mcp(Burn = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
plot(tuk.cld.T2, col=c("cornflowerblue", "tomato1"))


##separate burn vs unburned and test separately

Burnplots <- fungidata[which(fungidata$Burn=="Burned"), ]
unburnplot <- fungidata[which(fungidata$Burn=="Unburned"), ]

#do anova to test effect of fire on fungal richness in burned plots
model1_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = Burnplots)
anova(model1_burn)

model2_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = unburnplot)
anova(model2_burn)

#Tukey test to visualize which are sig different from each other
tuk_T2 <- glht(model1_burn , linfct = mcp(Fire = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
plot(tuk.cld.T2, ylab2="Burned Plots: Observed number of fungal species (OTUs)", xlab2="" ,col=c("cornflowerblue","tomato1")) 

#plot(tuk.cld.T2 ,col=c("cornflowerblue","tomato1")) 

tuk_T3 <- glht(model2_burn, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T3 <- cld(tuk_T3);tuk.cld.T3
plot(tuk.cld.T3, ylab2="Unburned Plots: Observed number of fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

pdf("fungalrichness_prepostfire_Tukey_Miseq.pdf", width=8, height=5)
par(mfrow=c(1,2))
#plot(tuk.cld.T3, ylab="Unburned Plots: Observed no. fungal species (OTUs)", xlab="", col=c("cornflowerblue","tomato1")) 
plot.new()
plot(tuk.cld.T2, ylab2="Burned Plots: Observed no. fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T3, col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T2, col=c("cornflowerblue","tomato1")) 
dev.off()

#set par bac to normal
par(mfrow=c(1,1))

fungibyplot_burns <- ggplot(Burnplots, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of fungal species (OTUs) in Burned Plots") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=14),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=0.05), #Original Size 18, original hjust 0.05
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1"));fungibyplot_burns #add in manual colors for points/lines

fungibyplot_burns 



fungibyplot_unburned <- ggplot(unburnplot, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of fungal OTUs in Unburned plot") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=16),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=0.05),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) #add in manual colors for points/lines


fungibyplot_unburned 

####Atempting to Use Fabi's Code to Make Tukey Plots Better#####

fabi_plot <- ggplot(fungidata, aes(x=Plot, y=S.obs, fill=Burn))+
  geom_boxplot()+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
               label = c("b","a","a"), vjust = -0.5, col="black")+
  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    legend.position = "bottom", text = element_text(size=18),
                    axis.text.y = element_text(size=18, angle=90,hjust=1),
                    axis.text.x = element_text(size=18))+ 
  xlab("Treatment")+ylab("Observed number of fungal species (OTUs)")+
  scale_fill_manual(labels=c("Prefire","Postfire"),
                    values=c("cornflowerblue", "tomato1")); fabi_plot

fabi_plot %>%
  ggexport(filename = "fabi_plot.pdf")
#################### WORKING HERE
attach(metadata)
fabi_plot_burn<-ggplot(Burnplots, aes(x=Plot, y=S.obs, fill=Fire))+
  geom_boxplot()+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
               label = c("b","a","b","a"), vjust = c(1,4,2,1), col="black")+
    theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    legend.position = "bottom", text = element_text(size=18),
                    axis.text.y = element_text(size=18, angle=90,hjust=1),
                    axis.text.x = element_text(size=18))+ 
  xlab("Treatment")+ylab("Observed number of fungal species (OTUs)")+
  scale_fill_manual(labels=c("Prefire","Postfire"),
                    values=c("cornflowerblue", "tomato1"));fabi_plot_burn



fabi_plot_burn %>%
  ggexport(filename = "fabi_plot_burn.pdf")

fabi_plot_unburn <- ggplot(unburnplot, aes(x=Plot, y=S.obs, fill=Burn))+
  geom_boxplot()+
  #  scale_fill_manual(values=c("cornflowerblue", "tomato1","tomato1"))+ 
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
               label = c("b","a","a"), vjust = -0.5, col="black")+
  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    legend.position = "bottom", text = element_text(size=18),
                    axis.text.y = element_text(size=18, angle=90,hjust=1),
                    axis.text.x = element_text(size=18, angle=45, hjust=1))+ 
  xlab("Treatment")+ylab("Observed number of fungal species (OTUs)")+
  scale_fill_manual(labels=c("Prefire","Postfire"),
                    values=c("cornflowerblue", "tomato1"))

fabi_plot %>%
  ggexport(filename = "fabi_plot.pdf")

           
#arrange using ggpubr package  (for example to put fungi and fungal plots side by side
richnessfigures <- ggarrange(fungibyplot_unburned,fungibyplot_burns, labels = c("A", "B"), ncol = 2, nrow = 1)
richnessfigures
richnessfigures %>%
  ggexport(filename = "fungalrichnesspreandpostfire_burnedvunburned_Miseq.pdf")

#########################################################################
#some preliminary beta diversity analysis

#########################################################################



#use vegdist function from vegan and make bray-curtis dissimilarity matrix from rarefied filtered 16S OTU table (1000 seq/sample)
braydist_1TS <- vegdist(otu_ITS.e12089, "bray", upper=TRUE, diag=TRUE)

#use avg dist function rarefy to 1000, make bray-curtis dissimilarity, do this 100 times, get the median, then square root transform
avgdist_e12089_sqrt <- avgdist(otu_ITS_filt_notax_trans, 12089, iterations=100, meanfun=median, transf= sqrt, dmethod="bray" )

mantel(braydist_1TS,avgdist_e12089_sqrt) #89.6% correlated

#check that row names are same as with meta data 
row.names(fungidata)== row.names(as.matrix(braydist_1TS ))



#make the nmds 
#https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation
plotnmds1 <- metaMDS(braydist_1TS, k=2, trymax=100)
plotnmds1 #stress=0.077
stressplot(plotnmds1)

#make dataframe of NMDS scores
scores1 <- as.data.frame(scores(plotnmds1))

#add scores to metadata
betadiversitydata <- cbind(fungidata, scores1)

betadiversitydata <- read.csv("ITS_betadiversity.csv", row.names=1)

#Do adonis to test if there are significant differences of treatments on bray curtis dissimilarity
adonis(braydist_1TS ~ Burn*Fire, data = betadiversitydata, permutations=999)

capture.output(adonis(braydist_1TS ~ Burn*Fire, data = betadiversitydata, permutations=999), file="fungi_adonis.doc")
#Post hoc test to tell which pairs are sig different, install_github("GuillemSalazar/EcolUtils")
adonis.pair(braydist_1TS, betadiversitydata$Fire, nper = 1000, corr.method = "fdr")
adonis.pair(braydist_1TS, betadiversitydata$Burn, nper = 1000, corr.method = "fdr")

#try running adonis on separate braydist16S files

#######

#separate unburned v burned
betadiversitydata

betadiversitydata_Burnplots <- betadiversitydata[which(betadiversitydata$Burn=="Burned"), ]
betadiversitydata_unburnplot <- betadiversitydata[which(betadiversitydata$Burn=="Unburned"), ]



#NMDS for burn plots
NMDS1_burn <- ggplot(betadiversitydata_Burnplots, aes(x=NMDS1, y=NMDS2, col=Fire, shape=Plot, group=Fire)) +
              geom_point(size=2) +
              theme_bw() +
              labs(title = "Fungal Bray-Curtis Dissimilarity")+
              #xlim(-0.2,0.7)+
              #ylim(-0.3,0.5)+
    #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
    #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
            theme(legend.position = "bottom", #put legend under graph,
                  legend.title = element_blank(),
              strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
               axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
              axis.title=element_text(size=18),#change size of x and y axis labels
              legend.spacing = unit(0,"cm"), #change spacing between legends
              legend.text=element_text(size=10), #change size of legend text
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
    stat_ellipse() + #include ellipse around the group
    scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                  values=c("cornflowerblue", "tomato1")) #add in manual colors for 

NMDS1_burn


#NMDS for unburned
NMDS2_unburn <- ggplot(betadiversitydata_unburnplot, aes(x=NMDS1, y=NMDS2, col=Fire, shape=Plot, group=Fire)) +
  geom_point(size=2) +
  theme_bw() +
  labs(title = "Fungal Bray-Curtis Dissimilarity")+
  #xlim(-0.2,0.7)+
  #ylim(-0.3,0.5)+
  #scale_color_manual(values=sitecolors, labels=sitenames) + #add in manual colors and change the legend names
  #scale_shape_manual(values=shapes, labels=sitenames) +  #add in manual shapes and change the legend name
  theme(legend.position = "bottom", #put legend under graph
      legend.title = element_blank(),
      strip.text.x = element_text(size = 18, colour = "black"), #make T1, T2, T3 labels bigger
        axis.text=element_text(size=18,angle=18, hjust=1),  #change size of  x and y axis tick labels
        axis.title=element_text(size=18),#change size of x and y axis labels
        legend.spacing = unit(0,"cm"), #change spacing between legends
        legend.text=element_text(size=10), #change size of legend text
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + #remove major and minor grid panels#+
  stat_ellipse() + #include ellipse around the group
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                   values=c("cornflowerblue", "tomato1")) #add in manual colors

NMDS2_unburn

#arrange using ggpubr package  (for example to put fungi and fungal plots side by side
NMDSfigures <- ggarrange(NMDS2_unburn,NMDS1_burn, labels = c("A", "B"), ncol = 2, nrow = 1)

NMDSfigures

pdf("fungalcommunitycompositionpreandpostfire_burnedvunburned_Miseq.pdf", height=5, width=8)
NMDSfigures
dev.off()


##Attempting to Use Fabi's Script for Richness Boxplots with Tukey Letters#####

#Original Paste in from Fabi's File
install.packages("ggfortify")
install.packages("phyloseq")

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(phyloseq)
library(tidyr)#manipulate data
library(dplyr)
library(multcomp)

#A<-ggplot(data, aes(x=TimeSinceFire, y=pH, group=Treatment, col=Treatment,shape=Treatment))+
#  stat_summary(fun.y=mean,geom="point", size=3)+
#  stat_summary(fun.data = mean_se,geom = "errorbar", size=0.5,
#               alpha=0.7,position = position_dodge(0.01))+
#  scale_colour_manual(values = c("#a50f15", "#ef3b2c","#f4a582","#0571b0"))+
#  geom_smooth(method = lm, formula=y ~ poly(x, 3), se=F, alpha=0.05)+
#  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                    axis.text.x = element_text(angle = 45, hjust=1, size=13), 
#                    axis.title=element_text(size=13,face="bold"))+
#  xlab("Time Since Fire")+ ylab("pH");A


#ggplot(T1, aes(x=Treatment, y=Ash, fill=Treatment))+
#  geom_boxplot()+
#  scale_fill_manual(values=c("#a50f15", "#ef3b2c","#f4a582"))+ 
#  stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
#               label = c("b","a","a"), vjust = -0.5, col="black")+
#  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
#                    axis.text.x = element_text(angle = 45,hjust=1, size=13), 
#                    axis.title=element_text(size=13,face="bold"))+
#  xlab("Treatment")+ylab("Ash Depth (cm)")
#dev.off() 


fungibyplot_burns <- ggplot(Burnplots, aes(x=Plot, y=S.obs, fill=Fire))+
  geom_boxplot()+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
               label = c("a","a","a","b","a","b"), col="black")+
  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    axis.text.x = element_text(angle = 45,hjust=1, size=13), 
                    axis.title=element_text(size=13,face="bold"))+
  ylab("Observed number of fungal species (OTUs) in Burned Plots") +
  theme(legend.position = "bottom", #put legend under graph
  text = element_text(size=18),
axis.title.x=element_blank(), #remove Plot title
axis.text.y = element_text(size=18, angle=90,hjust=0.05),
axis.text.x = element_text(size=18, angle=45, hjust=1)) +
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1"))

fungibyplot_burns

par(mfrow=c(1,1))

fungibyplot_burns <- ggplot(Burnplots, aes(x=Plot, y=S.obs)) + 
  geom_boxplot(aes(fill=Fire))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of fungal species (OTUs) in Burned Plots") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=0.05),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  stat_summary(fun.y=mean,geom="point", size=3)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=0.5,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(.7), 
             label = c("a","a","a","b","a","b"), vjust = -0.5, col="black")

fungibyplot_burns 



fungibyplot_unburned <- ggplot(unburnplot, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of fungal OTUs in Unburned plot") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=20),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=0.05),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) #add in manual colors for points/lines


fungibyplot_unburned 


#arrange using ggpubr package  (for example to put fungi and fungal plots side by side
richnessfigures <- ggarrange(fungibyplot_unburned,fungibyplot_burns, labels = c("A", "B"), ncol = 2, nrow = 1)
richnessfigures
richnessfigures %>%
  ggexport(filename = "fungalrichnesspreandpostfire_burnedvunburned_Miseq.pdf")
