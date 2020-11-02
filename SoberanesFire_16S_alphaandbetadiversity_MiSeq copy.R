#April 30, 2019

#this is the OTU table of 16S from the Illumina MiSeq

#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Dropbox/StatsandProgramming/SoberanesFire/")

#source in functions
source('~/Dropbox/StatsandProgramming/source/gettaxondd.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/getrowsums.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/pvalueadonis.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/pvaluelegend.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/avgdist2.R', chdir = TRUE)
source('~/Dropbox/StatsandProgramming/source/getbraydistfunction.R', chdir = TRUE)

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
otu_16S <- read.csv("16Sotu_table_filtered_2_MiSeq.csv", row.names=1)

####Filter Taxonomy
#make last column into a data frame
lastcolumn <- ncol(otu_16S)
taxon <- otu_16S[ ,lastcolumn]
id <- row.names(otu_16S)
dd <- data.frame(id,taxon)

#divide last row into subsets
otu_taxa1 <- ldply(str_split(string = dd$taxon, pattern=";"), rbind) # Divide a column using ";"and convert list to data frame
  names(otu_taxa1) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  otu_taxa2 <- as.data.frame(lapply(otu_taxa1, gsub, pattern=" ", replacement=""))
  dd <- cbind(dd[,1:2 ],otu_taxa2)

  rownames(dd) <- id

dim(dd)
dim(otu_16S)

#attach to the otu table
otu_16S_2 <- cbind(otu_16S, dd)

write.csv(otu_16S_2, "16S_otu_with_tax.csv")
#remove everything that is not bacteria, or "unassigned", remove chloroplasts and mitochondria
#unassigned <- which(otu_16S_2$Consensus.Lineage %in%c("Unassigned"))
#chloroplast <- which(otu_16S_2$Class %in%c("c__Chloroplast"))
#mitochondria <- which(otu_16S_2$Family %in%c("f__mitochondria"))

#combine into a vector  
#tofilter <- c(unassigned,chloroplast,mitochondria)

#remove those row from OTU table
#otu_16S_filtered <- otu_16S_2 [- tofilter, ]

#check dimensions of otu table before and after filtering
#dim(otu_16S)
#dim(otu_16S_filtered)

#write as csv 
#write.csv(otu_16S_filtered, "otu_16S_filtered.csv")

#remove taxonomy columns
#otu_16S_filt_notax  <- otu_16S_filtered[ ,1:74]

#write as csv 
#write.csv(otu_16S_filt_notax, "otu_16S_filtered_notax.csv")

#look at no DNA controls
noDNA <- which(colnames(otu_16S) %in%c("X16SPCRNC","X16SNODNA"))

#make otu table of no DNA controls
OTU_noDNA <-otu_16S[ ,noDNA]
colSums(OTU_noDNA)

#for example, extarctg 11 from OTU4184 if it's 0 or < 11 zero it, otherwise substract 11
#just some minimal spillover of high abundance OTUs, no big deal

########alpha diversity analysis
#transform table first
dim(otu_16S)
otu_16S_filt_notax_trans <- t(otu_16S[ ,1:74])

#employ getrowsums
getrowsums(otu_16S_filt_notax_trans)
#now to move on, need ot remove taxonomy files


#remove samples with not enough sequences
otu_16S_filt_notax_trans_2 <- otu_16S_filt_notax_trans[rowSums(otu_16S_filt_notax_trans)>10366,]

#rarefy OTU table to remove low abundance sequences but keep rest, normalize to 1,000 seq/sample, normalize 1x
otu_16S.e10367 <- rrarefy(otu_16S_filt_notax_trans_2, 10367)

otu_16S.e10367 <- read.csv("otu_16S.e10367.csv", row.names =1, check.names = FALSE)
#use EcoUtils function to normalize 100 times and get mean
otu_16S.e10367_mean <- rrarefy.perm(otu_16S_filt_notax_trans_2, sample =10367, n = 100, round.out = T)

#test if this matters - 97% correlated, it doesnt matter
mantel(otu_16S.e10367_mean, otu_16S.e10367)
#plot(otu_16S.e10367_mean ~ otu_16S.e10367)


#remove OTUs that got zero reads after normalization
zeroes <- which(colSums(otu_16S.e10367)==0)
#renove these columns 
otu_16S.e1000_nozeroes <- otu_16S.e10367[ ,-zeroes]
dim(otu_16S.e1000_nozeroes)

otu_16S.e1000_nozeroes <- read.csv("otu_16S.e10367_nozeroes.csv", row.names=1, check.names=FALSE)

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(otu_16S.e1000_nozeroes, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,13000))

#figure out how to plot number of OTUs against number of sequences


#now calculate richness with BioDiversityR package
#run function once to see what it does 
OTU_16_richness <- estimateR(otu_16S.e1000_nozeroes)
#run function on transformed table to switch rows and columns
OTU_16_richness<- t(estimateR(otu_16S.e1000_nozeroes))

#make a function to get estimates and make plots
#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_Bacterialrichness_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao1", col=alpha("red", 0.5),pch=16)
  #perform correlation test
  cortest1 <- cor.test(estimates2[,2],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest1$estimate, cortest1$p.value)
  mtext("A",side=3,adj=0)
  #plot S. Ace vs S.obs to see if they are correlated
  plot(estimates2[,4] ~estimates2[,1], xlab="S.obs",ylab="ACE",col=alpha("black", 0.5),pch=16)
  #perform correlation test
  cortest2 <- cor.test(estimates2[,4],estimates2[,1] )
  #invoke pvaluelegend to put pvalue on figure
  pvaluelegend(cortest2$estimate, cortest2$p.value)
  mtext("B",side=3,adj=0)
  dev.off()
  
}


#run function
estimates_plot_function(otu_16S.e1000_nozeroes,"otu_e1000")




##############################
#other ways to calculate species richness
##############################

# get species richness fo EMF not rarefied
otu.H <- diversity(otu_16S.e1000_nozeroes) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(otu_16S.e1000_nozeroes, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
OTU_16S_e1000_richness <- cbind(otu.richness,OTU_16_richness)

write.csv(OTU_16S_e1000_richness, "OTU_16S_e10367_richness.csv")

OTU_16S_e1000_richness <- read.csv("OTU_16S_e10367_richness.csv", row.names=1)

#test if data are normal distirbuted, if significant not normal
shapiro.test(OTU_16S_e1000_richness$S.obs)
histogram(OTU_16S_e1000_richness$S.obs)

########################################################################################
########Make some figures
##################################################################
#read in dataframe with METADATA - add in info about cardinal directions and meters
metadata <- read.csv("SoberanesFiremetadata_16S.csv", row.names=1)

#check that names match, if they dont match, use a matching function
row.names(OTU_16S_e1000_richness) == row.names(metadata)


bacteriadata <- cbind(OTU_16S_e1000_richness,metadata)

write.csv(bacteriadata, "bacteriarichnessdata.csv")

bacteriadata$Fire <- factor(bacteriadata$Fire, levels=c("Prefire","Postfire"))
bacteriadata$Burn <- factor(bacteriadata$Burn, levels=c("Unburned","Burned"))


bacteriabyplot <- ggplot(bacteriadata, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire, alpha=Burn))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of bacterial species (OTUs)") +  #change yaxis label
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
bacteriabyplot


bacteriabyplot%>%
  ggexport(filename = "Bacterialrichnesspreandpostfire_MiSeq.pdf")


######Adding Tukey Letters to the Previous Figure Using Fabi's Codes#########

attach(metadata)
richness_w_tukey<-ggplot(bacteriadata, aes(x=Plot, y=S.obs, fill=Fire, alpha=Burn))+
  geom_boxplot()+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(1), 
               label = c("a","a","b","a","b","a"), vjust = c(1,1,1,1,17,0.1), col="black")+
  theme_bw()+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                    legend.position = "bottom", text = element_text(size=17),
                    axis.text.y = element_text(size=18, angle=90,hjust=1),
                    axis.text.x = element_text(size=18))+ 
  ylab("Observed number of bacterial species (OTUs)")+
  scale_fill_manual(labels=c("Prefire","Postfire"),
                    values=c("cornflowerblue", "tomato1"));fabi_plot_burn

richness_w_tukey%>%
  ggexport(filename = "Bacterialrichnesspreandpostfire_MiSeq_Tukey.pdf")

#what's the % reduction from pre to post fire in the 2 plots?
#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.

bac_plots <- ddply(bacteriadata, c("Plot","Fire"), summarise,
                T1_mean = mean(S.obs, na.rm=TRUE),
                T1_sd = sd(S.obs, na.rm=TRUE),
                T1_n = sum(!is.na( S.obs)),
                T1_se = T1_sd/sqrt(T1_n)
)

head(bac_plots)

#this is based on means but try to figure out with standard error
plot58percentreduction <- (100*(2183.0833-2187.6667)/2183.0833) #0.2099508% growth
plot601percentreduction <- (100*(2082.25-1254.5833)/2082.25) #39.74867% reduction
plot603percentreduction <- (100*(2188.8333-843.8333)/2188.8333) #61.44826 % percent reduction


plot58percentreduction
plot601percentreduction
plot603percentreduction

########Just do a quick and dirty ANOVA to test for signfiicant differences, doesn't account for data not being normal, prob want to do a GLM

#load library for nonlinear mixed effect models
library(nlme)

# species richness response variable
lmeY <- lme(S.obs ~ Fire*Burn, random = ~ 1 | Plot, data = bacteriadata)
model1 <- anova(lmeY)
model1 
summary(model1)

Fit1 <- (aov(S.obs ~ Fire*Burn, data = bacteriadata))
summary(Fit1 )

#Perform post-hoc test to determine which comparisons are sig different
TukeyHSD(Fit1, ordered = FALSE, conf.level = 0.95)
#Another method of visualizing Tukey differences
######T1
model_1 <-aov(S.obs ~ Fire,  data=bacteriadata)
tuk_T1 <- glht(model_1, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T1 <- cld(tuk_T1) 
tuk.cld.T1
plot(tuk.cld.T1, ylab="S obs bacteria",xlab="" ) 

##separate burn vs unburned and test separately

Burnplots <- bacteriadata[which(bacteriadata$Burn=="Burned"), ]
unburnplot <- bacteriadata[which(bacteriadata$Burn=="Unburned"), ]

#do anova to test effect of fire on bacterial richness in burned plots
Fit2 <- (aov(S.obs ~ Fire, data = Burnplots))
summary(Fit2 )

#Tukey test to visualize which are sig different from each other
tuk_T2 <- glht(Fit2, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
plot(tuk.cld.T2, ylab="Burned Plots: Observed number of bacterial species (OTUs)", xlab="" ,col=c("cornflowerblue","tomato1")) 

Fit3 <- (aov(S.obs ~ Fire, data = unburnplot))
summary(Fit3 )

tuk_T3 <- glht(Fit3, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T3 <- cld(tuk_T3) 
tuk.cld.T3
plot(tuk.cld.T3, ylab="Unburned Plots: Observed number of bacterial species (OTUs)",xlab="", col=c("cornflowerblue","tomato1")) 

pdf("Bacterialrichness_prepostfire_Tukey_Miseq.pdf", width=8, height=5)
par(mfrow=c(1,2))
plot(tuk.cld.T3, ylab="Unburned Plots: Observed no. bacterial species (OTUs)", xlab="", col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T2, ylab="Burned Plots: Observed no. bacterial species (OTUs)", xlab="", col=c("cornflowerblue","tomato1")) 
dev.off()

#set par bac to normal
par(mfrow=c(1,1))

bacteriabyplot_burns <- ggplot(Burnplots, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of bacterial species (OTUs) in Burned Plots") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=0.05),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) #add in manual colors for points/lines


bacteriabyplot_burns 



bacteriabyplot_unburned <- ggplot(unburnplot, aes(x=Plot, y=S.obs)) +   
  geom_boxplot(aes(fill=Fire))+ #fill colors by burned v unburned
  theme_bw() + #make black and white
  ylab("Observed number of bacterial OTUs in Unburned plot") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=20),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=18, angle=90,hjust=0.05),
        axis.text.x = element_text(size=18, angle=45, hjust=1)) +  #make x axis text larger and angled
  scale_fill_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                    values=c("cornflowerblue", "tomato1")) #add in manual colors for points/lines


bacteriabyplot_unburned 

           
#arrange using ggpubr package  (for example to put bacteria and fungal plots side by side
richnessfigures <- ggarrange(bacteriabyplot_unburned,bacteriabyplot_burns, labels = c("A", "B"), ncol = 2, nrow = 1)
richnessfigures
richnessfigures %>%
  ggexport(filename = "Bacterialrichnesspreandpostfire_burnedvunburned_Miseq.pdf")

#########################################################################
#some preliminary beta diversity analysis

#########################################################################



#use vegdist function from vegan and make bray-curtis dissimilarity matrix from rarefied filtered 16S OTU table (1000 seq/sample)
braydist_16S <- vegdist(otu_16S.e10367, "bray", upper=TRUE, diag=TRUE)

#use avg dist function rarefy to 1000, make bray-curtis dissimilarity, do this 100 times, get the median, then square root transform
avgdist_e1000_sqrt <- avgdist(otu_16S_filt_notax_trans, 10367, iterations=100, meanfun=median, transf= sqrt, dmethod="bray" )

mantel(braydist_16S,avgdist_e1000_sqrt) #98% correlated

#check that row names are same as with meta data 
row.names(bacteriadata)== row.names(as.matrix(braydist_16S ))



#make the nmds 
#https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
#A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions, < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation
plotnmds1 <- metaMDS(braydist_16S, k=3, trymax=100)
plotnmds1 #stress=0.077
stressplot(plotnmds1)

#make dataframe of NMDS scores
scores <- as.data.frame(scores(plotnmds1))

#add scores to metadata
betadiversitydata <- cbind(bacteriadata, scores)

betadiversitydata <- read.csv("16S_betadiversity.csv", row.names=1, check.names=FALSE)

#Do adonis to test if there are significant differences of treatments on bray curtis dissimilarity
adonis(braydist_16S ~ Burn*Fire, data = betadiversitydata, permutations=999)

capture.output(adonis(braydist_16S ~ Burn*Fire, data = betadiversitydata, permutations=999), file="bacteria_adonis.doc")

#Post hoc test to tell which pairs are sig different, install_github("GuillemSalazar/EcolUtils")
adonis.pair(braydist_16S, betadiversitydata$Fire, nper = 1000, corr.method = "fdr")
adonis.pair(braydist_16S, betadiversitydata$Burn, nper = 1000, corr.method = "fdr")

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
              labs(title = "Bacterial Bray-Curtis Dissimilarity")+
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
  labs(title = "Bacterial Bray-Curtis Dissimilarity")+
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

#arrange using ggpubr package  (for example to put bacteria and fungal plots side by side
NMDSfigures <- ggarrange(NMDS2_unburn,NMDS1_burn, labels = c("A", "B"), ncol = 2, nrow = 1)

NMDSfigures

pdf("Bacterialcommunitycompositionpreandpostfire_burnedvunburned_Miseq.pdf", height=5, width=8)
NMDSfigures
dev.off()


######Making NMDS of total samples, NOT seperated by Burn#####




