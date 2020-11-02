#April 22 2020
#Soberanes Fire FunGuild Data

#Reset R's Brain
rm(list=ls())

#set working directory
setwd("~/Desktop/Big_Sur_Project/Big_Sur_Full_Run/FunGuild")

getwd()

#Make Sure R is up to date

#install.packages("installr")
#library(installr)

#updateR()

#Packages to load 
#install.packages("devtools")
library(devtools)
#install.packages("tidyverse")
library(tidyverse)
#install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)
#install.packages("SPECIES")
library(SPECIES)
#install.packages("vegan")
library(vegan)
#install.packages("BiodiversityR")
#install.packages("vegan3d")
#install.packages("shiny")
library(BiodiversityR)
library(scales)
#install.packages("ggpubr")
library(ggpubr)
#load library for nonlinear mixed effect models
library(nlme)
library(multcomp)
library("stringr")
library("plyr")


#read in matched OTU data

funguild <- read.delim("funguilds_fullset.txt")
dim(funguild)
head(funguild)

class(funguild)
class(funguild$Confidence)
str(funguild$Confidence["Highly_Probable"])

#Pull out only Highly Probable

fun_high <- subset(funguild, Confidence=="Highly_Probable")

dim(fun_high)

head(fun_high)

#save this so I have the whole highly probable set in case I need to come back to it.

write.csv(fun_high, "funguild_highly_probable.csv")

fun_high <- read.csv("funguild_highly_probable.csv", row.names =1)

head(fun_high)

#Now subset to only keep Ectomycorrhizal and Arbuscular Mycorrhizal taxa

fun_mf <- subset(fun_high, Guild=="Ectomycorrhizal" | Guild=="Arbuscular_Mycorrhizal")

#Write this to a CSV to be able to examine change in Excel

write.csv(fun_mf, "amf_emf.csv")

#Need to look at richness change in Ectomycorrhizae pre-post fire only

fun_emf <- subset(fun_high, Guild=="Ectomycorrhizal")
dim(fun_emf)
head(fun_emf)


#Going to use the already rarified otu table from the FULL dataset
#then going to join it the emf data using an innerjoin function

#First I am going to take the emf subset and save it

write.csv(fun_emf, "fun_emf.csv")

dim(fun_emf)
row.names(fun_emf)

#now read in the full rarified OTU table

rarefied_ITS <- read.csv("otu_ITS.e12089_nozeroes.csv", row.names =1)

#transpose it so that row names are OTU IDs

trans_12089 = t(rarefied_ITS)

dim(trans_12089)

row.names(trans_12089)

write.csv(trans_12089, "trans_12089.csv")

#Changed both fun_emf and trans_12089 in Excel so that OTU column has header "OTUID"
#Now I can read these in and use that to join them

trans_12089_2 <- read.csv("trans_12089_2.csv")
row.names(trans_12089_2)
trans_12089_2$OTUID

fun_emf_2 <- read.csv("fun_emf_2.csv")
row.names(fun_emf_2)
fun_emf_2$OTUID

#Now use the semi_join function to keep just the OTUs from the emf set

emf_12089 <- semi_join(trans_12089_2, fun_emf_2, by="OTUID", copy=FALSE)
dim(emf_12089)
row.names(emf_12089)
emf_12089$OTUID

#We lose 26 emf taxa through rarefaction. Now save this and strip the row name column in Excel

write.csv(emf_12089, "emf_12089.csv")

emf_12089<- read.csv("emf_12089.csv", row.names =1)

dim(emf_12089)

row.names(emf_12089)

#Nice! that worked. Now to go back through and do all the richness again from this point.

#First transpose the table again

emf_12089_t <- t(emf_12089)

#Don't need to getrowsums or to rarefy since this table is already rarefied. Just Dive in.

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(emf_12089_t, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,300))

#now calculate richness with BioDiversityR package
#run function once to see what it does 
emf_richness <- estimateR(emf_12089_t)
#run function on transformed table to switch rows and columns
emf_richness<- t(estimateR(emf_12089_t))

#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_emf_richness_fullset_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
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
estimates_plot_function(emf_12089_t,"emf_fullset")

##############################
#other ways to calculate species richness
##############################

# get species richness fo emf not rarefied
otu.H <- diversity(emf_12089_t) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(emf_12089_t, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
emf_12089_richness <- cbind(otu.richness,emf_richness)

write.csv(emf_12089_richness, "emf_fullset_richness.csv")

emf_12089_richness <- read.csv("emf_fullset_richness.csv", row.names=1)

#test if data are normal distirbuted, if significant not normal
shapiro.test(emf_12089_richness$S.obs) # less significantly different but still p < 0.5 (p=0.01472)
histogram(emf_12089_richness$S.obs) #skews left of center
qqnorm(emf_12089_richness$S.obs)

########################################################################################
########Make some figures
##################################################################

dim(emf_12089_richness)

#read in dataframe with METADATA - add in info about cardinal directions and meters
metadata <- read.csv("SoberanesFiremetadata_ITS.csv", row.names =1)

dim(metadata)

#check that names match, if they dont match, use a matching function
row.names(emf_12089_richness) == row.names(metadata)

#bind it all together now
emfdata <- cbind(emf_12089_richness,metadata)

write.csv(emfdata, "emfdata_fullset.csv")
emfdata <- read.csv("emfdata_fullset.csv")

emfdata$Fire <- factor(emfdata$Fire, levels=c("Prefire","Postfire"))
emfdata$Burn <- factor(emfdata$Burn, levels=c("Unburned","Burned"))

#make a figure
emfbyplot <- ggplot(emfdata, aes(x=Plot, y=S.obs)) +   
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
emfbyplot


emfbyplot%>%
  ggexport(filename = "fullset_emfrichnesspreandpostfire_MiSeq.pdf")

#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.

emf_plots <- ddply(emfdata, c("Plot","Fire"), summarise,
                   T1_mean = mean(S.obs, na.rm=TRUE),
                   T1_sd = sd(S.obs, na.rm=TRUE),
                   T1_n = sum(!is.na( S.obs)),
                   T1_se = T1_sd/sqrt(T1_n)
)

head(emf_plots)

#this is based on means but try to figure out with standard error
plot58percentreduction <- (100*(13.42-18.08)/18.08) #25.77% Increase
plot601percentreduction <- (100*(8.92-7.08)/8.92) #20.63% Decrease
plot603percentreduction <- (100*(13.27-4.3)/13.27) #67.60% Decrease

#Lets do the stats in the richness

#used linear mixed effect model from packaage nlme to run nested anova
#load library for nonlinear mixed effect models
library(nlme)
names(emfdata)

# species richness response variable
lmeY <- lme(S.obs ~ Fire*Burn, random = ~ 1 | Plot, data = emfdata)
model1 <- anova(lmeY)
model1 
summary(model1)
########################################################################################################
#---------------------------------Do stats on richness: Tukey ------------------------------------------
########################################################################################################

#lme on deadtanoak 2010 #sig effects of dead tanoak and plot by dead tanoak interaction
model1<- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = emfdata)
anova(model1)
tuk_T1 <- glht(model1, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T1 <- cld(tuk_T1) 
tuk.cld.T1
#plot(tuk.cld.T1, ylab="S obs fungi", col=c("cornflowerblue", "tomato1"))
plot(tuk.cld.T1, col=c("cornflowerblue", "tomato1"))


model2<- lme(S.obs ~ Burn, random = ~ 1 | Plot, data = emfdata)
anova(model2)
tuk_T2 <- glht(model2, linfct = mcp(Burn = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
plot(tuk.cld.T2, col=c("cornflowerblue", "tomato1"))


##separate burn vs unburned and test separately

Burnplots <- emfdata[which(emfdata$Burn=="Burned"), ]
unburnplot <- emfdata[which(emfdata$Burn=="Unburned"), ]

#do anova to test effect of fire on fungal richness in burned plots
model1_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = Burnplots)
burn_anova <- anova(model1_burn)
summary(burn_anova)

model2_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = unburnplot)
anova(model2_burn)

#Tukey test to visualize which are sig different from each other
tuk_T2 <- glht(model1_burn , linfct = mcp(Fire = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
a<- plot(tuk.cld.T2, ylab2="Burned Plots: Observed number of fungal species (OTUs)", xlab2="" ,col=c("cornflowerblue","tomato1")) 

a + facet_wrap(~Plot)
#plot(tuk.cld.T2 ,col=c("cornflowerblue","tomato1")) 

tuk_T3 <- glht(model2_burn, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T3 <- cld(tuk_T3);tuk.cld.T3
plot(tuk.cld.T3, ylab2="Unburned Plots: Observed number of fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

pdf("emfrichness_prepostfire_Tukey_Miseq.pdf", width=8, height=5)
par(mfrow=c(1,2))
#plot(tuk.cld.T3, ylab="Unburned Plots: Observed no. fungal species (OTUs)", xlab="", col=c("cornflowerblue","tomato1")) 
plot.new()
plot(tuk.cld.T2, ylab2="Burned Plots: Observed no. fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T3, col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T2, col=c("cornflowerblue","tomato1")) 
dev.off()

#set par bac to normal
par(mfrow=c(1,1))

#Now let's make the mean and SE figure like the other panels

C<-ggplot(emfdata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
 stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,-1,0),
               vjust=c(7,3,6,6.8,11.3,4),
               label = c("a","a","a","b","b","b"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Ectomycorrhizal Fungi") +
  theme_bw() + #make black and white
  ylab("Mean Number of Soil Fungal OTUs per Sample") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=20, angle=90,hjust=1),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled
#C

C %>%
  ggexport(filename = "emf_richnesspreandpostfire_MiSeq_meanandSE_one_tukey.pdf")

####Time to examine the change in emfrobes


#Now subset to only keep saprotroph taxa

fun_high <- read.csv("funguild_highly_probable_sap.csv", row.names =1)

fun_sap <- subset(fun_high, Guild==%in%"saprotroph"%in%)

write.csv(fun_sap, "funguild_saprotroph.csv")

#Need to go in Excel and make the first column say "OTUID"

#Now read in the transposed full OTU table and also read in the fixed saprotroph table

trans_12089_2 <- read.csv("trans_12089_2.csv")
row.names(trans_12089_2)
trans_12089_2$OTUID

fun_sap_2 <- read.csv("fun_sap_2.csv")
row.names(fun_sap_2)
fun_sap_2$OTUID

#Now use the semi_join function to keep just the OTUs from the sap set

sap_12089 <- semi_join(trans_12089_2, fun_sap_2, by="OTUID", copy=FALSE)
dim(sap_12089)
#We lose 10 OTUs to rarefaction
row.names(sap_12089)
sap_12089$OTUID

#We lose 10 sap taxa through rarefaction. Now save this and strip the row name column in Excel

write.csv(sap_12089, "sap_12089.csv")

sap_12089<- read.csv("sap_12089.csv", row.names =1)

dim(sap_12089)

row.names(sap_12089)

#First transpose the table again

sap_12089_t <- t(sap_12089)

#Don't need to getrowsums or to rarefy since this table is already rarefied. Just Dive in.

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(sap_12089_t, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,300))

#now calculate richness with BioDiversityR package
#run function once to see what it does 
sap_richness <- estimateR(sap_12089_t)
#run function on transformed table to switch rows and columns
sap_richness<- t(estimateR(sap_12089_t))

#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_sap_richness_fullset_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
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
estimates_plot_function(sap_12089_t,"sap_fullset")

##############################
#other ways to calculate species richness
##############################

# get species richness fo sap not rarefied
otu.H <- diversity(sap_12089_t) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(sap_12089_t, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
sap_12089_richness <- cbind(otu.richness,sap_richness)

write.csv(sap_12089_richness, "sap_fullset_richness.csv")

sap_12089_richness <- read.csv("sap_fullset_richness.csv", row.names=1)

#test if data are normal distirbuted, if significant not normal
shapiro.test(sap_12089_richness$S.obs) # less significantly different but still p < 0.5 (p=0.001296)
histogram(sap_12089_richness$S.obs) #skews left of center
qqnorm(sap_12089_richness$S.obs)

########################################################################################
########Make some figures
##################################################################

dim(sap_12089_richness)

#read in dataframe with METADATA - add in info about cardinal directions and meters
metadata <- read.csv("SoberanesFiremetadata_ITS.csv", row.names =1)

dim(metadata)

#check that names match, if they dont match, use a matching function
row.names(sap_12089_richness) == row.names(metadata)

#bind it all together now
sapdata <- cbind(sap_12089_richness,metadata)

write.csv(sapdata, "sapdata_fullset.csv")
sapdata <- read.csv("sapdata_fullset.csv")

sapdata$Fire <- factor(sapdata$Fire, levels=c("Prefire","Postfire"))
sapdata$Burn <- factor(sapdata$Burn, levels=c("Unburned","Burned"))

#make a figure
sapbyplot <- ggplot(sapdata, aes(x=Plot, y=S.obs)) +   
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
sapbyplot


sapbyplot%>%
  ggexport(filename = "fullset_saprichnesspreandpostfire_MiSeq.pdf")

#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.

sap_plots <- ddply(sapdata, c("Plot","Fire"), summarise,
                   T1_mean = mean(S.obs, na.rm=TRUE),
                   T1_sd = sd(S.obs, na.rm=TRUE),
                   T1_n = sum(!is.na( S.obs)),
                   T1_se = T1_sd/sqrt(T1_n)
)

head(sap_plots)

#this is based on means but try to figure out with standard error
plot58percentreduction <- (100*(4.33-5.6667)/5.6667) #23.59% Increase
plot601percentreduction <- (100*(4.5-1.4167)/4.5) #68.52% Decrease
plot603percentreduction <- (100*(4.636-1.1)/4.636) #76.27% Decrease

#Lets do the stats in the richness

#used linear mixed effect model from packaage nlme to run nested anova
#load library for nonlinear mixed effect models
library(nlme)
names(sapdata)

# species richness response variable
lmeY <- lme(S.obs ~ Fire*Burn, random = ~ 1 | Plot, data = sapdata)
model1 <- anova(lmeY)
model1 
summary(model1)
########################################################################################################
#---------------------------------Do stats on richness: Tukey ------------------------------------------
########################################################################################################

#lme on deadtanoak 2010 #sig effects of dead tanoak and plot by dead tanoak interaction
model1<- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = sapdata)
anova(model1)
tuk_T1 <- glht(model1, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T1 <- cld(tuk_T1) 
tuk.cld.T1
#plot(tuk.cld.T1, ylab="S obs fungi", col=c("cornflowerblue", "tomato1"))
plot(tuk.cld.T1, col=c("cornflowerblue", "tomato1"))


model2<- lme(S.obs ~ Burn, random = ~ 1 | Plot, data = sapdata)
anova(model2)
tuk_T2 <- glht(model2, linfct = mcp(Burn = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
plot(tuk.cld.T2, col=c("cornflowerblue", "tomato1"))


##separate burn vs unburned and test separately

Burnplots <- sapdata[which(sapdata$Burn=="Burned"), ]
unburnplot <- sapdata[which(sapdata$Burn=="Unburned"), ]

#do anova to test effect of fire on fungal richness in burned plots
model1_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = Burnplots)
burn_anova <- anova(model1_burn)
summary(burn_anova)

model2_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = unburnplot)
anova(model2_burn)

#Tukey test to visualize which are sig different from each other
tuk_T2 <- glht(model1_burn , linfct = mcp(Fire = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
a<- plot(tuk.cld.T2, ylab2="Burned Plots: Observed number of fungal species (OTUs)", xlab2="" ,col=c("cornflowerblue","tomato1")) 

a + facet_wrap(~Plot)
#plot(tuk.cld.T2 ,col=c("cornflowerblue","tomato1")) 

tuk_T3 <- glht(model2_burn, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T3 <- cld(tuk_T3);tuk.cld.T3
plot(tuk.cld.T3, ylab2="Unburned Plots: Observed number of fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

pdf("saprichness_prepostfire_Tukey_Miseq.pdf", width=8, height=5)
par(mfrow=c(1,2))
#plot(tuk.cld.T3, ylab="Unburned Plots: Observed no. fungal species (OTUs)", xlab="", col=c("cornflowerblue","tomato1")) 
plot.new()
plot(tuk.cld.T2, ylab2="Burned Plots: Observed no. fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T3, col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T2, col=c("cornflowerblue","tomato1")) 
dev.off()

#set par bac to normal
par(mfrow=c(1,1))

#Now let's make the mean and SE figure like the other panels

D<-ggplot(sapdata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,-1,0),
               vjust=c(2.5,11,10.5,6.8,11.3,2.5),
               label = c("a","a","a","b","b","b"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Saprotrophic Fungi") +
  theme_bw() + #make black and white
  ylab("Mean Number of Soil Fungal OTUs per Sample") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=20, angle=90,hjust=1),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled


D %>%
  ggexport(filename = "sap_richnesspreandpostfire_MiSeq_meanandSE_one_tukey.pdf")


####Time to examine the change in AMF


#Now subset to only keep AMF taxa

fun_high <- read.csv("funguild_highly_probable_sap.csv", row.names =1)

fun_amf <- subset(fun_high, Guild=="Arbuscular_Mycorrhizal")

write.csv(fun_amf, "funguild_amf.csv")

#Need to go in Excel and make the first column say "OTUID"

#Now read in the transposed full OTU table and also read in the fixed amf table

trans_12089_2 <- read.csv("trans_12089_2.csv")
row.names(trans_12089_2)
trans_12089_2$OTUID

fun_amf_2 <- read.csv("funguild_amf_2.csv")
row.names(fun_amf_2)
fun_amf_2$OTUID

#Now use the semi_join function to keep just the OTUs from the amf set

amf_12089 <- semi_join(trans_12089_2, fun_amf_2, by="OTUID", copy=FALSE)
dim(amf_12089)
#We lose no OTUs to rarefaction
row.names(amf_12089)
amf_12089$OTUID

#We lose no amf taxa through rarefaction. Now save this and strip the row name column in Excel

write.csv(amf_12089, "amf_12089.csv")

amf_12089<- read.csv("amf_12089.csv", row.names =1)

dim(amf_12089)

row.names(amf_12089)

#First transpose the table again

amf_12089_t <- t(amf_12089)

#Don't need to getrowsums or to rarefy since this table is already rarefied. Just Dive in.

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(amf_12089_t, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", ylim=c(0,300))

#now calculate richness with BioDiversityR package
#run function once to see what it does 
amf_richness <- estimateR(amf_12089_t)
#run function on transformed table to switch rows and columns
amf_richness<- t(estimateR(amf_12089_t))

#make a function to get estimates and make plots
estimates_plot_function <- function(datatrans,name){
  #plot S.chao vs S.obs to see if they are correlated
  estimates2<- t(estimateR(datatrans))
  pdf(paste("richnesscores_amf_richness_fullset_correlations_MiSeq",name,".pdf"))
  par(mfrow=c(1,2))
  plot(estimates2[,2] ~estimates2[,1], xlab="S.obs", ylab="chao", col=alpha("red", 0.5),pch=16)
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
estimates_plot_function(amf_12089_t,"amf_fullset")

##############################
#other ways to calculate species richness
##############################

# get species richness fo amf not rarefied
otu.H <- diversity(amf_12089_t) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(amf_12089_t, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
amf_12089_richness <- cbind(otu.richness,amf_richness)

write.csv(amf_12089_richness, "amf_fullset_richness.csv")

amf_12089_richness <- read.csv("amf_fullset_richness.csv", row.names=1)

#test if data are normal distirbuted, if significant not normal
shapiro.test(amf_12089_richness$S.obs) # horrible. (p=6.02e-10)
histogram(amf_12089_richness$S.obs) #skews left of center, but also splits in the middle (two dist, one left, smaller one right)
qqnorm(amf_12089_richness$S.obs)


#lets try transforming this to see if it improves

log_amf <- log(amf_12089_t)

#log doesnt work I don't think

sqrt_amf <- sqrt(amf_12089_t)

sqrt_amf_richness <- sqrt(amf_richness)

write.csv(sqrt_amf, "sqrt_amf.csv")

sqrt_amf <- read.csv("sqrt_amf.csv", row.names=1)

#run function
estimates_plot_function(sqrt_amf,"sqrt_amf")

##############################
#other ways to calculate species richness
##############################

# get species richness fo amf not rarefied
otu.H <- diversity(sqrt_amf) # Shannon entropy
otu.N1 <- exp(otu.H ) ## Shannon number of diversity
otu.N2 <- diversity(sqrt_amf, "inv") ## Simpson Diversity

#make data frane of shannon entropy, shannon diversity, simpson diversity
otu.richness <- data.frame(otu.H,otu.N1,otu.N2)

#add these to S obs, chao1, ACE
sqrt_amf_12089_richness <- cbind(otu.richness,sqrt_amf_richness)

write.csv(sqrt_amf_12089_richness, "sqrt_amf_fullset_richness.csv")

sqrt_amf_12089_richness <- read.csv("sqrt_amf_fullset_richness.csv", row.names=1)

#test if data are normal distirbuted, if significant not normal
shapiro.test(sqrt_amf_12089_richness$S.obs) # still horrible. (p=6.978e-09)
histogram(sqrt_amf_12089_richness$S.obs) #skews left of center, but also splits in the middle (two dist, one left, smaller one right)
qqnorm(sqrt_amf_12089_richness$S.obs)

########################################################################################
########Make some figures
##################################################################

dim(amf_12089_richness)

#read in dataframe with METADATA - add in info about cardinal directions and meters
metadata <- read.csv("SoberanesFiremetadata_ITS.csv", row.names =1)

dim(metadata)

#check that names match, if they dont match, use a matching function
row.names(amf_12089_richness) == row.names(metadata)

#bind it all together now
amfdata <- cbind(amf_12089_richness,metadata)

write.csv(amfdata, "amfdata_fullset.csv")
amfdata <- read.csv("amfdata_fullset.csv")

amfdata$Fire <- factor(amfdata$Fire, levels=c("Prefire","Postfire"))
amfdata$Burn <- factor(amfdata$Burn, levels=c("Unburned","Burned"))

#make a figure
amfbyplot <- ggplot(amfdata, aes(x=Plot, y=S.obs)) +   
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
amfbyplot


amfbyplot%>%
  ggexport(filename = "fullset_amfrichnesspreandpostfire_MiSeq.pdf")

#get mean, SD, SE of T1 and T2 by site by inoculum
# Calculate the means, sd, n, and se.

amf_plots <- ddply(amfdata, c("Plot","Fire"), summarise,
                   T1_mean = mean(S.obs, na.rm=TRUE),
                   T1_sd = sd(S.obs, na.rm=TRUE),
                   T1_n = sum(!is.na( S.obs)),
                   T1_se = T1_sd/sqrt(T1_n)
)

head(amf_plots)

#this is based on means but try to figure out with standard error
plot58percentreduction <- (100*(0.333-1.91667)/1.91667) #82.63% Increase
plot601percentreduction <- (100*(1.667-0.667)/1.667) #59.99% Decrease
plot603percentreduction <- (100*(1-0.2)/1) #80% Decrease

#Lets do the stats in the richness

#used linear mixed effect model from packaage nlme to run nested anova
#load library for nonlinear mixed effect models
library(nlme)
names(amfdata)

# species richness response variable
lmeY <- lme(S.obs ~ Fire*Burn, random = ~ 1 | Plot, data = amfdata)
model1 <- anova(lmeY)
model1 
summary(model1)
########################################################################################################
#---------------------------------Do stats on richness: Tukey ------------------------------------------
########################################################################################################

#lme on deadtanoak 2010 #sig effects of dead tanoak and plot by dead tanoak interaction
model1<- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = amfdata)
anova(model1)
tuk_T1 <- glht(model1, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T1 <- cld(tuk_T1) 
tuk.cld.T1
#plot(tuk.cld.T1, ylab="S obs fungi", col=c("cornflowerblue", "tomato1"))
plot(tuk.cld.T1, col=c("cornflowerblue", "tomato1"))


model2<- lme(S.obs ~ Burn, random = ~ 1 | Plot, data = amfdata)
anova(model2)
tuk_T2 <- glht(model2, linfct = mcp(Burn = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
plot(tuk.cld.T2, col=c("cornflowerblue", "tomato1"))


##separate burn vs unburned and test separately

Burnplots <- amfdata[which(amfdata$Burn=="Burned"), ]
unburnplot <- amfdata[which(amfdata$Burn=="Unburned"), ]

#do anova to test effect of fire on fungal richness in burned plots
model1_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = Burnplots)
burn_anova <- anova(model1_burn)
summary(burn_anova)

model2_burn <- lme(S.obs ~ Fire, random = ~ 1 | Plot, data = unburnplot)
anova(model2_burn)

#Tukey test to visualize which are sig different from each other
tuk_T2 <- glht(model1_burn , linfct = mcp(Fire = "Tukey")) 
tuk.cld.T2 <- cld(tuk_T2) 
tuk.cld.T2
a<- plot(tuk.cld.T2, ylab2="Burned Plots: Observed number of fungal species (OTUs)", xlab2="" ,col=c("cornflowerblue","tomato1")) 

a + facet_wrap(~Plot)
#plot(tuk.cld.T2 ,col=c("cornflowerblue","tomato1")) 

tuk_T3 <- glht(model2_burn, linfct = mcp(Fire = "Tukey")) 
tuk.cld.T3 <- cld(tuk_T3);tuk.cld.T3
plot(tuk.cld.T3, ylab2="Unburned Plots: Observed number of fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

pdf("amfrichness_prepostfire_Tukey_Miseq.pdf", width=8, height=5)
par(mfrow=c(1,2))
#plot(tuk.cld.T3, ylab="Unburned Plots: Observed no. fungal species (OTUs)", xlab="", col=c("cornflowerblue","tomato1")) 
plot.new()
plot(tuk.cld.T2, ylab2="Burned Plots: Observed no. fungal species (OTUs)", xlab2="", col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T3, col=c("cornflowerblue","tomato1")) 

plot(tuk.cld.T2, col=c("cornflowerblue","tomato1")) 
dev.off()

#set par bac to normal
par(mfrow=c(1,1))

#Now let's make the mean and SE figure like the other panels

E<-ggplot(amfdata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,-1,0),
               vjust=c(2.5,11,10.5,6.8,11.3,2.5),
               label = c("a","a","a","b","b","b"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Arbuscular Mycorrhizal Fungi") +
  theme_bw() + #make black and white
  ylab("Mean Number of Soil Fungal OTUs per Sample") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.text.y = element_text(size=20, angle=90,hjust=1),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled


E %>%
  ggexport(filename = "amf_richnesspreandpostfire_MiSeq_meanandSE_one_tukey.pdf")





#arrange using ggpubr package  (for example to put bacteria and fungal plots side by side

#Space to edit figures without having to run the whole code:

C<-ggplot(emfdata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,-1,0),
               vjust=c(18,8,16,17.8,24.5,11),
               label = c("a","b","a","c","b","d"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Ectomycorrhizal Fungi") +
  theme_bw() + #make black and white
  ylab("Mean Number of Soil Fungal OTUs per Sample") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.title.y = element_text(size=25),
        legend.text = element_text(size=22),
        axis.text.y = element_text(size=20, angle=90,hjust=1),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled


D<-ggplot(sapdata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,-1,0),
               vjust=c(6,28,27,18,25,7),
               label = c("a","a","a","b","c","c"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Saprotrophic Fungi") +
  theme_bw() + #make black and white
  #ylab("Mean Number of Soil Fungal OTUs per Sample") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.title.y=element_blank(),
        axis.text.y = element_text(size=20, angle=90,hjust=1),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled


E<-ggplot(amfdata, aes(x=Plot, y=S.obs, group=Fire, col=Fire,shape=Burn))+
  stat_summary(fun.y=mean,geom="point", size=6)+
  stat_summary(fun.data = mean_se,geom = "errorbar", size=1,
               alpha=0.7,position = position_dodge(0.01))+
  stat_summary(geom = 'text', fun.y = max, position = position_dodge(0), size = 8, hjust=c(0,0,0,0,-1,0),
               vjust=c(4,29.5,17,26,30.5,6.5),
               label = c("a","b","c","b","a","a"), col="black")+
  scale_color_manual(labels=c("Prefire","Postfire"), #manual labels for legend
                     values=c("cornflowerblue", "tomato1")) + #add in manual colors for points/lines
  ggtitle("Arbuscular Mycorrhizal Fungi") +
  theme_bw() + #make black and white
 # ylab("Mean Number of Soil Fungal OTUs per Sample") +  #change yaxis label
  theme(legend.position = "bottom", #put legend under graph
        text = element_text(size=18),
        axis.title.x=element_blank(), #remove Plot title
        axis.title.y=element_blank(),
        axis.text.y = element_text(size=20, angle=90,hjust=1),
        axis.text.x = element_text(size=20, angle=45, hjust=1))   #make x axis text larger and angled


#Then Bind them together

richnessfigures <- ggarrange(C,D,E, ncol = 3, nrow = 1, common.legend=TRUE,legend = "bottom" )

richnessfigures

#other option is to bind together fungi and bacteria together, add a column for bacteria v fungi, then make the plot in ggplot 2 and facet wrap by fungi vs bacteria

pdf("EMF_Sap_AMF_richness_meanandSE_3_column_2.pdf", height=15, width=20)
richnessfigures
dev.off()
