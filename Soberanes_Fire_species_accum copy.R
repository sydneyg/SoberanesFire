########Sepecies Accumulation Curves for Soberanes Fire#################

#set working directory
setwd("~/Desktop/Big_Sur_Project/R_Analysis")

#Reset R's Brain
rm(list=ls())


library(vegan)


##########Curves for Bacteria##################

bacteria_otuprefire <- read.csv("16S_e.10367_nozeroes_prefire.csv", row.names=1, check.names=FALSE)
bacteria_otuprefire <- as.matrix(bacteria_otuprefire)

bacteria_otupostfire<- read.csv("16S_e.10367_nozeroes_postfire.csv", row.names=1)
bacteria_otupostfire<- as.matrix(bacteria_otupostfire)

#make species accumulation curve to test if you are saturating community
specaccum1 <- specaccum(bacteria_otuprefire, method="exact")
plot(specaccum1 ,xlab="number of samples",ylab="number of OTUs", col="blue") #,ylim=c(0,2000))

specaccum1$richness

#?plot

#make species accumulation curve to test if you are saturating community
specaccum2 <- specaccum(bacteria_otupostfire , method="exact")
plot(specaccum2 ,xlab="number of samples",ylab="number of OTUs", col="red") #,ylim=c(0,2000))
specaccum2$richness
pdf("bacteria_species_accumulation_curves.pdf", height=5, width=6)
plot(specaccum1, 
     xlab="Number of samples", #make x axis label
     ylab="Bacterial OTUs", #make y axis label
     ylim=c(0,10000), #make y limits
     col="cornflowerblue", #color the line of the species accum curve
     ci.type = c("polygon"),  #change the SE to polygon
     ci.col=adjustcolor("cornflowerblue", alpha=0.5))#make the color of the polygon transparent
par(new=TRUE) #make it so next plot adds on to current plot
plot(specaccum2, #plot next thing
     col="tomato1", #make it a new color
     ylim=c(0,10000),  #use same y limits
     axes=FALSE, #supress axis ticks, numbers
     xlab="", ylab="", #suppress axis labels
     ci.type = c("polygon"), #make a polygon
     ci.col=adjustcolor("tomato1", alpha=0.3)) #make polygon transparent
legend("topleft", fill=c("cornflowerblue","tomato1"), bty="n",legend=c("Pre-fire Soil","Post-fire Soil"))
dev.off()

#############Curves for Fungi#####################

otuprefire <- read.csv("otu_ITS.e12089_nozeroes_prefire.csv", row.names=1)
otuprefire <- as.matrix(otuprefire)

otupostfire<- read.csv("otu_ITS.e12089_nozeroes_postfire.csv", row.names=1)
otupostfire <- as.matrix(otupostfire)

#make species accumulation curve to test if you are saturating community
specaccum3 <- specaccum(otuprefire, method="exact")
plot(specaccum3 ,xlab="number of samples",ylab="number of OTUs", col="blue" ,ylim=c(0,2000))
specaccum3$richness
#make species accumulation curve to test if you are saturating community
specaccum4 <- specaccum(otupostfire , method="exact")
plot(specaccum4 ,xlab="number of samples",ylab="number of OTUs", col="red",ylim=c(0,2000))

specaccum4$richness

pdf("fungi_species_accumulation_curves.pdf", height=5, width=6)
plot(specaccum3, 
     xlab="Number of samples", #make x axis label
     ylab="Fungal OTUs", #make y axis label
     ylim=c(0,2000), #make y limits
     col="cornflowerblue", #color the line of the species accum curve
     ci.type = c("polygon"),  #change the SE to polygon
     ci.col=adjustcolor("cornflowerblue", alpha=0.5))#make the color of the polygon transparent
par(new=TRUE) #make it so next plot adds on to current plot
plot(specaccum4, #plot next thing
     col="tomato1", #make it a new color
     ylim=c(0,4000),  #use same y limits
     axes=FALSE, #supress axis ticks, numbers
     xlab="", ylab="", #suppress axis labels
     ci.type = c("polygon"), #make a polygon
     ci.col=adjustcolor("tomato1", alpha=0.3)) #make polygon transparent
legend("topleft", fill=c("cornflowerblue","tomato1"), bty="n",legend=c("Pre-fire Soil","Post-fire Soil"))
dev.off()

  