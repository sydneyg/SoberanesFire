#Soberanes Indicator Species Analysis

#12/6/2019

#Clear R

rm(list=ls())

#set working directory
setwd("C:/Users/enrig/Desktop/R_Analysis")


#Load in Libraries

library(vegan)
library(cluster)
library(ggplot2)
library(dplyr)
library(ggdendro)
library(dendextend)
library(indicspecies)


#Load in Rarified 16S Dataset

OTU_16S = read.csv("otu_16S.e10367.csv", row.names = 1, check.names = FALSE)

#check dim of OTU table

dim(OTU_16S)

rownames(OTU_16S)

head(OTU_16S)
#Load in Metadata File

bacteria_meta = read.csv("SoberanesFiremetadata_16S.csv", row.names = 1, check.names = FALSE)

#check dim of metadata

dim(bacteria_meta)

#Attempt indicator species analysis, checking for Burned

burn.ind = multipatt(OTU_16S, bacteria_meta~Burn, func="IndVal.g", control = how(nperm=999))
summary(burn.ind)

#Checking the Dune dataset (which works) to figure out why mine isn't working

rownames(dune)
rownames(dune.env)

rownames(bacteria_meta)

#attempting to rerun the function with sample names removed

otu_post_bac = read.csv("otu_16s.e10367_post_fire.csv", row.names = 1, check.names = FALSE)
dim(otu_post)

otu_post_no_names = otu_post[ ,2:14883]
dim(otu_post_no_names)
rownames(otu_post_no_names)
head(otu_post_no_names)
colnames(otu_post_no_names)

otu_post_no_names[36,1]

meta_post_bac = read.csv("SoberanesFiremetadata_16s_post_fire.csv", row.names =1, check.names = FALSE)
dim(meta_post)
meta_post_no_names = meta_post[ ,2:5]
head(meta_post_no_names)


bac_post_fire.ind = multipatt(otu_post_bac, meta_post_bac$Burn, func="IndVal.g", control = how(nperm=999))

summary(bac_post_fire.ind, alpha = 0.001)

sum1 <- summary(bac_post_fire.ind, alpha = 0.001)

capture.output(summary(bac_post_fire.ind, alpha = 0.001), file= "bacteria_postfire_ISA.doc")
#Now let's repeat with the bacteria prefire dataset

otu_pre_bac = read.csv("otu_16S.e10367_pre_fire.csv", row.names = 1, check.names = FALSE)
dim(otu_pre_bac)

meta_pre_bac = read.csv("SoberanesFiremetadata_16s_pre_fire.csv", row.names =1, check.names = FALSE)
dim(meta_pre_bac)

bac_pre_fire.ind = multipatt(otu_pre_bac, meta_pre_bac$Burn, func="IndVal.g", control = how(nperm=999))
summary(bac_pre_fire.ind, alpha = 0.001)

capture.output(summary(bac_pre_fire.ind, alpha = 0.001), file= "bacteria_prefire_ISA.doc")
#Time for post-fire Fungi

otu_post_fun = read.csv("otu_ITS_rarified_postfire.csv", row.names = 1)
dim(otu_post_fun)

meta_post_fun = read.csv("ITS_meta_postfire.csv", row.names = 1)
dim(meta_post_fun)

fun_post_fire.ind = multipatt(otu_post_fun, meta_post_fun$Burn, func="IndVal.g", control = how(nperm=999))
summary(fun_post_fire.ind, alpha = 0.001)

capture.output(summary(fun_post_fire.ind, alpha = 0.001), file= "fungi_postfire_ISA.doc")

#Repeat for pre-fire fungi

otu_pre_fun = read.csv("otu_ITS_rarified_pre_fire.csv", row.names = 1)
dim(otu_pre_fun)

meta_pre_fun = read.csv("ITS_meta_prefire.csv", row.names = 1)
dim(meta_pre_fun)

fun_pre_fire.ind = multipatt(otu_pre_fun, meta_pre_fun$Burn, func ="IndVal.g", control = how(nperm=999))
summary(fun_pre_fire.ind, alpha = 0.001)

capture.output(summary(fun_pre_fire.ind, alpha = 0.001), file = "fungi_prefire_ISA.doc")
