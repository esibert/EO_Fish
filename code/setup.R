#############################################################
#               Setup Code for analysis of                  #
#               Eocene/Oligocene Fish Teeth                 #
#       For "No state change in pelagic fish production     #
#           and biodiversity during the Eocene-Oligocene    #
#           Transition"                                     #
#       by Elizabeth C. Sibert, Michelle E. Zill,           #
#           Ella T. Frigyik, and Ricahrd Norris             #
#                                                           #
#       Code written and maintained by Elizabeth Sibert     #
#           Contact: esibert@fas.harvard.edu                #
#                                                           #
#       Compiled: 5/6/2019                                  #
#       Last updated; 12/27/2019                            #
#############################################################

## This file contains the 'read-in' for all datasets used in this manuscript: 
# 1. The Ichthyolith Accumulation Rate datasets, alongside Barium accumulation and Silica accumulation compiled into a single CSV file, for sites DSDP 522, DSDP 596, ODP 689, ODP 748, ODP 886, ODP 1217, and IODP 1406
# 2. The tooth morphology codes and IDs for DSDP Site 596 
# 3. The tooth morphology codes and IDs for ODP 689 (coming soon)
# 4. A "unique morphotypes" CSV with character codes for each individual morphotype in this dataset.

###### Libraries, Function, RData files  #####
##if you don't have the ichthyoliths package yet
# # install and load the package: devtools:
# install.package('devtools')
# library(devtools)
# 
# # install ichthyoliths package from github (this requires the devtools package):
# install_github('esibert/ichthyoliths')

## Install stratigraph from binary
# install.packages('../abandoned_package_binaries/stratigraph_0.66.tar.gz', repos = NULL)

library(ichthyoliths)
library(doParallel)
library(viridis)
library(vegan)
library(RMark)
library(Hmisc)

# library(stringr)


## In-house functions
source('code/functions.R')

## Save/load R workspace image to work with 
# save.image('eo_fish.RData')
# load('eo_fish.RData')


## Restore graphical parameters
# par(mfrow=c(1,1))
# par(mar=c(5.1, 4.1, 4.1, 2.1))

##### Define variables used across figures #####

# Define x-axis plot age constraints (for all figures)
xlim.ages<-c(29,43)

# Define sites
sites<-c('IODP1406', 'ODP1217', 'DSDP522', 'DSDP596', 'ODP748', 'ODP689', 'ODP886')

# moving average filter
ma5<-c(1,1,1,1,1)/5
ma3 <- c(1,1,1)/3
ma4 <- c(1,1,1,1)/4

##### Read in and clean up datasets #####

### 1. IAR, BAR, SiAR datasets
all_AR_data<-read.csv('data/EOT_dataset_comp.csv', header=T)

iar.all<-subset(all_AR_data, Proxy=='IAR')


### 2. Morphotypes from DSDP Site 596

## call in and clean up dataset
dat.596 <- read.csv('data/DSDP_596_EO_fish.csv', header = TRUE)
morphdat.596 <- toothdat.cleanup(dat.596, sortby = 'age')

# Combine pairs of samples for larger sample sizes. These are the unique Age values that are being combined 
c01 <- c(28.49, 28.74)
c02 <- c(29.44, 29.75)
c03 <- c(30.52, 30.80)
c04 <- c(32.42, 32.68)
c05 <- c(33.52, 33.78)
c06 <- c(35.09, 35.34)
c07 <- c(36.34, 36.59)
c08 <- c(38.99, 39.42)
combines.596 <- list(c01, c02, c03, c04, c05, c06, c07, c08)
rm(c01, c02, c03, c04, c05, c06, c07, c08) #clean up!

morphdat.combined.596 <- combine.samples(morphdat = morphdat.596, combines = combines.596)
morphdat.all.596 <- morphdat.596
rm(morphdat.596) #forces me to choose which morphdat (combined or all) to use for everything


### 3. Morphotypes from ODP 689

## call in and clean up dataset
dat.689 <- read.csv('data/ODP_689_EO_fish_csv.csv', header = TRUE)
morphdat.689 <- dat.689[order(dat.689[,2]), ] 

AgeID.unique.689 <- unique(morphdat.689$AgeID)

c01 <- AgeID.unique.689[1:10]
c02 <- AgeID.unique.689[11:20]
c03 <- AgeID.unique.689[21:30]
c04 <- AgeID.unique.689[31:40]
c05 <- AgeID.unique.689[41:50]
c06 <- AgeID.unique.689[51:60]
c07 <- AgeID.unique.689[61:70]
c08 <- AgeID.unique.689[71:80]
c09 <- AgeID.unique.689[81:90]
c10 <- AgeID.unique.689[91:100]

combines.689 <- list(c01, c02, c03, c04, c05, c06, c07, c08, c09,c10)
rm(c01, c02, c03, c04, c05, c06, c07, c08, c09, c10) #clean up!

morphdat.combined.689 <- combine.samples(morphdat = morphdat.689, combines = combines.689)
morphdat.all.689 <- morphdat.689
rm(morphdat.689)


### 4. Sea Level and Temperature from Cramer 2011
# A note on age models: These datasets are taken from Cramer et al 2011, tables T7 (0-9 Ma) and T4 (9.2-108 Ma)
# https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JC007255 is the reference for this

# There is a break in ages from 83.45 to 83.3 because of age model conversion between GTS 2004 (what Cramer reports) and GTS 2012 (what we are using in this study). Used ODSN Timescale converter http://www.odsn.de/odsn/services/conv_ts/conv_ts.html to convert ages. 

temp_cramer <- read.csv('data/temp-cramer.csv', header = TRUE, check.names = TRUE)
names(temp_cramer)[1] <- 'age.temp.cramer.2012'
sealevel_cramer <- read.csv('data/sealevel-cramer.csv', header = TRUE)
names(sealevel_cramer)[1] <- 'age.sealevel.cramer.2012'

##### Disparity calculations #####
# ichthyolith disparity weights/traits for disparity analyses
traitset <- ichthyoliths::traits
weights <- ichthyoliths::weights
##### DSDP 596 Disparity #####

## select dataset
# morphdat.596 <- morphdat.combined.596
# morphdat.596 <- morphdat.all.596

## Calculate disparity using ichthyoliths package
toothdist.596 <- distances_clust(morph = morphdat.596, traits = traitset, weights = weights, 
                                 morphCols = c(8:33), traitsStartCol = 8, subsetWeights = TRUE, IDCol = 4, 
                                 contTraits = FALSE)

## test for NAs
#if there are any NA values, look at which teeth are problems, 
# check your original spreadsheets, and fix them, then re-run the toothdist line above
df.na <- subset(toothdist.596, is.na(toothdist.596$dist.sum))
unique(df.na$ObjectA) # list of teeth that have a typo of some sort

## save the dist-pairs output as a csv file for later
write.csv(toothdist.596, 'data/pairwisedist_morph_596.csv')
## to call this in: 
# saved.toothdist.596 <- read.csv('data/pairwisedist_morph_596.csv', header = TRUE)
# saved.toothdist.596 <- saved.toothdist.596[,-1] #remove first column, was rownames


## once there are no more NA values, generate distance matrix # 
distmat.596 <- distmat(toothdist.596)

## save the toothdistmat output as a csv file for later as well
write.csv(distmat.596, 'data/distmat_596.csv')
## to call this back in:
# saved.toothdistmat.596
# saved.distmat.596 <- read.csv('data/distmat_596.csv', header = FALSE)
# saved.distmat.596 <- saved.distmat.596[,-1] #remove first column
# saved.distmat.596 <- saved.distmat.596[-1,] # remove first row

# NOTE: The order that the toothdist.596 object is in is based on *age* and not on original row from the input csv file. This is due to the fact that Ella coded the samples out of order (this is a good thing). Because of this, the rownames are numerically strange, but can easily be traced back to the dat596 object. 


##### ODP 689 Disparity #####

## select dataset
# morphdat.689 <- morphdat.combined.689
# morphdat.689 <- morphdat.all.689


## Calculate disparity using ichthyoliths package
toothdist.689 <- distances_clust(morph = morphdat.689, traits = traitset, weights = weights, 
                                 morphCols = c(8:33), traitsStartCol = 8, subsetWeights = TRUE, IDCol = 4, 
                                 contTraits = FALSE)

## test for NAs
#if there are any NA values, look at which teeth are problems, 
# check your original spreadsheets, and fix them, then re-run the toothdist line above
df.na <- subset(toothdist.689, is.na(toothdist.689$dist.sum))
unique(df.na$ObjectA) # list of teeth that have a typo of some sort

## save the dist-pairs output as a csv file for later
write.csv(toothdist.689, 'data/pairwisedist_morph_689.csv')
## to call this in: 
# saved.toothdist.689 <- read.csv('data/pairwisedist_morph_689.csv', header = TRUE)
# saved.toothdist.689 <- saved.toothdist.689[,-1] #remove first column, was rownames


## once there are no more NA values, generate distance matrix # 
distmat.689 <- distmat(toothdist.689)

## save the toothdistmat output as a csv file for later as well
write.csv(distmat.689, 'data/distmat_689.csv')
## to call this back in:
# saved.toothdistmat.689
# saved.distmat.689 <- read.csv('data/distmat_689.csv', header = FALSE)
# saved.distmat.689 <- saved.distmat.689[,-1] #remove first column
# saved.distmat.689 <- saved.distmat.689[-1,] # remove first row


