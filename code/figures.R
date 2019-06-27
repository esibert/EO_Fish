#############################################################
#                       Figures For                         #
#        "Fish populations and diversity unaffected         #
#               by the Eocene-Oligocene Transition"         #
#                Sibert, Zill, Frigyik, Norris              #
#                                                           #
#       Compiled: 6/3/2019                                  #
#       Last updated; 6/3/2019                              #
#############################################################

library(stratigraph)
library(ichthyoliths)
library(doParallel)
library(viridis)
library(vegan)
library(RMark)
library(Hmisc)

source('code/functions.R')

load('eo_fish.RData')

##### Figure 2: The Fish Tooth Accumulation Figure ##### 

pdf(file = 'figures/Figure2_Fish.pdf', height = 7.5, width = 16, useDingbats = F)

par(mfrow = c(1, 9), 
    oma = c(8, 4, 4.5, 0.5) + 0.1,  #set outer margins for axis labels
    mar = c(0, 0, 0, 0) + 0.1)  # set plot martins to be very squished together
yaxis.age <- c(40, 28)
axis.scale <- 0.77
text.scale <- 1.2
pt.scale <- 1.8
iar.axis.text <- expression(paste('IAR (ich >38 ', mu, 'm/cm'^'2','/Myr)'))
eo.col <- adjustcolor('lightskyblue', alpha.f = 0.3)

# 1. 1406 
sub <- subset(iar.all, Site == 'IODP1406', select = c('Value', 'Age'))
sub <- sub[order(sub$Age),]
plot(sub, type = 'o', pch = 16, col = 'brown', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(1) # IAR values
mtext (text = iar.axis.text, side = 1, line = 2.5, cex = axis.scale)
axis(2) #age axis
mtext(text = 'Age (Ma)', side = 2, line = 2.5, cex = axis.scale)
text('IODP 1406', x = 8000, y = 40.2, cex = text.scale, font = 2)
abline(h = 33.9, col = eo.col, lwd = 5)

# 2. 886  
sub <- subset(iar.all, Site == 'ODP886', select = c('Value', 'Age'))
plot(sub, type = 'o', pch = 16, col = 'red', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(3) # IAR values
mtext (text = iar.axis.text, side = 3, line = 2.5, cex = axis.scale)
text('ODP 886', x = 1200, y = 27.8, cex = text.scale, font = 2)
abline(h = 33.9, col = eo.col, lwd = 5)


# 3. 1217  
sub <- subset(iar.all, Site == 'ODP1217', select = c('Value', 'Age'))
plot(sub, type = 'o', pch = 16, col = 'chocolate', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(1) # IAR values
mtext (text = iar.axis.text, side = 1, line = 2.5, cex = axis.scale)
text('ODP 1217', x = 3000, y = 40.2, cex = text.scale, font = 2)
abline(h = 33.9, col = eo.col, lwd = 5)

# 4. 522  
sub <- subset(iar.all, Site == 'DSDP522', select = c('Value', 'Age'))
plot(sub, type = 'o', pch = 16, col = 'plum3', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(3) # IAR values
mtext (text = iar.axis.text, side = 3, line = 2.5, cex = axis.scale)
text('DSDP 522', x = 6000, y = 27.8, cex = text.scale, font = 2)
abline(h = 33.9, col = eo.col, lwd = 5)

# 5. 596  
sub <- subset(iar.all, Site == 'DSDP596', select = c('Value', 'Age'))
plot(sub, type = 'o', pch = 16, col = 'mediumorchid4', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(1) # IAR values
mtext (text = iar.axis.text, side = 1, line = 2.5, cex = axis.scale)
text('DSDP 596', x = 410, y = 40.2, cex = text.scale, font = 2)
abline(h = 33.9, col = eo.col, lwd = 5)

# 6. 689  
sub <- subset(iar.all, Site == 'ODP689', select = c('Value', 'Age'))
plot(sub, type = 'o', pch = 16, col = 'royalblue3', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(3) # IAR values
mtext (text = iar.axis.text, side = 3, line = 2.5, cex = axis.scale)
text('ODP 689', x = 2700, y = 27.8, cex = text.scale, font = 2)
abline(h = 33.9, col = eo.col, lwd = 5)

# 7. 748  
sub <- subset(iar.all, Site == 'ODP748', select = c('Value', 'Age'))
plot(sub, type = 'o', pch = 16, col = 'dodgerblue', ylim = yaxis.age, cex = pt.scale, 
     xlim = c(0, max(sub$Value)), 
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis(1) # IAR values
mtext (text = iar.axis.text, side = 1, line = 2.5, cex = axis.scale, font = 2)
text('ODP 748', x = 400, y = 40.2, cex = text.scale)
abline(h = 33.9, col = eo.col, lwd = 5)

# 8. Deep ocean temperature
plot(temp_cramer$temp.cramer, temp_cramer$age.temp.cramer.2012, 
     type = 'l', lwd = 3, ylim = yaxis.age,
     xlim = c(4, 11),
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis (3)
mtext (text = 'Deep Ocean Temp (C)', side = 3, line = 2.5, cex = axis.scale)
abline(h = 33.9, col = eo.col, lwd = 5)

# 9. Relative Sea Level
plot(sealevel_cramer$sealevel.cramer, sealevel_cramer$age.sealevel.cramer.2012, 
     type = 'l', lwd = 3, ylim = yaxis.age,
     xlim = c(-10, 60),
     bty = 'n', axes = FALSE, xlab = '', ylab = '')
axis (1)
mtext (text = 'Sea Level (m)', side = 1, line = 2.5, cex = axis.scale)
abline(h = 33.9, col = eo.col, lwd = 5)


dev.off()



##### Figure 3: Range Charts #####
### normalized range chart plot 

#### Parameters and code for *percentage* figures
colors.vector <- viridis(5)
colors.vector[1] <- 'gray70' 
splits <- c(3, 6, 9, 12) #percentage splits 
largesize <- 1


############# Parameters and code for *counts* based figures
splits <- c(1, 3, 5, 10) # Counts splits - use this set! 
## define sizes for legend cex
sizes<-c(0:length(splits)) 
sizes<-((sizes/(max(sizes)/largesize)) + largesize) #evenly distributed betwen cex = largesize and cex = largesize*2

### 689
pdf(file='figures/RangeChart689_counts.pdf', height = 4, width = 10, useDingbats = FALSE)
par(mar=c(5, 4, 2, 2))

rangechart(counts.689, reorder = 'fad.by.lad', xaxis.labels = 'alphanum', print.xaxis = FALSE, 
           yaxis.ticks = FALSE, normalize.counts = FALSE, count.breaks = c(0,splits), 
           col.points = 'by.count', cols.vec = colors.vector,
           pch.points = 16, cex.points = 'by.count', largesize = largesize, main = 'ODP Site 689')

# Counts legend
legend('bottomright', legend=c('1', '2-3', '4-5', '6-10', '11+'), pch=16, 
       col=colors.vector, pt.cex=sizes, cex=0.8,
       horiz = TRUE, bg = 'n', bty = 'n') 

## Add EO line
abline(h=33.9, col='gray')

dev.off()

### 596
pdf(file='figures/RangeChart596_counts.pdf', height = 5, width = 10, useDingbats = FALSE)
par(mar=c(5, 4, 2, 2))

rangechart(counts.596, reorder = 'fad.by.lad', xaxis.labels = 'alphanum', print.xaxis = FALSE, 
           yaxis.ticks = FALSE, normalize.counts = FALSE, count.breaks = c(0,splits),
           col.points = 'by.count', cols.vec = colors.vector,
           pch.points = 16, cex.points = 'by.count', largesize = largesize, main = 'DSDP Site 596')

# Counts legend
legend('bottomright', legend=c('1', '2-3', '4-5', '6-10', '11+'), pch=16, 
       col=colors.vector, pt.cex=sizes, cex=0.8,
       horiz = TRUE, bg = 'n', bty = 'n') 

## Add EO line
abline(h=33.9, col='gray')

dev.off()

##### Percentage-based figures #####

## define legend values
splits.legend<-paste(c(0,splits), '-', c(splits,100), '%', sep='')
splits.legend[1]<-paste('<', splits[1], '%', sep='')
splits.legend[length(splits.legend)] <- paste('>', splits[length(splits)], '%', sep='')

## define sizes for legend cex
sizes<-c(0:length(splits)) 
sizes<-((sizes/(max(sizes)/largesize)) + largesize) #evenly distributed betwen cex = largesize and cex = largesize*2


##### DSDP 596
morphdat.596 <- morphdat.combined.596
# morphdat.596 <- morphdat.all.596
sub.596 <- data.frame(morphdat.596$AgeID, morphdat.596$Morphotype_Name)
counts.596 <- table(sub.596, exclude = '')

### Make the actual figure 
pdf(file='figures/RangeChart596_cleaned.pdf', height = 5, width = 10, useDingbats = FALSE)
par(mar=c(5, 4, 2, 2))

## Make the chart
rangechart(counts.596, reorder = 'fad.by.lad', xaxis.labels = 'alphanum', print.xaxis = FALSE, 
            yaxis.ticks = FALSE, normalize.counts = TRUE, count.breaks = c(0,splits), #count.breaks = NULL, 
            col.points = 'by.count', cols.vec = colors.vector,
            pch.points = 16, cex.points = 'by.count', largesize = largesize, main = 'DSDP Site 596')

## Add the legend
legend('bottomright', legend=splits.legend, pch=16, 
       col=colors.vector, pt.cex=sizes, cex=0.8,
       horiz = TRUE, bg = 'n', bty = 'n') 

## Add EO line
abline(h=33.9, col='gray')

dev.off()

## Get the number to type list in the console
typelist.596 <- rangechart(counts.596, reorder = 'fad.by.lad', xaxis.labels = 'alphanum', print.xaxis = TRUE, 
           yaxis.ticks = FALSE, normalize.counts = TRUE, count.breaks = c(0,splits), #count.breaks = NULL, 
           col.points = 'by.count', cols.vec = colors.vector,
           pch.points = 16, cex.points = 'by.count', largesize = largesize)


##### ODP Site 689
### normalized range chart plot 

## select dataset
morphdat.689 <- morphdat.combined.689
# morphdat.596 <- morphdat.all.596
sub.689 <- data.frame(morphdat.689$AgeID, morphdat.689$Morphotype_Name)
counts.689 <- table(sub.689, exclude = '')


### Make the actual figure 
pdf(file='figures/RangeChart689_cleaned.pdf', height = 5, width = 10, useDingbats = FALSE)
par(mar=c(5, 4, 2, 2))

## Make the chart
rangechart(counts.689, reorder = 'fad.by.lad', xaxis.labels = 'alphanum', print.xaxis = FALSE, 
           yaxis.ticks = FALSE, normalize.counts = TRUE, count.breaks = c(0,splits), #count.breaks = NULL, 
           col.points = 'by.count', cols.vec = colors.vector,
           pch.points = 16, cex.points = 'by.count', largesize = largesize, main = 'ODP Site 689')

## Add the legend
legend('bottomright', legend=splits.legend, pch=16, 
       col=colors.vector, pt.cex=sizes, cex=0.8,
       horiz = TRUE, bg = 'n', bty = 'n') 

## Add EO line
abline(h=33.9, col='gray')

dev.off()



## Get the number to type list in the console
typelist.689 <- rangechart(counts.689, reorder = 'fad.by.lad', xaxis.labels = 'alphanum', print.xaxis = TRUE, 
                       yaxis.ticks = FALSE, normalize.counts = TRUE, count.breaks = c(0,splits), #count.breaks = NULL, 
                       col.points = 'by.count', cols.vec = colors.vector,
                       pch.points = 16, cex.points = 'by.count', largesize = largesize)



