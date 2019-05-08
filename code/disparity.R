#########################################################
#                                                       #
#           Disparity analyses for 596 and 689          #
#               last updated 5/8/2019                   #
#                                                       #
#########################################################



##### 1. Select the distmat to run all of these calculations/figures/etc. on.
## Note that this whole script is *not* site-specific, and analyses are not saved as 596 or 689 until the end. 

toothdistmat <- distmat.596
ages.combined <- morphdat.combined.596$AgeID
ages.all <- morphdat.all.596$AgeID

##### 2.  Calculate NMDS values for plotting 
NMDS3 <- metaMDS(toothdistmat, k=3, try = 40)
stressplot(NMDS3)
ordiplot(NMDS3, type = 'text')

pts <- as.data.frame(NMDS3$points)
pts$ages.all <- ages.all
pts$ages.combined <- ages.combined

##### 3. Make timeseries NMDS figures #####

### a. plot with each time slice in morphospace
# define ages
ages <- unique(pts$ages.combined)
# make plot
plot(pts[,1], pts[,2], type='n', xlab='MDS1', ylab='MDS2')
cols<-rainbow(length(ages))
# cols<-colfunc(length(AgeID.unique))
for(k in 1:length(ages)) {
    sub<-subset(pts, ages.combined==ages[k])
    points(sub[,1], sub[,2], pch=16, col=cols[k])
    Plot_ConvexHull(sub[,1], sub[,2], lcol=cols[k], lwd=2)
}

### b. plot with just before and after EOT in morphospace
# define ages
ages <- unique(pts$ages.combined)
# make plot
plot(pts[,1], pts[,2], type='n', xlab='MDS1', ylab='MDS2')
# add before E/O
sub <- subset(pts, ages.combined >=33.9)
points(sub[,1], sub[,2], pch=16, col='red', cex = 2)
Plot_ConvexHull(sub[,1], sub[,2], lcol='red', lwd=2)
# add after E/O
sub <- subset(pts, ages.combined <=33.9)
points(sub[,1], sub[,2], pch=16, col='blue')
Plot_ConvexHull(sub[,1], sub[,2], lcol='blue', lwd=2)

### c. plot with each time bin into its own plot (18 time points)
ages <- unique(ages.all)
par(mfrow = c(3, 6))
cols = c(rep('blue', 10), rep('red', 8))
for(k in 1:length(ages)) {
    sub<-subset(pts, ages.all==ages[k])
    plot(sub[,1], sub[,2], pch=16, col=cols[k], main = ages[k])
    Plot_ConvexHull(sub[,1], sub[,2], lcol=cols[k], lwd=2)
}
#reset the workspace to have only one plotting window
par(mfrow = c(1,1))


# separate each time bin into its own plot (10 time points)
ages <- unique(ages.combined)
par(mfrow = c(2, 5))
cols = c(rep('blue', 5), rep('red', 5))
for(k in 1:length(ages)) {
    sub<-subset(pts, ages.combined==ages[k])
    plot(sub[,1], sub[,2], pch=16, col=cols[k], main = ages[k])
    Plot_ConvexHull(sub[,1], sub[,2], lcol=cols[k], lwd=2)
}
#reset the workspace to have only one plotting window
par(mfrow = c(1,1))
