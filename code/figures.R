#############################################################
#                       Figures For                         #
#        "Fish populations and diversity unaffected         #
#               by the Eocene-Oligocene Transition"         #
#                Sibert, Zill, Frigyik, Norris              #
#                                                           #
#       Compiled: 5/6/2019                                  #
#       Last updated; 5/6/2019                              #
#############################################################



##### Figure 1: Ichthyolith Accumulation rates (only) #####

## make sites vector with the appropriate order for the sites from N to S, for running the loop 
sites2 <- c("IODP1406", "ODP886", "ODP1217",  "DSDP522",  "DSDP596", "ODP689",  "ODP748") 
sites2.names <- c("Subarctic Atlantic (IODP 1406)", "North Pacific (ODP 886)", "Equatorial Pacific (ODP 1217)",  "South Atlantic (DSDP 522)",  "South Pacific (DSDP 596)", "Antarcic (ODP 689)",  "Antarctic (ODP 748)")

# not sure why I have a different age vector here but going with it for the moment... 
xlim.ages2 <- c(28, 40.5)

pdf('figures/fish_only_EOT.pdf', height = 10, width = 7, useDingbats = FALSE)
par(mfrow = c(7,1))
par(mar = c(2, 5, 1.5, 2))
for(i in 1:length(sites2)) {
    dat <- subset(iar.all, Site == sites2[i], select = c('Age', 'Value'))
    plot(dat, pch=16, col='black', cex=1, bty = 'n', xlab = '', ylab = '', axes = F, xlim = xlim.ages2)
    axis(2)
    mtext(paste(sites2.names[i]), side = 3, line = 0, font = 2, cex = 0.9) 
    mtext('IAR', side=2, line=3.2, cex=0.8)
    mtext(expression(paste('ich/cm'^2,'/kyr')), side=2, line=1.8, cex=0.6)
    lines(dat$Age, filter(dat$Value, ma5), lwd=3, col='black')
}
axis(1)
box(bty = 'l')
abline (v = 33.9)
dev.off()



##### Figure 2: Orders of Magnitude IAR values #####

#Variables for plotting and parameters
fillvec<-c('blue', 'magenta', 'gray', 'red', 'green3', 'yellow', 'purple')
outvec<-c('black', 'black', 'black', 'black', 'black', 'black', 'black')
pchvec<-c(24, 22, 21, 25, 21, 23, 25)

#plot
pdf('figures/magnitudes_log.pdf', height=7.5, width=7.5, useDingbats = F)
par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 2.1, 3.1))
plot(iar.all$Age, iar.all$Value, type='n', log='y', 
     xlim=xlim.ages, ylim = c(50, max(iar.all$Value)), 
     xlab='Age (Ma)', ylab='', 
     cex=1.3, axes=F)
axis(1)
#setup log axis
at.y <- outer(1:9, 10^(0:4))
lab.y <- ifelse(log10(at.y) %% 1 == 0, at.y, NA)
axis(2, at=at.y, labels=lab.y, las=1)
abline(h=800, lty=2, lwd=2)
abline(v=33.9, lty = 1, lwd = 6, col = 'lightblue')
box()
mtext('IAR (log scale)', side=2, line=3, cex=1)
mtext(expression(paste('ich/cm'^2,'/kyr')), side=2.2, line=2, cex=0.8)

# #add 5 pt moving averages (scaled 886)
# for(i in 1:length(sites)) {
#    lines(subset(iar.all, Site==sites[i])$Age, filter(subset(iar.all, Site==sites[i])$Scaled, ma5), 
#       pch=pchvec[i], col=adjustcolor(fillvec[i], alpha.f=1), lwd=5)
# }

#add points
for(i in 1:length(sites)) {
    points(x = subset(iar.all, Site==sites[i])$Age, y = subset(iar.all, Site==sites[i])$Value, 
           pch=pchvec[i], col=outvec[i], bg=fillvec[i], cex=1.3)
}
# legend('topright',legend=c('N. Atlantic', 'Eq. Pacific', 'S. Atlantic', 'S. Pacific', 'Antarctic 
#    (Kerguelen)', 'Antarctic 
#    (Maude)', 'N. Pacific'),
#    col=outvec, pt.bg=fillvec, pch=pchvec, bg='white', ncol=2)
legend('topright',legend=sites,col=outvec, pt.bg=fillvec, pch=pchvec, bg='white', ncol=2)
dev.off()



##### Figure 3: Range Chart for 596 #####
### normalized range chart plot 

## select dataset
morphdat.596 <- morphdat.combined.596
# morphdat.596 <- morphdat.all.596

## make SC objects necessary for plotting
sc.596<-build.strat.obj(morphdat.596)
ad.596<-a.datums(sc.596, depths=sc.596$absolute.ages, increasing.down = TRUE)


## Make splits for plot, normalize the range chart for plotting
splits<-c(3,6,9,12)
sc.plot <- breaks.fn (sc.596, splits)

## define colors
cols.vec <- rev(viridis(5))
cols.vec[1] <- 'gray70'

## define legend input
splits.legend<-paste(c(0,splits), '-', c(splits,100), '%', sep='')
splits.legend[1]<-paste('<', splits[1], '%', sep='')
splits.legend[length(splits.legend)] <- paste('>', splits[length(splits)], '%', sep='')

## Make the actual figure
pdf(file='figures/RangeChart596.pdf', height = 7, width = 10, useDingbats = FALSE)
par(mar=c(12, 4, 2, 4))
rangechart3(sc.plot, reorder='lad.by.fad', depths=sc.plot$absolute.ages, cex.xaxis=0.5, 
            cex.yaxis=1, cex.points="by.count", llwd=1, llcol = 'lightgray', llty = 3,
            col.points="by.count", colors.vec = cols.vec, xaxis.labels = 'names', 
            #xaxis.labels = c('numbers', 'names', 'alphanum')
            pch.points=16, baselines=FALSE, large.size=1, count.group=FALSE, legend.loc = 'bottomright', 
            legend.values = splits.legend, legend.horiz = TRUE, return.xaxis = FALSE) 
abline(h=33.9, col='gray')
dev.off()




##### Figure S1: Barium/Silica/IAR stack figure #####
pdf('figures/eot_IAR_ba_si.pdf', height=10, width=7, useDingbats = F)

    #set up window:
    par(mfrow=c(7,1))
    
    # Make plots
    plot.ba.si(all_data = all_AR_data, dataset='IODP1406', xlim.ages=xlim.ages, 
               plottitle = 'Subarctic Atlantic IODP U1403', plottitle.loc = c(30.7, 9700), age.axis=F, 
               #iar.ylim=c(0,2000),
               plot.margins=c(0,5,1,8))
    plot.ba.si(all_data = all_AR_data, dataset='ODP886', xlim.ages=xlim.ages, 
               plottitle = 'North Pacific ODP 886', plottitle.loc = c(30.3, 1300), age.axis=F, 
               use.scaled.IAR = FALSE, 
               plot.margins=c(0,5,1,8))
    plot.ba.si(all_data = all_AR_data, dataset='ODP1217', xlim.ages=xlim.ages, 
               plottitle = 'Equatorial Pacific ODP 1217', plottitle.loc = c(30.7, 4000), 
               si.units='Opal accumulation (g/cm2/myr)', ba.units='Barium accumulation (mg/cm2/kyr)', 
               age.axis=F)
    plot.ba.si(all_data = all_AR_data, dataset='DSDP522', xlim.ages=xlim.ages, 
               plottitle = 'Subtropical Atlantic DSDP 522', plottitle.loc = c(30.9, 7600), 
               si.units='Opal accumulation (g/cm2/myr)', 
               age.axis=F)
    plot.ba.si(all_data = all_AR_data, dataset='DSDP596', xlim.ages=xlim.ages, 
               plottitle = 'Subtropical Pacific DSDP 596', plottitle.loc = c(30.9, 430), 
               si.units='SiO2 accumulation (g/cm2/myr)', ba.units='Barium accumulation (mg/cm2/myr)', 
               age.axis=F)
    plot.ba.si(all_data = all_AR_data, dataset='ODP689', xlim.ages=xlim.ages, 
               plottitle = 'Antarctic ODP 689', plottitle.loc = c(30, 2600), 
               si.units='SiO2 accumulation (g/cm2/myr)', ba.units='Barium accumulation (mmol/cm2/myr)', 
               age.axis=F)
    plot.ba.si(all_data = all_AR_data, dataset='ODP748', xlim.ages=xlim.ages, 
               plottitle = 'Antarctic ODP 748', plottitle.loc = c(30, 500), 
               si.units='Si Accumulation (g/cm2/myr)',  
               age.axis=T, 
               plot.margins=c(4,5,0,8))
    box(bty='l')
    abline(v = 33.9)
    #par(mar=c(4,5,0,8))
    #plot(zachos$age_zachos, filter(zachos$O18_adj, ma5), type='l', xlim=xlim.ages, ylim=c(3,0.4), 
    #   bty='l', xlab='', ylab='')
    mtext('Age (Ma)', side=1, line=2)
    #mtext(expression(paste(delta^18,'O')), side=2, line=2.2, font=2)
    #text(30.3, 0.9, labels=expression(paste(delta^18,'O (Zachos et al 2008)')), font=2)
    #close the PDF device
dev.off()


##### Figure S2: Morphometrics for 596 #####

##### Figure: Range chart for 689 #####
### normalized range chart plot 

## select dataset
# morphdat.689 <- morphdat.combined.689
# morphdat.689 <- morphdat.all.689

## make SC objects necessary for plotting
sc.689<-build.strat.obj(morphdat.689)
ad.689<-a.datums(sc.689, depths=sc.689$absolute.ages, increasing.down = TRUE)


## Make splits for plot, normalize the range chart for plotting
splits<-c(3,6,9,12)
sc.plot <- breaks.fn (sc.689, splits)

## define colors
cols.vec <- rev(viridis(5))
cols.vec[1] <- 'gray70'

## define legend input
splits.legend<-paste(c(0,splits), '-', c(splits,100), '%', sep='')
splits.legend[1]<-paste('<', splits[1], '%', sep='')
splits.legend[length(splits.legend)] <- paste('>', splits[length(splits)], '%', sep='')

## Make the actual figure
pdf(file='figures/RangeChart689.pdf', height = 7, width = 10, useDingbats = FALSE)
par(mar=c(12, 4, 2, 4))
rangechart3(sc.plot, reorder='lad.by.fad', depths=sc.plot$absolute.ages, cex.xaxis=0.5, 
            cex.yaxis=1, cex.points="by.count", llwd=1, llcol = 'lightgray', llty = 3,
            col.points="by.count", colors.vec = cols.vec, xaxis.labels = 'names', 
            #xaxis.labels = c('numbers', 'names', 'alphanum')
            pch.points=16, baselines=FALSE, large.size=1, count.group=FALSE, legend.loc = 'bottomright', 
            legend.values = splits.legend, legend.horiz = TRUE, return.xaxis = FALSE) 
abline(h=33.9, col='gray')
dev.off()

