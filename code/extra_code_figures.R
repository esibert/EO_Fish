#################################################
#    Old figures and  code                      #
#   compiled 5/6/2019                           #
#################################################


##### non-log plot of all sites #####
pdf('figures/magnitudes_nonlog.pdf', height=7.5, width=7.5, useDingbats = F)

par(mfrow=c(1,1))
par(mar=c(5.1, 4.1, 2.1, 3.1))
plot(iar.all$Age, iar.all$Value, type='n', xlim=xlim.ages, xlab='Age (Ma)', ylab='')
#rect(36.8, 1, 37.3, 30000, col='gray85', border=NA)
#rect(33, 1, 34, 30000, col='azure', border=NA)
#rect(27, 1, 44, 800, col=adjustcolor('lightgray', alpha.f = 0.5), border=NA)
#rect(27, 800, 37.5, 30000, col=adjustcolor('lightgreen', alpha.f = 0.5), border=NA)
abline(h=800, lty=2, lwd=2)
mtext('IAR (log scale)', side=2, line=3, cex=1)
mtext(expression(paste('ich/cm'^2,'/kyr')), side=2.2, line=2, cex=0.8)

#add 5 pt moving averages (scaled 886)
for(i in 1:length(sites)) {
    lines(subset(iar.all, Site==sites[i])$Age, filter(subset(iar.all, Site==sites[i])$Value, ma5), 
          pch=pchvec[i], col=adjustcolor(fillvec[i], alpha.f=1), lwd=5)
}

#add points
for(i in 1:length(sites)) {
    points(subset(iar.all, Site==sites[i])$Age, subset(iar.all, Site==sites[i])$Value, 
           pch=pchvec[i], col=outvec[i], bg=fillvec[i], cex=1)
}
legend('topright',legend=sites,col=outvec, pt.bg=fillvec, pch=pchvec, bg='white', ncol=2)

dev.off()



##### t-test for DSDP Site 522 - is it higher pre-EO or post-EO? #####
d522<-subset(all_data, Site=='DSDP522')
le<-subset(d522, Age >=33.7, select='Value')
eo<-subset(d522, Age <=33.7 & Age >=31.5, select='Value')
t.test(le, eo)
