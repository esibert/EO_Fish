## Scale factor for converting ODP Site 886's counts of ichthyoliths >106 um to compare to the other sites, which are all counts >38 um. 

## For each site in the study, we calculated individual ratios of ich > 106um to ich > 38 um, and found that
## the scale factor for converting ich>106um to ich>38um is approximately 4.25. 
## in other words, ich >106um * 4.25 = ich >38um
## this was used to convert the IAR >106 um values from ODP 886 to Approx IAR >38um, 
## for comparison between ODP 886 and the other sites. 



# 596 % 106
xx<-readClipboard() #percentage of teeth >106um from DSDP 596
xx<-as.numeric(xx)
hist(xx)
mean(xx)
0.2499268
median(xx)
0.2297011
sd(xx)
0.08269078
#scaling factor
4.25

# 689 %106
yy<-readClipboard()
yy<-as.numeric(yy)
hist(yy)
mean(yy)
0.2230788
median(yy)
0.2222222
sd(yy)
0.1348631
#scaling factor
4.5

# 1217 %106
zz<-readClipboard()
zz<-as.numeric(zz)
hist(zz)
mean(zz)
0.263095
median(zz)
0.2666667
sd(zz)
0.1141313
#scaling factor
3.8

# 1406 %106
ww<-readClipboard()
ww<-as.numeric(ww)
hist(ww)
mean(ww)
0.2812705
median(ww)
0.2804487
sd(ww)
0.08461654
#scaling factor
3.57

#748
vv<-readClipboard()
vv<-as.numeric(vv)
hist(vv)
mean(vv)
0.1837183
median(vv)
0.2
sd(vv)
0.1418174
#scaling factor
5

#522 is not possible to calculate this factor with, as it was counted as >63 and <63, with no 106 designation


# all together: 
all106<-c(vv, ww, xx, yy, zz)
mean(c(vv, ww, xx, yy, zz))
0.2347096
median(c(vv, ww, xx, yy, zz))
0.25
sd(c(vv,ww,xx,yy,zz))
0.1250832
#Scaling factor (1/.235)
4.25


