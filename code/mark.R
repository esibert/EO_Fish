#################################################
#                                               #
#   Code for running Mark/Recapture analyses    #
#           Copy-pasted withou testing          #
#               5/8/2019                        #
#                                               #
#################################################




library(RMark)
MarkPath='C:/Software/MARK'   #Tell it where MARK is located!
setwd('mark/') #need to put the mark files in a separate file, otherwise everything is a mess
setwd('..') # resets the working directory to the main folder for later

# sample size per time bin
samplesize<-rev(as.vector(apply(sc$counts, 1, sum))) #number of teeth considered in each time bin, oldest to youngest for MARK to work

#time intervals between time bins
dt.rev <- c()
for (i in 1:length(sc$absolute.ages)-1) {
    ages <- rev(sc$absolute.ages)
    dt <- ages[i] - ages[i+1]
    dt.rev <- c(dt.rev, dt)
}

#set up the sample set for making .inp file
range.counts<-sc$counts #order by default is Eocene to Cretaceous
range.counts.rev<-range.counts[length(range.counts[,1]):1,] #reverse it so that order is Cretaceous to Eocene
# make .inp file
make.inp(x=range.counts.rev, filename='rangesrev_con.inp', header='ages old to young, full dataset, 1- age points')
# call in .inp file
rangesrev<-convert.inp('rangesrev_con')
# cleanup
rm(range.counts, range.counts.rev)


# Pradel recruitment
dp.pradrec<-process.data(rangesrev, model='Pradrec', time.intervals=dt.rev)  
ddl.pradrec<-make.design.data(dp.pradrec)
ddl.pradrec$p$samplesize<-samplesize #add sample size variable to test for variable p

#write wrapper function
Pradel.recruit<-function() {
    # setwd('C:/Elizabeth_Files/Personal/Research/R_Files/Manuscripts/ella-eo-fish/mark/') #need to put the mark files in here... 
    #Create formulas for Phi, p, and f
    Phi.time<-list(formula=~time)
    Phi.dot<-list(formula=~1)
    p.time<-list(formula=~time)
    p.dot<-list(formula=~1)
    #p.effort<-list(formula=~effort)
    p.samplesize <- list(formula =~ samplesize)
    #p.effortplussamplesize<-list(formula =~ effort + samplesize)
    f.time<-list(formula=~time)
    f.dot<-list(formula=~1)
    #create models to run by combining the above formulas, looking for objects with Phi. p. and f. at beginning
    cml<-create.model.list('Pradrec')
    #run all the models
    results<-mark.wrapper(cml,data=dp.pradrec,ddl=ddl.pradrec,output=FALSE,silent=TRUE)
    return(results)
}

# testing
# mark(dp.pradrec, ddl= ddl.pradrec, model = 'Pradrec', time.intervals = dt.rev)

# run MARK
Pradrec<-Pradel.recruit()    #run models
Pradrec  #displays table

pradrec.avg <- model.average(Pradrec, vcv = TRUE)

# Extinction estimates - inverse of Phi
pr.ext <- model.average(Pradrec, "Phi", vcv = TRUE)
pr.ext$estimates$estimate <- 1-pr.ext$estimates$estimate
pr.ext$estimates$lcl <- 1-pr.ext$estimates$lcl
pr.ext$estimates$ucl <- 1-pr.ext$estimates$ucl
names(pr.ext$estimates)[c(4, 5)] <- c('ucl', 'lcl')

# origination estimates - f
pr.orig <- model.average(Pradrec, "f", vcv = TRUE)

pr.p <- model.average(Pradrec, 'p', vcv = TRUE)



### plot MARK output
# calculate average ages for plotting values between time points
ages.int <- c()
for (i in 1:(length(ages)-1)) {
    avg <- mean(c(ages[i], ages[i+1]))
    ages.int <- c(ages.int, avg)
}

errbar(x=ages.int, y=pr.ext$estimates$estimate, yplus = pr.ext$estimates$ucl, yminus = pr.ext$estimates$lcl, xlim = c(42, 28), ylim = c(0, 0.1), 
       type = 'o', col = 'red', errbar.col = 'red', pch = 16, ylab = 'Estimate', xlab = 'Age (Ma)')
errbar(x=ages.int, y=pr.orig$estimates$estimate, yplus = pr.orig$estimates$ucl, yminus = pr.orig$estimates$lcl, xlim = c(42,28),
       type = 'o', col = 'blue', errbar.col = 'blue', pch = 16, ylab = '', add = "TRUE")
abline (v = 33.9, col = 'gray')
legend('topleft', legend = c('Origination', 'Extinction'), col = c('blue', 'red'), pch = 16, lty = 1)
# par(new = TRUE)
# errbar(x = ages, y = pr.p$estimates$estimate, yplus = pr.p$estimates$ucl, yminus = pr.p$estimates$lcl, xlim = c(28,42), 
#    type = 'o', col = 'black', errbar.col = 'black', pch = 16, ylab = '', axes = FALSE)
# axis(4)


# plot without error bars #

plot(ages.int, pr.ext$estimate, col = 'red', type = 'o', pch = 16)
points(ages.int, pr.orig$estimate, col = 'blue', type = 'o', pch = 16)
