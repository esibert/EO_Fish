#################################################
#                                               #
#   Code for running Mark/Recapture analyses    #
#               5/8/2019                        #
#    Running all on 596 for now, will make      #
#         new code for 689 when needed          #
#                                               #
#################################################


##### 1. Setup and select datasets #####
## Set MarkPath
MarkPath='C:/Software/MARK'   #Tell it where MARK is located!

## go into the MARK working directory
setwd('mark/') #need to put the mark files in a separate file, otherwise everything is a mess
setwd('..') # resets the working directory to the main folder for later

## Select dataset
morphdat.596 <- morphdat.combined.596
# morphdat.596 <- morphdat.all.596

## make SC objects necessary for plotting
sc<-build.strat.obj(morphdat.596)

##### 2. Set up MARK files #####
# sample size per time bin
samplesize<-rev(as.vector(apply(sc$counts, 1, sum))) #number of teeth considered in each time bin, oldest to youngest for MARK to work

#time intervals between time bins and make new ages vector
dt.rev <- c()
for (i in 1:length(sc$absolute.ages)-1) {
    ages <- rev(sc$absolute.ages)
    dt <- ages[i] - ages[i+1]
    dt.rev <- c(dt.rev, dt)
}

#set up the sample set for making .inp file
range.counts<-sc$counts #order by default is young to old
range.counts.rev<-range.counts[length(range.counts[,1]):1,] #reverse it so that order is old to young

# make .inp file
make.inp(x=range.counts.rev, filename='ranges_596.inp', header='ages old to young, full dataset, 1- age points')
# call in .inp file
rangesrev<-convert.inp('ranges_596')
# cleanup
rm(range.counts, range.counts.rev)

##### 3. Run Pradel models #####
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



##### 4. plot MARK output #####
# calculate average ages for plotting values between time points
ages.int <- c()
for (i in 1:(length(ages)-1)) {
    avg <- mean(c(ages[i], ages[i+1]))
    ages.int <- c(ages.int, avg)
}

### Plot model average values with "dropped estimates" 
errbar(x=ages.int, y=pr.ext$estimates$estimate, yplus = pr.ext$estimates$ucl, yminus = pr.ext$estimates$lcl, xlim = c(42, 28), ylim = c(0, 0.1), 
       type = 'o', col = 'red', errbar.col = 'red', pch = 16, ylab = 'Estimate', xlab = 'Age (Ma)')
errbar(x=ages.int, y=pr.orig$estimates$estimate, yplus = pr.orig$estimates$ucl, yminus = pr.orig$estimates$lcl, xlim = c(42,28),
       type = 'o', col = 'blue', errbar.col = 'blue', pch = 16, ylab = '', add = "TRUE")
abline (v = 33.9, col = 'gray')
legend('topleft', legend = c('Origination', 'Extinction'), col = c('blue', 'red'), pch = 16, lty = 1)
## Add probability of detection
# par(new = TRUE)
# errbar(x = ages, y = pr.p$estimates$estimate, yplus = pr.p$estimates$ucl, yminus = pr.p$estimates$lcl, xlim = c(42, 28),
#    type = 'o', col = 'black', errbar.col = 'black', pch = 16, ylab = '', axes = FALSE)
# axis(4)


### plot without error bars #

plot(ages.int, pr.ext$estimates$estimate, col = 'red', type = 'o', pch = 16)
points(ages.int, pr.orig$estimate$estimate, col = 'blue', type = 'o', pch = 16)

### Plot based on model average without using MARK to pull out indivual parameters (pull them out manually)

pradrec.avg <- model.average(Pradrec, vcv = TRUE)

#extract extinction = 1-survival
pr.ext.avg <- pradrec.avg$estimates[1:9,]
pr.ext.avg$estimate <- 1-pr.ext.avg$estimate
pr.ext.avg$lcl <- 1-pr.ext.avg$lcl
pr.ext.avg$ucl <- 1-pr.ext.avg$ucl
# extract probability of detection
pr.p.avg <- pradrec.avg$estimates[10:19,]
# extract probability of origination
pr.orig.avg <- pradrec.avg$estimates[20:28,]

# plot
errbar(x=ages.int, y=pr.ext.avg$estimate, yplus = pr.ext.avg$lcl, yminus = pr.ext.avg$ucl, xlim = c(42, 28), ylim = c(0, 1), 
       type = 'o', col = 'red', errbar.col = 'red', pch = 16, ylab = 'Estimate', xlab = 'Age (Ma)')
errbar(x=ages.int, y=pr.orig$estimates$estimate, yplus = pr.orig$estimates$ucl, yminus = pr.orig$estimates$lcl, xlim = c(42,28),
       type = 'o', col = 'blue', errbar.col = 'blue', pch = 16, ylab = '', add = "TRUE")
abline (v = 33.9, col = 'gray')
legend('topleft', legend = c('Origination', 'Extinction'), col = c('blue', 'red'), pch = 16, lty = 1)




### Plot only the best-fit model
model.to.plot <- Pradrec$Phi.dot.p.samplesize.f.dot

pr.ext.bestfit <- model.to.plot$results$real[1,]
pr.ext.bestfit$estimate <- 1-pr.ext.bestfit$estimate
pr.ext.bestfit$ucl <- 1- pr.ext.bestfit$ucl
pr.ext.bestfit$lcl <- 1 - pr.ext.bestfit$lcl
pr.p.bestfit <- model.to.plot$results$real[2:10,]
pr.orig.bestfit <- model.to.plot$results$real[11,]

# plot
errbar(x=ages.int, y=rep(pr.ext.bestfit$estimate, length(ages.int)), yplus = pr.ext.bestfit$lcl, yminus = pr.ext.bestfit$ucl, xlim = c(42, 28), ylim = c(0, 0.1), 
       type = 'o', col = 'red', errbar.col = 'red', pch = 16, ylab = 'Estimate', xlab = 'Age (Ma)')
errbar(x=ages.int, y=rep(pr.orig.bestfit$estimate, length(ages.int)), yplus = pr.orig.bestfit$ucl, yminus = pr.orig.bestfit$lcl, xlim = c(42,28),
       type = 'o', col = 'blue', errbar.col = 'blue', pch = 16, ylab = '', add = "TRUE")
abline (v = 33.9, col = 'gray')
legend('topleft', legend = c('Origination', 'Extinction'), col = c('blue', 'red'), pch = 16, lty = 1)
# add probability of detection
par(new = TRUE)
errbar(x = ages.int, y = pr.p.bestfit$estimate, yplus = pr.p.bestfit$ucl, yminus = pr.p.bestfit$lcl, xlim = c(42, 28),
    type = 'o', col = 'black', errbar.col = 'black', pch = 16, ylab = '', axes = FALSE)
axis(4)
