#########################################################
#                                                       #
#           functions.R for EO-fish project             #
#        Contains the functions used in the analyses    #
#   Created: 5/6/2019                                   #
#   Last updated: 6/5/2019                              #
#                                                       #
#########################################################


##### plot.ba.si function #####
plot.ba.si<-function(all_data, dataset, xlim.ages,
                     si.color='lightseagreen', ba.color='darkred', fish.color='black', 
                     si.pch=0, ba.pch=0, fish.pch=16,
                     si.lty=1, ba.lty=1, fish.lty=1,
                     si.ylim=NULL, iar.ylim=NULL,
                     plottitle=NULL, plottitle.loc=NULL, si.units=NULL, ba.units=NULL,
                     age.axis=FALSE, plot.margins=c(1,5,0,8), ...) { #
    
    if(dataset=='ODP748') {
        dat.site<-subset(all_data, Site=='ODP748' | Site=='ODP744')
    }
    else {
        dat.site<-subset(all_data, Site==dataset)
    }
    
    par(mar=plot.margins)
    ## Plot IAR
    
    #define dataset
    dat.fish<-subset(dat.site, Proxy=='IAR', select=c('Age', 'Value'))
    dat.fish <- dat.fish[order(dat.fish$Age),]
    
    #save dat.fish for later re-plotting
    dat<-dat.fish
    
    #if no value for the ylim for the plot in the parameters, make one
    if(is.null(iar.ylim)) {
        iar.ylim=c(0,max(dat$Value))
    }
    
    #make the plot
    plot(dat$Age, dat$Value, type='n', xlab='', ylab='', 
         xlim=xlim.ages, ylim=iar.ylim, bty='n', axes=F)
    axis(2, col=fish.color)
    # 5 point moving average
    lines(dat$Age, filter(dat$Value, ma5), lwd=1.5, col=fish.color)
    
    #ylabel
    mtext('IAR', side=2, line=3.2, cex=0.7)
    mtext(expression(paste('ich/cm'^2,'/kyr')), side=2, line=1.8, cex=0.5)
    
    #add annotations
    if(is.null(plottitle) == FALSE) {
        text(plottitle.loc[1], plottitle.loc[2], labels = plottitle, font=2)
    }
    if(age.axis==TRUE) {
        axis(side=1, at=seq(xlim.ages[1],xlim.ages[2], by=2))
        mtext('Age (Ma)', side=1, line=2, font=2, cex=1)
    }
    
    ## Add silica dataset
    if ('Si' %in% dat.site$Proxy) {
        dat<-subset(dat.site, Proxy=='Si', select=c('Age', 'Value'))
        dat <- dat[order(dat$Age),]
        if(is.null(si.ylim)==T) {
            par(new=T)
            plot(dat, pch=si.pch, lty=si.lty, type='o', col=si.color, xlab='', ylab='', 
                 xlim=xlim.ages, ylim=c(0, max(dat$Value)), bty='n', axes=F, lwd=2, cex=0.5)
            axis(4, line=-2, col=si.color, col.axis=si.color)
            mtext(si.units, side=4, line=0.5, cex=0.7, col=si.color)
        }
        if(is.null(si.ylim)!=T){
            par(new=T)
            plot(dat, pch=si.pch, lty=si.lty, type='o', col=si.color, xlab='', ylab='', 
                 xlim=xlim.ages, ylim=si.ylim, bty='n', axes=F, lwd=2, cex=0.5)
            axis(4, line=-2, col=si.color, col.axis=si.color)
            mtext(si.units, side=4, line=0.5, cex=0.7, col=si.color)
        }
    }
    
    ## add Barium dataset
    if ('Ba' %in% dat.site$Proxy) {
        dat<-subset(dat.site, Proxy=='Ba', select=c('Age', 'Value'))
        dat <- dat[order(dat$Age),]
        par(new=T)
        plot(dat, pch=ba.pch, lty=ba.lty, type='o', col=ba.color, xlab='', ylab='', 
             xlim=xlim.ages, ylim=c(0, max(dat$Value)), bty='n', axes=F, lwd=2, cex=0.5)
        axis(4, line=2.5, col=ba.color, col.axis=ba.color)
        mtext(ba.units, side=4, line=5, cex=0.7, col=ba.color)
    }
    
    ## re-add fish dataset on top of everything: 
    dat<-dat.fish
    par(new=T) 
    plot(dat$Age, dat$Value, pch=fish.pch, col=fish.color, xlab='', ylab='', 
         xlim=xlim.ages, ylim=iar.ylim, bty='n', axes=F)
    # 5 point moving average
    lines(dat$Age, filter(dat$Value, ma5), lwd=2, lty=fish.lty)
    
}


##### toothdat.cleanup function #####
toothdat.cleanup <- function(toothdat, sortby='age-obj') {
    teeth_all<-subset(toothdat, A==1)  #Clean up data to be just the teeth (TraitA = 1)
    teeth<-teeth_all[complete.cases(teeth_all),]  #clean up data to include just teeth that have been described
    teeth<-subset(teeth, B!=4) #Get rid of poor quality teeth if anything is still there... 
    teeth.dat <- teeth
    #sort the samples... 
    if(sortby=='age-obj') {
        teeth.dat<-teeth.dat[order(teeth.dat[,4]), ]  #age/object order (!) 
    }
    else if(sortby == 'age') {
        teeth.dat<-teeth.dat[order(teeth.dat[,2]), ] } #age-only order, objects may be scrambled
    else if(sortby == 'morph') {
        teeth.dat<-teeth.dat[order(teeth.dat[,3]), ] } #sort by morphotype number
    else if(sortby == 'original') {
        teeth.dat <- teeth.dat[order(as.numeric(row.names(teeth.dat))),] }  #put them in original (csv) order again
    else teeth.dat<-teeth.dat  #no ordering
    
    return(teeth.dat)
}


##### combine.samples(morphdat, combines) #####
# inputs: morphdat; combines - list of vectors defining samples in each bin 

combine.samples<-function(morphdat, combines) {
    
    #substitutions
    for(i in 1:length(combines)) {
        c.sub<-combines[[i]]
        morphdat$AgeID<-ifelse(morphdat$AgeID %in% c.sub, mean(c.sub), morphdat$AgeID)
    }
    
    return(morphdat)
}



##### build.strat.obj function #####
build.strat.obj<-function(morphdat) {
    #define unique age bins for the samples
    AgeID.unique<-unique(morphdat$AgeID)
    #call each unique taxa name
    taxa<-levels(morphdat$Morphotype_name)[1:length(levels(morphdat$Morphotype_Name))] 
    #build the strat column object for use with package stratigraph, and as storage for ranges
    sub<-data.frame(morphdat$AgeID, morphdat$Morphotype_Name)
    ranges<-table(sub, exclude='')  #make occurrance table#
    range.sc<-strat.column(ranges, absolute.ages=AgeID.unique, taxa=taxa)  #make "strat column" object
    #return the strat column object for manipulation
    return(range.sc)
}


##### normalize.sc #####
normalize.sc<-function(sc) {  #returns a normalized "counts" matrix. All rows (ages) should add up to 1. 
    norm.row<-function(row) {
        row/sum(row)
    }
    norm.sc<-apply(sc$counts, 1, norm.row)  #normalize the matrix by rows
    norm.sc<-t(norm.sc)                    #for some reason I have to transpose the output back to the normal "counts" form
    norm.sc
}


##### breaks.fn #####
# note that breaks.fn() uses normalize.sc() in it, so by running breaks.fn(), 
# it will automatically normalize the dataset. If it is already normalized, this 
# is fine - as everything will just be divided by 1, and therefore unchanged in value

breaks.fn <- function(sc, splits) {
    norm.sc<-100*normalize.sc(sc)
    splits<-c(0, splits, 100) #make breaks have full length
    for(i in 1:length(splits)-1) {
        norm.sc[norm.sc > splits[i] & norm.sc <= splits[i+1] ] = i
    }
    sc$counts <- norm.sc 
    return(sc)
}


##### rangechart3 #####
rangechart3 <- function (x = NULL, counts = NULL, depths = NULL, sample.labels = NULL, 
                         taxa = NULL, short.names = NULL, higher.grp = NULL, tax.cat = NULL, 
                         reorder = NULL, plot.points = FALSE, plot.depths.increasing.down = TRUE, 
                         llwd = 2, llcol = gray(0.5), llty = 1, cex.xaxis = 0.8, cex.yaxis = 1, 
                         cex.points = 1, pch.points = 1, y.axis.ticks = FALSE, 
                         col.points = 'black', colors.vec = NULL, baselines = FALSE, 
                         xaxis.labels = c('names', 'numbers', 'alphanum'), return.xaxis = FALSE, 
                         legend = TRUE, legend.values = NULL, legend.loc = NULL, legend.horiz = FALSE, 
                         legend.title = NULL, legend.bg = 'white', large.size = 1, count.group = FALSE, ...) 
{
    if (is.strat.column(x)) {
        counts <- x$counts
        if (is.null(depths)) 
            depths <- x$depths
        if (is.null(taxa)) 
            taxa <- x$taxa
        if (is.null(short.names)) 
            short.names <- x$short.names
        if (is.null(higher.grp)) 
            higher.grp <- x$higher.grp
        if (is.null(tax.cat)) 
            tax.cat <- x$tax.cat
    }
    if (is.data.frame(counts) || is.matrix(counts)) {
        counts <- apply(counts, 2, "as.numeric")
    }
    else if (is.character(counts) && length(counts) == 1) {
        counts <- read.csv(counts, header = TRUE, skip = 0, colClasses = "")
        counts <- apply(counts, 2, "as.numeric")
    }
    else {
        stop("argument to counts not understood")
    }
    if (is.null(depths)) {
        warning("no depths provided; plotting samples at regular intervals")
        depths <- 1:nrow(counts)
    }
    depths <- as.numeric(depths)
    if (length(depths) != nrow(counts)) {
        stop(paste(length(depths), " depths, and ", nrow(counts), 
                   " rows in the count matrix.", sep = ""))
    }
    if (sum(is.na(counts)) > 0) {
        warning(paste(sum(is.na(counts)), "missing values in count matrix replaced with zeros"))
        counts[is.na(counts)] <- 0
    }
    emptycols <- !(colSums(counts, na.rm = TRUE) > 0)
    if (any(emptycols)) {
        counts <- counts[, !emptycols]
        tax.cat <- tax.cat[!emptycols]
        taxa <- taxa[!emptycols]
        warning(paste(sum(emptycols), "columns with zero counts at all levels removed"))
    }
    
    ######alphanumeric rewrite
    if(xaxis.labels == 'alphanum') {
        colnames(counts) <- 1:length(colnames(counts))
    }
    
    if (!is.null(reorder)) {
        if (length(reorder) == 1 && is.character(reorder)) {
            funny <- function(x) return((1:length(x))[x > 0])
            if (pmatch(reorder, "fad.by.category", nomatch = FALSE)) {
                fads <- depths[as.numeric(lapply(apply(counts, 
                                                       2, funny), max))]
                reorder.vect <- sort(fads, decreasing = TRUE, 
                                     index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
                reorder.vect <- sort(as.character(taxa), index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
            }
            else if (pmatch(reorder, "lad.by.category", nomatch = FALSE)) {
                lads <- depths[as.numeric(lapply(apply(counts, 
                                                       2, funny), min))]
                reorder.vect <- sort(lads, decreasing = TRUE, 
                                     index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
                reorder.vect <- sort(as.character(taxa), index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
            }
            else if (pmatch(reorder, "lad.by.fad", nomatch = FALSE)) {
                lads <- depths[as.numeric(lapply(apply(counts, 
                                                       2, funny), min))]
                reorder.vect <- sort(lads, decreasing = TRUE, 
                                     index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
                fads <- depths[as.numeric(lapply(apply(counts, 
                                                       2, funny), max))]
                reorder.vect <- sort(fads, decreasing = TRUE, 
                                     index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
                reorder.vect <- sort(as.character(taxa), index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
            }
            else if (pmatch(reorder, "by.count", nomatch = FALSE)) {
                reorder.vect <- sort(colSums(counts), decreasing = TRUE, 
                                     index.return = TRUE)$ix
                counts <- counts[, reorder.vect]
            }
        }
        else if (length(reorder) == ncol(counts)) {
            reorder.vect <- as.numeric(reorder)
            counts <- counts[, reorder.vect]
        }
        else {
            stop("argument to reorder not understood")
        }
    }
    if (!is.null(tax.cat)) {
        if (length(tax.cat) != ncol(counts)) {
            warning("taxon category labels seem to be the wrong length")
            tax.cat <- NULL
        }
    }
    if (is.null(tax.cat)) {
        tax.cat <- rep("", ncol(counts))
    }
    
    if (xaxis.labels == 'names') { 
        xaxis.labels <- colnames(counts)
    }
    if (xaxis.labels == 'numbers') {
        colnum <- dim(counts)[2]
        xaxis.labels <- as.character(c(1:colnum))
    }
    
    if(xaxis.labels == 'alphanum') {
        xaxis.labels <- colnames(counts)
    }
    
    ad <- a.datums(strat.column(counts = counts, depths = depths))
    plot(1:ncol(counts), ylim = c(max(depths), min(depths)), 
         type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "Age (Ma)")
    segments(1:ncol(counts), ad[, 1], 1:ncol(counts), ad[, 2], 
             lwd = llwd, col = llcol, lty = llty)
    if (baselines == TRUE) { segments(1:ncol(counts), ad[, 2], 1:ncol(counts), rep(par()$usr[3], 
                                                                                   ncol(counts)), col = 'lightblue', lty = 3, lwd = 0.5) } #col=grey(0.5)
    axis(1, at = 1:ncol(counts), labels = xaxis.labels, cex.axis = cex.xaxis, ##### labels = colnames(counts)
         las = 3)
    axis(2, las = 1, cex.axis = cex.yaxis)
    if(y.axis.ticks == TRUE) {axis(2, at = depths, labels = FALSE, tck = -0.01)} 
    
    #added code: 
    if (is.null(colors.vec)) colors<-rainbow(max(counts), end=5/6) else colors<-colors.vec
    sizes<-seq(0,(max(counts)-1), 1)
    sizes<-((sizes/(max(sizes)/large.size)) + large.size)
    #groups
    group.fn<-function(x) {
        for (i in 1:length(x)) {
            if(x[i]>=2 & x[i]<=5) x[i]<-4
            if(x[i]>=6 & x[i]<=10) x[i]<-7
            if(x[i]>10) x[i] <- 15
        }
        return(x)
    }
    if (count.group == TRUE) counts.plt<-group.fn(counts) else counts.plt<-counts
    
    # make the plots! 
    for (i in 1:ncol(counts.plt)) {
        num.val<-c(counts.plt[,i][counts.plt[,i]>0])
        if (pch.points == 'by.count') points<-num.val else points<-pch.points
        if (col.points == 'by.count') cols<-c(colors[num.val]) else cols<-col.points
        if (cex.points == 'by.count') points.cex<-c(sizes[num.val]) else points.cex<-cex.points
        plocs <- depths[(counts.plt > 0)[, i]]
        points(rep(i, length(plocs)), plocs, cex=points.cex, pch=points, col=cols,
               ...)
    }
    
    # add a legend
    if(is.null(legend.values)) {legend.values <- seq(1:max(counts.plt))}
    if(is.null(legend.loc)) {legend.loc='topleft'} 
    if (legend == TRUE & count.group == FALSE) { 
        if (pch.points == 'by.count') leg.pch<-c(1:max(counts.plt)) else leg.pch<-pch.points
        if (col.points == 'by.count') leg.col<-colors else leg.col<-col.points
        if (cex.points == 'by.count') leg.cex<-sizes else leg.cex<-cex.points
        legend(legend.loc, legend=legend.values, pch=leg.pch, 
               col=leg.col, pt.cex=leg.cex, cex=large.size, title=legend.title,
               horiz = legend.horiz, bg = legend.bg) }
    if (legend == TRUE & count.group == TRUE) {
        if (pch.points == 'by.count') leg.pch<-sort(unique(as.vector(counts.plt)))[2:5] else leg.pch<-pch.points
        if (col.points == 'by.count') leg.col<-colors[sort(unique(as.vector(counts.plt)))[2:5]] else leg.col<-col.points
        if (cex.points == 'by.count') leg.cex<-sizes[sort(unique(as.vector(counts.plt)))[2:5]] else leg.cex<-cex.points
        legend(legend.loc, legend=c('1', '2-5', '6-10', '10+'), pch=leg.pch, 
               col=leg.col, pt.cex=leg.cex, cex=large.size, title=legend.title,
               horiz = legend.horiz, bg = legend.bg) }
    
    if(return.xaxis == TRUE) {
        if(xaxis.labels == 'alphanum') {
            return(taxa)}
        else {
            return(colnames(counts))
        }
    }
}



##### make.inp(x, filename, header=NULL) #####
make.inp<-function(x, filename, header=NULL) {  #make input file for RMark, x is table sc$counts 
    mat<-as.vector(x)
    mat<-matrix(mat, nrow=length(x[1,]), ncol=length(x[,1]), byrow=T)
    rownames(mat)<-colnames(x)
    mat<-ifelse(mat>0, 1, 0)  #change to 1's and 0's
    if(is.null(header)) header=filename
    filefoo<-file(paste(filename, sep=''), 'w')
    writeLines(paste('/*', header, '*/', sep=''), filefoo)
    for(i in 1:length(mat[,1])) {
        mm<-mat[i,]
        writeLines(paste('/* ', rownames(mat)[i], ' */ ', paste(mm, sep='', collapse=''), ' 1;', sep=''), filefoo) }
    close(filefoo)
}

# ##### rangechart #####
# 
# # Range chart runciton rewritten to *not* require stratigraph internal functions. This will be added to the ichthyoliths package at a later date. 
# 
# 
# # 0. counts is a table of occurrances of each taxa (columns) by age (rows)
# # 1. Grab ages and taxa from counts table as needed
# # 2. reorder options: (will include categories at some point but not right now)
# #       a) FAD (first occurrance datum) then LAD (fad.by.lad)
# #       b) LAD (last occurrance datum) then FAD (lad.by.fad)
# # 3. xaxis.labels = c('numbers', 'names', 'alphanum') and return.xaxis = c(TRUE, FALSE)
# #       a) 'numeric' means that they will be numbered 1-n on the x-axis. 
# #       b) 'names' means that the names of the types will be the x-axis labels
# #       c) 'alphanum' means that the original taxa order is preserved, and the numbers on the x-axis refer to the original order of taxa in the counts matrix. 
# #       NOTE: return.xaxis = TRUE will give a list of corresponding morphotype names for each x-axis numeric value. 
# # 4. normalize.counts (true/false) - uses internal function to normalize the counts 
# # 5. count.breaks = NULL; Input vector defining the 'breaks' for plotting numbers or normalized percentages of assemblages. 'default' value is count.breaks = c(0, 2, 5, 10, 100)
# #       a) if normalize.counts is TRUE and count.breaks = 'default', the breaks are at percentage values 1-2%, 2-5%, 5-10%, and >10% (count.breaks = c(0, 2, 5, 10, 100))
# #       b) if normalize.counts is FALSE and , the breaks are at 1-2, 2-5, 5-10, and 10+. (count.breaks = c(0,2,5,10,100))
# #       NOTE: can take inputs of any length, doesn't have to be 4 bins. Percentages and counts are both given in whole integers; 
# # 6. LINE Graphical Parameters notes: 
# #       Line segments connecting points are llwd (width), llcol (color), and llty (line type)
# #       Baseline (connecting x-axis with first point) are blwd (width), blcol (color), and blty (line type)
# # 7. POINT Graphical Parameters notes
# #       cols.vec is a vector of colors given for each point value or bin. Defaults to rainbow if not given
# #       pch.points, col.points, can be single values or 'by.count' or 'by.category'. 
# #           col.vec and pch.vec are inputs for the 'by.count' and 'by.category' options. 
# #       cex.points can be single values
# #       largesize is a scaling factor for making the points/plot larger.
# # 8. tax.cat is a vector matching the taxa (in order of the counts table) to a category. Must be numeric.
# # 9. LEGEND information: make your own legend after, placing it where you'd like, with the colors you'd like, etc... 
# #       
# 
# 
# # Note that if you just use the ages and taxa generated from the counts table, there's no worrying about reordering or mismatched ages. Just use the ages and taxa from the counts table...
# 
# rangechart <- function(counts, ages = NULL, taxa = NULL, tax.cat = NULL, reorder = NULL, 
#                         normalize.counts = FALSE, count.breaks = c(0, 2, 5, 10, 100), 
#                         cex.xaxis = 1, cex.yaxis = 1, yaxis.ticks = FALSE, 
#                         llwd = 1, llcol = 'gray70', llty = 3, 
#                         baselines = FALSE, blwd = 0.5, blcol = 'lightblue', blty = 3, 
#                         cols.vec = NULL, col.points = 'by.count', 
#                         pch.points = 16, pch.vec = NULL, 
#                         cex.points = 1, largesize = 1, 
#                         xaxis.labels = c('names', 'numeric', 'alphanum'), print.xaxis = FALSE, ...) {
#     
#     ##### set up the dataset #####
#     # Ages should be rownames of the counts table, and in increasing order
#     if(missing(ages)) {
#         ages <- as.numeric(rownames(counts))
#     }
#     else {
#         rownames(counts) <- ages
#     }
#     
#     # if the ages are not in increasing order, sort them and the counts table to be so
#     if(is.unsorted(ages) == TRUE) {
#         age.increasing <- sort(ages, index.return = TRUE)$ix
#         counts <- counts[age.increasing, ]
#         ages <- as.numeric(rownames(counts))
#     }
#     
#     # taxa should be column-names of the counts table, order doesn't matter. 
#     if(missing(taxa)) {
#         taxa <- as.character(colnames(counts))
#     }
#     else {
#         colnames(counts) <- taxa
#     }
#     
#     original.taxa <- taxa #useful for matching tax-cat later too. 
#     
#     
#     # clear NA values, if any, by replacing with zeros
#     if (sum(is.na(counts)) > 0) {
#         warning(paste(sum(is.na(counts)), "missing values in count matrix replaced with zeros"))
#         counts[is.na(counts)] <- 0
#     }
#     
#     # FAD: First (oldest) occurance datum calls the maximum index (mapped to the ages values) 
#     # of a non-zero count value for each taxa column in the counts matrix
#     fad <- ages[apply(counts, 2, function(x) {max(which (x!=0))})] 
#     
#     # LAD: Last (youngest) occurance datum calls the minimum index (mapped to the ages values) 
#     # of a non-zero count value for each taxa column in the counts matrix
#     lad <- ages[apply(counts, 2, function(x) {min(which (x!=0))})] 
#     
#     
#     ### Normalize the counts if they want to be normalized 
#     if(normalize.counts == TRUE) {
#         norm.row<-function(row) {
#             row/sum(row)
#         }
#         norm.counts<-apply(counts, 1, norm.row)  #normalize the matrix by rows
#         counts<-t(norm.counts)   #for some reason I have to transpose the output back to the normal "counts" form
#         counts <- 100 * counts #make this a percentage instead of a decimal value. 
#     }
#     
#     ### Group the counts into bins if you'd like them to be binned. 
#     if(!is.null(count.breaks)) {
#         count.breaks <- count.breaks
#     }
#     
#     for(i in 1:length(count.breaks)-1) {
#         counts[counts > count.breaks[i] & counts <= count.breaks[i+1] ] = i
#     }
#     
#     ##### reorder counts #####
#     if(!is.null(reorder)) {
#         
#         # fad.by.lad (origination)
#         if(reorder == 'fad.by.lad') {
#             # First reorder the counts by LAD and recalculate FAD
#             reorder.vect <- sort(lad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#             fad <- ages[apply(counts, 2, function(x) {max(which (x!=0))})] 
#             
#             # Next, reorder the counts by FAD
#             reorder.vect <- sort(fad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#             
#         }
#         
#         # lad.by.fad (extinction)
#         else if(reorder == 'lad.by.fad') {
#             # First reorder the counts by FAD and recalculate LAD
#             reorder.vect <- sort(fad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#             lad <- ages[apply(counts, 2, function(x) {min(which (x!=0))})] 
#             
#             # Next, reorder the counts by LAD
#             reorder.vect <- sort(lad, decreasing = TRUE, index.return = TRUE)$ix #pulls index of order by fads
#             counts <- counts[, reorder.vect]
#         }
#     }
#     
#     
#     ### re-generate 'taxa', 'ages', 'fad' and 'lad' and 'tax.cat' from the updated counts table
#     
#     taxa <- as.character(colnames(counts))
#     ages <- as.numeric(rownames(counts))
#     fad <- ages[apply(counts, 2, function(x) {max(which (x!=0))})] 
#     lad <- ages[apply(counts, 2, function(x) {min(which (x!=0))})] 
#     
#     if(!is.null(tax.cat)) { tax.cat <- tax.cat[match(taxa, original.taxa)] }
#     
#     ##### set up the graphical parameters #####
#     ## xaxis.labels
#     if(missing(xaxis.labels)) { xaxis.labels <- 'names' }
#     if(xaxis.labels == 'names') { xaxis.lab <- taxa }
#     if(xaxis.labels == 'numeric') { xaxis.lab <- 1:length(taxa) }
#     if(xaxis.labels == 'alphanum') { xaxis.lab <- match(taxa, original.taxa) }
#     
#     ## colors of points
#     if(missing(cols.vec)) {cols.vec = 'gray'}
#     
#     if(cols.vec[1] == 'rainbow') {
#         colors<-rainbow(max(counts), end=5/6) 
#     }
#     else if (cols.vec[1] == 'viridis') {
#         colors <- viridis::viridis(max(counts))
#     }
#     else colors<-cols.vec
#     
#     ## sizes of points (and scale of the whole thing...)
#     sizes<-seq(0,(max(counts)-1), 1)
#     sizes<-((sizes/(max(sizes)/largesize)) + largesize)
#     
#     
#     ##### Actually make the plot #####
#     ### make blank plot with appropriate dimensions, suppress x- and y- axes
#     plot(1:ncol(counts), ylim = c(max(ages), min(ages)), 
#          type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "Age (Ma)")
#     
#     ### add segments to the plot
#     segments(1:ncol(counts), lad, 1:ncol(counts), fad, 
#              lwd = llwd, col = llcol, lty = llty)
#     
#     if (baselines == TRUE) { 
#         segments(1:ncol(counts), fad, 1:ncol(counts), rep(par()$usr[3], ncol(counts)), 
#                  col = blcol, lty = blty, lwd = blwd) 
#     } 
#     
#     ### Add points to the plot
#     for (i in 1:ncol(counts)) {
#         num.val<-c(counts[,i][counts[,i]>0])
#         
#         # point characters
#         if (pch.points == 'by.count') {
#             pts.pch<-num.val
#         }
#         else if (pch.points == 'by.category') {
#             pts.pch <- pch.vec[tax.cat[i]] 
#         }
#         else { pts.pch<-pch.points }
#         
#         # point colors
#         if (col.points == 'by.count') {
#             pts.cols<-c(colors[num.val]) 
#         }
#         else if (col.points == 'by.category') {
#             pts.cols <- colors[tax.cat[i]] 
#         }
#         else { pts.cols<-col.points }
#         
#         # point size
#         if (cex.points == 'by.count') pts.cex<-c(sizes[num.val]) else pts.cex<-cex.points
#         
#         # point y-values
#         plocs <- ages[(counts > 0)[, i]]
#         
#         #actually add the points
#         points(rep(i, length(plocs)), plocs, cex=pts.cex, pch=pts.pch, col=pts.cols,
#                ...)
#     }
#     
#     ### add axes
#     axis(1, at = 1:ncol(counts), cex.axis = cex.xaxis, labels = xaxis.lab, 
#          las = 3)
#     axis(2, las = 1, cex.axis = cex.yaxis)
#     if(yaxis.ticks == TRUE) {axis(2, at = ages, labels = FALSE, tck = -0.01)} 
#     
#     
#     ##### print taxa in list #####
#     # this has to be the last thing that the function does, because R stops after a return value
#     if(print.xaxis == TRUE) {
#         if(xaxis.labels == 'alphanum') {
#             return(original.taxa)
#         } 
#         else {
#             return(taxa)
#         }
#     }
#     
# }