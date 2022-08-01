

################################################################################################################
#################### This Code is modified from https://github.com/PoisonAlien/maftools/blob/master/R/tcgacompare.R plot
################### Modifications are done, becuase of the different input type..... 
#############################################################################################################

## Example of input data
#We downloaded the ID9/mb data from cosmic3   ### https://cancer.sanger.ac.uk/signatures/documents/568/v3.2_ID9_TISSUE_yYzLQXR.txt
#https://cancer.sanger.ac.uk/signatures/id/id9/

ID9 <- read.csv("mutation_load_comparison_plot/example.csv", sep = '\t')  ## Read data

colnames(ID9) <- c("Cancer.Types","ID9mb")
cancers <- unique(ID9$Cancer.Types)

library(dplyr)
library(data.table)
median <- as_tibble(ID9) %>% 
  group_by(Cancer.Types) %>% 
  summarise(median(ID9mb))

median <- median[match(cancers, median$Cancer.Types),]
colnames(median) <- c("Cancer.Types","median1")

####
tt <- as.data.frame(table(ID9[1]))  ##Number of samples for each cancer types

#library(plyr)
#tt <- count(ID9$Cancer.Types)
tt<- tt[match(cancers, tt$Cancer.Types),]

tcga.cohort.med = merge(tt,median, by.x = 'Cancer.Types', by.y = 'Cancer.Types')
tcga.cohort.med <- tcga.cohort.med[match(cancers, tcga.cohort.med$Cancer.Types),]

#Median mutations
colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')
 
###############################

tcga.cohort = split(ID9, as.factor(ID9$Cancer.Types))

## reorder the plot.dat list
tcga.cohort <-tcga.cohort[cancers]

plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
  x = tcga.cohort[[i]]
  x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                             x[order(x$ID9mb, decreasing = TRUE), 'ID9mb'  ],
                             x[,'Cancer.Types']
                             )
  #colnames(x) <- c("V1","Cancer.Types","ID9mb")
})
names(plot.dat) = names(tcga.cohort)
 
logscale = TRUE

if(logscale){
  y_lims = range(log10(data.table::rbindlist(l = plot.dat)[V2 != 0][,V2]))
}else{
  y_lims = range(data.table::rbindlist(l = plot.dat)[,V2])
}

#####################

col = c('gray70', 'black')
bg_col = c('#EDF8B1', '#2C7FB8')
medianCol = 'red'
cohortFontSize = 0.6
axisFontSize = 0.8
#par(mar = c(4, 3, 2.5, 1.5))
#y_lims <- c(-0.01,1.4)
y_min = floor(min(y_lims))
y_max = ceiling(max(y_lims))
y_lims = c(y_min, y_max)
y_at = pretty(y_lims)


plot(NA, NA, xlim = c(0,length(cancers)), ylim = y_lims , axes = FALSE, xlab = NA, ylab = NA, frame.plot = TRUE)
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = grDevices::adjustcolor(col = "gray", alpha.f = 0.1))
rect(xleft = seq(0, length(cancers)-1, 1), ybottom = min(y_lims), xright = seq(1,length(cancers), 1),
     ytop = max(y_lims), col = grDevices::adjustcolor(col = bg_col, alpha.f = 0.2),
     border = NA)
#rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4], col = grDevices::adjustcolor(col = bg_col, alpha.f = 0.2))

abline(h = pretty(y_lims), lty = 2, col = "gray70")
#abline(v = seq(1, length(plot.dat)), lty = 1, col = "gray70")

########################################################3
 
lapply(seq_len(length(plot.dat)), function(i){
  x = plot.dat[[i]]
  if(x[1, 'V3'] == "Bladder-TCC"){   ###Cancer type you want to highlight
    if(logscale){
      points(x$V1, log10(x$V2), pch = 16, cex = 0.4, col = col[2])
    }else{
      points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[2])
    }
  }else{
    if(logscale){
      points(x$V1, log10(x$V2), pch = 16, cex = 0.4, col = col[1])
      #print("Haan log haiga va")
    }else{
      points(x$V1, x$V2, pch = 16, cex = 0.4, col = col[1])
    }
  }
})

samp_sizes = lapply(plot.dat, nrow)
axis(side = 1, at = seq(0.5, length(plot.dat)-0.5, 1), labels = names(plot.dat),
     las = 2, tick = FALSE, line = -0.8, cex.axis = cohortFontSize)
axis(side = 3, at = seq(0.5, length(plot.dat)-0.5, 1), labels = unlist(samp_sizes),
     las = 2, tick = FALSE, line = -0.8, font = 3, cex.axis = cohortFontSize)

 
tcga.cohort.med$Median_Mutations_log10 <-  log10(tcga.cohort.med$Median_Mutations)

sidePos = 2 # 4 is bad
linePos = 2


if(logscale){
axis(side = 2, at = y_at, las = 2, line = -0.6, tick = FALSE, labels = round(10^(y_at),5), cex.axis = axisFontSize)
mtext(text = "ID9-mutations-per-mb", side = sidePos, line = linePos)
}else{

  axis(side = 2, at = y_at, las = 2, line = -0.6, tick = FALSE, cex.axis = axisFontSize)
  mtext(text = "ID9-mutations-per-mb", side = sidePos, line = linePos)
  
  
}

if(logscale){
  lapply(seq_len(nrow(tcga.cohort.med)), function(i){
    segments(x0 = i-1, x1 = i, y0 = tcga.cohort.med[i, 'Median_Mutations_log10'],
             y1 = tcga.cohort.med[i, 'Median_Mutations_log10'], col = medianCol)
  })
}else{
  lapply(seq_len(nrow(tcga.cohort.med)), function(i){
    segments(x0 = i-1, x1 = i, y0 = tcga.cohort.med[i, 'Median_Mutations'],
             y1 = tcga.cohort.med[i, 'Median_Mutations'], col = medianCol)
  })
}

