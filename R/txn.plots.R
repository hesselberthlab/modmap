# exp.plots.R 
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 552 $'
#
# modmap pipeline for plotting signals relative to
# transcribed regions 

library(ggplot2)
library(Cairo)
library(plyr)

# get the filename
output = commandArgs(trailingOnly=TRUE)
if (length(output) != 3) {
   stop("usage: Rscript txn.plots.R infile sample.name output.dir")
}

infile = output[1]
sample.name = output[2]
output.dir = output[3]

COLNAMES <- c('region.name', 'region.score', 'region.strand',
              'region.type', 'signal.strand', 'operation',
              'signal', 'signal.type')

df <- read.table(infile, col.names=COLNAMES)
if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

txn.box.plot <- function(df, sample.name, ... ) {
    
    # add label columns
    compartment <- sapply(as.character(df$region.type), 
                          function (x) strsplit(x, '-')[[1]][1])
    df$compartment <- compartment
    
    region <- sapply(as.character(df$region.type), 
                     function (x) strsplit(x, '-')[[1]][2])
    df$region <- region
    
    stats <- ddply(df, operation ~ compartment,
                   function (x) {
                   t <- t.test(signal ~ region, data = x)
                   with(t, data.frame(statistic, p.value))})
    
    gp <- ggplot(data = df,
                 aes(factor(region),
                     log10(signal)))

    gp <- gp + geom_boxplot(aes(fill = factor(region)))
    gp <- gp + scale_fill_brewer(palette="Set1", guide=FALSE)
    
    median.val <- log10(mean(df$signal[df$signal > 0]))
    labeler <- function(x) {
      return(c(y = median.val, label = length(x)))
    }
    
    gp <- gp + stat_summary(fun.data = labeler, geom = "text")
                            
    gp <- gp + facet_grid(operation ~ compartment)
    gp <- gp + theme_bw()

    gp <- gp + geom_text(data = stats,
                         aes(label = paste("t.test stat=", signif(statistic,5), "\n", 
                                           "p=", signif(p.value,5), "\n", sep=''),
                         x=-Inf, y=Inf, hjust=0, vjust=1))
    # axis labels 
    gp <- gp + ylab('log10(rNMP / bp)')
    gp <- gp + xlab('')

    # add title
    title <- paste('modmap transcription-analysis \n',
                      'sample = ', sample.name, sep='')
    gp <- gp + ggtitle(title)

    return(gp)
}

gp.txn.plot <- txn.box.plot(df, sample.name)

# write the files
txn.pdf.filename <- paste(output.dir, '/', 'modmap.transcription',
                      '.', sample.name, '.pdf', sep='')
ggsave(filename = txn.pdf.filename, plot = gp.txn.plot, 
       height = 8.5, width = 11, device = CairoPDF)
