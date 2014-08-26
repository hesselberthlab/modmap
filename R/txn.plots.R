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
    df$compartment <- as.factor(compartment)
    
    region <- sapply(as.character(df$region.type), 
                     function (x) strsplit(x, '-')[[1]][2])
    df$region <- as.factor(region)
   
    # calculate t.tests among groups
    stats <- ddply(df, operation ~ compartment + signal.strand,
                   function (x) {
                   t <- t.test(signal ~ region, data = x)
                   with(t, data.frame(statistic, p.value))})
    
    gp <- ggplot(data = df,
                 aes(factor(region),
                     log10(signal)))

    gp <- gp + geom_boxplot(aes(fill = factor(region)),
                            show_guide = FALSE)
    gp <- gp + scale_fill_brewer(palette="Set1")
          
    gp <- gp + facet_grid(operation + compartment ~ signal.strand)
    gp <- gp + theme_bw()

    gp <- gp + geom_text(data = stats,
                         aes(label = paste("t.test stat=", signif(statistic,5), "\n", 
                                           "p=", signif(p.value,5), "\n", sep=''),
                         x=-Inf, y=Inf, hjust=0, vjust=1, size=4))
   
    # XXX val is log10 to be placed correctly on y-axis
    text.pos <- log10(mean(df$signal[df$signal > 0]))

    # XXX label num is divided by 2 to account for number of facets
    labeler <- function(x) {
      return(c(y = text.pos, label = length(x) / 2))
    }
    
    gp <- gp + stat_summary(fun.data = labeler, geom = "text")
    
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
