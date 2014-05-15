# origin.plots.R
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 552 $'
#
# modmap pipeline for plotting signals relative to origins

library(ggplot2)
library(Cairo)
library(gridExtra)

# get the filename
output = commandArgs(trailingOnly=TRUE)
if (length(output) != 3) {
   stop("usage: Rscript origin.plots.R infile sample.name output.dir")
}

infile = output[1]
sample.name = output[2]
output.dir = output[3]

COLNAMES <- c('nuc','offset','count','freq',
              'total.sites','max.timing',
              'flank.size','strand')

df <- read.table(infile, col.names=COLNAMES)
if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

origin.nuc.count.ggplot <- function(df, sample.name, ... ) {

    gp <- ggplot(data = df, 
                 aes(x=offset, y=count, fill=nuc))

    gp <- gp + geom_bar(stat='identity', position='dodge')
    gp <- gp + facet_grid(strand ~ max.timing + flank.size)
    gp <- gp + theme(legend.position = 'none')
    gp <- gp + scale_fill_brewer(palette="Set1")

    gp <- gp + theme_bw()

    # axis labels 
    gp <- gp + xlab('Offset')
    gp <- gp + ylab('Nucleotide counts')

    # add title
    title.top = paste('modmap origin-analysis (per-nuc counts)\n',
                      '(sample = ', sample.name, 
                      ')', sep='')
    title.bottom = "top row = max.timing; bottom row = flank.size"
    title = paste(title.top, title.bottom, sep='\n')
    gp <- gp + ggtitle(title)

    return(gp)
}

origin.agg.count.ggplot <- function(df, sample.name, ...) {
  
    gp <- ggplot(data = df,
                 aes(nuc, freq, offset, direction))

    gp <- gp + geom_bar(stat='identity')

    gp <- gp + theme_bw()
    # axis labels 
    gp <- gp + xlab('Direction')
    gp <- gp + ylab('Count')

    return(gp)
}

# gp.origin.freq <- origin.freq.ggplot(df, sample.name)
gp.origin.count <- origin.nuc.count.ggplot(df, sample.name)

# write the files
#freq.pdf.filename <- paste(output.dir, '/', 'modmap.origin.freq',
#                      '.', sample.name, '.pdf', sep='')
#ggsave(filename = freq.pdf.filename, plot = gp.nuc.freq, device = CairoPDF)

count.pdf.filename <- paste(output.dir, '/', 'modmap.origin.counts',
                           '.', sample.name, '.pdf', sep='')
ggsave(filename = count.pdf.filename, plot = gp.origin.count,
       height = 8.5, width = 11, device = CairoPDF)

