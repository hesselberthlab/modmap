# origin.plots.R
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 552 $'
#
# modmap pipeline for plotting signals relative to origins

library(ggplot2)
library(plyr)
library(RColorBrewer)
library(Cairo)

# get the filename
output = commandArgs(trailingOnly=TRUE)
if (length(output) != 3) {
   stop("usage: Rscript origin.plots.R infile sample.name output.dir")
}

infile = output[1]
sample.name = output[2]
output.dir = output[3]

COLNAMES <- c('nuc','offset','size','count','freq','total.sites','direction','sample.id')
df <- read.table(infile, col.names=COLNAMES)
if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

ggplot.origin.plot <- function(df, sample.name, ... ) {

    gp <- ggplot(data = df,
                 aes(nuc, freq, offset, direction))

    gp <- gp + geom_bar(stat = 'identity', aes(fill = factor(nuc)))
    gp <- gp + facet_grid(direction ~ offset)
    gp <- gp + theme(legend.position = 'none')

    gp <- gp + theme_bw()

    # axis labels 
    gp <- gp + xlab('Nucleotides')
    gp <- gp + ylab('Frequency')

    # add title
    title.top = paste('modmap origin analysis (sample ',
                      sample.name, ')', sep='')
    title.bottom = ""
    title = paste(title.top, title.bottom, sep='\n')
    gp <- gp + ggtitle(title)

    return(gp)
}


gp <- ggplot.origin.plot(df, sample.name)

# write the file
pdf.filename <- paste(output.dir, '/', 'modmap.origin.analysis',
                      '.', sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, plot = gp, device = CairoPDF)
