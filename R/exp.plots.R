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

COLNAMES <- ('region.name','region.score','region.strand',
             'signal.strand','operation','signal')

df <- read.table(infile, col.names=COLNAMES)
if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

# subset data
df <- subset(df, region.score > 1 & signal > 1)

exp.corr.plot <- function(df, sample.name, ... ) {

    gp <- ggplot(data = df, 
                 aes(x=log2(region.score),
                     y=log2(signal)))

    gp <- gp + geom_point()
    gp <- gp + geom_smooth(method = "lm") 
    gp <- gp + facet_grid(operation ~ region.strand + signal.strand,
                          margins = TRUE)
    gp <- gp + theme(legend.position = 'none')
    # gp <- gp + scale_fill_brewer(palette="Set1")

    gp <- gp + theme_bw()

    # axis labels 
    gp <- gp + xlab('log2(FPKM)')
    gp <- gp + ylab('log2(calculated signal)')

    # add title
    title.top = paste('modmap expression-analysis \n',
                      'sample = ', sample.name, sep='')
    title.bottom = "top row = region.strand; bottom row = signal.strand"
    title = paste(title.top, title.bottom, sep='\n')
    gp <- gp + ggtitle(title)

    return(gp)
}

# gp.origin.freq <- origin.freq.ggplot(df, sample.name)
gp.exp.plot <- exp.corr.plot(df, sample.name)

# write the files
exp.pdf.filename <- paste(output.dir, '/', 'modmap.expression',
                      '.', sample.name, '.pdf', sep='')
ggsave(filename = exp.pdf.filename, plot = gp.exp.plot, 
       height = 8.5, width = 11, device = CairoPDF)

