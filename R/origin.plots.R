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

COLNAMES <- c('nuc','offset','size','count','freq',
              'total.sites','direction','sample.id',
              'trep', 'flank.size')
df <- read.table(infile, col.names=COLNAMES)
if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

# stats for the table
max.trep <- df$trep[1]
flank.size <- df$flank.size[1]

nuc.freq.ggplot <- function(df, sample.name, ... ) {

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
    title.top = paste('modmap origin-analysis\n',
                      '(sample = ', sample.name, 
                      ' Trep max = ', max.trep,
                      ' Flank size = ', flank.size,
                      ')', sep='')
    title.bottom = ""
    title = paste(title.top, title.bottom, sep='\n')
    gp <- gp + ggtitle(title)

    return(gp)
}

origin.count.ggplot <- function(df, sample.name, ...) {
  
  gp <- ggplot(data = df, aes(direction, count))
  gp <- gp + geom_bar(stat='identity')
  
  gp <- gp + theme_bw()
  # axis labels 
  gp <- gp + xlab('Direction')
  gp <- gp + ylab('Count')
  
  return(gp)
}

gp.nuc.freq <- nuc.freq.ggplot(df, sample.name)
gp.origin.count <- origin.count.ggplot(df, sample.name)

# write the files
freq.pdf.filename <- paste(output.dir, '/', 'modmap.origin.freq',
                      '.', sample.name, '.pdf', sep='')
ggsave(filename = freq.pdf.filename, plot = gp.nuc.freq, device = CairoPDF)

count.pdf.filename <- paste(output.dir, '/', 'modmap.origin.counts',
                           '.', sample.name, '.pdf', sep='')
ggsave(filename = freq.pdf.filename, plot = gp.origin.count, device = CairoPDF)
