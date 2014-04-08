# hist.R
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 559 $'
#
# modmap pipeline for plotting nucleotide frequencies 
# makes histogram of signals

library(ggplot2)
library(RColorBrewer)
library(Cairo)

# get the filename
output = commandArgs(trailingOnly=TRUE)
if (length(output) != 4) {
   stop("usage: Rscript hist.R bedgraph.pos bedgraph.neg sample.name output.dir")
}

bg.pos = output[1]
bg.neg = output[2]
sample.name = output[3]
output.dir = output[4]

COLNAMES <- c('chrom','start','end','count')
df.pos <- read.table(bg.pos, col.names=COLNAMES)
df.neg <- read.table(bg.neg, col.names=COLNAMES)
df.all <- rbind(df.pos, df.neg)

if (nrow(df.pos) == 0 || nrow(df.neg) == 0) {
    warning("empty data frames?")
    quit(status=0)
}
head(df)

gp <- ggplot(df.all, aes(chrom, count))
gp <- 

pdf.filename <- paste(output.dir, '/', 'modmap.histogram.',
                      sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, 
       plot = gp.nuc.freq,
       device = CairoPDF)
