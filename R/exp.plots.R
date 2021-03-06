# exp.plots.R 
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 552 $'
#
# modmap pipeline for plotting signals relative to  expression levels

library(ggplot2)
library(Cairo)
library(plyr)
library(gridExtra)

# get the filename
output = commandArgs(trailingOnly=TRUE)
if (length(output) != 3) {
   stop("usage: Rscript exp.plots.R infile sample.name output.dir")
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


exp.corr.plot <- function(df, sample.name, ... ) {

    # subset data for plotting on log scale - extreme small numbers
    df <- subset(df, region.score >= 1 & signal >= 1)

    # remove outliers
    box.stats <- boxplot.stats(df$region.score)$stats
    lower.adj <- box.stats[1]
    upper.adj <- box.stats[5]

    df <- subset(df, region.score > lower.adj & region.score < upper.adj)

    # calculate correlations
    corrs <- ddply(df, operation ~ region.strand + signal.strand,
                   summarise,
                   cor = signif(cor.test(region.score, signal, method = 'pears')$estimate, 3),
                   pvalue = signif(cor.test(region.score, signal, method = 'pears')$p.value, 4),
                   num.regions = length(signal),
                   signal.sum = sum(signal))

    gp <- ggplot(data = df,
                 aes(x=log2(region.score),
                     y=log2(signal)))

    gp <- gp + geom_point(show_guide = FALSE)
    gp <- gp + geom_smooth(method = "lm") 
    gp <- gp + facet_grid(operation ~ region.strand + signal.strand)

    # x = -Inf and y = Inf put the label in the top left, see
    # https://groups.google.com/forum/#!topic/ggplot2/rFckCiNoU7U
    gp <- gp + geom_text(data = corrs,
                         aes(label = paste("pearson.r=", cor, "\n", 
                                           "p=", pvalue, "\n",
                                           "n.regions=", num.regions, "\n",
                                           "sum=", signal.sum, sep="")),
                         x=-Inf, y=Inf, hjust=0, vjust=1)

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

gp.exp.plot <- exp.corr.plot(df, sample.name)

# write the files
exp.pdf.filename <- paste(output.dir, '/', 'modmap.expression',
                      '.', sample.name, '.pdf', sep='')
ggsave(filename = exp.pdf.filename, plot = gp.exp.plot, 
       height = 8.5, width = 11, device = CairoPDF)

