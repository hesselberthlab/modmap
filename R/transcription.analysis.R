# transcription.plots.R 
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
library(dplyr)

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

add.missing.data <- function(df) {
  
  df.nonone <- df %>% filter(region.strand != 'region-none')
  df.neg <- df %>% filter(region.strand == 'region-none')
  df.neg$region.strand <- 'region-neg'
  df.pos <- df %>% filter(region.strand == 'region-none')
  df.pos$region.strand <- 'region-pos'
  
  dfx <- rbind(df.nonone, df.neg, df.pos)

  missing.row <- c('.','.','region-neg','mito-genic','signal-neg',
                         'count',0,'norm','mito','genic')
  dfx <- rbind(dfx, missing.row)
  
  missing.row <- c('.','.','region-neg','mito-genic','signal-neg',
                         'sum',0,'norm','mito','genic')
  dfx <- rbind(dfx, missing.row)
  
  missing.row <- c('.','.','region-neg','mito-genic','signal-pos',
                         'count',0,'norm','mito','genic')
  dfx <- rbind(dfx, missing.row)
  
  missing.row <- c('.','.','region-neg','mito-genic','signal-pos',
                         'sum',0,'norm','mito','genic')
  dfx <- rbind(dfx, missing.row)
             
  dfx$signal <- as.numeric(dfx$signal)
  return(dfx)
}

txn.box.plot <- function(df, sample.name, ... ) {
  
  # add label columns
  compartment <- sapply(as.character(df$region.type), 
                        function (x) strsplit(x, '-')[[1]][1])
  df$compartment <- as.factor(compartment)
  
  region <- sapply(as.character(df$region.type), 
                   function (x) strsplit(x, '-')[[1]][2])
  df$region <- as.factor(region)
  
  df <- droplevels(add.missing.data(df))

  # calculate t.tests among groups
  stats <- ddply(df,
                operation + compartment ~ 
                  signal.strand,
                 function (x) {
                   obj <- try(t.test(signal ~ region, data = x))
                   if (is(obj, "try-error")) {
                     pval <- NA
                   } else {
                     pval <- obj$p.value}},
                 .drop = FALSE)
  
  gp <- ggplot(data = df,
               aes(factor(region),
                   log10(signal)))
  
  gp <- gp + geom_boxplot(aes(fill = factor(region)),
                          show_guide = FALSE)
  gp <- gp + scale_fill_brewer(palette="Set1")
  
  gp <- gp + facet_grid(operation + compartment ~ 
                          signal.strand)
  gp <- gp + theme_bw()
  
  gp <- gp + geom_text(data = stats,
                       aes(label = paste("p=", signif(V1,5), "\n", sep=''),
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
