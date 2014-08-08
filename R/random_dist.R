t # random_dist.R
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 552 $'
#
# modmap pipeline for determining random distributions

library(ggplot2)
library(Cairo)
library(plyr)

output = commandArgs(trailingOnly=TRUE)
if (length(output) != 3) {
   stop("usage: Rscript random_dist.R infile sample.name output.dir")
}

infile = output[1]
sample.name = output[2]
output.dir = output[3]

COLNAMES <- c('obs', 'chrom')

df <- read.table(infile, col.names=COLNAMES)
if (nrow(df) == 0) {
    warning("empty data frame")
    quit(status=0)
}
head(df)

# generate observations from poisson using observed data

process.observations <- function(df) {

    obs.vals <- df$vals
    obs.lambda <- mean(obs.vals)
    num.obs <- length(obs.vals)

    # generate null distribution based on observed values
    null.vals <- rpois(num.obs, obs.lambda)

    combined.df <- data.frame(null.vals,
                          rep('pois'))

    return(combined.df)
}

ggplot.dist <- functions(df, ...) {

    df <- process.observations(df)

    pvals <- ddply(df, "chrom", 
                   function(x) {
                        res <- t.test(obs ~ null, data = x)
                        with(res, data.frame(statistic, p.value))
                   })

    gp <- ggplot(data = df, 
                 fill = factor(dist))

    gp <- gp + facet_grid(. ~ chrom)

    gp <- geom_histgram(position = 'dodge',
                        binwidth = 1,
                        color = 'black')

    gp <- scale_fill_brewer(palette="Set1")

    gp <- geom_text(data = pvals, 
                    aes(label = paste("P(t.test)=", p.value, sep=''),
                        x = -Inf, y = Inf,
                        hjust = 0, vjust = 0))

    gp <- gp + xlab('Bin') + ylab('Count')

    title = paste('modmap distribution \n(sample ',
                      sample.name, sep='')

    gp <- gp + ggtitle(title)

    return(gp)
}                       

gp.dist <- ggplot.dist(df)

# write the file
pdf.filename <- paste(output.dir, '/', 'modmap.dist',
                      '.', sample.name, '.pdf', sep='')

ggsave(filename = pdf.filename, plot = gp.dist,
       height = 8.5, width = 11,
       device = CairoPDF)
