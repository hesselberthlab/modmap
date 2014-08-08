# random_dist.R
#
# __author__ = 'Jay Hesselberth'
# __contact__ = 'jay.hesselberth@gmail.com'
# __version__ = '$Revision: 552 $'
#
# modmap pipeline for determining random distributions

library(ggplot2)
library(Cairo)

# generate observations from poisson using observed data

process.observations <- function(df) {

    obs.vals <- df$vals
    obs.lambda <- mean(obs.vals)
    num.obs <- length(obs.vals)

    # generate null distribution based on observed values
    null.vals <- rpois(num.obs, obs.lambda)

    combined.df <- data.frame(null.vals,
                          rep('pois'))

    result <- t.test(obs.vals, null.vals)
    pval <- result$p.value

    return(combined.df, pval)
}

ggplot.dist <- functions(df, ...) {

    df, pval <- process.observations(df)

    gp <- ggplot(data = df, 
                 fill = factor(dist))

    gp <- gp + facet_grid(. ~ chrom)

    gp <- geom_histgram(position = 'dodge',
                        binwidth = 1,
                        color = 'black')

    gp <- scale_fill_brewer(palette="Set1")

    #gp <- geom_text(data = pvals, 
    #                aes(label = paste("P(t.test)=", pval, sep=''),
    #                    x = -Inf, y = Inf,
    #                    hjust = 0, vjust = 0))

    gp <- gp + xlab('Bin') + ylab('Count')

    title.top = paste('modmap distribution \n(sample ',
                      sample.name, sep='')

    title = paste(title.top, title.bottom, sep='\n')
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
