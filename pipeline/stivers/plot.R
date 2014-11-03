library(ggplot2)
library(dplyr)

colnames <- c('pos','count','sample','strand',
              'align.mode','treatment','hyb','exp.name',
              'scale')

df <- read.table('combined.tab.gz', col.names=colnames, sep='\t')

df$sample.num <- as.integer(str_sub(df$sample,3,4))

df <- subset(df, sample.num > 40 & scale == "raw")
gp <- ggplot(data=df, aes(x = pos, y = count,
                          color=factor(treatment)))
gp <- gp + theme_bw()

gp <- gp + geom_line(size = 0.5, alpha = 0.5) + facet_grid(exp.name ~ strand + align.mode, 
                                    scales="free_y") + scale_color_brewer(palette="Set1")

gp <- gp + xlab("Position (bp)") + ylab("Coverage (raw)")

title <- "Uracil Excision-seq / Lockdown HIV enrichment"
sub.title <- "Coverage of pNL4-3-∆E-GFP sequence"

gp <- gp + ggtitle(paste(title, sub.title, sep="\n"))

# gp <- gp + xim(7500,10000)

gp
