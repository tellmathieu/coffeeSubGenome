#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("You need to have two arguments.", call.=FALSE)
} else if (length(args)==2) {
  buscoCSV = args[1]
  buscoPDF = args[2]
}

library("ggplot2")
library("reshape2")
library("tidyverse")

# get busco csv
busco_df <- read.csv(buscoCSV, header = TRUE)
busco_df$Strain <- as.factor(busco_df$Strain)
busco_df.melted <- melt(busco_df, id.vars = "Strain")
busco_df.melted$variable <-relevel(busco_df.melted$variable, "Missing")

# stacked bar plot
busco_plot <- ggplot(busco_df.melted, aes(x=Strain, fill=fct_rev(variable), y=value)) +
 geom_bar(position= "stack", width = 0.7, stat="identity") +
 labs(x = "Strain", y = "BUSCO", fill = "Type") +
 theme_bw() +
 theme(axis.text.x = element_text(angle=45, hjust=1, size = 12), axis.text.y = element_text(size = 12), axis.title=element_text(size=12))

# save plot
pdf(buscoPDF,width=8,height=5,paper='special')
print(busco_plot)
dev.off()
