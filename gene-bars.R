#!/usr/bin/env Rscript

library("ggplot2")
library("reshape2")

t <- read.table("../genes-table", header=TRUE)
tm <- melt(t, id.vars="asm", variable.name="type", value.name="genes")

p <- ggplot(tm, aes(x=asm, y=genes, fill=type)) +
    xlab("Assembly") +
    ylab("Number of genes") +
    scale_fill_discrete(name="Type",
                      labels=c("cegcom"="Cegma complete",
                          "cegpar"="Cegma partial",
                          "buscom"="Busco complete",
                          "buspar"="Busco partial",
                          "busmis"="Busco missing")) +
    geom_bar(stat="identity", position="dodge", colour="grey")

pdf(file="genes-bar.pdf", width=4, height=4)
print(p)
dev.off()
