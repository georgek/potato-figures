#!/usr/bin/env Rscript

library("ggplot2")
library("reshape2")
library("grid")

t <- read.table("../genes-table", header=TRUE)
tm <- melt(t, id.vars="asm", variable.name="type", value.name="genes")

p <- ggplot(tm, aes(x=asm, y=genes, fill=type)) +
    scale_fill_discrete(name="Type",
                      labels=c("cegcom"="Cegma complete",
                          "cegpar"="Cegma partial",
                          "buscom"="Busco complete",
                          "buspar"="Busco partial",
                          "busmis"="Busco missing")) +
    scale_x_discrete(name="Assembly") +
    scale_y_continuous(name="Number of genes") +
    geom_bar(stat="identity", position="dodge") +
    geom_bar(stat="identity", position="dodge", colour="black", show_guide=FALSE) +
    theme(legend.position="bottom", legend.key.size=unit(3, "mm"),
          legend.text=element_text(size=8),legend.title=element_text(size=8),
          plot.margin=unit(c(2,2,0,0), "mm")) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))

pdf(file="genes-bar.pdf", width=4, height=4)
print(p)
dev.off()
