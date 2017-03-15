#!/usr/bin/env Rscript

library("ggplot2")

t <- read.table("output/transcript-asm-cov-id-count",
                col.names=c("asm","cov","id","count"),
                colClasses=c("factor", "numeric", "numeric", "numeric"))

ggplot(subset(t, (cov*100)%%5 == 0 & id%%5 == 0), aes(x=id, y=count, fill=asm)) +
    geom_bar(stat="identity", position="dodge") +
    geom_bar(stat="identity", position="dodge", colour="black", show_guide=FALSE) +
    facet_grid(. ~ cov)
