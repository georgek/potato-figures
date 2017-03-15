#!/usr/bin/env Rscript

library("ggplot2")

args <- commandArgs(TRUE)

if (length(args) < 3) {
    print("usage: prog output coords ticks")
    quit()
}

output <- args[1]
coords <- args[2]
ticks <- args[3]

t <- read.table(coords,
                col.names=c("x1", "y1", "x2", "y2", "pid", "dir"))
t$dir <- factor(t$dir, levels=c("+","-"))

ticks <- read.table(ticks,
                    sep=",",
                    col.names=c("labels", "breaks"))

p <- ggplot(t, aes(x=x1,xend=x2,y=y1,yend=y2,alpha=pid,colour=dir)) +
    geom_segment() +
    geom_point() +
    geom_point(aes(x2,y2)) +
    scale_x_continuous(name="Reference") +
    scale_y_continuous(name="Assembly",
                       breaks=ticks$breaks,
                       labels=ticks$labels,
                       minor_breaks=NULL) +
    scale_colour_discrete(guide=FALSE) +
    scale_alpha_continuous(guide=FALSE) +
    theme(axis.text.y=element_text(size=6))

pdf(file=sprintf("%s-gg.pdf",output), width=7, height=6.5)
print(p)
dev.off()

png(file=sprintf("%s-gg.png",output), width=800, height=800)
print(p)
dev.off()
