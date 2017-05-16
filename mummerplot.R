#!/usr/bin/env Rscript

library("ggplot2")

args <- commandArgs(TRUE)

if (length(args) < 3) {
    print("usage: prog output coords xticks yticks")
    quit()
}

output <- args[1]
coords <- args[2]
xticks <- args[3]
yticks <- args[4]

t <- read.table(coords,
                col.names=c("x1", "y1", "x2", "y2", "pid", "dir"))
t$dir <- factor(t$dir, levels=c("+","-"))

xticks <- read.table(xticks,
                     sep=",",
                     col.names=c("labels", "breaks"))

yticks <- read.table(yticks,
                     sep=",",
                     col.names=c("labels", "breaks"))

p <- ggplot(t, aes(x=x1,xend=x2,y=y1,yend=y2,alpha=pid,colour=dir)) +
    geom_segment() +
    geom_point() +
    geom_point(aes(x2,y2)) +
    scale_x_continuous(name="Reference",
                       breaks=xticks$breaks,
                       labels=xticks$labels,
                       minor_breaks=NULL) +
    scale_y_continuous(name="Assembly",
                       breaks=yticks$breaks,
                       labels=yticks$labels,
                       minor_breaks=NULL) +
    scale_colour_discrete(guide=FALSE) +
    scale_alpha_continuous(guide=FALSE) +
    theme(axis.text.y=element_text(size=6,angle=45),
          axis.text.x=element_text(size=6,angle=45,hjust=1))

pdf(file=sprintf("%s-gg.pdf",output), width=7, height=6.5)
print(p)
dev.off()

png(file=sprintf("%s-gg.png",output), width=800, height=800)
print(p)
dev.off()
