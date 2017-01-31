#!/usr/bin/env Rscript

library("ggplot2")
library("grid")

args <- commandArgs(TRUE)

if (length(args) < 6) {
    print("usage: prog name dvc mpc dtc gc hom")
    quit()
}

name <- args[1]
dvcf <- args[2]
mpcf <- args[3]
dtcf <- args[4]
gcf  <- args[5]
homf <- args[6]

dvc <- read.table(dvcf, col.names=c("chr", "pos", "cov"))
mpc <- read.table(mpcf, col.names=c("chr", "pos", "cov"))
dtc <- read.table(dtcf, col.names=c("chr", "pos", "cov"))
gc  <- read.table(gcf,  col.names=c("chr", "pos", "gc"))
hom <- read.table(homf, col.names=c("chr", "base", "length", "beg", "end"))

plotmargin <- unit(c(0,5,0,5), "mm")

dvcmax <- 1000
## dvc[dvc$cov > dvcmax,]$cov <- dvcmax
dvcp <- ggplot(dvc, aes(x=pos, ymin=0, ymax=cov)) +
    xlab(NULL) +
    ylab("Discovar") +
    scale_x_continuous(labels=NULL) +
    scale_y_log10() +
    theme(axis.ticks.x = element_blank(), plot.margin=plotmargin) +
    geom_ribbon() +
    ggtitle(name)

mpcp <- ggplot(mpc, aes(x=pos, ymin=0, ymax=cov)) +
    xlab(NULL) +
    ylab("Mate pair") +
    scale_x_continuous(labels=NULL) +
    theme(axis.ticks.x = element_blank(), plot.margin=plotmargin) +
    geom_ribbon()

dtcp <- ggplot(dtc, aes(x=pos, ymin=0, ymax=cov)) +
    xlab(NULL) +
    ylab("Dovetail") +
    scale_x_continuous(labels=NULL) +
    theme(axis.ticks.x = element_blank(), plot.margin=plotmargin) +
    geom_ribbon()

pad <- 50
gcp <- ggplot(gc) +
    geom_rect(data=subset(hom,length > 4 & base != "N"),
              aes(ymin=0,ymax=1, xmin=beg-pad, xmax=end+pad,
                  fill=base, alpha=length)) +
    geom_line(aes(x=pos, y=gc)) +
    guides(alpha=FALSE, fill=FALSE) +
    scale_fill_manual(name="Base",
                      values=c(
                          "A" = "red",
                          "C" = "blue",
                          "G" = "yellow",
                          "T" = "green",
                          "N" = "black")) +
    xlab("Position") +
    ylab("GC Cont") +
    theme(plot.margin=plotmargin)

dvcg <- ggplotGrob(dvcp)
mpcg <- ggplotGrob(mpcp)
dtcg <- ggplotGrob(dtcp)
gcg <- ggplotGrob(gcp)

pdf(file=paste(name,".pdf"), height=3.5, width=12)
grid.draw(rbind(dvcg, mpcg, dtcg, gcg, size="last"))
dev.off()
