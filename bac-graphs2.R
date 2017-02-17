#!/usr/bin/env Rscript

gg_colour_hue <- function(n) {
    hues = seq(15, 375, length=n+1)[1:n]
    c(hcl(h=hues, l=65, c=100))
}


library("ggplot2")
library("grid")
library("Cairo")

args <- commandArgs(TRUE)

if (length(args) < 9) {
    print("usage: prog title output dvc mpc dtc gc hom con conorder")
    quit()
}

title <- args[1]
output <- args[2]
dvcf <- args[3]
mpcf <- args[4]
dtcf <- args[5]
gcf  <- args[6]
homf <- args[7]
conf <- args[8]
conorder <- args[9]

dvc <- read.table(dvcf, col.names=c("chr", "pos", "cov"))
mpc <- read.table(mpcf, col.names=c("chr", "pos", "cov"))
dtc <- read.table(dtcf, col.names=c("chr", "pos", "cov"))
gc  <- read.table(gcf,  col.names=c("chr", "pos", "gc"))
hom <- read.table(homf, col.names=c("chr", "base", "length",
                            "beg", "end"))
con <- read.table(conf, header=TRUE)
conorder <- unlist(strsplit(conorder, split=","))
con$file <- factor(con$file, levels=c(conorder, "space"))

plotmargin <- unit(c(0,0,0,4), "mm")

conp <- ggplot(con, aes(ymin=track+0.1, ymax=track+0.9,
                        xmin=beg, xmax=end, fill=file)) +
    guides(fill=FALSE) +
    xlab(NULL) +
    ylab(NULL) +
    scale_x_continuous(labels=NULL) +
    scale_y_continuous(breaks=c(0.5,2.0,4.5),
                       labels=c("Discovar", "Falcon", "Supernova")) +
    scale_fill_manual(values=c(gg_colour_hue(3),"light grey")) +
    theme(axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size=8),
          plot.margin=plotmargin) +
    geom_rect() +
    ggtitle(expression(paste(italic("S. verrucosum"), " BAC 22")))

dvcmax <- 1000
dvcp <- ggplot(dvc, aes(x=pos, ymin=0, ymax=cov)) +
    xlab(NULL) +
    ylab("Paired end") +
    scale_x_continuous(labels=NULL) +
    scale_y_log10(breaks=c(1,10,100,1000)) +
    theme(axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=9,vjust=0.1),
          plot.margin=plotmargin) +
    geom_ribbon()

mpcp <- ggplot(mpc, aes(x=pos, ymin=0, ymax=cov)) +
    xlab(NULL) +
    ylab("Mate pair") +
    scale_x_continuous(labels=NULL) +
    theme(axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=9,vjust=0.1),
          plot.margin=plotmargin) +
    geom_ribbon()

dtcp <- ggplot(dtc, aes(x=pos, ymin=0, ymax=cov)) +
    xlab(NULL) +
    ylab("Dovetail") +
    scale_x_continuous(labels=NULL) +
    theme(axis.ticks.x = element_blank(),
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=9,vjust=0.1),
          plot.margin=plotmargin) +
    geom_ribbon()

pad <- 50
gcp <- ggplot(gc) +
    geom_rect(data=subset(hom,length > 4 & base != "N"),
              aes(ymin=0,ymax=100, xmin=beg-pad, xmax=end+pad,
                  fill=base, alpha=length)) +
    geom_line(aes(x=pos, y=gc*100)) +
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
    theme(plot.margin=plotmargin,
          axis.text.y = element_text(size=8),
          axis.title.y = element_text(size=9,vjust=0.1))

dvcg <- ggplotGrob(dvcp)
mpcg <- ggplotGrob(mpcp)
dtcg <- ggplotGrob(dtcp)
gcg <- ggplotGrob(gcp)
cong <- ggplotGrob(conp)

pdf(file=sprintf("%s.pdf",output), height=4.5, width=12)
grid.draw(rbind(cong, dvcg, mpcg, dtcg, gcg, size="last"))
dev.off()

CairoPNG(file=sprintf("%s.png",output), height=4.5, width=12,
    units="in", dpi=80)
grid.draw(rbind(cong, dvcg, mpcg, dtcg, gcg, size="last"))
dev.off()
