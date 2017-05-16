#!/usr/bin/env Rscript

library("ggplot2")
library("grid")
library("gridExtra")
library("scales")
library("plyr")

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}

input <- read.table("output/all-assemblies.minlen-cumulative",
                    sep=",",
                    col.names=c("minlen", "cumulative", "fname", "name"))

input$fname <- revalue(input$fname, c("discovar-contig"="discovar"))

## "bigcomp" with the large assemblies
bigcomp <- subset(input,
                  fname %in% c("10x-bn", "falcon-dt-bn",
                               "discovar-mp-dt-bn"))
bigcomp$fname <- factor(bigcomp$fname,
                        levels=c("10x-bn",
                                 "falcon-dt-bn",
                                 "discovar-mp-dt-bn"))

bigcompplot <-
    ggplot(bigcomp, aes(x=minlen, y=cumulative, colour=fname)) +
    scale_x_continuous(trans=reverselog_trans(10),
                       breaks=c(10000000,1000000,100000,10000,1000),
                       labels=c("10Mbp", "1MBp", "100kbp", "10kbp", "1kbp")) +
    scale_y_log10(limits=c(400000,800000000),
                  breaks=c(800000000,100000000,10000000,1000000),
                  labels=c("800Mbp", "100Mbp", "10Mbp", "1Mbp")) +
    xlab("Minimal scaffold length") +
    ylab("Covered assembly length\n") +
    ggtitle(expression(paste("Comparison of ", italic("S. verrucosum"), " assemblies"))) +
    labs(colour="Assembly") +
    geom_line() +
    theme(legend.justification=c(1,0),legend.position=c(1,0))

## discovar comp
dvcomp <- subset(input,
                 fname %in% c("discovar", "discovar-mp",
                              "discovar-mp-bn", "discovar-mp-dt",
                              ## "discovar-mp-dt-jelly",
                              "discovar-mp-dt-bn"))

dvcompplot <-
    ggplot(dvcomp, aes(x=minlen, y=cumulative, colour=fname)) +
    scale_x_continuous(trans=reverselog_trans(10),
                       breaks=c(10000000,1000000,100000,10000,1000),
                       labels=c("10Mbp", "1MBp", "100kbp", "10kbp", "1kbp")) +
    scale_y_log10(limits=c(400000,800000000),
                  breaks=c(800000000,100000000,10000000,1000000),
                  labels=NULL) +
    xlab("Minimal scaffold length") +
    ylab(NULL) +
    ggtitle(expression(paste("Comparison of ", italic("S. verrucosum"), " Discovar assemblies"))) +
    labs(colour="Assembly") +
    geom_line() +
    theme(legend.justification=c(1,0),legend.position=c(1,0))

dvcompplotwaxis <- dvcompplot +
    scale_y_log10(limits=c(400000,800000000),
                  breaks=c(800000000,100000000,10000000,1000000),
                  labels=c("800Mbp", "100Mbp", "10Mbp", "1Mbp")) +
    ylab("Covered assembly length\n")

## pacbio comp
pbcomp <- subset(input,
                 fname %in% c("falcon", "falcon-dt",
                              "falcon-bn", "falcon-dt-bn"))

pbcompplot <-
    ggplot(pbcomp, aes(x=minlen, y=cumulative, colour=fname)) +
    scale_x_continuous(trans=reverselog_trans(10),
                       breaks=c(10000000,1000000,100000,10000,1000),
                       labels=c("10Mbp", "1MBp", "100kbp", "10kbp", "1kbp")) +
    scale_y_log10(limits=c(400000,800000000),
                  breaks=c(800000000,100000000,10000000,1000000),
                  labels=NULL) +
    xlab("Minimal scaffold length") +
    ylab(NULL) +
    ggtitle(expression(paste("Comparison of ", italic("S. verrucosum"), " Falcon assemblies"))) +
    labs(colour="Assembly") +
    geom_line() +
    theme(legend.justification=c(1,0),legend.position=c(1,0))

pbcompplotwaxis <- pbcompplot +
    scale_y_log10(limits=c(400000,800000000),
                  breaks=c(800000000,100000000,10000000,1000000),
                  labels=c("800Mbp", "100Mbp", "10Mbp", "1Mbp")) +
    ylab("Covered assembly length\n")

## pdf(file="cumulative-content.pdf", width=15, height=5)
## grid.arrange(bigcompplot, dvcompplot, pbcompplot, ncol=3)
## dev.off()

bigcompplotg <- ggplotGrob(bigcompplot)
dvcompplotg <- ggplotGrob(dvcompplot)
pbcompplotg <- ggplotGrob(pbcompplot)

pdf(file="cumulative-content.pdf", width=15, height=5)
grid.draw(cbind(bigcompplotg, dvcompplotg, pbcompplotg, size="last"))
dev.off()
