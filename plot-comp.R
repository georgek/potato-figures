#!/usr/bin/env Rscript

args <- commandArgs()

matrixFile <- args[1]
maxMul <- as.integer(args[2])
minCov <- as.integer(args[3])
covBands <- as.integer(args[4])
readsName <- args[5]
genName <- args[6]

transpose <- FALSE
isReads <- TRUE

library("reshape2")
library("ggplot2")

## returns list of n colours equally spaced in HSL, like ggplot default
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)[1:n]
    hcl(h=hues, l=60, c=120)
}

matrix <- read.table(matrixFile)
if (transpose) {
    matrix <- t(matrix)
}
matrix <- matrix[1:(maxMul+1),]
lastcolumn <- rowSums(matrix[,-seq_len(minCov+covBands-2)])
matrix <- matrix[,(minCov+1):(minCov+covBands-1)]
matrix <- cbind(matrix, lastcolumn)

counts <- data.frame(cbind(as.matrix(seq(0,maxMul)), matrix))
names(counts) <- c("multiplicity",
                   sapply(seq(minCov,minCov+covBands-1),
                          function (n) sprintf("%dx", n)))
names(counts)[length(names(counts))] <-
    sprintf("%s+", names(counts)[length(names(counts))])
counts <- melt(counts, id.vars=c("multiplicity"),
               variable.name="coverage", value.name="count")

## a big peak at 1 is often expected and not interesting, but any other peak is
totals <- rowSums(matrix)
peaksx <- which(diff(sign(diff(totals))) == -2)
peaksy <- totals[peaksx]
if (!isReads) {
    peaky <- max(totals)
} else if (peaksx[1] == 1) {
    peaky <- max(peaksy[-1])
} else {
    peaky <- max(peaksy)
}

p <- ggplot(counts, aes(x=multiplicity, y=count, fill=coverage))
p <- p + geom_bar(stat="identity")
if (minCov == 0) {
    p <- p + scale_fill_manual(values=c("#444444", gg_color_hue(covBands-1)))
} else {
    p <- p + scale_fill_manual(values=c(gg_color_hue(covBands)))
}
p <- p + coord_cartesian(xlim=c(0,maxMul))
p <- p + coord_cartesian(ylim=c(0,peaky*1.1))
p <- p + labs(title=sprintf("%s", genName),
              x="k-mer multiplicity", y="Number of distinct k-mers",
              fill="Coverage")
p <- p + theme(legend.justification=c(1,1),legend.position=c(1,1))
