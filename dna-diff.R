#!/usr/bin/env Rscript

library("ggplot2")
library("reshape2")
library("grid")

args <- commandArgs(TRUE)

if (length(args) < 2) {
    print("usage: dna-diff.R output indels snps")
    quit()
}

output <- args[1]
indelsf <- args[2]
snpsf <- args[3]

indels <- read.table(indelsf, col.names=c("asm", "base", "count", "bcount"))
snps <- read.table(snpsf, col.names=c("asm", "type", "count"))

indelsscp <- ggplot(indels, aes(x=base, y=(count/bcount)*10000, fill=base)) +
    geom_bar(colour="black", stat="identity", position="dodge") +
    facet_grid(. ~ asm) +
    scale_y_continuous(name="Insertions/deletions per 10kbp",
                       limits=c(-1,2.2),
                       breaks=c(-1,0,1,2)) +
    xlab("Base") +
    ggtitle("Indels") +
    guides(fill=FALSE, colour=FALSE)

indelsp <- ggplot(indels, aes(x=base, y=count, fill=base)) +
    geom_bar(colour="black", stat="identity", position="dodge") +
    facet_grid(. ~ asm) +
    xlab("Base") +
    ylab("Insertion/deletion count") +
    ggtitle("Indels") +
    guides(fill=FALSE, colour=FALSE)

snpsp <- ggplot(snps, aes(x=type, y=count, fill=type)) +
    geom_bar(colour="black", stat="identity", position="dodge") +
    facet_grid(. ~ asm) +
    xlab("Type of SNP") +
    ylab("Number of SNPs") +
    ggtitle("SNPs") +
    guides(fill=FALSE, colour=FALSE) +
    theme(axis.text.x = element_text(size=8, angle=90, hjust=1, vjust=0.5))

pdf(file=sprintf("%s/indels.pdf",output), height=4, width=4)
print(indelsp)
dev.off()

pdf(file=sprintf("%s/indels-scaled.pdf",output), height=4, width=4)
print(indelsscp)
dev.off()

pdf(file=sprintf("%s/snps.pdf",output), height=6, width=6)
print(snpsp)
dev.off()
