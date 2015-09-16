#!/bin/env Rscript
library(data.table)
library(ggplot2)
args <- commandArgs()
file_path <- file.path(args[6])
baseFileName <- args[7]

db <- fread(file_path)

p <- ggplot(db[db$type=="P",],aes(x=len))
p + geom_histogram(fill="white", colour="black", binwidth=10)+xlim(1,500)+ggtitle(paste(baseFileName,"Paired-End fragment length distribution"))+xlab("fragment length")
ggsave(paste(baseFileName,"PE_fraglen_distribution.png",sep = "_"))

p <- ggplot(db[db$type=="S",],aes(x=len))
p + geom_histogram(fill="white", colour="black", binwidth=10)+xlim(1,500)+ggtitle(paste(baseFileName,"Single-End fragment length distribution"))+xlab("fragment length")
ggsave(paste(baseFileName,"SE_fraglen_distribution.png",sep = "_"))

p <- ggplot(db,aes(x=len))
p + geom_histogram(fill="white", colour="black", binwidth=10)+xlim(1,500)+ggtitle(paste(baseFileName,"fragment length distribution"))+xlab("fragment length")
ggsave(paste(baseFileName,"fraglen_distribution.png",sep = "_"))

