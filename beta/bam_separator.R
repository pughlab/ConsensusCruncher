#!/usr/bin/env Rscript

### Create index file to separate bam file ###
## 1) Identify chr regions not part of bed -> these regions are separated by chr arm
## 2) Identify cytoband regions not part of bed -> these regions are separated by cytobands
## 3) Split cytoband regions part of bed based on bedfile coordinates

options(scipen=999)

args <- commandArgs(trailingOnly = TRUE) #input bedfile
args <- 'epic.bed'
bedfile <- read.table(args, stringsAsFactors = FALSE)
cytobands <- read.table('cytoband.txt', as.is = TRUE, stringsAsFactors = FALSE)
chr.arm <- read.table('hg19.armsizes.txt', stringsAsFactors = FALSE)

# Modify bedfile to include 2kb on each end (region where coverage fails)
bedfile$V2 <- bedfile$V2 - 2000
bedfile$V3 <- bedfile$V3 + 2000

## 1) Identify chr regions not part of bed -> these regions are separated by chr arm ##
bedfile.chr <- as.vector(unique(bedfile$V1)) # chr in bedfile
chr.arm.not.in.bed <- chr.arm[!chr.arm$V1 %in% bedfile.chr,]
chr.arm.not.in.bed$V2[which(chr.arm.not.in.bed$V2 == 1)] <- 0

## 2) Identify cytoband regions not part of bed -> these regions are separated by cytobands
# CYTOBAND REGIONS
chr.row <- match(bedfile.chr, bedfile$V1) # Rows where next chr begin
chr.row <- c(chr.row, nrow(bedfile))

# Identify probe regions 
probe.regions <- data.frame(matrix(nrow = length(chr.row) -1, ncol = 3))
for (i in seq(1, length(chr.row)-1)){
  start <- as.integer(bedfile[chr.row[i],][2])
  start <- max(0, start) # incase probe really close to start coordinate / prevent negative start 
  end <- as.integer(bedfile[chr.row[i+1] - 1,][3])
  probe.regions[i, ] <- unlist(c(bedfile[chr.row[i],][1], start, end))
}

# Overlapping probe and cytoband regions
probe.start.overlap <- unlist(apply(probe.regions, 1, function(x) which(cytobands$V1 == x[1] & cytobands$V2 < as.integer(x[2]) & cytobands$V3 > as.integer(x[2]))))
probe.stop.overlap <- unlist(apply(probe.regions, 1, function(x) which(cytobands$V1 == x[1] & cytobands$V2 < as.integer(x[3]) & cytobands$V3 > as.integer(x[3]))))
probe.overlap <- unique(c(probe.start.overlap, probe.stop.overlap))
# cytobands not in bedfile
not.bed.cytobands <- cytobands[-probe.overlap,]
not.bed.cytobands <- not.bed.cytobands[which(not.bed.cytobands$V1 %in% bedfile.chr),]
stopifnot(nrow(not.bed.cytobands) + nrow(probe.overlap) == nrow(cytobands))

# stopifnot(nrow(not.bed.cytobands) )

## 3) Split cytoband regions part of bed based on bedfile coordinates
bed.cytobands <- cytobands[probe.overlap,]
cytoband.seg <- data.frame()
for (i in seq(1, length(chr.row)-1)){
  chr.region <- bedfile[chr.row[i]:(chr.row[i+1]-1),]
  start.pos <- chr.region$V2
  end.pos <- chr.region$V3
  # find regions that don't overlap
  no.overlap <- sapply(seq(1, length(start.pos) - 1), function(x) start.pos[x + 1] - end.pos[x] > 0)
  segmentation <- chr.region[no.overlap, 3] # end position of non-overlapping region
  bed.region <- bed.cytobands[bed.cytobands$V1 == chr.region$V1[1],] # cytoband coordinates
  segmentation <- sort(as.vector(c(segmentation, bed.region$V2, bed.region$V3)))
  cytoband.seg <- rbind(cytoband.seg, t(sapply(seq(1, length(segmentation)-1), function(x) c(chr.region$V1[1], segmentation[x], segmentation[x+1]))))
}

## 4) combine chr arm, cytoband, and segmented cytoband regions
bam_sep <- rbind(chr.arm.not.in.bed[,1:3], not.bed.cytobands[,1:3], cytoband.seg)
chr_order <- unique(chr.arm$V1)
bam_sep$V1 <- factor(bam_sep$V1, levels = chr_order)
# bam_sep$V2 <- mixedsort(as.vector(bam_sep$V2))
library(gtools)
bam_sep$V2 <- factor(bam_sep$V2, levels = unique(mixedsort(as.vector(bam_sep$V2))))
bam_separator <- bam_sep[order(bam_sep$V1, bam_sep$V2),]

bam_separator <- paste(bam_separator$V1, bam_separator$V2, bam_separator$V3, sep = "_")

write.table(bam_separator, paste0(strsplit(args, '.bed')[[1]], '_separator_bed.txt'), sep = '\n', quote = FALSE, row.names = FALSE, col.names = FALSE)
