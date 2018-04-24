#!/usr/bin/env Rscript

### Create index file to separate bam file ###
# 1) Use cytoband coordinates for regions not overlapping with bedfile
# 2) Use bedfile coordinates to split overlapping cytoband regions

options(scipen=999)

args <- commandArgs(trailingOnly = TRUE) #input bedfile
bedfile <- read.table(args, stringsAsFactors = FALSE)
cytobands <- read.table('cytoBand.txt', as.is = TRUE, stringsAsFactors = FALSE)  # Need cytoband text file in same folder

# 1)  Use cytoband coordinates for regions not overlapping with bedfile
overlap.region.index <- unique(unlist(apply(bedfile, 1, function(x) which(cytobands$V1 %in% x[1] & cytobands$V2 < as.integer(x[2]) & cytobands$V3 > as.integer(x[2])))))
overlap.regions <- cytobands[overlap.region.index,]
non.overlap.regions <- cytobands[-overlap.region.index,]

# 2) Use bedfile coordinates to split overlapping cytoband regions
bed.split.cyto.regions <- data.frame(V1=NULL, V2=NULL, V3=NULL)
for (chr in unique(overlap.regions$V1)){
  cyto.overlap.chr <- overlap.regions[overlap.regions$V1 == chr,]
  bed.overlap.chr <- bedfile[bedfile$V1 == chr,]
  
  for (i in 1:nrow(cyto.overlap.chr)) {
    bed.overlap.cyto <- bed.overlap.chr[bed.overlap.chr$V2 >= cyto.overlap.chr[i, "V2"] & bed.overlap.chr$V3 <= cyto.overlap.chr[i, "V3"],]
    if (nrow(bed.overlap.cyto) > 5){
      # Divide bed coors overlapping cytoband region into 5 sections
      coor.overlap.section <- floor(nrow(bed.overlap.cyto)/5)  # number of bed coors in overlap section
      bed.split.cyto <- data.frame(V1=chr, V2=cyto.overlap.chr[i, 2], V3=bed.overlap.cyto[coor.overlap.section, 2])
      for (j in 2:4){
        bed.split.cyto <- rbind(bed.split.cyto,
                                data.frame(V1=chr, V2=bed.overlap.cyto[coor.overlap.section * (j-1), 2], V3=bed.overlap.cyto[coor.overlap.section * j, 2]))
      }
      bed.split.cyto <- rbind(bed.split.cyto,
                              data.frame(V1=chr, V2=bed.overlap.cyto[coor.overlap.section * 4, 2], V3=cyto.overlap.chr[i,3]))
    }
    else {
      bed.split.cyto <- cyto.overlap.chr[i, 1:3]
    }
    bed.split.cyto.regions <- rbind(bed.split.cyto.regions, bed.split.cyto)
  }
}

new.bedfile <- rbind(bed.split.cyto.regions, non.overlap.regions[,1:3])
new.bedfile <- new.bedfile[order(new.bedfile$V1, new.bedfile$V2),]
colnames(new.bedfile) <- c("chr", "start", "end")
# new.bedfile <- within(new.bedfile, strand <- rep('+', nrow(new.bedfile)))
new.bedfile <- within(new.bedfile, coor <- apply(new.bedfile, 1, function(x) paste(trimws(x[1:3]), collapse = '_')))

write.table(x = new.bedfile, file = paste0(strsplit(args, '[.]bed')[[1]], '_cytoband_bed.txt'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
