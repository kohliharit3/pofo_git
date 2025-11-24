dmr <- read.table("/media/harit/biggie/fast5_ont/dmr_analysis/dmr_result.bed", quote="", sep = "", header = FALSE, stringsAsFactors = FALSE)
View(dmr)
summary(dmr)

colnames(dmr) <- c("chrom", "start", "end", "name", "score", "strand", 
                   "h1_counts", "h1_total", "h2_counts", "h2_total",
                   "h1_percent", "h2_percent", "h1_frac_mod", "h2_frac_mod")


dmr <- dmr[, c(1:4, 13, 14, 5:12)]

View(dmr)

dmr$frac_diff <- dmr$h1_frac_mod-dmr$h2_frac_mod

View(dmr)
dmr <- dmr[, c(1:4, 15, 5:14)]


write.table(dmr, file = "~/Desktop/parent_of_origin/dmr")

x <- read.clipboard(header = FALSE)
dmr$methylated_parent <- y$V1

View(dmr)

write.table(dmr, file = "./dmr.tsv", sep = "\t", col.names = TRUE)

highly_dmr <- dmr[(dmr$frac_diff > 0.6 | dmr$frac_diff < -0.6), ]
View(highly_dmr)
highly_dmr <- na.omit(highly_dmr)
highly_dmr <- sort_by(x = highly_dmr, y = highly_dmr$chrom)
View(highly_dmr)

p_of_o_assignment <- function(frac_diff, methyl_parent) {
  if (methyl_parent=="Maternal") {
    
    if (frac_diff > 0) {
      return("h1_maternal")
    }
    else {
      return("h2_maternal")
    }
  }
  
  else {
    if (frac_diff > 0) {
      return("h2_maternal")
    }
    else {
      return("h1_maternal")
    }
  }
  
}

write.table(highly_dmr, file = "./highly_dmr.tsv", sep = "\t", col.names = TRUE)


View(highly_dmr)