dmr <- read.table("/media/harit/biggie/fast5_ont/dmr_analysis/dmr_result.bed", quote="", sep = "", header = FALSE, stringsAsFactors = FALSE)
View(dmr)
summary(dmr)

colnames(dmr) <- c("chrom", "start", "end", "name", "score", "strand", 
                   "h1_counts", "h1_total", "h2_counts", "h2_total",
                   "h1_percent", "h2_percent", "h1_frac_mod", "h2_frac_mod")


#View(dmr)
# Creating a new column which tells us the difference in methylation
dmr$frac_diff <- dmr$h1_frac_mod-dmr$h2_frac_mod

#View(dmr)
# Selecting relevant columns
dmr <- dmr[, c(1:4, 15, 5:14)]


#write.table(dmr, file = "~/Desktop/parent_of_origin/dmr")

# We copied the "methylated parent" variable from the list of iDMRs (new_icr.bed)
# The iDMR list file was the same one used to run "modkit pair"
methyl_pofo <- read.clipboard(header = FALSE)
dmr$methylated_parent <- methyl_pofo

#View(dmr)

#write.table(dmr, file = "./dmr.tsv", sep = "\t", col.names = TRUE)

# Identifying the actually differentially methylated iDMRs from the list
highly_dmr <- dmr[(dmr$frac_diff > 0.6 | dmr$frac_diff < -0.6), ]
#View(highly_dmr)

# Removing all NA values (all iDMRs which were not actually differentially methylated)
highly_dmr <- na.omit(highly_dmr)
highly_dmr <- sort_by(x = highly_dmr, y = highly_dmr$chrom)
#View(highly_dmr)

# Writing the output to be processed in the Python file which identifies PofO of variant alleles
write.table(highly_dmr, file = "./highly_dmr.tsv", sep = "\t", col.names = TRUE)


View(highly_dmr)
