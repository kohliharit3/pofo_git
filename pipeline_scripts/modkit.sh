modkit dmr pair -a h1_pileup.bed.gz -b h2_pileup.bed.gz -o dmr_result.bed \ # Two haplotype-split pileups created from split haplotagged BAMs
 -r /media/harit/biggie/fast5_ont/dmr_analysis/new_icr.bed \  # List of iDMRs
 --ref /media/harit/biggie/fast5_ont/ref_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \ # Reference genome
 --base C --threads 16 --log-filepath dmr.log
