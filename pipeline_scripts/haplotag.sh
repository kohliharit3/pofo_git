whatshap haplotag /
	-o vcf/haplotagged.bam /
	--reference ../GCA_000001405.15_GRCh38_no_alt_analysis_set.fna vcf/phased.vcf.gz big.bam /
	--output-threads=16 /
	--output-haplotag-list vcf/haplotag.txt /
	--ignore-read-groups
