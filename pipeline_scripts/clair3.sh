#!/bin/bash
#SBATCH --gres=gpu

singularity pull docker://hkubal/clair3:latest

singularity exec clair3_latest.sif opt/bin/run_clair3.sh \
  --bam_fn=/scratch/prj/mmg_human_chrm_seq/ultra/ultra_big.bam \
  --ref_fn=/scratch/prj/mmg_human_chrm_seq/ultra/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
  --threads=16 \               
  --platform="ont" \
  --model_path="/opt/models/r941_prom_sup_g5014/" \
  --output="/scratch/prj/mmg_human_chrm_seq/ultra/vcf/"
