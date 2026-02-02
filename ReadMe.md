---
title: "PofO_Git ReadMe"
author: "Harit Kohli"
format: html
jupyter: python3
---

# Parent of Origin Phasing Pipeline and Predictor

Parent of origin phasing (PofO-phasing) refers to the process of assigning parent of origin to alleles of heterozygous variants. The variant alleles are only inherited from one parent, and using the methodology of read-based phasing, can be assigned a parent-of-origin (PofO) with the help of parentally imprinted regions in the genome. The imprints, which are actually just certain differentially methylated regions, indicate from which parent a sequencing read (that spans such a region) came from.


## Pipeline

The pipeline scripts can be used to assign parent of origin to all identified heterozygous variants, provided they can be phased with an imprinted region.


## Predictor

The predictor consists of 3 files:
  1. chrom_dmrs.csv: This contains information on all imprinted regions in the genome. If any changes need      to be made, it must be ensured the format of the file is kept the same.
  2. ideal_chromosome_simulator.py: This is the backend of the tool. It reads the chrom_dmrs.csv file, and      runs simulations based on an "ideal" chromosome of 100 Megabases. These parameters are at the top of       the file, so can be changed.
  3. pofo_predictor.py: This is the executable file. It takes as input all sequencing parameters, and           outputs estimated regions that would be PofO-phased for each chromosome. 
