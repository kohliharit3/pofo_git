#!/bin/bash

module load dorado/0.7.2-gcc-13.2.0
dorado basecaller ../ont_models/dna_r9.4.1_e8_sup@v3.3 ../data/pod5_files/output.pod5 -x cpu -v > calls.bam
