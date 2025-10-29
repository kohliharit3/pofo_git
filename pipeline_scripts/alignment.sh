#!/bin/bash

module load dorado
dorado aligner ../reference_genome ../data/basecalled_files -o ../data/alignment -v --emit-summary
