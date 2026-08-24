---
title: "PofO_Git ReadMe"
author: "Harit Kohli"
format: html
jupyter: python3
---

# Parent of Origin Phasing Pipeline and Predictor

Parent of origin phasing (PofO-phasing) refers to the process of assigning parent of origin to alleles of heterozygous variants. The variant alleles are only inherited from one parent, and using the methodology of read-based phasing, can be assigned a parent-of-origin (PofO) with the help of parentally imprinted regions in the genome. The imprints, which are actually just certain differentially methylated regions, indicate from which parent a sequencing read (that spans such a region) came from.


## Pipeline

The pipeline scripts can be used to assign parent of origin to all identified heterozygous variants, provided they can be phased with an imprinted locus (iDMR). The workflow is in the diagram below.

<img width="911" height="641" alt="ont_pipeline" src="https://github.com/user-attachments/assets/b6c5907c-a777-4845-9efe-a205cefb31a7" />




## Predictor

The predictor consists of 3 files:
  1. chrom_dmrs.csv: This contains information on all imprinted regions in the genome. If any changes need      to be made, it must be ensured the format of the file is kept the same.
  2. ideal_chromosome_simulator.py: This is the backend of the tool. It reads the chrom_dmrs.csv file, and      runs simulations based on an "ideal" chromosome of 100 Megabases. These parameters are at the top of       the file, so can be changed.
  3. pofo_predictor.py: This is the executable file. It takes as input all sequencing parameters, and           outputs estimated regions that would be PofO-phased for each chromosome. 




## Algorithm

Details of each function are in the code (ideal_pofo_simulator.py). But overarching algorithm details are here.


# Algorithm for “find_best_read” 

1. Obtain starting position and a direction (“upper” or “lower”). 

2. Obtain the dataframe containing simulated reads. 

3. Obtain the known basecalling error rate, as a float between 0 and 1. 

4. Find all reads that cover the given starting position. These are called the “chosen reads”. 

5.  If the direction given is “upper”, find the read among the chosen reads which has the lowest starting position. This is the “best read”. Else, if the direction given is “lower”, find the read among the chosen reads which has the highest ending position. This is the “best read” (Fig. 14).

6.  Generate a random number between 0 and 1. If this number is less than the basecalling error rate, discard the “best read”, remove it from the “chosen reads”, and repeat step 5. Else, move on to step 7.

7.  Return the index of the “best read”. 

 

# Algorithm for “build_phase_blocks” 

1. Obtain list of iDMRs. 

2. For each iDMR, obtain the starting point and ending point of the DMR. 

3.  Initialize a dataframe “phase blocks”. 

4.  Run “find_best_read” in the upper direction for the ending point of the DMR. Now we have the best read for phasing, which spans the entire DMR. 

5.  Find the topmost variant in range, by using numpy.argwhere command to determine the variant in the given distribution that is furthest from the iDMR while being within the best read. 

6.  Set the topmost variant as the starting point. 

7.  Rerun the “find_best_read” function, this time for the new starting point, and change the topmost variant position. 

8.  Iterate until the topmost variant does not change. This means that no further variants are available within a single read. 

9.  The final determined variant is the upper limit for the phase block. 

10. Repeat steps 3 to 7 for the “lower” direction and the starting point of the DMR.  

11. Obtain the “bottommost variant”. This is the lower limit of the phase block. 

12. Add the obtained limits of the phase block to the phase block dataframe. 

The resultant phase blocks are checked for duplicates, and any duplicates are removed. The remaining phase blocks define the regions which can be PofO phased for the given chromosome and with the given iDMR distribution. 
