import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Path to the "phase_sets.tsv" file output by WhatsHap
phase_sets = pd.read_table("/media/harit/seagate/fast5_ont/ultra/true_phased/phase_sets.tsv")
phase_sets = phase_sets.iloc[:, 1:6]

# Path to the phased VCF file
phased = pd.read_table("/media/harit/seagate/fast5_ont/ultra/true_phased/copy_phased.vcf")

# Extracting phase set information for each heterozygous variant
phased[["GT", "GQ", "DP", "AD", "AF", "PS"]] = phased["SAMPLE"].str.split(":", expand=True)
het_phased = phased[(phased["GT"] == "1|0") | (phased["GT"] == "0|1")]

het_phased["PS"] = het_phased["PS"].astype(int)

# Initialising new column
het_phased["pofo"] = pd.Series(dtype="str")


# Loading the DMR results (the output) from the dmr_analysis.R script
dmrs = pd.read_table("/home/harit/Desktop/parent_of_origin/highly_dmr.tsv")

dmrs["h1_maternal"] = pd.Series(dtype="boolean")


# Using the difference in methylation ("frac_diff") between haplotypes (h1 and h2) and "methylated_parent" to determine whether h1 is maternal

i = 0
indices = dmrs.index
while i < len(indices):
    if ((dmrs.loc[indices[i], "methylated_parent"] == "Maternal") & (dmrs.loc[indices[i], "frac_diff"] > 0)):
        dmrs.loc[indices[i], "h1_maternal"] = True

    elif ((dmrs.loc[indices[i], "methylated_parent"] == "Paternal") & (dmrs.loc[indices[i], "frac_diff"] < 0)):
        dmrs.loc[indices[i], "h1_maternal"] = True

    else:
        dmrs.loc[indices[i], "h1_maternal"] = False

    i += 1


def pofo_assignment(chrom, block, h1_maternal):
    block = int(block)
    chr_phased = het_phased[het_phased["CHROM"] == chrom]
    print("Currently working on phase set: ", block, "in chromosome ", chrom)
    
    variants_in_block = chr_phased[chr_phased["PS"] == block]
    print("Total variants in current block: ", len(variants_in_block["PS"]))
    
    indices = variants_in_block.index
    
    j = 0
    print("Assigning parent of origin...")
    while (j < len(variants_in_block["PS"])):

        if ((variants_in_block.loc[indices[j], "GT"] == "1|0") & (h1_maternal == True)):
            variants_in_block.loc[indices[j], "pofo"] = "maternal"

        elif ((variants_in_block.loc[indices[j], "GT"] == "0|1") & (h1_maternal == False)):
            variants_in_block.loc[indices[j], "pofo"] = "maternal"

        else:
            variants_in_block.loc[indices[j], "pofo"] = "paternal"

        j += 1
    print("Parent of origin classification done on all variants in block: ", block)

    print("Returning parent of origin labels to the heterozygous phased variant table...")

    return variants_in_block["pofo"]
    
    
    i += 1



def merge_pofo_vcf(blocks):

    indices = blocks.index

    i = 0
    while i < len(indices):
        
        print("Extracting haplotype information...")
        chrom = blocks.loc[indices[i], "chrom"]
        block = int(blocks.loc[indices[i], "phase_block"])
        h1_maternal = blocks.loc[indices[i], "h1_maternal"]
        
        labels = pofo_assignment(chrom, block, h1_maternal)
        print("Parental labels obtained")
        print(labels)
        
        variants = labels.index
        print("Total variants in current phase block: ", len(variants))

        j = 0
        while j < len(variants):
            het_phased.loc[variants[j], "pofo"] = labels[variants[j]]
            j += 1
        
        print("Added parental information to VCF file. Moving on to the next block...")
        i += 1

    print("Done!")



#Assiging iDMRs to phase sets
i = 0
indices = dmrs.index

while i < len(indices):
    #selected a DMR, now to loop over the chromosomes looking for the right one

    j = 0
    while j < len(blocks.index):
        
        #comparison
        if dmrs.loc[indices[i], "chrom"] != blocks.loc[j, "chromosome"]:
            j += 1

        #now chromosomes are matched
        elif blocks.loc[j, "from"] > dmrs.loc[indices[i], "start"]:
            j += 1

        #now the selected phase set starts BEFORE the dmr
        elif blocks.loc[j, "to"] < dmrs.loc[indices[i], "end"]:
            j += 1

        #now the block completely covers the DMR
        else:
            dmrs.loc[indices[i], "phase_block"] = blocks.loc[j, "phase_set"]
            print("Found ya! Moving on to next dmr!")
            print("Just done: ", dmrs.loc[indices[i], ["chrom", "name", "phase_block"]])

            break

    i += 1 



            
# Adding parental information to the VCF

merge_pofo_vcf(dmrs)




