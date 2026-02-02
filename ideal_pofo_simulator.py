import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt
import warnings
import seaborn as sn
from functools import reduce
import time
from scipy.optimize import curve_fit
from sklearn import preprocessing
import bisect
from bisect import bisect_left
import matplotlib.patches as mpatches
import progressbar

def sigmoid(x, L, k, x0):
    y = np.float64(L / (1 + np.exp(-k*(x-x0))))
    return y

def inverse (x, m=1, x0= 0, c=0):
    y = np.exp(-m*(x-x0)) + c
    return y

def reduced_inverse (x, m=1, x0=0):
    y = np.exp(-m*(x-x0))
    return y


test_size = 100000000

chrom_dmrs = pd.read_csv("chrom_dmrs.csv")



test_icrs = pd.DataFrame(columns=["chrom", "start", "end", "name"])
test_icrs.loc[len(test_icrs)] = ["test", (test_size/2)-10, (test_size/2)+10, "ICR"]


#chrom_sizes = pd.read_csv("chrom_size.csv")
#lcrs = pd.read_table("/home/harit/Desktop/app_michael/hg38_gaps.bed", header=None)
#lcrs.columns = ["chrom", "start", "end"]
#lcrs["size"] = lcrs["end"]-lcrs["start"]
#new_dmrs = pd.read_csv("dmr_list.csv")
#icrs = pd.read_csv("germline_known_icrs.csv")
#het_phased = pd.read_csv("het_phased_all.csv")
#hq_variants = het_phased[het_phased["QUAL"]>19.45]

def create_simulated_reads(chrom='test', read_mu=8.62, read_sigma=1.11, mean_depth=35, show_plot=False, test_size=test_size):
    
    #First, we find the total size of the chromosome, giving us the upper range of "from"
    if chrom == 'test':
        chrom_size = test_size
        
    else:
        chrom_size = chrom_sizes.loc[chrom_sizes["chrom"]==chrom, "length"]
    
    #Now, for estimating number of reads
    mean_length = np.exp(read_mu + (read_sigma*read_sigma/2))
    num_reads = int(chrom_size * mean_depth / mean_length)
    #print("Mean length: ", mean_length)
    
    simulated_reads = pd.DataFrame(columns=["from", "to", "length"])
    
    #Now randomize the locations
    simulated_reads["from"] = np.rint(np.random.uniform(low=1, high=chrom_size, size=num_reads))    
    simulated_reads["length"] = np.random.lognormal(mean=read_mu, sigma=read_sigma, size=num_reads)
    
    #Rounding the lengths to integers
    simulated_reads["length"] = np.rint(simulated_reads["length"])
    simulated_reads["to"] = simulated_reads["from"] + simulated_reads["length"]

  
    #print("Created simulated reads for chromosome: IDEAL with mean length = ", mean_length, " and depth = ", mean_depth)
    #print("Number of reads created = ", num_reads)

    if show_plot==True:
        plt.hist(simulated_reads["from"], bins = int(chrom_size/100000))
        plt.suptitle("Location of reads in: "+chrom)
        plt.title("Each bin represents 100 kilobases (kb)", fontsize=10)
        plt.xlabel("Genomic Position (100 Megabases)")
        plt.ylabel("Frequency")
        plt.show()
        #plt.savefig("../results/chr8_read_distribution.jpg", dpi=1200)
        
        plt.hist(simulated_reads["length"], bins = 500)
        plt.suptitle("Read Length Distribution for: "+chrom)
        plt.title("Log mean length =  "+str(read_mu)+", SD = "+str(read_sigma)+", Depth = "+str(mean_depth))
        plt.xlabel("Read Lengths (bases)")
        plt.ylabel("Frequency")
        plt.show()
        #plt.savefig("../results/chr8_read_length_distribution.jpg", dpi=1200)
        #print(simulated_reads)
        
    return simulated_reads

def phasing_simulator (icrs, mean_depth, chrom="all", simulate_variants=True, pop_vcf=None, variant_per=1000, 
                       var_mu=5.82, var_sigma=1.528, read_mu = 8.51, read_sigma=1.158, error_rate=0.10):
    
    
    simulated_reads = create_simulated_reads(mean_depth=mean_depth, read_mu=read_mu, read_sigma=read_sigma)

    part_icrs = test_icrs
    
    #The following expects a VCF file. We will be using simulated variants from now on.
    #part_vcf = vcf_partition(variant_distribution, chromosome)

    if simulate_variants == False:
        part_vcf = pop_vcf

    else:
        part_vcf = create_variant_distribution(chrom, variant_per=variant_per, var_mu=var_mu, var_sigma=var_sigma)
        #The above creates a list of positions of variants
        
        #part_vcf = variants
        #warnings.warn("NOT SIMULATING VARIANTS. THIS IS THE TEST VERSION. COMMENT THIS OUT.")
        #This is to test if there is still randomness if variants are fixed.

    
    #Now, we add the ICR starting and ending points to the variant list, to avoid overlaps in the phase sets
    #part_vcf = np.append(part_vcf, part_icrs["start"])
    #part_vcf = np.append(part_vcf, part_icrs["end"])
    #part_vcf = np.unique(part_vcf)
    #part_vcf = np.sort(part_vcf)
    
    blocks = build_phase_blocks(dmrs=part_icrs, variants=part_vcf, reads=simulated_reads, error_rate=error_rate)
    #print("Basecalling error rate: ", error_rate)
    #print(blocks)
    frac = percentage_phased(blocks, chrom)
    return frac, blocks


def percentage_phased (blocks, chrom="test"):

    blocks = blocks.sort_values("size")

    #The next line was the default way of removing duplicates, but doesn't work when phase sets end at ICRs themselves
    #We need a better function to remove overlaps
    #blocks.drop_duplicates(keep="first", inplace=True)
    blocks.drop_duplicates(subset="start", keep="last", inplace=True)
    blocks.drop_duplicates(subset="end", keep="first", inplace=True)
    
    total = 0
    for index in blocks.index:
        total = total + blocks.loc[index, "size"]

    #print("Total bases which can be assigned parent of origin: ", total)

    #We find size of the current chromosome to calculate fraction
    chrom_size = test_size
    
    frac = total * 100 / chrom_size
    #print("Percent of chromosome which can be assigned PofO: ", frac)
    #print(frac)
    return frac


def find_best_read (pos, direction, reads, base_error=0.1):    
    chosen_reads = reads[(reads["from"] <= pos) & (reads["to"] >= pos)]
    #print("Number of chosen_reads = ", len(chosen_reads))
    
    #First, we check if the chosen_indices list isn't empty
    if chosen_reads.empty:
        #print("No read covers this ICR! Rerun the simulation.")
        #warnings.warn("An ICR was not spanned by any read in the simulation, so it was skipped. Rerun recommended for better results.")
        error_text = "no_read"
        return np.float64(pos)

    else:    
        #Now for determining the best for our purposes:
        if (direction == "up"):
            chosen_reads = chosen_reads.sort_values(by="from", ascending=True)
            #best_index = chosen_reads.iloc[0].name #the left-most "from" value is picked
                
        elif (direction == "down"):
            chosen_reads = chosen_reads.sort_values(by="to", ascending=False)
            #best_index = chosen_reads.iloc[0].name   #the right-most "to" value selected
    
        else:
            print("Couldn't understand direction, sorry!")

        #Error rate based random removal
        i = 0        
        roulette_value = np.random.uniform(0, 1, 1)                   #Randomly picks a number from 0 to 1 (pulling the trigger)
        while(roulette_value < base_error):                           #Checking if we're lucky (did we get the bullet?) 
            roulette_value = np.random.uniform(0, 1, 1)               #Respinning the wheel
            i += 1                                                    #Removing the best read (if we got the bullet)

        if i >= len(chosen_reads):
            limit_pos = np.float64(pos)
        
        elif (direction == "up"):
            limit_pos = np.float64(chosen_reads.iloc[i]["from"])
                
        elif (direction == "down"):
            limit_pos = np.float64(chosen_reads.iloc[i]["to"])
    
        else:
            print("Couldn't understand direction, sorry!")

        return limit_pos     


def build_phase_blocks(dmrs, variants, reads, error_rate=0.10):

    blocks = pd.DataFrame(columns = ["start", "end", "size"])

    #print("Starting to build phase blocks...")

    for index in dmrs.index:

        #print("Currently working on: ", dmrs.loc[index, "name"])

        dmr_start = dmrs.loc[index, "start"]
        dmr_end = dmrs.loc[index, "end"]

        #Checking that a read covers our ICR:
        if find_best_read(dmr_end, "up", reads, base_error=error_rate) == "no_read":
            #print("This ICR not covered by an read. Moving on to the next. Recommend restarting simulation.")
            #warnings.warn("This ICR not covered by an read. Moving on to the next. Recommend restarting simulation.")
            continue

        elif find_best_read(dmr_start, "down", reads, base_error=error_rate) == "no_read":
            #print("This ICR not covered by an read. Moving on to the next. Recommend restarting simulation.")
            #warnings.warn("This ICR not covered by an read. Moving on to the next. Recommend restarting simulation.")
            continue       
    
        else:
            #upper_limit = reads.loc[find_best_read(dmr_end, "up", reads, base_error=error_rate), "from"]
            #lower_limit = reads.loc[find_best_read(dmr_start, "down", reads, base_error=error_rate), "to"]
            upper_limit = find_best_read(dmr_end, "up", reads, base_error=error_rate)
            lower_limit = find_best_read(dmr_start, "down", reads, base_error=error_rate)
    
            topmost_variant = topmost_variant_in_range (variants, upper_limit)
            bottommost_variant = bottommost_variant_in_range (variants, lower_limit)
    
            no_upper = False
            if (topmost_variant >= dmr_start):
                no_upper = True
                #print("No variants in upper range for this ICR.")
                #warnings.warn("No variants were detected in the upper range. If you see this warning too much, something is wrong.")
    
            else:
                #print("Working on determining upper range...")
                current_pos = 0
                while ((upper_limit > 0) & (current_pos != topmost_variant)):
                    current_pos = topmost_variant
                    upper_limit = find_best_read(topmost_variant, "up", reads, base_error=error_rate)
                    topmost_variant = topmost_variant_in_range(variants, upper_limit)
    
                #print("Determined interval upper limit: ", topmost_variant)
    
    
            no_lower = False        
            #checking that we do find a variant after the DMR within range
            if (bottommost_variant <= dmr_end):
                no_lower = True
                #print("No variants detected in the lower range for this ICR.")
                #warnings.warn("No variants were detected in the lower range. If you see this warning too much, something is wrong.")
    
            else:
                #print("Working on determining lower range...")
                last_variant = variants[-1]
                #print("The last variant's position is: ", last_variant) 
                current_pos = last_variant
            
                while ((lower_limit < last_variant) & (current_pos != bottommost_variant)):
                    current_pos = bottommost_variant
                    lower_limit = find_best_read(bottommost_variant, "down", reads, base_error=error_rate)
                    bottommost_variant = bottommost_variant_in_range(variants, lower_limit)
    
                #print("Determined interval lower limit: ", bottommost_variant)


            if no_upper == True:
                phase_set_start = int(dmr_start)
    
            else:
                phase_set_start = int(topmost_variant)
    
    
            if no_lower == True:
                phase_set_end = int(dmr_end)
                
            else:
                phase_set_end = int(bottommost_variant)
    
    
            block_size = phase_set_end - phase_set_start 
            #print("This DMR can phase a region of: ", block_size)
    
            blocks.loc[index, "start"] = int(phase_set_start)
            blocks.loc[index, "end"] = int(phase_set_end)
            blocks.loc[index, "size"] = int(block_size)

            blocks.drop_duplicates(subset="start", keep="last", inplace=True)
            blocks.drop_duplicates(subset="end", keep="first", inplace=True)
    
    return blocks


def create_variant_distribution (chrom='test', show_plot=False, variant_per=1000, var_mu=6.2, var_sigma=1.4, lcr_exclusion=False):
    #Algo is simple: we know the number of variants in the whole genome, and we know the size of each chrom
    #We use the above to generate the number of variants for given chromosome (say: n)
    #We use log-normal distribution to generate positions of the n variants lying between 1 and chrom_size
    #For each variant, we substitute it for the closest number in the high complexity regions (HCRs).
    #This substitution does not change the underlying distance distribution.

    variant_pos = []

    size = test_size
    num_variants = int(size/variant_per) #assuming 1 variant every 1 kilobases
    #print("The number of variants to be generated here are: ", num_variants)

    ## The following section simulates inter-variant distances based on log-normal distribution
    variant_distances = np.random.lognormal(mean=var_mu, sigma=var_sigma, size=num_variants-1)
    
    #Initialising first variant
    variant_pos = [np.float64(10001)]

    #For each successive variant, we add the distance to the
    i = 0
    for dist in variant_distances:
        variant_pos.append(variant_pos[i] + int(dist))
        i += 1
    
    #print("Elements in variant_pos: ", len(variant_pos))
    actual_variant_pos = np.unique(variant_pos)
    
    #print("Total number of variants expected: ", num_variants)
    actual_variant_pos = [pos for pos in variant_pos if pos < size]
    #print("Actual number of variants: ", len(actual_variant_pos))
    

    if show_plot==True:
        plt.hist(actual_variant_pos, bins=10000)
        plt.suptitle("Simulated variant locations in: "+chrom)
        plt.title("Each bin represents 10 kilobases (kb)", fontsize=10)
        plt.ylabel("Frequency")
        plt.xlabel("Genomic Position (100 Megabases)")
        plt.show()
        #plt.savefig("../results/chr8_variant_distribution.jpg", dpi=1200)

        plt.hist(sorted(variant_distances), bins = 100)
        plt.suptitle("Intervariant distance distribution for: "+chrom)
        plt.title("Log mean distance =  "+str(var_mu)+", SD = "+str(var_sigma))
        plt.xlabel("Distance (bases)")
        plt.ylabel("Frequency")
        plt.show()

        plt.hist(np.log(sorted(variant_distances)), bins = 100)
        plt.suptitle("Log distance distribution for: "+chrom)
        plt.title("Log mean distance =  "+str(var_mu)+", SD = "+str(var_sigma))
        plt.xlabel("Distance (bases)")
        plt.ylabel("Frequency")
        plt.show()
    #print("Generated variants for: ", chrom)

    #actual_distances = [dist for dist in actual_distances if dist != 0]
    #return actual_distances
    return sorted(actual_variant_pos)

def topmost_variant_in_range (vcf, upper):
    #Checking if there are any variants in the upper range
    if len(np.argwhere(vcf > upper)) == 0:
        return upper
    else:
        j = np.argwhere(vcf > upper).min()
        return vcf[j]

def bottommost_variant_in_range (vcf, lower):
    if len(np.argwhere(vcf < lower)) == 0:
        return lower
    else:
        k = np.argwhere(vcf < lower).max()
        return vcf[k]


def all_sim_results (depths=[35], icrs=test_icrs, chrom="test", n_sim=20, simulate_variants=True, input_pop_vcf=None, 
                     variant_per=1000, filepath=None, output_blocks=False, overwrite=True, 
                     var_mu=6.2, var_sigma=1.4, read_mu=8.62, read_sigma=[1.11], error_rate=0.10):
#This is just like plotting_results, except it returns the entire resultant dataset, which can then be plotted with confidence intervals
    
    if chrom != "all":
        
        #We need to keep track of whether this is the first length-depth pair for that chromosome or not
        simulation_counter = 0
        
        #Initialising the dataframe containing pofo-phasing results for a certain length-depth pair
        phasing_results = pd.DataFrame(columns=["read_sigma", "depth", "error_rate", "mean_pct_phasing"])
        index=0
        
        for depth in depths:
            
            for sigma in read_sigma:

                #Initialising list containing dataframes of phase blocks
                phase_block_list = []

               
                
                print("Current depth: ", depth)
                print("Read distribution log-mean: ", read_mu, ", SD: ", sigma) 
                #print("Given read lengths: ", read_lengths)

                #Now, we run the simulation a pre-defined number of times, and append each result to the results table
                i = 0
                while i < n_sim:
                    print("Working on chromosome: ", chrom)
                    print("At simulation number: ", i+1)
                    #Running simulation, assigning total fraction to "random_result", and phase blocks to "blocks"
                    random_result, blocks = phasing_simulator(icrs=test_icrs, mean_depth=depth, chrom=chrom, 
                                                     simulate_variants=simulate_variants, pop_vcf=input_pop_vcf, variant_per=variant_per, 
                                                             var_mu = var_mu, var_sigma=var_sigma, read_mu=read_mu, read_sigma=sigma, 
                                                             error_rate=error_rate)

                    #Adding read_length and depth to the blocks dataframe
                    blocks["read_sigma"] = sigma
                    blocks["read_depth"] = depth
                    
                    #Appending fraction to the phasing_results dataframe
                    phasing_results.loc[index] = [sigma, depth, error_rate, random_result]

                    #Appending blocks to a temporary list
                    phase_block_list.append(blocks)
                    index += 1
                    i += 1
                    
                    '''
                    answer = input("Continue? (y/n)")
                    if answer == "y":
                        continue
                    else:
                        i = n_sim
                    '''

                    

                print("Done with the simulations for the parameters; chrom: ", chrom, " and depth: ", depth)

                if output_blocks==False:
                    data = phasing_results

                else:
                    data = pd.concat(phase_block_list)
                
                #If the code is just supposed to add results to a preexisiting file, we do that
                if overwrite ==  False:
                    data.to_csv(filepath, mode='a', index=False, header=False)

                #Else, we create a new file for the first simulation, then add the next results to the same file
                else:
                    #If this is the first simulation, we add the header and overwrite the old file
                    if simulation_counter == 0:
                        data.to_csv(filepath, index=False, header=True)
    
                    #For all following simulations, we add the results to the file just created for the first simulation
                    else:
                        data.to_csv(filepath, mode='a', index=False, header=False)               
                
                
                simulation_counter += 1

        #Concatenating all phase blocks in the list into one dataframe
        phase_blocks = pd.concat(phase_block_list)
        print(phase_blocks)
        #return phase_blocks

        #To get mean % phasing
        print("Mean phasing %age: ", np.mean(phasing_results["mean_pct_phasing"]))
        print("Median phasing %age: ", np.median(phasing_results["mean_pct_phasing"]))
        print("Minimum phasing %age: ", np.min(phasing_results["mean_pct_phasing"]))
        print("Maximum phasing %age: ", np.max(phasing_results["mean_pct_phasing"]))
        return data
        
        
    else:
        print("not yet!")
        
def test_sim_results (depth=35, icrs=test_icrs, chrom="test", n_sim=20, simulate_variants=True, input_pop_vcf=None, 
                     variant_per=1000, filepath=None, output_blocks=False, overwrite=True, 
                     var_mu=6.2, var_sigma=1.4, read_mu=8.62, read_sigma=1.11, error_rate=0.10):

    #Initialising the dataframe containing pofo-phasing results for a certain length-depth pair
    phasing_results = pd.DataFrame(columns=["read_sigma", "depth", "error_rate", "mean_pct_phasing"])
    index=0

    #Initialising list containing dataframes of phase blocks
    phase_block_list = []

    #print("Current depth: ", depth)
    #print("Read distribution log-mean: ", read_mu, ", SD: ", read_sigma) 
    #print("Given read lengths: ", read_lengths)

    #Now, we run the simulation a pre-defined number of times, and append each result to the results table
    print("Simulating phase blocks...")
    b = progressbar.ProgressBar(maxval=n_sim)
    b.start()
    i = 0
    sim_duration = 1
    while i < n_sim:
        due_remaining_min = sim_duration * (n_sim - i) / 60
        start_time = time.time()
        #print("Estimated time remaining (minutes): %d" % due_remaining_min, end='\r')
        b.update(i+1)
        
        #print("Working on chromosome: ", chrom)
        #print("At simulation number: ", i+1)
        #Running simulation, assigning total fraction to "random_result", and phase blocks to "blocks"
        random_result, blocks = phasing_simulator(icrs=test_icrs, mean_depth=depth, chrom=chrom,
                                                  simulate_variants=simulate_variants, pop_vcf=input_pop_vcf, variant_per=variant_per, 
                                                  var_mu = var_mu, var_sigma=var_sigma, read_mu=read_mu, read_sigma=read_sigma, 
                                                  error_rate=error_rate)

        #Adding read_length and depth to the blocks dataframe
        blocks["read_sigma"] = read_sigma
        blocks["read_depth"] = depth
        
        #Appending fraction to the phasing_results dataframe
        phasing_results.loc[index] = [read_sigma, depth, error_rate, random_result]

        #Appending blocks to a temporary list
        phase_block_list.append(blocks)
        index += 1

        sim_duration = time.time() - start_time  
        
        i += 1
        
    #print("Done with the simulations for the parameters; chrom: ", chrom, " and depth: ", depth)

    if output_blocks==False:
        data = phasing_results

    else:
        data = pd.concat(phase_block_list)
    
    
    return data
    
    
def non_ideal_chrom(chrom, chrom_dmrs=chrom_dmrs, bases_phased=15026166):
    #bases_phased = function(read_sigma, depth)

    #Here, looking at the mean:
    #bases_phased = 15585014
    
    half_length = np.int64(bases_phased/2)

    # Our DMRs:
    our_dmrs = chrom_dmrs[chrom_dmrs["chrom"]==chrom]
    our_dmrs.reset_index(inplace=True)
    our_dmrs = our_dmrs.iloc[:, 1:]

    # Now we make the prediction blocks
    predicted_blocks = pd.DataFrame(columns=["chrom", "dmr_name", "left_end", "right_end", "size"])

    centstart = our_dmrs.loc[0, "centstart"]
    centend = our_dmrs.loc[0, "centend"]
    length = our_dmrs.loc[0, "length"]
    
    i = 0
    while i < len(our_dmrs):
        # We need to do a bunch of things.
        # 1) We need to calculate the phase block for this DMR if it was in an ideal chromosome.
        # 2) We need to trim the phase block.
        # 3) Trim limits will be imposed by the chromosome size (0, length) and centromere (centstart, centend)
        # 4) This trimming depends on whether the DMR is on the left or right arm of the chromosome.

        dmr_start = our_dmrs.loc[i, "start"]
        dmr_end = our_dmrs.loc[i, "end"]
        
        untrimmed_left = our_dmrs.loc[i, "end"] - half_length
        untrimmed_right = our_dmrs.loc[i, "start"] + half_length

        # Setting final values as untrimmed values. Will change them during the trimming process.
        left_limit = untrimmed_left
        right_limit = untrimmed_right

        # Determining arm
        if dmr_start > centend:
            arm = 'right'

        elif dmr_end < centstart:
            arm = 'left'

        else:
            warnings.warn("Can't figure out which arm the DMR lies on!")
            continue


        # Trimming based on arm
        if arm == 'left':
            if untrimmed_left < 1:
                left_limit = 1

            if untrimmed_right > centstart:
                right_limit = centstart

        elif arm == 'right':
            if untrimmed_left < centend:
                left_limit = centend

            if untrimmed_right > length:
                right_limit = length

        else:
            warnings.warn("No arm detected! What is this, Spain in the 50s?")
        
        size = right_limit - left_limit
        
        predicted_blocks.loc[len(predicted_blocks)] = [chrom, our_dmrs.loc[i, "name"], left_limit, right_limit, size]

        i += 1

    total_phased_bases = remove_overlaps(predicted_blocks)
    print("Phasing %age: ", total_phased_bases*100/length)
    
    return predicted_blocks

def remove_overlaps (blocks):
    
    range_df = blocks.copy()
    #print(range_df)
    indices = list(range_df.index)
    #print(indices)
    encapsulated_indices = []
    
    i = 1
    while i < len(indices):
        if range_df.loc[indices[i], "left_end"] <= range_df.loc[indices[i-1], "right_end"]:
            #This means our second range is overlapping with the first

            if range_df.loc[indices[i], "right_end"] <= range_df.loc[indices[i-1], "right_end"]:
                #This means the second range is completely encapsulated by the first
                #print(indices[i])
                #print(encapsulated_indices)
                encapsulated_indices.append(i)
                
            else:
                range_df.loc[indices[i], "left_end"] = range_df.loc[indices[i-1], "right_end"] + 1
                

        i += 1

    range_df.drop(encapsulated_indices, axis=0, inplace=True)

    range_df["size"] = range_df["right_end"]-range_df["left_end"]

    total_bases_phased = np.sum(range_df["size"])    
    return total_bases_phased
