import ideal_pofo_simulator as ideal
import pandas as pd
import numpy as np

def phased_bases (test_sim_results):

    results = test_sim_results
    
    frac = results["mean_pct_phasing"].mean()
    if frac > 99:
        ideal.warnings.warn("Max phasing reached. Needs tweaking in the ideal chromosome size.")

    bases_phased = np.round(frac * 1000000)

    return np.int64(bases_phased)

def chrom_pofo_predictor (chrom='all', bases_phased=None, read_mu=8.62, read_sigma=1.11, 
                          var_mu=6.2, var_sigma=1.4, depth=35, 
                          error_rate=0.1, variant_per=1000, n_sim=None):

    start_time = ideal.time.time()
    

    if bases_phased == None:
        print("How accurate do you want the predictions to be? Better accuracy means longer wait times!")
        code = input("Enter code (1-5): \n1. Low (nsim=2) \n2. Medium (nsim=10) \n3. High (nsim=100) \n4. Super (nsim=1000) \n5. Custom\n")
        nsim_codes = {1:2, 2:10, 3:100, 4:1000}
        code = int(code)
        
        while (code not in [1,2,3,4,5]):
            print("Sorry, could not understand. Accepted inputs are: 1, 2, 3, 4 and 5.")
            code = input("Enter code (1-5): \n1. Low (nsim=2) \n2. Medium (nsim=10) \n3. High (nsim=100) \n4. Super (nsim=1000) \n5. Custom\n")
            code = int(code)
            
        if code == 5:
            n_sim = int(input("Enter custom n_sim: "))

        else:
            n_sim = nsim_codes[code]
        
        results = ideal.test_sim_results(read_mu=read_mu, read_sigma=read_sigma, var_mu=var_mu, 
                                         var_sigma=var_sigma, error_rate=error_rate, 
                                         variant_per=variant_per, n_sim=n_sim)
        
        bases_phased = phased_bases(results)
        print("\nIn an ideal chromosome with these params, we phase: ", round(bases_phased), " bases.")

    else:
        bases_phased =  bases_phased

        print("We are given that in an ideal chromosome, we phase: ", round(bases_phased), " bases.")

    if chrom=="all":
        for chrom in ideal.chrom_dmrs["chrom"].unique():
            print("\n\nChromosome: ", chrom)
            print(ideal.non_ideal_chrom(chrom=chrom, bases_phased=bases_phased))

    else:
        print(ideal.non_ideal_chrom(chrom=chrom, bases_phased=bases_phased))

    print("\n\nTime taken (seconds): ----%s----" % (ideal.time.time()-start_time))

def take_inputs():
    print("Hello! Welcome to the PofO Phasing Predictor.")
    bases_phased = input("If you already know how many bases you're phasing (e.g. for re-checking results), enter: \n")

    if bases_phased != "":
        bases_phased = int(bases_phased)
        chrom_pofo_predictor(bases_phased=bases_phased)

    else:
        print("We will have to find how many bases we phase with our parameters in an ideal chromosome.")
        print("WARNING: If you leave an input blank, the default will be used.")

        var_mu = float(input("Variant Distribution - Log-mean (default 6.2): "))
        var_sigma = float(input("Variant Distribution - Log-SD (default 1.4): "))
        variant_per = float(input("You have one variant per ---- bases (default 1,000): "))

        print("Moving on to sequencing parameters...")

        depth = float(input("Sequencing depth of coverage: "))

        calc_log = input("Do you have log-parameters of read length distribution (y/n)? If not, you must enter mean and median read lengths. ")

        if calc_log == 'y':
            read_mu = float(input("Read Distribution - Log-mean (default 8.62): "))
            read_sigma = float(input("Read Distribution - Log-SD (default 1.11): "))

        elif calc_log == 'n':
            mean_read = float(input("Read Distribution - Mean read length (default 10,260 bp): \n"))
            median_read = float(input("Read Distribution - Median read length (default 5,541 bp): \n"))

            read_mu = np.log(median_read)
            read_sigma = np.sqrt(2*(np.log(mean_read) - np.log(median_read)))

        error_rate = float(input("Basecalling error rate (between 0 to 1): "))

        chrom_pofo_predictor (read_mu=read_mu, read_sigma=read_sigma, var_mu=var_mu, 
                              var_sigma=var_sigma, depth=depth, 
                              error_rate=error_rate, variant_per=variant_per, n_sim=None)

take_inputs()
        
    


