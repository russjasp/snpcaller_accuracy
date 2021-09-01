

### Author Contact Information

Russ Jasper
rjjasper@ucalgary.ca

### Usage and license information

If you use or are inspired by code in this repository please cite the following work or contact me about how to cite. 

<link to come here> #### FINISH THIS ####

---

# Evaluating the accuracy of variant calling methods using the frequency of parent-offspring genotype mismatch

After calling snps with FreeBayes, HaplotypeCaller, SAMtools, UnifiedGenotyper, and VarScan and performing a base level of filtering, the resulting snp sets varied over orders of magnitude in the number of sites called, and in the error rates.

I used R to iteratively filter these 5 data sets to the same number of snps called. Upon filtering each data set to the same number of snps I could compare the error rates between variant callers in a standardized manner.

---

# Environment

R version 1.4.1106


# Usage

### snpcaller_accuracy.R

Call snps with the variant caller of your choice, convert  vcf file to a table (eg, VariantsToTable), and read into R (object "snp_df").

Choose the different filtering metrics you would like to use to filter (eg, I ended up using QUAL, DP, GQ). Metrics by site: "SITE_tot_filters" and metrics by genotype: "GT_filters".

Choose the different values for the number of snps you would like to filter all variant callers to (vector "target_num_snps").

Code will loop through each target number of snps and iteratively filter the snp set down using the filtering metrics of your choice until it gets as close as possible to the target number of snps.

Code will then loop through each target number of snps (and the resulting filtered data set for each target number of snps) and calculate the results. Results include, number of snps/genotypes remaining, error rates (by site and by genotype), filtering quantile and value for each filtering metric.

Code outputs a lot of messages as it goes which I was just using to error check things as I went, can just ignore the verbosity.


# Repo Table of Contents

## Code

### snpcaller_accuracy.R
As described above

## Data

### Example_Data.table
An example snp set. Use this as the input data (object "snp_df").

### Example_Results.csv
I have filtered the Example_Data.table using QUAL, DP, GQ.
I used the target snp values: 5000, 2500, 1000.

Example_Results.csv contains the results at each level of filtering. Column names are as described:  
Sites - number of sites remaining after filtering  
GT - number of genotype calls remaining after filtering  
NA.gt - number of NA genotype calls remaining after filtering  
Match.abs.site (gt) - absolute number of sites (genotypes) matching between P1 & F2  
Miss.abs.site (gt) - absolute number of sites (genotypes) mismatching between P1 & F2  
Match.prop.site (gt) - proportion of sites (genotypes) matching between P1 & F2  
Miss.prop.site (gt) - proportion of sites (genotypes) mismatching between P1 & F2  
quantile - quantile of the distribution of each metric the data set was filtered at  
QUAL - value of QUAL corresponding to the given quantile of the QUAL distribution in F2  
P1.QUAL - value of QUAL corresponding to the given quantile of the QUAL distribution in P1  



