
## Evaluating the accuracy of variant calling methods using the frequency of parent-offspring genotype mismatch

## Code to iteratively filter down snp sets from different variant caller programs
## At each level of filtering, calculate:
  # snps & genotypes remaining
  # error rate by site & by genotype
  # value of each filtering metric

## Russ Jasper
## rjjasper@ucalgary.ca



library(data.table)
options(scipen=999)

######################################################################################################

## Work around because UG data has different distributions for by site metrics for P1 vs F2
## ie, for a given site UG has a QUAL value for the P1 and a QUAL value for the F2
## this is in contrast to data from FB, HC, SAM, etc where a given site has one single QUAL value for the entire cohort

## if UG data (or data with different by site metrics between P1 & F2)
UGHACK <- TRUE
## if normal data (a single value for the entire site, not different between P1 & F2)
UGHACK <- FALSE


######################################################################################################
### FUNCTIONS ###

## find the value of a metric that corresponds to the quantile value
## ie, A percentile of 10% in HC sample data gives a QUAL of 70
find.quant.val <- function(Filters,DF,quant,filter.type){
  df.var0 <- DF[grepl(Filters, names(DF))]
  if(filter.type=="SITE"){df.var <- df.var0[!grepl("LP",names(df.var0))]}
  if(filter.type=="GT"){df.var <- df.var0[grepl("LP",names(df.var0))]}
  quant.val <- quantile(df.var, quant, na.rm = TRUE)
  return(quant.val)
}


## find the labels of all sites to filter out with - by SITE - metrics
## (by SITE metrics = metrics with a value for the entire site/snp, eg, QUAL, FS, MQ. In contrast to by GT metrics, eg, GQ, DP)
site.filter.fun <- function(site.filter.list, quant.vals, HACK){
  
  metric_var <- site.filter.list
  
  filter_var <- as.numeric(unname(quant.vals[grepl(metric_var, names(quant.vals))]))
  
  snp_df_var0 <- snp_df_no_P1_down_filters_REVERSE[grepl(paste0("label|",metric_var), names(snp_df_no_P1_down_filters_REVERSE))]
  
  snp_df_var <- snp_df_var0[!grepl("LP",names(snp_df_var0))]
  
  filtered_df <- snp_df_var[snp_df_var[,metric_var] < filter_var | is.na(snp_df_var[,metric_var]),]
  
  
  ## have to reverse engineer the P1 by site metric
  ## (because I didn't realize that UG would have different by site metrics for P1 and for F2 when I originally wrote this)
  if(paste0(P1name,".",metric_var) %in% colnames(snp_df_P1)){
    
    if(HACK){
      quantile_var2 <- gsub("%","",names(quant.vals[grepl(metric_var, names(quant.vals))]))
      quantile_var1 <- gsub(paste0(metric_var,"."),"",quantile_var2)
      quantile_var <- as.numeric(quantile_var1)/100
      prop_offspring <- sum(snp_df_no_P1_down_filters_REVERSE$QUAL < filter_var, na.rm=TRUE) / dim(snp_df_no_P1_down_filters_REVERSE)[1]
      if(quantile_var*1.3 < prop_offspring | quantile_var/1.3 > prop_offspring){message("ERROR UG");break}
    }
    if(HACK==FALSE){
      quantile_var <- quant.vals[,"quantile"]
    }
    
    message("P1 quantile: ",quantile_var)
    
    P1_df0 <- snp_df_P1[grepl(paste0(metric_var,"|label"), names(snp_df_P1))]
    
    P1_df <- P1_df0[!grepl("label", names(P1_df0))]
    
    P1_quant <- quantile(P1_df, quantile_var, na.rm = TRUE)
    
    message("P1 value: ",P1_quant)
    
    P1_Filter <- P1_df0[P1_df < P1_quant,"label"]

    sites_to_remove_list <- unique(c(filtered_df$label,P1_Filter))
       
     
  } else { sites_to_remove_list <- filtered_df$label }
  
  
  cat(length(sites_to_remove_list), "absolute sites filtered out -",metric_var,"\n")
  cat(round(length(sites_to_remove_list)/dim(snp_df_var)[1],3), "proportion sites filtered out\n")
  
  return(sites_to_remove_list)
}




## aggregate the values of all filter metrics for all quantile values
## count number of snps that remain at each quantile filtering level
## ie, 10% quantile level results in 600k snps remaining in HC sample data
quantile.matrix.function <- function(sequence.list){
  
  ## DOWN filter metrics are reported with OPPOSITE sign ##
  ## (DOWN filter metrics, eg FS, SOR, metrics where decreasing the value of the metric increases the stringency)
  snp_results <- array(NA, dim = c(1,tot_num_filters+2), dimnames = list(c(),c("quantile",SITE_tot_filters,GT_filters,"Sites")))
  
  
  snp_results[1,1] <- sequence.list
  
  
  sites_quant_vals <- sapply(SITE_tot_filters, FUN = find.quant.val, DF = snp_df_no_P1_down_filters_REVERSE, quant = sequence.list, filter.type = "SITE")
  
  snp_results[1,c(2:c(length(SITE_tot_filters)+1))] <- sites_quant_vals
  
  if(identical(paste0(names(snp_results[1,c(2:c(length(SITE_tot_filters)+1))]),".",sequence.list*100,"%"),names(sites_quant_vals)) == FALSE){message("error -- fill array sites");break}
  
  
  gt_quant_vals <- sapply(GT_filters, FUN = find.quant.val, DF = snp_df_no_P1_down_filters_REVERSE, quant = sequence.list, filter.type = "GT")
  
  snp_results[1,c(c(length(SITE_tot_filters)+2):c(length(GT_filters)+length(SITE_tot_filters)+1))] <- gt_quant_vals
  
  if(identical(paste0(names(snp_results[1,c(c(length(SITE_tot_filters)+2):c(length(GT_filters)+length(SITE_tot_filters)+1))]),".",sequence.list*100,"%"),names(gt_quant_vals)) == FALSE){message("error -- fill array gt");break}
  
  
  
  ## keep a list of the sites to filter out from each metric
  ## then take unique of list from each metric
  
  message(paste0("\nQUANTILE: ", sequence.list))
  
  sites_to_filter_out_full_list <- lapply(SITE_tot_filters, FUN = site.filter.fun, quant.vals = sites_quant_vals, HACK=UGHACK)
  names(sites_to_filter_out_full_list) <- SITE_tot_filters
  
  unique_sites_to_remove_by_SITE <- NULL
  
  for(l in 1: length(sites_to_filter_out_full_list)){
    
    unique_sites_to_remove_by_SITE <- unique(c(unique_sites_to_remove_by_SITE, sites_to_filter_out_full_list[[l]]))
    
  }
  
  length(unique_sites_to_remove_by_SITE)
  
  
  
  
  ## filter GQ and then DP iteratively
  
  GQ_filt1 <- as.numeric(unname(gt_quant_vals[grepl("GQ", names(gt_quant_vals))]))
  DP_filt1 <- as.numeric(unname(gt_quant_vals[grepl("DP", names(gt_quant_vals))]))
  
  
  snp_df_GT <- snp_df_no_P1_down_filters_REVERSE[grepl("GT|label", names(snp_df_no_P1_down_filters_REVERSE))]
  rownames(snp_df_GT) <- snp_df_GT[,"label"]
  
  
  ## filter GQ
  snp_df_GQ <- snp_df_no_P1_down_filters_REVERSE[grepl("GQ|label", names(snp_df_no_P1_down_filters_REVERSE))]
  rownames(snp_df_GQ) <- snp_df_GQ[,"label"]
  
  snp_df_GQ_filt <- snp_df_GQ[grepl("GQ", names(snp_df_GQ))] < GQ_filt1
  
  snp_df_GT_filt1 <- snp_df_GT[grepl("GT", names(snp_df_GT))]
  pre_GT_missing <- rowSums(snp_df_GT_filt1 == "./." | snp_df_GT_filt1 == "." | snp_df_GT_filt1 == ".|.", na.rm=TRUE)
  cat("GT missing PRE filtering:",sum(pre_GT_missing),"\n")
  cat("Proportion:",sum(pre_GT_missing)/(dim(snp_df_GT_filt1)[1]*dim(snp_df_GT_filt1)[2]),"\n")
  
  snp_df_GT_filt1[snp_df_GQ_filt] <- "./."
  
  GQ_filt_GT_missing <- rowSums(snp_df_GT_filt1 == "./." | snp_df_GT_filt1 == "." | snp_df_GT_filt1 == ".|.", na.rm=TRUE)
  cat("GT missing GQ filtering:",sum(GQ_filt_GT_missing),"\n")
  cat("Proportion:",sum(GQ_filt_GT_missing)/(dim(snp_df_GT_filt1)[1]*dim(snp_df_GT_filt1)[2]),"\n")
  
  
  ## filter DP
  snp_df_DP <- snp_df_no_P1_down_filters_REVERSE[grepl("DP|label", names(snp_df_no_P1_down_filters_REVERSE))]
  rownames(snp_df_DP) <- snp_df_DP[,"label"]
  
  snp_df_DP_filt <- snp_df_DP[grepl("DP", names(snp_df_DP))] < DP_filt1
  
  snp_df_GT_filt2 <- snp_df_GT_filt1
  
  snp_df_GT_filt2[snp_df_DP_filt] <- "./."
  
  DP_filt_GT_missing <- rowSums(snp_df_GT_filt2 == "./." | snp_df_GT_filt2 == "." | snp_df_GT_filt2 == ".|.", na.rm=TRUE)
  cat("GT missing DP filtering:",sum(DP_filt_GT_missing),"\n")
  cat("Proportion:",sum(DP_filt_GT_missing)/(dim(snp_df_GT_filt1)[1]*dim(snp_df_GT_filt1)[2]),"\n")
  
  
  snp_df_GT_filt3 <- snp_df_GT_filt2
  FINAL_filt_GT_missing <- rowSums(snp_df_GT_filt3 == "./." | snp_df_GT_filt3 == "." | snp_df_GT_filt3 == ".|.", na.rm=TRUE)
  
  snp_df_GT_filt3$label <- rownames(snp_df_GT_filt3)
  
  unique_sites_to_remove_by_GT <- snp_df_GT_filt3[FINAL_filt_GT_missing > 53, "label"]
  
  cat("Sites removed by GT filters:", length(unique_sites_to_remove_by_GT),"\n")
  
  
  FULL_sites_to_remove0 <- unique(c(unique_sites_to_remove_by_GT,unique_sites_to_remove_by_SITE))
  length(FULL_sites_to_remove0)
  
  
  if(identical(rownames(snp_df_GT),rownames(snp_df_GQ))==FALSE){message("ERROR A");break}
  if(identical(rownames(snp_df_GT),rownames(snp_df_DP))==FALSE){message("ERROR B");break}
  if(identical(rownames(snp_df_GQ),rownames(snp_df_DP))==FALSE){message("ERROR C");break}
  
  
  ## also remove sites where Parent is below GQ/DP threshold
  
  P1_GQ_filt1 <- quantile(snp_df_P1[grepl("GQ", names(snp_df_P1))], sequence.list, na.rm=TRUE)
  P1_DP_filt1 <- quantile(snp_df_P1[grepl("DP", names(snp_df_P1))], sequence.list, na.rm=TRUE)
  
  filt_P1_GQ <- snp_df_P1[snp_df_P1[grepl("GQ", names(snp_df_P1))] < P1_GQ_filt1,]
  
  filt_P1_DP <- snp_df_P1[snp_df_P1[grepl("DP", names(snp_df_P1))] < P1_DP_filt1,]
  
  FULL_sites_to_remove <- unique(c(FULL_sites_to_remove0, filt_P1_GQ$label, filt_P1_DP$label))
  
  length(FULL_sites_to_remove)
  
  
  
  num_snps_remain <- dim(snp_df_no_P1_orig)[1]-length(FULL_sites_to_remove)
  
  snp_results[1,c(2+tot_num_filters)] <- num_snps_remain
  
  
  return(snp_results)
  
  # if(!is.na(snps.tar) & num_snps_remain <= snps.tar){message("finished quantile matrix early!");break}
  
  
} ## function end


## finds a new, finer-scale sequence of quantiles to filter by
find.new.quantile.seq <- function(target,DF){
  
  all_quants_more_than_target_snps <- DF[DF[,"Sites"] >= target,]
  lower_range_quant <- all_quants_more_than_target_snps[all_quants_more_than_target_snps[,"Sites"] == min(all_quants_more_than_target_snps[,"Sites"]),]
  
  all_quants_less_than_target_snps <- DF[DF[,"Sites"] <= target,]
  upper_range_quant <- all_quants_less_than_target_snps[all_quants_less_than_target_snps[,"Sites"] == max(all_quants_less_than_target_snps[,"Sites"]),]
  
  finer_scale_range <- c(lower_range_quant[,"quantile"],upper_range_quant[,"quantile"])
  
  pre_new_seq <- seq(finer_scale_range[1]*100, finer_scale_range[2]*100, length.out = 11)/100
  
  new_seq <- pre_new_seq[-c(1,11)]
  
  results_list <- vector(mode = "list", length = 3)
  results_list[[1]] <- new_seq
  results_list[[2]] <- lower_range_quant
  results_list[[3]] <- upper_range_quant
  
  return(results_list)
  
}


### End of all the functions



## total SNPs in each dataset
## input total SNPs so I can compare the error rates at each value here
UG_tot_snp  <- 3489299
HC_tot_snp <- 2238403
SAM_tot_snp <- 459477
FB_tot_snp <- 119318
VAR_tot_snp <- 87868

All_tot_snp <- c(UG_tot_snp,HC_tot_snp,SAM_tot_snp,FB_tot_snp,VAR_tot_snp)



## runs one snp set at a time
setwd("")
snp_df <- data.frame(fread("Example_Data.table", header = TRUE, stringsAsFactors = FALSE))


## input snp caller name
SNPCALLER_NAME <- "Example"

## input parent name in the dataset
## make sure parent name is distinct from offspring columns
P1name <- "LP.mg0P1"
cat("does P1 name exist:",sum(grepl(P1name,colnames(snp_df)))>1) # just making sure


snp_df_P1 <- snp_df[grepl(paste0(P1name,"|label"), names(snp_df))]

snp_df_no_P1_orig <- snp_df[!grepl(P1name, names(snp_df))]


## choose which filters to use
SITE_UP_filters <- c("QUAL") ## metrics where INCREASING the value of the metric increases stringency
SITE_DOWN_filters <- NULL ## metric where DECREASING the value of the metric increases stringency
SITE_tot_filters <- c(SITE_UP_filters, SITE_DOWN_filters)
GT_filters <- c("GQ","DP")
tot_num_filters <- sum(length(SITE_tot_filters),length(GT_filters))


## FLIP THE ORIENTATION OF THE DISTRIBUTION FOR DOWN FILTERS ##
# ie, filter FS downwards so flip the distribution and then filter it upwards in SAME manner as the up metrics

snp_df_no_P1_down_filters_REVERSE <- snp_df_no_P1_orig
for(dd in SITE_DOWN_filters){ snp_df_no_P1_down_filters_REVERSE[,dd] <- snp_df_no_P1_down_filters_REVERSE[,dd] * -1 }



### VARSCAN: add dummy variable if you don't want to use a by site metric ###
# snp_df_no_P1_down_filters_REVERSE$varscan <- 1
###


## TEST SEQUENCE START ##
# no need to change
# this just filters from 0 to 100% in 10% increments as a starting point
test_seq <- seq(0,100, by=10)/100


######################################################################################################

# FIRST ITERATION #

## remember-- DOWN filter metrics are reported with OPPOSITE sign
snp_results_mat <- lapply(test_seq, FUN = quantile.matrix.function)
snp_results_DF1 <- data.frame(matrix(unlist(snp_results_mat), nrow=length(snp_results_mat), byrow=T, dimnames=list(NULL, colnames(snp_results_mat[[1]]))),stringsAsFactors=FALSE)

Total_num_snps <- dim(snp_df)[1]


# Choose Target Number of SNPs #
# target snp == 0   means the max number of snps after filtering out NA's at all by site metrics



## INPUT THE DIFFERENT VALUES FOR THE NUMBERS OF SNPS YOU WANT TO FILTER TO ##
target_num_snps_standard <- c(5e3, 2.5e3, 1e3) 
## INPUT THE DIFFERENT VALUES FOR THE NUMBERS OF SNPS YOU WANT TO FILTER TO ##



## This vector contains all of the "target snp" values less than or equal to the total snps in the current dataset
target_num_snps0 <- sort(c(All_tot_snp,target_num_snps_standard), decreasing = TRUE)
target_num_snps <- target_num_snps0[target_num_snps0 <= Total_num_snps]

quant_df_list <- vector(mode = "list", length = length(target_num_snps));names(quant_df_list) <- paste0("SNP_target",target_num_snps)

# Choose How Close Realized Number of SNPs must be to Target Number
# Proportion
snp_num_stringency <- 0.0005


for(tt in 1:length(target_num_snps)){ ## loop through each target number of snps, filter as close as possible
  
  target_snps <- target_num_snps[tt]
  
  if(target_snps == 0 | target_snps == Total_num_snps){quant_df_list[[tt]] <- snp_results_DF1;next}
  
  stringency_target <- target_snps * snp_num_stringency
  message("### ### ### ### ### ### ### ### ###\nTarget SNPs: ",target_snps,"\n",Sys.time(),"\n### ### ### ### ### ### ### ### ###")
  
  closest_to_target <- snp_results_DF1[min(abs(snp_results_DF1[,"Sites"] - target_snps)) == abs(snp_results_DF1[,"Sites"] - target_snps),]
  
  diff_in_num_snps1 <- abs(target_snps - closest_to_target[,"Sites"])
  diff_in_num_snps2 <- diff_in_num_snps1
  diff_in_num_snps_repeat <- diff_in_num_snps1
  
  snp_results_DF_repeat <- snp_results_DF1
  
  while(diff_in_num_snps2 > stringency_target){
    
    # WHILE ITERATION #
    
    new_quant2 <- find.new.quantile.seq(target_snps, DF = snp_results_DF_repeat)
    test_seq2 <- new_quant2[[1]]
    snp_results_mat2 <- lapply(test_seq2, FUN = quantile.matrix.function)
    
    ## quantile dataframe ##
    snp_results_DF2_pre <- data.frame(matrix(unlist(snp_results_mat2), nrow=length(snp_results_mat2), byrow=T, dimnames=list(NULL, colnames(snp_results_mat2[[1]]))),stringsAsFactors=FALSE)
    snp_results_DF2A <- rbind(new_quant2[[2]],snp_results_DF2_pre)
    
    if( dim(snp_results_DF2_pre)[1] == length(test_seq2) ){ snp_results_DF2 <- rbind(snp_results_DF2A,new_quant2[[3]]) } else { snp_results_DF2 <- snp_results_DF2A }
    
    
    
    ## test number of snps ##
    closest_to_target2 <- snp_results_DF2[min(abs(snp_results_DF2[,"Sites"] - target_snps)) == abs(snp_results_DF2[,"Sites"] - target_snps),]
    
    diff_in_num_snps2 <- abs(target_snps - closest_to_target2[1,"Sites"])
    
    
    
    ## (if number of snps does not improve then stop) ##
    if(diff_in_num_snps2 == diff_in_num_snps_repeat){message("INCREASING STRINGENCY DID NOT IMPROVE NUMBER OF SNPS -- DONE");break}
    ## (two different quantile levels have the same number of snps, filtering within this range will not improve) ##
    if(dim(closest_to_target2)[1] >= 2){message("Difference in number of SNPs: ",diff_in_num_snps2);message("FURTHER FILTERING WILL NOT IMPROVE FIT -- DONE");break}
    
    
    ## save results for next iteration
    snp_results_DF_repeat <- snp_results_DF2
    diff_in_num_snps_repeat <- diff_in_num_snps2
    
    message("Difference in number of SNPs: ",diff_in_num_snps2)
    
    
  }
  
  
  quant_df_list[[tt]] <- snp_results_DF2  ## save quantile, value for each metric, and number of sites remaining
  
} ## end for loop



## error rates, by site and by gt, sites remaining, gt remaining, value for each metric, etc
RESULTS <- array(NA, dim = c(length(target_num_snps), length(quant_df_list[[1]])+13), dimnames = list(c(),c("Sites","GT","NA.gt","Match.abs.site","Miss.abs.site","Match.prop.site","Miss.prop.site","Match.abs.gt","Miss.abs.gt","Match.prop.gt","Miss.prop.gt",names(quant_df_list[[1]][-dim(quant_df_list[[1]])[2]]),"P1.QUAL","P1.GQ","P1.DP")))
correct_size_df_list <- vector(mode = "list", length = length(target_num_snps));names(quant_df_list) <- paste0("SNP_target",target_num_snps)
error_check_df <- array(NA, dim = c(length(target_num_snps), 9), dimnames = list(c(),c("Sites","quantile","QUAL","GQ","DP","QUAL.sites","GQ.sites","DP.sites","P1.sites")))

## Quantify the error rates, number of gt called, etc for each data set
for(tt in 1:length(target_num_snps)){
  
  target_snps2 <- target_num_snps[tt]
    
  if(target_snps2 == Total_num_snps){ FINAL_filter_df_correct_size <- snp_df_no_P1_down_filters_REVERSE
    first_closest_to_target <- quant_df_list[[tt]][1,]
    first_closest_to_target[,-dim(quant_df_list[[tt]][1,])[2]] <- NA } else {
    
    quant_df_var <- quant_df_list[[tt]]
    
    first_closest_to_target <- quant_df_var[min(abs(quant_df_var[,"Sites"] - target_snps2)) == abs(quant_df_var[,"Sites"] - target_snps2),]
    
    if(dim(first_closest_to_target)[1] > 1){ first_closest_to_target <- first_closest_to_target[1,] }
    
    realized_num_snps <- first_closest_to_target[,"Sites"]
    
    
    message("### ### ### ### ### ### ### ### ###\nTarget SNPs: ",target_snps2,"\nFilter at quantile: ",first_closest_to_target[,"quantile"],"\n",Sys.time(),"\n### ### ### ### ### ### ### ### ###")
    
    sites_to_filter_out_full_list2 <- lapply(SITE_tot_filters, FUN = site.filter.fun, quant.vals = first_closest_to_target, HACK=UGHACK)
    names(sites_to_filter_out_full_list2) <- SITE_tot_filters
    
    P1_QUAL_filt2 <- quantile(snp_df_P1$LP.mg0P1.QUAL, first_closest_to_target$quantile, na.rm = TRUE)

    
    unique_sites_to_remove_by_SITE2_pre <- NULL
    
    for(l in 1: length(sites_to_filter_out_full_list2)){
      
      unique_sites_to_remove_by_SITE2_pre <- unique(c(unique_sites_to_remove_by_SITE2_pre, sites_to_filter_out_full_list2[[l]]))
      
    }
    
    length(unique_sites_to_remove_by_SITE2_pre)
    cat("by SITE filters - Sites removed:",length(unique_sites_to_remove_by_SITE2_pre),"\n")
    
    
    keep_these_sites_df_GT0 <- snp_df_no_P1_down_filters_REVERSE[!snp_df_no_P1_down_filters_REVERSE$label %in% unique_sites_to_remove_by_SITE2_pre,]
    
    if(dim(keep_these_sites_df_GT0)[1] < realized_num_snps){ message("ERROR - by SITES filter - filtered more sites than expected");break}
    
    
    keep_these_sites_df_GT1 <- keep_these_sites_df_GT0 ## DF after filtering for by SITE metrics
    error_check_df[tt,6] <- dim(keep_these_sites_df_GT1)[1]
    
    sites_filtered_labels2 <- NULL
    
    ### fixed so it filters GQ and then DP iteratively and not 1 at a time ###
    
    GQ_filt2 <- as.numeric(first_closest_to_target$GQ)
    DP_filt2 <- as.numeric(first_closest_to_target$DP)
    
    
    ## filter GQ
    df_GQ1 <- keep_these_sites_df_GT1[grepl("GQ|label", names(keep_these_sites_df_GT1))]
    rownames(df_GQ1) <- df_GQ1[,"label"]
    
    filt_GQ1 <- df_GQ1[grepl("GQ", names(df_GQ1))] < GQ_filt2
    
    
    ## filter DP
    df_DP1 <- keep_these_sites_df_GT1[grepl("DP|label", names(keep_these_sites_df_GT1))]
    rownames(df_DP1) <- df_DP1[,"label"]
    
    filt_DP1 <- df_DP1[grepl("DP", names(df_DP1))] < DP_filt2
    
    
    if(identical(keep_these_sites_df_GT1$label,rownames(df_GQ1))==FALSE){message("ERROR A2");break}
    if(identical(keep_these_sites_df_GT1$label,rownames(df_DP1))==FALSE){message("ERROR B2");break}
    if(identical(rownames(df_GQ1),rownames(df_DP1))==FALSE){message("ERROR C2");break}
    
    
    keep_these_sites_df_GT2 <- keep_these_sites_df_GT1
    if("P1" %in% names(keep_these_sites_df_GT2)){message("ERROR D");break}
    
    GT_tab1 <- keep_these_sites_df_GT1[grepl("GT", names(keep_these_sites_df_GT1))]
    pre_GT_missingA <- rowSums(GT_tab1 == "./." | GT_tab1 == "." | GT_tab1 == ".|.", na.rm=TRUE)
    cat("GT missing PRE filtering:",sum(pre_GT_missingA),"\n")
    cat("Proportion:",sum(pre_GT_missingA)/(dim(GT_tab1)[1]*dim(GT_tab1)[2]),"\n")
    
    
    keep_these_sites_df_GT2[grepl("GT", names(keep_these_sites_df_GT2))][filt_GQ1] <- "./."
    keep_these_sites_df_GT3 <- keep_these_sites_df_GT2
    
    GT_tab2 <- keep_these_sites_df_GT2[grepl("GT", names(keep_these_sites_df_GT2))]
    pre_GT_missingA1 <- rowSums(GT_tab2 == "./." | GT_tab2 == "." | GT_tab2 == ".|.", na.rm=TRUE)
    cat("GT missing GQ filtering:",sum(pre_GT_missingA1),"\n")
    cat("Proportion:",sum(pre_GT_missingA1)/(dim(GT_tab2)[1]*dim(GT_tab2)[2]),"\n")
    error_check_df[tt,7] <- sum(pre_GT_missingA1 <= 53)
    
    keep_these_sites_df_GT3[grepl("GT", names(keep_these_sites_df_GT3))][filt_DP1] <- "./."
    keep_these_sites_df_GT4 <- keep_these_sites_df_GT3
    
    GT_tab3 <- keep_these_sites_df_GT3[grepl("GT", names(keep_these_sites_df_GT3))]
    pre_GT_missingA2 <- rowSums(GT_tab3 == "./." | GT_tab3 == "." | GT_tab3 == ".|.", na.rm=TRUE)
    cat("GT missing DP filtering:",sum(pre_GT_missingA2),"\n")
    cat("Proportion:",sum(pre_GT_missingA2)/(dim(GT_tab3)[1]*dim(GT_tab3)[2]),"\n")
    error_check_df[tt,8] <- sum(pre_GT_missingA2 <= 53)
    
    
    GT_tab4 <- keep_these_sites_df_GT4[grepl("GT", names(keep_these_sites_df_GT4))]
    GT_missingA <- rowSums(GT_tab4 == "./." | GT_tab4 == "." | GT_tab4 == ".|.", na.rm=TRUE)
    
    keep_these_sites_df_GT5A <- keep_these_sites_df_GT4[GT_missingA <= 53,]
    
    
    ## also remove sites where Parent is below GQ/DP threshold ##
    
    P1_GQ_filt2 <- quantile(snp_df_P1[grepl("GQ", names(snp_df_P1))], first_closest_to_target[,"quantile"], na.rm=TRUE)
    P1_DP_filt2 <- quantile(snp_df_P1[grepl("DP", names(snp_df_P1))], first_closest_to_target[,"quantile"], na.rm=TRUE)
    
    filt_P1_GQ2 <- snp_df_P1[snp_df_P1[grepl("GQ", names(snp_df_P1))] < P1_GQ_filt2,]
    
    filt_P1_DP2 <- snp_df_P1[snp_df_P1[grepl("DP", names(snp_df_P1))] < P1_DP_filt2,]
    
    P1_sites_to_remove2 <- unique(c(filt_P1_GQ2$label, filt_P1_DP2$label))
    
    
    keep_these_sites_df_GT5 <- keep_these_sites_df_GT5A[!keep_these_sites_df_GT5A$label %in% P1_sites_to_remove2,]
    
    
    realized_sites_removed_due_to_P1 <- sum(duplicated(c(P1_sites_to_remove2,keep_these_sites_df_GT5A$label)))
    error_check_df[tt,9] <- realized_sites_removed_due_to_P1
    
    
    realized_sites <- dim(keep_these_sites_df_GT5)[1]
    
    if(realized_sites != first_closest_to_target[,"Sites"]){message("ERROR - realized number of sites");break}
    
    
    sites_filtered_labels2 <- keep_these_sites_df_GT4[GT_missingA > 53,"label"]
    
    
    if(sum(length(sites_filtered_labels2), length(unique_sites_to_remove_by_SITE2_pre), first_closest_to_target[,"Sites"], realized_sites_removed_due_to_P1) != dim(snp_df)[1]){message("big error\n");break}
    
    
    
    cat("by GT filters - Sites removed: ",length(sites_filtered_labels2),"\n\n")  
    
    unique_sites_to_remove_by_SITE2 <- unique(c(unique_sites_to_remove_by_SITE2_pre,sites_filtered_labels2,P1_sites_to_remove2))
    
    cat("Total Sites Removed: ",length(unique_sites_to_remove_by_SITE2),"\n\n")  
    
    
    FINAL_filter_df <- keep_these_sites_df_GT5
    
    if(dim(FINAL_filter_df)[1] != realized_num_snps){message("ERROR -- number of snps after filtering does not match");break}
    
    FINAL_filter_df_correct_size <- FINAL_filter_df
    
    } ## if else -- target number of snps is the max in the data frame
  
  correct_size_df_list[[tt]] <- FINAL_filter_df_correct_size
  
  ### END OF CORRECTLY SIZING DATA ###
  
  target_snps2 <- dim(FINAL_filter_df_correct_size)[1]
  
  ## one missingness check / GQ DP filter check ##
  
  GT_tab_fin <- FINAL_filter_df_correct_size[!grepl("P1", names(FINAL_filter_df_correct_size)) & grepl("GT", names(FINAL_filter_df_correct_size))]
  
  missing_per_row0 <- rowSums(GT_tab_fin=="./." | GT_tab_fin=="." | GT_tab_fin==".|.", na.rm = TRUE)
  ## missingness based on NA data
  if(max(missing_per_row0) > 53){message("ERROR - GQ DP FILTER");break}
  ## missingness based on NA and Hetero data
  missing_per_row <- rowSums(GT_tab_fin=="./." | GT_tab_fin=="." | GT_tab_fin==".|." | GT_tab_fin=="H", na.rm = TRUE)
  
  good_gt_per_row <- rowSums(GT_tab_fin!="./." & GT_tab_fin!="." & GT_tab_fin!=".|." & GT_tab_fin!="H", na.rm = TRUE)
  
  full_num_gt <- sum(good_gt_per_row)
  full_na_gt <- sum(missing_per_row)
  
  
  
  
  ### BEGIN ASSESSING ERROR RATES ###
  
  ####################################
  ## VARSCAN HETERO ##
  check_for_heteros <- FINAL_filter_df_correct_size[(grepl("GT", names(FINAL_filter_df_correct_size)))]
  if( sum(grepl("H", check_for_heteros))>0 ){
    FINAL_filter_df_correct_size1 <- data.frame(lapply(FINAL_filter_df_correct_size[(grepl("GT", names(FINAL_filter_df_correct_size)))], function(x) {gsub("H", "./.", x)}))
    FINAL_filter_df_correct_size[(grepl("GT", names(FINAL_filter_df_correct_size)))] <- FINAL_filter_df_correct_size1
  }
  ####################################
  
  offspring_labels <- FINAL_filter_df_correct_size$label
  
  GT_tab0 <- FINAL_filter_df_correct_size[grepl("label|GT", names(FINAL_filter_df_correct_size))]
  GT_tab0_nolab <- GT_tab0[!grepl("label", names(GT_tab0))]
  
  num_GT_na <- sum(is.na(GT_tab0) | GT_tab0 == "." | GT_tab0 == ".|." | GT_tab0 == "./.")
  num_GT <- sum(!is.na(GT_tab0_nolab) & GT_tab0_nolab != "." & GT_tab0_nolab != ".|." & GT_tab0_nolab != "./.")
  if( (dim(GT_tab0)[1]*dim(GT_tab0_nolab)[2]) != (num_GT_na+num_GT) ){message("ERROR -- Number of GT does not add up");break}
  
  
  subset_parent_to_same_sites <- snp_df[snp_df$label %in% offspring_labels,]
  parental_GT00 <- subset_parent_to_same_sites[grepl(paste0("label|",P1name), names(subset_parent_to_same_sites))]
  parental_GT0 <- parental_GT00[grepl("label|GT", names(parental_GT00))]
  
  
  ## sort label names so that SNPs line up horizontally across both datasets
  GT_tab <- GT_tab0[order(GT_tab0$label),]
  parental_GT <- parental_GT0[order(parental_GT0$label),]
  
  if(identical(parental_GT$label,GT_tab$label)==FALSE){message("ERROR -- parent/offspring dataframes");break}
  if(identical(dim(parental_GT)[1],dim(GT_tab)[1])==FALSE){message("ERROR -- parent/offspring dataframes");break}
  # check for NA in parent
  if( sum(is.na(parental_GT[grepl("GT", names(parental_GT))]) | parental_GT[grepl("GT", names(parental_GT))] == "." | parental_GT[grepl("GT", names(parental_GT))] == "./." | parental_GT[grepl("GT", names(parental_GT))] == ".|.") >0 ){message("ERROR -- parental NA");break}
  
  # convert parent GT into single alleles
  P1_GT <- as.matrix(parental_GT[grepl("GT", names(parental_GT))])
  rownames(P1_GT) <- parental_GT$label
  P1_GT_A1 <- gsub("/.*|\\|.*","",P1_GT)
  P1_GT_A2 <- gsub(".*/|.*\\|","",P1_GT)
  
  # offspring GT
  off_GT <- as.matrix(GT_tab[!grepl("label", names(GT_tab))])
  rownames(off_GT) <- GT_tab$label
  
  
  # convert offspring GT into single alleles if needed (Varscan)
  off_GT_A1 <- gsub("/.*","",off_GT)
  off_GT_A2 <- gsub(".*/","",off_GT)
  off_GT_A1 <- gsub(".", NA, off_GT_A1, fixed=TRUE)
  off_GT_A2 <- gsub(".", NA, off_GT_A2, fixed=TRUE)
  
  
  if( identical(rownames(P1_GT_A1),rownames(P1_GT_A2))==FALSE ){message("ERROR label - P1 alleles");break}
  if( identical(rownames(P1_GT_A1),rownames(off_GT_A1))==FALSE ){message("ERROR label - P1/OFF alleles");break}
  if( identical(rownames(P1_GT_A1),rownames(off_GT_A2))==FALSE ){message("ERROR label - P1/OFF alleles");break}
  if( identical(rownames(off_GT_A1),rownames(off_GT_A2))==FALSE ){message("ERROR label - OFF alleles");break}
  
  
  # Test if haploid or diploid offspring data #
  if( identical(off_GT_A1,off_GT_A2) ){
    ## Haploid
    Test_GT <- (data.frame(off_GT_A1) == P1_GT_A1 | data.frame(off_GT_A1) == P1_GT_A2);cat("Haploid offspring data\n") 
  } else {
    ## Diploid
    Test_GT <- (data.frame(off_GT_A1) == P1_GT_A1 | data.frame(off_GT_A1) == P1_GT_A2) & (data.frame(off_GT_A2) == P1_GT_A1 | data.frame(off_GT_A2) == P1_GT_A2);cat("Diploid offspring data\n")
  }
  
  # BY GT MATCH ABS
  MATCHES_GT <- sum(Test_GT, na.rm = TRUE)
  MISMATCHES_GT <- sum(Test_GT==FALSE, na.rm = TRUE)
  NA_gt <- sum(is.na(Test_GT))
  
  if(num_GT_na != NA_gt){message("Different number of NA GT?");break}
  if(sum(MATCHES_GT,MISMATCHES_GT) != num_GT){message("Number of non NA GT don't add up?");break}
  if(full_num_gt != num_GT){message("ERROR F1");break}
  if(full_na_gt != NA_gt){message("ERROR F2");break}
  
  
  ## BY GT MATCH RATE
  GT_Match_Prop <- MATCHES_GT / num_GT
  GT_Mismatch_Prop <- MISMATCHES_GT / num_GT
  
  ## BY SITE MATCH 
  Mismatch_by_Site <- rowSums(Test_GT==FALSE, na.rm= TRUE)
  MATCH_SITE <- sum(Mismatch_by_Site == 0)
  MISATCH_SITE <- sum(Mismatch_by_Site > 0)
  if( sum(MATCH_SITE,MISATCH_SITE) != target_snps2){message("ERROR -- number of matches and mismatches doesn't equal number of SNPs");break}
  
  SITE_Match_Prop <- MATCH_SITE / target_snps2
  SITE_Mismatch_Prop <- MISATCH_SITE / target_snps2
  
  ######
  
  RESULTS[tt,1] <- dim(FINAL_filter_df_correct_size)[1]
  RESULTS[tt,2] <- num_GT
  RESULTS[tt,3] <- NA_gt
  RESULTS[tt,4] <- MATCH_SITE
  RESULTS[tt,5] <- MISATCH_SITE
  RESULTS[tt,6] <- SITE_Match_Prop
  RESULTS[tt,7] <- SITE_Mismatch_Prop
  RESULTS[tt,8] <- MATCHES_GT
  RESULTS[tt,9] <- MISMATCHES_GT
  RESULTS[tt,10] <- GT_Match_Prop
  RESULTS[tt,11] <- GT_Mismatch_Prop
  RESULTS[tt,c(12:c(10+length(quant_df_list[[1]])))] <- as.matrix(first_closest_to_target[1,-dim(first_closest_to_target)[2]])
  if(target_snps2 != Total_num_snps){ RESULTS[tt,c(11+length(quant_df_list[[1]]))] <- P1_QUAL_filt2
    RESULTS[tt,c(12+length(quant_df_list[[1]]))] <- P1_GQ_filt2
    RESULTS[tt,c(13+length(quant_df_list[[1]]))] <- P1_DP_filt2
    }
  
  ## Reverse SIGN on DOWN filters
  if( length(SITE_DOWN_filters)>0 & target_snps2 != Total_num_snps ){
    RESULTS[tt,colnames(RESULTS) %in% SITE_DOWN_filters] <- RESULTS[tt,colnames(RESULTS) %in% SITE_DOWN_filters] * -1
    # Make sure the values of the DOWN filters in my results exist in the ORIGINAL dataframe
    if(round(RESULTS[tt,colnames(RESULTS) %in% SITE_DOWN_filters[1]],3) %in% round(snp_df_no_P1_orig[,colnames(snp_df_no_P1_orig) %in% SITE_DOWN_filters[1]],3)){print("Down filter results exist in original (1)")} else {message("ERROR -- DOWN filter results do not exist in ORIGINAL dataset (1)");break}
    if(length(SITE_DOWN_filters)>1){if(round(RESULTS[tt,colnames(RESULTS) %in% SITE_DOWN_filters[2]],3) %in% round(snp_df_no_P1_orig[,colnames(snp_df_no_P1_orig) %in% SITE_DOWN_filters[2]],3)){print("Down filter results exist in original (2)")} else {message("ERROR -- DOWN filter results do not exist in ORIGINAL dataset (2)");break}}
  }
  
  
  error_check_df[tt,1] <- target_snps2
  error_check_df[tt,2] <- first_closest_to_target[1,]$quantile
  if("QUAL" %in% names(first_closest_to_target)){ error_check_df[tt,3] <- first_closest_to_target[1,]$QUAL }
  error_check_df[tt,4] <- first_closest_to_target[1,]$GQ
  error_check_df[tt,5] <- first_closest_to_target[1,]$DP
  # setwd("")
  # write.csv(error_check_df, paste0(SNPCALLER_NAME,"_error_check.csv"), row.names = FALSE)
  
  message("Finished SNPs RESULTS: ", target_snps2)

  
} ## end for loop

setwd("")
write.csv(RESULTS, paste0(SNPCALLER_NAME,"_Results.csv"), row.names = FALSE)

# Sites - number of sites remaining after filtering
# GT - number of genotype calls remaining after filtering
# NA.gt - number of NA genotype calls remaining after filtering
# Match.abs.site (gt) - absolute number of sites (genotypes) matching between P1 & F2
# Miss.abs.site (gt) - absolute number of sites (genotypes) mismatching between P1 & F2
# Match.prop.site (gt) - proportion of sites (genotypes) matching between P1 & F2
# Miss.prop.site (gt) - proportion of sites (genotypes) mismatching between P1 & F2
# quantile - quantile of the distribution of each metric the data set was filtered at
# QUAL - value of QUAL corresponding to the given quantile of the QUAL distribution in F2
# P1.QUAL - value of QUAL corresponding to the given quantile of the QUAL distribution in P1

