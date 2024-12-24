#Multivariable Mendelian Randomization-------------------------------------------


#MR_getTableFromFunctionR_mv_local_local-----------------------------------------

rDataPath <- "" #GBM_GSE141982_count.rData//NSCLC_GSE127471_count.rData//NSCLC_GSE131907_count.rData
source(rDataPath) 

library(MendelianRandomization)

#perform local multivariable Mendelian Randomization analysis
mvdat <- mrm_local_local(exposure_file,
                         exposure_name,
                         outcome_file,
                         outcome_name,
                         p_threshold = pval)

#check whether "mvdat" contains results
if(is.list(mvdat) && length(mvdat$exposure_beta) <= 3 || is.data.frame(mvdat) && nrow(mvdat) == 0){
  print(TRUE)
} else { 
    print(FALSE)
  }


#perform MR-PRESSO
library(MRPRESSO)
library(tidyverse)

pressoA.toString()
df <- data.frame(df)
pressoB.toString()
mr_presso_res <- pressoC.toString()

#check whether "mr_presso_res" contains results
if(is.data.frame(mr_presso_res) && nrow(mr_presso_res) == 0){
  print(TRUE)
} else { 
  print(FALSE)
}

mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]]
mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]

res_mv <- method(mvdat,
                 intercept = intercept, 
                 instrument_specific = instrument_specific, 
                 pval_threshold = pval)

res_mv <- mv_subset_lq(mvdat,
                       intercept = intercept, 
                       instrument_specific = instrument_specific, 
                       pval_threshold = pval,
                       method = method)

#check whether "res_mv" contains results
if(is.null(res_mv$result) || nrow(res_mv$result) <= 0){
  print(TRUE)
} else { 
  print(FALSE)
}

res_mv_odds <- generate_odds_ratios(res_mv$result) #calculate odds ratio
paste0(apply(res_mv_odds, 1, paste0, collapse = "@"), collapse = "lqnb")






#MR_getTableFromFunctionR_mv_local_online----------------------------------------

rDataPath <- "" #GBM_GSE141982_count.rData//NSCLC_GSE127471_count.rData//NSCLC_GSE131907_count.rData
source(rDataPath) 

library(MendelianRandomization)

mvdat <- mrm_local_online(exposure_file,
                         exposure_name,
                         outcome_file,
                         p_threshold = pval)

#check whether "mvdat" contains results
if(is.list(mvdat) && length(mvdat$exposure_beta) <= 3 || is.data.frame(mvdat) && nrow(mvdat) == 0){
  print(TRUE)
} else { 
    print(FALSE)
  }


#perform MR-PRESSO
library(MRPRESSO)
library(tidyverse)

pressoA.toString()
df <- data.frame(df)
pressoB.toString()
mr_presso_res <- pressoC.toString()

#check whether "mr_presso_res" contains results
if(is.data.frame(mr_presso_res) && nrow(mr_presso_res) == 0){
  print(TRUE)
} else { 
  print(FALSE)
}

mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["RSSobs"]]
mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]

res_mv <- method(mvdat,
                 intercept = intercept, 
                 instrument_specific = instrument_specific, 
                 pval_threshold = pval)

res_mv <- mv_subset_lq(mvdat,
                       intercept = intercept, 
                       instrument_specific = instrument_specific, 
                       pval_threshold = pval,
                       method = method)

res_mv_odds <- generate_odds_ratios(res_mv$result) #calculate odds ratio
paste0(apply(res_mv_odds, 1, paste0, collapse = "@"), collapse = "lqnb")



#Multivariable Mendelian Randomization Functions---------------------------------

mrm_local_local <- function(filepathlist_exposure, phenotype_list, filepath_outcome, phnotype_outcome, p_threshold = 5e-5){
  
  filepathlist_exposure <- unlist(strsplit(filepathlist_exposure,"@@@"))
  phenotype_list <- unlist(strsplit(phenotype_list,"@@@"))
  
  
  library(TwoSampleMR)
  library(data.table)
  library(MendelianRandomization)
  
  
  #-step1----exposure----------------------------
  
  
  SNPs <- c()
  for(filepath_exposure in filepathlist_exposure){
    
    dat_exposure_temp <- fread(filepath_exposure,sep="\t",header=T)
    colnames(dat_exposure_temp) <- c("SNP","alt","ref","beta","se","p","freq")
    SNPs <- c(SNPs,subset(dat_exposure_temp,p<p_threshold)$SNP)
  }
  SNPs <- unique(SNPs)
  
  
  dat_exposure_all <- c()
  for(i in 1:length(filepathlist_exposure)){
    
    filepath_exposure = filepathlist_exposure[i]
    dat_exposure_temp <- fread(filepath_exposure,sep="\t",header=T)
    colnames(dat_exposure_temp) <- c("SNP","alt","ref","beta","se","p","freq")
    
    dat_exposure_temp <- dat_exposure_temp[dat_exposure_temp$SNP %in% SNPs,]
    
    dat_exposure_temp$phenotype = phenotype_list[i]
    
    dat_exposure <- TwoSampleMR::format_data( dat_exposure_temp,
                                              type = "exposure", #"outcome",
                                              phenotype_col ="phenotype",
                                              snp_col ="SNP",
                                              beta_col = "beta",
                                              se_col ="se",
                                              effect_allele_col ="alt",
                                              other_allele_col = "ref",
                                              pval_col = "p",
    )
    
    dat_exposure_all <- rbind(dat_exposure_all,dat_exposure)			
  }
  
  
  #-step2----outcome-----------------------------
  
  
  dat_outcome_temp <- fread(filepath_outcome,sep="\t",header=T)
  
  colnames(dat_outcome_temp) <- c("SNP","alt","ref","beta","se","p","freq")
  
  dat_outcome_temp$phenotype = phnotype_outcome
  
  dat_outcome <- TwoSampleMR::format_data( dat_outcome_temp, 
                                           type = "outcome", #"exposure",
                                           phenotype_col ="phenotype",
                                           snp_col ="SNP",
                                           beta_col = "beta",
                                           se_col = "se",
                                           effect_allele_col = "alt",
                                           other_allele_col = "ref",
                                           pval_col = "p",
  )
  dat_outcome <- dat_outcome[dat_outcome$SNP %in% SNPs,]
  
  
  #-step3----harmonsise--------------------------
  
  
  mvdat <- tryCatch({
    mv_harmonise_data(dat_exposure_all, dat_outcome)
  }, error = function(err) {
    data.frame()  #If an error occurs, return an empty data frame.
  })
  
  
  mvdat
  
}



mrm_local_online <- function(filepathlist_exposure, phenotype_list, data_online_id, p_threshold=5e-5){
  
  filepathlist_exposure <- unlist(strsplit(filepathlist_exposure,"@@@"))
  phenotype_list <- unlist(strsplit(phenotype_list,"@@@"))
  
  

  library(TwoSampleMR)
  library(data.table)
  
  p_threshold = 5e-5 #Set the P_value threshold
  
  
  #-step1----exposure----------------------------
  
  
  SNPs <- c()
  for(filepath_exposure in filepathlist_exposure){
    
    dat_exposure_temp <- fread(filepath_exposure,sep = "\t",header = T)
    colnames(dat_exposure_temp) <- c("SNP","alt","ref","beta","se","p","freq")
    SNPs <- c(SNPs,subset(dat_exposure_temp,p < p_threshold)$SNP)
  }
  SNPs <- unique(SNPs)
  
  
  dat_exposure_all <- c()
  for(i in 1:length(filepathlist_exposure)){
    
    filepath_exposure = filepathlist_exposure[i]
    dat_exposure_temp <- fread(filepath_exposure,sep="\t",header=T)
    colnames(dat_exposure_temp) <- c("SNP","alt","ref","beta","se","p","freq")
    
    dat_exposure_temp <- dat_exposure_temp[dat_exposure_temp$SNP %in% SNPs,]
    
    dat_exposure_temp$phenotype = phenotype_list[i]
    
    dat_exposure <- TwoSampleMR::format_data( dat_exposure_temp, 
                                              type = "exposure", #"outcome",
                                              phenotype_col ="phenotype",
                                              snp_col ="SNP",
                                              beta_col = "beta",
                                              se_col = "se",
                                              effect_allele_col = "alt",
                                              other_allele_col = "ref",
                                              pval_col = "p",
    )
    
    dat_exposure_all <- rbind(dat_exposure_all,dat_exposure)			
  }
  
  
  #-step2----outcome-----------------------------
  

  dat_outcome <- extract_outcome_data(snps = dat_exposure_all$SNP,outcomes = data_online_id,proxies = F)
  
  
  #-step3----harmonsise--------------------------
  
  
  mvdat <- mv_harmonise_data(dat_exposure_all, dat_outcome) 
  
  mvdat <- tryCatch({
    mv_harmonise_data(dat_exposure_all, dat_outcome)
  }, error = function(err) {
    data.frame()  #If an error occurs, return an empty data frame.
  })
  
  mvdat
  
}



mv_subset_lq <- function(mvdat, intercept = FALSE, instrument_specific = FALSE, pval_threshold = 5e-8, method = 'mv_multiple', plots = FALSE)
{
  features = mv_lasso_feature_selection(mvdat)
  
  if(nrow(features) == 0 ){
    res_mv <- list(result = data.frame())  #Initialize an empty 'result' data frame.
    return(res_mv)
  }else{
    #Update mvdat object
    mvdat$exposure_beta <- mvdat$exposure_beta[, features$exposure, drop=FALSE]
    mvdat$exposure_se <- mvdat$exposure_se[, features$exposure, drop=FALSE]
    mvdat$exposure_pval <- mvdat$exposure_pval[, features$exposure, drop=FALSE]
    
    #Find relevant instruments
    instruments <- apply(mvdat$exposure_pval, 1, function(x) any(x < pval_threshold))
    stopifnot(sum(instruments) > nrow(features))
    
    mvdat$exposure_beta <- mvdat$exposure_beta[instruments,drop=FALSE]
    mvdat$exposure_se <- mvdat$exposure_se[instruments,drop=FALSE]
    mvdat$exposure_pval <- mvdat$exposure_pval[instruments,drop=FALSE]	
    mvdat$outcome_beta <- mvdat$outcome_beta[instruments]
    mvdat$outcome_se <- mvdat$outcome_se[instruments]
    mvdat$outcome_pval <- mvdat$outcome_pval[instruments]
    
    if(method == 'mv_multiple'){
      res_mv <- mv_multiple(mvdat, intercept=intercept, instrument_specific=instrument_specific, pval_threshold=pval_threshold, plots=plots)
    }else{
      res_mv <- mv_residual(mvdat, intercept=intercept, instrument_specific=instrument_specific, pval_threshold=pval_threshold, plots=plots)
    }
    
    if (is.null(res_mv$result) || nrow(res_mv$result) == 0) {
      res_mv$result <- data.frame()  #Ensure that the returned result always has a 'result' attribute.
    }
    
    return (res_mv)
    
  }
  
}



generate_odds_ratios <- function(mr_res)
{
  mr_res$lo_ci <- mr_res$b - 1.96 * mr_res$se
  mr_res$up_ci <- mr_res$b + 1.96 * mr_res$se
  mr_res$or <- exp(mr_res$b)
  mr_res$or_lci95 <- exp(mr_res$lo_ci)
  mr_res$or_uci95 <- exp(mr_res$up_ci)
  return(mr_res)
}








