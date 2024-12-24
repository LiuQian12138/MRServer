#Univariable Mendelian Randomization---------------------------------------------


#MR_getTableFromFunctionR_local_local--------------------------------------------

rDataPath <- "" #GBM_GSE141982_count.rData//NSCLC_GSE127471_count.rData//NSCLC_GSE131907_count.rData
source(rDataPath)

#perform local Univariable Mendelian Randomization analysis
dat <- mr1_local_local(exposure_file, 
                       outcome_file, 
                       p_threshold = pval) 

#check whether "dat" contains results
dim(dat)[1] <= 3


#perform MR-PRESSO
library(MRPRESSO)

mr_presso_res <- mr_Presso(dat, num = presso_num)

#check whether "mr_presso_res" contains results
if(is.data.frame(mr_presso_res) && nrow(mr_presso_res) == 0){
  print(TRUE)
} else {
  print(FALSE)
}

mr_presso_res$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_res$`MR-PRESSO results`$`Global Test`$Pvalue

dat <- mr_presso_snp(mr_presso_res, mr_presso_main, dat, type = "data")

#Obtain the table of harmonise results 
dat$exposure <- exposure_name
dat$outcome <-  outcome_name
paste0(apply(dat, 1, paste0, collapse = "@"), collapse = "lqnb")

#Obtain the table of MR results
res <- mr(dat,method_list = c(method))
res <- generate_odds_ratios(res) #calculate odds ratio
res$exposure <- exposure_name
res$outcome <- outcome_name
paste0(apply(res, 1, paste0, collapse = "@"), collapse = "lqnb")



#MR_getTableFromFunctionR_local_online-------------------------------------------

rDataPath <- "" #GBM_GSE141982_count.rData//NSCLC_GSE127471_count.rData//NSCLC_GSE131907_count.rData
source(rDataPath)

#perform online Univariable Mendelian Randomization analysis
dat <- mr1_local_online(exposure_file,
                        outcome_id,
                        p_threshold = pval)

#check whether "dat" contains results
dim(dat)[1] <= 3


#perform MR-PRESSO
library(MRPRESSO)

mr_presso_res <- mr_Presso(dat, num = presso_num)

#check whether "mr_presso_res" contains results
if(is.data.frame(mr_presso_res) && nrow(mr_presso_res) == 0){
  print(TRUE)
} else {
  print(FALSE)
}

mr_presso_res$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_res$`MR-PRESSO results`$`Global Test`$Pvalue

dat <- mr_presso_snp(mr_presso_res, mr_presso_main, dat, type = "data")

#Obtain the table of harmonise results 
dat$exposure <- exposure_name
dat$outcome <-outcome_name
paste0(apply(dat, 1, paste0, collapse = "@"), collapse = "lqnb")

#Obtain the table of MR results
res <- mr(dat,method_list = c(method))
res <- generate_odds_ratios(res) #calculate odds ratio
res$exposure <- exposure_name
res$outcome <- outcome_name
paste0(apply(res, 1, paste0, collapse = "@"), collapse = "lqnb")



#Univariable Mendelian Randomization Functions-----------------------------------

mr1_local_local <- function(filepath_exposure, filepath_outcome, p_threshold = 5e-5){
  
  library(TwoSampleMR)
  library(data.table)
  

  #-step1----exposure----------------------------
  
  
  dat_exposure_temp <- fread(filepath_exposure,sep="\t",header=T)
  colnames(dat_exposure_temp) <- c("SNP","alt","ref","beta","se","p","freq")
  dat_exposure_temp1 <- subset(dat_exposure_temp,p<p_threshold)
  
  
  if(nrow(dat_exposure_temp1) == 0){ #The threshold is too low, resulting in a limited number of SNPs.
    dat <- data.frame()
    return (dat)
  }
  
  
  dat_exposure <- TwoSampleMR::format_data( dat_exposure_temp1, 
                                            type = "exposure", #"outcome",
                                            snp_col = "SNP",
                                            beta_col = "beta",
                                            se_col = "se",
                                            effect_allele_col = "alt",
                                            other_allele_col = "ref",
                                            pval_col = "p",
  )
  
  
  #-step2----outcome-----------------------------
  
  
  dat_outcome_temp <- fread(filepath_outcome,sep="\t",header=T)
  
  colnames(dat_outcome_temp) <- c("SNP","alt","ref","beta","se","p","freq")
  
  dat_outcome <- TwoSampleMR::format_data( dat_outcome_temp, 
                                           type = "outcome",#"exposure",
                                           snp_col = "SNP",
                                           beta_col = "beta",
                                           se_col = "se",
                                           effect_allele_col = "alt",
                                           other_allele_col = "ref",
                                           pval_col = "p",
  )
  
  
  #-step3----harmonsise--------------------------
  
  
  dat <- harmonise_data(
    exposure_dat = dat_exposure, 
    outcome_dat = dat_outcome
  )
  
  if(nrow(dat) == 0){    #If there are no results from harmonsise(), proceed to this step.
    dat <- data.frame()
    return (dat)
  }else{
    dat <- subset(dat,select = c("SNP",
                                 "effect_allele.exposure",
                                 "other_allele.exposure",
                                 "effect_allele.outcome",
                                 "other_allele.outcome",
                                 "mr_keep",
                                 "remove",
                                 "palindromic" ,
                                 "ambiguous",
                                 
                                 "id.exposure",
                                 "exposure",
                                 "beta.exposure",
                                 "se.exposure",
                                 "pval.exposure",
                                 
                                 "id.outcome",
                                 "outcome",
                                 "beta.outcome",
                                 "se.outcome",
                                 "pval.outcome"
    ))
    
  }
  
  dat
  
}



mr1_local_online <- function(filepath_exposure, data_online_id, p_threshold = 5e-5){
  
  library(TwoSampleMR)
  library(data.table)
  
  
  #-step1----exposure----------------------------
  
  
  dat_exposure_temp <- fread(filepath_exposure,sep="\t",header=T)
  colnames(dat_exposure_temp) <- c("SNP","alt","ref","beta","se","p","freq")
  dat_exposure_temp1 <- subset(dat_exposure_temp,p < p_threshold)
  
  if(nrow(dat_exposure_temp1) == 0){ #The threshold is too low, resulting in a limited number of SNPs.
    dat <- data.frame()
    return (dat)
  }
  
  colnames(dat_exposure_temp1) <- c("SNP","alt","ref","beta","se","p","freq")
  
  
  dat_exposure <- TwoSampleMR::format_data( dat_exposure_temp1, 
                                            type = "exposure", #"outcome",
                                            snp_col = "SNP",
                                            beta_col = "beta",
                                            se_col = "se",
                                            effect_allele_col = "alt",
                                            other_allele_col = "ref",
                                            pval_col = "p",
  )
  
  
  #-step2----outcome-----------------------------
  
  
  dat_outcome <- extract_outcome_data(snps = dat_exposure$SNP,outcomes = data_online_id,proxies = F)
  
  if(is.null(dat_outcome) ){ #If the SNPs for the exposure do not match the outcome, resulting in no results, proceed to this step.
    dat <- data.frame()
    return (dat)
  }
  
  
  #-step3----harmonsise--------------------------
  
  
  dat <- harmonise_data(
    exposure_dat = dat_exposure, 
    outcome_dat = dat_outcome
  )
  
  if(nrow(dat) == 0){ #If there are no results from harmonsise(), proceed to this step.
    dat <- data.frame()
    return (dat)
  }else{
    dat <- subset(dat,select = c("SNP",
                                 "effect_allele.exposure",
                                 "other_allele.exposure",
                                 "effect_allele.outcome",
                                 "other_allele.outcome",
                                 "mr_keep",
                                 "remove",
                                 "palindromic" ,
                                 "ambiguous",
                                 
                                 "id.exposure",
                                 "exposure",
                                 "beta.exposure",
                                 "se.exposure",
                                 "pval.exposure",
                                 
                                 "id.outcome",
                                 "outcome",
                                 "beta.outcome",
                                 "se.outcome",
                                 "pval.outcome"
    ))
    
  }
  
  
  dat
  
}



mr <- function(dat, parameters = default_parameters(), method_list = subset(mr_method_list(), use_by_default)$obj)
{
  mr_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
  {
    x <- subset(x1, mr_keep)
    if(nrow(x) == 0)
    {
      message("No SNPs available for MR analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    } else {
      message("Analysing '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
    }
    res <- lapply(method_list, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    })
    methl <- mr_method_list()
    mr_tab <- data.frame(
      outcome = x$outcome[1],
      exposure = x$exposure[1],
      method = methl$name[match(method_list, methl$obj)],
      nsnp = sapply(res, function(x) x$nsnp),
      b = sapply(res, function(x) x$b),
      se = sapply(res, function(x) x$se),
      pval = sapply(res, function(x) x$pval)
    )
    mr_tab <- subset(mr_tab, !(is.na(b) & is.na(se) & is.na(pval)))
    return(mr_tab)
  })
  
  return(mr_tab)
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



mr_Presso <- function(dat, num = 10000){
  library(TwoSampleMR)
  library(MRPRESSO)
  library(dplyr)
  set.seed(123)
  try (mr_presso_res <- mr_presso(BetaOutcome ="beta.outcome", BetaExposure = "beta.exposure", SdOutcome ="se.outcome", SdExposure = "se.exposure", 
                                  OUTLIERtest = TRUE,DISTORTIONtest = TRUE, data = dat,  
                                  SignifThreshold = 0.05,NbDistribution = num))
  return(mr_presso_res)
  
}



mr_presso_snp <- function(mr_presso_res, mr_presso_main, dat, type = "list"){
  data_re <- list()
  if(type == "list"){
    for(i in 1:length(mr_presso_res)){
      res <- mr_presso_res$`MR-PRESSO results`[[i]]
      main <- mr_presso_main[[i]]
      data <- dat[[i]]
      try(if(is.na(main[2,6])==FALSE){
        outliers <- which(res$`Outlier Test`$Pvalue<0.05)
        data$mr_keep[outliers]<-FALSE
      })
      data_re[[i]] <- data
      names(data_re)[[i]] <- names(dat)[[i]]
    }
    return(data_re)
  }
  
  if(type=="data"){
    res <- mr_presso_res$`MR-PRESSO results`
    main <- mr_presso_main
    data <- dat
    try(if(is.na(main[2,6])==FALSE){
      outliers <- which(res$`Outlier Test`$Pvalue<0.05)
      data$mr_keep[outliers]<-FALSE
    })
    return(data)
  }
}



