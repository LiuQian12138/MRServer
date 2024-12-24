#Circulating Proteins------------------------------------------------------------


#Load packages
library(TwoSampleMR)
library(ieugwasr)
library(data.table) 

#Set the P_value threshold
p_threshold = 1e-5

#Set the working directory
workpath <- "G:/mengdeer"
setwd(workpath)

#Read data
data_info <- read.table("./Circulating_Proteins/PMID37563310_studies_export.tsv",header=T,sep="\t")


start = 1
end = 91

write.table(c(),paste0(workpath,"/Circulating_Proteins/res_significant_snp_clump.txt",start,"_",end),row.names=F,col.names=F,quote=F,sep="\t")



#----------------------------------Data Process----------------------------------

for (i in start:end) {
  
  file_name <- data_info$accessionId[i]

  dat1 <- read_exposure_data(
    filename = paste0(workpath,"/Circulating_Proteins/",file_name,".tsv.gz"),
    clump = F, #After filtering by P-values, proceed to clumping.
    sep = "\t",
    snp_col = "rsid",
    beta_col = "beta",
    se_col = "standard_error",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "effect_allele_frequency",
    pval_col = "p_value"
  )
  
  
  print(paste0("read rs_",file_name," ",i," ok"))
  
  
  if(min(dat1$pval.exposure) > p_threshold){next} #If the P_value threshold is not satisfied, skip to the next file.
  dat1_res <- subset(dat1,pval.exposure < p_threshold)
  dat1_res$dataid = i
  

  
#----------------------------------Local Clump-----------------------------------

    b <- ld_clump(
    dplyr::tibble(rsid=dat1_res$SNP, pval=dat1_res$pval.exposure, id=dat1_res$dataid),
    
    #get_plink_exe()
    plink_bin = "D:/R/R-4.2.1/library/plinkbinr/bin/plink_Windows.exe",
    
    #European population reference genome file path
    bfile = paste0(workpath,"/data/ref/EUR")
  )
  dat1_res_clump <- dat1_res[dat1_res$SNP %in% b$rsid,]
  
  
  
  res <- dat1_res                    #The results were statistically significant (P < P_threshold).
  if(dim(dat1_res_clump)[1] > 0){
    res <- dat1_res_clump            #The results were statistically significant and clumped.
  }
  
  write.table(res[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure","pval.exposure","dataid")],paste0(workpath,"/Circulating_Proteins/res_significant_snp_clump.txt",start,"_",end),row.names=F,col.names=F,quote=F,sep="\t",append=T) 
}



#----------------------------------Output The Results----------------------------

#Combine the calculated results from the above sections into a single file: res_significant_snp_all_clump.txt.
snp_all_clump <- read.table("./Circulating_Proteins/res_significant_snp_clump.txt",header=F)
names(snp_all_clump) <- c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure","pval.exposure","dataid")


#bmi t2d smoke123
snp_5_clump <- read.table("./Circulating_Proteins/significant_clump_all_files.txt",header=T)


#Extract SNPs from all data in the snp_all file for each individual file, and save them separately.
data_info <- read.table("./Circulating_Proteins/PMID37563310_studies_export.tsv",header=T,sep="\t")


start = 1
end = 91


#Iteratively output result files
for (i in start:end) {
  
  file_name <- data_info$accessionId[i]
  
  dat1 <- read.table(paste0(workpath,"/Circulating_Proteins/",file_name,".tsv.gz"),header=T,sep="\t")
  
  print(paste0("read rs_",file_name," ",i," ok"))
  

  index_i <- dat1$rsid %in% unique(c(snp_all_clump$SNP,snp_5_clump$SNP))
  dat1_res <- dat1[index_i,c("variant_id","effect_allele","other_allele","beta","standard_error","p_value","effect_allele_frequency")]  #2 alt,3 ref
  
  
  print(dim(dat1_res))
  
  write.table(dat1_res,paste0(workpath,"/Circulating_Proteins/Processed_data/",i,"_res_",file_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
  R.utils::gzip(paste0(workpath,"/Circulating_Proteins/Processed_data/",i,"_res_",file_name,".txt"))
}



