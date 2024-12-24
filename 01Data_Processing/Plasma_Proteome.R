#Plasma Proteome-----------------------------------------------------------------


#Load packages
library(TwoSampleMR)
library(ieugwasr)
library(data.table) 

#Set the P_value threshold
p_threshold = 1e-5

#Set the working directory
workpath = "G:/mengdeer"
setwd(workpath)

#Read data
data_info <- read.table("./Plasma_Proteome/a.txt",header=F,sep=" ")

start = 4501
end = 4910

write.table(c(),paste0(workpath,"/","Plasma_Proteome/res_significant_snp_clump.txt",start,"_",end),row.names=F,col.names=F,quote=F,sep="\t")



#----------------------------------Data Process----------------------------------

for (i in start:end) {
  
  file_name <- data_info$V11[i]
  
  dat1 <- fread(paste0(workpath, "/", file_name),
                header = T,
                sep = "\t")
 
  
  #Filter values above the threshold
  if (min(dat1$Pval) > p_threshold) {
    next
  } #If the P_value threshold is not satisfied, skip to the next file.
  
  dat1_res <- subset(dat1, Pval < p_threshold)

  
  dat1_res <- TwoSampleMR::format_data(
    dat1_res, #The data that was just loaded
    type = "exposure",
    snp_col = "rsids",
    beta_col = "Beta",
    se_col = "SE",
    eaf_col = "ImpMAF",
    effect_allele_col = "effectAllele",
    other_allele_col = "otherAllele",
    pval_col = "Pval",
    chr_col = "Chrom",
    pos_col = "Pos",
    samplesize_col = "N"
  )
  
  
  dat1_res$dataid = i
  
  print(paste0("read rs_", file_name, " ", i, " ok 4910"))
  
  
  
#----------------------------------Local Clump-----------------------------------
  
  #Clumping aims to remove redundant variants that are highly correlated due to linkage disequilibrium, retaining statistically independent and significant variants.
  b <- ld_clump(
    dplyr::tibble(
      rsid = dat1_res$SNP,
      pval = dat1_res$pval.exposure,
      id = dat1_res$dataid
    ),
    
    #get_plink_exe()
    plink_bin = "D:/R/R-4.2.1/library/plinkbinr/bin/plink_Windows.exe",
    
    #European population reference genome file path
    bfile = paste0(workpath, "/data/ref/EUR")
  )
  dat1_res_clump <- dat1_res[dat1_res$SNP %in% b$rsid, ]
  
  
  res <- dat1_res            #The results were statistically significant (P < P_threshold).
  if (dim(dat1_res_clump)[1] > 0) {
    res <- dat1_res_clump    #The results were statistically significant and clumped.
  }
  
  
  write.table(
    res[, c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure","pval.exposure","dataid")],
    paste0(workpath,"/Plasma_Proteome/res_significant_snp_clump.txt",start,"_",end),row.names = F,col.names = F,quote = F,sep = "\t",append = T)
  
}



#----------------------------------Output The Results----------------------------

#Combine the calculated results from the above sections into a single file: res_significant_snp_all_clump.txt.
snp_all_clump <- read.table("./Plasma_Proteome/res_significant_snp_clump.txt",header=F)
names(snp_all_clump) <- c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure","pval.exposure","dataid")


#bmi t2d smoke123
snp_5_clump <- read.table("./Plasma_Proteome/significant_clump_all_files.txt",header=T)


#Extract SNPs from all data in the snp_all file for each individual file, and save them separately.
data_info <- read.table("./Plasma_Proteome/data_info.txt",header=F,sep=" ")

#Export data_info
#write.table(data_info,paste0(workpath,"/Plasma_Proteome/data_info.txt"),quote=F,col.names=T,row.names=F,sep="\t")


start = 4910
end = 4001


#Iteratively output result files
for (i in start:end) {
  
  file_name <- data_info$V11[i]
  
  if(file.exists(paste0(workpath,"/Plasma_Proteome/Processed_data/",i,"_res_",file_name,".txt.gz"))|file.exists(paste0(workpath,"/Plasma_Proteome/Processed_data/",i,"_res_",file_name))){
    
    print(paste0("ok 4910 ",workpath,"/Plasma_Proteome/Processed_data/",i,"_res_",file_name,".txt.gz allready ok"))
    next
  }
  
  dat1 <- fread(paste0(workpath,"/Plasma_Proteome/",file_name),header=T,sep="\t")
  
  print(paste0("read rs_",file_name," ",i," ok"))
  
  index_i <- dat1$rsids %in% unique(snp_all_clump$SNP,snp_5_clump$SNP)
  dat1_res <- dat1[index_i,c("rsids","effectAllele","otherAllele","Beta","SE","Pval","ImpMAF")]  #2 alt,3 ref
  rm(dat)
  gc()
  
  
  print(dim(dat1_res))
  
  file_name <- unlist(strsplit(file_name,"\\."))[1]
  
  write.table(dat1_res,paste0(workpath,"/Plasma_Proteome/Processed_data/",i,"_res_",file_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
  R.utils::gzip(paste0(workpath,"/Plasma_Proteome/Processed_data/",i,"_res_",file_name,".txt"))
}





