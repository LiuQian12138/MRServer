#Immune Cells--------------------------------------------------------------------


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
data_info <- read.table("./Immune_Cells/data_info.txt",header=F,sep=" ")

start = 1
end = 731

write.table(c(),paste0(workpath,"/Immune_Cells/res_significant_snp_clump.txt",start,"_",end),row.names=F,col.names=F,quote=F,sep="\t")



#----------------------------------Data Process----------------------------------

for (i in start:end) {
  
  file_name <- data_info$GWAS.Catalog.Accession.Number[i]
  
  dat1 <- read_exposure_data( #Note: After importing the data, the column names "effect_allele_col" and "other_allele_col" might switch order."
    filename = paste0(workpath,"/Immune_Cells/rs_",file_name,"_buildGRCh37.tsv.gz"),
    clump = F, #After filtering by P-values, proceed to clumping.
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    effect_allele_col ="A2",
    other_allele_col = "A1",
    eaf_col = "effect_allele_frequency",
    pval_col = "P"
  )
  
  
  print(paste0("read rs_",file_name," ",i," ok"))
  

  if(min(dat1$pval.exposure) > p_threshold){
    next
    } #If the P_value threshold is not satisfied, skip to the next file.
  
  dat1_res <- subset(dat1,pval.exposure < p_threshold)
  dat1_res$dataid = i
  
  print(dim(dat1_res))
  dat1_res_clump <- clump_data(dat1_res) #clump data 
  print(dim(dat1_res_clump))
  
  #The clumping process might fail due to network issues, resulting in zero rows of output.
  
  for(n in 1:30){
    
    if(dim(dat1_res_clump)[1] > 0){	
      break
    }else{
      print(paste0("Retry clumping data ",n))
      dat1_res_clump <- clump_data(dat1_res) #Attempt clumping the data again to resolve the issue.
    }
  }
  
  res <- dat1_res                    #The results were statistically significant (P < P_threshold).
  if(dim(dat1_res_clump)[1] > 0){
    res <- dat1_res_clump            #The results were statistically significant and clumped.
  }
  
  
  write.table(res[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure","pval.exposure","dataid")],paste0(workpath,"/Immune_Cells/res_significant_snp_clump.txt",start,"_",end),row.names=F,col.names=F,quote=F,sep="\t",append=T)
}



#----------------------------------Output The Results----------------------------

#Combine the calculated results from the above sections into a single file: res_significant_snp_all_clump.txt.
snp_all_clump <- read.table("./Immune_Cells/res_significant_snp_clump.txt",header=F)
names(snp_all_clump) <- c("SNP", "other_allele.exposure" ,"effect_allele.exposure", "beta.exposure", "pval.exposure", "se.exposure",'dataid')

#bmi t2d smoke123
snp_5_clump <- read.table("./Immune_Cells/significant_clump_all_files.txt.txt",header=T)


#Extract SNPs from all data in the snp_all file for each individual file, and save them separately.
data_info <- read.table("./Immune_Cells/data_info.txt",header=T,sep="\t")


start = 1
end = 731


#Iteratively output result files
for (i in start:end) {
  file_name <- data_info$GWAS.Catalog.Accession.Number[i]
  
  dat1<-read.table(paste0(workpath,"/Immune Cells/rs_",file_name,"_buildGRCh37.tsv.gz"),header=T)
  
  print(paste0(workpath,"/read rs_",file_name," ",i," ok"))
  
  index_i <- dat1$SNP %in% unique(c(snp_all_clump$SNP,snp_5_clump$SNP))
  dat1_res <- dat1[index_i,c("SNP","A2","A1","BETA","SE","P","EFFECT_ALLELE_FREQUENCY")]  #A2 alt,A1 ref
  
  
  print(dim(dat1_res))
  
  write.table(dat1_res,paste0(workpath,"/Immune_Cells/Processed_data/",i,"_res_",file_name,".txt"),row.names=F,col.names=T,quote=F,sep="\t")
  R.utils::gzip(paste0(workpath,"/Immune_Cells/Processed_data/",i,"_res_",file_name,".txt"))
}






