#Sensitivity Analysis Of Univariable Mendelian Randomization---------------------


#MR_getTableFromFunctionR_local_local--------------------------------------------

#perform Univariable Heterogeneity Test
het <- mr_heterogeneity(dat,method_list = c(method_heterogeneity)) 
paste0(apply(het, 1, paste0, collapse = "@"), collapse = "lqnb") #Obtain the table of heterogeneity test results


#perform Horizontal Pleiotropy Test
ple <- mr_pleiotropy_test(dat)
paste0(apply(ple, 1, paste0, collapse = "@"), collapse = "lqnb") #Obtain the table of pleiotropy test results



#MR_getTableFromFunctionR_local_online-------------------------------------------

#perform Univariable Heterogeneity Test
het <- mr_heterogeneity(dat,method_list = c(method_heterogeneity))
paste0(apply(het, 1, paste0, collapse = "@"), collapse = "lqnb") #Obtain the table of heterogeneity test results


#perform Horizontal Pleiotropy Test
ple <- mr_pleiotropy_test(dat)
paste0(apply(ple, 1, paste0, collapse = "@"), collapse = "lqnb") #Obtain the table of pleiotropy test results



