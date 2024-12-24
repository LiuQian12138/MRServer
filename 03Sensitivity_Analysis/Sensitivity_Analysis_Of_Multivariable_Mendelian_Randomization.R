#Sensitivity Analysis Of Multivariable Mendelian Randomization-------------------


#MR_getTableFromFunctionR_mv_local_local-----------------------------------------

#perform Multivariable Heterogeneity Test
X.toString()
SummaryStats <- data.frame(SummaryStats)
Y.toString()

#Inverse Variance Weighting method(IVW)
MRMVObject_mrivw <- mr_mvivw(MRMVInputObject,
                             model = ivw_model, 
                             robust = ivw_robust,
                             distribution = ivw_distribution)

#get the result_mrivw table
result_mrivw <- data.frame(exposure = unlist(strsplit(exposure_name,"@@@")), 
                           outcome = outcome_name, 
                           estimate = unname(MRMVObject_mrivw@Estimate),
                           StdError = unname(MRMVObject_mrivw@StdError),
                           CILower = unname(MRMVObject_mrivw@CILower),
                           CIUpper = unname(MRMVObject_mrivw@CIUpper),
                           Pvalue = unname(MRMVObject_mrivw@Pvalue)) 

paste0(apply(result_mrivw, 1, paste0, collapse = "@"), collapse = "lqnb")


MRMVObject_mrivw@SNPs #parameter 1 of the result_mrivw table
round(MRMVObject_mrivw@RSE,4) #parameter 2 of the result_mrivw table
MRMVObject_mrivw@Heter.Stat[1] #parameter 3 of the result_mrivw table
round(MRMVObject_mrivw@Heter.Stat[2],4) #parameter 4 of the result_mrivw table


#MR-Egger method
MRMVObject_mregger <- mr_mvegger(MRMVInputObject, 
                                 orientate = egger_orientate, 
                                 distribution = egger_distribution)
#get the result_mregger table
result_mregger <- data.frame(exposure = unlist(strsplit(exposure_name,"@@@")), 
                             outcome = outcome_name, 
                             estimate = MRMVObject_mregger@Estimate,
                             StdError = MRMVObject_mregger@StdError.Est,
                             CILower = MRMVObject_mregger@CILower.Est,
                             CIUpper = MRMVObject_mregger@CIUpper.Est,
                             Pvalue = MRMVObject_mregger@Pvalue.Est)

result_mregger <- rbind(result_mregger, 
                        c(intercept), 
                        c(outcome_name), 
                        MRMVObject_mregger@Intercept,
                        MRMVObject_mregger@StdError.Int,
                        MRMVObject_mregger@CILower.Int,
                        MRMVObject_mregger@CIUpper.Int,
                        MRMVObject_mregger@Pvalue.Int)

paste0(apply(result_mregger, 1, paste0, collapse = "@"), collapse = "lqnb")

MRMVObject_mregger@SNPs #parameter 1 of the result_mregger table
round(MRMVObject_mregger@RSE,4) #parameter 2 of the result_mregger table
MRMVObject_mregger@Heter.Stat[1] #parameter 3 of the result_mregger table
round(MRMVObject_mregger@Heter.Stat[2],4) #parameter 4 of the result_mregger table
MRMVObject_mregger@Orientate #parameter 5 of the result_mregger table


#perform the F Test
library(MVMR)

Z.toString()
sres <- strength_mvmr(r_input = F.data, gencov = 0)
paste0(apply(round(sres,4), 1, paste0, collapse = "@"), collapse = "lqnb")


#perform the Q Test
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

pres$Qstat #Q1
round(pres$Qpval,4) #Q2
length(F.data[,2])-num-1 #Q3

#get the mvmr_ivw table
res <- cbind(Exposure= result_mrivw$exposure,ivw_mvmr(r_input = F.data, gencov = 0))
paste0(apply(res, 1, paste0, collapse = "@"), collapse = "lqnb")
A_sum = ivw_mvmr_lq(F.data, gencov = 0)

round(A_sum$sigma,4) #parameter 1 of the mvmr_ivw table
A_sum$df[2] #parameter 2 of the mvmr_ivw table
MRMVObject_mrivw@Model #model of the mvmr_ivw table



#MR_getTableFromFunctionR_mv_local_online----------------------------------------

#perform Multivariable Heterogeneity Test
X.toString()
SummaryStats <- data.frame(SummaryStats)
Y.toString()

#Inverse Variance Weighting method(IVW)
MRMVObject_mrivw <- mr_mvivw(MRMVInputObject,
                             model = ivw_model, 
                             robust = ivw_robust, 
                             distribution = ivw_distribution )

#get the result_mrivw table
result_mrivw <- data.frame(exposure = unlist(strsplit(exposure_name,"@@@")), 
                           outcome = outcome_name, 
                           estimate = unname(MRMVObject_mrivw@Estimate),
                           StdError = unname(MRMVObject_mrivw@StdError),
                           CILower = unname(MRMVObject_mrivw@CILower),
                           CIUpper = unname(MRMVObject_mrivw@CIUpper),
                           Pvalue = unname(MRMVObject_mrivw@Pvalue))

paste0(apply(result_mrivw, 1, paste0, collapse = "@"), collapse = "lqnb")

MRMVObject_mrivw@SNPs #parameter 1 of the result_mrivw table
round(MRMVObject_mrivw@RSE,4) #parameter 2 of the result_mrivw table
MRMVObject_mrivw@Heter.Stat[1] #parameter 3 of the result_mrivw table
round(MRMVObject_mrivw@Heter.Stat[2],4) #parameter 4 of the result_mrivw table

#MR-Egger method
MRMVObject_mregger <- mr_mvegger(MRMVInputObject, 
                                 orientate  = egger_orientate, 
                                 distribution = egger_distribution)


#get the result_mregger table
result_mregger <- data.frame(exposure = unlist(strsplit(exposure_name,"@@@")) , 
                             outcome = outcome_name, 
                             estimate = MRMVObject_mregger@Estimate,
                             StdError = MRMVObject_mregger@StdError.Est,
                             CILower = MRMVObject_mregger@CILower.Est,
                             CIUpper = MRMVObject_mregger@CIUpper.Est,
                             Pvalue = MRMVObject_mregger@Pvalue.Est)

result_mregger <- rbind(result_mregger, 
                        c(intercept, 
                          c(outcome_name), 
                          MRMVObject_mregger@Intercept,
                          MRMVObject_mregger@StdError.Int,
                          MRMVObject_mregger@CILower.Int,
                          MRMVObject_mregger@CIUpper.Int,
                          MRMVObject_mregger@Pvalue.Int))

paste0(apply(result_mregger, 1, paste0, collapse = "@"), collapse = "lqnb")

MRMVObject_mregger@SNPs #parameter 1 of the result_mregger table
round(MRMVObject_mregger@RSE,4) #parameter 2 of the result_mregger table
MRMVObject_mregger@Heter.Stat[1] #parameter 3 of the result_mregger table
round(MRMVObject_mregger@Heter.Stat[2],4) #parameter 4 of the result_mregger table
MRMVObject_mregger@Orientate #parameter 5 of the result_mregger table


#perform the F Test
library(MVMR)

Z.toString()
sres <- strength_mvmr(r_input = F.data, gencov = 0)
paste0(apply(round(sres,4), 1, paste0, collapse = "@"), collapse = "lqnb")


#perform the Q Test
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0)

pres$Qstat #Q1
round(pres$Qpval,4) #Q2
length(F.data[,2])-num-1 #Q3

#get the mvmr_ivw table
res <- cbind(Exposure = result_mrivw$exposure,
             ivw_mvmr(r_input = F.data, gencov = 0))
paste0(apply(res, 1, paste0, collapse = "@"), collapse = "lqnb")
A_sum = ivw_mvmr_lq(F.data, gencov = 0)

round(A_sum$sigma,4) #parameter 1 of the mvmr_ivw table
A_sum$df[2] #parameter 2 of the mvmr_ivw table
MRMVObject_mrivw@Model #model of the mvmr_ivw table















