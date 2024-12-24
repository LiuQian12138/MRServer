#Sensitivity Analysis Functions--------------------------------------------------


#Univariable Mendelian Randomization Sensitivity Analysis Functions--------------

mr_heterogeneity <- function(dat, parameters = default_parameters(), method_list = subset(mr_method_list(), heterogeneity_test & use_by_default)$obj)
{
  het_tab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
  {
    x <- subset(x1, mr_keep)
    if(nrow(x) < 2)
    {
      message("Not enough SNPs available for Heterogeneity analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    res <- lapply(method_list, function(meth)
    {
      get(meth)(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, parameters)
    })
    methl <- mr_method_list()
    het_tab <- data.frame(
      outcome = x$outcome[1],
      exposure = x$exposure[1],
      method = methl$name[match(method_list, methl$obj)],
      Q = sapply(res, function(x) x$Q),
      Q_df = sapply(res, function(x) x$Q_df),
      Q_pval = sapply(res, function(x) x$Q_pval)
    )
    het_tab <- subset(het_tab, !(is.na(Q) & is.na(Q_df) & is.na(Q_pval)))
    return(het_tab)
  })
  
  return(het_tab)
}



mr_pleiotropy_test <- function(dat)
{
  ptab <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(x1)
  {
    x <- subset(x1, mr_keep)
    if(nrow(x) < 2)
    {
      message("Not enough SNPs available for pleiotropy analysis of '", x1$id.exposure[1], "' on '", x1$id.outcome[1], "'")
      return(NULL)
    }
    res <- mr_egger_regression(x$beta.exposure, x$beta.outcome, x$se.exposure, x$se.outcome, default_parameters())
    out <- data.frame(
      outcome = x$outcome[1],
      exposure = x$exposure[1],
      egger_intercept = res$b_i,
      se = res$se_i,
      pval = res$pval_i
    )
    return(out)
  })
  return(ptab)
}



mr_leaveoneout <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
  if(!"samplesize.outcome" %in% names(dat))
  {
    dat$samplesize.outcome <- NA
  }
  
  stopifnot("outcome" %in% names(dat))
  stopifnot("exposure" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))
  
  
  res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(X)
  {
    x <- subset(X, mr_keep)
    nsnp <- nrow(x)
    if(nsnp == 0)
    {
      x <- X[1,]
      d <- data.frame(
        SNP = "All",
        b = NA,
        se = NA,
        p = NA,
        samplesize = NA,
        outcome = x$outcome[1],
        exposure = x$exposure[1]
      )
      return(d)
    }
    if(nsnp > 2)
    {
      l <- lapply(1:nsnp, function(i)
      {
        with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i], parameters))
      })
      l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      d <- data.frame(
        SNP = c(as.character(x$SNP), "All"),
        b = sapply(l, function(y) y$b),
        se = sapply(l, function(y) y$se),
        p = sapply(l, function(y) y$pval),
        samplesize = x$samplesize.outcome[1]
      )
      d$outcome <- x$outcome[1]
      d$exposure <- x$exposure[1]
      
    } else {
      a <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      d <- data.frame(
        SNP = "All",
        b = a$b,
        se = a$se,
        p = a$pval,
        samplesize = x$samplesize.outcome[1]
      )
      d$outcome <- x$outcome[1]
      d$exposure <- x$exposure[1]
    }
    return(d)
  })
  res <- subset(res, select=c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p))
  return(res)
}



mr_singlesnp <- function(dat, parameters=default_parameters(), single_method="mr_wald_ratio", all_method=c("mr_ivw", "mr_egger_regression"))
{
  
  if(!"samplesize.outcome" %in% names(dat))
  {
    dat$samplesize.outcome <- NA
  }
  
  stopifnot("outcome" %in% names(dat))
  stopifnot("exposure" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))
  
  res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(X)
  {
    x <- subset(X, mr_keep)
    nsnp <- nrow(x)
    if(nsnp == 0)
    {
      x <- X[1,]
      d <- data.frame(
        SNP = "No available data",
        b = NA,
        se = NA,
        p = NA,
        samplesize = NA,
        outcome = x$outcome[1],
        exposure = x$exposure[1]
      )
      return(d)
    }
    l <- lapply(1:nsnp, function(i)
    {
      with(x, get(single_method)(beta.exposure[i], beta.outcome[i], se.exposure[i], se.outcome[i], parameters))
    })
    nom <- c()
    for(i in seq_along(all_method))
    {
      l[[nsnp+i]] <- with(x, get(all_method[i])(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      
      nom <- c(nom, paste0("All - ", subset(mr_method_list(), obj==all_method[i])$name))
    }
    
    d <- data.frame(
      SNP = c(as.character(x$SNP), nom),
      b = sapply(l, function(y) y$b),
      se = sapply(l, function(y) y$se),
      p = sapply(l, function(y) y$pval),
      samplesize = x$samplesize.outcome[1]
    )
    d$outcome <- x$outcome[1]
    d$exposure <- x$exposure[1]
    return(d)
  })
  res <- subset(res, select=c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p))
  return(res)
}



#Multivariable Mendelian Randomization Sensitivity Analysis Functions------------

ivw_mvmr_lq <- function(r_input, gencov = 0){
  
  #convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }
  
  #Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }
  
  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  
  #If weights is missing, first order weights are used by default.
  
  
  #Inverse variance weighting is used.
  
  Wj <- 1/r_input[,3]^2
  
  #Determine the number of exposures included in the model
  
  exp.number <- length(names(r_input)[-c(1,2,3)])/2
  
  #Fit the IVW MVMR model
  
  A_sum <- summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))
  
  return(A_sum)
  
}



ivw_mvmr <- function(r_input,gencov = 0){
  
  #convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }
  
  #Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }
  
  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  
  #If weights is missing, first order weights are used by default.
  
  
  #Inverse variance weighting is used.
  
  Wj <- 1/r_input[,3]^2
  
  #Determine the number of exposures included in the model
  
  exp.number <- length(names(r_input)[-c(1,2,3)])/2
  
  #Fit the IVW MVMR model
  
  A_sum <- summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))
  
  A <- summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))$coef
  
  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i] <- paste0("exposure",i,collapse="")
  }
  
  
  ##########
  # Output #
  ##########
  
  #Print a few summary elements that are common to both lm and plm model summary objects
  cat("\n")
  
  cat("Multivariable MR\n")
  
  cat("\n")
  
  print(A)
  
  cat("\nResidual standard error:", round(A_sum$sigma,3), "on", A_sum$df[2], "degrees of freedom")
  
  cat("\n")
  
  cat("\n")
  
  cat("\n")
  
  
  return(A)
  
}



strength_mvmr <- function(r_input, gencov = 0){
  
  #convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }
  
  #Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }
  
  #gencov is the covariance between the effect of the genetic variants on each exposure.
  #By default it is set to 0.
  
  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  
  
  #Inverse variance weighting is used.
  
  Wj <- 1/r_input[,3]^2
  
  #Determine the number of exposures included in the model
  
  exp.number <- length(names(r_input)[-c(1,2,3)])/2
  
  A <- summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,data=r_input))$coef
  
  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i] <- paste0("exposure",i,collapse="")
  }
  
  #############################################
  # Generalised instrument strength het.stats #
  #############################################
  
  #Create an empty matrix for delta value (coefficients regressing each set of exposure effects upon other
  #exposure effects)
  
  delta_mat <- matrix(0,ncol=exp.number,nrow=exp.number-1)
  
  
  #Obtain delta values fitting regression models for each set of exposure effects upon other exposure effects
  for(i in 1:exp.number){
    regressand <- names(r_input[3 + i])
    regressors <- names(r_input)[-c(1,2,3,
                                  4+exp.number:length(names(r_input)))]
    C<-paste(regressand, "~", "-1 +", paste(regressors[-i], collapse="+"))
    D.reg<-stats::lm(C,data=r_input)
    delta_mat[,i] <- D.reg$coefficients
  }
  
  sigma2xj_dat <- matrix(ncol=exp.number,nrow=length(r_input[,1]),0)
  
  if(length(gencov) < 2){
    
    #Create a subset containing only standard errors for exposure effect estimates
    sebetas <- r_input[,(exp.number + 4):length(r_input)]
    
    #Generates the sigma2xj values for each exposure
    for(i in 1:exp.number){
      se.temp <- as.matrix(sebetas[,-i])
      for(j in 1:(exp.number-1)){
        sigma2xj_dat[,i] <- sigma2xj_dat[,i] + (se.temp[,j]^2 * delta_mat[j,i]^2)
      }
      sigma2xj_dat[,i] <- sigma2xj_dat[,i] + sebetas[,i]^2
      
    }
    
    
  }
  
  if(length(gencov) > 2){
    sigma2xj_dat <- matrix(ncol=exp.number,nrow=length(r_input[,1]),0)
    delta.temp <- matrix(0,ncol=exp.number,nrow=exp.number)
    
    
    #Generates the sigma2xj values for each exposure
    for(i in 1:exp.number){
      
      if(i == 1) {
        delta.temp[,i] <- c(-1, delta_mat[,i])
      }
      
      if(i>1 && i<exp.number){
        delta.temp[,i] <- c(delta_mat[1:(i-1),i],-1,delta_mat[i:(exp.number-1),i])
      }
      
      if(i == exp.number){
        delta.temp[,i] <- c(delta_mat[,i],-1)
      }
      
      
      for(l in seq_along(r_input[,1])){
        
        sigma2xj_dat[l,i] <- sigma2xj_dat[l,i] +  t(delta.temp[,i])%*%gencov[[l]]%*%delta.temp[,i]
        
      }
      
      
    }
  }
  
  #Create an empty matrix for instrument strength Q statistics
  Q_strength <- matrix(ncol=exp.number,nrow=1,0)
  
  #Generates the component of the Q statistic to be subtracted from the exposure estimates
  for(i in 1:exp.number){
    betas <- r_input[,c(4:(3+exp.number))]
    betas <- data.frame(betas[,-i])
    temp.sub <- 0
    for(j in 1:(exp.number-1)){
      temp.sub <- temp.sub + (delta_mat[j,i] * betas[,j])
    }
    
    #Populates matrix of Q statistics with respect to instrument strength
    Q_strength[i] <- sum((1/sigma2xj_dat[,i]) * ((r_input[,3+i] - temp.sub)^2))
    Q_strength[i] <- Q_strength[i]/nrow(r_input)
    
  }
  
  Q_strength <- data.frame(Q_strength)
  names(Q_strength) <- dimnames(A)[[1]]
  rownames(Q_strength) <- "F-statistic"
  
  
  ##########
  # Output #
  ##########
  
  #Print a few summary elements that are common to both lm and plm model summary objects
  cat("\n")
  
  cat("Conditional F-statistics for instrument strength\n")
  
  cat("\n")
  
  print(Q_strength)
  
  return(Q_strength)
  
}



pleiotropy_mvmr <- function(r_input, gencov = 0){
  
  #convert MRMVInput object to mvmr_format
  if ("MRMVInput" %in% class(r_input)) {
    r_input <- mrmvinput_to_mvmr_format(r_input)
  }
  
  #Perform check that r_input has been formatted using format_mvmr function
  if(!("mvmr_format" %in%
       class(r_input))) {
    stop('The class of the data object must be "mvmr_format", please resave the object with the output of format_mvmr().')
  }
  
  #gencov is the covariance between the effect of the genetic variants on each exposure.
  #By default it is set to 0.
  
  if(!is.list(gencov) && gencov == 0) {
    warning("Covariance between effect of genetic variants on each exposure not specified. Fixing covariance at 0.")
  }
  
  #Inverse variance weighting is used.
  
  Wj <- 1/r_input[,3]^2
  
  #Determine the number of exposures included in the model
  
  exp.number <- length(names(r_input)[-c(1,2,3)])/2
  
  #Fit the IVW MVMR model
  
  A_sum <- summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))
  
  A <- summary(stats::lm(stats::as.formula(paste("betaYG~ -1 +", paste(names(r_input)[
    seq(4,3+exp.number,by=1)], collapse="+")))
    ,weights=Wj,data=r_input))$coef
  
  #Rename the regressors for ease of interpretation
  for(i in 1:exp.number){
    dimnames(A)[[1]][i] <- paste0("exposure",i,collapse="")
  }
  
  
  #Create a subset containing only standard errors for exposure effect estimates
  sebetas <- r_input[,(exp.number + 4):length(r_input)]
  
  
  ########################
  ## Instrument Validity #
  ########################
  
  if(length(gencov) < 2){
    
    #Generate Sigma^2_A values
    sigma2A <- r_input[,3]^2
    for(i in 1:exp.number){
      sigma2A <- sigma2A + (A[i]^2 * sebetas[,i]^2)
    }
    
    #Create a subset of exposure effect estimates
    betas <- r_input[,c(4:(3+exp.number))]
    
    #Generates the component of the Q statistic to be subtracted from the outcome estimates
    temp.sub2 <- 0
    for(i in 1:exp.number){
      temp.sub2 <- temp.sub2 + (betas[,i] * A[i])
    }
    
    #Calculates Q statistic for instrument validity
    Q_valid <- sum((1/sigma2A)*(r_input[,2]-temp.sub2)^2)
    
    #Calculates p_value for instrument validity
    Q_chiValid <-stats::pchisq(Q_valid,length(r_input[,2])-exp.number-1,lower.tail = FALSE)
    
    
  }
  
  if(length(gencov) > 2){
    
    #Generate Sigma^2_A values
    sigma2A <- r_input[,3]^2
    for(i in seq_along(r_input[,3])){
      sigma2A[i] <- sigma2A[i] + (t(as.matrix(A[,1])) %*% gencov[[i]]%*% as.matrix(A[,1]))
    }
    
    #Create a subset of exposure effect estimates
    betas <- r_input[,c(4:(3+exp.number))]
    
    #Generates the component of the Q statistic to be subtracted from the outcome estimates
    temp.sub2 <- 0
    for(i in 1:exp.number){
      temp.sub2 <- temp.sub2 + (betas[,i] * A[i])
    }
    
    #Calculates Q statistic for instrument validity
    Q_valid <- sum((1/sigma2A)*(r_input[,2]-temp.sub2)^2)
    
    #Calculates p_value for instrument validity
    Q_chiValid <- stats::pchisq(Q_valid,length(r_input[,2])-exp.number-1,lower.tail = FALSE)
    
    
  }
  
  
  ##########
  # Output #
  ##########
  
  cat("Q-Statistic for instrument validity:")
  
  cat("\n")
  
  cat(Q_valid, "on", length(r_input[,2])-exp.number-1, "DF",",", "p-value:" , Q_chiValid)
  
  cat("\n")
  
  multi_return <- function() {
    Out_list <- list("Qstat" = Q_valid, "Qpval"=Q_chiValid)
    
    #Defines class of output object
    
    return(Out_list)
  }
  OUT <- multi_return()
}





