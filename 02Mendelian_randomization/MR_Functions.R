#MR Functions--------------------------------------------------------------------


format_data <- function(dat, type="exposure", snps=NULL, header=TRUE,
                        phenotype_col="Phenotype", snp_col="SNP",
                        beta_col="beta", se_col="se", eaf_col="eaf",
                        effect_allele_col="effect_allele",
                        other_allele_col="other_allele", pval_col="pval",
                        units_col="units", ncase_col="ncase",
                        ncontrol_col="ncontrol", samplesize_col="samplesize",
                        gene_col="gene", id_col="id", min_pval=1e-200,
                        z_col="z", info_col="info", chr_col="chr",
                        pos_col="pos", log_pval=FALSE)
{
  
  if (inherits(dat, "data.table")) {
    datname <- deparse(substitute(dat))
    stop(paste0(
      "Your ", datname, " data.frame is also of class 'data.table', ",
      "please reformat as simply a data.frame with ", datname, " <- data.frame(",
      datname, ") and then rerun your format_data() call."
    ))
  }
  
  all_cols <- c(phenotype_col, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, units_col, ncase_col, ncontrol_col, samplesize_col, gene_col, id_col, z_col, info_col, chr_col, pos_col)
  
  i <- names(dat) %in% all_cols
  if(sum(i) == 0)
  {
    stop("None of the specified columns present")
  }
  dat <- dat[,i]
  
  if(! snp_col %in% names(dat))
  {
    stop("SNP column not found")
  }
  
  names(dat)[names(dat) == snp_col] <- "SNP"
  snp_col <- "SNP"
  dat$SNP <- tolower(dat$SNP)
  dat$SNP <- gsub("[[:space:]]", "", dat$SNP)
  dat <- subset(dat, !is.na(SNP))
  
  if(!is.null(snps))
  {
    dat <- subset(dat, SNP %in% snps)
  }
  
  if(! phenotype_col %in% names(dat))
  {
    message("No phenotype name specified, defaulting to '", type, "'.")
    dat[[type]] <- type
  } else {
    dat[[type]] <- dat[[phenotype_col]]
    if(phenotype_col != type)
    {
      dat <- dat[,-which(names(dat)==phenotype_col)]
    }
  }
  
  if(log_pval)
  {
    dat$pval <- 10^-dat[[pval_col]]
  }
  
  #Remove duplicated SNPs
  dat <- plyr::ddply(dat, type, function(x) {
    x <- plyr::mutate(x)
    dup <- duplicated(x$SNP)
    if(any(dup))
    {
      warning("Duplicated SNPs present in exposure data for phenotype '", x[[type]][1], ". Just keeping the first instance:\n", paste(x$SNP[dup], collapse="\n"))
      x <- x[!dup,]
    }
    return(x)
  })
  
  #Check if columns required for MR are present
  mr_cols_required <- c(snp_col, beta_col, se_col, effect_allele_col)
  mr_cols_desired <- c(other_allele_col, eaf_col)
  if(! all(mr_cols_required %in% names(dat)))
  {
    warning("The following columns are not present and are required for MR analysis\n", paste(mr_cols_required[!mr_cols_required %in% names(dat)]), collapse="\n")
    dat$mr_keep.outcome <- FALSE
  } else {
    dat$mr_keep.outcome <- TRUE
  }
  
  if(! all(mr_cols_desired %in% names(dat)))
  {
    warning("The following columns are not present but are helpful for harmonisation\n", paste(mr_cols_desired[!mr_cols_desired %in% names(dat)]), collapse="\n")
  }
  
  #Check beta
  i <- which(names(dat) == beta_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "beta.outcome"
    if(!is.numeric(dat$beta.outcome))
    {
      warning("beta column is not numeric. Coercing...")
      dat$beta.outcome <- as.numeric(dat$beta.outcome)
    }
    index <- !is.finite(dat$beta.outcome)
    index[is.na(index)] <- TRUE
    dat$beta.outcome[index] <- NA
  }
  
  #Check se
  i <- which(names(dat) == se_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "se.outcome"
    if(!is.numeric(dat$se.outcome))
    {
      warning("se column is not numeric. Coercing...")
      dat$se.outcome <- as.numeric(dat$se.outcome)
    }
    index <- !is.finite(dat$se.outcome) | dat$se.outcome <= 0
    index[is.na(index)] <- TRUE
    dat$se.outcome[index] <- NA
  }
  
  #Check eaf
  i <- which(names(dat) == eaf_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "eaf.outcome"
    if(!is.numeric(dat$eaf.outcome))
    {
      warning("eaf column is not numeric. Coercing...")
      dat$eaf.outcome <- as.numeric(dat$eaf.outcome)
    }
    index <- !is.finite(dat$eaf.outcome) | dat$eaf.outcome <= 0 | dat$eaf.outcome >= 1
    index[is.na(index)] <- TRUE
    dat$eaf.outcome[index] <- NA
  }
  
  #Check effect_allele
  i <- which(names(dat) == effect_allele_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "effect_allele.outcome"
    if(is.logical(dat$effect_allele.outcome))
    {
      dat$effect_allele.outcome <- substr(as.character(dat$effect_allele.outcome), 1, 1)
    }
    if(!is.character(dat$effect_allele.outcome))
    {
      warning("effect_allele column is not character data. Coercing...")
      dat$effect_allele.outcome <- as.character(dat$effect_allele.outcome)
    }
    
    dat$effect_allele.outcome <- toupper(dat$effect_allele.outcome)
    index <- ! (grepl("^[ACTG]+$", dat$effect_allele.outcome) | dat$effect_allele.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if(any(index))
    {
      warning("effect_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded.")
      dat$effect_allele.outcome[index] <- NA
      dat$mr_keep.outcome[index] <- FALSE
    }
  }
  
  
  #Check other_allele
  i <- which(names(dat) == other_allele_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "other_allele.outcome"
    if(is.logical(dat$other_allele.outcome))
    {
      dat$other_allele.outcome <- substr(as.character(dat$other_allele.outcome), 1, 1)
    }
    if(!is.character(dat$other_allele.outcome))
    {
      warning("other_allele column is not character data. Coercing...")
      dat$other_allele.outcome <- as.character(dat$other_allele.outcome)
    }
    
    dat$other_allele.outcome <- toupper(dat$other_allele.outcome)
    index <- ! (grepl("^[ACTG]+$", dat$other_allele.outcome) | dat$other_allele.outcome %in% c("D", "I"))
    index[is.na(index)] <- TRUE
    if(any(index))
    {
      warning("other_allele column has some values that are not A/C/T/G or an indel comprising only these characters or D/I. These SNPs will be excluded")
      dat$other_allele.outcome[index] <- NA
      dat$mr_keep.outcome[index] <- FALSE
    }
  }
  
  
  #Check pval
  i <- which(names(dat) == pval_col)[1]
  if(!is.na(i))
  {
    names(dat)[i] <- "pval.outcome"
    if(!is.numeric(dat$pval.outcome))
    {
      warning("pval column is not numeric. Coercing...")
      dat$pval.outcome <- as.numeric(dat$pval.outcome)
    }
    index <- !is.finite(dat$pval.outcome) | dat$pval.outcome < 0 | dat$pval.outcome > 1
    index[is.na(index)] <- TRUE
    dat$pval.outcome[index] <- NA
    index <- dat$pval.outcome < min_pval
    index[is.na(index)] <- FALSE
    dat$pval.outcome[index] <- min_pval
    
    dat$pval_origin.outcome <- "reported"
    if(any(is.na(dat$pval.outcome)))
    {
      if("beta.outcome" %in% names(dat) && "se.outcome" %in% names(dat))
      {
        index <- is.na(dat$pval.outcome)
        dat$pval.outcome[index] <- stats::pnorm(abs(dat$beta.outcome[index])/dat$se.outcome[index], lower.tail=FALSE)
        dat$pval_origin.outcome[index] <- "inferred"
      }
    }
  }
  
  #If no pval column then create it from beta and se if available
  if("beta.outcome" %in% names(dat) && "se.outcome" %in% names(dat) && ! "pval.outcome" %in% names(dat))
  {
    message("Inferring p-values")
    dat$pval.outcome <- stats::pnorm(abs(dat$beta.outcome)/dat$se.outcome, lower.tail=FALSE) * 2
    dat$pval_origin.outcome <- "inferred"
  }
  
  
  if(ncase_col %in% names(dat))
  {
    names(dat)[which(names(dat) == ncase_col)[1]] <- "ncase.outcome"
    if(!is.numeric(dat$ncase.outcome))
    {
      warning(ncase_col, " column is not numeric")
      dat$ncase.outcome <- as.numeric(dat$ncase.outcome)
    }
  }
  if(ncontrol_col %in% names(dat))
  {
    names(dat)[which(names(dat) == ncontrol_col)[1]] <- "ncontrol.outcome"
    if(!is.numeric(dat$ncontrol.outcome))
    {
      warning(ncontrol_col, " column is not numeric")
      dat$ncontrol.outcome <- as.numeric(dat$ncontrol.outcome)
    }
  }
  
  
  
  if(samplesize_col %in% names(dat))
  {
    names(dat)[which(names(dat) == samplesize_col)[1]] <- "samplesize.outcome"
    if(!is.numeric(dat$samplesize.outcome))
    {
      warning(samplesize_col, " column is not numeric")
      dat$samplesize.outcome <- as.numeric(dat$samplesize.outcome)
    }
    
    if("ncontrol.outcome" %in% names(dat) && "ncase.outcome" %in% names(dat))
    {
      index <- is.na(dat$samplesize.outcome) & !is.na(dat$ncase.outcome) & !is.na(dat$ncontrol.outcome)
      if(any(index))
      {
        message("Generating sample size from ncase and ncontrol")
        dat$samplesize.outcome[index] <- dat$ncase.outcome[index] + dat$ncontrol.outcome[index]
      }
    }
  } else if("ncontrol.outcome" %in% names(dat) && "ncase.outcome" %in% names(dat))
  {
    message("Generating sample size from ncase and ncontrol")
    dat$samplesize.outcome <- dat$ncase.outcome + dat$ncontrol.outcome
  }
  
  if(gene_col %in% names(dat))
  {
    names(dat)[which(names(dat) == gene_col)[1]] <- "gene.outcome"
  }
  
  if(info_col %in% names(dat))
  {
    names(dat)[which(names(dat) == info_col)[1]] <- "info.outcome"
  }
  
  if(z_col %in% names(dat))
  {
    names(dat)[which(names(dat) == z_col)[1]] <- "z.outcome"
  }
  
  if(chr_col %in% names(dat))
  {
    names(dat)[which(names(dat) == chr_col)[1]] <- "chr.outcome"
  }
  
  if(pos_col %in% names(dat))
  {
    names(dat)[which(names(dat) == pos_col)[1]] <- "pos.outcome"
  }
  
  if(units_col %in% names(dat))
  {
    names(dat)[which(names(dat) == units_col)[1]] <- "units.outcome"
    dat$units.outcome_dat <- as.character(dat$units.outcome)
    temp <- check_units(dat, type, "units.outcome")
    if(any(temp$ph))
    {
      dat[[type]] <- paste0(dat[[type]], " (", dat$units.outcome, ")")
    }
  }
  
  #Create id column
  if(id_col %in% names(dat))
  {
    names(dat)[which(names(dat) == id_col)[1]] <- "id.outcome"
    dat$id.outcome <- as.character(dat$id.outcome)
  } else {
    dat$id.outcome <- create_ids(dat[[type]])
  }
  
  if(any(dat$mr_keep.outcome))
  {
    mrcols <- c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome")
    mrcols_present <- mrcols[mrcols %in% names(dat)]
    dat$mr_keep.outcome <- dat$mr_keep.outcome & apply(dat[, mrcols_present], 1, function(x) !any(is.na(x)))
    if(any(!dat$mr_keep.outcome))
    {
      warning("The following SNP(s) are missing required information for the MR tests and will be excluded\n", paste(subset(dat, !mr_keep.outcome)$SNP, collapse="\n"))
    }
  }
  if(all(!dat$mr_keep.outcome))
  {
    warning("None of the provided SNPs can be used for MR analysis, they are missing required information.")
  }
  
  #Add in missing MR cols
  for(col in c("SNP", "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"))
  {
    if(! col %in% names(dat))
    {
      dat[[col]] <- NA
    }
  }
  
  names(dat) <- gsub("outcome", type, names(dat))
  rownames(dat) <- NULL
  return(dat)
}



extract_outcome_data <- function(snps, outcomes, proxies = TRUE, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3, opengwas_jwt=ieugwasr::get_opengwas_jwt(), splitsize=10000, proxy_splitsize=500)
{
  outcomes <- ieugwasr::legacy_ids(unique(outcomes))
  
  snps <- unique(snps)
  firstpass <- extract_outcome_data_internal(snps, outcomes, proxies = FALSE, opengwas_jwt=opengwas_jwt, splitsize = splitsize)
  
  if(proxies)
  {
    for(i in seq_along(outcomes))
    {
      if(is.null(firstpass))
      {
        missedsnps <- snps
      } else {
        missedsnps <- snps[!snps %in% subset(firstpass, id.outcome == outcomes[i])$SNP]
      }
      if(length(missedsnps)>0)
      {
        message("Finding proxies for ", length(missedsnps), " SNPs in outcome ", outcomes[i])
        temp <- extract_outcome_data_internal(missedsnps, outcomes[i], proxies = TRUE, rsq, align_alleles, palindromes, maf_threshold, opengwas_jwt = opengwas_jwt, splitsize = proxy_splitsize)
        if(!is.null(temp))
        {
          firstpass <- plyr::rbind.fill(firstpass, temp)
        }
      }
    }
  }
  
  return(firstpass)
}



harmonise_data <- function(exposure_dat, outcome_dat, action=2)
{
  stopifnot(all(action %in% 1:3))
  check_required_columns(exposure_dat, "exposure")
  check_required_columns(outcome_dat, "outcome")
  res.tab <- merge(outcome_dat, exposure_dat, by="SNP")
  ncombinations <- length(unique(res.tab$id.outcome))
  if(length(action) == 1)
  {
    action <- rep(action, ncombinations)
  } else if(length(action) != ncombinations) {
    stop("Action argument must be of length 1 (where the same action will be used for all outcomes), or number of unique id.outcome values (where each outcome will use a different action value)")
  }
  
  res.tab <- harmonise_cleanup_variables(res.tab)
  if (nrow(res.tab) == 0) {
    # NOTE: currently the columns are not in the same order as they would
    # be should there be overlapping variants
    return(res.tab)
  }
  
  d <- data.frame(id.outcome=unique(res.tab$id.outcome), action=action)
  res.tab <- merge(res.tab, d, by="id.outcome")
  
  combs <- subset(res.tab, !duplicated(paste(id.exposure, id.outcome)), select=c(id.exposure, id.outcome))
  
  fix.tab <- list()
  mr_cols <- c("beta.exposure", "beta.outcome", "se.exposure", "se.outcome")
  for(i in seq_len(nrow(combs)))
  {
    x <- subset(res.tab, id.exposure == combs$id.exposure[i] & id.outcome == combs$id.outcome[i])
    message("Harmonising ", x$exposure[1], " (", x$id.exposure[1], ") and ", x$outcome[1], " (", x$id.outcome[1], ")")
    x <- harmonise(x, 0.08, x$action[1])
    attr(x, "log")[['candidate_variants']] <- sum(exposure_dat$id.exposure == x$id.exposure[1])
    attr(x, "log")[['variants_absent_from_reference']] <- sum(exposure_dat$id.exposure == x$id.exposure[1]) - nrow(x)
    
    x$mr_keep[apply(x[, mr_cols], 1, function(y) any(is.na(y)))] <- FALSE
    attr(x, "log")[["total_variants"]] <- nrow(x)
    attr(x, "log")[["total_variants_for_mr"]] <- sum(x$mr_keep)
    attr(x, "log")[["proxy_variants"]] <- ifelse(is.null(x$proxy.outcome), 0, sum(x$proxy.outcome, na.rm=TRUE))
    fix.tab[[i]] <- x
  }

  
  jlog <- plyr::rbind.fill(lapply(fix.tab, function(x) attr(x, "log")))
  fix.tab <- plyr::rbind.fill(fix.tab)
  attr(fix.tab, "log") <- jlog
  
  if(!"samplesize.outcome" %in% names(fix.tab))
  {
    fix.tab$samplesize.outcome <- NA
  }
  
  
  return(fix.tab)
}



mv_harmonise_data <- function(exposure_dat, outcome_dat, harmonise_strictness=2)
{
  
  stopifnot(all(c("SNP", "id.exposure", "exposure", "effect_allele.exposure", "beta.exposure", "se.exposure", "pval.exposure") %in% names(exposure_dat)))
  nexp <- length(unique(exposure_dat$id.exposure))
  stopifnot(nexp > 1)
  tab <- table(exposure_dat$SNP)
  keepsnp <- names(tab)[tab == nexp]
  exposure_dat <- subset(exposure_dat, SNP %in% keepsnp)
  
  
  exposure_mat <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="beta.exposure")
  
  
  #Get outcome data
  dat <- harmonise_data(subset(exposure_dat, id.exposure == exposure_dat$id.exposure[1]), outcome_dat, action=harmonise_strictness)
  dat <- subset(dat, mr_keep)
  dat$SNP <- as.character(dat$SNP)
  
  exposure_beta <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="beta.exposure")
  exposure_beta <- subset(exposure_beta, SNP %in% dat$SNP)
  exposure_beta$SNP <- as.character(exposure_beta$SNP)
  
  exposure_pval <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="pval.exposure")
  exposure_pval <- subset(exposure_pval, SNP %in% dat$SNP)
  exposure_pval$SNP <- as.character(exposure_pval$SNP)
  
  exposure_se <- reshape2::dcast(exposure_dat, SNP ~ id.exposure, value.var="se.exposure")
  exposure_se <- subset(exposure_se, SNP %in% dat$SNP)
  exposure_se$SNP <- as.character(exposure_se$SNP)
  
  index <- match(exposure_beta$SNP, dat$SNP)
  dat <- dat[index, ]
  stopifnot(all(dat$SNP == exposure_beta$SNP))
  
  exposure_beta <- as.matrix(exposure_beta[,-1])
  exposure_pval <- as.matrix(exposure_pval[,-1])
  exposure_se <- as.matrix(exposure_se[,-1])
  
  rownames(exposure_beta) <- dat$SNP
  rownames(exposure_pval) <- dat$SNP
  rownames(exposure_se) <- dat$SNP
  
  outcome_beta <- dat$beta.outcome
  outcome_se <- dat$se.outcome
  outcome_pval <- dat$pval.outcome
  
  expname <- subset(exposure_dat, !duplicated(id.exposure), select=c(id.exposure, exposure))
  outname <- subset(outcome_dat, !duplicated(id.outcome), select=c(id.outcome, outcome))
  
  
  return(list(exposure_beta=exposure_beta, exposure_pval=exposure_pval, exposure_se=exposure_se, outcome_beta=outcome_beta, outcome_pval=outcome_pval, outcome_se=outcome_se, expname=expname, outname=outname))
}















