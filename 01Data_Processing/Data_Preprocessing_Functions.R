#Data Preprocessing Functions----------------------------------------------------


read_exposure_data <- function(filename, clump=FALSE, sep=" ", phenotype_col="Phenotype", snp_col="SNP", beta_col="beta", se_col="se", eaf_col="eaf", effect_allele_col="effect_allele", other_allele_col="other_allele", pval_col="pval", units_col="units", ncase_col="ncase", ncontrol_col="ncontrol", samplesize_col="samplesize", gene_col="gene", id_col="id", min_pval=1e-200, log_pval=FALSE, chr_col="chr", pos_col="pos")
{
  exposure_dat <- data.table::fread(filename, header=TRUE, sep=sep)
  exposure_dat <- format_data(
    as.data.frame(exposure_dat),
    type="exposure",
    snps=NULL,
    phenotype_col=phenotype_col,
    snp_col=snp_col,
    beta_col=beta_col,
    se_col=se_col,
    eaf_col=eaf_col,
    effect_allele_col=effect_allele_col,
    other_allele_col=other_allele_col,
    pval_col=pval_col,
    units_col=units_col,
    ncase_col=ncase_col,
    ncontrol_col=ncontrol_col,
    samplesize_col=samplesize_col,
    gene_col=gene_col,
    id_col=id_col,
    min_pval=min_pval,
    log_pval=log_pval,
    chr_col=chr_col,
    pos_col=pos_col
  )
  exposure_dat$data_source.exposure <- "textfile"
  if(clump)
  {
    exposure_dat <- clump_data(exposure_dat)
  }
  return(exposure_dat)
}



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
  
  #If no pval column then create it from beta and se if available.
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



ld_clump <- function(dat = NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99, pop = "EUR", opengwas_jwt = get_opengwas_jwt(), bfile = NULL, plink_bin = NULL)
{
  
  stopifnot("rsid" %in% names(dat))
  stopifnot(is.data.frame(dat))
  
  if(is.null(bfile))
  {
    message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
  }
  
  if(! "pval" %in% names(dat))
  {
    if( "p" %in% names(dat))
    {
      warning("No 'pval' column found in dat object. Using 'p' column.")
      dat[["pval"]] <- dat[["p"]]
    } else {
      warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
      dat[["pval"]] <- clump_p
    }
  }
  
  if(! "id" %in% names(dat))
  {
    dat$id <- random_string(1)
  }
  
  ids <- unique(dat[["id"]])
  res <- list()
  for(i in 1:length(ids))
  {
    x <- subset(dat, dat[["id"]] == ids[i])
    if(nrow(x) == 1)
    {
      message("Only one SNP for ", ids[i])
      res[[i]] <- x
    } else {
      message("Clumping ", ids[i], ", ", nrow(x), " variants, using ", pop, " population reference")
      if(is.null(bfile))
      {
        res[[i]] <- ld_clump_api(x, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, pop=pop, opengwas_jwt=opengwas_jwt)
      } else {
        res[[i]] <- ld_clump_local(x, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin=plink_bin)
      }
    }
  }
  res <- dplyr::bind_rows(res)
  return(res)
}



get_plink_exe <- function()
{
  os <- Sys.info()['sysname']
  a <- paste0("bin/plink_", os)
  if(os == "Windows") a <- paste0(a, ".exe")
  plink_bin <- system.file(a, package="plinkbinr")
  if(!file.exists(plink_bin))
  {
    stop("No plink2 executable available for OS '", os, "'. Please provide your own plink2 executable file using the plink_bin argument.")
  }
  return(plink_bin)
}




