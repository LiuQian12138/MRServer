#Result Visualization------------------------------------------------------------


#Forest Plot Visualization-------------------------------------------------------

ForsetPlot_lq <- function(singlesnp_results, exponentiate=FALSE)
{
  
  d <- plyr::mutate(singlesnp_results)
  if(sum(!grepl("All", d$SNP)) < 2) {
    return(
      blank_plot("Insufficient number of SNPs")
    )
  }
  levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
  levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
  am <- grep("All", d$SNP, value=TRUE)
  d$up <- d$b + 1.96 * d$se
  d$lo <- d$b - 1.96 * d$se
  d$tot <- 0.01
  d$tot[d$SNP %in% am] <- 1
  d$SNP <- as.character(d$SNP)
  nom <- d$SNP[! d$SNP %in% am]
  nom <- nom[order(d$b)]
  d <- rbind(d, d[nrow(d),])
  d$SNP[nrow(d)-1] <- ""
  d$b[nrow(d)-1] <- NA
  d$up[nrow(d)-1] <- NA
  d$lo[nrow(d)-1] <- NA
  d$SNP <- ordered(d$SNP, levels=c(am, "", nom))
  
  xint <- 0
  if(exponentiate)
  {
    d$b <- exp(d$b)
    d$up <- exp(d$up)
    d$lo <- exp(d$lo)
    xint <- 1
  }
  
  d
}



mr_forest_plot <- function(singlesnp_results, exponentiate=FALSE)
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    levels(d$SNP)[levels(d$SNP) == "All - Inverse variance weighted"] <- "All - IVW"
    levels(d$SNP)[levels(d$SNP) == "All - MR Egger"] <- "All - Egger"
    am <- grep("All", d$SNP, value=TRUE)
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 0.01
    d$tot[d$SNP %in% am] <- 1
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[! d$SNP %in% am]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c(am, "", nom))
    
    xint <- 0
    if(exponentiate)
    {
      d$b <- exp(d$b)
      d$up <- exp(d$up)
      d$lo <- exp(d$lo)
      xint <- 1
    }
    
    ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      ggplot2::geom_vline(xintercept=xint, linetype="dotted") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, linewidth=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
      ggplot2::scale_colour_manual(values=c("black", "red")) +
      ggplot2::scale_linewidth_manual(values=c(0.3, 1)) +
      ggplot2::theme(
        legend.position="none",
        axis.text.y=ggplot2::element_text(size=8),
        axis.ticks.y=ggplot2::element_line(linewidth=0),
        axis.title.x=ggplot2::element_text(size=8)) +
      ggplot2::labs(y="", x=paste0("MR effect size for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))
  })
  res
}
