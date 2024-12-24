#Result Visualization------------------------------------------------------------


#Density Plot visualization------------------------------------------------------

mr_density_plot <- function(singlesnp_results, mr_results, exponentiate=FALSE, bandwidth="nrd0")
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    d$SNP <- as.character(d$SNP)
    
    d2 <- subset(d, !grepl("All - ", SNP))
    d1 <- subset(mr_results, id.exposure == d2$id.exposure[1] & id.outcome == d2$id.outcome[1])
    
    xint <- 0
    if(exponentiate)
    {
      d$b <- exp(d$b)
      d$up <- exp(d$up)
      d$lo <- exp(d$lo)
      xint <- 1
    }
    
    ggplot2::ggplot(d2, ggplot2::aes(x=b)) +
      ggplot2::geom_vline(xintercept=xint, linetype="dotted") +
      ggplot2::geom_density(ggplot2::aes(weight=1/se), bw=bandwidth) +
      ggplot2::geom_point(y=0, colour="red", ggplot2::aes(size=1/se)) +
      ggplot2::geom_vline(data=mr_results, ggplot2::aes(xintercept=b, colour=method)) +
      ggplot2::scale_colour_brewer(type="qual") +
      ggplot2::labs(x="Per SNP MR estimate")
  })
  res
}
