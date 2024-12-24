#Result Visualization------------------------------------------------------------


#Funnel Plot Visualization-------------------------------------------------------

mr_funnel_plot <- function(singlesnp_results)
{
  res <- plyr::dlply(singlesnp_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(sum(!grepl("All", d$SNP)) < 2) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    am <- grep("All", d$SNP, value=TRUE)
    d$SNP <- gsub("All - ", "", d$SNP)
    am <- gsub("All - ", "", am)
    ggplot2::ggplot(subset(d, ! SNP %in% am), ggplot2::aes(y = 1/se, x=b)) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(data=subset(d, SNP %in% am), ggplot2::aes(xintercept=b, colour = SNP)) +
      ggplot2::scale_colour_manual(values = c("#a6cee3",
                                              "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
                                              "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",
                                              "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(y=expression(1/SE[IV]), x=expression(beta[IV]), colour="MR Method") +
      ggplot2::theme(legend.position="top", legend.direction="vertical")
  })
  res
}

