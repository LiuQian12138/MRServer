#Result Visualization------------------------------------------------------------


#Scatter Plot Visualization------------------------------------------------------

perform_mr_analysis <- function(har_mod, mr_mod) {
  d <- har_mod
  mr_results <- mr_mod
  d <- plyr::mutate(d)
  
  if (nrow(d) < 2 | sum(d$mr_keep) == 0) {
    print("Insufficient number of SNPs")
    return(list(d = d, mrres = NULL))
  } else {
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    
    if ("MR Egger" %in% mrres$method) {
      temp <- mr_egger_regression(
        d$beta.exposure,
        d$beta.outcome,
        d$se.exposure,
        d$se.outcome,
        default_parameters()
      )
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    mrres$x1 <- 0
    mrres$y1 <- mrres$a
    mrres$x2 <- max(d$beta.exposure)+0.01
    mrres$y2 <- mrres$b*mrres$x2+mrres$a
    
    #Generate a data frame to store the results of the scatter plot.
    result_scatter1 <- data.frame(SNP = paste("'",d$SNP,"'",sep = ""),beta_outcome=d$beta.outcome,beta_exposure=d$beta.exposure,left=d$beta.exposure-d$se.exposure,right=d$beta.exposure+d$se.exposure,low=d$beta.outcome-d$se.outcome,high=d$beta.outcome+d$se.outcome)
    
    result_scatter1 <- mutate_all(result_scatter1, as.character)
    
    #Use the dplyr package to manipulate the data frame.
    result1 <- result_scatter1 %>%rowwise() %>% mutate(NewCol = paste(c("[",SNP,",",beta_outcome,",",beta_exposure,",",left,",",right,",",low,",",high , "]"), collapse = "")) %>% select(NewCol)
    
    
    result_scatter2 <- data.frame(method = mrres$method,x1=mrres$x1,x2=mrres$x2,y1=mrres$y1,y2=mrres$y2)
    result2 <-result_scatter2 %>% rowwise() %>% mutate(NewCol = paste("{type: 'line', name: '", method, "', symbolSize:0, data: [[", x1, ", ", y1, "],[", x2, ", ", y2, "]]}", collapse = "",sep = "")) %>% pull(NewCol) %>% paste(collapse = ",")
    
    
    return(list(result1 = result1, result2 = result2))
  }
}



mr_scatter_plot <- function(mr_results, dat)
{
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position="top", legend.direction="vertical") +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
  })
  mrres
}


