#' gen.tab
#'
#' Creates a summary table of the distribution of the genetic features by a specific outcome/covariate of interest. The difference
#' is tested using Fisher's exact test and further adjusted for multiple comparisons. Note that continuous genetic factors
#' are dichotomized at their median.
#'
#' @param gen.dat A matrix or dataframe, with patients as rows and features as columns.
#' @param outcome A leveled vector of length equal to the number of rows in gen.dat.
#' @param filter a numeric value between 0 and 1 (1 not included) that is the lower bound for the proportion of patients
#' having a genetic event (only for binary features). All features with an event rate lower than that value will be removed.
#' Default is 0 (all features included).
#' @param paired Boolean if the data are paired. Default is FALSE.
#' @param cont Should the outcome be treated as a continuous value. Default is FALSE treated as categorical.
#' @param rank Should the table returned be oredered by Pvalue. Boolean, default is T
#' @return fits : a table of odds ratio and pvalues.
#' @return forest.plot : A forest plot of the top 10 hits.
#'
#' @export
#'
#' @examples library(gnomeR)
#' patients <- unique(clin$DMPID)[1:500]
#' mut.only <- create.bin.matrix(patients = patients,maf = mut)
#' gen.dat <- mut.only$mut
#' outcome = as.factor(as.character(clin$Sex[match(patients,clin$DMPID)]))
#' test <- gen.tab(gen.dat,outcome)
#' head(test$fits)
#' test$forest.plot

gen.tab <- function(gen.dat,outcome,filter=0,paired = F,cont=F,rank = T){

  # remove all columns that are constant #
  if(length(which(apply(gen.dat, 2, function(x){length(unique(x[!is.na(x)])) == 1} || all(is.na(x))))) > 0)
    gen.dat <- gen.dat[,-which(apply(gen.dat, 2, function(x){length(unique(x[!is.na(x)])) <= 1 || all(is.na(x))}))]

  if(filter > 0){
    # get binary cases #
    temp <- apply(gen.dat, 2, function(x){length(unique(x[!is.na(x)])) == 2})
    genes.bin <- names(temp[which(temp)])
    if(length(genes.bin) == ncol(gen.dat)) rm <- apply(gen.dat, 2, function(x){sum(x, na.rm = T)/length(x) < filter})
    else rm <- apply(gen.dat[,genes.bin], 2, function(x){sum(x, na.rm = T)/length(x) < filter})
    genes.rm <- names(rm[which(rm)])
    # print(genes.rm)
    gen.dat <- gen.dat %>%
      select(-one_of(genes.rm))
  }
  if(is.null(dim(gen.dat)) )
    stop("Only one or fewer genes are left after filtering. We need a minimum of two. Please relax the filter argument.")


  if(!cont){
    if(is.character(outcome)) outcome <- as.factor(outcome)

    fits <- as.data.frame(t(apply(gen.dat, 2, function(x) {
      if (length(unique(x))/length(x) > 0.5) {
        x <- ifelse(x > median(x, na.rm = T), 1, 0)
      }
      if(paired == F) test <- fisher.test(x, outcome)
      if(paired == T){
        # tt = with(as.data.frame(cbind(x,outcome)), table(x,outcome))
        test <- exact2x2::mcnemar.exact(x = x[1:(length(x)/2)],y = x[(length(x)/2+1):length(x)])
      }
      out <- c(sum(x, na.rm = T)/length(x))
      for (i in 1:length(levels(outcome))) {
        out <- c(out, sum(x[which(outcome == levels(outcome)[i])], na.rm = T)/length(which(outcome == levels(outcome)[i])))
      }
      out <- paste0(round(as.numeric(out)*100,digits = 2),"%")
      if (!is.null(test$estimate)) {
        out <- c(out, round(as.numeric(test$estimate),digits = 2), formatC(test$p.value, format = "e", digits = 2), round(as.numeric(test$conf.int),digits =2))
        names(out) <- c("Overall",levels(outcome)[1:length(levels(outcome))],
                        "OddsRatio", "Pvalue", "Lower", "Upper")
      }
      else {
        out <- c(out, test$p.value, round(as.numeric(test$conf.int),digits =2))
        names(out) <- c("Overall",levels(outcome)[1:length(levels(outcome))],
                        "Pvalue")
      }
      return(out)
    })))

    colnames(fits)[2:(length(levels(outcome))+1)] <- paste0(colnames(fits)[2:(length(levels(outcome))+1)],
                                                            "(N=",as.numeric(summary(outcome)),")")
    fits$FDR <- formatC(p.adjust(as.numeric(as.character(fits$Pvalue)),method="fdr"), format = "e", digits = 2)
    if(rank) fits <- fits[order(as.numeric(as.character(fits$Pvalue))),]

    if (!is.null(fits$OddsRatio)){
      # forest plot #
      f.dat <- fits
      f.dat$Gene <- rownames(f.dat)
      f.dat <- f.dat %>%
        filter(!is.infinite(OddsRatio)) %>%
        mutate(
          OddsRatio = as.numeric(as.character(OddsRatio)),
          Lower = as.numeric(as.character(Lower)),
          Upper = as.numeric(as.character(Upper)),
        )
      f.dat <- f.dat[1:min(10,nrow(f.dat)),]
      forest.plot <- f.dat %>%
        ggplot(aes(x = Gene,y = OddsRatio, ymin = Lower, ymax = Upper ))+
        geom_pointrange(aes(col=Gene))+
        geom_hline(aes(fill=Gene),yintercept =1, linetype=2)+
        xlab('Gene')+ ylab("Risk Ratio (95% Confidence Interval)")+
        geom_errorbar(aes(ymin=Lower, ymax=Upper,col=Gene),width=0.5,cex=1)+
        coord_flip()
    }
    else{
      forest.plot <- NULL
    }
    return(list("fits"=fits,"forest.plot"=forest.plot))
  }


  ############


  if(cont){
    fits <- as.data.frame(do.call('rbind',apply(gen.dat, 2, function(x) {
      if(length(unique(x))/length(x) > 0.25 || length(unique(x)) == 2){
        fit <- lm(outcome ~ x)
        temp <- as.data.frame(summary(fit)$coefficient)
        colnames(temp) <-  c("Estimate","SD", "tvalue","Pvalue")
        if(is.numeric(x)) temp$MutationFreq <- sum(x,na.rm = T)/length(x[!is.na(x)])
        else temp$MutationFreq <- 0
        out <- temp[2,c(1,2,4,5)]
      }
      else{
        fit <- lm(outcome ~ as.factor(x))
        temp <- as.data.frame(summary(fit)$coefficient)
        colnames(temp) <-  c("Estimate","SD", "tvalue","Pvalue")
        # if(is.numeric(x)) temp$MutationFreq <- sum(x,na.rm = T)/length(x[!is.na(x)])
        # else
        temp$MutationFreq <- 0#rep(0,nrow(temp))
        out <- as.data.frame(temp[2:nrow(temp),c(1,2,4,5)])
        rownames(out) <- gsub("as.factor\\(x\\)","",rownames(out))
      }
      return(out)
    })))


    fits$FDR <- p.adjust(fits$Pvalue,method = "fdr")
    fits$GeneName <- rownames(fits)

    if(all(apply(gen.dat,2,is.numeric))){
      vPlot <- try(plot_ly(data = fits %>% filter(!is.na(Pvalue),is.numeric(Pvalue)), x = ~Estimate, y = ~-log10(Pvalue),
                       text = ~paste('Gene :',GeneName,
                                     '</br> Estimate :',round(Estimate,digits=2)),
                       mode = "markers",size = ~MutationFreq,color = ~Estimate) %>%
        layout(title ="Volcano Plot"),silent = T)
    }
    else{
      vPlot <- try(plot_ly(data = fits %>% filter(!is.na(Pvalue),is.numeric(Pvalue)), x = ~Estimate, y = ~-log10(Pvalue),
                       text = ~paste('Gene :',GeneName,
                                     '</br> Estimate :',round(Estimate,digits=2)),
                       mode = "markers") %>%
        layout(title ="Volcano Plot"),silent =T)
    }
    # fits <- fits[,-match("GeneName",colnames(fits))]
    if(rank) fits <- fits[order(fits$Pvalue),-ncol(fits)]
    return(list("fits"=fits,"vPlot"=vPlot))
  }
}
