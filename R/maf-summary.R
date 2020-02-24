#' maf.summary
#'
#' Creates a set of plot summarising a maf file.
#'
#' @param maf the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param mut.type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param spe.plat boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is FALSE.
#' @return p.class Barplot of counts of each variant classification
#' @return p.type Barplot of counts of each variant type
#' @return p.SNV Histogram of counts of each SNV class
#' @return p.patient.variant Histogram of counts variants per patient
#' @return p.variant.bp Boxplot of the distribution of variant classification per patient
#' @return p.variant.dist Boxplot of the relative proportion of each variant class in individual patients
#' @return p.variant.dist.bar Stacked barplot of the variants distribution in all patients
#' @return p.SNV.dist Boxplot of the relative proportion of each SNV class in individual patients
#' @return p.corr Correlation heatmap of the top 10 genes
#' @return p.comut Heatmap of the comutation of the top 10 genes
#' @export
#'
#' @examples library(gnomeR)
#' all.plots <- maf.summary(maf=mut)
#' @import
#' dplyr
#' stringr
#' ggplot2
#' reshape2
#' RColorBrewer
#' GGally


maf.summary <- function(maf,mut.type = "SOMATIC", spe.plat = F){

  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
  if(length(match("Variant_Classification",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
  if(length(match("Mutation_Status",colnames(maf))) == 0)
    stop("The MAF file inputted is missing a mutation status column. (Mutation_Status)")

  # recode gene names that have been changed between panel versions to make sure they are consistent and counted as the same gene
  if(!is.character(maf$Hugo_Symbol)) maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)
  if (sum(grepl("KMT2D", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "KMT2D" ~ "MLL2",
        TRUE ~ Hugo_Symbol
      ))

    warning("KMT2D has been recoded to MLL2")
  }

  if (sum(grepl("KMT2C", maf$Hugo_Symbol)) > 1) {
    maf <- maf %>%
      mutate(Hugo_Symbol = case_when(
        Hugo_Symbol == "KMT2C" ~ "MLL3",
        TRUE ~ Hugo_Symbol
      ))

    warning("KMT2C has been recoded to MLL3")
  }

  substrRight <- function(x, n){
    x <- as.character(x)
    substr(x, nchar(x)-n+1, nchar(x))
  }

  # maf <- as.data.frame(maf)
  if(mut.type == "ALL") Mut.filt = unique(maf$Mutation_Status) else Mut.filt = mut.type
  maf <- maf %>% filter(tolower(Mutation_Status) %in% tolower(Mut.filt))
  nb.cols <- 20 #length(unique(maf$Variant_Classification))
  ## summarise variant wise ##

  # variants call summary plot #
  maf$Variant_Classification <- factor(as.character(maf$Variant_Classification),
                                       levels =  as.character(unlist(maf %>%
                                                                       group_by(Variant_Classification) %>%
                                                                       summarise(N = n()) %>%
                                                                       arrange(-N) %>% select(Variant_Classification))))

  p.class <- maf %>% ggplot(aes(x = Variant_Classification,color=Variant_Classification,fill = Variant_Classification)) +
    geom_bar() + coord_flip() + theme(legend.position="none") +
    ggtitle("Variant Classification count") + xlab("Variant Classification") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) )


  # variants type summary plot #
  maf$Variant_Type <- factor(as.character(maf$Variant_Type),
                             levels =  as.character(unlist(maf %>%
                                                             group_by(Variant_Type) %>%
                                                             summarise(N = n()) %>%
                                                             arrange(-N) %>% select(Variant_Type))))
  p.type <- maf %>% ggplot(aes(x = Variant_Type,color=Variant_Type,fill = Variant_Type)) +
    geom_bar() + coord_flip() + theme(legend.position="none") +
    ggtitle("Variant Type count") + xlab("Variant Type")

  # SNV #
  maf <- maf %>%
    filter(Variant_Type == "SNP",
           HGVSc != "") %>%
    mutate(SNV_Class = substrRight(HGVSc,3))
  p.SNV <- maf %>%
    filter(!grepl("N",SNV_Class)) %>%
    ggplot(aes(x = SNV_Class,color=SNV_Class,fill = SNV_Class)) +
    geom_bar() + coord_flip() + theme(legend.position="none") +
    ggtitle("SNV Class count") + xlab("SNV Class") +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) )



  ## summarise patient wise ##
  # distribution of variant per sample #
  maf$Tumor_Sample_Barcode <- factor(as.character(maf$Tumor_Sample_Barcode),levels = as.character(unlist(maf %>%
                                                                                                           group_by(Tumor_Sample_Barcode) %>%
                                                                                                           summarise(N = n()) %>%
                                                                                                           arrange(-N) %>% select(Tumor_Sample_Barcode))))

  p.patient.variant <- maf %>%
    group_by(Tumor_Sample_Barcode,Variant_Classification) %>%
    ggplot(aes(x = Tumor_Sample_Barcode)) + geom_bar(aes(fill = Variant_Classification)) +
    ggtitle("Variants per sample") + ylab("Variant Count") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) )

  # same but in boxplot format # --> simpler to read
  p.variant.bp <- maf %>%
    group_by(Tumor_Sample_Barcode,Variant_Classification) %>%
    summarise(N = n()) %>%
    ggplot(aes(x = Variant_Classification, y = N,color = Variant_Classification)) +
    ggtitle("Variant classification distribution") + xlab("Variant Classification") +
    geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) )


  # most mutated genes and class distrib #
  top.genes <- as.character(unlist(maf %>%
                                     group_by(Hugo_Symbol) %>%
                                     summarise(N = n()) %>%
                                     arrange(-N) %>%
                                     select(Hugo_Symbol)))
  maf$Hugo_Symbol <- factor(as.character(maf$Hugo_Symbol),levels = top.genes)

  p.genes <-  maf %>%
    filter(Hugo_Symbol %in% top.genes[1:10]) %>%
    ggplot(aes(x = Hugo_Symbol)) + geom_bar(aes(fill = Variant_Classification)) +
    coord_flip() +
    ggtitle("Top genes variants classification") + xlab("Gene Name")+
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) )


  # distrib of variants per patient #
  varprop.maf <- maf %>% tbl_df() %>%
    group_by(Tumor_Sample_Barcode) %>%
    mutate(totalMut = n()) %>%
    ungroup() %>%
    group_by(Tumor_Sample_Barcode,Variant_Classification) %>%
    summarise(N=n(),
              varProp = N/unique(totalMut))

  p.variant.dist <- varprop.maf %>%
    ggplot(aes(x = Variant_Classification,y = varProp)) + geom_boxplot(aes(fill = Variant_Classification))+
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="none") + ylab("% Variant")

  p.variant.dist.bar <- varprop.maf %>%
    ggplot(aes(x = Tumor_Sample_Barcode, y=varProp)) +
    geom_bar(aes(fill = Variant_Classification),stat = "identity") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) ) +
    ylab("% Variants")


  # distrib of mutations per patient #
  p.SNV.dist <- maf %>% tbl_df() %>%
    group_by(Tumor_Sample_Barcode) %>%
    mutate(totalMut = n()) %>%
    ungroup() %>%
    group_by(Tumor_Sample_Barcode,SNV_Class) %>%
    summarise(N=n(),
              varProp = N/unique(totalMut)) %>%
    ggplot(aes(x = SNV_Class,y = varProp)) + geom_boxplot(aes(fill = SNV_Class))+
    scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Accent"))(nb.cols) ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(legend.position="none") + ylab("% SNV")


  # comutation patterns #
  bin.maf <- binmat(maf = maf,mut.type = mut.type, spe.plat = spe.plat)
  bin.maf <- bin.maf
  keep <- names(sort(apply(bin.maf,2,function(x){sum(x)}),decreasing = T))[1:10]
  bin.maf <- bin.maf[,keep]
  p.corr <- ggcorr(bin.maf,limits = NULL)


  co.mut <- apply(bin.maf,2,function(x){
    apply(bin.maf,2,function(y){
      sum(y == 1 & x == 1,na.rm = T)/length(x)
    })
  })
  p.comut <- ggcorr(co.mut,limits = NULL)

  return(list("p.class"=p.class,"p.type"=p.type,"p.SNV" = p.SNV,
              "p.patient.variant"=p.patient.variant,"p.variant.bp" = p.variant.bp,"p.genes" = p.genes,
              "p.variant.dist"=p.variant.dist,"p.variant.dist.bar"=p.variant.dist.bar,"p.SNV.dist"=p.SNV.dist,
              "p.corr"=p.corr,"p.comut"=p.comut))
}
