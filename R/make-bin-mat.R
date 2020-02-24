#' binmat
#'
#' Enables creation of a binary matrix from a maf file with
#' a predifined list of patients (rows are patients and columns are genes)
#'
#' @param patients a character vector that let's the user specify the patients to be used to create the matrix.
#' Default is NULL is which case all patients in the MAF file will be used.
#' @param maf A MAF file.
#' @param mut.type The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
#' keep all mutations regardless of status (not recommended). Default is SOMATIC.
#' @param SNP.only Boolean to rather the genetics events to be kept only to be SNPs (insertions and deletions will be removed).
#' Default is FALSE.
#' @param include.silent Boolean to keep or remove all silent mutations. TRUE keeps, FALSE removes. Default is FALSE.
#' @param fusion An optional MAF file for fusions. If inputed the outcome will be added to the matrix with columns ending in ".fus".
#' Default is NULL.
#' @param cna An optional CNA files. If inputed the outcome will be added to the matrix with columns ending in ".del" and ".amp".
#' Default is NULL. Note that this file must have patients as columns and genes as rows. binmat expects a matrix with
#' values between -2 and 2. Please do not use any other format. Other functions in the package are available to deal with more detailed
#' CNA data.
#' @param cna.relax for cna data only enables to count both gains and shallow deletions as amplifications and deletions respectively.
#' @param spe.plat boolean specifying if specific IMPACT platforms should be considered. When TRUE NAs will fill the cells for genes
#' of patients that were not sequenced on that plaform. Default is TRUE.
#' @return mut : a binary matrix of mutation data
#' @export
#' @examples library(gnomeR)
#' mut.only <- binmat(maf = mut)
#' all.platforms <- binmat(patients = unique(mut$Tumor_Sample_Barcode)[1:100],maf = mut,fusion = fusion,cna = cna)
#' @import dplyr
#' @import stringr


###############################################
###### MAIN FUNCTION GROUPING EVERYTHING ######
###############################################

binmat <- function(patients=NULL, maf = NULL, mut.type = "SOMATIC",SNP.only = F,include.silent = F,
                   fusion = NULL,cna = NULL,cna.relax = F, spe.plat = F){

  if(is.null(maf) && is.null(fusion) && is.null(cna)) stop("You must provided one of the three following files: MAF, fusion or CNA.")

  mut <- NULL

  if(!is.null(maf)){
    # quick data checks #
    if(length(match("Tumor_Sample_Barcode",colnames(maf))) == 0)
      stop("The MAF file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
    if(length(match("Hugo_Symbol",colnames(maf))) == 0)
      stop("The MAF file inputted is missing a gene name column. (Hugo_Symbol)")
    if(length(match("Variant_Classification",colnames(maf))) == 0)
      stop("The MAF file inputted is missing a variant classification column. (Variant_Classification)")
    if(length(match("Mutation_Status",colnames(maf))) == 0)
      warning("The MAF file inputted is missing a mutation status column (Mutation_Status). It will be assumed that
            all variants are of the same type (SOMATIC/GERMLINE).")

    # set maf to maf class #
    maf <- structure(maf,class = c("data.frame","maf"))
    # getting mutation binary matrix #
    mut <- createbin(obj = maf, patients = patients, mut.type = mut.type,cna.relax = cna.relax,
                     SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
  }

  # fusions #
  if(!is.null(fusion)){
    fusion <- as.data.frame(fusion)
    fusion <- structure(fusion,class = c("data.frame","fusion"))
    fusion <- createbin(obj = fusion, patients = patients, mut.type = mut.type,
                        SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
    if(!is.null(mut)){
      mut <- as.data.frame(cbind(mut,fusion))
      rownames(mut) <- patients}
    else mut <- fusion
  }

  # cna #
  if(!is.null(cna)){
    cna <- as.data.frame(cna)
    cna <- structure(cna,class = c("data.frame","cna"))
    cna <- createbin(obj = cna, patients = patients, mut.type = mut.type,cna.relax = cna.relax,
                     SNP.only = SNP.only, include.silent = include.silent, spe.plat = spe.plat)
    if(!is.null(mut)){
      mut <- as.data.frame(cbind(mut,cna))
      rownames(mut) <- patients}
    else mut <- cna
  }

  # specific platform for IMPACT #
  if(spe.plat){
    v=strsplit(patients, "-IM")
    if(!all(lapply(v, length) == 2))
      warning("All patients were not sequenced on the IMPACT platform or some were mispecified. '-IM' requiered in sample ID.
              The spe.plat argument has been overwritten to FALSE.")
    else{
      v=unlist(lapply(1:length(v), function(x)v[[x]][2]))
      if(length(unique(v)) == 1)
        warning("All patients were not sequenced on the IMPACT platform or some were mispecified. '-IM' requiered in sample ID.
              The spe.plat argument has been overwritten to FALSE.")
      else{
        # remove 410 platform patients #
        missing <- setdiff(c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp")),
                           c(g.impact$g410, paste0(g.impact$g410,".fus"),paste0(g.impact$g410,".Del"),paste0(g.impact$g410,".Amp")))
        if(sum(v == "5") > 0 && sum(missing %in% colnames(mut)) > 0)
          mut[which(v == "5"), na.omit(match(missing, colnames(mut)))] <- NA

        # remove 341 platform patients #
        missing <- setdiff(c(g.impact$g468, paste0(g.impact$g468,".fus"),paste0(g.impact$g468,".Del"),paste0(g.impact$g468,".Amp")),
                           c(g.impact$g341, paste0(g.impact$g341,".fus"),paste0(g.impact$g341,".Del"),paste0(g.impact$g341,".Amp")))
        if(sum(v == "3") > 0 && sum(missing %in% colnames(mut)) > 0)
          mut[which(v == "3"), na.omit(match(missing, colnames(mut)))] <- NA
      }
    }
  }

  return(mut)
}


##############################################
###### CREATE BINARIES FOR DIFF CLASSES ######
##############################################


createbin <- function(obj, patients, mut.type, SNP.only,include.silent, cna.relax, spe.plat){
  UseMethod("createbin")
}

createbin.default <- function(obj) {
  cat("The data did not match any known data type. Please review it and make sure it is correctly specified.")
}


##############################################
############# MUTATION MATRIX ################
##############################################

createbin.maf <- function(obj, patients, mut.type, SNP.only, include.silent, cna.relax, spe.plat){
  maf <- obj
  maf$Hugo_Symbol <- as.character(maf$Hugo_Symbol)
  # recode gene names that have been changed between panel versions to make sure they are consistent and counted as the same gene
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


  # filter/define patients #
  if(!is.null(patients)) maf <- maf[maf$Tumor_Sample_Barcode %in% patients,]
  else patients <- as.character(unique(maf$Tumor_Sample_Barcode))

  # # clean gen dat #
  if(SNP.only) SNP.filt = "SNP"
  else SNP.filt = unique(maf$Variant_Type)

  if(!include.silent) Variant.filt = "Silent"
  else Variant.filt = ""

  if(tolower(mut.type) == "all") Mut.filt = unique(maf$Mutation_Status)
  else Mut.filt = mut.type

  maf <- maf %>% filter(Variant_Classification != Variant.filt,
                        Variant_Type %in% SNP.filt,
                        tolower(Mutation_Status) %in% tolower(Mut.filt))


  #### out frame
  mut <- as.data.frame(matrix(0L,nrow=length(patients),ncol=length(unique(maf$Hugo_Symbol))))
  colnames(mut) <- unique(maf$Hugo_Symbol)
  rownames(mut) <- patients

  for(i in patients){
    genes <- maf$Hugo_Symbol[maf$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0){mut[match(i,rownames(mut)),match(unique(as.character(genes)),colnames(mut))] <- 1}
  }

  missing.mut <- apply(mut,1,function(x){sum(x)==0})
  if(sum(missing.mut) > 0)
    warning(paste0("Some patients did not have any mutations found in the MAF file.", rownames(mut)[missing.mut], collapse = ","))

  return(mut)
}


###########################################
############# FUSION MATRIX ###############
###########################################

createbin.fusion <- function(obj, patients, mut.type, SNP.only,include.silent, cna.relax, spe.plat){
  fusion <- obj
  # quick data checks #
  if(length(match("Tumor_Sample_Barcode",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a patient name column. (Tumor_Sample_Barcode)")
  if(length(match("Hugo_Symbol",colnames(fusion))) == 0)
    stop("The fusion file inputted is missing a gene name column. (Hugo_Symbol)")


  fusion <- fusion %>%
    filter(Tumor_Sample_Barcode %in% patients)

  #### out frame
  fusion.out <- as.data.frame(matrix(0L,nrow=length(patients),ncol=length(unique(fusion$Hugo_Symbol))))
  colnames(fusion.out) <- unique(fusion$Hugo_Symbol)
  rownames(fusion.out) <- patients

  for(i in patients){
    genes <- fusion$Hugo_Symbol[fusion$Tumor_Sample_Barcode %in% i]
    if(length(genes) != 0){fusion.out[match(i,rownames(fusion.out)),
                                      match(unique(as.character(genes)),colnames(fusion.out))] <- 1}
  }
  colnames(fusion.out) <- paste0(colnames(fusion.out),".fus")
  return(fusion.out)
}


################################################
############# COPY NUMBER MATRIX ###############
################################################

createbin.cna <- function(obj, patients, mut.type, SNP.only,include.silent, cna.relax, spe.plat){
  cna <- obj
  rownames(cna) <- cna[,1]
  cna <- cna[,-1]
  cna <- as.data.frame(t(cna))
  rownames(cna) <- gsub("\\.","-",rownames(cna))
  cna <- cna[rownames(cna) %in% patients,]

  temp <- do.call("cbind",apply(cna,2,function(x){
    if(cna.relax){
      yA <- ifelse(x>=0.9,1,0)
      yD <- ifelse(x<=-0.9,1,0)
    }
    if(!cna.relax){
      yA <- ifelse(x==2,1,0)
      yD <- ifelse(x==-2,1,0)
    }
    out <- as.data.frame(cbind(yA,yD))
    colnames(out) <- c("Amp","Del")
    return(out)
  }))

  cna <- temp[,apply(temp,2,function(x){sum(x,na.rm=T) > 0})]

  # add missing
  missing <- patients[which(is.na(match(patients,rownames(cna))))]
  add <- as.data.frame(matrix(0L,nrow = length(missing),ncol = ncol(cna)))
  rownames(add )  <- missing
  colnames(add )<- colnames(cna)
  cna <- as.data.frame(rbind(cna,add))
  cna <- cna[match(patients,rownames(cna)),]

  return(cna)
}
