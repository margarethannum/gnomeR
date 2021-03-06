#' facets.dat
#'
#' Creates a copy number alteration matrix from segment files.
#'
#' @param seg a segmentation file containing the segmentation information of multiple patients
#' @param filenames the names of the segment files to be loaded and processed (Note must end in ".Rdata").
#' @param path the relative path to the files folder from your current directory
#' @param patients the names of the patients of the respective filenames. Default simply 1 to number of files.
#' @param min.purity the minimum purity of the sample required to be kept in the final dataset. Default is 0.3.
#' @param epsilon level of unions when aggregating segments between
#' @param adaptive
#' @return out.cn : a matrix of the union of all segment files with patients as rows and segments as columns
#' @return ploidy : a vector of the ploidy values for the patients in out.cn
#' @return purity : a vector of the purity values for the patients in out.cn
#' @return FGAs : a vector of the fragment of genome altered values for the patients in out.cn
#' @export
#'
#' @examples library(gnomeR)
#'
#' @import
#' iClusterPlus
#' dplyr
#' cluster


facets.dat <- function(seg = NULL,filenames = NULL, path=NULL,
                       patients=NULL, min.purity = 0.3, epsilon = 0.005,
                       adaptive = F){

  if(is.null(seg) && is.null(filenames))
    stop("You must provide either a complete segmentation file
         or a list of files to be loaded with their corresponding path")

  if(!is.null(seg) && !is.null(filenames))
    stop("Please provide either a complete segmentation file or a
         list of segmentation files to be loaded")

  ######################################


  if(!is.null(filenames)){
    if (!file.exists(path)) stop("The path provided cannot be found")
    if(!is.null(patients))
      if(length(patients) != length(filenames)) stop("Length of patients differs from length of filenames")
    if(is.null(patients)) patients <- as.character(1:length(filenames))
    if(min.purity < 0 || min.purity > 1) stop("Please select a purity between 0 and 1, included.")

    ### segment files ###
    for(i in 1:length(filenames)){
      ## set up ##
      all.dat <- data.frame()
      FGAs <- c()
      dipLogR <- c()
      ploidy <- c()
      purity <- c()
      missing <- c()
      s <- 0
      ##
      fit <- NULL
      try(load(paste0(path,"/",filenames[i])),silent=T)

      if(!is.null(fit) && !is.na(fit$purity) && fit$purity >= min.purity){
        ## keep track ##
        s = s + 1
        dipLogR[s] <- fit$dipLogR
        ploidy[s] <- fit$ploidy
        purity[s] <- fit$purity
        #################

        cncf <- as.data.frame(fit$cncf %>% select(chrom,start,end,tcn.em,lcn.em,num.mark))
        cncf.comp <- cncf[complete.cases(cncf),]
        FGAs[s] <- as.numeric(unique(cncf.comp %>%
                                       mutate(
                                         numerator = end - start,
                                         seg.sum = sum(end-start),
                                         FGA = sum(numerator[!(tcn.em==2 & lcn.em==1)])/seg.sum
                                       ) %>%
                                       select(FGA)))

        cncf$sample <- rep(patients[i],nrow(cncf))
        cncf$seg.mean <- log2(cncf$tcn.em/fit$ploidy+ 1*10^(-6))
        cncf <- cncf[,c("sample", "chrom", "start", "end",
                        "num.mark", "seg.mean")]
        all.dat <- rbind(all.dat,cncf)
      }
      else{
        missing <- c(missing,filenames[i])
      }

    }


    out.cn <- CNregions.mod(seg = all.dat,epsilon = epsilon,adaptive = adaptive)
    names(ploidy) <- rownames(out.cn)
    names(purity) <- rownames(out.cn)

    if(!is.null(missing)){
      warning("Some files in the list were not found. You can see a full list in the 'missing' output.")
      return(list("out.cn"=out.cn,"ploidy"=ploidy,"purity"=purity,"FGA"=FGAs,"missing"=missing))
    }
    else{
      return(list("out.cn"=out.cn,"ploidy"=ploidy,"purity"=purity,"FGA"=FGAs))
    }
  }

  ####################################################

  if(!is.null(seg)){
    if(is.null(patients)) patients <- as.character(unique(seg[,1]))
    # else patients <- intersect(patients,unique(seg$ID))
    seg.filt <- seg %>%
      filter(ID %in% patients)
    all.dat <- data.frame()
    ### segment files ###
    for(i in 1:length(patients)){


      cncf <- as.data.frame(seg.filt %>%
                              filter(ID == patients[i]) %>%
                              rename(sample = ID,start = loc.start, end = loc.end) %>%
                              mutate(chrom = as.numeric(as.character(chrom)),
                                     start = as.numeric(as.character(start)),
                                     end = as.numeric(as.character(end)),
                                     num.mark = as.numeric(as.character(num.mark)),
                                     seg.mean = as.numeric(as.character(seg.mean)))
      )
      all.dat <- rbind(all.dat,cncf)
    }
    all.dat <- all.dat[complete.cases(all.dat),]
    out.cn <- CNregions.mod(seg = all.dat,epsilon = epsilon,adaptive = adaptive)
    # out.cn <- out.cn[match(rownames(out.cn),patients),]
    return(list("out.cn"=out.cn,"patients" = patients))
  }

}

