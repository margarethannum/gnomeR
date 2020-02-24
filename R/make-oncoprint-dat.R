#' dat.oncoPrint
#'
#' Enables creation of a matrix used to generate an OncoPrint heatmap.
#'
#' @param gen.dat A binary matrix or dataframe, with patients as rows and features as columns. Note that the names of the
#' columns must end in ".Del" or ".Amp" to recognize copy number alterations. (see create.bin.matrix for more details on this format).
#' @param clin.dat An optional clinical file, including only the features the user wishes to add to the plot. Default is NULL.
#' @return mat : a matrix ready to be plotted using plot.Oncoprint().
#' @export
#'
#' @import
#' tibble


dat.oncoPrint <- function(gen.dat,clin.dat=NULL){

  # would be best if genetics also had an unknown #
  if(anyNA(gen.dat))
    gen.dat[is.na(gen.dat)] <- 0

  if(!is.null(clin.dat)){
    # seet NA's to UNKNOWN #
    clin.dat <- clin.dat %>%
      tibble::rownames_to_column('sample') %>%
      mutate_all(as.character) %>%
      tibble::column_to_rownames('sample')
    if(anyNA(clin.dat)) clin.dat[is.na(clin.dat)] <- "Unknown"
    # subset data #
    patients <- intersect(rownames(gen.dat),rownames(clin.dat))
    gen.dat <- gen.dat[match(patients,rownames(gen.dat)),]
    clin.dat <- clin.dat %>%
      tibble::rownames_to_column('sample')
    clin.dat <- clin.dat[match(patients,clin.dat$sample),]
    rownames(clin.dat) <- clin.dat$sample
    clin.dat <- clin.dat %>%
      select(-one_of("sample"))
    ## names ##
    genes <- c(colnames(clin.dat),unique(gsub(".Amp|.Del|.fus","",colnames(gen.dat))))
  }
  else genes <- unique(gsub(".Amp|.Del|.fus","",colnames(gen.dat)))

  mat <- as.data.frame(matrix(nrow=length(genes),ncol=nrow(gen.dat)))
  colnames(mat) <- rownames(gen.dat)
  rownames(mat) <- genes

  # mutations #
  mut <- gen.dat[,na.omit(match(genes,colnames(gen.dat)))]

  for(i in 1:nrow(mut)){
    for(j in 1:ncol(mut)){
      if(!is.na(mut[i,j])){
        if(mut[i,j]==1) {
          mat[match(colnames(mut)[j],rownames(mat)),match(rownames(mut)[i],colnames(mat))] <- "MUT;"
        }
      }
    }
  }

  # del #
  if(length(na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))) > 0){
    del <- gen.dat[,na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))]
    if(is.null(dim(del))) {
      del <- as.data.frame(del)
      rownames(del) <- rownames(gen.dat)
    }
    colnames(del) <- gsub(".Del","",colnames(gen.dat)[na.omit(match(paste0(genes,".Del"),colnames(gen.dat)))])
    for(i in 1:nrow(del)){
      for(j in 1:ncol(del)){
        if(del[i,j]==1) {
          if(!is.na(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))])){
            mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <-
              paste0(mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))],"DEL;")}
          else{mat[match(colnames(del)[j],rownames(mat)),match(rownames(del)[i],colnames(mat))] <- "DEL;"}
        }
      }
    }
  }

  # amp #
  if(length(na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))) > 0){
    amp <- gen.dat[,na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))]
    if(is.null(dim(amp))) {
      amp <- as.data.frame(amp)
      rownames(amp) <- rownames(gen.dat)
    }
    colnames(amp) <- gsub(".Amp","",colnames(gen.dat)[na.omit(match(paste0(genes,".Amp"),colnames(gen.dat)))])
    for(i in 1:nrow(amp)){
      for(j in 1:ncol(amp)){
        if(amp[i,j]==1) {
          if(!is.na(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))])){
            mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <-
              paste0(mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))],"AMP;")}
          else{mat[match(colnames(amp)[j],rownames(mat)),match(rownames(amp)[i],colnames(mat))] <- "AMP;"}
        }
      }
    }
  }

  # fusions #
  if(length(na.omit(match(paste0(genes,".fus"),colnames(gen.dat)))) > 0){
    fusion <- gen.dat[,na.omit(match(paste0(genes,".fus"),colnames(gen.dat)))]
    if(is.null(dim(fusion))) {
      fusion <- as.data.frame(fusion)
      rownames(fusion) <- rownames(gen.dat)
    }
    colnames(fusion) <- gsub(".fus","",colnames(gen.dat)[na.omit(match(paste0(genes,".fus"),colnames(gen.dat)))])
    for(i in 1:nrow(fusion)){
      for(j in 1:ncol(fusion)){
        if(fusion[i,j]==1) {
          if(!is.na(mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))])){
            mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))] <-
              paste0(mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))],"FUS;")}
          else{mat[match(colnames(fusion)[j],rownames(mat)),match(rownames(fusion)[i],colnames(mat))] <- "FUS;"}
        }
      }
    }
  }
  ### clin ###
  if(!is.null(clin.dat)){
    temp <- clin.dat
    for(x in colnames(clin.dat)){
      y <- clin.dat[,x]

      # if binary #
      if(length(unique(y[!is.na(y)])) == 2){
        mat[match(x,rownames(mat)),] <- paste0(x,"_",y,";") #ifelse(y==1,"CLIN;",NA)
      }

      else if(length(unique(y[!is.na(y)])) %in% c(3:5)){
        mat[match(x,rownames(mat)),] <- paste0(x,"_",y,";") #abbreviate(x,1)
      }

      else{
        mat[match(x,rownames(mat)),] <- ifelse(y > median(y),"CLIN;",NA)
      }
    }
  }
  ###
  mat[is.na(mat)] <- "  "
  mat <- as.matrix(mat)

  return(mat)
}
