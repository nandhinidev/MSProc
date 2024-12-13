#' Missing value imputation and normalization
#'
#' This function fills missing values with a constant and normalizes intensities to TIC
#' @param dtable Feature table after processing
#' @param fvalue missing values are imputed with (min intensity)/fvalue 
#' @return Missing value imputed and normalized feature table
#'
imputeproc <- function(dtable,fvalue){
  # Remove zero variance features (single intensity value for a feature)
  filt = dtable[rowSums(!is.na(dtable))>1,]
  min.int = min(filt[filt>0],na.rm=TRUE)/fvalue
  filt[is.na(filt)] = min.int
  # Total intensity normalization (final values from 0 to 1000)
  data.norm = as.data.frame(filt) %>% purrr::map(function(x) 1000*x/sum(x,na.rm=TRUE)) %>% as.data.frame()
  final = cbind("mz"=row.names(filt),data.norm)
  return(cbind.data.frame("mz"=row.names(filt),filt))
  }
