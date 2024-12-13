#' Outlier Replicate Removal and replicate spectra averaging
#'
#' This function removes outlier replicates and averages the remaining replicates for each sample
#' @param datalist MALDIquant object (list of mass spectra)
#' @param outliers vector of replicate outliers
#' @return MALDIquant object (List of averaged mass spectrum/label after outlier repliactes removed)
#'
average_spectra <- function(datalist,outliers){
  fname = unname(unlist(lapply(datalist,function(x) metaData(x)$file)))
  nreps = lapply(unique(fname),function(x) sum(stringr::str_count(fname,paste0("^",x,"$"))))
  spread = lapply(nreps,function(x) seq(1,x,1))
  reps = lapply(spread,function(x) paste0("Rep",x))
  repsname = paste0(fname,unlist(reps))
  names(datalist) = repsname
  datalist[names(datalist) %in% unname(outliers)] = NULL
  names.trim = unname(unlist(lapply(datalist,function(x) metaData(x)$file)))
  msobj.avg = suppressWarnings(MALDIquant::averageMassSpectra(datalist,labels=factor(names.trim,unique(names.trim)),method="mean"))
  msobj.avg.write = MALDIquantForeign::exportMzMl(msobj.avg,path=file.path(wdir,"output",batch,"data_avg"))
  return(msobj.avg)
}
