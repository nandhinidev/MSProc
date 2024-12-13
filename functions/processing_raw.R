#' MALDIquant processing (peak detection,binning & filtering)
#' 
#' This function processes the smoothed and baselined spectra
#' @param datalist MALDIquant object (List of mass spectra)
#' @param hws MALDIquant halfWindowSize parameter
#' @param pdetect.method MALDIquant detectPeaks method
#' @param snr MALDIquant detectPeaks SNR parameter
#' @param bin.method MALDIquant binPeaks method
#' @param tol MALDIquant binPeaks tolerance parameter
#' @param dir directory to write processed spectra
#' @return Table of features (rows) and sample replicates (columns) extracted from MALDIquant list object
#'
process <- function(datalist,hws,pdetect.method,snr,bin.method,tol,dir,m1,m2){
  fname = unname(unlist(lapply(datalist,function(x) metaData(x)$file)))
  msobj.pl = MALDIquant::detectPeaks(datalist,halfWindowSize=hws,method=pdetect.method,SNR=snr)
  msobj.pl.tr = trim(msobj.pl,range=c(m1,m2))
  
  # Create directory to store
  if(!dir.exists(file.path(getwd(),"output",batch,dir))){dir.create(file.path(wdir,"output",batch,dir))}
  dir.cent = lapply(unique(fname),function(x) dir.create(file.path(wdir,"output",batch,dir,x)))
 
  # Write centroided peaks
  msobj.pl.cent = lapply(msobj.pl.tr,function(x) as.data.frame(cbind(mz=MALDIquant::mass(x),int=MALDIquant::intensity(x))))
  msobj.pl.aligned = MALDIquant::binPeaks(msobj.pl.tr,method=bin.method,tolerance=tol)
  #msobj.pl.aligned.fl = MALDIquant::filterPeaks(msobj.pl.aligned,minFrequency = 0.75,labels=factor(unlist(fname),unique(unlist(fname))))
  int = MALDIquant::intensityMatrix(msobj.pl.aligned)
  datatable = t(int)

  #count replicates
  nreps = lapply(unique(unlist(fname)),function(x) sum(stringr::str_count(unlist(fname),paste0("^",x,"$"))))
  spread = lapply(nreps,function(x) seq(1,x,1))
  reps = lapply(spread,function(x) paste0("Rep",x))
  repsname = paste0(unlist(fname),unlist(reps))
  colnames(datatable) = repsname
  msobj.pl.cent.write = mapply(function(x,y,z) write.table(x,file=file.path(wdir,"output",batch,dir,z,paste0(y,"_cent.txt")),sep="\t",row.names = FALSE,quote = FALSE,col.names = FALSE), x = msobj.pl.cent,y = repsname,z=fname)
  return(datatable)
}

