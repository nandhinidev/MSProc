#' MALDIquant pre-processing (smoothing, baseline removal & alignment)
#'
#' This function takes individual spectrum files (mzml) and pre-process them
#'
#' Individual mass spectrum of replicates within a group should be inside a folder for that group.
#' Top directory has to be "data"
#' Folder names are parsed for group ids, so name appropriately. See below for required dir structure
#' -data
#' --Group1
#' ---Group1MSFile.mzml
#' --Group2
#' ---Group2MSFile.mzml
#'
#' @param file path to input files
#' @param rpeaks reference peaks in the spectrum to be used for internal calibration if calibrate=TRUE
#' @param hws MALDIquant halfWindowSize parameter
#' @param smooth.method MALDIquant smoothIntensity method
#' @param br.method MALDIquant removeBaseline method
#' @param pdetect.method MALDIquant detectPeaks method
#' @param snr MALDIquant detectPeaks SNR parameter
#' @param tol MALDIquant alignSpectra tolerance parameter
#' @param calibrate logical,internally calibrate spectrum with reference peaks
#' @return MALDIquant list object with baselined,smoothed and aligned spectra from input files
#' 
data.preprocess.raw <- function(file,rpeaks,hws,smooth.method,br.method,pdetect.method,snr,tol,calibrate){
  # Read data
  files = list.files(file,pattern="mz*",recursive=TRUE)
  fname = list()
  msobj = MALDIquantForeign::importMzMl(file,verbose=FALSE)
  for (i in 1:length(files))
  {
    fname[[i]] = unlist(strsplit(files[i],"/"))[1]
    metaData(msobj[[i]])$file = fname[[i]]
  }
  msobj.sm = MALDIquant::smoothIntensity(msobj,method=smooth.method,halfWindowSize=hws)
  msobj.sm.br = MALDIquant::removeBaseline(msobj.sm,method=br.method,iterations=10)
  msobj.pl = MALDIquant::detectPeaks(msobj.sm.br,halfWindowSize=hws,method=pdetect.method,SNR=snr)
  if(calibrate){
    refpeaks = MALDIquant::createMassPeaks(mass=rpeaks,intensity=rep(1,length(rpeaks)))
    pdf(file.path(file,"../","Sample_Calibration.pdf"))
    w= MALDIquant::determineWarpingFunctions(msobj.pl,reference=refpeaks,tolerance=2000e-06,method="linear",plot=TRUE)
    dev.off()
    msobj.al = MALDIquant::alignSpectra(msobj.sm.br,halfWindowSize=hws,noiseMethod=pdetect.method,SNR=snr,tolerance=2000e-06,reference=refpeaks,warpingMethod = "linear")
    return(msobj.al)
  } else{
    return(msobj.sm.br)
  }
}

