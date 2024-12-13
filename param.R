#' Parameter file for MS data processing
#'
#' All MALDIquant and other parameters are set here
#' 
#' wdir: directory where spectra are located
#' hws: MALDIquant halfWindowSize parameter
#' smooth.method: MALDIquant smoothIntensity method
#' br.method: MALDIquant removeBaseline method
#' pdetect.method: MALDIquant detectPeaks method
#' snr: MALDIquant detectPeaks SNR parameter
#' bin.method: MALDIquant binPeaks method
#' tol: MALDIquant binPeaks tolerance parameter
#' thresh: Cosine similarity threshold to remove outlier technical replicates
#' m1,m2: mass range of data
#' rpeaks: reference peaks for internal recalibration
#' impute: logical,impute missing values and normalize to TIC of each spectrum
#' fvalue: missing values are imputed with (min intensity)/fvalue

# Check to see if required packages are installed and install them if not found
packages = c("tidyverse","DescTools","MALDIquant","MALDIquantForeign","reshape2","grDevices","matrixStats","coop")
check.pack = lapply(packages,function(x) if(!require(x,character.only=TRUE)){install.packages(x,dependancies=TRUE)
  library(x,character.only=TRUE)})

# Parameters for MALDIquant
wdir = "/Users/nandhini_sokkalingam/edibleoil/"
file=file.path(wdir,"data")
hws = 5 
smooth.method="SavitzkyGolay"
br.method="SNIP"
pdetect.method="MAD" 
snr = 5 
bin.method="strict"
tol=700e-06 
thresh = 0.9 
m1 = 100 
m2 = 1000 

# set internal calibration TRUE/FALSE
calibrate = FALSE

# set peaks (m/z) for internal recalibration if calibrate = TRUE
# Set at least 4-5 peaks for better calibration
#rpeaks = c(100,250,500,750,1000) # Example

# set impute and normalize to TRUE/FALSE
impute = TRUE
fvalue = 5


