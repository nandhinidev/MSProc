source("param.R")
source("functions/initialize_raw.R")
source("functions/processing_raw.R")
source("functions/outliers.R")
source("functions/replicate_correlation.R")
source("functions/average_spectra.R")
source("functions/postprocess_raw.R")
source("functions/imputeproc.R")

# Create output directory
batch = paste0(pdetect.method,"snr",snr,"hws",hws,"tol",tol,"cal",calibrate,"br",remove,"imp",impute,"thresh",thresh)

# Create directories for writing output
if(!dir.exists(file.path(wdir,"output"))){dir.create(file.path(wdir,"output"))}
if(!dir.exists(file.path(wdir,"output",batch))){dir.create(file.path(wdir,"output",batch))}

# Initialize processed data
mass.spectra = suppressWarnings(data.preprocess.raw(file,rpeaks,hws,smooth.method,br.method,pdetect.method,snr,tol,calibrate))

# Generate feature table from processed data for outlier detection
feature.table = process(mass.spectra,hws,pdetect.method,snr,bin.method,tol,dir="data_cent",m1,m2)
write.csv(feature.table,file=file.path(wdir,"output",batch,"preprocessed_data.csv"))

# Outlier replicate detection
out = outliers(feature.table,thresh,fvalue)

# Remove outlier replicates, average remaining spectra per group and generate feature table for statistical analysis
if(!dir.exists(file.path(wdir,"output",batch,"data_avg"))){dir.create(file.path(wdir,"output",batch,"data_avg"))}
spectra.avg = average_spectra(mass.spectra,out)

features.table.gq = postprocess(spectra.avg,hws,pdetect.method,snr,bin.method,tol,dir="data_avg_cent",m1,m2)
write.csv(features.table.gq,file=file.path(wdir,"output",batch,"preprocessed_data_avg.csv"))

# Missing value imputation and normalization
if(impute){
  features.table.gq.imp = imputeproc(dtable = features.table.gq,fvalue=fvalue)
  write.csv(features.table.gq.imp,file=file.path(wdir,"output",batch,"preprocessed_data_avg_imp.csv"),row.names=FALSE)
  }



