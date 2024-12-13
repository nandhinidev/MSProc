#' Calculates cosine similarity and correlation coefficient
#'
#' This function generates a correlation table for replicates within each group with four different correlation coefficients
#' Correlation coefficients available to use: Pearson, Spearman, Cosine similarity, Lin's concordance
#' @param dataset Feature table from MALDIquant processing
#' @return Correlation table for replicates within each group
#' @keywords internal
#' 
compute_cc <- function(dataset){
  cosine.sim = corr_coeff(mdata=dataset,"cos.sim")
  lcc = corr_coeff(mdata=dataset,"lcc")
  groups = unique(unlist(lapply(strsplit(colnames(dataset),"Rep"),function(x) x[1])))
  output = lapply(seq_along(groups),function(x) cbind(cosine.sim = cosine.sim[[x]],lcc = lcc[[x]]))
  names(output) = groups
  return(output)
}

compute_cv <- function(dataset){
  feature.table = as.data.frame(dataset)
  groups = unique(unlist(lapply(strsplit(colnames(feature.table),"Rep"),function(x) x[1])))
  data.split = lapply(seq_along(groups),function(x) dplyr::select(feature.table,starts_with(paste0(groups[x],"Rep"))))
  calc = lapply(data.split,function(x) x %>% mutate(av=rowMeans(as.data.frame(.),na.rm=TRUE),sd=matrixStats::rowSds(as.matrix(.),na.rm=TRUE)))
  cv.calc = lapply(calc,function(x) x %>% mutate(cv=100*sd/av))
  cv.intra = lapply(cv.calc,function(x) median(na.omit(as.data.frame(x)$cv)))
  names(cv.intra) = groups
  return(cv.intra)
}

# Functions for Lin's CCC
lc <- function(x,y){
  value = DescTools::CCC(x,y,na.rm=TRUE)
  value$rho.c[,1]
}

# Function for calculating Lin's CCC from data
lc_calc <- function(data){
  data = as.data.frame(data)
  ind = combn(seq(1,ncol(data),1),2)
  cc = lapply(seq(1,ncol(ind),1),function(x) lc(data[,(ind[,x][1])],data[,(ind[,x][2])]))
  new = diag(ncol(data))
  new[lower.tri(new)] = unlist(cc)
  new = t(new)
  new[lower.tri(new)] = unlist(cc)
  new
}

# Corr Coeff (Pearson, Spearman, cosine similarity, LCC)
corr_coeff <- function(mdata,type=c("cos.sim","lcc")){
    mdata.t = as.data.frame(mdata)
    groups = unique(unlist(lapply(strsplit(colnames(mdata.t),"Rep"),function(x) x[1])))
    data.split = lapply(seq_along(groups),function(x) dplyr::select(mdata.t,starts_with(paste0(groups[x],"Rep"))))
  if (!type == "lcc")
  {
    type = match.arg(type)
    corsim = lapply(data.split,function(x) coop::cosine(x,use="complete.obs"))
    avg.coeff = lapply(corsim,function(x) apply(x,1,function(y) mean(y[y!=1])))
    final = lapply(avg.coeff,data.frame)
    out = lapply(final,setNames,type)
    return(out)
  }
  else
  {
    ccc.values = lapply(data.split,lc_calc)
    ccc.avg = lapply(ccc.values,function(x) apply(x,1,function(y) mean(y[y!=1])))
    ccc.final = lapply(ccc.avg,data.frame)
    ccc.out = lapply(ccc.final,setNames,type)
    names(ccc.out) = groups
    return(ccc.out)
  }
}




