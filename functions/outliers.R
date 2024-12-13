#' Outlier replicate detection 
#'
#' This function detects outliers based on calculated Cosine similarity values (missing values imputed and normalized before similarity is calculated)
#' @param dataset Feature table after processing with MALDIquant
#' @param thresh Cosine similarity threshold to remove outlier technical replicates 
#' @return vector of replicate outliers
#'
outliers <- function(dataset,thresh,fvalue){
  feature.table = as.data.frame(dataset)
  mv.rep = min(feature.table[feature.table>0],na.rm=TRUE)/fvalue
  feature.table[is.na(feature.table)] = mv.rep
  data.norm = as.data.frame(feature.table) %>% purrr::map(function(x) x/sum(x,na.rm=TRUE)) %>% as.data.frame(optional=TRUE)
  colnames(data.norm) = colnames(feature.table)
  groups = unlist(lapply(strsplit(colnames(feature.table),"Rep"),function(x) x[1]))
  repsname = colnames(feature.table)

# Cosine similarity
  ccvalues = compute_cc(data.norm)
  out.id = lapply(ccvalues,function(x) which(as.data.frame(x)$cos.sim<thresh))
  nout.id = lapply(ccvalues,function(x) which(as.data.frame(x)$cos.sim>thresh))
  outliers1 = lapply(out.id,function(x) if(length(x)>=1){paste0("Rep",x)})
  # samples with < 3 non-outliers
  outliers2 = lapply(nout.id,function(x) if(length(x)<3 & length(x)>0){paste0("Rep",x)})
  rout1 = mapply(function(x,y) if(length(y)>0){paste0(x,y)},x=names(outliers1),y=outliers1)
  rout2 = mapply(function(x,y) if(length(y)>0){paste0(x,y)},x=names(outliers2),y=outliers2)
  rout = unique(c(unlist(rout1),unlist(rout2)))
  write.table(rout,file=file.path(wdir,"output",batch,"outliers.txt"),sep="\n",quote=FALSE,row.names=FALSE,col.names=FALSE)
  coeff = lapply(ccvalues,function(x) as.data.frame(x)$cos.sim)
  cplot = ggplot(reshape2::melt(coeff),aes(x=L1,y=value)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90),axis.text.y = element_text(size=8)) + ylim(0.2,1) + labs(x="Group",y="Cosine Similarity") + coord_flip() + ggtitle("All Replicates")
  ggsave(plot = cplot,device = "png",width = 5,height=8,units="in",file = file.path(wdir,"output",batch,"Cos-Sim_all-replicates.png"))
  if(!is_empty(rout)){
  data.norm.gq = data.norm[,!(names(data.norm) %in% rout)]
  ccvalues.gq = compute_cc(data.norm.gq)
  coeff.gq = lapply(ccvalues.gq,function(x) as.data.frame(x)$cos.sim)
  if(!is_empty(coeff.gq)){
  cplot.gq = ggplot(reshape2::melt(coeff.gq),aes(x=L1,y=value)) + geom_boxplot() + theme(axis.text.x=element_text(angle=90),axis.text.y = element_text(size=8)) + ylim(0.2,1) + labs(x="Group",y="Cosine Similarity") + coord_flip() + ggtitle("Outliers Removed")
  ggsave(plot = cplot.gq,device = "png",width = 5,height=8,units="in",file = file.path(wdir,"output",batch,"Cos-Sim_outliers-removed.png"))}
  else{
    stop(paste("No replicates passed the outlier test"))
  }
}
return(rout)
}
