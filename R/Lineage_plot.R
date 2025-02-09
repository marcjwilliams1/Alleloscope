#' Generate genotype plot (scatter plot) for each region and save in the plot directory.
#'
#' @param Obj_filtered An Alleloscope object with a n cell by (m region * 2) genotype_values matrix and seg_table_filtered matrix.
#' Every 2 columns in the genotype_values matrix are (rho_hat, theta_hat) of each region.
#' @param nSNP An integer for the minimum number of SNPs across segments. Segments with the number of SNPs < nSNP are excluded. 
#' @param clust_method Method for clustering. Please refer to the "pheatmap" function. 
#' @param plot_conf Logical (TRUE/FALSE). Whether or not to plot the confidence scores under the lineage tree.
#' @param nclust An integer for the number of subclones gapped in the plot.
#' @param plot_path The path for saving the plot.
#' @param all_chr Logical (TRUE/FALSE). Whether or not the analysis is at the whole-genome level. 
#' @param maxcp Integer. Setting the maximum number of copies for the analysis.
#'
#' @return Plot the lineage tree constructed using cell-level genotypes across all regions and return a vector indicating the cluster identity (of the nclust clusters) of each cell.
#'
#' @import ggplot2
#' @import pheatmap
#' @import cluster
#' @export

Lineage_plot=function(Obj_filtered=NULL, nSNP=0, clust_method='ward.D2', nclust=5, plot_conf=FALSE,plot_path=NULL, all_chr=FALSE, maxcp=6){
  
  # check parameters
  if(is.null(Obj_filtered)){
    stop("Please provide a valid Alleloscope object for Obj_filtered.")
  }else if(nSNP<0){
    stop("Please provide a positive integer for nSNP.")
  }else if(nclust<0 | nclust > length(Obj_filtered$barcodes)){
    stop("nclust should be in the range [0, ncell].")
  }
  
  # assign values
  samplename=Obj_filtered$samplename
  ref=Obj_filtered$ref
  theta_hat_cbn=Obj_filtered$genotype_values
  segmentation=Obj_filtered$seg_table_filtered
  #region_list=Obj_filtered$seg_table_filtered$chrr
  region_list=sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2)[!duplicated(sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2))]
  if(is.null(plot_path)){
    plot_path=paste0(Obj_filtered$dir_path,'/plots/', "lineage_ref_",ref,'.pdf')[1]}
  #plot_path=paste0(Obj_filtered$dir_path,'/plots/lineage_ref_',ref,'.pdf')}
  #dir.create(plot_path)
  
  if(!dir.exists(paste0(Obj_filtered$dir_path,'/plots'))){
    dir.create(plot_path)
  }
  
  pp_list=c()
  
  cluster_cbn=NULL  # category
  cluster_cbn2=NULL # values for plotting
  
  for(chrr in region_list){
    theta_N_sub=theta_hat_cbn[,which(sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2)==chrr)]
    df=data.frame("rho_hat"=theta_N_sub[,1], "theta_hat"=theta_N_sub[,2])
    
    cluster1=genotype_neighbor(df,cluster_name = TRUE, maxcp = maxcp)
    cluster2=genotype_neighbor(df, cluster_name = FALSE, maxcp = maxcp)
    
    cluster_cbn=cbind(cluster_cbn, cluster1)
    cluster_cbn2=cbind(cluster_cbn2, cluster2)
    
  }
  
  colnames(cluster_cbn)=region_list
  rownames(cluster_cbn)=rownames(theta_hat_cbn)
  cluster_cbn=data.frame(cluster_cbn, stringsAsFactors = T)
  
  colnames(cluster_cbn2)=region_list
  rownames(cluster_cbn2)=rownames(theta_hat_cbn)
  cluster_cbn2=data.frame(cluster_cbn2, stringsAsFactors = F)
  
  if(plot_conf){
    cluster_cbn_conf=Obj_filtered$genotypeConfidence###
    colnames(cluster_cbn_conf)=region_list
    rownames(cluster_cbn_conf)=rownames(theta_hat_cbn)
    cluster_cbn_conf=data.frame(cluster_cbn_conf, stringsAsFactors = F)
  }
  
  # select segments for plotting
  ind=which(region_list %in% segmentation$chrr[which(segmentation$Freq>nSNP)])
  segmentation=segmentation[match(region_list[ind],segmentation$chrr),]
  cluster_cbn=cluster_cbn[,ind]
  cluster_cbn2=cluster_cbn2[,ind]
  if(plot_conf){
    cluster_cbn_conf=cluster_cbn_conf[,ind]
  }
  
  # plot the genotypes on the genome (5000000 bins)
  nrep=round(as.numeric(segmentation[,5])/5000000)
  
  cnrep=c(0,cumsum(nrep))
  #cluster_cbn_all=matrix(ncol=sum(nrep), nrow=nrow(cluster_cbn))
  cluster_cbn2_all=matrix(ncol=sum(nrep), nrow=nrow(cluster_cbn))
  cluster_cbn_conf_all=matrix(ncol=sum(nrep), nrow=nrow(cluster_cbn))
  
  for(ii in 1:ncol(cluster_cbn)){
    if(nrep[ii]>0){
      #rr=(replicate(nrep[ii],cluster_cbn[,ii]))
      #cluster_cbn_all[,(cnrep[ii]+1):(cnrep[ii+1])]=rr
      rr=(replicate(nrep[ii],cluster_cbn2[,ii]))
      cluster_cbn2_all[,(cnrep[ii]+1):(cnrep[ii+1])]=rr
      
      if(plot_conf){
        rr=(replicate(nrep[ii],cluster_cbn_conf[,ii]))
        cluster_cbn_conf_all[,(cnrep[ii]+1):(cnrep[ii+1])]=rr}
    }
  }
  chrgap=c()
  for(ii in 1:22){
    chrgap=c(chrgap,sum(nrep[which(segmentation$chr==ii)]))
  }
  
  mu0=NULL
  for(ii in 1:maxcp){
    mu_tmp=c(rep(0.5*ii,ii+1),c(0,( (1:(ii))/ii)))
    mu0=rbind(mu0, matrix(mu_tmp, byrow=F, ncol=2))
  }
  
  if(nrow(mu0)>=21){
    col=c('#4d9efa','#0323a1','#9ecae1','#b0b0b0','#00d9ff',
          "#fff1ba","#ffb521","#DC7633","#BA4A00",
          "#fde0dd","#fcc5c0","#f768a1","#ae017e","#49006a",
          "#c7e9b4","#7fcdbb","#41b6c4","#41ab5d","#006d2c", "#000000",
          "#7B241C", rep("#7B241C", (nrow(mu0)-21)))
  }else{
    col=c('#4d9efa','#0323a1','#9ecae1','#b0b0b0','#00d9ff',
          "#fff1ba","#ffb521","#DC7633","#BA4A00",
          "#fde0dd","#fcc5c0","#f768a1","#ae017e","#49006a",
          "#c7e9b4","#7fcdbb","#41b6c4","#41ab5d","#006d2c", "#000000",
          "#7B241C")[1:nrow(mu0)]
  }
  
  #cluster_cbn_all=data.frame(cluster_cbn_all, stringsAsFactors = F)
  
  test=gower.dissimilarity.mtrx <- cluster::daisy(cluster_cbn, metric = c("gower"))
  hh=hclust(test, method=clust_method)
  hh$height=log(hh$height+1) # linear transform thee tree length
  
  
  plot_matrix<- cluster_cbn2_all
  
  pdf(paste0(plot_path ), width = 12,height = 6)
  
  if(all_chr==TRUE){
    col_lab=rep(" ", ncol(cluster_cbn2_all))
    col_lab[c(0, cumsum(chrgap)[1:(length(chrgap)-1)])+chrgap/2]=paste0("chr",as.character(1:22))
    
    tmp=pheatmap::pheatmap(plot_matrix,
                           cluster_cols = F, cluster_rows = hh,
                           show_rownames = F,
                           labels_col=col_lab,
                           clustering_distance_rows = "euclidean",
                           clustering_method = clust_method,
                           gaps_col=cumsum(chrgap),
                           cutree_rows = nclust,
                           color =col,
                           breaks = 0:26)
  }else{
    tmp=pheatmap::pheatmap(plot_matrix,
                           cluster_cols = F, cluster_rows = hh,
                           show_rownames = F,
                           clustering_distance_rows = "euclidean",
                           clustering_method = clust_method,
                           gaps_col=cumsum(chrgap),
                           cutree_rows = nclust,
                           color =col,
                           breaks = 0:26)}
  
  dev.off()
  
  message(paste0("Lineage plot is successfully saved in the path:",plot_path))
  
  if(plot_conf){
    #conf_plot(Obj_filtered = Obj_filtered, nclust = nclust)
    plot_matrix<- cluster_cbn_conf_all
    conf_path=paste0(Obj_filtered$dir_path,'/plots/', "conf_ref_",ref,'.pdf')
    
    pdf(conf_path, width = 12,height = 6)
    
    if(all_chr==TRUE){
      col_lab=rep(" ", ncol(cluster_cbn2_all))
      col_lab[c(0, cumsum(chrgap)[1:(length(chrgap)-1)])+chrgap/2]=paste0("chr",as.character(1:22))
      
      pheatmap::pheatmap(plot_matrix,
                         cluster_cols = F, cluster_rows = hh,
                         show_rownames = F,
                         labels_col=col_lab,
                         clustering_distance_rows = "euclidean",
                         clustering_method = clust_method,
                         gaps_col=cumsum(chrgap),
                         cutree_rows = nclust)
      #color =col,
      #breaks = 0:26)
    }else{
      pheatmap::pheatmap(plot_matrix,
                         cluster_cols = F, cluster_rows = hh,
                         show_rownames = F,
                         clustering_distance_rows = "euclidean",
                         clustering_method = clust_method,
                         gaps_col=cumsum(chrgap),
                         cutree_rows = nclust)
      #color =col,
      #breaks = 0:26)
    }
    
    
    
    dev.off()
    
    message(paste0("Confidence plot is successfully saved in the path:",conf_path))
    
  }
  
  clust=cutree(tmp$tree_row, k=nclust)
  clust=clust[hh$order]
  
  return(clust)
  
}
