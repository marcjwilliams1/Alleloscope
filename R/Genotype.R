#' Genotype each cell for each region and plot the genotypes.
#'
#' @param Obj_filtered An Alleloscope object with a n cell by (m region * 4) genotype_values matrix.
#' Every 2 columns in the genotype_table matrix are (rho_hat, theta_hat) of each region.
#' @param xmax An integer for the x-axis maximum limit.
#' @param plot_path The path for saving the plot.
#' @param rds_path The path for saving the rds object. 
#' @param ref_gt A reference "genotypes" (from scDNA-seq) to help with genotype estimation.
#' @param cell_type A matrix with two columns: COL1- cell barcodes; COL2- cell types ("tumor" and others)
#' @param maxcp Integer. Setting the maximum number of copies for the analysis.
#' @param legend Logical (TRUE/FALSE) Whether or not to show the figure legends.
#'
#' @return A list of ggplot objects of the genotyping results for all the regions.
#'
#' @import ggplot2
#' @import cowplot
#' @export
Genotype=function(Obj_filtered=NULL, xmax=NULL, plot_path=NULL,rds_path=NULL, ref_gt=NULL, cell_type=NULL,maxcp=6, legend=FALSE){
  # check parameters
  if(is.null(Obj_filtered)){
    stop("Please provide a valid Alleloscope object for Obj_filtered.")
  }
  
  
  samplename=Obj_filtered$samplename
  assay=Obj_filtered$assay
  ref=Obj_filtered$ref
  
  theta_hat_cbn=Obj_filtered$genotype_values
  #region_list=sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2)[(1:(ncol(theta_hat_cbn)/4))*4]
  region_list=sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2)[!duplicated(sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2))]
  
  if(is.null(plot_path)){
    plot_path=paste0(Obj_filtered$dir_path,'/plots/gtype_scatter_ref_',ref,'.pdf')[1]}
  
  if(is.null(rds_path)){
    rds_path=paste0(Obj_filtered$dir_path,'/rds/genotypes.rds')}
  
  
  if(!dir.exists(paste0(Obj_filtered$dir_path,'/plots'))){
    dir.create(plot_path)
  }
  
  if(!dir.exists(paste0(Obj_filtered$dir_path,"/rds"))){
    dir.create(paste0(Obj_filtered$dir_path,"/rds"))
  }
  
  genotype_table=matrix(nrow=nrow(theta_hat_cbn), ncol=length(region_list))
  genotypeProb=list()
  genotypeConfidence=matrix(nrow=nrow(theta_hat_cbn), ncol=length(region_list))
  
  pp_list=list()
  
  mu0=NULL
  for(ii in 1:maxcp){
    mu_tmp=c(rep(0.5*ii,ii+1),c(0,( (1:(ii))/ii)))
    mu0=rbind(mu0, matrix(mu_tmp, byrow=F, ncol=2))
  }
  
  rownames(mu0)=paste0("center",1:dim(mu0)[1])
  mu0=as.data.frame(mu0, stringsAsFactors = F)
  
  message("Start genotyping each cell for each region. ")
  ii=0
  for(chrr in region_list){
    ii=ii+1
    theta_N_sub=theta_hat_cbn[,which(sapply(strsplit(colnames(theta_hat_cbn),'_'),'[',2)==chrr)]
    df=data.frame("rho_hat"=theta_N_sub[,1], "theta_hat"=theta_N_sub[,2])
    
    if(is.null(xmax)){
      xm=max(df$rho_hat)
    }else{
      xm=xmax
    }
    
    if(is.null(cell_type)){
      
      if(is.null(ref_gt)){
        cluster=genotype_neighbor(df, maxcp=maxcp)
        genotypeConfidence[,ii]=genotype_conf(X=df, gt=cluster)
        
        
      }else{
        snpCoverage = as.numeric(theta_N_sub[,3])+ as.numeric(theta_N_sub[,4])
        priorProb=table(ref_gt[,chrr])/length(ref_gt[,chrr])                                                
        possibleGenotypes = as.numeric(names(priorProb))                                                   
        est_gt=genotype_ref(X = df, snpCoverage = snpCoverage,
                            possibleGenotypes = possibleGenotypes , priorProb = priorProb )
        cluster=est_gt$genotypes
        genotypeProb[[paste0('chr',chrr)]]=est_gt$genotypeProb
        genotypeConfidence[,ii]=est_gt$genotypeConfidence
        
      }
      
      df$cluster=factor(cluster, levels = 1:nrow(mu0))
      genotype_table[,ii]=cluster
      
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
      
      names(col) <- as.character(1:nrow(mu0))
      
      pp1=ggplot(df, aes(rho_hat, theta_hat)) +
        geom_hline(yintercept = 0.5)+
        geom_vline(xintercept = 1)+
        geom_point(aes(x=V1, y=V2), data=mu0 ,stroke=0.5, size=1, colour="black")+
        geom_point(aes(color = cluster),alpha = 0.5,stroke=0.5, size=0.5)+
        #geom_point(aes(color = cluster, shape=cluster),alpha = 1,stroke=0.5, size=0.8, show.legend = FALSE)+
        ggtitle(paste0(" (chr", as.character(chrr),")")) +
        scale_color_manual(values = col)+
        theme(plot.title = element_text(hjust = 0.5, size = 14))+
        ylab(NULL)+
        xlab(NULL)+
        xlim(0,xm)+ylim(0,1)+
        guides(color=legend)
      
    }else{
      barcodes_tumor=cell_type[which(cell_type[,2]=='tumor'),1]
      barcodes_normal=cell_type[which(cell_type[,2]!='tumor'),1]
      stype=rep("unknown", nrow(df))
      stype[which(rownames(df) %in% barcodes_tumor)]='tumor'
      stype[which(rownames(df) %in% barcodes_normal)]='normal'
      df$type = stype
      
      if(length(unique(df[,1]))!=1){
        pp1=ggplot(df, aes(rho_hat, theta_hat)) + 
          geom_point(alpha = 0.3, size=1,stroke=0, aes(color = type))+
          geom_point(aes(x=V1, y=V2), data=mu0 ,stroke=0.5, size=1, colour="black")+
          stat_density2d(data = df[df$type == "tumor",], fill = "#e41a1c", color = "#e41a1c", geom="polygon", alpha = .1, bins = 5) +
          stat_density2d(data = df[df$type == "normal",], fill = "#1f78b4", color = "#1f78b4", geom="polygon", alpha = .1, bins = 5)+
          scale_color_manual(values = c("tumor" = "#e41a1c","normal" = "#1f78b4", "unknown"='#9F9998')) + 
          geom_hline(yintercept = 0.5)+
          geom_vline(xintercept = 1)+
          ggtitle(paste0(" (chr", as.character(chrr),")")) +
          ylab(NULL) +
          xlab(NULL)+
          xlim(0,xm)+ylim(0,1)+
          guides(color=legend)+
          theme_bw()
      }else{
        pp1=ggplot(df, aes(rho_hat, theta_hat)) + 
          geom_point(alpha = 0.3, size=1,stroke=0, aes(color = type))+
          geom_point(aes(x=V1, y=V2), data=mu0 ,stroke=0.5, size=1, colour="black")+
          scale_color_manual(values = c("tumor" = "#e41a1c","normal" = "#1f78b4", "unknown"='#9F9998')) + 
          geom_hline(yintercept = 0.5)+
          geom_vline(xintercept = 1)+
          ggtitle(paste0(" (chr", as.character(chrr),")")) +
          ylab(NULL) +
          xlab(NULL)+
          xlim(0,xm)+ylim(0,1)+
          guides(color=legend)+
          theme_bw()
      }
      
    }
    
    
    pp_list[[paste0("pp_",as.character(chrr))]]=pp1
    
    cat(paste0(chrr," "))
    
    
  }
  cat("\n")
  
  colnames(genotype_table)=region_list
  rownames(genotype_table)=rownames(theta_hat_cbn)
  
  pdf(paste0(plot_path), width = 6)
  
  pp_list=c(pp_list, vector(mode = "list", length = ceiling(length(pp_list)/6) * 6 - length(pp_list)))
  
  for(ii in c(6*(0:(length(pp_list)/6-1))+1)){
    pp_combine=plot_grid(pp_list[[ii]], pp_list[[ii+1]], pp_list[[ii+2]], pp_list[[ii+3]], pp_list[[ii+4]], pp_list[[ii+5]] ,
                         #labels = region_label[ii:(ii+5)],
                         ncol = 2, nrow = 3)
    
    print(pp_combine)
  }
  
  dev.off()
  
  
  
  
  Obj_filtered$genotypes=genotype_table
  Obj_filtered$genotypeProb=genotypeProb
  Obj_filtered$genotypeConfidence=genotypeConfidence

  
  
  #Obj_filtered$region_plot=pp_list
  
  saveRDS(genotype_table, rds_path)
  
  message(paste0("Genotype plots succefully created and stored in the path:", plot_path))
  cat("\"genotypes\" is added to the Obj_filtered object.\n")
  cat(paste0("Matrix for cell specific genotypes for each region is stored as genotypes.rds in the path", Obj_filtered$dir_path,"\n"))
  
  return(Obj_filtered)
}

