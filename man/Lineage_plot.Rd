% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lineage_plot.R
\name{Lineage_plot}
\alias{Lineage_plot}
\title{Generate genotype plot (scatter plot) for each region and save in the plot directory.}
\usage{
Lineage_plot(
  Obj_filtered = NULL,
  nSNP = 0,
  clust_method = "ward.D2",
  nclust = 5,
  plot_conf = FALSE,
  plot_path = NULL,
  all_chr = FALSE,
  maxcp = 6
)
}
\arguments{
\item{Obj_filtered}{An Alleloscope object with a n cell by (m region * 2) genotype_values matrix and seg_table_filtered matrix.
Every 2 columns in the genotype_values matrix are (rho_hat, theta_hat) of each region.}

\item{nSNP}{An integer for the minimum number of SNPs across segments. Segments with the number of SNPs < nSNP are excluded.}

\item{clust_method}{Method for clustering. Please refer to the "pheatmap" function.}

\item{nclust}{An integer for the number of subclones gapped in the plot.}

\item{plot_conf}{Logical (TRUE/FALSE). Whether or not to plot the confidence scores under the lineage tree.}

\item{plot_path}{The path for saving the plot.}

\item{all_chr}{Logical (TRUE/FALSE). Whether or not the analysis is at the whole-genome level.}

\item{maxcp}{Integer. Setting the maximum number of copies for the analysis.}
}
\value{
Plot the lineage tree constructed using cell-level genotypes across all regions and return a vector indicating the cluster identity (of the nclust clusters) of each cell.
}
\description{
Generate genotype plot (scatter plot) for each region and save in the plot directory.
}
