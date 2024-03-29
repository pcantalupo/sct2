% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Self_scmapCluster.R
\name{Self_scmapCluster}
\alias{Self_scmapCluster}
\title{Self_scmapCluster}
\usage{
Self_scmapCluster(sce, clustercols = "RNA_snn_res.0.8")
}
\arguments{
\item{sce}{A SingleCellExperiment object}

\item{clustercols}{Metadata column(s) to use for clustering}
}
\value{
A list of scmapCluster() objects. A self Sankey plot will be generated
for each 'clustercols' and automatically opened in a browser
}
\description{
From a SingleCellExperiment object, this function performs a scmap cluster
analysis on itself. Useful for validating clusters.
}
\examples{

selfscmap = Self_scmapCluster(pbmc_small_sce)

}
\keyword{scmap}
