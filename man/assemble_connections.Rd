% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runCicero.R
\name{assemble_connections}
\alias{assemble_connections}
\title{Combine and reconcile cicero models}
\usage{
assemble_connections(cicero_model_list, silent = FALSE)
}
\arguments{
\item{cicero_model_list}{A list of cicero output objects, generally, the
output of \code{\link{generate_cicero_models}}.}

\item{silent}{Logical, should the function run silently?}
}
\value{
A data frame of connections with their cicero co-accessibility
  scores.
}
\description{
Function which takes the output of \code{\link{generate_cicero_models}} and
assembles the connections into a data frame with cicero co-accessibility
scores.
}
\details{
This function combines glasso models computed on overlapping windows of the
genome. Pairs of sites whose regularized correlation was calculated twice
are first checked for qualitative concordance (both zero, positive or
negative). If they not concordant, NA is returned. If they are concordant
the mean is returned.
}
\examples{
  data("cicero_data")
  data("human.hg19.genome")
  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  sample_genome$V2[1] <- 100000
  input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
  input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                               reduction_method = 'tSNE',
                               norm_method = "none")
  tsne_coords <- t(reducedDimA(input_cds))
  row.names(tsne_coords) <- row.names(pData(input_cds))
  cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)
  model_output <- generate_cicero_models(cicero_cds,
                                         distance_parameter = 0.3,
                                         genomic_coords = sample_genome)
  cicero_cons <- assemble_connections(model_output)

}
\seealso{
\code{\link{generate_cicero_models}}
}
