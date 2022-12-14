% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_bgsig.R
\name{RemoveBackground_vector_single}
\alias{RemoveBackground_vector_single}
\title{Return single vector}
\usage{
RemoveBackground_vector_single(
  background_profile,
  sig_profile,
  sampling_number,
  start_num,
  boundary = 2
)
}
\arguments{
\item{background_profile}{background mutational catalog}

\item{sig_profile}{experimentally generated mutational catalog}

\item{sampling_number}{number of bootstrapping samples}

\item{start_num}{mutation burden to start}

\item{boundary}{Range of signature around centroid (default = 2)}
}
\value{
data.frame including background signature and experiment signature
}
\description{
Return single vector
}
