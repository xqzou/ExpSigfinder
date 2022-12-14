% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_bgsig.R
\name{Wrap_KOSig}
\alias{Wrap_KOSig}
\title{Wrap for subtracting the background to estimate the experimentally generated signature}
\usage{
Wrap_KOSig(
  MutCatalogue,
  bg_column,
  ko_column,
  sampling_number,
  start_num,
  boundary,
  outputname
)
}
\arguments{
\item{MutCatalogue}{96-channel mutational catalog}

\item{bg_column}{the column name of the background catalog}

\item{ko_column}{the column names of the experiment catalog}

\item{sampling_number}{number of bootstrapping samples}

\item{start_num}{mutation burden to start}

\item{boundary}{Range of signature around centroid (default = 2)}

\item{outputname}{output file name}
}
\description{
Wrap for subtracting the background to estimate the experimentally generated signature
}
