% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct_natural_abundance.R
\name{flx.correct_natural_abundance}
\alias{flx.correct_natural_abundance}
\title{Correction of isotopic natural abundance.}
\usage{
flx.correct_natural_abundance(
  dataset,
  resolution = 120000,
  resolution_defined_at = 200
)
}
\arguments{
\item{dataset}{The input may be a dataset after background subtraction, or without subtraction. In the latter case, the blank samples, if any in the dataset input, will be processed just as a normal sample.}

\item{resolution}{Mass spectrometry scan resolution. Default value set to 120000}

\item{resolution_defined_at}{Molecular mass (Dalton) at which the mass spectrometry resolution is defined. Default value set to 200 Dalton.}
}
\value{
A list of four data frames: the original data set, corrected intensity, corrected and normalized intensity, and the pool.
}
\description{
This is a wrapper function of "natural_abundance_correction" in the accucor package. It functions to correct the 13-carbon natural abundance to unveil the net 13-C enrichment due to experimental 13-C tracer labeling.
}
