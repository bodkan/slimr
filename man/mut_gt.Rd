% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf.R
\name{mut_gt}
\alias{mut_gt}
\title{Load mutation genotypes.}
\usage{
mut_gt(vcf, mut_type = NULL, pop_origin = NULL, t_min = -Inf,
  t_max = Inf)
}
\arguments{
\item{vcf}{VCF object.}

\item{mut_type}{Mutation type as an integer.}

\item{pop_origin}{Integer ID of a population of origin.}

\item{t_min}{Lower bound on time of origin.}

\item{t_max}{Upper bound on time of origin.}
}
\value{
GRanges object.
}
\description{
Load mutation genotypes.
}
