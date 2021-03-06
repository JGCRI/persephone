% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cross_validate_LOO.R
\name{cross_validate_LOO}
\alias{cross_validate_LOO}
\title{A function for leave one out cross validation on a predefined partition of data}
\usage{
cross_validate_LOO(ensSiteData, testnmes, ctwtests, individSiteDataPath,
  functionalform, sdfunctionalform, ensemblename, testpartition)
}
\arguments{
\item{ensSiteData}{A list of sites for ensemble under consideration}

\item{testnmes}{The names of the CTW tests}

\item{ctwtests}{The actual ctw tests in full}

\item{functionalform}{The name functional form for the mean of the likelihood}

\item{sdfunctionalform}{The name of the functional form for the standard deviation of the likelihood}

\item{ensemblename}{The name of the ensemble under consideration}

\item{testpartition}{The names of the CTW tests that will be excluded from training to be used for testing.}
}
\description{
A function for leave one out cross validation on a predefined partition of data
}
