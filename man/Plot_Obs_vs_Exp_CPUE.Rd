% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SS_BDM.R
\name{Plot_Obs_vs_Exp_CPUE}
\alias{Plot_Obs_vs_Exp_CPUE}
\title{Plot observed vs expected CPUE (asymptotic errors)}
\usage{
Plot_Obs_vs_Exp_CPUE(DatFromCSVFile, model_outputs, xaxis_lab = "Year")
}
\arguments{
\item{DatFromCSVFile}{model data read in from csv file}

\item{model_outputs}{object containing outputs from fitted TMB model}

\item{xaxis_lab}{x axis label, to overide default}
}
\value{
plots of observed vs expected cpue
}
\description{
This function plots observed vs expected CPUE and associated 95 percent confidence limits
based on asymptotic errors
}
