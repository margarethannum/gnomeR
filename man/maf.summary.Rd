% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maf-summary.R
\name{maf.summary}
\alias{maf.summary}
\title{maf.summary}
\usage{
maf.summary(maf, mut.type = "SOMATIC")
}
\arguments{
\item{maf}{the names of the segment files to be loaded and processed (Note must end in ".Rdata").}

\item{mut.type}{The mutation type to be used. Options are "SOMATIC", "GERMLINE" or "ALL". Note "ALL" will
keep all mutations regardless of status (not recommended). Default is SOMATIC.}
}
\value{
p.class Barplot of counts of each variant classification

p.type Barplot of counts of each variant type

p.SNV Histogram of counts of each SNV class

p.patient.variant Histogram of counts variants per patient

p.variant.bp Boxplot of the distribution of variant classification per patient

p.variant.dist Boxplot of the relative proportion of each variant class in individual patients

p.variant.dist.bar Stacked barplot of the variants distribution in all patients

p.SNV.dist Boxplot of the relative proportion of each SNV class in individual patients

p.corr Correlation heatmap of the top 10 genes

p.comut Heatmap of the comutation of the top 10 genes
}
\description{
Creates a set of plot summarising a maf file.
}
\examples{
library(gnomeR)
all.plots <- maf.summary(maf=mut)
}
