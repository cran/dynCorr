\name{dynCorrData}
\docType{data}
\alias{dynCorrData}
\title{
An example dataset for use in the example calls in the 
help files for the dynamicCorrelation and bootstrapCI functions
}
\description{
This dataset contains three longitudinal responses for each subject, each measured every ten days over a varying length of follow-up depending upon subject (for each of 34 subjects). Each row contains all three responses for a given subject and a given follow-up time. The features are that the first and third responses have a moderate-sized positive dynamical correlation, while the first and second, and second and third, respectively, have lesser negative dynamical correlations. The dataset is called in the examples on the dynamicCorrelation and bootstrapCI help pages.
}
\usage{data(dynCorrData)}
\format{A 648 by 5 data frame.}
\keyword{datasets}
