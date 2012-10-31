\name{bootstrapCI}
\alias{bootstrapCI}
\title{Bootstrap Confidence Interval}
\description{ Computes percentile bootstrap (BS) confidence intervals for dynamical correlation for pairs of longitudinal responses, including consideration of lags and derivatives, following a local polynomial regression smoothing step.
}
\usage{ 
   bootstrapCI(dataFrame, depVar, indepVar, subjectVar, function.choice, width.range, width.place, boundary.trunc, byOrder, max.dynCorrLag, B, percentile, by.deriv.only, seed)
}
\arguments{
   \item{dataFrame}{ The data frame that contains the dependent variables/responses, the independent variable (often time), and the subject/individual identification; there should be one row entry for each combination of subject/individual and indepVar (often time).}
   \item{depVar}{ Dependent variables/responses; at least two are necessary for purposes of calculating at least one dynamical correlation estimate of interest; there should be a unique column within depVar for each response.}
   \item{indepVar}{ Independent variable, typically the discrete recorded time points at which the dependent variables were collected; note that this is the independent variable for purposes of curve creation leading into estimating the dynamical correlations between pairs of dependent variables; must be contained in a single column.}
   \item{subjectVar}{ Column name of the individuals; there should be one row entry for each combination of subject/individual and indepVar.}
   \item{function.choice}{ A vector of length 3 to indicate which derivatives desired for local polynomial curve smoothing; 1st entry is for 0th derivative (i.e., original function), 2nd entry is for 1st derivative, 3rd is for 2nd derivative; 1=yes, 0=no. e.g., c(1,0,1) would be specified if interest is in looking at derivative 0 and derivative 2, c(1,0,0) for looking at original function (0th derivative) only, etc. }
   \item{width.range}{ Bandwidth for local polynomial regression curve smoothing for each dependent variable/response; it can be a list that specifies a distinct bandwidth two-element vector for each response, or a single two-element vector that is used for each of the responses --- the program is currently set up to allow linearly increasing or decreasing bandwidths, specified by this two-element vector with the increase (or decrease) occurring from the first argument in width.place to its second argument; the lpepa function within the lpridge package is called, which uses Epanecknikov kernel weighting, and the specifications of bandwidth in width.range will be used there; the default bandwidth is the range of indepVar (usually time) divided by 4, i.e., a constant global bandwidth. }
   \item{width.place}{ Endpoints for width change assuming a non-constant bandwidth requested --- the program is currently set up to allow linearly increasing or decreasing bandwidths, specified in width.range, with the increase (or decrease) occurring from the first argument in width.place to its second argument; it can be a list that specifies different endpoints for each response, or a single vector of endpoints that is used for each of the responses; default is no endpoints specified which means constant global bandwidth throughout the range of indepVar.}
   \item{boundary.trunc}{ Indicate the boundary of indepVar that should be truncated after smoothing; this may be done in case of concerns about estimating dynamical correlation at the boundaries of indepVar; this is a two element vector, where the first argument is how much to truncate from the right of the minimum value of indep.var and the second argument is how much to truncate from the left of the maximum value of indep var, within an individual and specific response; default is no truncation, i.e., c(0,0).
}
   \item{byOrder}{ A vector that specifies the order of the variables/responses and derivatives (if any) to be the leading variable in the calculations; this will have an effect on how lag.input is to be interpreted; default is to use the order as specified in the depVar argument.}
   \item{max.dynCorrLag}{ A specified lag value at which you would like to obtain a percentile BS dynamical correlation interval estimate; only single lag values are allowed (due to computational considerations), and the lag value considered might be that which is output from a lag analysis using the dynamicCorrelation function 
(see Example 2 on the dynamicCorrelation help page).}
   \item{B}{ The number of samples used for the bootstrap. }
   \item{percentile}{ The percentile used to construct confidence intervals for dynamic correlations.}
   \item{by.deriv.only}{ If TRUE, the inter-dynamical correlations between different derivatives are not computed, which can save computation time (e.g., when function.choice=c(1,0,1) is specified and by.deriv.only=T, then only dynamical correlations will be calculated within (and not across) the 0th and 2nd derivative estimates, respectively); default is TRUE.  }
   \item{seed}{ The seed used in generating the BS samples.}
}
\details{ This function will provide smooth estimates (curves/functions or their derivatives) of responses of interest, then generate dynamical correlation estimates for each pair of responses.  Lags of interest can be specified, using the lag.input argument.  For smoothing, the function uses local polynomial smoothing using Epanecknikov kernel weighting (by calling the function lpepa within the lpridge package).  The default global bandwidth is generated by taking the range of indepVar (usually time) and dividing by 4.  This, by default, will be a constant global bandwidth, but proper specification of the width.range and width.place arguments can allow for a more flexible bandwidth choice, including different specification for each response in depVar.

Whereas the dynamicCorrelation program, which produces only point estimates, is fast, the bootstrapCI program is slow in its current form, as quite a bit of processing is required for each bootstrap sample. Details of the two-stage bootstrap algorithm can be found in Dubin and Muller (2005).  We will attempt to boost computing speed in future versions. 

In addition, as pointed out in Dubin and Muller (2005),
a downward shift of the two-stage bootstrap CI toward 0 is intentional, in order to account for error in the curve creation step.
However, it should be noted that greater than expected downward shifts may occur for dynamical correlations
that are high in magnitude.  A future version of this function will attempt to make this bootstrap approach more robust
in this situation.  
}

\author{ Joel A. Dubin,  Dandi Qiao, Hans-Georg Muller}
\seealso{\link[lpridge]{lpepa}
}
\examples{

      data(dynCorrData)

      ## Example 1: using default smoothing parameters, obtain bootstrap CI estimates
      ##            for all three pairs of responses, for original function only. 
      ##	    B=200 or greater should be considered, but not run here due to 	 
      ##	    computational time required.
      
      examp1.bs <- bootstrapCI(dataFrame=dynCorrData, 
                               depVar=c('resp1', 'resp2', 'resp3'),
                               indepVar='time',
                               subjectVar = 'subject',
                               function.choice = c(1,0,0),
                               B = 2, percentile=c(0.025, 0.975), seed = 5)
      examp1.bs

      ## Example 2: using default smoothing parameters, obtain bootstrap CI estimates
      ##             for all three pairs of responses, for original function only, at -10 
      ##            days lag, at .01 and .99 percentiles. 
      ##	    B=200 or greater should be considered, but not run here due to 	 
      ##	    computational time required.	
      
      examp2.bs <- bootstrapCI(dataFrame=dynCorrData, 
                               depVar=c('resp1', 'resp2', 'resp3'),
                               indepVar='time',
                               subjectVar = 'subject',
                               function.choice = c(1,0,0), max.dynCorrLag = -10,
                               B = 2, percentile=c(0.01, 0.99), seed = 7)
      examp2.bs
}
\keyword{multivariate 
}
