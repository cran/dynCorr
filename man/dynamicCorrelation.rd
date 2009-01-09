\name{dynamicCorrelation}
\alias{dynamicCorrelation}
\title{Dynamic Correlation}
\description{
Computes dynamical correlation estimates for pairs of longitudinal responses, including consideration of lags and derivatives, following a local polynomial regression smoothing step.
}
\usage{
   dynamicCorrelation(dataFrame, depVar, indepVar, subjectVar, function.choice, width.range, width.place, boundary.trunc, lag.input, byOrder, by.deriv.only, full.lag.output)
}
\arguments{
   \item{dataFrame}{ The data frame that contains the dependent variables/responses, the independent variable (often time), and the subject/individual identification;there should be one row entry for each combination of subject/individual and indepVar (often time)}
   \item{depVar}{ Dependent variables/responses; at least two are necessary for purposes of calculating at least one dynamical correlation estimate of interest; there should be a unique column within depVar for each response}
   \item{indepVar}{ Independent variable, typically the discrete recorded time points at which the dependent variables were collected; note that this is the independent variable for purposes of curve creation leading into estimating the dynamical correlations between pairs of dependent variables; must be contained in a single column}
   \item{subjectVar}{ Column name of the individuals; there should be one row entry for each combination of subject/individual and indepVar}
   \item{function.choice}{ A vector of length 3 to indicate which derivatives desired for local polynomial curve smoothing; 1st entry is for 0th derivative (i.e., original function), 2nd entry is for 1st derivative, 3rd is for 2nd derivative; 1=yes, 0=no. e.g., c(1,0,1) would be specified if interest is in looking at derivative 0 and derivative 2, c(1,0,0) for looking at original function (0th derivative) only, etc.}
   \item{width.range}{Bandwidth for local polynomial regression curve smoothing for each dependent variable/response; it can be a list that specifies a distinct bandwidth two-element vector for each response, or a single two-element vector that is used for each of the responses --- the program is currently set up to allow linearly increasing or decreasing bandwidths, specified by this two-element vector with the increase (or decrease) occurring from the first argument in width.place to its second argument; the lpepa function within the lpridge package is called, which uses Epanecknikov kernel weighting, and the specifications of bandwidth in width.range will be used there; the default bandwidth is the range of indepVar (usually time) divided by 4, i.e., a constant global bandwidth }
   \item{width.place}{ Endpoints for width change assuming a non-constant bandwidth requested --- the program is currently set up to allow linearly increasing or decreasing bandwidths, specified in width.range, with the increase (or decrease) occurring from the first argument in width.place to its second argument; it can be a list that specifies different endpoints for each response, or a single vector of endpoints that is used for each of the responses; default is no endpoints specified which means constant global bandwidth throughout the range of indepVar}
   \item{boundary.trunc}{ Indicate the boundary of indepVar that should be truncated after smoothing; this may be done in case of concerns about estimating dynamical correlation at the boundaries of indepVar; this is a two element vector, where the first argument is how much to truncate from the right of the minimum value of indep.var and the second argument is how much to truncate from the left of the maximum value of indep var, within an individual and specific response; default is no truncation, i.e., c(0,0)}
   \item{lag.input}{ Values of lag to be considered; can be a vector of requested lags, for which a dynamical correlation estimate will be produced for each pair of responses at each lag; a positive value of lag.input means that the first entry for the dynamical correlation leads (occurs before) the second entry --- conversely, a negative value means that the second entry for the correlation leads the first entry; default is no lag at all considered}
   \item{byOrder}{ A vector that specifies the order of the variables/responses and derivatives (if any) to be the leading variable in the calculations; this will have an effect on how lag.input is to be interpreted; default is to use the order as specified in the depVar argument}
   \item{by.deriv.only}{ If TRUE, the inter-dynamical correlations between different derivatives are not computed, which can save computation time (e.g., when function.choice=c(1,0,1) is specified and by.deriv.only=T, then only dynamical correlations will be calculated within (and not across) the 0th and 2nd derivative estimates, respectively); default is TRUE }
   \item{full.lag.output}{ If TRUE, the dynamical correlation values for each pair of responses and requested derivative will be stored in vectors corresponding to different lag values, which enables plotting the correlations as a function of lag values; all the vectors will be stored in the returned attribute resultMatrix; default is FALSE }
}
\details{ This function will provide smooth estimates (curves/functions or their derivatives) of longitudinal responses of interest per individual, then generate dynamical correlation estimates for each pair of responses.  Lags of interest can be specified, using the lag.input argument.  For smoothing, the function uses local polynomial smoothing using Epanecknikov kernel weighting (by calling the function lpepa within the lpridge package).  The default global bandwidth is generated by taking the range of indepVar (usually time) and dividing by 4.  This, by default, will be a constant global bandwidth, but proper specification of the width.range and width.place arguments can allow for a more flexible bandwidth choice, including different specification for each response in depVar. 

Details of the methodology for dynamical correlation can be found in Dubin and Muller (2005).
}
\author{Joel A. Dubin, Dandi Qiao, Hans-Georg Muller}
\seealso{\link[lpridge]{lpepa}
}
\examples{

      data(dynCorrData)

      ## Example 1: using default smoothing parameters, obtain dynamical correlation estimates
      ##            for all three pairs of responses, for both original function and the first derivative 
      examp1 <- dynamicCorrelation(dataFrame=dynCorrData, 
                                   depVar=c('resp1', 'resp2', 'resp3'),
                                   indepVar='time',
                                   subjectVar = 'subject',
                                   function.choice = c(1,1,0))
      examp1

      ## Example 2: using default smoothing parameters, obtain dynamical correlation estimates
      ##             for all three pairs of responses, looking at range of lags between -10 and +10,
      ##             for original functions only 
      examp2 <- dynamicCorrelation(dataFrame=dynCorrData, 
                                   depVar=c('resp1', 'resp2', 'resp3'),
                                   indepVar='time',
                                   subjectVar = 'subject',
                                   function.choice = c(1,0,0),
                                   lag.input=seq(-20,20, by=1))
      examp2
      
      ## note: output includes zero lag correlations, as well as maximum correlation (in absolute value)
      ##             in max.dynCorr and and its corresponding lag value in max.dynCorrLag

      ## Example 3: re-rerun example 2, but set up for plotting of specified lagged correlations
      
      examp3 <- dynamicCorrelation(dataFrame=dynCorrData, 
                                   depVar=c('resp1', 'resp2', 'resp3'),
                                   indepVar='time',
                                   subjectVar = 'subject',
                                   function.choice = c(1,0,0),
                                   lag.input=seq(-20,10, by=1),
                                   full.lag.output=TRUE)

      # conduct plotting, with one panel for each pair of responses considered; the ylim adjustment
      #  is made here for the different magnitude of the correlations between the two pairs
      par(mfrow=c(1,2))
      plot(seq(-20,10, by=1), examp3$lagResultMatrix[[1]][1,], type='b',
	   xlab = 'lag order (in days)', ylab = 'lagged correlations', ylim = c(-0.4, -0.2),
	   main = 'dyncorr b/t resp1 and resp2 as a function of lag')
      abline(v = examp3$max.dynCorrLag[[1]][1,2], lty = 2)
      plot(seq(-20,10, by=1), examp3$lagResultMatrix[[1]][2,], type='b',
	   xlab = 'lag order (in days)', ylab = 'lagged correlations', ylim = c(0.3, 0.5),
	   main = 'dyncorr b/t resp1 and resp3 as a function of lag')
      abline(v = examp3$max.dynCorrLag[[1]][1,3], lty = 2)

      ## Example 4: same as the original function piece of Example 1, except now
      ##             adjust the constant global bandwidth from the default to 40
      
      examp4 <- dynamicCorrelation(dataFrame=dynCorrData,    
                                   depVar=c('resp1', 'resp2', 'resp3'),
                                   indepVar='time',
                                   subjectVar = 'subject',
                                   function.choice = c(1,0,0),
                                   width.range = c(40, 40) ) 
     examp4
}
\keyword{multivariate}
