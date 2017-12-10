#' Bootstrap Confidence Interval
#'
#' Computes percentile bootstrap (BS) confidence intervals for dynamical
#' correlation for pairs of longitudinal responses, including consideration
#' of lags and derivatives, following a local polynomial regression smoothing
#' step.
#'
#' @usage
#' bootstrapCI(dataFrame, depVar, indepVar, subjectVar, function.choice,
#'             width.range, width.place, min.obs, points.length, points.by,
#'             boundary.trunc, byOrder, max.dynCorrLag, B, percentile,
#'             by.deriv.only, seed)
#'
#' @param dataFrame The data frame that contains the dependent
#'  variables/responses, the independent variable (often time),
#'  and the subject/individual identification; there should be
#'  one row entry for each combination of subject/individual and
#'  indepVar (often time).
#'
#' @param depVar Dependent variables/responses; at least two are
#'  necessary for purposes of calculating at least one dynamical
#'  correlation estimate of interest; there should be a unique
#'  column within depVar for each response.
#'
#' @param indepVar Independent variable, typically the discrete
#'  recorded time points at which the dependent variables were
#'  collected; note that this is the independent variable for
#'  purposes of curve creation leading into estimating the dynamical
#'  correlations between pairs of dependent variables; must be
#'  contained in a single column.
#'
#' @param subjectVar Column name of the individuals; there should be
#'  one row entry for each combination of subject/individual and indepVar.
#'
#' @param function.choice A vector of length 3 to indicate which
#'  derivatives desired for local polynomial curve smoothing;
#'  1st entry is for 0th derivative (i.e., original function),
#'  2nd entry is for 1st derivative, 3rd is for 2nd derivative;
#'  1=yes, 0=no. e.g., c(1,0,1) would be specified if interest
#'  is in looking at derivative 0 and derivative 2, c(1,0,0)
#'  for looking at original function (0th derivative) only, etc.
#'
#' @param width.range Bandwidth for local polynomial regression curve
#'  smoothing for each dependent variable/response; it can be a list
#'  that specifies a distinct bandwidth two-element vector for each
#'  response, or a single two-element vector that is used for each of
#'  the responses — the program is currently set up to allow linearly
#'  increasing or decreasing bandwidths, specified by this two-element
#'  vector with the increase (or decrease) occurring from the first argument
#'  in width.place to its second argument; the lpepa function within the
#'  lpridge package is called, which uses Epanecknikov kernel weighting,
#'  and the specifications of bandwidth in width.range will be used there;
#'  the default bandwidth is the range of indepVar (usually time) divided
#'  by 4, i.e., a constant global bandwidth.
#'
#' @param width.place Endpoints for width change assuming a non-constant
#'  bandwidth requested — the program is currently set up to allow linearly
#'  increasing or decreasing bandwidths, specified in width.range, with the
#'  increase (or decrease) occurring from the first argument in width.place
#'  to its second argument; it can be a list that specifies different
#'  endpoints for each response, or a single vector of endpoints that is
#'  used for each of the responses; default is no endpoints specified
#'  which means constant global bandwidth throughout the range of indepVar.
#'
#' @param min.obs Minimum oberservation (follow-up period) required.
#'  If specified, individuals whose follow-up period shorter than min.obs will
#'  be removed from calculation. Default = NA (use whole dataset provided).
#'
#' @param points.length Number of indep (time) points for dynamic correlation 
#'  calculation for each response of each individual. This is the number of 
#'  points between time 0 and minimum of maximum of individual follow-up 
#'  periods (max_common_obs). Note that each individual’s full follow-up 
#'  time span is used in the local polynomial regression curve smoothing 
#'  step, but only the first points.length number of time points (from time 0 
#'  to max_common_obs) is used in the following dynamic correlation calculation. 
#'  Default points.length value is set to 100; points.length takes precedence 
#'  unless points.by is specified.
#'
#' @param points.by Interval between indep (time) points for local polynomial
#'  regression curve smoothing for each response of each individual.
#'  Both integer and non-integer value could be specified, and time grid will
#'  be computed accordingly. Note that points.length takes precedence
#'  (default = 100) unless points.by is specified.
#'
#' @param boundary.trunc Indicate the boundary of indepVar that should be
#'  truncated after smoothing; this may be done in case of concerns about
#'  estimating dynamical correlation at the boundaries of indepVar; this
#'  is a two element vector, where the first argument is how much to
#'  truncate from the right of the minimum value of indep.var and the
#'  second argument is how much to truncate from the left of the maximum
#'  value of indep var, within an individual and specific response;
#'  default is no truncation, i.e., c(0,0).
#'
#' @param byOrder A vector that specifies the order of the
#'  variables/responses and derivatives (if any) to be the leading
#'  variable in the calculations; this will have an effect on how
#'  lag.input is to be interpreted; default is to use the order as
#'  specified in the depVar argument.
#'
#' @param max.dynCorrLag A specified lag value at which you would
#'  like to obtain a percentile BS dynamical correlation interval
#'  estimate; only single lag values are allowed (due to computational
#'  considerations), and the lag value considered might be that which
#'  is output from a lag analysis using the dynamicCorrelation function
#'  (see Example 2 on the dynamicCorrelation help page).
#'
#' @param B The number of samples used for the bootstrap.
#'
#' @param percentile The percentile used to construct confidence intervals
#'  for dynamic correlations.
#'
#' @param by.deriv.only If TRUE, the inter-dynamical correlations between
#'  different derivatives are not computed, which can save computation time
#'  (e.g., when function.choice=c(1,0,1) is specified and by.deriv.only=T,
#'  then only dynamical correlations will be calculated within
#'  (and not across) the 0th and 2nd derivative estimates, respectively);
#'  default is TRUE.
#'
#' @param seed The seed used in generating the BS samples.
#'
#' @details This function will provide smooth estimates
#'  (curves/functions or their derivatives) of responses
#'  of interest, then generate dynamical correlation estimates
#'  for each pair of responses. Lags of interest can be specified,
#'  using the lag.input argument. For smoothing, the function uses
#'  local polynomial smoothing using Epanecknikov kernel weighting
#'  (by calling the function lpepa within the lpridge package).
#'  The default global bandwidth is generated by taking the range
#'  of indepVar (usually time) and dividing by 4. This, by default,
#'  will be a constant global bandwidth, but proper specification of
#'  the width.range and width.place arguments can allow for a more
#'  flexible bandwidth choice, including different specification for
#'  each response in depVar. \cr \cr
#'  Whereas the dynamicCorrelation program, which produces only point
#'  estimates, is fast, the bootstrapCI program is slow in its current
#'  form, as quite a bit of processing is required for each bootstrap
#'  sample. Details of the two-stage bootstrap algorithm can be found
#'  in Dubin and Muller (2005). We will attempt to boost computing speed
#'  in future versions. \cr \cr
#'  In addition, as pointed out in Dubin and Muller (2005), a downward
#'  shift of the two-stage bootstrap CI toward 0 is intentional,
#'  in order to account for error in the curve creation step. However,
#'  it should be noted that greater than expected downward shifts may
#'  occur for dynamical correlations that are high in magnitude. A future
#'  version of this function will attempt to make this bootstrap approach
#'  more robust in this situation.
#'
#' @author Joel A. Dubin, Mike Li, Dandi Qiao, Hans-Georg Müller
#' @seealso \link[lpridge]{lpepa}
#' @examples
#'
#' ## Example 1: using default smoothing parameters, obtain bootstrap CI
#' ##            estimates for all three pairs of responses, for original
#' ##            function only. Note that B=200 or greater should be
#' ##            considered for real data analysis.
#'
#' examp1.bs <- bootstrapCI(dataFrame = dynCorrData,
#'                          depVar = c('resp1', 'resp2', 'resp3'),
#'                          indepVar = 'time',
#'                          subjectVar = 'subject',
#'                          points.by = 1,
#'                          function.choice = c(1,0,0),
#'                          B = 2, percentile = c(0.025, 0.975), seed = 5)
#' examp1.bs
#'
#' ## Example 2: using default smoothing parameters, obtain bootstrap CI
#' ##            estimates for all three pairs of responses, for original
#' ##            function only, at -10 days lag, at .01 and .99 percentiles.
#' ##            Note that B=200 or greater should be considered for real
#' ##            data analysis.
#'
#' examp2.bs <- bootstrapCI(dataFrame = dynCorrData,
#'                          depVar = c('resp1', 'resp2', 'resp3'),
#'                          indepVar = 'time',
#'                          subjectVar = 'subject',
#'                          points.by = 1,
#'                          function.choice = c(1,0,0), max.dynCorrLag = -10,
#'                          B = 2, percentile = c(0.01, 0.99), seed = 7)
#' examp2.bs
#'
#' @export
#'
#' @importFrom lpridge lpepa
#' @importFrom stats quantile var

bootstrapCI <- function (dataFrame,
                         depVar = c("resp1", "resp2", "resp3", "resp4", "resp5"),
                         indepVar = "time",
                         subjectVar = "subject",
                         function.choice = c(1),
                         width.range = c(((range(dataFrame[[indepVar]]))[2] -
                                            (range(dataFrame[[indepVar]]))[1])/4,
                                         ((range(dataFrame[[indepVar]]))[2] -
                                            (range(dataFrame[[indepVar]]))[1])/4),
                         width.place = c(NA, NA),
                         min.obs = NA,
                         points.length = 100,
                         points.by = NA,
                         boundary.trunc = c(0, 0),
                         byOrder = c(),
                         max.dynCorrLag = 0,
                         B = 100,
                         percentile = c(0.025, 0.975),
                         by.deriv.only = TRUE,
                         seed = 299)
{
  ## --------------------Set up -------------------

  # Nov-2017 Introduce new parametew min.obs (modify dataset if min.obs
  # is specified)
  if (!is.na(min.obs)) {
    vec <- unique(dataFrame[[subjectVar]]) # list of individuals id
    dep_var = dataFrame[depVar]
    subject_var = dataFrame[[subjectVar]]
    indep_var = dataFrame[[indepVar]]
    
    to_remove <- c()
    for (i in 1:length(vec)) {
      ind_time <-max(indep_var[subject_var == vec[i]])
      if (ind_time < min.obs) {
        to_remove <- c(to_remove, vec[i])
      }
    }
    
    for (i in 1:length(to_remove)) {
      dataFrame <- dataFrame[subject_var!=to_remove[i],]
      subject_var <-subject_var[subject_var!=to_remove[i]]
      indep_var <- indep_var[subject_var!=to_remove[i]]
    }
  }

  vec <- unique(dataFrame[[subjectVar]]) # list of individuals id
  dep_var = dataFrame[depVar]
  subject_var = dataFrame[[subjectVar]]
  indep_var = dataFrame[[indepVar]]

  # process function choice
  funcVar <- c()
  for (i in 1:length(function.choice)) {
    if (function.choice[i] == 1) funcVar <- c(funcVar, i)
  }
  num_funcVar <- length(funcVar) # number of function methods
  num_vec <- length(vec) # num of individuals
  num_depVar <- length(depVar) # number of dependent Variables
  l_trunc <- boundary.trunc[1] # lower boundary
  h_trunc <- boundary.trunc[2] # upper boundary

  # check input byOrder
  if (length(byOrder) != 0) {
    largest <- num_depVar * num_funcVar
    if (max(byOrder) > largest) {
      stop("The specified index order of leading variables is over range")
    }
    if (length(byOrder) != largest && length(byOrder) != 0) {
      stop("The index order of leading variables is not completely specified")
    }
  }

  # check inputs width.range and width.point.
  # If they are invalid, throw exception
  if (!is.list(width.range)) {
    temp_widthRange <- width.range
    width.range <- vector(mode = "list", num_depVar)
    for (i in 1:num_depVar) {
      width.range[[i]] <- temp_widthRange
    }
  }
  if (!is.list(width.place)) {
    temp_widthPlace <- width.place
    width.place <- vector(mode = "list", num_depVar)
    for (i in 1:num_depVar) {
      width.place[[i]] <- temp_widthPlace
    }
  }
  if (length(width.range) != num_depVar) {
    stop("If width.range is a list, it needs to have the same number of
         components as the number of responses")
  }
  if (length(width.place) != num_depVar) {
    stop("If width.place is a list, it needs to have the same number of
         components as the number of responses")
  }

  # find max common observation
  v_ob_time <- c()
  for (i in 1:num_vec) {
    v_ob_time <- c(v_ob_time, max(indep_var[subject_var == vec[i]]))
  }
  limit <- min(v_ob_time) # max common obs
  

  # Sep 2017 - calculate points.by if it is not specified
  if (is.na(points.by)) {
    points.by <- limit/points.length
  }

  # Nov 2017 - add output table
  max_limit <- max(v_ob_time)
  base <- min(indep_var[subject_var == vec[1]])
  max_limit <- max_limit - base
  min_max_limit <- limit - base
  n <- length(vec)
  data_summary <- matrix(c(n,min_max_limit,max_limit), ncol=3, byrow=TRUE)
  colnames(data_summary) <- c('sample.size','min.max.time','max.max.time')
  
  ## ------------ Smooth Curves -----------------------
  smoothedCurves <- vector(mode = "list", num_vec)
  for (k in 1:num_vec) {
    indep <- indep_var[subject_var == vec[k]]
    smoothedCurves[[k]] <- vector(mode = "list", num_depVar)
    max_indep <- max(indep[!is.na(indep)])

    for (i in 1:num_depVar) {
      cur_wrange <- width.range[[i]]
      cur_wplace <- width.range[[i]]

      if (is.na(width.place[[i]][1])) {
        cur_wplace[1] <- min(indep[!is.na(indep)])
      }
      if (is.na(width.place[[i]][2])) {
        cur_wplace[2] <- max_indep
      }

      # Aug-2017 Update: change points.use to accommodate non-integer time grid
      # prior to v1.0.0, points.use assumes integer values
      # points.use = seq(l_trunc, ceiling(max_indep - h_trunc), by = 1)
      # size <- ceiling(max_indep - h_trunc) - l_trunc + 1
      points.use <- seq(l_trunc, ceiling(max_indep - h_trunc), by=points.by)

      size <- length(points.use)
      band <- c()
      for (count in 1:size) {
        num_points <- points.use[count]
        widthrange <- cur_wrange[2] - cur_wrange[1]

        if (abs(widthrange) < 1e-05) width <- cur_wrange[1]
        else if (num_points < cur_wplace[1]) width <- cur_wrange[1]
        else if (num_points >= cur_wplace[1] & num_points <= cur_wplace[2]) {
          width <- cur_wrange[1] +
            ((num_points - cur_wplace[1]) / (cur_wplace[2] - cur_wplace[1])) *
            widthrange
        } else width <- cur_wrange[2]
        band <- c(band, width)
      }

      curve_x <- indep_var[subject_var == vec[k]][!is.na(
        dep_var[[i]][subject_var == vec[k]])]
      curve_y <- dep_var[[i]][subject_var == vec[k]][!is.na(
        dep_var[[i]][subject_var == vec[k]])]

      temp <- lpepa(x = curve_x, y = curve_y, bandwidth = band,
                    deriv = 0, n.out = size, x.out = points.use,
                    order = 1, var = FALSE)$est
      smoothedCurves[[k]][[i]] <- vector(mode = "list", length = 2)
      smoothedCurves[[k]][[i]][[1]] <- points.use
      smoothedCurves[[k]][[i]][[2]] <- temp
    }
  }
  days.out <- l_trunc:(ceiling(limit) - h_trunc)
  resid <- vector(mode = "list", num_vec)
  smooth <- vector(mode = "list", num_vec)

  for (k in 1:num_vec) {
    resid[[k]] <- vector(mode = "list", num_depVar)
    smooth[[k]] <- vector(mode = "list", num_depVar)
    for (i in 1:num_depVar) {
      indep <-
        indep_var[subject_var ==
                    vec[k]][!is.na(dep_var[[i]][subject_var == vec[k]])]

      # fixed code
      # prior to version 1.1, calculation uses
      # indep[indep >= l_trunc & indep <= ceiling(max(indep)) - h_trunc]
      temp.obs <-
        indep[indep >= l_trunc & indep <= ceiling(max(indep) - h_trunc)]

      dep.v <-
        dep_var[[i]][subject_var ==
                       vec[k]][!is.na(dep_var[[i]][subject_var == vec[k]])]

      depValue <- dep.v[indep >= l_trunc & indep <=
                          ceiling(max(indep[!is.na(indep)])) - h_trunc]
      temp.obs <- round(temp.obs, 0)
      smooth[[k]][[i]] <- vector(mode = "numeric", length(temp.obs))
      for (j in 1:length(temp.obs)) {
        smooth[[k]][[i]][[j]] <-
          smoothedCurves[[k]][[i]][[2]][
            is.element(smoothedCurves[[k]][[i]][[1]], temp.obs[j])]
      }
      resid[[k]][[i]] <- depValue - smooth[[k]][[i]]
    }
  }
  rm(temp.obs)
  bs.seq.samp <- vector(mode = "list", length = B)
  set.seed(seed)
  indexVec <- 1:num_vec
  for (b in 1:B) {
    bs.seq.samp[[b]] <- sort(sample(indexVec, replace = TRUE))
  }
  recreate.B.list <- vector(mode = "list", length(B))
  for (b in 1:B) {
    vecSample <- bs.seq.samp[[b]]
    recreate.B.list[[b]] <- vector(mode = "list", length(vecSample))
    for (j in 1:length(vecSample)) {
      recreate.B.list[[b]][[j]] <- vector(mode = "list",
                                          num_depVar)
      for (i in 1:num_depVar) {
        recreate.B.list[[b]][[j]][[i]] <-
          sample(resid[[vecSample[j]]][[i]], replace = TRUE) +
          smooth[[vecSample[j]]][[i]]
      }
    }
  }
  days.out <- l_trunc:(ceiling(limit) - h_trunc)
  temp.days <- tempData <- vector(mode = "list", length(B))
  for (b in 1:B) {
    vecSample <- bs.seq.samp[[b]]
    tempData[[b]] <- vector(mode = "list", length(vecSample))
    temp.days[[b]] <- vector(mode = "list", length(vecSample))
    for (j in 1:length(vecSample)) {
      tempData[[b]][[j]] <- vector(mode = "list", num_depVar)
      temp.days[[b]][[j]] <- vector(mode = "list", num_depVar)
      for (i in 1:num_depVar) {
        indep <-
          indep_var[subject_var ==
                      vec[vecSample[j]]][!is.na(dep_var[[i]][
                        subject_var == vec[vecSample[j]]])]
        temp.1 <- indep[indep >= 0 & indep < l_trunc]
        temp.2 <- indep[indep > ceiling(max(indep[!is.na(indep)])) -
                          h_trunc & indep <= max(indep[!is.na(indep)])]
        temp.days[[b]][[j]][[i]] <-
          indep[indep >= 0 & indep <= max(indep[!is.na(indep)])]
        depValue <- dep_var[[i]][subject_var == vec[vecSample[j]]]
        depV1 <- depValue[is.element(indep_var[subject_var ==
                                                 vec[vecSample[j]]], temp.1)]
        depV2 <- depValue[is.element(indep_var[subject_var ==
                                                 vec[vecSample[j]]], temp.2)]
        if (length(temp.2) > 0) {
          tempData[[b]][[j]][[i]] <- c(depV1, recreate.B.list[[b]][[j]][[i]],
                                       depV2)
        }
        else if (length(temp.1) > 0) {
          tempData[[b]][[j]][[i]] <- c(depV1, recreate.B.list[[b]][[j]][[i]])
        }
        else {
          tempData[[b]][[j]][[i]] <- recreate.B.list[[b]][[j]][[i]]
        }
      }
    }
  }
  len <- num_vec
  cov.wgt.mtx.B <- vector(mode = "list", length = B)
  for (b in 1:B) {
    cat(" B =", b, "\n")
    vecSample <- bs.seq.samp[[b]]
    smoothedCurves2 <- vector(mode = "list", length(vecSample))
    for (k in 1:length(vecSample)) {
      smoothedCurves2[[k]] <- vector(mode = "list", num_depVar)
      for (i in 1:num_depVar) {
        indep <- temp.days[[b]][[k]][[i]]
        max_indep <- max(indep[!is.na(indep)])
        if (is.na(width.place[[i]][1])) {
          width.place[[i]][1] <- min(indep[!is.na(indep)])
        }
        if (is.na(width.place[[i]][2])) {
          width.place[[i]][2] <- max_indep
        }

        # Sep-2017 Update: change points.use to accommodate non-integer time grid
        # prior to v1.0.0, points.use assumes integer values
        # points.use = seq(l_trunc, ceiling(max_indep - h_trunc), by = 1)
        # size <- ceiling(max_indep - h_trunc) - l_trunc + 1
        points.use = seq(l_trunc, ceiling(max_indep - h_trunc), by=points.by)
        size <- length(points.use)
        band <- c()
        for (count in 1:size) {
          if (abs(width.range[[i]][1] - width.range[[i]][2]) <
              1e-05)
            width <- width.range[[i]][1]
          else {
            widthrange <- width.range[[i]][2] - width.range[[i]][1]
            if (points.use[count] < width.place[[i]][1])
              width <- width.range[[i]][1]
            else if (points.use[count] >= width.place[[i]][1] &
                     points.use[count] <= width.place[[i]][2])
              width <- width.range[[i]][1] +
                ((points.use[count] - width.place[[i]][1])/
                   (width.place[[i]][2] - width.place[[i]][1])) * widthrange
            else width <- width.range[[i]][2]
          }
          band <- c(band, width)
        }
        smoothedCurves2[[k]][[i]] <- vector(mode = "list", num_funcVar)
        for (j in 1:num_funcVar) {
          tempFunc <- lpepa(x = indep,
                            y = tempData[[b]][[k]][[i]],
                            bandwidth = band,
                            deriv = (funcVar[j] - 1),
                            n.out = size,
                            x.out = points.use,
                            order = funcVar[j], var = FALSE)$est
          smoothedCurves2[[k]][[i]][[j]] <- tempFunc
        }
      }
    }
    # Sep-2017 update: max_len is always calculated by points.by
    max_len <- length(seq(l_trunc, ceiling(limit) - h_trunc, by=points.by))

    meanMatrix <- list()
    for (i in 1:num_depVar) {
      forEachDepMean <- list()
      for (j in 1:num_funcVar) {
        forEachDerivMean <- rep(NA, max_len)
        for (m in 1:max_len) {
          for (n in 1:length(vecSample)) {
            if (n == 1) {
              total <- 1
              forEachDerivMean[m] <- smoothedCurves2[[n]][[i]][[j]][m]
            }
            else if (n > 1) {
              total <- total + 1
              forEachDerivMean[m] <-
                (1/total) * ((total - 1) * forEachDerivMean[m]
                             + smoothedCurves2[[n]][[i]][[j]][m])
            }
          }
        }
        forEachDepMean[[j]] <- forEachDerivMean
      }
      meanMatrix[[i]] <- forEachDepMean
    }
    if (by.deriv.only == FALSE) {
      dim <- num_depVar * num_funcVar
      cov.lag.mtx.listz <- list()
      weights.vec <- rep(NA, length(vecSample))
      time.extend <- (l_trunc + max_len - 1) +
        l_trunc
      for (i in 1:length(vecSample)) {
        cov.lag.mtx.listz[[i]] <- list()
        weights.vec[i] <-
          length(subject_var[(subject_var == vec[vecSample[i]]) &
                               (indep_var >= 0) &
                               (indep_var <= time.extend)])
        dep.correct <- list()
        for (m in 1:num_depVar) {
          dep.correct.function <- list()
          for (n in 1:num_funcVar) {
            dep.correct.function[[n]] <-
              smoothedCurves2[[i]][[m]][[n]][1:max_len] - meanMatrix[[m]][[n]]
          }
          dep.correct[[m]] <- dep.correct.function
        }
        cov.lag.mtx.listz[[i]] <- matrix(nrow = dim, ncol = dim)
        diag(cov.lag.mtx.listz[[i]]) <- 1
        if (max.dynCorrLag >= 0) {
          lag.end <- max_len - max.dynCorrLag
          lag.beg <- 1 + max.dynCorrLag
          m.support <- lag.end
          if (length(byOrder) == 0) {
            correst.standz <- list()
            index <- 1
            mean <- mean(dep.correct[[1]][[1]][1:lag.end])
            sd <- sqrt(var(dep.correct[[1]][[1]][1:lag.end]))
            correst.standz[[index]] <-
              (dep.correct[[1]][[1]][1:lag.end] - mean)/
              (sqrt((m.support - 1)/m.support) * sd)
            for (index in 2:dim) {
              if (dim >= 2) {
                dep <- (index - 1)%/%num_funcVar + 1
                deriv <- (index - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
                correst.standz[[index]] <-
                  (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                cov.lag.mtx.listz[[i]][1, index] <-
                  cov.lag.mtx.listz[[i]][index, 1] <-
                  (1/m.support) *
                  sum(correst.standz[[1]] * correst.standz[[index]])
              }
            }
            for (d in 2:(dim - 1)) {
              if (dim >= 3) {
                dep <- (d - 1)%/%num_funcVar + 1
                deriv <- (d - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                correst.standz[[d]] <-
                  (dep.correct[[dep]][[deriv]][1:lag.end] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                for (e in (d + 1):dim) {
                  cov.lag.mtx.listz[[i]][d, e] <-
                    cov.lag.mtx.listz[[i]][e, d] <-
                    (1/m.support) *
                    sum(correst.standz[[d]] * correst.standz[[e]])
                }
              }
            }
          }
          else {
            correst.standz <- list()
            index <- byOrder[1]
            dep <- (index - 1)%/%num_funcVar + 1
            deriv <- (index - 1)%%num_funcVar + 1
            mean <- mean(dep.correct[[dep]][[deriv]][1:lag.end])
            sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
            correst.standz[[index]] <-
              (dep.correct[[dep]][[deriv]][1:lag.end] - mean)/
              (sqrt((m.support - 1)/m.support) * sd)
            for (ord in 2:dim) {
              if (dim >= 2) {
                index <- byOrder[ord]
                dep <- (index - 1)%/%num_funcVar +
                  1
                deriv <- (index - 1)%%num_funcVar +
                  1
                mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
                correst.standz[[index]] <-
                  (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                cov.lag.mtx.listz[[i]][byOrder[1], index] <-
                  cov.lag.mtx.listz[[i]][index, byOrder[1]] <-
                  (1/m.support) *
                  sum(correst.standz[[byOrder[1]]] * correst.standz[[index]])
              }
            }
            for (ord in 2:(dim - 1)) {
              if (dim >= 3) {
                d <- byOrder[ord]
                dep <- (d - 1)%/%num_funcVar + 1
                deriv <- (d - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                correst.standz[[d]] <-
                  (dep.correct[[dep]][[deriv]][1:lag.end] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                for (ord2 in (ord + 1):dim) {
                  e <- byOrder[ord2]
                  cov.lag.mtx.listz[[i]][d, e] <-
                    cov.lag.mtx.listz[[i]][e, d] <-
                    (1/m.support) *
                    sum(correst.standz[[d]] * correst.standz[[e]])
                }
              }
            }
          }
        }
        else if (max.dynCorrLag < 0) {
          lag.end <- max_len + max.dynCorrLag
          lag.beg <- 1 - max.dynCorrLag
          m.support <- lag.end
          if (length(byOrder) == 0) {
            correst.standz <- list()
            index <- 1
            mean <- mean(dep.correct[[1]][[1]][lag.beg:max_len])
            sd <- sqrt(var(dep.correct[[1]][[1]][lag.beg:max_len]))
            correst.standz[[index]] <-
              (dep.correct[[1]][[1]][lag.beg:max_len] - mean)/
              (sqrt((m.support - 1)/m.support) * sd)
            for (index in 2:dim) {
              if (dim >= 2) {
                dep <- (index - 1)%/%num_funcVar + 1
                deriv <- (index - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                correst.standz[[index]] <-
                  (dep.correct[[dep]][[deriv]][1:lag.end] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                cov.lag.mtx.listz[[i]][1, index] <-
                  cov.lag.mtx.listz[[i]][index, 1] <-
                  (1/m.support) *
                  sum(correst.standz[[1]] * correst.standz[[index]])
              }
            }
            for (d in 2:(dim - 1)) {
              if (dim >= 3) {
                dep <- (d - 1)%/%num_funcVar + 1
                deriv <- (d - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
                correst.standz[[d]] <-
                  (dep.correct[[dep]][[deriv]][lag.beg:max_len] -  mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                for (e in (d + 1):dim) {
                  cov.lag.mtx.listz[[i]][d, e] <-
                    cov.lag.mtx.listz[[i]][e, d] <-
                    (1/m.support) *
                    sum(correst.standz[[d]] * correst.standz[[e]])
                }
              }
            }
          }
          else {
            correst.standz <- list()
            index <- byOrder[1]
            dep <- (index - 1)%/%num_funcVar + 1
            deriv <- (index - 1)%%num_funcVar + 1
            mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
            sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
            correst.standz[[index]] <-
              (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean)/
              (sqrt((m.support - 1)/m.support) * sd)
            for (ord in 2:dim) {
              index <- byOrder[ord]
              if (dim >= 2) {
                dep <- (index - 1)%/%num_funcVar +
                  1
                deriv <- (index - 1)%%num_funcVar +
                  1
                mean <- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                correst.standz[[index]] <-
                  (dep.correct[[dep]][[deriv]][1:lag.end] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                cov.lag.mtx.listz[[i]][1, index] <-
                  cov.lag.mtx.listz[[i]][index, 1] <-
                  (1/m.support) *
                  sum(correst.standz[[1]] * correst.standz[[index]])
              }
            }
            for (ord in 2:(dim - 1)) {
              if (dim >= 3) {
                d <- byOrder[ord]
                dep <- (d - 1)%/%num_funcVar + 1
                deriv <- (d - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
                correst.standz[[d]] <-
                  (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean)/
                  (sqrt((m.support - 1)/m.support) * sd)
                for (ord2 in (ord + 1):dim) {
                  e <- byOrder[ord2]
                  cov.lag.mtx.listz[[i]][d, e] <-
                    cov.lag.mtx.listz[[i]][e, d] <-
                    (1/m.support) *
                    sum(correst.standz[[d]] * correst.standz[[e]])
                }
              }
            }
          }
        }
      }
      cov.lag.wgt.mtx <- matrix(0, nrow = dim, ncol = dim)
      for (i in 1:length(vecSample)) {
        cov.lag.wgt.mtx <-
          cov.lag.wgt.mtx + (weights.vec[i] * cov.lag.mtx.listz[[i]])
      }
      cov.lag.wgt.mtx <- (1/sum(weights.vec)) * cov.lag.wgt.mtx
      diag(cov.lag.wgt.mtx) <- 1
      names <- c()
      for (i in 1:length(dep_var)) {
        name <- depVar[i]
        for (j in 1:num_funcVar) {
          name2 <- paste(name, (funcVar[j] - 1))
          names <- c(names, name2)
        }
      }
      dimnames(cov.lag.wgt.mtx) <- list(names, names)
      cov.wgt.mtx.B[[b]] <- cov.lag.wgt.mtx
    }
    else {
      weights.vec <- rep(NA, length(vecSample))
      time.extend <- (l_trunc + max_len - 1) +
        l_trunc
      cov.lag.mtx.listz <- vector(mode = "list", num_vec)
      dim <- num_depVar
      for (i in 1:length(vecSample)) {
        weights.vec[i] <-
          length(subject_var[(subject_var == vec[vecSample[i]]) &
                               (indep_var >= 0) &
                               (indep_var <= time.extend)])
        cov.lag.mtx.listz[[i]] <- vector(mode = "list", num_funcVar)
        dep.correct <- list()
        for (m in 1:num_depVar) {
          dep.correct.function <- list()
          for (n in 1:num_funcVar) {
            dep.correct.function[[n]] <-
              smoothedCurves2[[i]][[m]][[n]][1:max_len] -  meanMatrix[[m]][[n]]
          }
          dep.correct[[m]] <- dep.correct.function
        }
        for (deriv in 1:num_funcVar) {
          cov.lag.mtx.listz[[i]][[deriv]] <- matrix(nrow = dim, ncol = dim)
          diag(cov.lag.mtx.listz[[i]][[deriv]]) <- 1
          if (max.dynCorrLag >= 0) {
            lag.end <- max_len - max.dynCorrLag
            lag.beg <- 1 + max.dynCorrLag
            m.support <- lag.end
            if (length(byOrder) == 0) {
              correst.standz <- list()
              index <- 1
              mean <- mean(dep.correct[[index]][[deriv]][1:lag.end])
              sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
              correst.standz[[index]] <-
                (dep.correct[[index]][[deriv]][1:lag.end] - mean)/
                (sqrt((m.support - 1)/m.support) * sd)
              for (index in 2:dim) {
                if (dim >= 2) {
                  mean <- mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
                  sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
                  correst.standz[[index]] <-
                    (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  cov.lag.mtx.listz[[i]][[deriv]][1, index] <-
                    cov.lag.mtx.listz[[i]][[deriv]][index, 1] <-
                    (1/m.support) *
                    sum(correst.standz[[1]] * correst.standz[[index]])
                }
              }
              for (d in 2:(dim - 1)) {
                if (dim >= 3) {
                  mean <- mean(dep.correct[[d]][[deriv]][1:lag.end])
                  sd <- sqrt(var(dep.correct[[d]][[deriv]][1:lag.end]))
                  correst.standz[[d]] <-
                    (dep.correct[[d]][[deriv]][1:lag.end] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  for (e in (d + 1):dim) {
                    cov.lag.mtx.listz[[i]][[deriv]][d, e] <-
                      cov.lag.mtx.listz[[i]][[deriv]][e, d] <-
                      (1/m.support) *
                      sum(correst.standz[[d]] * correst.standz[[e]])
                  }
                }
              }
            }
            else {
              byOrderPerD <- c()
              for (ite in 1:length(byOrder)) {
                dep <- (byOrder[ite] - 1)%/%num_funcVar + 1
                deriv2 <- (byOrder[ite] - 1)%%num_funcVar + 1
                if (deriv2 == deriv) byOrderPerD <- c(byOrderPerD, dep)
              }
              if (length(byOrderPerD) != num_depVar) {
                stop("**********Error for the length of the order of
                     variables for each derivatives")
              }
              correst.standz <- list()
              index <- byOrderPerD[1]
              mean <- mean(dep.correct[[index]][[deriv]][1:lag.end])
              sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
              correst.standz[[index]] <-
                (dep.correct[[index]][[deriv]][1:lag.end] - mean)/
                (sqrt((m.support - 1)/m.support) * sd)
              for (ord in 2:dim) {
                if (dim >= 2) {
                  index <- byOrderPerD[ord]
                  mean <- mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
                  sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
                  correst.standz[[index]] <-
                    (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  cov.lag.mtx.listz[[i]][[deriv]][byOrderPerD[1], index] <-
                    cov.lag.mtx.listz[[i]][[deriv]][index, byOrderPerD[1]] <-
                    (1/m.support) *
                    sum(correst.standz[[byOrderPerD[1]]] *
                          correst.standz[[index]])
                }
              }
              for (ord in 2:(dim - 1)) {
                if (dim >= 3) {
                  d <- byOrderPerD[ord]
                  mean <- mean(dep.correct[[d]][[deriv]][1:lag.end])
                  sd <- sqrt(var(dep.correct[[d]][[deriv]][1:lag.end]))
                  correst.standz[[d]] <-
                    (dep.correct[[d]][[deriv]][1:lag.end] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  for (ord2 in (ord + 1):dim) {
                    e <- byOrderPerD[ord2]
                    cov.lag.mtx.listz[[i]][[deriv]][d, e] <-
                      cov.lag.mtx.listz[[i]][[deriv]][e, d] <-
                      (1/m.support) *
                      sum(correst.standz[[d]] * correst.standz[[e]])
                  }
                }
              }
              }
          }
          else if (max.dynCorrLag < 0) {
            lag.end <- max_len + max.dynCorrLag
            lag.beg <- 1 - max.dynCorrLag
            m.support <- lag.end
            if (length(byOrder) == 0) {
              correst.standz <- list()
              index <- 1
              mean <- mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
              sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
              correst.standz[[index]] <-
                (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean)/
                (sqrt((m.support - 1)/m.support) * sd)
              for (index in 2:dim) {
                if (dim >= 2) {
                  mean <- mean(dep.correct[[index]][[deriv]][1:lag.end])
                  sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
                  correst.standz[[index]] <-
                    (dep.correct[[index]][[deriv]][1:lag.end] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  cov.lag.mtx.listz[[i]][[deriv]][1, index] <-
                    cov.lag.mtx.listz[[i]][[deriv]][index, 1] <-
                    (1/m.support) *
                    sum(correst.standz[[1]] * correst.standz[[index]])
                }
              }
              for (d in 2:(dim - 1)) {
                if (dim >= 3) {
                  mean <- mean(dep.correct[[d]][[deriv]][lag.beg:max_len])
                  sd <- sqrt(var(dep.correct[[d]][[deriv]][lag.beg:max_len]))
                  correst.standz[[d]] <-
                    (dep.correct[[d]][[deriv]][lag.beg:max_len] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  for (e in (d + 1):dim) {
                    cov.lag.mtx.listz[[i]][[deriv]][d, e] <-
                      cov.lag.mtx.listz[[i]][[deriv]][e, d] <-
                      (1/m.support) *
                      sum(correst.standz[[d]] * correst.standz[[e]])
                  }
                }
              }
            }
            else {
              byOrderPerD <- c()
              for (ite in 1:length(byOrder)) {
                dep <- (byOrder[ite] - 1)%/%num_funcVar +
                  1
                deriv2 <- (byOrder[ite] - 1)%%num_funcVar +
                  1
                if (deriv2 == deriv) {
                  byOrderPerD <- c(byOrderPerD, dep)
                }
              }
              if (length(byOrderPerD) != num_depVar) {
                stop("**********Error for the length of the order
                     of variables for each derivatives")
              }
              correst.standz <- list()
              index <- byOrderPerD[1]
              mean <- mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
              sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
              correst.standz[[index]] <-
                (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean)/
                (sqrt((m.support - 1)/m.support) * sd)
              for (ord in 2:dim) {
                index <- byOrderPerD[ord]
                if (dim >= 2) {
                  mean <- mean(dep.correct[[index]][[deriv]][1:lag.end])
                  sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
                  correst.standz[[index]] <-
                    (dep.correct[[index]][[deriv]][1:lag.end] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  cov.lag.mtx.listz[[i]][[deriv]][byOrderPerD[1], index] <-
                    cov.lag.mtx.listz[[i]][[deriv]][index, byOrderPerD[1]] <-
                    (1/m.support) *
                    sum(correst.standz[[1]] * correst.standz[[index]])
                }
              }
              for (ord in 2:(dim - 1)) {
                if (dim >= 3) {
                  d <- byOrderPerD[ord]
                  mean <- mean(dep.correct[[d]][[deriv]][lag.beg:max_len])
                  sd <- sqrt(var(dep.correct[[d]][[deriv]][lag.beg:max_len]))
                  correst.standz[[d]] <-
                    (dep.correct[[d]][[deriv]][lag.beg:max_len] - mean)/
                    (sqrt((m.support - 1)/m.support) * sd)
                  for (ord2 in (ord + 1):dim) {
                    e <- byOrderPerD[ord2]
                    cov.lag.mtx.listz[[i]][[deriv]][d, e] <-
                      cov.lag.mtx.listz[[i]][[deriv]][e, d] <-
                      (1/m.support) *
                      sum(correst.standz[[d]] * correst.standz[[e]])
                  }
                }
              }
              }
          }
        }
      }
      cov.lag.wgt.mtx <- vector(mode = "list", num_funcVar)
      for (deriv in 1:num_funcVar) {
        cov.lag.wgt.mtx[[deriv]] <- matrix(0, nrow = dim,
                                           ncol = dim)
        for (i in 1:length(vecSample)) {
          cov.lag.wgt.mtx[[deriv]] <-
            cov.lag.wgt.mtx[[deriv]] +
            (weights.vec[i] * cov.lag.mtx.listz[[i]][[deriv]])
        }
        cov.lag.wgt.mtx[[deriv]] <-
          (1/sum(weights.vec)) * cov.lag.wgt.mtx[[deriv]]
      }
      for (deriv in 1:num_funcVar) {
        names <- c()
        for (i in 1:length(dep_var)) {
          name <- depVar[i]
          name2 <- paste(name, (funcVar[deriv] - 1))
          names <- c(names, name2)
        }
        dimnames(cov.lag.wgt.mtx[[deriv]]) <- list(names, names)
      }
      cov.wgt.mtx.B[[b]] <- cov.lag.wgt.mtx
    }
  }
  if (by.deriv.only == FALSE) {
    dim <- num_depVar * num_funcVar
    mean.cov.wgt.mtx.B <- matrix(0, nrow = dim, ncol = dim)
    for (b in 1:B) {
      mean.cov.wgt.mtx.B <- mean.cov.wgt.mtx.B + cov.wgt.mtx.B[[b]]
    }
    mean.cov.wgt.mtx.B <- mean.cov.wgt.mtx.B/B
  }
  else {
    mean.cov.wgt.mtx.B <- vector(mode = "list", num_funcVar)
    dim <- num_depVar
    for (deriv in 1:num_funcVar) {
      mean.cov.wgt.mtx.B[[deriv]] <- matrix(0, nrow = dim, ncol = dim)
      for (b in 1:B) {
        mean.cov.wgt.mtx.B[[deriv]] <-
          mean.cov.wgt.mtx.B[[deriv]] + cov.wgt.mtx.B[[b]][[deriv]]
      }
      mean.cov.wgt.mtx.B[[deriv]] <- mean.cov.wgt.mtx.B[[deriv]]/B
    }
  }
  mean.cov.wgt.mtx.B

  # add Sep-2017 result for cor est
  corr.est <- dynamicCorrelation(dataFrame = dataFrame,
                                 depVar = depVar,
                                 indepVar = indepVar,
                                 subjectVar = subjectVar,
                                 function.choice = function.choice,
                                 width.range = width.range,
                                 width.place = width.place,
                                 min.obs = min.obs,
                                 points.by = points.by,
                                 points.length = points.length,
                                 boundary.trunc = boundary.trunc,
                                 lag.input = max.dynCorrLag,
                                 byOrder = byOrder,
                                 by.deriv.only = by.deriv.only)

  if (by.deriv.only == FALSE) {
    dim <- num_depVar * num_funcVar
    corr <- vector(mode = "list", dim)
    for (i in 1:dim) {
      corr[[i]] <- vector(mode = "list", length(dim))
      for (j in 1:dim) {
        corr[[i]][[j]] <- c(cov.wgt.mtx.B[[1]][i, j])
        for (b in 2:B) {
          corr[[i]][[j]] <- c(corr[[i]][[j]], cov.wgt.mtx.B[[b]][i, j])
        }
      }
    }
    CI <- matrix(0, dim * (dim - 1)/2, 6)
    count <- 1
    namesRow <- c()
    for (i in 1:dim) {
      for (j in (i + 1):dim) {
        if (dim >= (i + 1)) {
          rowVar <- (i - 1)%/%num_funcVar + 1
          rowDer <- (i - 1)%%num_funcVar + 1
          colVar <- (j - 1)%/%num_funcVar + 1
          colDer <- (j - 1)%%num_funcVar + 1
          temp <- quantile(corr[[i]][[j]], percentile)
          CI[count, 4] <- round(temp[[1]], 4)
          CI[count, 5] <- round(temp[[2]], 4)

          # add Sep-2017 cor.est
          CI[count, 3] <- round(corr.est$dynCorrMatrix[i,j], 4)
          CI[count, 2] <- max.dynCorrLag
          CI[count, 1] <- " "
          if (mean.cov.wgt.mtx.B[i, j] > 0) {
            CI[count, 6] <- round(2 * sum(corr[[i]][[j]] <= 0)/B, 5)
          }
          else {
            CI[count, 6] <- round(2 * sum(corr[[i]][[j]] >= 0)/B, 5)
          }
          count <- count + 1
          name <- paste(depVar[rowVar], funcVar[rowDer] - 1,
                        "x", depVar[colVar], funcVar[colDer] - 1)
          namesRow <- c(namesRow, name)
        }
      }
    }
    perc1 <- paste(percentile[1] * 100, "%")
    perc2 <- paste(percentile[2] * 100, "%")
    namesCol <- c(" ", "lag", "dCorr.est", perc1, perc2, "BS p-value")
    dimnames(CI) <- list(namesRow, namesCol)
    CI <- as.data.frame(CI)
    return(list(quantilesMatrix = CI,
                # added Nov-2017 include data summary
                data.summary = data_summary))
  }
  else {
    dim <- num_depVar
    corr <- vector(mode = "list", num_funcVar)

    for (deriv in 1:num_funcVar) {
      corr[[deriv]] <- vector(mode = "list", dim)
      for (i in 1:dim) {
        corr[[deriv]][[i]] <- vector(mode = "list", length(dim))
        for (j in 1:dim) {
          corr[[deriv]][[i]][[j]] <- c(cov.wgt.mtx.B[[1]][[deriv]][i, j])
          for (b in 2:B) {
            corr[[deriv]][[i]][[j]] <- c(corr[[deriv]][[i]][[j]],
                                         cov.wgt.mtx.B[[b]][[deriv]][i, j])
          }
        }
      }
    }
    CI <- vector(mode = "list", num_funcVar)
    perc1 <- paste(percentile[1] * 100, "%")
    perc2 <- paste(percentile[2] * 100, "%")
    namesCol <- c(" ", "lag", "dCorr.est", perc1, perc2, "BS p-value")
    for (deriv in 1:num_funcVar) {
      CI[[deriv]] <- matrix(0, dim * (dim - 1)/2, 6)
      count <- 1
      namesRow <- c()
      for (i in 1:dim) {
        for (j in (i + 1):dim) {
          if (dim >= (i + 1)) {
            temp <- quantile(corr[[deriv]][[i]][[j]],
                             percentile)
            CI[[deriv]][count, 4] <- round(temp[[1]], 4)
            CI[[deriv]][count, 5] <- round(temp[[2]], 4)
            # add dyn Corr result - Sep 2017
            CI[[deriv]][count, 3] <- round(corr.est$dynCorrMatrix[[deriv]][i, j],4)
            CI[[deriv]][count, 2] <- max.dynCorrLag
            CI[[deriv]][count, 1] <- " "
            if (mean.cov.wgt.mtx.B[[deriv]][i, j] > 0) {
              CI[[deriv]][count, 6] <-
                round(2 * sum(corr[[deriv]][[i]][[j]] <= 0)/B, 5)
            }
            else {
              CI[[deriv]][count, 6] <-
                round(2 * sum(corr[[deriv]][[i]][[j]] >= 0)/B, 5)
            }
            count <- count + 1
            name <- paste(depVar[i], funcVar[deriv] - 1, "x",
                          depVar[j], funcVar[deriv] - 1)
            namesRow <- c(namesRow, name)
          }
        }
      }
      dimnames(CI[[deriv]]) <- list(namesRow, namesCol)
      CI[[deriv]] <- as.data.frame(CI[[deriv]])
    }
    return(list(quantilesMatrix = CI,
                # added Nov-2017 include data summary
                data.summary = data_summary))
  }
  }
