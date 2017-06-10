#' Dynamic Correlation
#'
#' Computes dynamical correlation estimates for pairs of longitudinal
#' responses, including consideration of lags and derivatives,
#' following a local polynomial regression smoothing step.
#'
#' @usage
#' dynamicCorrelation(dataFrame, depVar, indepVar, subjectVar,
#'                    function.choice, width.range, width.place,
#'                    boundary.trunc, lag.input, byOrder,
#'                    by.deriv.only, full.lag.output)
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
#' @param boundary.trunc Indicate the boundary of indepVar that should be
#'  truncated after smoothing; this may be done in case of concerns about
#'  estimating dynamical correlation at the boundaries of indepVar; this
#'  is a two element vector, where the first argument is how much to
#'  truncate from the right of the minimum value of indep.var and the
#'  second argument is how much to truncate from the left of the maximum
#'  value of indep var, within an individual and specific response;
#'  default is no truncation, i.e., c(0,0).
#'
#' @param lag.input Values of lag to be considered; can be a vector of
#' requested lags, for which a dynamical correlation estimate will be
#' produced for each pair of responses at each lag; a positive value
#' of lag.input means that the first entry for the dynamical correlation
#' leads (occurs before) the second entry — conversely, a negative value
#' means that the second entry for the correlation leads the first entry;
#' default is no lag at all considered
#'
#' @param byOrder A vector that specifies the order of the
#'  variables/responses and derivatives (if any) to be the leading
#'  variable in the calculations; this will have an effect on how
#'  lag.input is to be interpreted; default is to use the order as
#'  specified in the depVar argument.
#'
#' @param by.deriv.only If TRUE, the inter-dynamical correlations
#' between different derivatives are not computed, which can save
#' computation time (e.g., when function.choice=c(1,0,1) is specified
#' and by.deriv.only=T, then only dynamical correlations will be
#' calculated within (and not across) the 0th and 2nd derivative
#' estimates, respectively); default is TRUE
#'
#' @param full.lag.output If TRUE, the dynamical correlation values
#' for each pair of responses and requested derivative will be stored
#' in vectors corresponding to different lag values, which enables
#' plotting the correlations as a function of lag values; all the vectors
#' will be stored in the returned attribute resultMatrix; default is FALSE
#'
#' @details This function will provide smooth estimates (curves/functions
#'  or their derivatives) of longitudinal responses of interest per
#'  individual, then generate dynamical correlation estimates for each pair
#'  of responses. Lags of interest can be specified, using the lag.input
#'  argument. For smoothing, the function uses local polynomial smoothing
#'  using Epanecknikov kernel weighting (by calling the function lpepa within
#'  the lpridge package). The default global bandwidth is generated by taking
#'  the range of indepVar (usually time) and dividing by 4. This, by default,
#'  will be a constant global bandwidth, but proper specification of the
#'  width.range and width.place arguments can allow for a more flexible
#'  bandwidth choice, including different specification for each response
#'  in depVar. \cr \cr
#'  Details of the methodology for dynamical correlation can be found in
#'  Dubin and Muller (2005).
#'
#' @author Joel A. Dubin, Mike Li, Dandi Qiao, Hans-Georg Müller
#' @seealso \link[lpridge]{lpepa}
#' @examples
#'
#' ## Example 1: using default smoothing parameters, obtain dynamical
#' ##            correlation estimates for all three pairs of responses,
#' ##            for both original function and the first derivative
#'
#' examp1 <- dynamicCorrelation(dataFrame=dynCorrData,
#'                              depVar=c('resp1', 'resp2', 'resp3'),
#'                              indepVar='time',
#'                              subjectVar = 'subject',
#'                              function.choice = c(1,1,0))
#' examp1
#'
#' ## Example 2: using default smoothing parameters, obtain dynamical
#' ##            correlation estimates for all three pairs of responses,
#' ##            looking at range of lags between -10 and +10, for original
#' ##            functions only
#'
#' examp2 <- dynamicCorrelation(dataFrame=dynCorrData,
#'                              depVar=c('resp1', 'resp2', 'resp3'),
#'                              indepVar='time',
#'                              subjectVar = 'subject',
#'                              function.choice = c(1,0,0),
#'                              lag.input=seq(-20,20, by=1))
#' examp2
#'
#' ## note: output includes zero lag correlations, as well as maximum
#' ##       correlation (in absolute value) in max.dynCorr and and its
#' ##       corresponding lag value in max.dynCorrLag
#'
#' ## Example 3: re-rerun example 2, but set up for plotting of specified
#' ##            lagged correlations
#'
#' examp3 <- dynamicCorrelation(dataFrame=dynCorrData,
#'                              depVar=c('resp1', 'resp2', 'resp3'),
#'                              indepVar='time',
#'                              subjectVar = 'subject',
#'                              function.choice = c(1,0,0),
#'                              lag.input=seq(-20,10, by=1),
#'                              full.lag.output=TRUE)
#'
#' # conduct plotting, with one panel for each pair of responses considered;
#' # the ylim adjustment is made here for the different magnitude of the
#' # correlations between the two pairs
#'
#' par(mfrow=c(1,2))
#' plot(seq(-20,10, by=1),
#'      examp3$lagResultMatrix[[1]][1,],
#'      type='b',
#'      xlab = 'lag order (in days)',
#'      ylab = 'lagged correlations',
#'      ylim = c(-0.4, -0.2),
#'      main = 'dyncorr b/t resp1 and resp2 as a function of lag')
#'
#' abline(v = examp3$max.dynCorrLag[[1]][1,2], lty = 2)
#'
#' plot(seq(-20,10, by=1),
#'      examp3$lagResultMatrix[[1]][2,],
#'      type='b',
#'      xlab = 'lag order (in days)',
#'      ylab = 'lagged correlations',
#'      ylim = c(0.3, 0.5),
#'      main = 'dyncorr b/t resp1 and resp3 as a function of lag')
#'
#' abline(v = examp3$max.dynCorrLag[[1]][1,3], lty = 2)
#'
#' ## Example 4: same as the original function piece of Example 1,
#' ##            except now adjust the constant global bandwidth
#' ##            from the default to 40
#'
#' examp4 <- dynamicCorrelation(dataFrame=dynCorrData,
#'                              depVar=c('resp1', 'resp2', 'resp3'),
#'                              indepVar='time',
#'                              subjectVar = 'subject',
#'                              function.choice = c(1,0,0),
#'                              width.range = c(40, 40))
#' examp4
#'
#' @export
#'
#' @importFrom lpridge lpepa
#' @importFrom stats var

# ---------------------------------------------
dynamicCorrelation <- function (dataFrame,
                                depVar = c("resp1", "resp2", "resp3", "resp4", "resp5"),
                                indepVar = "time",
                                subjectVar = "subject",
                                function.choice = c(1),
                                width.range = c(((range(dataFrame[[indepVar]]))[2] -
                                                   (range(dataFrame[[indepVar]]))[1])/4,
                                                ((range(dataFrame[[indepVar]]))[2] -
                                                   (range(dataFrame[[indepVar]]))[1])/4),
                                width.place = c(NA, NA),
                                boundary.trunc = c(0, 0),
                                lag.input = c(),
                                byOrder = c(),
                                by.deriv.only = TRUE,
                                full.lag.output = FALSE)
{

  ## --------------------Set up -------------------

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
  for (i in 1:length(vec)) {
    v_ob_time <- c(v_ob_time, max(indep_var[subject_var == vec[i]]))
  }
  limit <- min(v_ob_time) # max common obs



  ## ------------ Smooth Curves -----------------------

  # create a list to store curves
  smoothedCurves <- vector(mode = "list", length(vec))

  for (k in 1:length(vec)) {
    indep <- indep_var[subject_var == vec[k]]
    smoothedCurves[[k]] <- vector(mode = "list", num_depVar)
    max_indep <- max(indep[!is.na(indep)])

    for (i in 1:num_depVar) {
      # Set width.range and width.place for current depVar
      cur_wrange <- width.range[[i]]
      cur_wplace <- width.range[[i]]

      if (is.na(width.place[[i]][1])) {
        cur_wplace[1] <- min(indep[!is.na(indep)])
      }
      if (is.na(width.place[[i]][2])) {
        cur_wplace[2] <- max_indep
      }

      points.use = seq(l_trunc, ceiling(max_indep - h_trunc), by = 1)
      size <- ceiling(max_indep - h_trunc) - l_trunc + 1
      band <- c()

      # calculate band width
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

      smoothedCurves[[k]][[i]] <- vector(mode = "list", num_funcVar)

      curve_x <- indep_var[subject_var == vec[k]][!is.na(
        dep_var[[i]][subject_var == vec[k]])]
      curve_y <- dep_var[[i]][subject_var == vec[k]][!is.na(
        dep_var[[i]][subject_var == vec[k]])]

      # produce curves
      for (j in 1:num_funcVar) {
        der <- funcVar[j] - 1
        temp <- lpepa(x = curve_x, y = curve_y, bandwidth = band, deriv = der,
                      n.out = size, x.out = points.use, order = funcVar[j],
                      var = FALSE)$est
        smoothedCurves[[k]][[i]][[j]] <- temp
      }
    }
  }

  max_len <- (ceiling(limit) - h_trunc) - l_trunc + 1
  meanMatrix <- vector(mode = "list", num_depVar)

  # calculate mean matrix
  for (i in 1:num_depVar) {
    meanMatrix[[i]] <- vector(mode = "list", num_funcVar)
    for (j in 1:num_funcVar) {
      forEachDerivMean <- rep(NA, max_len)
      for (m in 1:max_len) {
        for (n in 1:length(vec)) {
          if (n == 1) {
            total <- 1
            forEachDerivMean[m] <- smoothedCurves[[n]][[i]][[j]][m]
          }
          else if (n > 1) {
            total <- total + 1
            forEachDerivMean[m] <-
              (1/total) * ((total - 1) * forEachDerivMean[m]
                           + smoothedCurves[[n]][[i]][[j]][m])
          }
        }
      }
      meanMatrix[[i]][[j]] <- forEachDerivMean
    }
  }

  # correct terms
  correst_all_stand <- vector(mode = "list", length(vec))
  for (k in 1:length(vec)) {
    correst_all_stand[[k]] <- list()
    correst_all_stand[[k]][[1]] <- list(l_trunc:(l_trunc + max_len - 1))
    for (i in 2:(num_depVar + 1)) {
      correst_all_stand[[k]][[i]] <- vector(mode = "list", num_funcVar)
      for (j in 1:num_funcVar) {
        temp <- smoothedCurves[[k]][[i-1]][[j]][1:max_len]
        correct <- temp - meanMatrix[[i-1]][[j]]
        mean_correct <- mean(correct)
        sd_correct <- sqrt(var(correct))
        correst_all_stand[[k]][[i]][[j]] <-
          list((correct - mean_correct) /
                 (sqrt((max_len - 1)/max_len) * sd_correct))
      }
    }
  }
  weights_vec <- rep(NA, length(vec))
  time_extend <- (l_trunc + max_len - 1) + l_trunc

  if (by.deriv.only == TRUE) {
    cov.mtx.listz <- vector(mode = "list", num_funcVar)
    cov.wgt.mtxz <- vector(mode = "list", num_funcVar)
    dim <- num_depVar
    for (deriv in 1:num_funcVar) {
      cov.mtx.listz[[deriv]] <- vector(mode = "list", length(vec))
      for (i in 1:length(vec)) {
        weights_vec[i] <- length(subject_var[(subject_var == vec[i]) &
                                               (indep_var >= 0) &
                                               (indep_var <= time_extend)])
        cov.mtx.listz[[deriv]][[i]] <- matrix(nrow = dim, ncol = dim)
        m <- max_len
        diag(cov.mtx.listz[[deriv]][[i]]) <- 1
        for (j in 1:dim) {
          for (k in (j + 1):dim) {
            if (k <= dim) {
              cov.mtx.listz[[deriv]][[i]][j, k] <-
                cov.mtx.listz[[deriv]][[i]][k, j] <-
                ((1/m) *
                   sum(correst_all_stand[[i]][[(j + 1)]][[deriv]][[1]] *
                         correst_all_stand[[i]][[(k +  1)]][[deriv]][[1]]))
            }
          }
        }
      }
      cov.wgt.mtxz[[deriv]] <- matrix(0, dim, dim)
      for (i in 1:length(vec)) {
        cov.wgt.mtxz[[deriv]] <- cov.wgt.mtxz[[deriv]] +
          (weights_vec[i] * cov.mtx.listz[[deriv]][[i]])
      }
      cov.wgt.mtxz[[deriv]] <- (1/sum(weights_vec)) * cov.wgt.mtxz[[deriv]]
      names <- c()
      for (i in 1:length(dep_var)) {
        name <- depVar[i]
        name2 <- paste(name, (deriv - 1))
        names <- c(names, name2)
      }
      dimnames(cov.wgt.mtxz[[deriv]]) <- list(names, names)
    }
  }
  else {
    cov.mtx.listz <- list()
    dim <- length(dep_var) * num_funcVar
    for (i in 1:length(vec)) {
      weights_vec[i] <- length(subject_var[(subject_var == vec[i]) &
                                             (indep_var >= 0) &
                                             (indep_var <= time_extend)])
      cov.mtx.listz[[i]] <- matrix(nrow = dim, ncol = dim)
      m <- max_len
      diag(cov.mtx.listz[[i]]) <- 1
      for (j in 1:dim) {
        for (k in (j + 1):dim) {
          if (k <= dim) {
            m1 <- (j - 1)%/%num_funcVar + 2
            m2 <- (j - 1)%%num_funcVar + 1
            s1 <- (k - 1)%/%num_funcVar + 2
            s2 <- (k - 1)%%num_funcVar + 1
            cov.mtx.listz[[i]][j, k] <-
              cov.mtx.listz[[i]][k,j] <-
              ((1/m) * sum(correst_all_stand[[i]][[m1]][[m2]][[1]] *
                             correst_all_stand[[i]][[s1]][[s2]][[1]]))
          }
        }
      }
    }
    cov.wgt.mtxz <- matrix(0, dim, dim)
    for (i in 1:length(vec)) {
      cov.wgt.mtxz <- cov.wgt.mtxz + (weights_vec[i] *
                                        cov.mtx.listz[[i]])
    }
    cov.wgt.mtxz <- (1/sum(weights_vec)) * cov.wgt.mtxz
    names <- c()
    for (i in 1:num_depVar) {
      name <- depVar[i]
      for (j in 1:num_funcVar) {
        name2 <- paste(name, (funcVar[j] - 1))
        names <- c(names, name2)
      }
    }
    dimnames(cov.wgt.mtxz) <- list(names, names)
  }

  cov.wgt.mtxz

  if (length(lag.input) != 0) {
    if (by.deriv.only == FALSE) {
      dim <- num_depVar * num_funcVar
      cov.lag.mtx.listz <- list()
      for (i in 1:length(vec)) {
        cov.lag.mtx.listz[[i]] <- list()
        dep.correct <- list()
        for (m in 1:num_depVar) {
          dep.correct.function <- list()
          for (n in 1:num_funcVar) {
            dep.correct.function[[n]] <-
              smoothedCurves[[i]][[m]][[n]][1:max_len] - meanMatrix[[m]][[n]]
          }
          dep.correct[[m]] <- dep.correct.function
        }
        for (j in 1:length(lag.input)) {
          cov.lag.mtx.listz[[i]][[j]] <- matrix(nrow = dim, ncol = dim)
          diag(cov.lag.mtx.listz[[i]][[j]]) <- 1

          if (lag.input[j] >= 0) {
            lag_end <- max_len - lag.input[j]
            lag.beg <- 1 + lag.input[j]
            m.support <- lag_end #

            correst.standz <- list()
            if (length(byOrder) == 0) index0 <- 1
            else index0 <- byOrder[1]

            dep <- (index0 - 1)%/%num_funcVar + 1
            deriv <- (index0 - 1)%%num_funcVar + 1
            mean <- mean(dep.correct[[dep]][[deriv]][1:lag_end])
            sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag_end]))
            correst.standz[[index0]] <-
              (dep.correct[[dep]][[deriv]][1:lag_end] - mean)/
              (sqrt((m.support - 1)/m.support) * sd)

            for (ord in 2:dim) {
              if (dim >= 2) {
                index <- byOrder[ord]
                dep <- (index - 1)%/%num_funcVar + 1
                deriv <- (index - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
                correst.standz[[index]] <-
                  (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)
                cov.lag.mtx.listz[[i]][[j]][index0, index] <-
                  cov.lag.mtx.listz[[i]][[j]][index, index0] <-
                  (1/m.support) *
                  sum(correst.standz[[index0]] * correst.standz[[index]])
              }
            }

            for (ord in 2:(dim - 1)) {
              if (dim >= 3) {
                if (length(byOrder) == 0) d <- ord
                else d <- byOrder[ord]
                dep <- (d - 1)%/%num_funcVar + 1
                deriv <- (d - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][1:lag_end])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag_end]))
                correst.standz[[d]] <-
                  (dep.correct[[dep]][[deriv]][1:lag_end] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)
                for (ord2 in (ord + 1):dim) {
                  if (length(byOrder) == 0) e <- ord2
                  else e <- byOrder[ord2]
                  cov.lag.mtx.listz[[i]][[j]][d, e] <-
                    cov.lag.mtx.listz[[i]][[j]][e, d] <-
                    (1/m.support) *
                    sum(correst.standz[[d]] * correst.standz[[e]])
                }
              }
            }
          }

          else if (lag.input[j] < 0) {
            lag_end <- max_len + lag.input[j]
            lag.beg <- 1 - lag.input[j]
            m.support <- lag_end

            correst.standz <- list()
            if (length(byOrder) == 0) index0 <- 1
            else index0 <- byOrder[1]
            dep <- (index0 - 1)%/%num_funcVar + 1
            deriv <- (index0 - 1)%%num_funcVar + 1
            mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
            sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
            correst.standz[[index0]] <-
              (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean) /
              (sqrt((m.support - 1)/m.support) * sd)

            for (ord in 2:dim) {
              if (length(byOrder) == 0) index <- ord
              else index <- byOrder[ord]

              if (dim >= 2) {
                dep <- (index - 1)%/%num_funcVar + 1
                deriv <- (index - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][1:lag_end])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag_end]))
                correst.standz[[index]] <-
                  (dep.correct[[dep]][[deriv]][1:lag_end] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)

                # code fixed prior to version 1.1. change from 1 to index0
                cov.lag.mtx.listz[[i]][[j]][index0, index] <-
                  cov.lag.mtx.listz[[i]][[j]][index, index0] <-
                  (1/m.support) *
                  sum(correst.standz[[index0]] * correst.standz[[index]])
              }
            }

            for (ord in 2:(dim - 1)) {
              if (dim >= 3) {
                if (length(byOrder) == 0) d <- ord
                else d <- byOrder[ord]

                dep <- (d - 1)%/%num_funcVar + 1
                deriv <- (d - 1)%%num_funcVar + 1
                mean <- mean(dep.correct[[dep]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max_len]))
                correst.standz[[d]] <-
                  (dep.correct[[dep]][[deriv]][lag.beg:max_len] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)

                for (ord2 in (ord + 1):dim) {
                  if (length(byOrder) == 0) e <- ord2
                  else e <- byOrder[ord2]

                  cov.lag.mtx.listz[[i]][[j]][d, e] <-
                    cov.lag.mtx.listz[[i]][[j]][e, d] <-
                    (1/m.support) *
                    sum(correst.standz[[d]] * correst.standz[[e]])
                }
              }
            }
          }
        }
      }

      cov.lag.wgt.mtx.listz <- list()

      for (j in 1:length(lag.input)) {
        cov.lag.wgt.mtx.listz[[j]] <- matrix(0, nrow = dim, ncol = dim)
        for (i in 1:length(vec)) {
          cov.lag.wgt.mtx.listz[[j]] <- cov.lag.wgt.mtx.listz[[j]] +
            (weights_vec[i] * cov.lag.mtx.listz[[i]][[j]])
        }
        cov.lag.wgt.mtx.listz[[j]] <- (1/sum(weights_vec)) *
          cov.lag.wgt.mtx.listz[[j]]
      }

      cov.lag.wgt.mtxz <- matrix(0, nrow = dim, ncol = dim)
      lag.max.cov.mtxz <- matrix(NA, nrow = dim, ncol = dim)

      for (j in 1:length(lag.input)) {
        for (k in 1:(dim - 1)) {
          for (l in (k + 1):dim) {
            if (abs(cov.lag.wgt.mtx.listz[[j]][k, l]) >
                abs(cov.lag.wgt.mtxz[k, l])) {
              cov.lag.wgt.mtxz[k, l] <-
                cov.lag.wgt.mtxz[l, k] <- cov.lag.wgt.mtx.listz[[j]][k, l]
              lag.max.cov.mtxz[k, l] <-
                lag.max.cov.mtxz[l, k] <- lag.input[j]
            }
          }
        }
      }
      diag(cov.lag.wgt.mtxz) <- 1
      names <- c()
      for (i in 1:length(dep_var)) {
        name <- depVar[i]
        for (j in 1:num_funcVar) {
          name2 <- paste(name, (funcVar[j] - 1))
          names <- c(names, name2)
        }
      }
      dimnames(cov.lag.wgt.mtxz) <-
        dimnames(lag.max.cov.mtxz) <- list(names, names)
      round(cov.lag.wgt.mtxz, 3)
      lag.max.cov.mtxz
    }
    else {
      cov.lag.mtx.listz <- vector(mode = "list", length(vec))
      dim <- num_depVar
      for (i in 1:length(vec)) {
        cov.lag.mtx.listz[[i]] <- vector(mode = "list", num_funcVar)
        dep.correct <- list()
        for (m in 1:num_depVar) {
          dep.correct.function <- list()
          for (n in 1:num_funcVar) {
            dep.correct.function[[n]] <-
              smoothedCurves[[i]][[m]][[n]][1:max_len] - meanMatrix[[m]][[n]]
          }
          dep.correct[[m]] <- dep.correct.function
        }
        for (deriv in 1:num_funcVar) {
          cov.lag.mtx.listz[[i]][[deriv]] <-
            vector(mode = "list", length(lag.input))
          for (j in 1:length(lag.input)) {
            cov.lag.mtx.listz[[i]][[deriv]][[j]] <-
              matrix(nrow = dim, ncol = dim)
            diag(cov.lag.mtx.listz[[i]][[deriv]][[j]]) <- 1

            if (lag.input[j] >= 0) {
              lag_end <- max_len - lag.input[j]
              lag.beg <- 1 + lag.input[j]
              m.support <- lag_end
              if (length(byOrder) == 0) {
                correst.standz <- list()
                index <- 1
                mean <- mean(dep.correct[[index]][[deriv]][1:lag_end])
                sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag_end]))
                correst.standz[[index]] <-
                  (dep.correct[[index]][[deriv]][1:lag_end] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)
                for (index in 2:dim) {
                  if (dim >= 2) {
                    mean <-
                      mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
                    sd <-
                      sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
                    correst.standz[[index]] <-
                      (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    cov.lag.mtx.listz[[i]][[deriv]][[j]][1, index] <-
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][index, 1] <-
                      (1/m.support) *
                      sum(correst.standz[[1]] * correst.standz[[index]])
                  }
                }
                for (d in 2:(dim - 1)) {
                  if (dim >= 3) {
                    mean <- mean(dep.correct[[d]][[deriv]][1:lag_end])
                    sd <- sqrt(var(dep.correct[[d]][[deriv]][1:lag_end]))
                    correst.standz[[d]] <-
                      (dep.correct[[d]][[deriv]][1:lag_end] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    for (e in (d + 1):dim) {
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][d, e] <-
                        cov.lag.mtx.listz[[i]][[deriv]][[j]][e, d] <-
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
                mean <- mean(dep.correct[[index]][[deriv]][1:lag_end])
                sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag_end]))
                correst.standz[[index]] <-
                  (dep.correct[[index]][[deriv]][1:lag_end] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)
                for (ord in 2:dim) {
                  if (dim >= 2) {
                    index <- byOrderPerD[ord]
                    mean <-
                      mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
                    sd <-
                      sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
                    correst.standz[[index]] <-
                      (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    cov.lag.mtx.listz[[i]][[deriv]][[j]][byOrderPerD[1], index] <-
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][index, byOrderPerD[1]] <-
                      (1/m.support) *
                      sum(correst.standz[[byOrderPerD[1]]] *correst.standz[[index]])
                  }
                }
                for (ord in 2:(dim - 1)) {
                  if (dim >= 3) {
                    d <- byOrderPerD[ord]
                    mean <- mean(dep.correct[[d]][[deriv]][1:lag_end])
                    sd <- sqrt(var(dep.correct[[d]][[deriv]][1:lag_end]))
                    correst.standz[[d]] <-
                      (dep.correct[[d]][[deriv]][1:lag_end] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    for (ord2 in (ord + 1):dim) {
                      e <- byOrderPerD[ord2]
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][d, e] <-
                        cov.lag.mtx.listz[[i]][[deriv]][[j]][e, d] <-
                        (1/m.support) *
                        sum(correst.standz[[d]] * correst.standz[[e]])
                    }
                  }
                }
              }
            }
            else if (lag.input[j] < 0) {
              lag_end <- max_len + lag.input[j]
              lag.beg <- 1 - lag.input[j]
              m.support <- lag_end
              if (length(byOrder) == 0) {
                correst.standz <- list()
                index <- 1
                mean <- mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
                correst.standz[[index]] <-
                  (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)
                for (index in 2:dim) {
                  if (dim >= 2) {
                    mean <- mean(dep.correct[[index]][[deriv]][1:lag_end])
                    sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag_end]))
                    correst.standz[[index]] <-
                      (dep.correct[[index]][[deriv]][1:lag_end] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    cov.lag.mtx.listz[[i]][[deriv]][[j]][1, index] <-
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][index, 1] <-
                      (1/m.support) *
                      sum(correst.standz[[1]] * correst.standz[[index]])
                  }
                }
                for (d in 2:(dim - 1)) {
                  if (dim >= 3) {
                    mean <- mean(dep.correct[[d]][[deriv]][lag.beg:max_len])
                    sd <- sqrt(var(dep.correct[[d]][[deriv]][lag.beg:max_len]))
                    correst.standz[[d]] <-
                      (dep.correct[[d]][[deriv]][lag.beg:max_len] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    for (e in (d + 1):dim) {
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][d, e] <-
                        cov.lag.mtx.listz[[i]][[deriv]][[j]][e, d] <-
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
                mean <- mean(dep.correct[[index]][[deriv]][lag.beg:max_len])
                sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max_len]))
                correst.standz[[index]] <-
                  (dep.correct[[index]][[deriv]][lag.beg:max_len] - mean) /
                  (sqrt((m.support - 1)/m.support) * sd)
                for (ord in 2:dim) {
                  index <- byOrderPerD[ord]
                  if (dim >= 2) {
                    mean <- mean(dep.correct[[index]][[deriv]][1:lag_end])
                    sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag_end]))
                    correst.standz[[index]] <-
                      (dep.correct[[index]][[deriv]][1:lag_end] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    cov.lag.mtx.listz[[i]][[deriv]][[j]][byOrderPerD[1], index] <-
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][index, byOrderPerD[1]] <-
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
                      (dep.correct[[d]][[deriv]][lag.beg:max_len] - mean) /
                      (sqrt((m.support - 1)/m.support) * sd)
                    for (ord2 in (ord + 1):dim) {
                      e <- byOrderPerD[ord2]
                      cov.lag.mtx.listz[[i]][[deriv]][[j]][d, e] <-
                        cov.lag.mtx.listz[[i]][[deriv]][[j]][e, d] <-
                        (1/m.support) *
                        sum(correst.standz[[d]] * correst.standz[[e]])
                    }
                  }
                }
              }
            }
          }
        }
      }
      cov.lag.wgt.mtx.listz <- vector(mode = "list", num_funcVar)
      for (deriv in 1:num_funcVar) {
        cov.lag.wgt.mtx.listz[[deriv]] <-
          vector(mode = "list", length(lag.input))
        for (j in 1:length(lag.input)) {
          cov.lag.wgt.mtx.listz[[deriv]][[j]] <-
            matrix(0, nrow = dim, ncol = dim)
          for (i in 1:length(vec)) {
            cov.lag.wgt.mtx.listz[[deriv]][[j]] <-
              cov.lag.wgt.mtx.listz[[deriv]][[j]] +
              (weights_vec[i] * cov.lag.mtx.listz[[i]][[deriv]][[j]])
          }
          cov.lag.wgt.mtx.listz[[deriv]][[j]] <-
            (1/sum(weights_vec)) *
            cov.lag.wgt.mtx.listz[[deriv]][[j]]
        }
      }
      cov.lag.wgt.mtxz <- vector(mode = "list", num_funcVar)
      lag.max.cov.mtxz <- vector(mode = "list", num_funcVar)
      for (deriv in 1:num_funcVar) {
        cov.lag.wgt.mtxz[[deriv]] <- matrix(0, nrow = dim, ncol = dim)
        lag.max.cov.mtxz[[deriv]] <- matrix(NA, nrow = dim, ncol = dim)
        for (j in 1:length(lag.input)) {
          for (k in 1:(dim - 1)) {
            for (l in (k + 1):dim) {
              if (abs(cov.lag.wgt.mtx.listz[[deriv]][[j]][k, l]) >
                  abs(cov.lag.wgt.mtxz[[deriv]][k, l])) {
                cov.lag.wgt.mtxz[[deriv]][k, l] <-
                  cov.lag.wgt.mtxz[[deriv]][l, k] <-
                  cov.lag.wgt.mtx.listz[[deriv]][[j]][k, l]
                lag.max.cov.mtxz[[deriv]][k, l] <-
                  lag.max.cov.mtxz[[deriv]][l, k] <- lag.input[j]
              }
            }
          }
        }
        diag(cov.lag.wgt.mtxz[[deriv]]) <- 1
        names <- c()
        for (i in 1:length(dep_var)) {
          name <- depVar[i]
          name2 <- paste(name, (funcVar[deriv] - 1))
          names <- c(names, name2)
        }
        dimnames(cov.lag.wgt.mtxz[[deriv]]) <-
          dimnames(lag.max.cov.mtxz[[deriv]]) <- list(names, names)
      }
      cov.lag.wgt.mtxz
      lag.max.cov.mtxz
    }
  }

  if (by.deriv.only == FALSE) {
    if (full.lag.output == FALSE) {
      if (length(lag.input) == 0) {
        return(list(dynCorrMatrix = cov.wgt.mtxz))
      }
      else {
        return(list(dynCorrMatrix = cov.wgt.mtxz,
                    lag.input = lag.input,
                    max.dynCorr = cov.lag.wgt.mtxz,
                    max.dynCorrLag = lag.max.cov.mtxz))
      }
    }
    else {
      if (length(lag.input) == 0) {
        return(list(dynCorrMatrix = cov.wgt.mtxz))
      }
      else {
        dim <- num_depVar * num_funcVar
        numRow <- dim * (dim - 1)/2
        result <- matrix(0, numRow, length(lag.input))
        namesRow <- c()
        namesCol <- c()
        for (m in 1:length(lag.input)) {
          lagName <- paste(lag.input[m])
          namesCol <- c(namesCol, lagName)
        }
        row <- 1
        for (i in 1:dim) {
          for (j in (i + 1):dim) {
            if (dim >= (i + 1)) {
              dep1 <- (i - 1)%/%num_funcVar + 1
              deriv1 <- (i - 1)%%num_funcVar + 1
              dep2 <- (j - 1)%/%num_funcVar + 1
              deriv2 <- (j - 1)%%num_funcVar + 1
              name <- depVar[dep1]
              name2 <- depVar[dep2]
              nameT <- paste(name, (funcVar[deriv1] - 1), " x ",
                             name2, (funcVar[deriv2] - 1))
              namesRow <- c(namesRow, nameT)
              for (m in 1:length(lag.input)) {
                result[row, m] <- cov.lag.wgt.mtx.listz[[m]][i, j]
              }
              row <- row + 1
            }
          }
        }
        dimnames(result) <- list(namesRow, namesCol)
        return(list(dynCorrMatrix = cov.wgt.mtxz,
                    lag.input = lag.input,
                    lagResultMatrix = result,
                    max.dynCorr = cov.lag.wgt.mtxz,
                    max.dynCorrLag = lag.max.cov.mtxz))
      }
    }
  }
  else {
    if (full.lag.output == FALSE) {
      if (length(lag.input) == 0) {
        return(list(dynCorrMatrix = cov.wgt.mtxz))
      }
      else {
        return(list(dynCorrMatrix = cov.wgt.mtxz,
                    lag.input = lag.input,
                    max.dynCorr = cov.lag.wgt.mtxz,
                    max.dynCorrLag = lag.max.cov.mtxz))
      }
    }
    else {
      if (length(lag.input) == 0) {
        return(list(dynCorrMatrix = cov.wgt.mtxz))
      }
      else {
        dim <- num_depVar
        numRows <- dim * (dim - 1)/2
        namesCol <- c()
        for (m in 1:length(lag.input)) {
          lagName <- paste(lag.input[m])
          namesCol <- c(namesCol, lagName)
        }
        result <- vector(mode = "list", num_funcVar)
        for (k in 1:num_funcVar) {
          result[[k]] <- matrix(0, numRows, length(lag.input))
          namesRow <- c()
          row <- 1
          for (i in 1:dim) {
            for (j in (i + 1):dim) {
              if (dim >= i + 1) {
                name <- paste(depVar[i], funcVar[k] - 1, "x",
                              depVar[j], funcVar[k] - 1)
                namesRow <- c(namesRow, name)
                for (m in 1:length(lag.input)) {
                  result[[k]][row, m] <- cov.lag.wgt.mtx.listz[[k]][[m]][i, j]
                }
                row <- row + 1
              }
            }
          }
          dimnames(result[[k]]) <- list(namesRow, namesCol)
        }
        return(list(dynCorrMatrix = cov.wgt.mtxz,
                    lag.input = lag.input,
                    lagResultMatrix = result,
                    max.dynCorr = cov.lag.wgt.mtxz,
                    max.dynCorrLag = lag.max.cov.mtxz))
      }
    }
  }
}
