# --------------------------------------------------------

# Functions computing the Sobol indices (package sensobol)

# --------------------------------------------------------

sobol_boot <- function(d, i, N, params, matrices, R, first, total, boot) {
  
  # Stopping rule to check concordance between estimators and sample matrix
  # -------------------------------------------------------------------
  
  ms <- "Revise the correspondence between the matrices and the estimators"
  
  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    
    if (!first == "saltelli" & !first == "jansen" |
        !total == "jansen" & !total == "sobol" & !total == "homma" &
        !total == "janon" & !total == "glen") {
      
      stop(ms)
      
    }
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    
    if (!first == "sobol"| !total == "saltelli") {
      
      stop(ms)
      
    }
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {
    
    if (!first == "azzini" | !total == "azzini" &
        !total == "jansen" & !total == "sobol" & !total == "homma" &
        !total == "janon" & !total == "glen" & !total == "saltelli") {
      
      if (!total == "azzini" | !first == "saltelli" & !first == "jansen" &
          !first == "azzini" & !first == "sobol") {
        
        stop(ms)
        
      }
      
    }
    
  }
  
  # -------------------------------------
  
  k <- length(params)
  
  if (boot == TRUE) {
    
    m <- d[i, ]
    
  } else if (boot == FALSE) {
    
    m <- d
    
  }

  # Define vectors based on sample design
  # ------------------------------------------------------------------
  
  if (isTRUE(all.equal(matrices, c("A", "B", "AB")))) {
    
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, -c(1, 2)]
    Y_AB1<-m[,3]
    Y_ABn<-m[,ncol(m)]
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_BA <- m[,-c(1, 2)]
    Y_BA1<-m[,3]
    Y_BAn<-m[,ncol(m)]
    
  } else if (isTRUE(all.equal(matrices, c("A", "B", "AB", "BA")))) {
    
    Y_A <- m[, 1]
    Y_B <- m[, 2]
    Y_AB <- m[, 3:(k + 2)]
    Y_BA <- m[, (k + 3):ncol(m)]
    Y_AB1<-m[,3]
    Y_ABn<-m[,(k+2)]
    Y_BA1<-m[,(k + 3)]
    Y_BAn<-m[,ncol(m)]
    
  } # A warning might be needed here
  if (isTRUE(all.equal(matrices, c("A", "B", "AB"))) |
      isTRUE(all.equal(matrices, c("A", "B", "BA")))) {
    
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
    
  }
  
  # Define first-order estimators
  # --------------------------------------------------------------------
  
  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (first == "saltelli" | first == "jansen" | first == "sobol") {
    
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
  }
  
  # ----------------------------------
  if (first == "sobol") {
    
    Vifull <- 1 / N * sum(Y_A * Y_BA1 - f0^2)
    
    Viminoneind <- 1 / N * sum(Y_A * Y_BAn - f0^2) 
    
  } else if (first == "saltelli") {
    
    Vifull <- 1 / N * sum(Y_B * (Y_AB1 - Y_A))
    
    Viminoneind <- 1 / N * sum(Y_B * (Y_ABn - Y_A)) 
    
  } else if (first == "jansen") {
    
    Vifull <- VY - 1 / (2 * N) * sum((Y_B - Y_AB1)^2)
    
    Viminoneind <- VY - 1 / (2 * N) * sum((Y_B - Y_ABn)^2) 
    
  } else if (first == "azzini") {
    
    VYfull <- sum((Y_A - Y_B)^2 + (Y_BA1 - Y_AB1)^2)
    Vifull <- (2 * sum((Y_BA1 - Y_B) * (Y_A - Y_AB1)))
    
    VYminoneind <- sum((Y_A - Y_B)^2 + (Y_BAn - Y_ABn)^2)
    Viminoneind <- (2 * sum((Y_BAn - Y_B) * (Y_A - Y_ABn))) 
    
  } else {
    stop("first should be sobol, saltelli, jansen or azzini")
  }
  
  if (first == "azzini") {
    
    Sifull <- Vifull[1:length(params)] / VYfull[1:length(params)]
    
    Siminoneind <- Viminoneind[1:length(params)] / VYminoneind[1:length(params)] 
    
  } else {
    
    Sifull <- Vifull / VY
    
    Siminoneind <- Viminoneind / VY
    
  }
  
  # Define total-order estimators
  # --------------------------------------------------------------------
  
  # Define variance for estimators with A, B, AB; or A, B, BA matrices
  if (total == "azzini" | total == "jansen" | total == "sobol" |
      total == "homma" | total == "janon" | total == "glen" | total == "saltelli") {
    
    f0 <- 1 / (2 * N) * sum(Y_A + Y_B)
    VY <- 1 / (2 * N - 1) * sum((Y_A - f0)^2 + (Y_B - f0)^2)
    
  }
  
  # ----------------------------------
  if (total == "jansen") {
    
    Tifull <- (1 / (2 * N) * sum((Y_A - Y_AB1)^2)) / VY
    
    Timinoneind <- (1 / (2 * N) * sum((Y_A - Y_ABn)^2)) / VY
    
  } else if (total == "sobol") {
    
    Tifull <- ((1 / N) * sum(Y_A * (Y_A - Y_AB1))) / VY
    
    Timinoneind <- ((1 / N) * sum(Y_A * (Y_A - Y_ABn))) / VY
    
  } else if (total == "homma") {
    
    Tifull <- (VY - (1 / N) * sum(Y_A * Y_AB1) + f0^2) / VY
    
    Timinoneind <- (VY - (1 / N) * sum(Y_A * Y_ABn) + f0^2) / VY
    
  } else if (total == "saltelli") {
    
    Tifull <- 1 - ((1 / N * sum(Y_B * Y_BA1 - f0^2)) / VY)
    
    Timinoneind <- 1 - ((1 / N * sum(Y_B * Y_BAn - f0^2)) / VY)
    
  } else if (total == "janon") {
    
    Tifull <- 1 - (1 / N * sum(Y_A * Y_AB1) -
                     (1/ N * sum((Y_A + Y_AB1) / 2))^2) /
      (1 / N * sum((Y_A ^ 2 + Y_AB1^2) / 2) -
         (1/ N * sum((Y_A + Y_AB1) / 2))^2)
    
    Timinoneind <- 1 - (1 / N * sum(Y_A * Y_ABn) -
                          (1/ N * sum((Y_A + Y_ABn) / 2))^2) /
      (1 / N * sum((Y_A ^ 2 + Y_ABn^2) / 2) -
         (1/ N * sum((Y_A + Y_ABn) / 2))^2)
    
  } else if (total == "glen") {
    
    Tifull <- 1 - (1 / (N - 1) *
                     sum(((Y_A - mean(Y_A)) * (Y_AB1 - mean(Y_AB1))) /
                                      sqrt(stats::var(Y_A) * stats::var(Y_AB1))))
    
    Timinoneind <- 1 - (1 / (N - 1) *
                          sum(((Y_A - mean(Y_A)) * (Y_ABn - mean(Y_ABn))) /
                                sqrt(stats::var(Y_A) * stats::var(Y_ABn))))
    
  } else if (total == "azzini") {
    
    Tifull <- sum((Y_B - Y_BA1)^2 + (Y_A - Y_AB1)^2) /
      sum((Y_A - Y_B)^2 + (Y_BA1 - Y_AB1)^2)
    
    Timinoneind <- sum((Y_B - Y_BAn)^2 + (Y_A - Y_ABn)^2) /
      sum((Y_A - Y_B)^2 + (Y_BAn - Y_ABn)^2)
    
  } else {
    
    stop("total should be jansen, sobol, homma saltelli, janon, glen or azzini")
    
  }

  return(c(Sifull,Siminoneind, Tifull, Timinoneind))
  
}


bootstats <- function(b, conf = conf, type = type) {
  p <- length(b$t0)
  lab <- c("original", "bias", "std.error", "low.ci", "high.ci")
  tmp <- as.data.frame(matrix(nrow = p,
                              ncol = length(lab),
                              dimnames = list(NULL, lab)))
  for (i in 1:p) {
    # original estimation, bias, standard deviation
    tmp[i, "original"] <- b$t0[i]
    tmp[i, "bias"] <- mean(b$t[, i]) - b$t0[i]
    tmp[i, "std.error"] <- stats::sd(b$t[, i])
    # confidence interval
    
    if (type == "norm") {
      ci <- boot::boot.ci(b, index = i, type = "norm", conf = conf)
      
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$norm[2]
        tmp[i, "high.ci"] <- ci$norm[3]
      }
      
    } else if (type == "basic") {
      ci <- boot::boot.ci(b, index = i, type = "basic", conf = conf)
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$basic[4]
        tmp[i, "high.ci"] <- ci$basic[5]
      }
      
    } else if (type == "percent") {
      ci <- boot::boot.ci(b, index = i, type = "perc", conf = conf)
      
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$percent[4]
        tmp[i, "high.ci"] <- ci$percent[5]
      }
      
    } else if (type == "bca") {
      ci <- boot::boot.ci(b, index = i, conf = conf)
      
      if (!is.null(ci)) {
        tmp[i, "low.ci"] <- ci$bca[4]
        tmp[i, "high.ci"] <- ci$bca[5]
      }
    }
  }
  return(tmp)
}


#' Computation of Sobol' indices
#'
#' It allows to compute Sobol' indices up to the fourth-order using state-of-the-art estimators.
#'
#'@param matrices Character vector with the required matrices. The default is \code{matrices = c("A", "B", "AB")}.
#' See \code{\link{sobol_matrices}}.
#' @param Y  numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param params Character vector with the name of the model inputs.
#' @param first Estimator to compute first-order indices. Options are:
#' * \code{first = "saltelli"} \insertCite{Saltelli2010a}{sensobol}.
#' * \code{first = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{first = "sobol"}  \insertCite{Sobol1993}{sensobol}.
#' * \code{first = "azzini"} \insertCite{Azzini2020}{sensobol}.
#' @param total Estimator to compute total-order indices. Options are:
#' * \code{total = "jansen"} \insertCite{Jansen1999}{sensobol}.
#' * \code{total = "sobol"} \insertCite{Sobol2001}{sensobol}.
#' * \code{total = "homma"} \insertCite{Homma1996}{sensobol}.
#' * \code{total = "janon"} \insertCite{Janon2014}{sensobol}.
#' * \code{total = "glen"} \insertCite{Glen2012}{sensobol}.
#' * \code{total = "azzini"} \insertCite{Azzini2020}{sensobol}.
#' * \code{total = "saltelli"} \insertCite{Saltelli2008}{sensobol}.
#' @param order Whether to compute "first", "second", "third" or fourth-order Sobol' indices. Default
#' is \code{order = "first"}.
#' @param boot Logical. If TRUE, the function bootstraps the Sobol' indices. If FALSE, it provides point
#' estimates. Default is \code{boot = FALSE}.
#' @param R Positive integer, number of bootstrap replicas. Default is NULL.
#' @param parallel The type of parallel operation to be used (if any).
#' If missing, the default is taken from the option "boot.parallel"
#' (and if that is not set, "no"). For more information, check the
#' \code{parallel} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param ncpus Positive integer: number of processes to be used in parallel operation:
#' typically one would chose this to the number of available CPUs.
#' Check the \code{ncpus} option in the \code{boot} function of the \code{\link{boot}} package.
#' @param conf Confidence interval if \code{boot = TRUE}. Number between 0 and 1. Default is \code{conf = 0.95}.
#' @param type Method to compute the confidence interval if \code{boot = TRUE}. Default is "norm".
#' Check the \code{type} option in the \code{boot} function of the \code{\link{boot}} package.
#' @importFrom rlang ":="
#' @importFrom Rdpack reprompt
#' @importFrom stats var
#' @references
#' \insertAllCited{}
#'
#' @return A \code{sensobol} object.
#' @seealso Check the function \code{\link{boot}} for further details on the bootstrapping
#' with regards to the methods available for the computation of confidence intervals in the \code{type} argument.
#' @export
#'
#' @details Any first and total-order estimator can be combined with the appropriate sampling design.
#' Check Table 3 of the vignette for a summary of all possible combinations, and Tables 1 and 2 for a
#' mathematical description of the estimators. If the analyst mismatches estimators and sampling designs,
#' the function will generate an error and urge to redefine the sample matrices or the estimators.
#'
#' For all estimators except \insertCite{Azzini2020;textual}{sensobol}'s and \insertCite{Janon2014;textual}{sensobol}'s,
#' \code{sobol_indices()} calculates the sample mean as \deqn{\hat{f}_0=\frac{1}{2N} \sum_{v=1}^{N}(f(\mathbf{A})_v + f(\mathbf{B})_v)\,,}
#' where \eqn{N} is the row dimension of the base sample matrix, and the unconditional sample variance as
#'
#' \deqn{\hat{V}(y) = \frac{1}{2N-1} \sum{v=1}^{N} ((f(\mathbf{A})_v - \hat{f})^2 + (f(\mathbf{B})_v - \hat{f})^2)\,,}
#' where \eqn{f(\mathbf{A})_v} (\eqn{f(\mathbf{B})_v}) indicates the model output \eqn{y} obtained after running the model \eqn{f}
#' in the \eqn{v}-th row of the \eqn{\mathbf{A}} (\eqn{\mathbf{B}}) matrix.
#'
#' For the Azzini estimator,
#' \deqn{\hat{V}(y) = \sum_{v=1}^{N} (f(\mathbf{A})_v - f(\mathbf{B})_v)^2 + (f(\mathbf{B}_A^{(i)})_v - f(\mathbf{A}_B^{(i)})_v) ^ 2}
#'
#' and for the Janon estimator,
#' \deqn{\hat{V}(y)=\frac{1}{N} \sum_{v=1}^{N} \frac{f(\mathbf{A})_v^2 + f(\mathbf{A}_B^{(i)})_v^2}{2}-f_0^2}
#'
#'where \eqn{f(\mathbf{A}_B^{(i)})_v} (\eqn{f(\mathbf{B}_A^{(i)})_v}) is the model output obtained after running the model \eqn{f} in
#'the \eqn{v}-th row of an \eqn{\mathbf{A}_B^{(i)})_v} (\eqn{\mathbf{B}_A^{(i)})_v}) matrix, where all columns come from \eqn{\mathbf{A}} (\eqn{\mathbf{B}})
#'except the \eqn{i}-th, which comes from \eqn{\mathbf{B}} (\eqn{\mathbf{A}}).
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Compute and bootstrap Sobol' indices
#' ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = R)
sobol_indices <- function(matrices = c("A","B","BA"), Y, N, params,
                          first = "sobol", total = "saltelli",
                          boot = FALSE, R = NULL,
                          parallel = "no", ncpus = 1, conf = 0.95, type = "norm",
                          full,indep) {
  
  # Check concordance between boot and R arguments
  # ---------------------------------------------------------------------
  
  if (boot == FALSE & is.null(R) == FALSE | boot == TRUE & is.null(R) == TRUE) {
    stop("Bootstrapping requires boot = TRUE and an integer in R")
  }
  
  # Define parameters
  # ----------------------------------------------------------------------
  
  sensitivity <- parameters <- NULL
  k <- length(params)
  d <- matrix(Y, nrow = N)
  
  # Function when boot = FALSE
  # -----------------------------------------------------------------------
  
  if (boot == FALSE) {
    tmp <- sobol_boot(d = d, N = N, params = params, first = first, total = total, boot = FALSE, matrices = matrices)
    out <- data.table::data.table(tmp)
    #data.table::setnames(out, "tmp", "original")
    out <- list(Sobolindices=data.table::setnames(out, "tmp", "original"))
    
    # Function when boot = TRUE
    # -----------------------------------------------------------------------
    
  } else if (boot == TRUE) {
    tmp <- boot(data = d, statistic = sobol_boot, R = R, N = N, params = params,
                first = first, total = total, matrices = matrices,
                parallel = parallel, ncpus = ncpus, boot = TRUE)
    out <- list(Sobolindices=data.table::data.table(bootstats(tmp$boot0, conf = conf, type = type)),bootstrapsamples=tmp$bootstrapsamples)
    
  } else {
    stop("boot has to be TRUE or FALSE")
  }
  
  # Vectors of parameters and sensitivity indices
  # -------------------------------------------------------------------------
  
  parameters <- rep(c(full,indep),2)
  sensitivity <- c("Sifull","Si-1ind", "Tifull","Ti-1ind")
    
  # Create class and output
  # ----------------------------------------------------------------------
  
  ind <- structure(list(), class = "sensobol") # Create class
  ind$results <- cbind(Sobolindices=out$Sobolindices, sensitivity, parameters) # Add Sobol' indices
  original <- NULL
  ind$first <- first # First-order estimator
  ind$total <- total # Total-order estimator
  ind$C <- length(Y) # Total number of model runs
  ind$bootstrapsamples<-out$bootstrapsamples
  
  return(ind)
  
}

#' Display the results obtained with the \code{sobol_indices} function.
#'
#' @param x A \code{sensobol} object produced by \code{sobol_indices}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @return The function \code{print.sensobol} informs on the first and total-order
#' estimators used in the computations, the total number of model runs and
#' the sum of first-order index. It also plots the estimated results.
#' @export
#'
print.sensobol <- function(x, ...) {
  cat("\nFirst-order estimator:", x$first, "| Total-order estimator:", x$total, "\n")
  cat("\nTotal number of model runs:", x$C, "\n")
  cat("\nSum of first order indices:", x$si.sum, "\n")
  print(data.table::data.table(x$results))
}


# PERSONALIZED GGPLOT2 THEME
##################################################################################

theme_AP <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent",
                                                             color = NA),
                   legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
                   strip.background = ggplot2::element_rect(fill = "white"),
                   legend.position = "top")
}

# PLOT SOBOL' FIRST AND TOTAL-ORDER INDICES
##################################################################################
#' Visualization of first, total, second, third and fourth-order Sobol' indices.
#'
#' It plots first, total, second, third and fourth-order Sobol' indices.
#'
#' @param x The output of \code{\link{sobol_indices}}.
#' @param order If \code{order = "first"}, it plots first and total-order effects.
#' If \code{order = "second"}, it plots second-order effects. If \code{order = "third"}, it plots
#' third-order effects. If \code{order = "fourth"}, it plots
#' third-order effects. Default is \code{order = "first"}.
#' @param dummy The output of \code{\link{sobol_dummy}}. Default is NULL.
#' @param ... Other graphical parameters to plot.
#'
#' @return A \code{ggplot} object.
#' @rdname plot.sensobol
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Compute and bootstrap Sobol' indices
#' ind <- sobol_indices(Y = Y, N = N, params = params, boot = TRUE, R = R)
#'
#' # Plot Sobol' indices
#' plot(ind)

plot.sensobol <- function(x, order = "first", dummy = NULL, ...) {
  sensitivity <- parameters <- original <- low.ci <- high.ci <- NULL
  data <- x$results
  colNames <- colnames(data)
  
  # Plot only first-order indices
  # -----------------------------------------
  
  if (order == "first") {
    dt <- data[sensitivity %in% c("Si", "Ti")]
    gg <- ggplot2::ggplot(dt, ggplot2::aes(parameters, original, fill = sensitivity)) +
      ggplot2::geom_bar(stat = "identity",
                        position = ggplot2::position_dodge(0.6),
                        color = "black") +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
      ggplot2::labs(x = "",
                    y = "Sobol' index") +
      ggplot2::scale_fill_manual(values=c("black","grey"),
                                 name = "Sobol' indices",
                                 labels = c(expression(S[italic(i)]),
                                            expression(T[italic(i)])))+
      theme_AP()
    
    # Check if there are confidence intervals
    # -----------------------------------------
    
    if (any(grepl("high.ci", colNames)) == TRUE) {
      gg <- gg +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = low.ci,
                                            ymax = high.ci),
                               position = ggplot2::position_dodge(0.6))
    }
    
    # Check if there are indices for the dummy parameter
    # -----------------------------------------
    
    if (is.null(dummy) == FALSE) {
      col_names <- colnames(dummy)
      
      if(any(grepl("high.ci", col_names)) == TRUE) {
        lmt <- dummy$high.ci
        
      } else {
        lmt <- dummy$original
      }
      gg <- gg +
        ggplot2::geom_hline(data = dummy,
                            ggplot2::aes(yintercept = lmt, color = sensitivity),
                            lty = 2,lwd=2) +
        scale_color_manual(values=c("black","grey"))+
        ggplot2::guides(linetype = FALSE, color = FALSE)
    }
    
  } else if (!order == "first") {
    
    # Define for second and third-order indices
    # -----------------------------------------
    
    if (order == "second") {
      dt <- data[sensitivity %in% "Sij"][low.ci > 0]
      
    } else if (order == "third") {
      dt <- data[sensitivity %in% "Sijl"][low.ci > 0]
      
    } else if (order == "fourth") {
      dt <- data[sensitivity %in% "Sijlm"][low.ci > 0]
      
    } else {
      stop("Order should be first, second or third")
    }
    gg <- ggplot2::ggplot(dt, ggplot2::aes(stats::reorder(parameters, original),
                                           original)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = low.ci,
                                          ymax = high.ci)) +
      ggplot2::labs(x = "",
                    y = "Sobol' index") +
      ggplot2::geom_hline(yintercept = 0,
                          lty = 2,
                          color = "red") +
      theme_AP()
  }
  return(gg)
}


# PLOT MODEL OUTPUT UNCERTAINTY
##################################################################################

#' Visualization of the model output uncertainty
#'
#' It creates an histogram with the model output distribution.
#'
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Plot uncertainty
#' plot_uncertainty(Y = Y, N = N)

plot_uncertainty <- function(Y, N = NULL) {
  
  # Ensure that Y is a vector
  # -----------------------------------------
  
  if (is.vector(Y) == FALSE) {
    stop("Y should be a vector")
  }
  
  # Ensure that N is defined
  # -----------------------------------------
  
  if (is.null(N) == TRUE) {
    stop("The size of the base sample matrix N should be specified")
  }
  
  Y <- Y[1:N]
  df <- data.frame(Y)
  gg <- ggplot2::ggplot(df, ggplot2::aes(Y)) +
    ggplot2::geom_histogram(color = "black",
                            fill = "white") +
    ggplot2::labs(x = "y",
                  y = "Count") +
    theme_AP()
  return(gg)
}

# PLOT SCATTERPLOTS OF MODEL OUTPUT AGAINST MODEL INPUTS
#################################################################################

#' Scatter plots of the model output against the model inputs.
#'
#' It creates scatter plots of the model output against the model inputs.
#'
#' @param data The matrix created with \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{sobol_matrices}.
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param params Character vector with the name of the model inputs.
#' @param method The type of plot. If \code{method = "point"} (the default), each simulation is a point.
#' If \code{method = "bin"}, bins are used to aggregate simulations.
#' @param size Number between 0 and 1, argument of \code{geom_point()}. Default is 0.7.
#' @param alpha Number between 0 and 1, transparency scale of \code{geom_point()}. Default is 0.2.
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Plot scatter
#' plot_scatter(data = mat, Y = Y, N = N, params = params)
plot_scatter <- function(data, N, Y, params, method = "point", size = 0.7, alpha = 0.2) {
  value <- y <- NULL
  
  dt <- data.table::data.table(cbind(data, Y))[1:N]
  colnames(dt)[length(colnames(dt))] <- "y"
  out <- data.table::melt(dt, measure.vars = params)
  
  # Define the plot skeleton
  # -----------------------------------------
  
  gg <- ggplot2::ggplot(out, ggplot2::aes(value, y)) +
    ggplot2::facet_wrap(~variable, scales = "free_x") +
    ggplot2::labs(x = "Value", y = "y") +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent",
                                                             color = NA),
                   legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
                   strip.background = ggplot2::element_rect(fill = "white"),
                   legend.position = "top",
                   strip.text = element_text(size=25))
  
  # Precise for geom_point
  # -----------------------------------------
  
  if (method == "point") {
    gg <- gg + ggplot2::geom_point(size = size, alpha = alpha) +
      ggplot2::stat_summary_bin(fun = "mean", geom = "point", colour = "black", size = 6,shape=18)
    
    # Precise for geom_hex
    # -----------------------------------------
    
  } else if (method == "bin") {
    gg <- gg + ggplot2::geom_hex() +
      ggplot2::stat_summary_bin(fun = "mean", geom = "point", colour = "red", size = 3)
    
  } else {
    stop("Method should be either point or bin")
  }
  return(gg)
}

# PLOT SCATTERPLOT MATRIX OF PAIRS OF PARAMETERS
##################################################################################

#' Pairwise combinations of model inputs with the colour
#' proportional the model output value.
#'
#' It plots all pairwise combinations of model inputs with the colour
#' proportional the model output value.
#'
#' @param data The matrix created with \code{\link{sobol_matrices}}.
#' @param N Positive integer, the initial sample size of the base sample matrix created with \code{\link{sobol_matrices}}.
#' @param Y A numeric vector with the model output obtained from the matrix created with
#' \code{\link{sobol_matrices}}.
#' @param params Character vector with the name of the model inputs.
#' @param smpl The number of simulations to plot.
#' The default is NULL.
#'
#' @importFrom data.table .SD .N
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Define settings
#' N <- 1000; params <- paste("X", 1:3, sep = ""); R <- 10
#'
#' # Create sample matrix
#' mat <- sobol_matrices(N = N, params = params)
#'
#' # Compute Ishigami function
#' Y <- ishigami_Fun(mat)
#'
#' # Plot scatterplot matrix
#' plot_multiscatter(data = mat, N = N, Y = Y, params = params)
plot_multiscatter <- function(data, N, Y, params, smpl = NULL) {
  xvar <- yvar <- x <- y <- NULL
  dt <- data.table::data.table(data)
  out <- t(utils::combn(params, 2))
  da <- list()
  
  # Define pairwise combinations
  # -----------------------------------------
  
  for (i in 1:nrow(out)) {
    cols <- out[i, ]
    da[[i]] <- cbind(dt[1:N, .SD, .SDcols = (cols)], cols[1], cols[2], Y[1:N])
    data.table::setnames(da[[i]], colnames(da[[i]]), c("xvar", "yvar", "x", "y", "output"))
  }
  
  output <- data.table::rbindlist(da)
  
  # Define option to plot just a fraction of the sample
  # -----------------------------------------
  
  if (is.null(smpl) == FALSE) {
    if (is.numeric(smpl) == FALSE) {
      stop("smpl should be a number")
      
    } else {
      output <- output[,.SD[sample(.N, min(smpl,.N))], by = list(x, y)]
    }
  }
  
  # Plot
  # -----------------------------------------
  
  gg <- ggplot2::ggplot(output, ggplot2::aes(xvar, yvar, color = output)) +
    ggplot2::geom_point(size = 0.5) +
    ggplot2::scale_colour_gradientn(colours = grDevices::terrain.colors(10), name = "y") +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
    ggplot2::facet_wrap(x~y, scales = "free") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.background = ggplot2::element_rect(fill = "transparent",
                                                             color = NA),
                   legend.key = ggplot2::element_rect(fill = "transparent", color = NA),
                   strip.background = ggplot2::element_rect(fill = "white"),
                   legend.position = "top")
  return(gg)
}