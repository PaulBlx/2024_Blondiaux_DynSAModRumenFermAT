# --------------------------------------------------------------------------------------------------------------------

# Functions computing the Shapley effects using the random permutation method (Song et al., 2016)(package sensitivity)

# --------------------------------------------------------------------------------------------------------------------

shapleyPermRand <- function(model = NULL, Xall, Xset, d, Nv, m, No = 1, Ni = 3, colnames = NULL, ...) {
  #################################################################################
  # Authors: 
  # First version by: Eunhye Song, Barry L. Nelson, Jeremy Staum (Northwestern University), 2015
  # Modified for 'sensitivity' package by: Bertrand Iooss (EDF R&D), 2017
  #################################################################################
  
  #################################################################################
  # This function implements Algorithm 1 to calculate the Shapley effects and
  # their standard errors by randomly sampling permutations of inputs.
  # It also estimates First order and total Sobol' indices 
  #
  # List of inputs to this function:
  # model: a function, or a model with a predict method, defining the model to analyze
  # Xall(n): a function to generate a n-sample of a d-dimensional input vector
  # Xset(n, Sj, Sjc, xjc): a function to generate a n- sample an input vector corresponding 
  #   to the indices in Sj conditional on the input values xjc with the index set Sjc
  # d: number of inputs
  # Nv: Monte Carlo (MC) sample size to estimate the output variance
  # m: number of randomly sampled permutations
  # No: outer MC sample size to estimate the cost function
  # Ni: inner MC sample size to estimate the cost function
  #
  # This function requires R package "gtools" 
  #
  #################################################################################
  
  # Generate m permutations
  perms <- matrix(NA, ncol=d, nrow=m)
  #for (p in 1:m) perms[p,] <- gtools::permute(1:d)
  for (p in 1:m) perms[p,] <- gtools::permute(colnames)
  
  
  X <- matrix(NA, ncol=d, nrow=Nv+m*(d-1)*No*Ni)
  if (is.null(colnames)) for (i in 1:d) colnames <- c(colnames,paste("X",i,sep=""))
  colnames(X) <- colnames
  
  X[1:Nv,] <- Xall(Nv)
  
  for (p in 1:m)
  {
    pi <- perms[p,]
    #pi_s <- sort(pi,index.return=TRUE)$ix
    pi_s <- match(colnames(X),pi)
    
    for (j in (1:(d-1)))
    {
      Sj <- pi[c(1:j)] # set of the 1st-jth elements in pi 
      Sjc <- pi[-c(1:j)] # set of the (j+1)th-dth elements in pi
      
      xjcM <- matrix(Xset(No, Sjc, NULL, NULL),nrow=No) # sampled values of the inputs in Sjc
      for (l in 1:No)
      {
        xjc <- xjcM[l,]
        
        # sample values of inputs in Sj conditional on xjc
        xj <- Xset(Ni, Sj, Sjc, xjc)
        xx <- cbind(xj, matrix(xjc,nrow=Ni,ncol=length(xjc),byrow=T))
        X[(Nv+(p-1)*(d-1)*No*Ni+(j-1)*No*Ni+(l-1)*Ni+1):(Nv+(p-1)*(d-1)*No*Ni+(j-1)*No*Ni+l*Ni),] <- xx[,pi_s]
      }
    }
  }
  
  sh <- list(model = model, Xall = Xall, xset = Xset, d = d, Nv = Nv, m = m, No = No, Ni = Ni, X = X, 
             colnames = colnames, perms = perms, call = match.call())
  class(sh) <- "shapleyPermRand"
  
  if (!is.null(sh$model)) {
    response(sh, ...)
    tell(sh)
  }
  
  return(list(SamplingInputparameters=X,sh=sh))
  
}


tell.shapleyPermRand <- function(x, y = NULL, return.var = NULL, ...) {
  id <- deparse(substitute(x))
  
  if (! is.null(y)) {
    x$y <- y
  } else{ 
    y <- x$y
    if (is.null(x$y)) stop("y not found")
  }
  
  d <- x$d ; m <- x$m
  
  # Initialize Shapley value for all players
  Sh <- rep(0, d) ; Sh2 <- rep(0, d)
  
  # Initialize main and total (Sobol) effects for all players
  Vsob <- rep(0, d) ; Vsob2 <- rep(0, d) ; nV <- rep(0, d)
  Tsob <- rep(0, d) ; Tsob2 <- rep(0, d) ; nT <- rep(0, d)
  
  # Estimate Var[Y] 
  Y <- y[1:x$Nv] ; y <- y[-(1:x$Nv)]
  EY <- mean(Y)
  VarY <- var(Y)
  
  # Estimate Shapley effects
  for (p in 1:m)
  {
    pi <- match(x$perms[p,],colnames(x$X))
    prevC <- 0
    for (j in 1:d)
    {
      if (j == d)
      {    
        Chat <- VarY
        del <- Chat - prevC
        Vsob[pi[j]] <- Vsob[pi[j]] + prevC # first order effect
        Vsob2[pi[j]] <- Vsob2[pi[j]] + prevC^2
        nV[pi[j]] <- nV[pi[j]] + 1
      }
      else
      {
        cVar <- NULL
        for (l in 1:x$No)
        {
          Y <- y[1:x$Ni] ; y <- y[-(1:x$Ni)]
          cVar <- c(cVar,var(Y))
        }
        Chat <- mean(cVar)
        del <- Chat - prevC
      }
      
      Sh[pi[j]] <- Sh[pi[j]] + del
      Sh2[pi[j]] <- Sh2[pi[j]] + del^2
      
      prevC <- Chat
      
      if (j == 1){
        Tsob[pi[j]] <- Tsob[pi[j]] + Chat # Total effect
        Tsob2[pi[j]] <- Tsob2[pi[j]] + Chat^2
        nT[pi[j]] = nT[pi[j]] + 1
      }
    }
  }
  Sh <- Sh / m / VarY
  Sh2 <- Sh2 / m / VarY^2
  ShSE <- sqrt((Sh2 - Sh^2) / m)
  
  Vsob <- Vsob / nV / VarY # averaging by number of permutations with j=d-1
  Vsob2 <- Vsob2 / nV / VarY^2
  VsobSE <- sqrt((Vsob2 - Vsob^2) / nV)
  Vsob <- 1 - Vsob 
  Vsob2 <- 1 - Vsob2 
  
  Tsob <- Tsob / nT / VarY # averaging by number of permutations with j=1
  Tsob2 <- Tsob2 / nT / VarY^2
  TsobSE <- sqrt((Tsob2 - Tsob^2) / nT)
  
  Shapley <- data.frame(cbind(Sh,ShSE,Sh-2*ShSE,Sh+2*ShSE),row.names=x$colnames)
  names(Shapley) <- c("original","std. error", "min. c.i.", "max. c.i.")
  Vsobol <- data.frame(cbind(Vsob,VsobSE,Vsob-1.96*VsobSE,Vsob+1.96*VsobSE),row.names=x$colnames)
  names(Vsobol) <- c("original","std. error", "min. c.i.", "max. c.i.")
  Tsobol <- data.frame(cbind(Tsob,TsobSE,Tsob-1.96*TsobSE,Tsob+1.96*TsobSE),row.names=x$colnames)
  names(Tsobol) <- c("original","std. error", "min. c.i.", "max. c.i.")
  
  x$Shapley <- Shapley
  x$SobolS <- Vsobol
  x$SobolT <- Tsobol
  x$V <- VarY
  x$E <- EY
  
  for (i in return.var) {
    x[[i]] <- get(i)
  }
  
  assign(id, x, parent.frame())
  
  return(x)
}


print.shapleyPermRand <- function(x, ...) {
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  if (!is.null(x$y)) {
    cat("\nModel runs:", length(x$y), "\n")
    cat("\nShapley' effects:\n")
    print(x$Shapley)
    cat("\nFull first order Sobol' indices:\n")
    print(x$SobolS)
    cat("\nIndependent total Sobol' indices:\n")
    print(x$SobolT)
  }
}


plot.shapleyPermRand <- function(x, ylim = c(0, 1), ...) {
  if (!is.null(x$y)) {
    pch = c(21, 24, 25)
    nodeplot(x$Shapley, xlim = c(1, x$d + 1), ylim = ylim, pch = pch[1])
    nodeplot(x$SobolS, xlim = c(1, x$d + 1), ylim = ylim, labels = FALSE,
             pch = pch[2], at = (1:x$d)+.2, add = TRUE)
    nodeplot(x$SobolT, xlim = c(1, x$d + 1), ylim = ylim, labels = FALSE,
             pch = pch[3], at = (1:x$d)+.4, add = TRUE)
    legend(x = "topright", legend = c("Shapley effect","Full first Sobol'", "Independent total Sobol'"), pch = pch)
  }
}

ggplot.shapleyPermRand <- function(x, ylim = c(0, 1), title = NULL, ...) {
  if (!is.null(x$y)) {
    pch = c(21, 24, 25)
    nodeggplot(list(x$Shapley,x$SobolS,x$SobolT), xname=c("Shapley effect","Full first Sobol'", "Independent total Sobol'"), ylim = ylim, title = title, pch = pch)
  }
}

#                       Nodeplot: anti-boxplot
#                         Gilles Pujol 2006
# Modified by Frank Weber (2016): support model functions
# returning a matrix or a 3-dimensional array.

nodeplot <- function(x, xlim = NULL, ylim = NULL, labels = TRUE,
                     col = par("col"), pch = 21, bg = "white",
                     add = FALSE, at = NULL, y_col = NULL, y_dim3 = NULL, ...) {
  n <- nrow(x)
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (add) {
    par(new = TRUE)
  }
  
  # axes
  
  plot(0, xlim = xlim, ylim = ylim, axes = FALSE,
       xlab = "", ylab = "", type = "n", ...)
  #  if (class(labels) == "logical") {
  if (inherits(labels, "logical")){
    if (labels) {
      axis(side = 1, at = at, labels = rownames(x))
    } else {
      axis(side = 1, at = at, labels = FALSE, tick = FALSE)
    }
    #  } else if (class(labels) == "character") {
  } else if (inherits(labels, "character")){
    axis(side = 1, at = at, labels = labels)
  }
  axis(side = 2)
  box()
  
  # bias
  
  if ("bias" %in% dimnames(x)[[2]]) {
    if(is.null(y_col) && is.null(y_dim3)){
      xx <- x[["original"]] - x[["bias"]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      xx <- x[, "original", y_col] - x[, "bias", y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      xx <- x[, "original", y_col, y_dim3] - x[, "bias", y_col, y_dim3]
    }
  } else {
    if(is.null(y_col) && is.null(y_dim3)){
      xx <- x[["original"]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      xx <- x[, y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      xx <- x[, y_col, y_dim3]
    }
  }
  
  # confidence intervals
  
  if (("min. c.i." %in% dimnames(x)[[2]]) & "max. c.i." %in% dimnames(x)[[2]]) {
    if(is.null(y_col) && is.null(y_dim3)){
      min_ci <- x[["min. c.i."]]
      max_ci <- x[["max. c.i."]]
    } else if(!is.null(y_col) && is.null(y_dim3)){
      min_ci <- x[, "min. c.i.", y_col]
      max_ci <- x[, "max. c.i.", y_col]
    } else if(!is.null(y_col) && !is.null(y_dim3)){
      min_ci <- x[, "min. c.i.", y_col, y_dim3]
      max_ci <- x[, "max. c.i.", y_col, y_dim3]
    }
    for (i in 1 : n) {
      lines(c(at[i], at[i]), c(min_ci[i], max_ci[i]),
            col = col)
    }
  }
  
  # points
  
  points(at, xx, col = col, pch = pch, bg = bg)
}

#                       Nodeggplot: anti-boxplot in ggplot
#                         Sebastien Da Veiga (June 2019, Seignosse)

if (getRversion() >= "2.15.1") utils::globalVariables(c("x","y","id"))

nodeggplot <- function(listx, xname, xlim = NULL, ylim = NULL, labels = TRUE, title = NULL, 
                       col = par("col"), pch = 21, at = NULL, y_col = NULL, y_dim3 = NULL, ...) {
  
  ngraphs <- length(listx)
  x <- unlist(listx)
  n <- nrow(listx[[1]])
  if (is.null(xlim)) {
    xlim <- c(1, n)
  }
  if (n<=10){
    angle <- 0
    hjust <- 0
  }
  if (n>10 & n <=20){
    angle <- 45
    hjust <- 1
  }
  if (n>20){
    angle <- 90
    hjust <- 1
  }
  if (is.null(ylim)) {
    ylim <- c(min(x), max(x))
  }
  if (is.null(at)) {
    at <- 1 : n
  }
  if (is.null(title)){
    title <- title
  }
  
  #  if (class(labels) == "logical") {
  if (inherits(labels, "logical")){
    if (labels) {
      l <- rownames(listx[[1]])
    } else {
      l <- NULL
    }
    #  } else if (class(labels) == "character") {
  } else if (inherits(labels, "character")){
    l <- labels
  }
  
  # bias
  
  d <- NULL
  for (i in 1:ngraphs){
    x <- listx[[i]]
    if ("bias" %in% dimnames(x)[[2]]) {
      if(is.null(y_col) && is.null(y_dim3)){
        xx <- x[["original"]] - x[["bias"]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        xx <- x[, "original", y_col] - x[, "bias", y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        xx <- x[, "original", y_col, y_dim3] - x[, "bias", y_col, y_dim3]
      }
    } else {
      if(is.null(y_col) && is.null(y_dim3)){
        xx <- x[["original"]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        xx <- x[, y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        xx <- x[, y_col, y_dim3]
      }
    }
    d <- rbind(d,data.frame(x=at,y=xx,id=xname[i]))
  }
  
  # confidence intervals
  d2 <- NULL
  n2 <- 0
  for (i in 1:ngraphs){
    x <- listx[[i]]
    if (("min. c.i." %in% dimnames(x)[[2]]) & "max. c.i." %in% dimnames(x)[[2]]) {
      if(is.null(y_col) && is.null(y_dim3)){
        min_ci <- x[["min. c.i."]]
        max_ci <- x[["max. c.i."]]
      } else if(!is.null(y_col) && is.null(y_dim3)){
        min_ci <- x[, "min. c.i.", y_col]
        max_ci <- x[, "max. c.i.", y_col]
      } else if(!is.null(y_col) && !is.null(y_dim3)){
        min_ci <- x[, "min. c.i.", y_col, y_dim3]
        max_ci <- x[, "max. c.i.", y_col, y_dim3]
      }
      n2 <- n2 +1
    }else{
      min_ci <- rep(NA,n)
      max_ci <- rep(NA,n)
    }
    d2 <- rbind(d2,data.frame(min_ci=min_ci,max_ci=max_ci))
  }
  
  d <- cbind(d,d2)
  
  if (ngraphs>1){
    pd <- position_dodge(0.3)
    if (n2>0){
      g <- ggplot(d, aes(x=x, y=y, shape=id)) + 
        geom_point(size=3, colour=col, position = pd) +
        geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width=.1, position = pd) + 
        scale_shape_manual(values=pch) +
        coord_cartesian(ylim=ylim) +
        labs(y="", x = "", shape = title) +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), legend.position = c(0.8, 0.9))
    }else{
      g <- ggplot(d, aes(x=x, y=y, shape=id)) + 
        geom_point(size=3, colour=col, position = pd) +
        scale_shape_manual(values=pch) +
        coord_cartesian(ylim=ylim) +
        labs(y="", x = "", shape = title) +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), legend.position = c(0.8, 0.9))
    }
  }else{
    if (n2>0){
      g <- ggplot(d, aes(x=x, y=y)) + 
        geom_point(size=3, shape=pch, colour=col) +
        geom_errorbar(aes(ymin=min_ci, ymax=max_ci), width=.1, colour=col) + 
        coord_cartesian(ylim=ylim)+
        labs(title= xname, y="", x = "") +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 15))
    }else{
      g <- ggplot(d, aes(x=x, y=y)) + 
        geom_point(size=3, shape=pch, colour=col) +
        coord_cartesian(ylim=ylim)+
        labs(title= xname, y="", x = "") +
        scale_x_discrete(limits=l) + 
        theme_bw() +
        theme(axis.text.x = element_text(face = "bold", size = 12, angle = angle, hjust = hjust), axis.text.y = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 15)) 
    }
  }
  return(g)
}