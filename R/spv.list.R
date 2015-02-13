#' @rdname spv
#' @method spv list
#' @export
spv.list <- function(n, design, type = "spherical", formula, at = FALSE, 
                     keepfun, sample, unscaled = FALSE, ...){
  cll <- match.call()
  type <- match.arg(type, c("spherical", "cuboidal", "lhs", "mlhs", "slhs", 
                            "rslhs"))
  nr <- length(design)
  if(missing(sample)){
    sample <- sampler(n = n, design = design[[1]], type = type, at = at, ...)
    if(!missing(keepfun)) {
      repeat{
        keep <- keepfun(sample)
        cnt <- sum(keep)
        sample <- sample[keep, ]
        if(cnt >= n) break
        rate <- cnt/n
        cat("Retained samples:",  round(cnt, digits = 2), 
            "-- Adding some more...\n")
        addsample <- sampler(n = max(ceiling((n - cnt)/rate), ceiling(n/10)), 
                             design = design[[1]], type = type, at = at, ...)
        sample <- rbind(sample, addsample)
      }
      cat("Final sample of size", nrow(sample))
    }
  }
  spvdesign <- function(design, sample, formula, call, unscaled){
    ndes <- nrow(design)
    n <- nrow(sample)
    mat <- model.matrix(formula, data = as.data.frame(sample))
    m <- ncol(design)
    mod.mat <- model.matrix(formula, data = as.data.frame(design))
    p <- ncol(mod.mat)
    FtF.inv <- solve(crossprod(mod.mat))
    tmp <- .Fortran("fds", as.integer(p), as.integer(n), as.integer(ndes), 
                    as.double(FtF.inv), as.double(mat), double(n), 
                    PACKAGE = "vdg")
    spv <- tmp[[6]]
    if(unscaled) spv <- spv / ndes
    out <- list(spv = spv, sample = sample, type = type, call = call, at = at, 
                FtF.inv = FtF.inv, formula = formula, ndes = ndes, 
                unscaled = unscaled)
    class(out) <- c("spv", "list")
    out
  }
  cl <- makeCluster(getOption("cl.cores", min(detectCores(), nr)))  
  clusterEvalQ(cl, library(vdg))
  if(is(formula, "formula")){
    out <- parLapply(cl, design, spvdesign, sample = sample, formula = formula, 
                     call = cll, unscaled = unscaled)
    stopCluster(cl)
    nms <- names(design)
    if(is.null(nms)) names(out) <- paste0("Des", seq_along(design))
    else names(out) <- nms
    class(out) <- c("spvlist", "list")
    return(out)
  }
  if(is.list(formula)){
    nf <- length(formula)
    out <- lapply(formula, function(y) {
      out <- parLapply(cl, design, spvdesign, sample = sample, formula = y, 
                       call = cll, unscaled = unscaled)
      class(out) <- c("spvlist", "list")
      return(out)
      })
    stopCluster(cl)
    class(out) <- c("spvlistforlist", "list")
    return(out)
  }
}