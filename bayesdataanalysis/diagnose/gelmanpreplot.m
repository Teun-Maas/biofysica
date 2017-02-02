
  function gelmanpreplot(x,varargin) 

			binwidth	= keyval('binwidth',varargin,10);
			maxbins		= keyval('maxbins',varargin,50);
            confidence	= keyval('confidence',varargin,0.95);
			transform	= keyval('transform',varargin,false);
            autoburnin	= keyval('autoburnin',varargin,true);
			
			


		
  x <- as.mcmc.list(x)
  if (niter(x) <= 50) 
    stop("Less than 50 iterations in chain")
  nbin <- min(floor((niter(x) - 50)/thin(x)), max.bins)
  binw <- floor((niter(x) - 50)/nbin)
  last.iter <- c(seq(from = start(x) + 50 * thin(x), by = binw * 
                     thin(x), length = nbin), end(x))
  shrink <- array(dim = c(nbin + 1, nvar(x), 2))
  dimnames(shrink) <- list(last.iter, varnames(x),
                           c("median", paste(50 * (confidence + 1), "%",
                                             sep = ""))
                           )
  for (i in 1:(nbin + 1)) {
    shrink[i, , ] <- gelman.diag(window(x, end = last.iter[i]), 
                                 confidence = confidence,
                                 transform = transform,
                                 autoburnin = autoburnin)$psrf
  }
  all.na <- apply(is.na(shrink[, , 1, drop = FALSE]), 2, all)
  if (any(all.na)) {
    cat("\n******* Error: *******\n")
    cat("Cannot compute Gelman & Rubin's diagnostic for any chain \n")
    cat("segments for variables", varnames(x)[all.na], "\n")
    cat("This indicates convergence failure\n")
  }
  return(list(shrink = shrink, last.iter = last.iter))
