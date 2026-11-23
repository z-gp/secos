a <- commandArgs(trailingOnly = TRUE); compito <- a[1]; x <- readRDS(a[2])
tempo <- function(f) { invisible(gc()); t1 <- proc.time()[["elapsed"]]; f(); proc.time()[["elapsed"]] - t1 }
ris <- switch(compito,
  secos = {
    source("R/secos.R")
    z <- secos_map(x$E, n_neighbors = x$k, min_dist = x$md, replicates = x$R, n_threads = 1L)
    c(z$timing, fidelity = tempo(function() secos_fidelity(z)), components = z$settings$components)
  },
  mca = { loadNamespace("ca"); tempo(function() ca::mjca(data.frame(apply(x$B, 2, factor)), nd = 2, lambda = "indicator", ps = "_")) },
  ca  = { loadNamespace("ca"); tempo(function() ca::ca(x$B, nd = 2)) },
  cs  = {
    suppressPackageStartupMessages(library(bibliometrix))
    arg <- list(M = x$M, field = "DE", method = "MCA", minDegree = x$s, clust = "auto", k.max = 3, stemming = FALSE,
                remove.terms = x$stop, documents = 5, graph = FALSE)
    tempo(function() do.call(conceptualStructure, arg[names(arg) %in% names(formals(conceptualStructure))]))
  },
  stop("unknown task"))
saveRDS(ris, a[3])
