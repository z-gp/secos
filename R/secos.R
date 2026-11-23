# =============================================================================
# secos.R - Semantic Conceptual Structure (SeCoS)
#
# Maps the author keywords of a bibliographic collection from the meaning of the
# terms, instead of from how often they co-occur. The functions work on the data
# frame produced by bibliometrix::convert2df(), and the language model is a
# parameter: the default is a general model, which a user working in a specific
# field can replace.
#
#   secos()              whole procedure, from the bibliometrix data frame to the map
#   secos_terms()        terms above a minimum document frequency
#   secos_embed()        embeddings of the isolated terms (needs Python)
#   secos_cooccurrence() co-occurrence distance, for the hybrid version
#   secos_map()          consensus map from a matrix of embeddings
#   secos_stability()    neighborhood stability across replicates
#   secos_fidelity()     trustworthiness, continuity and related measures
#
# Requires: uwot, Matrix, igraph. Embeddings also require reticulate and a Python
# environment with torch and transformers; they can be computed elsewhere and
# passed as a matrix instead.
#
# Example
#   library(bibliometrix)
#   M <- convert2df("savedrecs.txt", dbsource = "wos", format = "plaintext")
#   reticulate::use_condaenv("my-env")
#   x <- secos(M, field = "DE", min_docs = 50, remove_terms = c("sport", "injury"))
#   x                                # summary
#   plot(x, labels = 25)             # map
#   head(x$map)                      # coordinates, ellipses, stability, uncertainty
#   secos_fidelity(x)                # fidelity with respect to the embedding space
#
# Another model, here SciBERT with the mean over the tokens of the term:
#   x <- secos(M, model = "allenai/scibert_scivocab_uncased", layer = -2L, pooling = "term")
#
# Hybrid version, half semantic distance and half co-occurrence distance:
#   x <- secos(M, lambda = 0.5)
# =============================================================================

SECOS_MODEL <- "sentence-transformers/all-mpnet-base-v2"
SECOS_REVISION <- "e8c3b32edf5434bc2275fc9bab85f82640a19130"

# ---- Internal helpers --------------------------------------------------------

.secos_keywords <- function(M, field, sep) {
  if (!is.data.frame(M) || !field %in% names(M))
    stop("M must be a data frame with a column named '", field, "'")
  x <- as.character(M[[field]])
  x[is.na(x)] <- ""
  lapply(strsplit(tolower(x), sep, fixed = TRUE), trimws)
}

.secos_cosine <- function(X) {                     # cosine distance between the rows
  X <- as.matrix(X)
  Z <- X / sqrt(rowSums(X^2))
  D <- 1 - tcrossprod(Z)
  D[D < 0] <- 0
  diag(D) <- 0
  dimnames(D) <- list(rownames(X), rownames(X))
  D
}

.secos_neighbours <- function(D, k) {              # ties go to the lowest index
  n <- nrow(D)
  diag(D) <- Inf
  N <- t(vapply(seq_len(n), function(i) order(D[i, ], seq_len(n))[seq_len(k)], integer(k)))
  if (k == 1) N <- matrix(N, ncol = 1)
  N
}

.secos_standardise <- function(D) {                # mean 0 and SD 1 over the pairs of terms
  x <- D[lower.tri(D)]
  S <- (D - mean(x)) / stats::sd(x)
  diag(S) <- 0
  S
}

.secos_procrustes <- function(maps, iterations = 10) {   # generalised, without scaling
  Y <- lapply(maps, function(y) sweep(y, 2, colMeans(y)))
  reference <- Y[[1]]
  for (it in seq_len(iterations)) {
    Y <- lapply(Y, function(y) { s <- svd(crossprod(y, reference)); y %*% (s$u %*% t(s$v)) })
    new <- Reduce(`+`, Y) / length(Y)
    if (max(abs(new - reference)) < 1e-9) break
    reference <- new
  }
  Y
}

.secos_median <- function(Y) {
  M <- apply(simplify2array(Y), c(1, 2), stats::median)
  dimnames(M) <- list(rownames(Y[[1]]), NULL)
  M
}

# Do the positions of a term across the replicates form two separate zones?
.secos_zones <- function(P, seed) {
  km <- tryCatch({ set.seed(seed); stats::kmeans(P, centers = 2, nstart = 20, iter.max = 50) },
                 error = function(e) NULL)
  if (is.null(km)) return(c(D = 0, share = 0))
  u <- km$centers[2, ] - km$centers[1, ]
  p <- drop(P %*% (u / sqrt(sum(u^2))))
  g <- km$cluster
  D <- sqrt(2) * abs(diff(tapply(p, g, mean))) / sqrt(sum(tapply(p, g, stats::var)))
  c(D = unname(D), share = min(table(g)) / nrow(P))
}

.secos_ranks <- function(D) {                      # rank of j among the neighbors of i
  n <- nrow(D)
  R <- matrix(0L, n, n)
  for (i in seq_len(n)) {
    d <- D[i, ]
    d[i] <- -Inf
    R[i, order(d, seq_len(n))] <- seq_len(n) - 1L
  }
  R
}

# ---- Terms ---------------------------------------------------------------------

#' Terms above a minimum document frequency
#'
#' @param M Data frame from bibliometrix::convert2df().
#' @param field Column with the keywords, by default the author keywords.
#' @param min_docs Minimum number of documents in which a term appears.
#' @param remove_terms Terms to leave out, such as those used in the query.
#' @param sep Separator between keywords.
#' @return Data frame with columns term and documents, by decreasing frequency.
secos_terms <- function(M, field = "DE", min_docs = 50, remove_terms = NULL, sep = ";") {
  kw <- lapply(.secos_keywords(M, field, sep), unique)
  f <- table(unlist(kw, use.names = FALSE))
  f <- f[!names(f) %in% c("", "na", tolower(trimws(remove_terms)))]
  f <- f[f >= min_docs]
  if (length(f) == 0) stop("no term reaches min_docs = ", min_docs)
  out <- data.frame(term = names(f), documents = as.integer(f), stringsAsFactors = FALSE)
  out <- out[order(-out$documents, out$term, method = "radix"), ]
  rownames(out) <- NULL
  out
}

# ---- Embeddings ----------------------------------------------------------------

.secos_python <- '
import time
import numpy as np, torch
from transformers import AutoTokenizer, AutoModel

def secos_embed_py(terms, model, revision, layer, pooling, batch_size, threads):
    if threads is not None: torch.set_num_threads(int(threads))
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(model, revision=revision)
    mod = AutoModel.from_pretrained(model, revision=revision, output_hidden_states=True)
    mod.eval()
    t1 = time.time()
    out, n_tokens = [], []
    with torch.no_grad():
        for i in range(0, len(terms), int(batch_size)):
            enc = tok(list(terms[i:i + int(batch_size)]), padding=True, return_tensors="pt",
                      return_special_tokens_mask=True)
            special = enc.pop("special_tokens_mask").bool()
            attention = enc["attention_mask"].bool()
            h = mod(**enc).hidden_states[int(layer)]
            term = attention & ~special
            mask = term if pooling == "term" else attention
            n_tokens += term.sum(1).tolist()
            w = mask.unsqueeze(-1).float()
            out.append(((h * w).sum(1) / w.sum(1)).numpy())
    t2 = time.time()
    return dict(E=np.vstack(out).astype("float64"), load=t1 - t0, encode=t2 - t1, n_tokens=n_tokens,
                commit=str(getattr(mod.config, "_commit_hash", None)))
'

#' Embeddings of the isolated terms
#'
#' Calls transformers in Python through reticulate; choose the Python environment
#' first, for example with reticulate::use_condaenv(). Hyphens and underscores
#' become spaces before encoding.
#'
#' @param terms Character vector of terms.
#' @param model Hugging Face model. The default is a general sentence encoder.
#' @param revision Commit of the model, so that results can be repeated. NULL takes
#'   the latest version.
#' @param layer Layer to use: -1 is the last, -2 the last but one.
#' @param pooling "all" averages over all the tokens, as sentence encoders are
#'   trained to do; "term" averages only over the tokens of the term, leaving out
#'   [CLS] and [SEP], which suits models such as SciBERT.
#' @param batch_size Terms per batch.
#' @param threads Threads for torch; NULL leaves the default.
#' @return Matrix with one row per term. Attributes: commit, n_tokens, and timing,
#'   the seconds taken to load the model and to encode the terms.
secos_embed <- function(terms, model = SECOS_MODEL,
                        revision = if (identical(model, SECOS_MODEL)) SECOS_REVISION else NULL,
                        layer = -1L, pooling = c("all", "term"), batch_size = 32L, threads = NULL) {
  pooling <- match.arg(pooling)
  if (!requireNamespace("reticulate", quietly = TRUE)) stop("secos_embed() needs the reticulate package")
  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE", TOKENIZERS_PARALLELISM = "false")
  if (!reticulate::py_module_available("transformers"))
    stop("the Python in use has no transformers module: choose an environment with reticulate::use_condaenv()")
  if (!reticulate::py_has_attr(reticulate::import_main(), "secos_embed_py"))
    reticulate::py_run_string(.secos_python)
  r <- reticulate::py$secos_embed_py(as.list(gsub("[-_]+", " ", terms)), model, revision, as.integer(layer),
                                     pooling, as.integer(batch_size), threads)
  E <- as.matrix(r$E)
  dimnames(E) <- list(terms, NULL)
  attr(E, "commit") <- r$commit
  attr(E, "n_tokens") <- as.integer(unlist(r$n_tokens))
  attr(E, "timing") <- c(load = r$load, encode = r$encode)
  E
}

# ---- Co-occurrence -------------------------------------------------------------

#' Co-occurrence distance between terms
#'
#' Cosine distance between the profiles of association strength,
#' c_ij / (f_i f_j). Terms that never co-occur with another term are at distance 1
#' from all the others.
#'
#' @param M Data frame from bibliometrix::convert2df().
#' @param terms Character vector of terms.
#' @param field Column with the keywords.
#' @param sep Separator between keywords.
#' @return List: D (distances), counts (co-occurrences, with the document
#'   frequencies on the diagonal), documents (documents with at least one term).
secos_cooccurrence <- function(M, terms, field = "DE", sep = ";") {
  kw <- .secos_keywords(M, field, sep)
  pos <- lapply(kw, function(v) { m <- unique(match(v, terms)); m[!is.na(m)] })
  X <- Matrix::sparseMatrix(i = rep(seq_along(pos), lengths(pos)), j = unlist(pos), x = 1,
                            dims = c(length(pos), length(terms)), dimnames = list(NULL, terms))
  X <- X[Matrix::rowSums(X) > 0, , drop = FALSE]
  C <- as.matrix(Matrix::crossprod(X))
  f <- diag(C)
  A <- C / outer(pmax(f, 1), pmax(f, 1))
  diag(A) <- 0
  dimnames(A) <- list(terms, terms)
  D <- suppressWarnings(.secos_cosine(A))
  D[is.na(D)] <- 1
  diag(D) <- 0
  list(D = D, counts = C, documents = nrow(X))
}

# ---- Map ------------------------------------------------------------------------

#' Consensus map from a matrix of embeddings
#'
#' Runs UMAP several times with exact neighbors and a fixed spectral start, aligns
#' the replicates with generalised Procrustes (no scaling) and takes the median
#' coordinates. Each term also gets a 95% ellipse and a stability value. A term has
#' an uncertain position if its ellipse is large or its replicates fall in two
#' separate zones (Ashman's D and minimum share of the smaller zone); such terms
#' should not be interpreted by their position.
#'
#' With lambda < 1 the distance is the convex combination
#' lambda * D_sem + (1 - lambda) * D_cooc, after standardising the two distances
#' over the pairs of terms. With lambda = 1 only the embeddings are used.
#'
#' @param E Matrix of embeddings, with the terms as row names.
#' @param D_cooc Co-occurrence distances (from secos_cooccurrence()), needed if lambda < 1.
#' @param lambda Weight of the semantic distance, between 0 and 1.
#' @param n_neighbors,min_dist UMAP parameters.
#' @param replicates Number of UMAP runs.
#' @param seeds One seed per replicate. The function always uses R's default random
#'   number generator and leaves the state of the session as it was.
#' @param area_limit Share of the map above which an ellipse makes a position uncertain.
#' @param ashman,min_share Thresholds for two zones.
#' @param seed_zones Base seed of the k-means used to look for two zones.
#' @param n_threads Threads for uwot; NULL leaves the default.
#' @return Object of class "secos": map (data frame with one row per term),
#'   replicates (aligned maps), D (distances), settings, timing (seconds by step).
secos_map <- function(E, D_cooc = NULL, lambda = 1, n_neighbors = 10, min_dist = 0.1, replicates = 100L,
                      seeds = seq_len(replicates), area_limit = 0.10, ashman = 4, min_share = 0.10,
                      seed_zones = 2026L, n_threads = NULL) {
  E <- as.matrix(E)
  terms <- rownames(E)
  n <- nrow(E)
  if (is.null(terms) || anyDuplicated(terms)) stop("E must have unique terms as row names")
  if (lambda < 0 || lambda > 1) stop("lambda must lie between 0 and 1")
  if (lambda < 1 && is.null(D_cooc)) stop("lambda < 1 needs D_cooc, from secos_cooccurrence()")
  if (n <= n_neighbors) stop("more terms than n_neighbors are needed")
  if (length(seeds) != replicates || replicates < 2) stop("one seed per replicate, and at least two replicates")
  for (p in c("uwot", "igraph")) if (!requireNamespace(p, quietly = TRUE)) stop("secos_map() needs the ", p, " package")
  # R's default generator whatever the session uses, so that a seed always gives the same
  # map (some packages leave another generator on); the session is restored on exit.
  kind <- RNGkind()
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit({
    suppressWarnings(RNGkind(kind[1], kind[2], kind[3]))
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv())) rm(".Random.seed", envir = globalenv())
  })
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  timing <- numeric()
  tic <- function(name, expr) { t0 <- proc.time()[["elapsed"]]; value <- force(expr)
    timing[[name]] <<- proc.time()[["elapsed"]] - t0; value }

  # Distances and exact neighbors (uwot counts the term itself among them)
  nn <- tic("neighbours", {
    Z <- E / sqrt(rowSums(E^2))
    D <- .secos_cosine(E)
    if (lambda < 1) {
      Dc <- D_cooc[terms, terms]
      D <- lambda * .secos_standardise(D) + (1 - lambda) * .secos_standardise(Dc)
      D <- D - min(D[lower.tri(D)]) + 1e-6
      diag(D) <- 0
    }
    N <- .secos_neighbours(D, n_neighbors - 1)
    idx <- cbind(seq_len(n), N)
    dimnames(idx) <- NULL
    storage.mode(idx) <- "integer"
    list(idx = idx, dist = t(vapply(seq_len(n), function(i) unname(D[i, idx[i, ]]), numeric(n_neighbors))))
  })
  run <- function(s, ...) {
    set.seed(s)
    Y <- uwot::umap(Z, n_neighbors = n_neighbors, min_dist = min_dist, metric = "cosine", n_components = 2,
                    nn_method = nn, n_threads = n_threads, ...)
    rownames(Y) <- terms
    Y
  }

  # Fixed start: spectral if the graph is connected. Otherwise uwot would fall back on
  # the PCA of the embeddings, so with lambda < 1 the start comes from the combined distance.
  graph <- igraph::graph_from_edgelist(cbind(rep(seq_len(n), n_neighbors - 1), as.vector(nn$idx[, -1])), directed = FALSE)
  components <- igraph::components(graph)$no
  start <- tic("start", {
    if (components == 1) run(seeds[1], n_epochs = 0) else {
      Y <- stats::cmdscale(D, k = 2)
      10 * Y / max(abs(Y))
    }
  })

  maps <- tic("replicates", lapply(seeds, run, init = start))
  aligned <- tic("alignment", .secos_procrustes(maps))
  M <- tic("median", .secos_median(aligned))

  # Ellipses at 95% and search for two zones
  box <- prod(apply(M, 2, function(x) diff(range(x))))
  A <- simplify2array(aligned)
  ellipses <- tic("ellipses", do.call(rbind, lapply(seq_len(n), function(i) {
    e <- eigen(stats::cov(t(A[i, , ])), symmetric = TRUE)
    q <- stats::qchisq(0.95, 2)
    data.frame(axis1 = sqrt(q * max(e$values[1], 0)), axis2 = sqrt(q * max(e$values[2], 0)),
               angle = atan2(e$vectors[2, 1], e$vectors[1, 1]))
  })))
  zones <- tic("zones", t(vapply(seq_len(n), function(i) .secos_zones(t(A[i, , ]), seed_zones + i), c(D = 0, share = 0))))
  map <- data.frame(term = terms, x = M[, 1], y = M[, 2], ellipses, stringsAsFactors = FALSE, row.names = NULL)
  map$area_share <- pi * map$axis1 * map$axis2 / box
  map$ashman_D <- zones[, "D"]
  map$smaller_zone <- zones[, "share"]
  map$two_zones <- !is.na(map$ashman_D) & map$ashman_D >= ashman & map$smaller_zone >= min_share
  map$uncertain <- map$area_share >= area_limit | map$two_zones
  x <- structure(list(map = map, replicates = aligned, D = D,
                      settings = list(lambda = lambda, n_neighbors = n_neighbors, min_dist = min_dist,
                                      replicates = replicates, seeds = seeds, components = components),
                      timing = timing), class = "secos")
  map$stability <- tic("stability", secos_stability(x)$per_term)
  x$map <- map
  x$timing <- timing
  x
}

# ---- Evaluation -----------------------------------------------------------------

#' Neighborhood stability
#'
#' For each term, the Jaccard index between its neighbors on the map, averaged over
#' all the pairs of replicates.
#'
#' @param x Object from secos_map(), or a list of matrices with two columns.
#' @param k Number of neighbors.
#' @param replicates Replicates to use, by default all.
#' @return List: per_term, mean, median, chance (expected value with random neighbors).
secos_stability <- function(x, k = 10, replicates = NULL) {
  Y <- if (inherits(x, "secos")) x$replicates else x
  if (!is.null(replicates)) Y <- Y[seq_len(replicates)]
  R <- length(Y)
  if (R < 2) stop("at least two replicates are needed")
  n <- nrow(Y[[1]])
  keys <- lapply(Y, function(y) as.vector(.secos_neighbours(as.matrix(stats::dist(y)), k)) + rep((seq_len(n) - 1) * n, k))
  total <- numeric(n)
  for (a in seq_len(R - 1)) for (b in (a + 1):R) {
    shared <- rowSums(matrix(keys[[a]] %in% keys[[b]], n, k))
    total <- total + shared / (2 * k - shared)
  }
  per_term <- total / (R * (R - 1) / 2)
  names(per_term) <- rownames(Y[[1]])
  i <- 0:k
  list(per_term = per_term, mean = mean(per_term), median = stats::median(per_term),
       chance = sum(stats::dhyper(i, k, n - 1 - k, k) * i / (2 * k - i)))
}

#' Fidelity of the map to the original space
#'
#' Trustworthiness penalizes neighbors that exist only on the map, continuity those
#' that the map loses (Venna and Kaski 2001). Also reports the share of neighbors kept
#' and the Spearman correlation between all the distances. The value of k should
#' differ from the n_neighbors of UMAP.
#'
#' @param x Object from secos_map().
#' @param k Number of neighbors, below half of the terms.
#' @return Named vector.
secos_fidelity <- function(x, k = 15) {
  Dh <- x$D
  n <- nrow(Dh)
  if (k >= n / 2) stop("k must be below half of the terms")
  Dl <- as.matrix(stats::dist(as.matrix(x$map[, c("x", "y")])))
  Rh <- .secos_ranks(Dh); Rl <- .secos_ranks(Dl)
  Nh <- .secos_neighbours(Dh, k); Nl <- .secos_neighbours(Dl, k)
  trust <- cont <- 0
  for (i in seq_len(n)) {
    trust <- trust + sum(Rh[i, setdiff(Nl[i, ], Nh[i, ])] - k)
    cont <- cont + sum(Rl[i, setdiff(Nh[i, ], Nl[i, ])] - k)
  }
  g <- 2 / (n * k * (2 * n - 3 * k - 1))
  low <- lower.tri(Dh)
  c(trustworthiness = 1 - g * trust, continuity = 1 - g * cont,
    neighbours_kept = mean(vapply(seq_len(n), function(i) length(intersect(Nh[i, ], Nl[i, ])) / k, 0)),
    spearman = stats::cor(rank(Dh[low]), rank(Dl[low])))
}

# ---- Whole procedure --------------------------------------------------------------

#' SeCoS map from a bibliometrix data frame
#'
#' @param M Data frame from bibliometrix::convert2df().
#' @param field Column with the keywords, by default the author keywords.
#' @param min_docs Minimum number of documents in which a term appears.
#' @param remove_terms Terms to leave out, such as those used in the query.
#' @param terms Terms to map, instead of those chosen with min_docs.
#' @param E Matrix of embeddings with the terms as row names, instead of computing them.
#' @param model,revision,layer,pooling See secos_embed().
#' @param lambda,n_neighbors,min_dist,replicates,seeds See secos_map().
#' @param sep Separator between keywords.
#' @param ... Other arguments of secos_map().
#' @return Object of class "secos". The print method gives a summary and the plot
#'   method draws the map.
secos <- function(M, field = "DE", min_docs = 50, remove_terms = NULL, terms = NULL, E = NULL,
                  model = SECOS_MODEL, revision = if (identical(model, SECOS_MODEL)) SECOS_REVISION else NULL,
                  layer = -1L, pooling = c("all", "term"), lambda = 1, n_neighbors = 10, min_dist = 0.1,
                  replicates = 100L, seeds = seq_len(replicates), sep = ";", ...) {
  pooling <- match.arg(pooling)
  frequency <- secos_terms(M, field, 1, remove_terms, sep)
  if (is.null(terms)) terms <- frequency$term[frequency$documents >= min_docs]
  if (length(terms) == 0) stop("no term reaches min_docs = ", min_docs)
  if (is.null(E)) {
    E <- secos_embed(terms, model, revision, layer, pooling)
  } else {
    if (!all(terms %in% rownames(E))) stop("E has no row for some of the terms")
    E <- E[terms, , drop = FALSE]
  }
  co <- if (lambda < 1) secos_cooccurrence(M, terms, field, sep)$D
  x <- secos_map(E, D_cooc = co, lambda = lambda, n_neighbors = n_neighbors, min_dist = min_dist,
                 replicates = replicates, seeds = seeds, ...)
  x$map$documents <- frequency$documents[match(x$map$term, frequency$term)]
  x$settings$model <- model
  x$settings$revision <- revision
  x
}

#' @export
print.secos <- function(x, ...) {
  s <- x$settings
  cat(sprintf("SeCoS map: %d terms, %d replicates, lambda = %g, n_neighbors = %d, min_dist = %g\n",
              nrow(x$map), s$replicates, s$lambda, s$n_neighbors, s$min_dist))
  if (!is.null(s$model)) cat("Model:", s$model, "\n")
  cat(sprintf("Stability: mean %.3f, median %.3f\n", mean(x$map$stability), stats::median(x$map$stability)))
  cat(sprintf("Terms with an uncertain position: %d (not to be interpreted by position)\n", sum(x$map$uncertain)))
  if (s$components > 1) cat("The graph of neighbors has", s$components, "components: the start is not spectral\n")
  cat(sprintf("Elapsed: %.1f seconds\n", sum(x$timing)))
  invisible(x)
}

#' Plot of the consensus map
#'
#' Hollow markers and dashed ellipses mark the terms with an uncertain position. UMAP
#' does not preserve global distances: the separation between distant regions carries
#' no interpretation.
#'
#' @param x Object from secos_map() or secos().
#' @param labels Number of terms to label, the most frequent if frequencies are known.
#' @param col Color of the markers, one value or one per term.
#' @param ... Other arguments of plot().
#' @export
plot.secos <- function(x, labels = 20, col = "#2a78d6", ...) {
  m <- x$map
  col <- rep_len(col, nrow(m))
  plot(m$x, m$y, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "", ...)
  t <- seq(0, 2 * pi, length.out = 40)
  for (i in seq_len(nrow(m))) {
    a <- m$axis1[i] * cos(t); b <- m$axis2[i] * sin(t); w <- m$angle[i]
    graphics::polygon(m$x[i] + a * cos(w) - b * sin(w), m$y[i] + a * sin(w) + b * cos(w),
                      border = if (m$uncertain[i]) "grey55" else NA, lty = 2,
                      col = if (m$uncertain[i]) NA else grDevices::adjustcolor(col[i], 0.07))
  }
  size <- if (is.null(m$documents)) rep(0.7, nrow(m)) else 0.4 + sqrt(m$documents / max(m$documents)) * 1.4
  graphics::points(m$x, m$y, pch = ifelse(m$uncertain, 1, 19), col = col, cex = size)
  ordered <- if (is.null(m$documents)) seq_len(nrow(m)) else order(-m$documents)
  chosen <- ordered[seq_len(min(labels, nrow(m)))]
  graphics::text(m$x[chosen], m$y[chosen], m$term[chosen], cex = 0.65, pos = 3, offset = 0.3,
                 font = ifelse(m$uncertain[chosen], 3, 1))
  invisible(x)
}
