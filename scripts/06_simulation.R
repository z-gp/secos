source("R/functions.R")
R_REP           <- 10
G               <- 6
D_EMB           <- 768
DOC_PER_TERMINE <- 125
Q_DISTURBO      <- 20
KW_MEDIE        <- 4
MIN_DOC         <- 5
RUMORE  <- list(basso = c(sigma = 0.6, eps = 0.2), alto = c(sigma = 1.0, eps = 0.4))
SOVRAPP <- c(bassa = 0.3, alta = 0.5)
K_UMAP  <- 10
MIN_DIST <- 0.1
K_PUR   <- 10
K_VAL   <- 15
GRIGLIA <- expand.grid(T = c(100, 300), rumore = names(RUMORE),
                       sovrapposizione = names(SOVRAPP), stringsAsFactors = FALSE)
genera_dati <- function(T, sigma, eps, omega, seed) {
  set.seed(seed)
  gr <- sample(rep(seq_len(G), length.out = T))
  dirichlet <- function(n) { U <- matrix(rgamma(n * G, 1), n, G); U / rowSums(U) }
  Pi <- (1 - omega) * diag(G)[gr, ] + omega * dirichlet(T)
  M  <- matrix(rnorm(G * D_EMB), G, D_EMB); M <- M / sqrt(rowSums(M^2))
  W  <- matrix(rnorm(Q_DISTURBO * D_EMB), Q_DISTURBO, D_EMB); W <- W / sqrt(rowSums(W^2))
  E  <- Pi %*% M + sigma * matrix(rnorm(T * Q_DISTURBO), T, Q_DISTURBO) %*% W / sqrt(Q_DISTURBO)
  f  <- rlnorm(T, 0, 0.75)
  nd <- DOC_PER_TERMINE * T
  gd <- sample(seq_len(G), nd, replace = TRUE)
  P  <- sweep(t(Pi)[gd, , drop = FALSE], 2, f, "*")
  B  <- matrix(0L, nd, T)
  for (i in seq_len(nd)) {
    L  <- min(T, 1 + rpois(1, KW_MEDIE - 1))
    nr <- rbinom(1, L, eps)
    s1 <- if (L > nr) sample.int(T, L - nr, prob = P[i, ]) else integer(0)
    altri <- setdiff(seq_len(T), s1)
    s2 <- if (nr > 0) altri[sample.int(length(altri), nr, prob = f[altri])] else integer(0)
    B[i, c(s1, s2)] <- 1L
  }
  tieni <- colSums(B) >= MIN_DOC
  B <- B[, tieni, drop = FALSE]; B <- B[rowSums(B) > 0, , drop = FALSE]
  nomi <- sprintf("t%03d", which(tieni))
  colnames(B) <- nomi
  E <- E[tieni, , drop = FALSE]; rownames(E) <- nomi
  list(E = E, B = B, gruppo = gr[tieni], Pi = Pi[tieni, , drop = FALSE])
}
ca_termini <- function(B) {
  B <- as.matrix(B); storage.mode(B) <- "double"
  N <- sum(B); L <- rowSums(B); cm <- colSums(B) / N
  S <- (crossprod(B / sqrt(L)) / N - tcrossprod(cm)) / sqrt(tcrossprod(cm))
  e <- eigen(S, symmetric = TRUE)
  X <- e$vectors[, 1:2] / sqrt(cm)
  dimnames(X) <- list(colnames(B), c("Dim1", "Dim2")); X
}
mca_termini <- function(B) {
  B <- as.matrix(B); storage.mode(B) <- "double"
  n <- nrow(B); T <- ncol(B); C <- crossprod(B); f <- diag(C)
  i0 <- 2 * seq_len(T) - 1; i1 <- 2 * seq_len(T)
  Bu <- matrix(0, 2 * T, 2 * T)
  Bu[i1, i1] <- C
  Bu[i1, i0] <- f - C
  Bu[i0, i1] <- t(f - C)
  Bu[i0, i0] <- n - outer(f, f, "+") + C
  diag(Bu)[i0] <- n - f; Bu[cbind(i0, i1)] <- 0; Bu[cbind(i1, i0)] <- 0
  P <- Bu / sum(Bu); cm <- colSums(P)
  e <- eigen((P - tcrossprod(cm)) / sqrt(tcrossprod(cm)), symmetric = TRUE)
  X <- sweep(e$vectors[, 1:2] / sqrt(cm), 2, sqrt(e$values[1:2]), "*")[i1, ]
  dimnames(X) <- list(colnames(B), c("Dim1", "Dim2")); X
}
proietta <- function(d, seed) {
  Z <- norm_l2(d$E); T <- nrow(Z)
  set.seed(seed)
  secos <- uwot::umap(Z, n_neighbors = K_UMAP, min_dist = MIN_DIST, metric = "cosine",
                      n_components = 2)
  set.seed(seed)
  tsne <- Rtsne::Rtsne(Z, dims = 2, perplexity = min(30, floor((T - 1) / 3)),
                       check_duplicates = FALSE)$Y
  list("SeCoS (UMAP)" = secos, "t-SNE" = tsne,
       "PCA" = stats::prcomp(Z)$x[, 1:2],
       "MDS" = stats::cmdscale(stats::as.dist(dist_coseno(d$E)), k = 2),
       "CA"  = ca_termini(d$B),
       "MCA" = mca_termini(d$B))
}
misura <- function(X, Dist, d, seed) {
  set.seed(seed)
  km <- stats::kmeans(X, centers = G, nstart = 25, iter.max = 50)$cluster
  N <- vicini_knn(Dist, K_PUR)
  pur <- mean(vapply(seq_len(nrow(N)), function(i) mean(d$gruppo[N[i, ]] == d$gruppo[i]),
                     numeric(1)))
  tc <- trust_cont(dist_euclidea(d$Pi), Dist, K_VAL)
  c(ARI = ari(km, d$gruppo), purezza = pur, trust = unname(tc[1]), continuity = unname(tc[2]))
}
cooc_assoc <- function(B) {
  C <- crossprod(B); fr <- diag(C)
  A <- C / outer(fr, fr); diag(A) <- 0; A
}
FONTE <- c("SeCoS (UMAP)" = "embedding", "t-SNE" = "embedding", "PCA" = "embedding",
           "MDS" = "embedding", "CA" = "documenti", "MCA" = "documenti")
righe <- list()
for (c_i in seq_len(nrow(GRIGLIA))) {
  g <- GRIGLIA[c_i, ]
  for (r in seq_len(R_REP)) {
    seed <- SEED + 1000 * c_i + r
    d <- genera_dati(g$T, RUMORE[[g$rumore]]["sigma"], RUMORE[[g$rumore]]["eps"],
                     SOVRAPP[[g$sovrapposizione]], seed)
    Ys <- proietta(d, seed)
    m <- t(vapply(Ys, function(Y) misura(Y, dist_euclidea(Y), d, seed), numeric(4)))
    A <- cooc_assoc(d$B)
    rif <- rbind("Embeddings (original space)" = misura(norm_l2(d$E), dist_coseno(d$E), d, seed),
                 "Co-occurrence (original space)" = misura(A, dist_coseno(A), d, seed))
    m <- rbind(m, rif)
    righe[[length(righe) + 1]] <- data.frame(
      condizione = c_i, T = g$T, rumore = g$rumore, sovrapposizione = g$sovrapposizione,
      replica = r, seed = seed, T_eff = nrow(d$E), n_doc = nrow(d$B),
      metodo = rownames(m),
      fonte = c(FONTE[rownames(m)[1:6]], "embedding", "documenti"),
      tipo = c(rep("proiezione 2D", 6), rep("riferimento", 2)),
      m, row.names = NULL, check.names = FALSE)
  }
}
ris <- do.call(rbind, righe)
misure_nomi <- c("ARI", "purezza", "trust", "continuity")
chiave <- c("condizione", "T", "rumore", "sovrapposizione", "metodo", "fonte", "tipo")
media <- aggregate(ris[, misure_nomi], ris[, chiave], mean)
sdev  <- aggregate(ris[, misure_nomi], ris[, chiave], sd)
tab <- merge(media, sdev, by = chiave, suffixes = c("_media", "_sd"))
ordine_metodi <- c(names(FONTE), "Embeddings (original space)", "Co-occurrence (original space)")
tab <- tab[order(tab$condizione, match(tab$metodo, ordine_metodi)), ]
etich_cond <- function(rum, sov) {
  lab <- sprintf("Noise: %s\nOverlap: %s",
                 ifelse(rum == "basso", "low", "high"), ifelse(sov == "bassa", "low", "high"))
  livelli <- sprintf("Noise: %s\nOverlap: %s", c("low", "low", "high", "high"),
                     c("low", "high", "low", "high"))
  factor(lab, levels = livelli)
}
figura <- function(misura_nome, titolo_y, file) {
  x <- tab[tab$tipo == "proiezione 2D", ]
  x$metodo <- factor(x$metodo, levels = names(FONTE))
  x$source <- factor(ifelse(x$fonte == "embedding", "Embeddings", "Documents"),
                     levels = c("Embeddings", "Documents"))
  x$cond <- etich_cond(x$rumore, x$sovrapposizione)
  x$Tlab <- paste("T =", x$T)
  rf <- tab[tab$tipo == "riferimento", ]
  rf$source <- factor(ifelse(rf$fonte == "embedding", "Embeddings", "Documents"),
                      levels = c("Embeddings", "Documents"))
  rf$cond <- etich_cond(rf$rumore, rf$sovrapposizione); rf$Tlab <- paste("T =", rf$T)
  x$m <- x[[paste0(misura_nome, "_media")]]; x$s <- x[[paste0(misura_nome, "_sd")]]
  rf$m <- rf[[paste0(misura_nome, "_media")]]
  col <- c(Embeddings = "#2a78d6", Documents = "#eb6834")
  p <- ggplot2::ggplot(x, ggplot2::aes(metodo, m, colour = source, shape = source)) +
    ggplot2::geom_hline(data = rf, ggplot2::aes(yintercept = m, colour = source),
                        linetype = "dashed", linewidth = 0.4, show.legend = FALSE) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = pmax(m - s, 0), ymax = pmin(m + s, 1)),
                           width = 0.25, linewidth = 0.5) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::facet_grid(Tlab ~ cond) +
    ggplot2::scale_colour_manual(values = col, name = "Input data") +
    ggplot2::scale_shape_manual(values = c(Embeddings = 16, Documents = 17), name = "Input data") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    ggplot2::labs(x = NULL, y = titolo_y,
                  caption = paste0("Points: mean over ", R_REP, " replicates; bars: \u00b1 1 SD. ",
                                   "Dashed lines: same measure in the original space of each input ",
                                   "(no projection).")) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   legend.position = "bottom", plot.caption = ggplot2::element_text(hjust = 0))
  ggplot2::ggsave(paste0(file, ".pdf"), p, width = 8, height = 5.2)
}
inglese <- function(x) {
  x$rumore <- unname(c(basso = "low", alto = "high")[as.character(x$rumore)]); x$sovrapposizione <- unname(c(bassa = "low", alta = "high")[as.character(x$sovrapposizione)])
  x$fonte <- unname(c(embedding = "embeddings", documenti = "documents")[as.character(x$fonte)]); x$tipo <- unname(c("proiezione 2D" = "2D projection", riferimento = "reference")[as.character(x$tipo)])
  stopifnot(!anyNA(x[c("rumore", "sovrapposizione", "fonte", "tipo")]))
  nomi <- c(condizione = "condition", rumore = "noise", sovrapposizione = "overlap", replica = "replicate", metodo = "method", fonte = "input", tipo = "type")
  j <- names(x) %in% names(nomi); names(x)[j] <- nomi[names(x)[j]]
  names(x) <- sub("_media$", "_mean", sub("^purezza", "purity", names(x))); x
}
scrivi_csv(inglese(ris), "results/simulation_replicates.csv")
scrivi_csv(inglese(tab), "results/simulation_summary.csv")
figura("ARI", "Adjusted Rand index (k-means vs true groups)", "results/fig5_simulation_ari")
figura("purezza", sprintf("Neighborhood purity (k = %d)", K_PUR), "results/figS7_simulation_purity")
sessione("06")
