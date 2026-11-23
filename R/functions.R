SEED <- 2026
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
PKG <- c("bibliometrix", "uwot", "Rtsne", "ca", "ggplot2", "poweRlaw")
mancanti <- PKG[!vapply(PKG, requireNamespace, logical(1), quietly = TRUE)]
if (length(mancanti) > 0) stop("Install the packages: ", paste(mancanti, collapse = ", "))
invisible(lapply(PKG, library, character.only = TRUE))
for (d in c("results", "wos")) if (!dir.exists(d)) dir.create(d)
sessione <- function(nome) writeLines(capture.output(sessionInfo()), sprintf("results/sessionInfo_%s.txt", nome))
leggi_csv <- function(f, ...) read.csv(f, stringsAsFactors = FALSE, encoding = "UTF-8", ...)
scrivi_csv <- function(x, f) write.csv(x, f, row.names = FALSE, fileEncoding = "UTF-8")
norm_l2 <- function(X) {
  X <- as.matrix(X)
  X / sqrt(rowSums(X^2))
}
dist_coseno <- function(X) {
  Z <- norm_l2(X)
  D <- 1 - tcrossprod(Z)
  D[D < 0] <- 0
  diag(D) <- 0
  dimnames(D) <- list(rownames(X), rownames(X))
  D
}
dist_euclidea <- function(Y) {
  D <- as.matrix(stats::dist(as.matrix(Y)))
  dimnames(D) <- list(rownames(Y), rownames(Y))
  D
}
vicini_knn <- function(D, k) {
  D <- as.matrix(D); n <- nrow(D)
  stopifnot(k >= 1, k < n)
  diag(D) <- Inf
  N <- t(vapply(seq_len(n), function(i) order(D[i, ], seq_len(n))[seq_len(k)],
                integer(k)))
  if (k == 1) N <- matrix(N, ncol = 1)
  rownames(N) <- rownames(D)
  N
}
ranghi_vicini <- function(D) {
  D <- as.matrix(D); n <- nrow(D)
  R <- matrix(0L, n, n)
  for (i in seq_len(n)) {
    d <- D[i, ]; d[i] <- -Inf
    o <- order(d, seq_len(n))
    R[i, o] <- seq_len(n) - 1L
  }
  R
}
jaccard_vicini <- function(N1, N2) {
  vapply(seq_len(nrow(N1)), function(i) {
    a <- N1[i, ]; b <- N2[i, ]
    length(intersect(a, b)) / length(union(a, b))
  }, numeric(1))
}
stabilita_vicini <- function(lista_Y, k) {
  R <- length(lista_Y); stopifnot(R >= 2)
  N <- lapply(lista_Y, function(Y) vicini_knn(dist_euclidea(Y), k))
  coppie <- utils::combn(R, 2)
  J <- vapply(seq_len(ncol(coppie)),
              function(c) jaccard_vicini(N[[coppie[1, c]]], N[[coppie[2, c]]]),
              numeric(nrow(N[[1]])))
  J <- matrix(J, nrow = nrow(N[[1]]))
  S <- rowMeans(J)
  names(S) <- rownames(lista_Y[[1]])
  list(per_termine = S, media = mean(S), mediana = stats::median(S),
       R = R, k = k)
}
jaccard_casuale <- function(n, k) {
  x <- 0:k
  p <- stats::dhyper(x, m = k, n = n - 1 - k, k = k)
  sum(p * x / (2 * k - x))
}
trust_cont <- function(D_alto, D_basso, k) {
  n <- nrow(D_alto)
  if (k >= n / 2) stop("trust_cont needs k < n/2")
  Ra <- ranghi_vicini(D_alto); Rb <- ranghi_vicini(D_basso)
  Na <- vicini_knn(D_alto, k); Nb <- vicini_knn(D_basso, k)
  pen_t <- pen_c <- 0
  for (i in seq_len(n)) {
    U <- setdiff(Nb[i, ], Na[i, ])
    V <- setdiff(Na[i, ], Nb[i, ])
    pen_t <- pen_t + sum(Ra[i, U] - k)
    pen_c <- pen_c + sum(Rb[i, V] - k)
  }
  g <- 2 / (n * k * (2 * n - 3 * k - 1))
  c(trustworthiness = 1 - g * pen_t, continuity = 1 - g * pen_c)
}
conservazione_knn <- function(D_alto, D_basso, k, per_termine = FALSE) {
  Na <- vicini_knn(D_alto, k); Nb <- vicini_knn(D_basso, k)
  q <- vapply(seq_len(nrow(Na)),
              function(i) length(intersect(Na[i, ], Nb[i, ])) / k, numeric(1))
  if (per_termine) { names(q) <- rownames(D_alto); q } else mean(q)
}
correlazione_distanze <- function(D_alto, D_basso, permutazioni = 0, seed = 1) {
  lt <- lower.tri(D_alto)
  ra <- rank(D_alto[lt]); rb <- rank(D_basso[lt])
  r <- stats::cor(ra, rb)
  p <- NA_real_
  if (permutazioni > 0) {
    Ra <- matrix(0, nrow(D_alto), ncol(D_alto)); Ra[lt] <- ra; Ra <- Ra + t(Ra)
    set.seed(seed)
    rp <- vapply(seq_len(permutazioni), function(b) {
      o <- sample(nrow(Ra)); stats::cor(Ra[o, o][lt], rb)
    }, numeric(1))
    p <- (sum(rp >= r) + 1) / (permutazioni + 1)
  }
  c(spearman = r, p_mantel = p)
}
ari <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  s_ij <- sum(c2(tab)); s_a <- sum(c2(rowSums(tab))); s_b <- sum(c2(colSums(tab)))
  atteso <- s_a * s_b / c2(sum(tab))
  massimo <- (s_a + s_b) / 2
  if (massimo == atteso) return(1)
  (s_ij - atteso) / (massimo - atteso)
}
nmi <- function(a, b, variante = c("media", "sqrt", "max")) {
  variante <- match.arg(variante)
  p <- table(a, b) / length(a)
  pa <- rowSums(p); pb <- colSums(p)
  H <- function(q) { q <- q[q > 0]; -sum(q * log(q)) }
  nz <- p > 0
  I <- sum(p[nz] * log(p[nz] / outer(pa, pb)[nz]))
  Ha <- H(pa); Hb <- H(pb)
  if (Ha == 0 && Hb == 0) return(1)
  den <- switch(variante, media = (Ha + Hb) / 2, sqrt = sqrt(Ha * Hb), max = max(Ha, Hb))
  if (den == 0) return(0)
  I / den
}
sim_albero <- function(alberi_a, alberi_b) {
  if (length(alberi_a) == 0 || length(alberi_b) == 0) return(NA_real_)
  pa <- strsplit(alberi_a, ".", fixed = TRUE)
  pb <- strsplit(alberi_b, ".", fixed = TRUE)
  best <- 0
  for (x in pa) for (y in pb) {
    m <- min(length(x), length(y))
    lcp <- if (m == 0) 0 else { d <- which(x[seq_len(m)] != y[seq_len(m)])
                                 if (length(d) == 0) m else d[1] - 1 }
    best <- max(best, 2 * lcp / (length(x) + length(y)))
  }
  best
}
accordo_gerarchia <- function(D, alberi, k, permutazioni = 999, seed = 1) {
  alberi <- alberi[vapply(alberi, length, integer(1)) > 0]
  t_ok <- intersect(rownames(D), names(alberi))
  m <- length(t_ok)
  if (m <= k) stop("accordo_gerarchia: mapped terms (", m, ") not more than k")
  S <- matrix(0, m, m)
  for (i in seq_len(m - 1)) for (j in (i + 1):m)
    S[i, j] <- S[j, i] <- sim_albero(alberi[[t_ok[i]]], alberi[[t_ok[j]]])
  N <- vicini_knn(D[t_ok, t_ok], k)
  stat <- function(SS) mean(vapply(seq_len(m), function(i) mean(SS[i, N[i, ]]), numeric(1)))
  oss <- stat(S)
  per_termine <- vapply(seq_len(m), function(i) mean(S[i, N[i, ]]), numeric(1))
  names(per_termine) <- t_ok
  set.seed(seed)
  perm <- vapply(seq_len(permutazioni), function(b) { o <- sample(m); stat(S[o, o]) },
                 numeric(1))
  list(statistica = oss, media_permutazioni = mean(perm),
       p = (sum(perm >= oss) + 1) / (permutazioni + 1),
       termini_mappati = m, k = k, per_termine = per_termine)
}
valuta_proiezione <- function(D_alto, Y, k_val = 15, permutazioni = 0) {
  Db <- dist_euclidea(Y)
  tc <- trust_cont(D_alto, Db, k_val)
  cd <- correlazione_distanze(D_alto, Db, permutazioni)
  data.frame(T = nrow(D_alto), k_val = k_val,
             trustworthiness = unname(tc[1]), continuity = unname(tc[2]),
             conservazione_knn = conservazione_knn(D_alto, Db, k_val),
             spearman = unname(cd[1]), p_mantel = unname(cd[2]))
}
wilson <- function(x, n) {
  p <- x / n; z <- qnorm(0.975)
  c0 <- (p + z^2 / (2 * n)) / (1 + z^2 / n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(p, c0 - h, c0 + h)
}
matrice_sim <- function(alb, termini) {
  tt <- intersect(termini, names(alb)[lengths(alb) > 0]); m <- length(tt)
  S <- matrix(0, m, m, dimnames = list(tt, tt))
  for (i in seq_len(m - 1)) for (j in (i + 1):m) S[i, j] <- S[j, i] <- sim_albero(alb[[tt[i]]], alb[[tt[j]]])
  S
}
acc_S <- function(Dx, S, k) {
  tt <- rownames(S); N <- vicini_knn(Dx[tt, tt], k)
  setNames(vapply(seq_along(tt), function(i) mean(S[i, N[i, ]]), 0), tt)
}
ic_diff <- function(a, b, seme = SEED, B = 2000) {
  n <- intersect(names(a), names(b)); d <- b[n] - a[n]
  set.seed(seme)
  bs <- vapply(seq_len(B), function(i) mean(sample(d, replace = TRUE)), 0)
  c(differenza = mean(d), inf = unname(quantile(bs, 0.025)), sup = unname(quantile(bs, 0.975)))
}
gpa <- function(lista, iter = 10) {
  Y <- lapply(lista, function(y) sweep(y, 2, colMeans(y))); rif <- Y[[1]]
  for (it in seq_len(iter)) {
    Y <- lapply(Y, function(y) { s <- svd(crossprod(y, rif)); y %*% (s$u %*% t(s$v)) })
    nuovo <- Reduce(`+`, Y) / length(Y)
    if (max(abs(nuovo - rif)) < 1e-9) break
    rif <- nuovo
  }
  Y
}
mediana <- function(Y) { A <- simplify2array(Y); M <- apply(A, c(1, 2), median); dimnames(M) <- list(rownames(Y[[1]]), NULL); M }
quota_vicini <- function(l, k) {
  N <- lapply(l, function(Y) vicini_knn(dist_euclidea(Y), k)); cp <- combn(length(N), 2)
  mean(vapply(seq_len(ncol(cp)), function(j) { a <- N[[cp[1, j]]]; b <- N[[cp[2, j]]]
    mean(vapply(seq_len(nrow(a)), function(i) length(intersect(a[i, ], b[i, ])) / k, 0)) }, 0))
}
mesh_alberi <- function(f = "data/mesh_mapping.csv") {
  m <- leggi_csv(f, na.strings = "")
  a <- function(col) { k <- !is.na(m[[col]]); setNames(lapply(strsplit(m[[col]][k], "; ", fixed = TRUE), trimws), m$term[k]) }
  list(finale = a("tree_numbers"), autore = a("tree_numbers_before_review"), automatica = a("tree_numbers_automatic"))
}
CODING_EN <- c(finale = "final", autore = "before_review", automatica = "automatic")
matrice_emb <- function(f) { E <- as.matrix(read.csv(f, row.names = 1, check.names = FALSE, encoding = "UTF-8")); dimnames(E) <- list(rownames(E), NULL); E }
COL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INCHIOSTRO <- "#2b2b28"; GRIGIO <- "#8c8b83"; CHIARO <- "#d4d3cb"
tema <- theme_minimal(base_size = 8.5) + theme(
  text = element_text(colour = INCHIOSTRO), plot.title = element_text(size = 9, face = "bold"), plot.subtitle = element_text(size = 7, colour = "#6b6a63"),
  plot.background = element_rect(fill = "white", colour = NA), panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "#ecebe5", linewidth = 0.25), legend.box = "vertical", legend.spacing.y = unit(0.5, "mm"), strip.text = element_text(size = 8, face = "bold", hjust = 0), legend.position = "bottom",
  legend.title = element_blank(), legend.key.size = unit(3.2, "mm"), legend.text = element_text(size = 7.5))
salva <- function(p, nome, w, h) ggsave(file.path("results", paste0(nome, ".pdf")), p, width = w, height = h, units = "mm", device = grDevices::pdf, bg = "white")
