source("R/functions.R")
SEMI <- 1:20; K_B <- 10; MD_B <- 0.1; K_V <- 15; R_CONS <- 100; R_LOTTI <- c(5, 10, 20, 50); SOGLIE <- c(30, 75, 100)
VARIANTI <- list("k = 5" = list(k = 5), "k = 15" = list(k = 15), "k = 30" = list(k = 30), "k = 50" = list(k = 50),
                 "min_dist = 0.01" = list(md = 0.01), "min_dist = 0.05" = list(md = 0.05),
                 "min_dist = 0.25" = list(md = 0.25), "min_dist = 0.5" = list(md = 0.5),
                 "Euclidean" = list(metrica = "euclidean"), "correlation" = list(metrica = "correlation"))
TERMINI <- leggi_csv("results/terms.csv")$term; n_T <- length(TERMINI)
ENC <- readRDS("results/embeddings.rds"); um14 <- readRDS("results/umap_runs.rds")
freq <- leggi_csv("results/term_frequencies.csv")
COD <- mesh_alberi()
E0 <- ENC$sbert[TERMINI, ]; Z <- norm_l2(E0); D0 <- dist_coseno(E0)
SIM <- matrice_sim(COD$finale, TERMINI)
misura <- function(l, Dh) {
  tt <- rownames(l[[1]]); st10 <- stabilita_vicini(l, 10)
  f <- do.call(rbind, lapply(l, function(Y) valuta_proiezione(Dh[tt, tt], Y, k_val = K_V)))
  S <- SIM[rownames(SIM) %in% tt, rownames(SIM) %in% tt]
  c(T = length(tt), stability = st10$media, stability_k15 = stabilita_vicini(l, K_V)$media, shared = quota_vicini(l, 10),
    trust = mean(f$trustworthiness), continuity = mean(f$continuity), neighbours = mean(f$conservazione_knn), spearman = mean(f$spearman),
    mesh = mean(rowMeans(sapply(l, function(Y) acc_S(dist_euclidea(Y), S, 10)))))
}
um_run <- function(X, s, k = K_B, md = MD_B, metrica = "cosine", ...) {
  set.seed(s)
  Y <- uwot::umap(X, n_neighbors = k, min_dist = md, metric = metrica, n_components = 2, ...)
  rownames(Y) <- rownames(X); Y
}
E_all <- rbind(E0, ENC$sbert_extra); D_all <- dist_coseno(E_all)
G <- list()
for (v in names(VARIANTI)) {
  a <- VARIANTI[[v]]; m <- if (is.null(a$metrica)) "cosine" else a$metrica
  X <- if (m == "euclidean") E0 else Z
  G[[v]] <- lapply(SEMI, function(s) um_run(X, s, k = if (is.null(a$k)) K_B else a$k, md = if (is.null(a$md)) MD_B else a$md, metrica = m))
}
for (s in SOGLIE) G[[sprintf("threshold %d", s)]] <- lapply(SEMI, function(r) um_run(norm_l2(E_all[freq$term[freq$documents >= s], ]), r))
righe <- list("base" = misura(um14$sbert, D0))
for (v in names(VARIANTI)) righe[[v]] <- misura(G[[v]], D0)
for (s in SOGLIE) righe[[sprintf("threshold %d", s)]] <- misura(G[[sprintf("threshold %d", s)]], D_all)
NOMI_E <- c(scibert = "SciBERT isolated", scibert_context = "SciBERT contextualized", pubmed = "PubMedBERT", fasttext = "fastText")
for (e in names(NOMI_E)) righe[[paste("encoder", NOMI_E[[e]])]] <- misura(um14[[e]], dist_coseno(ENC[[e]][TERMINI, ]))
acc_base <- function(l) {
  tt <- intersect(rownames(l[[1]]), TERMINI); k <- min(10, length(tt) - 1)
  mean(vapply(seq_along(SEMI), function(s) { A <- vicini_knn(dist_euclidea(l[[s]][tt, ]), k); B <- vicini_knn(dist_euclidea(um14$sbert[[s]][tt, ]), k)
    mean(vapply(seq_along(tt), function(i) length(intersect(A[i, ], B[i, ])) / k, 0)) }, 0))
}
tab <- as.data.frame(do.call(rbind, righe)); tab$setting <- names(righe)
tab$vs_base <- vapply(tab$setting, function(v) if (is.null(G[[v]])) NA_real_ else acc_base(G[[v]]), 0)
tsne <- lapply(SEMI, function(s) { set.seed(s)
  Y <- Rtsne::Rtsne(Z, dims = 2, perplexity = min(30, floor((n_T - 1) / 3)), check_duplicates = FALSE)$Y; rownames(Y) <- TERMINI; Y })
tab <- rbind(tab, data.frame(as.list(misura(tsne, D0)), setting = "t-SNE", vs_base = NA, check.names = FALSE))
st_sb <- stabilita_vicini(um14$sbert, 10); st_ts <- stabilita_vicini(tsne, 10)
ic_ts <- ic_diff(st_ts$per_termine, st_sb$per_termine)
scrivi_csv(tab, "results/table_S9_sensitivity.csv")
scrivi_csv(data.frame(term = TERMINI, stability_secos = st_sb$per_termine, stability_tsne = st_ts$per_termine), "results/stability_terms.csv")
scrivi_csv(data.frame(secos_mean = st_sb$media, secos_median = st_sb$mediana, share_above_0.5 = mean(st_sb$per_termine > 0.5),
                      tsne_mean = st_ts$media, tsne_median = st_ts$mediana, tsne_min = min(st_ts$per_termine),
                      difference = ic_ts[[1]], ci_low = ic_ts[[2]], ci_high = ic_ts[[3]], chance = jaccard_casuale(n_T, 10)), "results/stability_summary.csv")

Nx <- vicini_knn(D0, K_B - 1); idx <- cbind(seq_len(n_T), Nx); dimnames(idx) <- NULL; storage.mode(idx) <- "integer"
nn <- list(idx = idx, dist = t(vapply(seq_len(n_T), function(i) unname(D0[i, idx[i, ]]), numeric(K_B))))
Y_sp <- um_run(Z, 1, nn_method = nn, n_epochs = 0)
Y_pca <- stats::prcomp(Z)$x[, 1:2]
Y_pca <- sweep(sweep(Y_pca, 2, apply(Y_pca, 2, sd), "/"), 2, apply(Y_sp, 2, sd), "*")
rep_c <- lapply(seq_len(R_CONS), function(s) um_run(Z, s, nn_method = nn, init = Y_sp))
rep_p <- lapply(SEMI, function(s) um_run(Z, s, nn_method = nn, init = Y_pca))
RNGkind("L'Ecuyer-CMRG")
Y_sp2 <- um_run(Z, 1, nn_method = nn, n_epochs = 0)
rep_c2 <- lapply(seq_len(R_CONS), function(s) um_run(Z, s, nn_method = nn, init = Y_sp2))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
saveRDS(list(main = rep_c, second = rep_c2), "results/consensus_replicates.rds")
ali <- gpa(rep_c); M_cons <- mediana(ali)
st_rep <- stabilita_vicini(rep_c, 10)
ripet <- do.call(rbind, lapply(R_LOTTI, function(r) {
  gr <- split(seq_len(R_CONS), ceiling(seq_len(R_CONS) / r)); gr <- gr[lengths(gr) == r]
  cm <- lapply(gr, function(i) mediana(gpa(rep_c[i])))
  data.frame(R = r, consensus_maps = length(cm), stability = stabilita_vicini(cm, 10)$media)
}))
f_rep <- do.call(rbind, lapply(rep_c, valuta_proiezione, D_alto = D0, k_val = K_V))
f_cons <- valuta_proiezione(D0, M_cons, k_val = K_V)
A <- simplify2array(ali)
ell <- do.call(rbind, lapply(seq_len(n_T), function(i) {
  e <- eigen(cov(t(A[i, , ])), symmetric = TRUE); q <- qchisq(0.95, 2)
  data.frame(asse1 = sqrt(q * max(e$values[1], 0)), asse2 = sqrt(q * max(e$values[2], 0)), angolo = atan2(e$vectors[2, 1], e$vectors[1, 1]))
}))
ell$area <- pi * ell$asse1 * ell$asse2
scatola <- prod(apply(M_cons, 2, function(x) diff(range(x))))
compress <- scatola / mean(vapply(ali, function(Y) prod(apply(Y, 2, function(x) diff(range(x)))), 0))
cons <- data.frame(termine = TERMINI, x = M_cons[, 1], y = M_cons[, 2], ell, stabilita = st_rep$per_termine)
scrivi_csv(data.frame(term = cons$termine, x = cons$x, y = cons$y, axis1 = cons$asse1, axis2 = cons$asse2, angle = cons$angolo, area = cons$area,
                      stability = cons$stabilita), "results/consensus_map.csv")
scrivi_csv(ripet, "results/consensus_repeatability.csv")
scrivi_csv(data.frame(item = c("single replicates", "consensus map"),
                      trust = c(mean(f_rep$trustworthiness), f_cons$trustworthiness), continuity = c(mean(f_rep$continuity), f_cons$continuity),
                      neighbours = c(mean(f_rep$conservazione_knn), f_cons$conservazione_knn), spearman = c(mean(f_rep$spearman), f_cons$spearman),
                      mesh = c(mean(rowMeans(sapply(rep_c, function(Y) acc_S(dist_euclidea(Y), SIM, 10)))), mean(acc_S(dist_euclidea(M_cons), SIM, 10))),
                      stability_first20_fixed_spectral = stabilita_vicini(rep_c[SEMI], 10)$media, stability_fixed_pca = stabilita_vicini(rep_p, 10)$media,
                      bounding_box_ratio = compress, ellipse_median_share = median(cons$area) / scatola, ellipse_max_share = max(cons$area) / scatola),
           "results/consensus_summary.csv")

crp <- readRDS("wos/corpus_clean.rds")
LAMBDA <- seq(0, 1, by = 0.1); N_META <- 10; B_RIC <- 20
et <- function(l) sprintf("%.1f", l)
kw <- lapply(strsplit(tolower(ifelse(is.na(crp$DE), "", as.character(crp$DE))), ";", fixed = TRUE), trimws)
pos <- lapply(kw, function(v) { m <- unique(match(v, TERMINI)); m[!is.na(m)] })
X <- Matrix::sparseMatrix(i = rep(seq_along(pos), lengths(pos)), j = unlist(pos), x = 1, dims = c(length(pos), n_T), dimnames = list(NULL, TERMINI))
tiene <- Matrix::rowSums(X) > 0; X <- X[tiene, ]; ANNI <- as.numeric(crp$PY)[tiene]; n_doc <- nrow(X)
rm(crp, kw, pos); invisible(gc())
prof_cooc <- function(C) {
  f <- diag(C); A <- C / outer(pmax(f, 1), pmax(f, 1)); diag(A) <- 0; dimnames(A) <- list(TERMINI, TERMINI)
  D <- suppressWarnings(dist_coseno(A)); D[is.na(D)] <- 1; diag(D) <- 0
  list(D = D, A = A, vuoti = sum(rowSums(A) == 0))
}
Cc <- as.matrix(Matrix::crossprod(X)); fr <- diag(Cc)
pc <- prof_cooc(Cc); D_cooc <- pc$D; AS <- pc$A; D_sem <- D0
D2 <- n_doc * (outer(1 / fr, 1 / fr, "+") - 2 * Cc / outer(fr, fr)); diag(D2) <- 0
D_chi <- sqrt(pmax(D2, 0)); dimnames(D_chi) <- list(TERMINI, TERMINI)
fuori <- function(D) D[lower.tri(D)]
std <- function(D) { x <- fuori(D); S <- (D - mean(x)) / sd(x); diag(S) <- 0; S }
rango <- function(D) { R <- matrix(0, nrow(D), ncol(D), dimnames = dimnames(D)); r <- rank(fuori(D)); R[lower.tri(R)] <- r; R <- R + t(R); R / length(r) }
combina <- function(A, B, l) { D <- l * A + (1 - l) * B; D <- D - min(fuori(D)) + 1e-6; diag(D) <- 0; D }
S_sem <- std(D_sem); S_cooc <- std(D_cooc)
DL <- setNames(lapply(LAMBDA, function(l) combina(S_sem, S_cooc, l)), et(LAMBDA))
SIM3 <- lapply(COD, matrice_sim, termini = TERMINI)
L <- readRDS("results/louvain_runs.rds")
k_L <- leggi_csv("results/benchmark_summary.csv")$louvain_communities
contro_L <- function(P) mean(vapply(L, function(l) ari(P, l), 0))
riga_spazio <- function(D, Ds = D_sem, Dc = D_cooc) c(
  neighbours_sem = conservazione_knn(D, Ds, K_B), neighbours_cooc = conservazione_knn(D, Dc, K_B),
  spearman_sem = correlazione_distanze(D, Ds)[["spearman"]], spearman_cooc = correlazione_distanze(D, Dc)[["spearman"]],
  mesh = mean(acc_S(D, SIM3$finale, K_B)), mesh_before_review = mean(acc_S(D, SIM3$autore, K_B)), mesh_automatic = mean(acc_S(D, SIM3$automatica, K_B)))
tab_sp <- data.frame(lambda = LAMBDA, t(vapply(DL, riga_spazio, numeric(7))), row.names = NULL)
tab_sp$ari_louvain <- vapply(DL, function(D) contro_L(cluster::pam(stats::as.dist(D), k = k_L, diss = TRUE, cluster.only = TRUE)), 0)
mesh_sp <- lapply(DL, acc_S, S = SIM3$finale, k = K_B)
ic_sp <- t(vapply(mesh_sp, function(a) ic_diff(mesh_sp[[et(1)]], a), numeric(3)))
tab_sp$mesh_diff <- ic_sp[, 1]; tab_sp$mesh_low <- ic_sp[, 2]; tab_sp$mesh_high <- ic_sp[, 3]
controllo <- function(nome, A, B, Dc = D_cooc) data.frame(variant = nome, lambda = LAMBDA, t(vapply(LAMBDA, function(l) riga_spazio(combina(A, B, l), Dc = Dc), numeric(7))), row.names = NULL)
tab_ctrl <- rbind(controllo("not standardized", D_sem, D_cooc), controllo("ranks", rango(D_sem), rango(D_cooc)), controllo("chi-square", S_sem, std(D_chi), Dc = D_chi))

da_vicini <- function(D) {
  Nx <- vicini_knn(D, K_B - 1); idx <- cbind(seq_len(n_T), Nx); dimnames(idx) <- NULL; storage.mode(idx) <- "integer"
  list(idx = idx, dist = t(vapply(seq_len(n_T), function(i) unname(D[i, idx[i, ]]), numeric(K_B))))
}
componenti <- function(nn) igraph::components(igraph::graph_from_edgelist(cbind(rep(seq_len(n_T), K_B - 1), as.vector(nn$idx[, -1])), directed = FALSE))$no
um_l <- function(s, nn, ...) { set.seed(s)
  Y <- uwot::umap(Z, n_neighbors = K_B, min_dist = MD_B, metric = "cosine", n_components = 2, nn_method = nn, ...); rownames(Y) <- TERMINI; Y }
init_sp <- function(nn, D) {
  if (componenti(nn) == 1) return(um_l(1, nn, n_epochs = 0))
  Y <- stats::cmdscale(D, k = 2); 10 * Y / max(abs(Y))
}
M <- list()
for (l in et(LAMBDA)) { nn <- da_vicini(DL[[l]]); M[[l]] <- lapply(SEMI, um_l, nn = nn, init = init_sp(nn, DL[[l]])) }
km <- function(Y) { set.seed(SEED); stats::kmeans(Y, centers = k_L, nstart = 100, iter.max = 100)$cluster }
fedelta <- function(l, Dr) { r <- do.call(rbind, lapply(l, function(Y) valuta_proiezione(Dr, Y, k_val = K_V)))
  c(trust = mean(r$trustworthiness), continuity = mean(r$continuity), neighbours = mean(r$conservazione_knn), spearman = mean(r$spearman)) }
stab <- lapply(M, stabilita_vicini, k = K_B)
mesh_m <- lapply(M, function(l) rowMeans(sapply(l, function(Y) acc_S(dist_euclidea(Y), SIM3$finale, K_B))))
tab_m <- do.call(rbind, lapply(et(LAMBDA), function(l) {
  a <- fedelta(M[[l]], DL[[l]]); b <- fedelta(M[[l]], D_sem); d <- fedelta(M[[l]], D_cooc)
  data.frame(lambda = as.numeric(l), stability = stab[[l]]$media, shared = quota_vicini(M[[l]], K_B), mesh = mean(mesh_m[[l]]),
             ari_louvain = mean(vapply(M[[l]], function(Y) contro_L(km(Y)), 0)),
             t(setNames(a, paste0(names(a), "_own"))), t(setNames(b, paste0(names(b), "_sem"))), t(setNames(d, paste0(names(d), "_cooc"))))
}))
ic_mesh <- t(vapply(mesh_m, function(a) ic_diff(mesh_m[[et(1)]], a), numeric(3)))
ic_stab <- t(vapply(stab, function(s) ic_diff(stab[[et(1)]]$per_termine, s$per_termine), numeric(3)))
tab_m$mesh_diff <- ic_mesh[, 1]; tab_m$mesh_low <- ic_mesh[, 2]; tab_m$mesh_high <- ic_mesh[, 3]
tab_m$stab_diff <- ic_stab[, 1]; tab_m$stab_low <- ic_stab[, 2]; tab_m$stab_high <- ic_stab[, 3]
tab_m$components <- vapply(DL, function(D) componenti(da_vicini(D)), 0)

ok_a <- !is.na(ANNI); anno_m <- floor(stats::median(ANNI[ok_a]))
parti <- function(A, B) {
  sa <- std(prof_cooc(as.matrix(Matrix::crossprod(X[A, ])))$D); sb <- std(prof_cooc(as.matrix(Matrix::crossprod(X[B, ])))$D)
  list(v = vapply(LAMBDA, function(l) conservazione_knn(combina(S_sem, sa, l), combina(S_sem, sb, l), K_B), 0),
       s = vapply(LAMBDA, function(l) correlazione_distanze(combina(S_sem, sa, l), combina(S_sem, sb, l))[["spearman"]], 0))
}
per <- parti(which(ok_a & ANNI <= anno_m), which(ok_a & ANNI > anno_m))
cas <- lapply(seq_len(N_META), function(r) { set.seed(SEED + r); i <- sample(n_doc); h <- floor(n_doc / 2); parti(i[seq_len(h)], i[-seq_len(h)]) })
cas_v <- sapply(cas, `[[`, "v"); cas_s <- sapply(cas, `[[`, "s")
tab_c <- data.frame(lambda = LAMBDA, periods_neighbours = per$v, halves_neighbours = rowMeans(cas_v), halves_neighbours_min = apply(cas_v, 1, min),
                    periods_spearman = per$s, halves_spearman = rowMeans(cas_s))

vic_l <- function(Sc, l) vicini_knn(l * S_sem + (1 - l) * Sc, K_B)
mca_piano <- function(C, n) {
  f <- diag(C); i0 <- 2 * seq_len(n_T) - 1; i1 <- 2 * seq_len(n_T)
  Bu <- matrix(0, 2 * n_T, 2 * n_T)
  Bu[i1, i1] <- C; Bu[i1, i0] <- f - C; Bu[i0, i1] <- t(f - C)
  Bu[i0, i0] <- n - outer(f, f, "+") + C
  diag(Bu)[i0] <- n - f; Bu[cbind(i0, i1)] <- 0; Bu[cbind(i1, i0)] <- 0
  P <- Bu / sum(Bu); cm <- colSums(P)
  e <- eigen((P - tcrossprod(cm)) / sqrt(tcrossprod(cm)), symmetric = TRUE)
  Y <- sweep(e$vectors[, 1:2] / sqrt(cm), 2, sqrt(e$values[1:2]), "*")[i1, ]; rownames(Y) <- TERMINI; Y
}
sov <- function(a, b) mean(vapply(seq_len(nrow(a)), function(i) length(intersect(a[i, ], b[i, ])) / K_B, 0))
L_R <- c(0, 0.5, 0.8, 1)
tutte <- function(Sc, Y) c(setNames(lapply(L_R, function(l) vic_l(Sc, l)), paste("lambda", L_R)), list("MCA plane" = vicini_knn(dist_euclidea(Y), K_B)))
intero <- tutte(std(prof_cooc(Cc)$D), mca_piano(Cc, n_doc))
ric <- lapply(seq_len(B_RIC), function(r) { set.seed(SEED + r); i <- sample(n_doc, replace = TRUE)
  C <- as.matrix(Matrix::crossprod(X[i, ])); tutte(std(prof_cooc(C)$D), mca_piano(C, n_doc)) })
cp <- combn(B_RIC, 2)
tab_r <- do.call(rbind, lapply(names(intero), function(v) {
  a <- vapply(ric, function(x) sov(x[[v]], intero[[v]]), 0)
  b <- vapply(seq_len(ncol(cp)), function(j) sov(ric[[cp[1, j]]][[v]], ric[[cp[2, j]]][[v]]), 0)
  data.frame(item = v, with_whole = mean(a), with_whole_min = min(a), between_resamples = mean(b), between_resamples_min = min(b))
}))
scrivi_csv(tab_sp, "results/hybrid_space.csv")
scrivi_csv(tab_ctrl, "results/hybrid_space_controls.csv")
scrivi_csv(tab_m, "results/hybrid_maps.csv")
scrivi_csv(tab_c, "results/hybrid_periods.csv")
scrivi_csv(tab_r, "results/resampling.csv")
scrivi_csv(data.frame(median_year = anno_m, documents_up_to = sum(ok_a & ANNI <= anno_m), documents_after = sum(ok_a & ANNI > anno_m)), "results/hybrid_periods_documents.csv")
sessione("09")
