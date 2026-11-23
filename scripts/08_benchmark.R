source("R/functions.R")
K_FED <- c(5, 10, 15, 30); K_MAPPA <- 10; KS <- c(5, 10, 20); PERM <- 999
N_LOUVAIN <- 200; N_COPPIE <- 2000; MIN_PTS <- 5; K_GRIGLIA <- 4:12; NSTART <- 100
TERMINI <- leggi_csv("results/terms.csv")$term; n_T <- length(TERMINI)
ENC <- readRDS("results/embeddings.rds"); um <- readRDS("results/umap_runs.rds")
mc <- leggi_csv("results/mca_coordinates.csv")
Y_mca <- as.matrix(mc[, c("Dim1", "Dim2")]); dimnames(Y_mca) <- list(mc$term, c("Dim1", "Dim2")); Y_mca <- Y_mca[TERMINI, ]
Cc <- as.matrix(read.csv("results/cooccurrence_matrix.csv", row.names = 1, check.names = FALSE, encoding = "UTF-8"))[TERMINI, TERMINI]
n_doc <- leggi_csv("results/mca_summary.csv")$cooccurrence_documents; fr <- diag(Cc)
D2 <- n_doc * (outer(1 / fr, 1 / fr, "+") - 2 * Cc / outer(fr, fr)); diag(D2) <- 0
D_chi <- sqrt(pmax(D2, 0)); dimnames(D_chi) <- list(TERMINI, TERMINI)
AS <- Cc / outer(fr, fr); diag(AS) <- 0; dimnames(AS) <- list(TERMINI, TERMINI)
D_prof <- dist_coseno(AS)
D_sb <- dist_coseno(ENC$sbert[TERMINI, ]); D_sc <- dist_coseno(ENC$scibert[TERMINI, ])
RIF <- list("chi-square" = D_chi, "co-occurrence" = D_prof, "SBERT" = D_sb, "SciBERT" = D_sc)
PIANI <- list("MCA" = list(Y_mca), "SeCoS SBERT" = um$sbert, "SeCoS SciBERT" = um$scibert)
PROPRIO <- c("MCA" = "chi-square", "SeCoS SBERT" = "SBERT", "SeCoS SciBERT" = "SciBERT")

fedelta <- function(lista, Dr, k) {
  r <- do.call(rbind, lapply(lista, function(Y) valuta_proiezione(Dr, Y, k_val = k)))
  c(trust = mean(r$trustworthiness), trust_sd = sd(r$trustworthiness), cont = mean(r$continuity), cont_sd = sd(r$continuity),
    neighbours = mean(r$conservazione_knn), neighbours_sd = sd(r$conservazione_knn), spearman = mean(r$spearman), spearman_sd = sd(r$spearman))
}
griglia <- expand.grid(map = names(PIANI), reference = names(RIF), k = K_FED, stringsAsFactors = FALSE)
griglia <- griglia[!(griglia$map == "SeCoS SBERT" & griglia$reference == "SciBERT") & !(griglia$map == "SeCoS SciBERT" & griglia$reference == "SBERT"), ]
tab_fed <- cbind(griglia, t(vapply(seq_len(nrow(griglia)), function(i) fedelta(PIANI[[griglia$map[i]]], RIF[[griglia$reference[i]]], griglia$k[i]), numeric(8))))
D_media <- lapply(PIANI, function(l) Reduce(`+`, lapply(l, dist_euclidea)) / length(l))
mantel <- vapply(names(PIANI), function(p) correlazione_distanze(RIF[[PROPRIO[[p]]]], D_media[[p]], permutazioni = PERM, seed = SEED)[["p_mantel"]], 0)

acc_lista <- function(l, S, k) rowMeans(sapply(l, function(Y) acc_S(dist_euclidea(Y), S, k)))
COD <- mesh_alberi()
SIM <- lapply(COD, matrice_sim, termini = TERMINI)
caso <- mean(SIM$finale[lower.tri(SIM$finale)])
massimo <- vapply(KS, function(k) mean(apply(`diag<-`(SIM$finale, -1), 1, function(x) mean(sort(x, TRUE)[seq_len(k)]))), 0)
tab_mesh <- do.call(rbind, lapply(names(COD), function(cd) do.call(rbind, lapply(KS, function(k) {
  S <- SIM[[cd]]; a0 <- acc_lista(PIANI$MCA, S, k)
  voci <- c(setNames(lapply(names(PIANI), function(p) acc_lista(PIANI[[p]], S, k)), names(PIANI)),
            setNames(lapply(RIF, function(D) acc_S(D, S, k)), paste("space", names(RIF))))
  do.call(rbind, lapply(names(voci), function(v) {
    ic <- if (v == "MCA") c(0, NA, NA) else ic_diff(a0, voci[[v]])
    data.frame(coding = cd, terms = nrow(S), k = k, item = v, agreement = mean(voci[[v]]), difference = ic[[1]], ci_low = ic[[2]], ci_high = ic[[3]])
  }))
}))))
tt <- rownames(SIM$finale)
tutte <- vapply(names(PIANI), function(p) correlazione_distanze(D_media[[p]][tt, tt], 1 - SIM$finale, permutazioni = PERM, seed = SEED), c(spearman = 0, p_mantel = 0))

g <- igraph::graph_from_adjacency_matrix(AS, mode = "undirected", weighted = TRUE, diag = FALSE)
LV <- lapply(seq_len(N_LOUVAIN), function(r) { set.seed(SEED + r - 1); igraph::cluster_louvain(g, resolution = 1) })
mod <- vapply(LV, igraph::modularity, 0)
L <- lapply(LV, function(x) as.integer(igraph::membership(x)))
L_best <- L[[which.max(mod)]]; k_L <- length(unique(L_best))
saveRDS(L, "results/louvain_runs.rds")
set.seed(SEED)
cp <- matrix(sample(N_LOUVAIN, 2 * N_COPPIE, replace = TRUE), ncol = 2); cp <- cp[cp[, 1] != cp[, 2], ]
rif_ari <- vapply(seq_len(nrow(cp)), function(i) ari(L[[cp[i, 1]]], L[[cp[i, 2]]]), 0)
rif_nmi <- vapply(seq_len(nrow(cp)), function(i) nmi(L[[cp[i, 1]]], L[[cp[i, 2]]]), 0)
km <- function(Mx, k) { set.seed(SEED); stats::kmeans(Mx, centers = k, nstart = NSTART, iter.max = 100)$cluster }
contro_L <- function(P, ok = rep(TRUE, n_T)) {
  a <- vapply(L, function(l) ari(P[ok], l[ok]), 0); v <- vapply(L, function(l) nmi(P[ok], l[ok]), 0)
  c(ari_mean = mean(a), ari_sd = sd(a), nmi_mean = mean(v), ari_best = ari(P[ok], L_best[ok]), nmi_best = nmi(P[ok], L_best[ok]))
}
Z_sb <- norm_l2(ENC$sbert[TERMINI, ]); Z_sc <- norm_l2(ENC$scibert[TERMINI, ])
P_sb <- km(Z_sb, k_L); P_sc <- km(Z_sc, k_L); P_mca <- km(Y_mca, k_L)
partiz <- rbind("k-means, SBERT embeddings" = contro_L(P_sb), "k-means, SciBERT embeddings" = contro_L(P_sc),
                "k-means, MCA plane" = contro_L(P_mca),
                "k-means, SeCoS SBERT maps" = rowMeans(sapply(um$sbert, function(Y) contro_L(km(Y, k_L)))),
                "k-means, SeCoS SciBERT maps" = rowMeans(sapply(um$scibert, function(Y) contro_L(km(Y, k_L)))))
griglia_k <- do.call(rbind, lapply(K_GRIGLIA, function(k) {
  Ps <- km(Z_sb, k); Pm <- km(Y_mca, k)
  data.frame(k = k, ari_embeddings = mean(vapply(L, function(l) ari(Ps, l), 0)), ari_mca = mean(vapply(L, function(l) ari(Pm, l), 0)))
}))
hdb <- dbscan::hdbscan(stats::as.dist(D_sb), minPts = MIN_PTS)
P_hdb <- as.integer(hdb$cluster); ok_h <- P_hdb != 0
h_ris <- if (length(unique(P_hdb[ok_h])) > 1) contro_L(P_hdb, ok_h) else rep(NA_real_, 5)

N_mca <- vicini_knn(dist_euclidea(Y_mca), K_MAPPA)
sov_piani <- rowMeans(sapply(um$sbert, function(Y) { N <- vicini_knn(dist_euclidea(Y), K_MAPPA)
  vapply(seq_len(n_T), function(i) length(intersect(N[i, ], N_mca[i, ])) / K_MAPPA, 0) }))
lt <- lower.tri(AS)
ij <- which(lt, arr.ind = TRUE); dif <- rank(-AS[lt], ties.method = "average") - rank(D_sb[lt])
coppia <- function(o, tipo) data.frame(type = tipo, a = TERMINI[ij[o, 1]], b = TERMINI[ij[o, 2]], documents = as.integer(Cc[ij[o, , drop = FALSE]]))

scrivi_csv(tab_fed, "results/benchmark_fidelity.csv")
scrivi_csv(data.frame(map = names(mantel), mantel_p = mantel, spearman_mesh_all_pairs = tutte["spearman", ], mantel_p_mesh = tutte["p_mantel", ]), "results/benchmark_mantel.csv")
tab_mesh$coding <- unname(CODING_EN[tab_mesh$coding]); stopifnot(!anyNA(tab_mesh$coding))
scrivi_csv(tab_mesh, "results/benchmark_mesh.csv")
scrivi_csv(data.frame(k = KS, chance = caso, maximum = massimo), "results/mesh_reference.csv")
scrivi_csv(rbind(data.frame(partition = rownames(partiz), partiz, row.names = NULL),
                 data.frame(partition = c("HDBSCAN, SBERT embeddings", "Louvain with itself"), rbind(h_ris, c(mean(rif_ari), sd(rif_ari), mean(rif_nmi), NA, NA)), row.names = NULL)),
           "results/benchmark_partitions.csv")
scrivi_csv(griglia_k, "results/benchmark_partitions_grid.csv")
scrivi_csv(data.frame(term = TERMINI, documents = fr, louvain = L_best, kmeans_sbert = P_sb, kmeans_mca = P_mca, hdbscan_sbert = P_hdb, shared_with_mca = sov_piani),
           "results/benchmark_terms.csv")
scrivi_csv(data.frame(louvain_communities = k_L, hdbscan_clusters = length(unique(P_hdb[ok_h])), hdbscan_unassigned = sum(!ok_h),
                      shared_neighbours_maps = mean(sov_piani), shared_neighbours_chance = K_MAPPA / (n_T - 1)), "results/benchmark_summary.csv")
scrivi_csv(rbind(coppia(head(order(-dif), 12), "close in meaning"), coppia(head(order(dif), 12), "close in co-occurrence")), "results/benchmark_pairs.csv")
sessione("08")
