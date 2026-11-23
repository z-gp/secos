source("R/functions.R")
RNGkind("L'Ecuyer-CMRG")
DEF <- readRDS("wos/corpus_clean.rds")
if (!"DB" %in% names(DEF) || all(is.na(DEF$DB))) DEF$DB <- "ISI"
STOPWORD <- readLines("data/stopwords.txt", encoding = "UTF-8")
PARAM_RETE <- list(n = 50, normalize = "association", cluster = "louvain", community.repulsion = 0.1,
                   remove.isolates = TRUE, remove.multiple = FALSE, noloops = TRUE, weighted = TRUE,
                   edges.min = 0, type = "auto", labelsize = 1, size = 3, Title = "")
PARAM_TM <- list(field = "DE", n = 250, minfreq = 5, ngrams = 1, stemming = FALSE, size = 0.5,
                 n.labels = 3, community.repulsion = 0.1, repel = TRUE, remove.terms = STOPWORD, cluster = "louvain")
chiama <- function(f, args) do.call(f, args[names(args) %in% names(formals(f))])
prendi <- function(df, nomi) { j <- intersect(nomi, names(df)); if (length(j)) df[[j[1]]] else NULL }
quadrante <- function(cl) {
  xc <- prendi(cl, c("rcentrality", "centrality")); yd <- prendi(cl, c("rdensity", "density"))
  ifelse(xc >= mean(xc), ifelse(yd >= mean(yd), "motor", "basic"), ifelse(yd >= mean(yd), "niche", "emerging or declining"))
}

NetMatrix <- biblioNetwork(DEF, analysis = "co-occurrences", network = "author_keywords", sep = ";", remove.terms = STOPWORD)
pdf("results/fig3_network.pdf", width = 10, height = 7.5)
net <- chiama(networkPlot, c(list(NetMatrix = NetMatrix), PARAM_RETE, list(seed = SEED, verbose = TRUE)))
dev.off()
cr <- net$cluster_res
mod_rete <- if (!is.null(net$modularity)) net$modularity else igraph::modularity(net$cluster_obj)
pr <- grep("pagerank", names(cr), value = TRUE)[1]
scrivi_csv(cr[order(cr$cluster, -cr[[pr]]), ], "results/table_S7_network.csv")

tm <- chiama(thematicMap, c(list(M = DEF), PARAM_TM, list(seed = SEED)))
ggsave("results/fig4_thematic_map.pdf", tm$map, width = 9, height = 7)
cl <- as.data.frame(tm$clusters); cl <- cl[, !vapply(cl, is.list, logical(1)), drop = FALSE]
cl$quadrant <- quadrante(cl)
wd <- as.data.frame(tm$words); wd <- wd[, !vapply(wd, is.list, logical(1)), drop = FALSE]
c_occ <- intersect(c("Occurrences", "sC", "freq"), names(wd))[1]; c_g <- intersect(c("Cluster", "groups", "cluster"), names(wd))[1]
wd <- wd[order(wd[[c_g]], -as.numeric(wd[[c_occ]])), ]
scrivi_csv(cl, "results/thematic_map_themes.csv")
scrivi_csv(wd, "results/thematic_map_terms.csv")
scrivi_csv(do.call(rbind, lapply(split(wd, wd[[c_g]]), head, 10)), "results/table_S8_thematic_map.csv")
scrivi_csv(data.frame(nodes = nrow(cr), communities = length(unique(cr$cluster)), modularity = mod_rete,
                      themes = nrow(cl), theme_terms = nrow(wd)), "results/network_summary.csv")

SEMI <- SEED + c(0, 1000, 2000, 3000, 4000)
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
acc_coppie <- function(lst) {
  cp <- utils::combn(length(lst), 2)
  vapply(seq_len(ncol(cp)), function(i) { p <- lst[[cp[1, i]]]; q <- lst[[cp[2, i]]]; t <- intersect(names(p), names(q)); ari(p[t], q[t]) }, numeric(1))
}
pdf(NULL)
reti <- lapply(SEMI, function(s) {
  net <- chiama(networkPlot, c(list(NetMatrix = NetMatrix), PARAM_RETE, list(seed = s, verbose = FALSE)))
  cr <- net$cluster_res; setNames(as.integer(factor(cr$cluster)), as.character(cr$vertex))
})
dev.off()
rm(NetMatrix); invisible(gc())
mappe <- lapply(SEMI, function(s) {
  tm <- chiama(thematicMap, c(list(M = DEF), PARAM_TM, list(seed = s)))
  cl <- as.data.frame(tm$clusters); cl <- cl[, !vapply(cl, is.list, logical(1)), drop = FALSE]
  wd <- as.data.frame(tm$words); wd <- wd[, !vapply(wd, is.list, logical(1)), drop = FALSE]
  gc_id <- prendi(cl, c("groups", "Cluster", "cluster"))
  w_term <- as.character(prendi(wd, c("words", "Words"))); w_grp <- prendi(wd, c("groups", "Cluster", "cluster"))
  invisible(gc())
  list(gruppi = setNames(as.integer(factor(w_grp)), w_term), grp = w_grp, term = w_term, occ = prendi(wd, c("sC", "Occurrences", "freq")),
       quad = setNames(quadrante(cl), as.character(gc_id)), n = nrow(cl))
})
base <- mappe[[1]]
temi <- do.call(rbind, lapply(unique(base$grp), function(g) {
  tb <- base$term[base$grp == g]; tb <- tb[order(-base$occ[base$grp == g])]
  js <- vapply(mappe[-1], function(m) max(vapply(unique(m$grp), function(h) jacc(tb, m$term[m$grp == h]), numeric(1))), numeric(1))
  qs <- vapply(mappe, function(m) { i <- match(tb[1], m$term); if (is.na(i)) NA_character_ else unname(m$quad[as.character(m$grp[i])]) }, "")
  data.frame(theme = tb[1], terms = length(tb), mean_jaccard = mean(js), min_jaccard = min(js), quadrant = qs[1],
             same_quadrant = sum(qs == qs[1], na.rm = TRUE), stringsAsFactors = FALSE)
}))
temi$stable <- temi$mean_jaccard >= 0.70 & temi$same_quadrant >= 4
ari_rete <- acc_coppie(reti); ari_tm <- acc_coppie(lapply(mappe, `[[`, "gruppi"))
scrivi_csv(temi[order(-temi$terms), ], "results/theme_stability.csv")
scrivi_csv(data.frame(analysis = c("network", "thematic map"),
                      clusters_per_seed = c(paste(vapply(reti, function(r) length(unique(r)), 1L), collapse = ", "),
                                            paste(vapply(mappe, `[[`, 1L, "n"), collapse = ", ")),
                      mean_ari = c(mean(ari_rete), mean(ari_tm)), min_ari = c(min(ari_rete), min(ari_tm))), "results/seed_stability.csv")
sessione("04")
