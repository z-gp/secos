source("R/functions.R")
source("R/secos.R")
suppressPackageStartupMessages(library(patchwork))
TERMINI <- leggi_csv("results/terms.csv")$term; n_T <- length(TERMINI)
ENC <- readRDS("results/embeddings.rds")
RC <- readRDS("results/consensus_replicates.rds"); rep_c <- RC$main
cons <- leggi_csv("results/consensus_map.csv"); cons <- cons[match(TERMINI, cons$term), ]
part <- leggi_csv("results/benchmark_terms.csv"); part <- part[match(TERMINI, part$term), ]
stab <- leggi_csv("results/stability_terms.csv"); stab <- stab[match(TERMINI, stab$term), ]
griglia <- leggi_csv("results/table_S9_sensitivity.csv")
sp <- leggi_csv("results/hybrid_space.csv"); mp <- leggi_csv("results/hybrid_maps.csv"); dc <- leggi_csv("results/hybrid_periods.csv")
mc <- leggi_csv("results/mca_coordinates.csv"); Y_mca <- as.matrix(mc[match(TERMINI, mc$term), c("Dim1", "Dim2")])
D_sem <- dist_coseno(ENC$sbert[TERMINI, ])
fai <- function(nome, codice) invisible(codice)
o <- order(-stab$stability_secos); u <- order(stab$stability_secos)
scrivi_csv(data.frame(rank = 1:10, most_stable = TERMINI[o[1:10]], stability = stab$stability_secos[o[1:10]],
                      least_stable = TERMINI[u[1:10]], stability_least = stab$stability_secos[u[1:10]]), "results/table6_stability.csv")
SOGLIA_AREA <- 0.10
D_MIN <- 4
QUOTA_MIN <- 0.10
N_PER_GRUPPO <- 2
N_TOP <- 15
R_LOTTO <- 20
GUIDE <- c("sports injury" = "Injuries and trauma", "resistance training" = "Training and exercise therapy",
           "athlete" = "Athletes and sports", "epidemiology" = "Epidemiology, health and psychology",
           "rehabilitation" = "Clinical care and research methods", "biomechanics" = "Biomechanics and movement",
           "anterior cruciate ligament" = "Joints and musculoskeletal conditions", "adolescent" = "Age and sex")
ABBR <- c("anterior cruciate ligament" = "ACL", "posterior cruciate ligament" = "PCL",
          "magnetic resonance imaging" = "MRI", "mild traumatic brain injury" = "mild TBI",
          "traumatic brain injury" = "TBI", "chronic traumatic encephalopathy" = "CTE")
accorcia <- function(v) { for (a in names(ABBR)) v <- sub(a, ABBR[[a]], v, fixed = TRUE); v }
vuoto <- theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(fill = NA, colour = CHIARO, linewidth = 0.4))
zone <- function(P, med, seme) {
  set.seed(seme); km <- stats::kmeans(P, centers = 2, nstart = 20, iter.max = 50)
  u <- km$centers[2, ] - km$centers[1, ]; L <- sqrt(sum(u^2))
  pr <- drop(P %*% (u / L)); g <- km$cluster
  D <- sqrt(2) * abs(diff(tapply(pr, g, mean))) / sqrt(sum(tapply(pr, g, var)))
  raggio <- sqrt(km$tot.withinss / nrow(P))
  list(D = unname(D), quota_min = min(table(g)) / nrow(P), distanza = L,
       presso_mediana = mean(sqrt(rowSums(sweep(P, 2, med)^2)) <= raggio), g = g)
}
dieci_in <- function(A, i, rr) {
  f <- table(unlist(lapply(rr, function(r) { d <- sqrt(rowSums(sweep(A[, , r], 2, A[i, , r])^2)); d[i] <- Inf; order(d)[1:10] })))
  as.integer(names(sort(f, decreasing = TRUE)[1:10]))
}
incerti <- function(rep) {
  ali <- gpa(rep); M <- mediana(ali)
  A <- simplify2array(ali)[match(TERMINI, rownames(ali[[1]])), , ]
  scatola <- prod(apply(M, 2, function(x) diff(range(x))))
  area <- vapply(seq_len(n_T), function(i) { e <- eigen(cov(t(A[i, , ])), symmetric = TRUE)$values; q <- qchisq(0.95, 2)
    pi * sqrt(q * max(e[1], 0)) * sqrt(q * max(e[2], 0)) }, 0)
  Z <- lapply(seq_len(n_T), function(i) zone(t(A[i, , ]), M[TERMINI[i], ], SEED + i))
  ctrl <- data.frame(termine = TERMINI, quota_area = area / scatola, D = vapply(Z, `[[`, 0, "D"), quota_min = vapply(Z, `[[`, 0, "quota_min"))
  ctrl$due_zone <- ctrl$D >= D_MIN & ctrl$quota_min >= QUOTA_MIN
  ctrl$incerto <- ctrl$quota_area >= SOGLIA_AREA | ctrl$due_zone
  ctrl$vicini_comuni <- NA_real_
  for (i in which(ctrl$due_zone)) ctrl$vicini_comuni[i] <- length(intersect(dieci_in(A, i, which(Z[[i]]$g == 1)), dieci_in(A, i, which(Z[[i]]$g == 2)))) / 10
  ctrl
}
ctrl <- incerti(rep_c)
ctrl2 <- incerti(RC$second)
cons$quota_area <- ctrl$quota_area
scrivi_csv(data.frame(term = TERMINI, documents = part$documents, ellipse_share = ctrl$quota_area, ashman_D = ctrl$D, smaller_zone = ctrl$quota_min,
                      two_zones = ctrl$due_zone, uncertain = ctrl$incerto, shared_neighbours_zones = ctrl$vicini_comuni, uncertain_second_set = ctrl2$incerto),
           "results/uncertain_terms.csv")
scrivi_csv(data.frame(two_zones = sum(ctrl$due_zone), large_ellipse = sum(ctrl$quota_area >= SOGLIA_AREA), uncertain = sum(ctrl$incerto),
                      shared_neighbours_median = median(ctrl$vicini_comuni, na.rm = TRUE), shared_neighbours_0.7 = sum(ctrl$vicini_comuni >= 0.7, na.rm = TRUE),
                      shared_neighbours_below_0.5 = sum(ctrl$vicini_comuni < 0.5, na.rm = TRUE), uncertain_second_set = sum(ctrl2$incerto),
                      in_common = sum(ctrl$incerto & ctrl2$incerto), added = sum(!ctrl$incerto & ctrl2$incerto)), "results/uncertain_summary.csv")
g_num <- part$kmeans_sbert
id_g <- sort(unique(g_num))
freq_g <- function(k, n = 2) { j <- which(g_num == k); accorcia(TERMINI[j[order(-part$documents[j])][seq_len(min(n, length(j)))]]) }
nome_g <- setNames(vapply(id_g, function(k) paste(freq_g(k), collapse = ", "), ""), id_g)
gd <- intersect(names(GUIDE), TERMINI)
k_gd <- g_num[match(gd, TERMINI)]
guide_ok <- length(gd) == length(GUIDE) && !anyDuplicated(k_gd) && length(id_g) == length(GUIDE)
if (guide_ok) nome_g[as.character(k_gd)] <- GUIDE[gd]
nome_g <- nome_g[order(-table(g_num)[names(nome_g)])]
GRUPPO <- factor(nome_g[as.character(g_num)], levels = unname(nome_g))
col_g <- setNames(COL[seq_along(nome_g)], unname(nome_g))
medoide <- vapply(id_g, function(k) { j <- which(g_num == k); j[which.min(rowMeans(D_sem[j, j, drop = FALSE]))] }, 0L)
sel <- unique(c(medoide, unlist(lapply(id_g, function(k) { j <- which(g_num == k); j[order(-part$documents[j])][seq_len(min(N_PER_GRUPPO, length(j)))] }))))
dati <- data.frame(termine = TERMINI, nome = accorcia(TERMINI), gruppo = GRUPPO, louvain = part$louvain, documenti = part$documents,
                   x = cons$x, y = cons$y, mx = Y_mca[, 1], my = Y_mca[, 2], asse1 = cons$axis1, asse2 = cons$axis2, angolo = cons$angle,
                   quota_area = cons$quota_area, stabilita = stab$stability_secos, incerto = ctrl$incerto, etichetta = seq_len(n_T) %in% sel,
                   stringsAsFactors = FALSE)
dati$stato <- factor(ifelse(dati$incerto, "Uncertain position", "Stable position"), levels = c("Stable position", "Uncertain position"))
elisse <- function(d, n = 40) do.call(rbind, lapply(seq_len(nrow(d)), function(i) { t <- seq(0, 2 * pi, length.out = n)
  a <- d$asse1[i] * cos(t); b <- d$asse2[i] * sin(t); w <- d$angolo[i]
  data.frame(termine = d$termine[i], gruppo = d[i, "colore"], incerto = d$incerto[i],
             x = d$x[i] + a * cos(w) - b * sin(w), y = d$y[i] + a * sin(w) + b * cos(w)) }))
rip <- function(d, s = 2.05) ggrepel::geom_text_repel(data = d, aes(x, y, label = nome), inherit.aes = FALSE, size = s, colour = INCHIOSTRO,
  fontface = ifelse(d$incerto, "italic", "plain"), seed = SEED, max.overlaps = Inf, box.padding = 0.22, point.padding = 0.12, force = 1.4, max.iter = 20000,
  min.segment.length = 0.1, segment.size = 0.18, segment.colour = GRIGIO)
finestra <- function(v) { q <- stats::quantile(v, c(0.25, 0.75)); c(max(min(v), q[[1]] - 3 * diff(q)), min(max(v), q[[2]] + 3 * diff(q))) }
coppia_mappe <- function(d, scala, titolo_b) {
  d$colore <- d[["cl"]]
  a <- transform(d, x = mx, y = my)
  fx <- finestra(a$x); fy <- finestra(a$y); inq <- diff(range(a$x)) > 1.6 * diff(fx) || diff(range(a$y)) > 1.6 * diff(fy)
  if (!inq) { fx <- range(a$x); fy <- range(a$y) }
  fx <- fx + c(-1, 1) * 0.05 * diff(fx); fy <- fy + c(-1, 1) * 0.05 * diff(fy)
  dentro <- a$x >= fx[1] & a$x <= fx[2] & a$y >= fy[1] & a$y <= fy[2]
  pa <- ggplot(a[dentro, ], aes(x, y)) +
    geom_hline(yintercept = 0, colour = CHIARO, linewidth = 0.2) + geom_vline(xintercept = 0, colour = CHIARO, linewidth = 0.2) +
    geom_point(aes(colour = colore, size = documenti), alpha = 0.9, stroke = 0) +
    rip(a[dentro & a$etichetta, ]) + scala + scale_size_area(max_size = 3, guide = "none") + guides(colour = "none") +
    coord_cartesian(xlim = fx, ylim = fy, expand = FALSE) + labs(title = "(a) MCA, first two dimensions", subtitle = if (inq) "Inset: the whole plane; the rectangle marks the area shown" else " ") + tema + vuoto
  if (inq) {
    cc <- table(factor(a$x[dentro] > mean(fx), c(TRUE, FALSE)), factor(a$y[dentro] > mean(fy), c(TRUE, FALSE)))
    w <- which(cc == min(cc), arr.ind = TRUE)[1, ]
    sx <- c(0.70, 0.02)[w[1]]; sy <- c(0.70, 0.02)[w[2]]
    ins <- ggplot(a, aes(x, y)) + geom_point(colour = GRIGIO, size = 0.25, stroke = 0) +
      annotate("rect", xmin = fx[1], xmax = fx[2], ymin = fy[1], ymax = fy[2], fill = NA, colour = INCHIOSTRO, linewidth = 0.25) +
      theme_void() + theme(panel.border = element_rect(fill = NA, colour = CHIARO, linewidth = 0.3), panel.background = element_rect(fill = "white", colour = NA))
    pa <- pa + inset_element(ins, sx, sy, sx + 0.28, sy + 0.28, align_to = "panel")
  }
  el <- elisse(d)
  pb <- ggplot(d, aes(x, y)) +
    geom_polygon(data = el[!el$incerto, ], aes(group = termine, fill = gruppo), alpha = 0.07, colour = NA) +
    geom_polygon(data = el[el$incerto, ], aes(group = termine), fill = NA, colour = GRIGIO, linewidth = 0.2, linetype = "22") +
    geom_point(aes(colour = colore, size = documenti, shape = stato), alpha = 0.9, stroke = 0.45) +
    rip(d[d$etichetta, ]) + scala + scale_fill_manual(values = attr(scala, "valori"), guide = "none") +
    scale_shape_manual(values = c("Stable position" = 19, "Uncertain position" = 1), guide = "none") +
    scale_size_area(max_size = 3, guide = "none") + coord_equal() + labs(title = titolo_b, subtitle = "Shaded: 95% ellipses. Hollow marker, dashed ellipse: uncertain position") + tema + vuoto
  list(p = (pa | pb) + plot_layout(guides = "collect") & theme(legend.position = "bottom"), fuori = a$termine[!dentro], inq = inq)
}
s_col <- function(val, n) { s <- scale_colour_manual(values = val, guide = guide_legend(nrow = n, override.aes = list(size = 2.2, alpha = 1, shape = 19))); attr(s, "valori") <- val; s }
titolo_b <- expression(bold("(b) SeCoS consensus map (" * lambda * " = 1), 100 replicates"))
fai("p18_fig9_mappe", { d <- dati; d$cl <- d$gruppo; fig9 <- coppia_mappe(d, s_col(col_g, 3), titolo_b); salva(fig9$p, "fig6_maps", 180, 104) })
fai("p18_figS_mappe_louvain", {
  o <- names(sort(table(dati$louvain), decreasing = TRUE)); d <- dati
  d$cl <- factor(paste("Community", match(d$louvain, o)), levels = paste("Community", seq_along(o)))
  v <- setNames(COL[seq_along(o)], levels(d$cl)); r <- coppia_mappe(d, s_col(v, 2), titolo_b); salva(r$p, "figS9_maps_louvain", 180, 104) })
zoom_fuori <- list()
fai("p18_figS_zoom", {
  pz <- do.call(rbind, lapply(levels(GRUPPO), function(g) { m <- dati[dati$gruppo == g, ]
    fx <- finestra(m$x); fy <- finestra(m$y); lon <- m$x < fx[1] | m$x > fx[2] | m$y < fy[1] | m$y > fy[2]
    zoom_fuori[[g]] <<- m$nome[lon]; m <- m[!lon, ]
    rx <- range(m$x) + c(-1, 1) * 0.08 * diff(range(m$x)); ry <- range(m$y) + c(-1, 1) * 0.08 * diff(range(m$y))
    f <- dati[dati$gruppo != g & dati$x >= rx[1] & dati$x <= rx[2] & dati$y >= ry[1] & dati$y <= ry[2], ]
    rbind(transform(m, pannello = g, sfondo = FALSE), if (nrow(f)) transform(f, pannello = g, sfondo = TRUE)) }))
  pz$pannello <- factor(pz$pannello, levels = levels(GRUPPO))
  p <- ggplot(pz, aes(x, y)) + geom_point(data = pz[pz$sfondo, ], colour = CHIARO, size = 0.7, stroke = 0) +
    geom_point(data = pz[!pz$sfondo, ], aes(colour = gruppo, shape = stato), size = 1.3, stroke = 0.4) +
    ggrepel::geom_text_repel(data = pz[!pz$sfondo, ], aes(label = nome), size = 2.05, colour = INCHIOSTRO, seed = SEED, max.overlaps = Inf,
                             box.padding = 0.22, point.padding = 0.15, force = 2.5, min.segment.length = 0.15, segment.size = 0.15, segment.colour = GRIGIO, max.iter = 20000,
                             fontface = ifelse(pz$incerto[!pz$sfondo], "italic", "plain")) +
    facet_wrap(~pannello, ncol = 2, scales = "free") + scale_colour_manual(values = col_g, guide = "none") +
    scale_shape_manual(values = c("Stable position" = 19, "Uncertain position" = 1)) + tema + vuoto
  salva(p, "figS8_zoom", 180, 235) })
lotti <- split(seq_along(rep_c), ceiling(seq_along(rep_c) / R_LOTTO))
lotti <- lotti[lengths(lotti) == R_LOTTO]
cm <- lapply(lotti, function(i) mediana(gpa(rep_c[i])))
st_cons <- stabilita_vicini(lapply(cm, function(Y) Y[TERMINI, ]), 10)
caso <- jaccard_casuale(n_T, 10)
s10 <- rbind(data.frame(metodo = "SeCoS, single run", valore = stab$stability_secos), data.frame(metodo = "t-SNE, single run", valore = stab$stability_tsne),
             data.frame(metodo = sprintf("SeCoS, consensus of %d runs", R_LOTTO), valore = st_cons$per_termine))
liv <- c("t-SNE, single run", "SeCoS, single run", sprintf("SeCoS, consensus of %d runs", R_LOTTO))
s10$metodo <- factor(s10$metodo, levels = liv)
med <- tapply(s10$valore, s10$metodo, stats::median)
fai("p18_fig10_stabilita", {
  v <- setNames(c(COL[2], COL[1], COL[3]), liv)
  dd <- data.frame(metodo = factor(liv, levels = liv), m = as.numeric(med[liv]), h = c(0.96, 0.84, 0.72))
  p <- ggplot(s10, aes(valore, colour = metodo, fill = metodo)) +
    geom_density(alpha = 0.16, linewidth = 0.5, bounds = c(0, 1)) +
    geom_vline(data = dd, aes(xintercept = m, colour = metodo), linewidth = 0.45, show.legend = FALSE) +
    do.call(geom_label, c(list(data = dd, mapping = aes(x = m, y = Inf, label = sprintf("median %.2f", m)), vjust = 1.2 + 1.5 * (seq_along(liv) - 1), size = 2.4, colour = INCHIOSTRO, fill = "white",
               label.padding = unit(0.6, "mm"), show.legend = FALSE), if (utils::packageVersion("ggplot2") >= "3.5.0") list(linewidth = 0) else list(label.size = 0))) +
    geom_vline(xintercept = c(caso, 1), colour = GRIGIO, linewidth = 0.3, linetype = "22") +
    annotate("text", x = caso, y = 0, label = sprintf("  random %.2f", caso), angle = 90, hjust = 0, vjust = 1.3, size = 2.4, colour = INCHIOSTRO) +
    annotate("text", x = 1, y = 0, label = "  MCA 1 (deterministic)", angle = 90, hjust = 0, vjust = -0.4, size = 2.4, colour = INCHIOSTRO) +
    scale_colour_manual(values = v) + scale_fill_manual(values = v) + scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = expansion(mult = 0.01)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + labs(x = "Neighborhood stability of a term (Jaccard, 10 neighbors)", y = "Density") + tema
  salva(p, "fig7_stability", 130, 85) })
PAN <- c("(a) Neighbors shared, original space", "(b) Map trustworthiness", "(c) Agreement with MeSH",
         "(d) Agreement with Louvain (ARI)", "(e) Stability of the map", "(f) Neighbors shared, subcorpora")
lin <- function(p, l, v, s) data.frame(pannello = PAN[p], lambda = l, valore = v, serie = s)
m1 <- mp$mesh[mp$lambda == 1]
s1 <- mp$stability[mp$lambda == 1]
L <- rbind(lin(1, sp$lambda, sp$neighbours_sem, "Semantic"), lin(1, sp$lambda, sp$neighbours_cooc, "Co-occurrence"),
           lin(2, mp$lambda, mp$trust_sem, "Semantic"), lin(2, mp$lambda, mp$trust_cooc, "Co-occurrence"),
           lin(3, mp$lambda, mp$mesh, "Map"), lin(3, sp$lambda, sp$mesh, "Original space"),
           lin(4, mp$lambda, mp$ari_louvain, "Map"), lin(4, sp$lambda, sp$ari_louvain, "Original space"),
           lin(5, mp$lambda, mp$stability, "Map"),
           lin(6, dc$lambda, dc$periods_neighbours, "Two periods"), lin(6, dc$lambda, dc$halves_neighbours, "Random halves"))
L$pannello <- factor(L$pannello, levels = PAN)
SER <- c("Semantic" = COL[1], "Co-occurrence" = COL[2], "Map" = INCHIOSTRO, "Original space" = GRIGIO, "Two periods" = COL[7], "Random halves" = GRIGIO)
TIP <- c("Semantic" = "solid", "Co-occurrence" = "solid", "Map" = "solid", "Original space" = "22", "Two periods" = "solid", "Random halves" = "22")
fai("p18_fig_lambda", {
  rb <- rbind(data.frame(pannello = PAN[3], lambda = mp$lambda, lo = m1 + mp$mesh_low, hi = m1 + mp$mesh_high),
              data.frame(pannello = PAN[5], lambda = mp$lambda, lo = s1 + mp$stab_low, hi = s1 + mp$stab_high))
  rb$pannello <- factor(rb$pannello, levels = PAN)
  et <- data.frame(p = c(1, 1, 2, 2, 3, 3, 4, 4, 6, 6), x = c(0.20, 0.98, 0.05, 0.97, 0.72, 0.12, 0.78, 0.40, 0.03, 0.97), y = c(0.97, 0.62, 0.87, 0.77, 0.19, 0.075, 0.20, 0.405, 0.70, 0.45),
                   h = c(0, 1, 0, 1, 0, 0, 0, 0, 0, 1), testo = c("Co-occurrence", "Semantic", "Co-occurrence", "Semantic", "Map (solid)", "Original space (dashed)", "Map (solid)", "Original space (dashed)", "Random halves (dashed)", "Two periods (solid)"))
  et$pannello <- factor(PAN[et$p], levels = PAN)
  p <- ggplot(L, aes(lambda, valore)) + geom_ribbon(data = rb, aes(x = lambda, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = CHIARO, alpha = 0.7) +
    geom_line(aes(colour = serie, linetype = serie), linewidth = 0.55) + geom_point(aes(colour = serie), size = 0.8, stroke = 0) +
    geom_text(data = et, aes(x, y, label = testo, hjust = h), inherit.aes = FALSE, size = 2.3, colour = INCHIOSTRO) +
    facet_wrap(~pannello, ncol = 3, scales = "free_y") + expand_limits(y = 0) + scale_colour_manual(values = SER, guide = "none") + scale_linetype_manual(values = TIP, guide = "none") +
    scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) + labs(x = expression(lambda), y = NULL) + tema
  salva(p, "fig9_lambda", 180, 105) })
fai("p18_fig_sensibilita", {
  b <- griglia[grepl("^base", griglia$setting), ]; v <- griglia$setting
  riga <- function(sel, gr, et) { x <- griglia[sel, ]; x$gruppo <- gr; x$livello <- et; x }
  etT <- function(x) sprintf("%s (T = %d)", sub("threshold ", "", x), griglia$T[match(x, v)])
  G <- rbind(riga(c(which(grepl("^k = ", v)), which(grepl("^base", v))), "Neighbors k", NA), riga(c(which(grepl("^min_dist", v)), which(grepl("^base", v))), "min_dist", NA),
             riga(c(which(v %in% c("Euclidean", "correlation")), which(grepl("^base", v))), "Metric", NA), riga(c(which(grepl("^threshold", v)), which(grepl("^base", v))), "Frequency threshold", NA),
             riga(c(which(grepl("^encoder", v)), which(grepl("^base", v))), "Encoder", NA))
  G$base <- grepl("^base", G$setting)
  G$livello <- ifelse(G$base, ifelse(G$gruppo == "Neighbors k", "10", ifelse(G$gruppo == "min_dist", "0.1", ifelse(G$gruppo == "Metric", "Cosine", ifelse(G$gruppo == "Encoder", "SBERT", sprintf("50 (T = %d)", b$T))))),
    ifelse(grepl("^threshold", G$setting), etT(G$setting), sub("^k = ", "", sub("^min_dist = ", "", sub("^encoder ", "", G$setting)))))
  G$livello <- sub("^correlation$", "Correlation", G$livello)
  ord <- function(x) { n <- suppressWarnings(as.numeric(sub(" .*", "", x))); x[order(if (all(!is.na(n))) -n else seq_along(x))] }
  G$lv <- paste(G$gruppo, G$livello, sep = "|"); lv <- unlist(lapply(unique(G$gruppo), function(g) paste(g, ord(G$livello[G$gruppo == g]), sep = "|")))
  G$lv <- factor(G$lv, levels = lv)
  MIS <- c(stability = "Stability", trust = "Trustworthiness", continuity = "Continuity", mesh = "MeSH agreement")
  Lg <- do.call(rbind, lapply(names(MIS), function(m) data.frame(G[, c("gruppo", "lv", "base")], valore = G[[m]], misura = MIS[[m]], base_v = b[[m]])))
  Lg$misura <- factor(Lg$misura, levels = MIS); Lg$gruppo <- factor(Lg$gruppo, levels = unique(G$gruppo))
  Lg$confrontabile <- !(Lg$misura == MIS[["mesh"]] & Lg$gruppo == "Frequency threshold" & !Lg$base)
  p <- ggplot(Lg, aes(valore, lv)) + geom_vline(aes(xintercept = base_v), colour = CHIARO, linewidth = 0.3) +
    geom_point(aes(colour = base, shape = confrontabile), size = 1.5, stroke = 0.45) + facet_grid(gruppo ~ misura, scales = "free", space = "free_y") +
    scale_colour_manual(values = c("FALSE" = GRIGIO, "TRUE" = COL[1]), guide = "none") + scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1), guide = "none") +
    scale_y_discrete(labels = function(x) sub(".*\\|", "", x)) + scale_x_continuous(breaks = function(l) pretty(l, 3)) + labs(x = NULL, y = NULL) + tema + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7), panel.spacing.x = unit(5, "mm"))
  salva(p, "fig8_sensitivity", 180, 125) })
fai("p18_mappa_interattiva.html", {
  D2 <- dist_euclidea(cbind(dati$x, dati$y)); Nv <- vicini_knn(D2, 5)
  js <- function(s) paste0('"', gsub('"', '\\\\"', s), '"')
  righe <- vapply(seq_len(n_T), function(i) with(dati[i, ], sprintf(
    '{"t":%s,"g":%s,"d":%d,"s":%.3f,"a":%.2f,"x":%.4f,"y":%.4f,"mx":%.5f,"my":%.5f,"e1":%.4f,"e2":%.4f,"w":%.4f,"u":%d,"n":[%s]}',
    js(nome), js(as.character(gruppo)), documenti, stabilita, 100 * quota_area, x, y, mx, my, asse1, asse2, angolo, as.integer(incerto),
    paste(js(dati$nome[Nv[i, ]]), collapse = ","))), "")
  cl <- paste(sprintf("%s:%s", js(names(col_g)), js(unname(col_g))), collapse = ",")
  tmpl <- r"---(<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>SeCoS and MCA maps</title>
<style>
:root{--bg:#fcfcfb;--ink:#2b2b28;--mut:#6b6a63;--line:#d4d3cb;--pop:#fff}
@media (prefers-color-scheme:dark){:root{--bg:#1a1a19;--ink:#f0efec;--mut:#b5b4aa;--line:#44443f;--pop:#262624}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:14px/1.45 system-ui,sans-serif}
main{max-width:1100px;margin:0 auto;padding:16px}h1{font-size:20px;margin:0 0 4px}p{margin:0 0 10px;color:var(--mut)}
.bar{display:flex;flex-wrap:wrap;gap:10px;align-items:center;margin:10px 0}
button,input[type=search]{font:inherit;color:inherit;background:var(--pop);border:1px solid var(--line);border-radius:6px;padding:5px 10px}
button[aria-pressed=true]{border-color:var(--ink);font-weight:600}input[type=search]{min-width:200px}label{color:var(--mut)}
</style></head><body><main>
<h1>SeCoS and MCA maps of sports injury keywords</h1>
<p>171 author keywords. Colors are k-means groups on the embeddings, the same in both maps. Hollow markers have an uncertain position on the SeCoS map and are not interpreted. UMAP does not preserve global distances: the separation between distant regions carries no interpretation. Hover over a term for its details, scroll to zoom, drag to move.</p>
<div class="bar"><button id="bS" aria-pressed="true">SeCoS consensus map</button><button id="bM" aria-pressed="false">MCA map</button>
<label><input type="checkbox" id="cE"> 95% ellipses</label><label><input type="checkbox" id="cL" checked> Labels of frequent terms</label>
<input id="q" type="search" placeholder="Find a term" aria-label="Find a term" list="ls"><datalist id="ls"></datalist><button id="bR">Reset view</button></div>
<div id="leg"></div><div id="w"><svg id="s" role="img" aria-label="Map of keywords"></svg><div id="tip"></div></div>
</main><script>
const D=[/*DATI*/],C={/*COLORI*/},T=/*NETICHETTE*/;
const s=document.getElementById('s'),tip=document.getElementById('tip'),w=document.getElementById('w'),NS='http://www.w3.org/2000/svg';
let map='S',view=null,off=new Set(),sel=null;const W=1000,H=700;
const ord=[...D].sort((a,b)=>b.d-a.d),lab=new Set(ord.slice(0,T).map(d=>d.t));
const px=d=>map=='S'?d.x:d.mx,py=d=>map=='S'?d.y:d.my;
function box(){const q=(a,p)=>{a=[...a].sort((x,y)=>x-y);return a[Math.floor(p*(a.length-1))]};
 const xs=D.map(px),ys=D.map(py);let x0=map=='S'?Math.min(...xs):q(xs,.01),x1=map=='S'?Math.max(...xs):q(xs,.99),y0=map=='S'?Math.min(...ys):q(ys,.01),y1=map=='S'?Math.max(...ys):q(ys,.99);
 const mx=(x1-x0)*.08,my=(y1-y0)*.08;x0-=mx;x1+=mx;y0-=my;y1+=my;const k=Math.max((x1-x0)/W,(y1-y0)/H),cx=(x0+x1)/2,cy=(y0+y1)/2;return{cx,cy,k}}
const fx=x=>W/2+(x-view.cx)/view.k,fy=y=>H/2-(y-view.cy)/view.k;
function el(n,a){const e=document.createElementNS(NS,n);for(const k in a)e.setAttribute(k,a[k]);return e}
function draw(){s.setAttribute('viewBox',`0 0 ${W} ${H}`);s.textContent='';const E=document.getElementById('cE').checked&&map=='S',L=document.getElementById('cL').checked;
 const ng=sel?new Set(sel.n):null,vis=D.filter(d=>!off.has(d.g));
 if(E)for(const d of vis){const a=d.e1/view.k,b=d.e2/view.k;s.append(el('ellipse',{cx:fx(d.x),cy:fy(d.y),rx:a,ry:b,transform:`rotate(${-d.w*180/Math.PI} ${fx(d.x)} ${fy(d.y)})`,fill:C[d.g],'fill-opacity':.08,stroke:d.u?'#8c8b83':'none','stroke-dasharray':'4 3'}))}
 for(const d of vis){const on=!sel||d==sel||ng.has(d.t),r=2.5+Math.sqrt(d.d)/9,hollow=d.u&&map=='S',c=el('circle',{cx:fx(px(d)),cy:fy(py(d)),r:r,fill:hollow?'none':C[d.g],stroke:C[d.g],'stroke-width':hollow?2:1,opacity:on?.95:.15,tabindex:0,'aria-label':d.t});
  c.style.pointerEvents='all';c.addEventListener('pointerenter',e=>show(d,e));c.addEventListener('pointerleave',()=>tip.style.display='none');c.addEventListener('focus',e=>show(d,e));c.addEventListener('click',()=>{sel=sel==d?null:d;draw()});s.append(c)}
 for(const d of vis){const on=sel?(d==sel||ng.has(d.t)):(L&&lab.has(d.t));if(on){const t=el('text',{x:fx(px(d))+6,y:fy(py(d))+4});t.textContent=d.t;t.style.pointerEvents='none';s.append(t)}}}
function show(d,e){tip.innerHTML=`<b></b><span></span>`;tip.firstChild.textContent=d.t;tip.lastChild.textContent=`${d.g}. ${d.d} documents. Stability ${d.s.toFixed(2)}. Ellipse ${d.a.toFixed(1)}% of the map${d.u&&map=='S'?' (uncertain position)':''}. Nearest on the SeCoS map: ${d.n.join(', ')}.`;
 const r=w.getBoundingClientRect(),b=(e.target||s).getBoundingClientRect();tip.style.display='block';let l=b.left-r.left+12,t=b.top-r.top+12;if(l+290>r.width)l=Math.max(4,l-300);if(t+130>r.height)t=Math.max(4,t-130);tip.style.left=l+'px';tip.style.top=t+'px'}
function setM(m){map=m;view=box();sel=null;document.getElementById('bS').setAttribute('aria-pressed',m=='S');document.getElementById('bM').setAttribute('aria-pressed',m=='M');draw()}
document.getElementById('bS').onclick=()=>setM('S');document.getElementById('bM').onclick=()=>setM('M');document.getElementById('bR').onclick=()=>{view=box();sel=null;draw()};
document.getElementById('cE').onchange=draw;document.getElementById('cL').onchange=draw;
const leg=document.getElementById('leg');for(const g in C){const b=document.createElement('button');b.setAttribute('aria-pressed','true');const i=document.createElement('i');i.style.background=C[g];b.append(i,g);b.onclick=()=>{off.has(g)?off.delete(g):off.add(g);b.setAttribute('aria-pressed',!off.has(g));draw()};leg.append(b)}
const dl=document.getElementById('ls');for(const d of [...D].sort((a,b)=>a.t.localeCompare(b.t))){const o=document.createElement('option');o.value=d.t;dl.append(o)}
document.getElementById('q').oninput=e=>{const d=D.find(x=>x.t==e.target.value);if(d){sel=d;draw()}};
s.addEventListener('wheel',e=>{e.preventDefault();const r=s.getBoundingClientRect(),f=e.deltaY>0?1.15:1/1.15,sc=Math.max(W/r.width,H/r.height),mx=(e.clientX-r.left-r.width/2)*sc,my=(e.clientY-r.top-r.height/2)*sc;
 view.cx+=mx*view.k*(1-f);view.cy-=my*view.k*(1-f);view.k*=f;draw()},{passive:false});
let dr=null;s.addEventListener('pointerdown',e=>{if(e.target==s)dr=[e.clientX,e.clientY]});window.addEventListener('pointerup',()=>dr=null);
window.addEventListener('pointermove',e=>{if(!dr)return;const r=s.getBoundingClientRect(),sc=Math.max(W/r.width,H/r.height);view.cx-=(e.clientX-dr[0])*sc*view.k;view.cy+=(e.clientY-dr[1])*sc*view.k;dr=[e.clientX,e.clientY];draw()});
setM('S');
</script></body></html>)---"
  html <- sub("/*DATI*/", paste(righe, collapse = ","), tmpl, fixed = TRUE)
  html <- sub("/*COLORI*/", cl, html, fixed = TRUE); html <- sub("/*NETICHETTE*/", sum(dati$etichetta), html, fixed = TRUE)
  con <- file("results/interactive_map.html", open = "wb"); writeLines(enc2utf8(html), con, useBytes = TRUE); close(con) })
DEF <- readRDS("wos/corpus_clean.rds")
if (!"DB" %in% names(DEF) || all(is.na(DEF$DB))) DEF$DB <- "ISI"
STOP <- readLines("data/stopwords.txt", encoding = "UTF-8")
E0 <- ENC$sbert[TERMINI, ]
x <- secos_map(E0, n_neighbors = 10, min_dist = 0.1, replicates = 100L)
m <- x$map[match(TERMINI, x$map$term), ]
scrivi_csv(data.frame(max_diff_coordinates = max(abs(m$x - cons$x), abs(m$y - cons$y)), max_diff_axes = max(abs(m$axis1 - cons$axis1), abs(m$axis2 - cons$axis2)),
                      max_diff_stability = max(abs(m$stability - cons$stability)), same_uncertain_terms = identical(as.logical(m$uncertain), ctrl$incerto)),
           "results/secos_check.csv")
rm(x); invisible(gc())

SOGLIE <- c(100, 50, 30, 15, 8, 4); R_MAPPA <- 100L; N_RIPETI <- 3L; SOTTO_S <- 60; LIMITE_MIN <- 30; LIMITE_SECOS <- 60; FATTORE_MEM <- 13
PASSI <- c("neighbours", "start", "replicates", "alignment", "median", "ellipses", "zones", "stability")
ram <- tryCatch(if (.Platform$OS.type == "windows") round(as.numeric(system2("powershell", c("-NoProfile", "-Command", "(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory"), stdout = TRUE)[1]) / 1024^3, 1) else
  round(as.numeric(gsub("[^0-9]", "", grep("MemTotal", readLines("/proc/meminfo"), value = TRUE))) / 1024^2, 1), error = function(e) NA_real_, warning = function(w) NA_real_)
cpu <- tryCatch(if (.Platform$OS.type == "windows") system2("powershell", c("-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"), stdout = TRUE)[1] else
  sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]), error = function(e) Sys.getenv("PROCESSOR_IDENTIFIER"), warning = function(w) Sys.getenv("PROCESSOR_IDENTIFIER"))
writeLines(c(cpu, paste(ram, "GB"), R.version.string, sessionInfo()$BLAS, sessionInfo()$LAPACK), "results/timing_machine.txt")
LIMITE_GB <- if (is.na(ram)) 6 else 0.7 * ram
if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) reticulate::use_condaenv("r-text", required = TRUE)
kw <- .secos_keywords(DEF, "DE", ";")
freq <- secos_terms(DEF, "DE", min(SOGLIE), STOP)
matrice_doc <- function(tt) {
  pos <- lapply(kw, function(v) { p <- unique(match(v, tt)); p[!is.na(p)] })
  X <- Matrix::sparseMatrix(i = rep(seq_along(pos), lengths(pos)), j = unlist(pos), x = 1, dims = c(length(pos), length(tt)), dimnames = list(NULL, tt))
  as.matrix(X[Matrix::rowSums(X) > 0, , drop = FALSE])
}
corpus_ridotto <- function(tt) {
  D2 <- DEF
  D2$DE <- vapply(kw, function(v) paste(unique(v[v %in% tt]), collapse = "; "), "")
  D2$DE <- toupper(D2$DE); D2$DE[!nzchar(D2$DE)] <- NA
  D2
}
Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
misura <- function(compito, dati) {
  f_in <- tempfile(fileext = ".rds"); f_out <- tempfile(fileext = ".rds"); saveRDS(dati, f_in)
  st <- system2(file.path(R.home("bin"), "Rscript"), c("--no-init-file", "R/timing_task.R", compito, f_in, f_out), stdout = FALSE, stderr = FALSE)
  if (st != 0 || !file.exists(f_out)) stop("timing task '", compito, "' failed")
  r <- readRDS(f_out); unlink(c(f_in, f_out)); r
}
ripeti <- function(f) {
  r <- list(f())
  if (sum(r[[1]]) < SOTTO_S) for (i in seq_len(N_RIPETI - 1)) r[[i + 1]] <- f()
  list(tempi = apply(do.call(cbind, r), 1, stats::median), n = length(r))
}
stima <- function(tt, ss) {
  ok <- is.finite(ss); tt <- tt[ok]; ss <- ss[ok]; n <- length(ss)
  if (n == 0) return(function(t) 0)
  e <- if (n >= 2) min(max(log(ss[n] / ss[n - 1]) / log(tt[n] / tt[n - 1]), 1), 3) else 2
  function(t) ss[n] * (t / tt[n])^e
}
RIS <- list()
for (s in SOGLIE) {
  tt <- freq$term[freq$documents >= s]; nT <- length(tt)
  fatti <- RIS[order(vapply(RIS, `[[`, 0, "T"))]; Tv <- vapply(fatti, `[[`, 0, "T")
  if (stima(Tv, vapply(fatti, function(z) z$s_tot, 0))(nT) / 60 > LIMITE_SECOS) next
  E_T <- NULL
  r_e <- ripeti(function() { z <- secos_embed(tt, threads = 1L); E_T <<- z; attr(z, "timing") })
  r_m <- ripeti(function() misura("secos", list(E = E_T, k = 10, md = 0.1, R = R_MAPPA)))
  B <- matrice_doc(tt); gb <- FATTORE_MEM * nrow(B) * 2 * nT * 8 / 1e9
  att <- function(k) stima(Tv, vapply(fatti, function(z) z$altri[[k]], 0))(nT)
  salta <- function(k) { if (gb > LIMITE_GB && k != "ca") return(sprintf("memory %.1f GB", gb))
    if (att(k) / 60 > LIMITE_MIN) return(sprintf("expected %.0f min", att(k) / 60)); "" }
  altri <- c(mca = NA_real_, cs = NA_real_, ca = NA_real_); nota <- c(mca = "", cs = "", ca = "")
  esegui <- function(k, dati) {
    nota[[k]] <<- salta(k); if (nzchar(nota[[k]])) return(invisible())
    altri[[k]] <<- tryCatch(ripeti(function() misura(k, dati))$tempi, error = function(e) { nota[[k]] <<- conditionMessage(e); NA_real_ })
  }
  esegui("ca", list(B = B)); esegui("mca", list(B = B))
  esegui("cs", list(M = corpus_ridotto(tt), s = s, stop = STOP))
  RIS[[as.character(s)]] <- list(soglia = s, T = nT, doc = nrow(B), s_tot = sum(r_e$tempi) + sum(r_m$tempi[c(PASSI, "fidelity")]), e = r_e$tempi, m = r_m$tempi, altri = altri, nota = nota)
  rm(B, E_T); invisible(gc())
}
RIS <- RIS[order(vapply(RIS, `[[`, 0, "T"))]
tab <- do.call(rbind, lapply(RIS, function(z) { m <- z$m; e <- z$e
  data.frame(threshold = z$soglia, T = z$T, documents = z$doc, load = e[["load"]], embedding = e[["encode"]], neighbours = m[["neighbours"]] + m[["start"]],
             umap_100 = m[["replicates"]], umap_one = m[["replicates"]] / R_MAPPA, alignment = m[["alignment"]] + m[["median"]],
             ellipses_zones = sum(m[c("ellipses", "zones")]), stability = m[["stability"]],
             secos_total = e[["load"]] + e[["encode"]] + sum(m[c("neighbours", "start", "replicates", "alignment", "median", "ellipses", "zones")]) + m[["stability"]],
             secos_one_map = e[["encode"]] + m[["neighbours"]] + m[["replicates"]] / R_MAPPA,
             mca_step = z$altri[["mca"]], mca_function = z$altri[["cs"]], ca = z$altri[["ca"]], note_mca = z$nota[["mca"]], note_function = z$nota[["cs"]]) }))
rownames(tab) <- NULL
pendenza <- function(v) { i <- tab$T >= 171 & is.finite(v); if (sum(i) < 3) return(NA_real_); unname(coef(lm(log(v[i]) ~ log(tab$T[i])))[2]) }
scrivi_csv(tab, "results/table7_timing.csv")
scrivi_csv(data.frame(measure = c("embedding", "umap_100", "secos_total", "mca_step", "mca_function", "ca"),
                      slope = vapply(c("embedding", "umap_100", "secos_total", "mca_step", "mca_function", "ca"), function(n) pendenza(tab[[n]]), 0)),
           "results/timing_slopes.csv")
sessione("10")
