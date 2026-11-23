source("R/functions.R")
RNGkind("L'Ecuyer-CMRG")
A <- readRDS("wos/corpus_broad.rds")
UT_RIS <- readLines("results/ut_proximity_corpus.txt")
UT_DEF <- readLines("results/ut_final_corpus.txt")
th <- leggi_csv("data/thesaurus.csv", na.strings = character(0), colClasses = c(author_keyword = "character", term = "character"))
sv <- leggi_csv("data/acronyms.csv", na.strings = "")
sin <- leggi_csv("data/synonyms.csv")
STOPWORD <- readLines("data/stopwords.txt", encoding = "UTF-8")
TERMINI_T <- readLines("results/terms_threshold_50.txt", encoding = "UTF-8")
SOGLIA_T <- 50; ANNI <- 1989:2025; INIZIO <- 1991; N_AVVII <- 10; K <- 10
RE_OMBRELLA <- "\\bSYSTEM\\b|\\bAFFILIATES\\b"
CORPUS <- c("Ampio", "Ristretto", "Definitivo")

if (!"DB" %in% names(A) || all(is.na(A$DB))) A$DB <- "ISI"
if (!"AU1_CO" %in% names(A)) A <- metaTagExtraction(A, Field = "AU1_CO", sep = ";")
au1 <- trimws(A$AU1_CO); au1[!nzchar(au1)] <- NA
A <- A[, intersect(c("UT", "PY", "DE", "AU", "C3", "TC", "TI", "AB", "SO", "J9", "DI"), names(A)), drop = FALSE]
invisible(gc())
G <- list(Ampio = rep(TRUE, nrow(A)), Ristretto = A$UT %in% UT_RIS, Definitivo = A$UT %in% UT_DEF)
dividi <- function(x) {
  p <- strsplit(ifelse(is.na(x), "", x), ";", fixed = TRUE)
  v <- trimws(unlist(p)); id <- rep(seq_along(p), lengths(p)); ok <- nzchar(v)
  list(id = id[ok], v = v[ok])
}
una_volta <- function(l) { d <- !duplicated(paste(l$id, l$v)); list(id = l$id[d], v = l$v[d]) }

UK_US <- c(paediatric = "pediatric", paediatrics = "pediatrics", orthopaedic = "orthopedic",
           orthopaedics = "orthopedics", behaviour = "behavior", behavioural = "behavioral",
           programme = "program", programmes = "program", centre = "center", oedema = "edema",
           anaesthesia = "anesthesia", haemorrhage = "hemorrhage", ischaemia = "ischemia",
           ischaemic = "ischemic", tumour = "tumor", fibre = "fiber", fibres = "fiber",
           randomised = "randomized", haematoma = "hematoma", oestrogen = "estrogen",
           aetiology = "etiology", haemoglobin = "hemoglobin", labour = "labor")
NON_PLURALI <- c("diabetes", "series", "species", "biceps", "triceps", "quadriceps", "lens",
                 "pilates", "achilles", "mumps", "measles", "rickets", "news", "means",
                 "sports medicine", "status", "pes", "physics", "herpes", "rabies")
base_norm <- function(x) {
  x <- tolower(x); x <- gsub("'", "", x, fixed = TRUE)
  x <- gsub("[^[:alnum:]()]+", " ", x); trimws(gsub(" +", " ", x))
}
estrai_sigla <- function(x) {
  m <- regmatches(x, regexec("^(.+?) ?\\(([[:alnum:] ]+)\\)$", x))
  t(vapply(m, function(v) if (length(v) == 3) c(trimws(v[2]), trimws(v[3])) else c(NA_character_, NA_character_), character(2)))
}
singolare <- function(w) {
  out <- w
  ok <- !(w %in% NON_PLURALI) & !grepl("(ss|us|is|ics|sis|ness)$", w) & grepl("[aeiouy]", w)
  i1 <- ok & grepl("ies$", w) & nchar(w) > 4;  out[i1] <- sub("ies$", "y", w[i1])
  i2 <- ok & !i1 & grepl("(ches|shes|xes)$", w); out[i2] <- sub("es$", "", w[i2])
  i3 <- ok & !i1 & !i2 & grepl("s$", w) & nchar(w) > 3; out[i3] <- sub("s$", "", w[i3])
  out
}
norm_parole <- function(x) {
  vapply(strsplit(x, " "), function(p) {
    if (length(p) == 0) return("")
    p <- ifelse(p %in% names(UK_US), UK_US[p], p)
    p <- gsub("isation$", "ization", p)
    frase <- paste(p, collapse = " ")
    if (frase %in% NON_PLURALI) return(frase)
    p[length(p)] <- singolare(p[length(p)])
    paste(p, collapse = " ")
  }, character(1))
}
sigla_plausibile <- function(s, e) {
  mapply(function(s, e) {
    if (is.na(s) || is.na(e) || !nzchar(e)) return(FALSE)
    s <- trimws(s); e <- trimws(e)
    if (!nzchar(s) || !nzchar(e)) return(FALSE)
    cs <- strsplit(s, "")[[1]]; ce <- strsplit(gsub(" ", "", e), "")[[1]]
    if (length(cs) == 0 || length(ce) == 0) return(FALSE)
    if (cs[1] != ce[1] || s %in% strsplit(e, " ")[[1]] || nchar(s) >= nchar(e)) return(FALSE)
    pos <- 0
    for (ch in cs) {
      if (pos >= length(ce)) return(FALSE)
      k <- which(ce[(pos + 1):length(ce)] == ch)
      if (length(k) == 0) return(FALSE); pos <- pos + k[1] }
    TRUE
  }, s, e, USE.NAMES = FALSE)
}
espandi <- function(x, mappa) {
  if (nrow(mappa) == 0) return(x)
  vapply(strsplit(x, " "), function(p) {
    j <- match(p, mappa$sigla); p[!is.na(j)] <- mappa$estesa[j[!is.na(j)]]
    paste(p, collapse = " ")
  }, character(1))
}
kw <- dividi(A$DE)
u_raw <- unique(kw$v)
b <- base_norm(u_raw)
sg <- estrai_sigla(b)
plaus <- sigla_plausibile(sg[, 2], ifelse(is.na(sg[, 1]), NA, norm_parole(sg[, 1])))
dup <- !is.na(sg[, 2]) & !is.na(sg[, 1]) & norm_parole(ifelse(is.na(sg[, 2]), "", sg[, 2])) == norm_parole(ifelse(is.na(sg[, 1]), "", sg[, 1]))
senza_par <- trimws(gsub(" +", " ", ifelse(plaus | dup, sg[, 1], gsub("\\(|\\)", " ", b))))
chiave1 <- norm_parole(senza_par)
sig_ok <- data.frame(sigla = sv$acronym[sv$accepted == 1], estesa = sv$final_expansion[sv$accepted == 1], stringsAsFactors = FALSE)
corr <- sv[!is.na(sv$corrected_expansion) & sv$corrected_expansion != sv$proposed_expansion, c("proposed_expansion", "corrected_expansion")]
ch <- espandi(chiave1, sig_ok)
j <- match(ch, corr$proposed_expansion); ch[!is.na(j)] <- corr$corrected_expansion[j[!is.na(j)]]
for (it in 1:5) { j <- match(ch, sin$term); if (all(is.na(j))) break; ch[!is.na(j)] <- sin$merged_into[j[!is.na(j)]] }
j <- match(u_raw, th$author_keyword); ch[!is.na(j)] <- th$term[j[!is.na(j)]]
term <- ch[match(kw$v, u_raw)]
ok <- nzchar(term)
kw_t <- una_volta(list(id = kw$id[ok], v = term[ok]))
u_t <- unique(kw_t$v)
B <- Matrix::sparseMatrix(i = kw_t$id, j = match(kw_t$v, u_t), x = 1, dims = c(nrow(A), length(u_t)), dimnames = list(NULL, u_t))
rm(kw, kw_t, term); invisible(gc())

tokens_campo <- function(x) {
  p <- strsplit(gsub("[^a-z]+", " ", tolower(ifelse(is.na(x), "", x))), " ", fixed = TRUE)
  t <- unlist(p); d <- rep(seq_along(p), lengths(p)); pos <- sequence(lengths(p))
  ok <- nzchar(t); list(t = t[ok], d = d[ok], pos = pos[ok])
}
SPORT_ESATTI <- c("football", "footballer", "footballers", "soccer", "rugby", "hockey", "basketball", "baseball", "softball",
                  "handball", "volleyball", "tennis", "cricket", "netball", "lacrosse", "futsal", "skiing", "skier", "skiers",
                  "snowboarding", "snowboarder", "snowboarders", "dancer", "dancers", "swimmer", "swimmers", "cyclist",
                  "cyclists", "boxer", "boxers", "boxing", "runner", "runners", "military", "soldier", "soldiers")
SPORT_PREF <- c("athlet", "player", "gymnast", "wrestl", "triathl", "marathon", "golf")
INF_ESATTI <- c("reinjury", "reinjuries", "preinjury", "strain", "strains", "tear", "tears", "acl", "tendinopathy", "tendinopathies")
INF_PREF <- c("concuss", "sprain", "ruptur", "dislocat", "fractur")
r3_campo <- function(x) {
  z <- tokens_campo(x); t <- z$t; u <- unique(t)
  pref <- function(v) Reduce(`|`, lapply(v, function(p) startsWith(u, p)))
  S1 <- pref(c("sport", "exercise")) | pref(SPORT_PREF) | u %in% SPORT_ESATTI
  I1 <- u %in% c("injury", "injuries") | u %in% INF_ESATTI | pref(INF_PREF)
  jj <- match(t, u); nxt <- c(t[-1], "")
  ia <- which(I1[jj] | (t == "anterior" & nxt == "cruciate")); ib <- which(S1[jj] | (t == "physical" & startsWith(nxt, "activit")))
  docs <- intersect(z$d[ia], z$d[ib]); out <- rep(FALSE, nrow(A))
  if (length(docs) == 0) return(out)
  pa <- split(z$pos[ia], z$d[ia]); pb <- split(z$pos[ib], z$d[ib]); chd <- as.character(docs)
  out[docs] <- vapply(chd, function(x) any(abs(outer(pa[[x]], pb[[x]], "-")) <= K), logical(1), USE.NAMES = FALSE)
  out
}
R3 <- r3_campo(A$TI) | r3_campo(A$AB) | r3_campo(gsub(";", strrep(" zz", K + 1), A$DE, fixed = TRUE))
B3 <- readLines("results/exclusion_terms.txt")
tiene_b3 <- !grepl(paste(B3, collapse = "|"), tolower(paste(A$TI, A$AB, A$DE)), perl = TRUE)
G$R3 <- (G$Ristretto | R3) & tiene_b3

freq <- lapply(G, function(g) { f <- Matrix::colSums(B[g, , drop = FALSE]); f[f > 0] })
con_de <- vapply(G, function(g) sum(Matrix::rowSums(B[g, , drop = FALSE]) > 0), numeric(1))
fr <- lapply(freq, function(f) { f <- f[!names(f) %in% STOPWORD]; f[order(-f, names(f))] })
T_soglia <- lapply(fr, function(f) names(f)[f >= SOGLIA_T]); T_soglia$Definitivo <- TERMINI_T
T_top <- lapply(fr, function(f) names(f)[seq_len(min(length(TERMINI_T), length(f)))]); T_top$Definitivo <- TERMINI_T
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))

py <- suppressWarnings(as.integer(A$PY))
crescita <- function(g) {
  n <- as.integer(table(factor(py[g], levels = ANNI)))
  d <- data.frame(anno = ANNI, documenti = n)
  m <- glm(documenti ~ anno, family = quasipoisson, data = d[d$anno >= INIZIO, ])
  b <- coef(m)[["anno"]]; se <- sqrt(vcov(m)["anno", "anno"]); q <- qt(0.975, df.residual(m))
  c(composto = (n[length(n)] / n[ANNI == INIZIO])^(1 / (max(ANNI) - INIZIO)) - 1,
    qp = exp(b) - 1, qp_inf = exp(b - q * se) - 1, qp_sup = exp(b + q * se) - 1, quota_2021_25 = sum(n[ANNI >= 2021]) / sum(n))
}
cr <- lapply(G, crescita)
au <- dividi(toupper(A$AU))
au <- una_volta(list(id = au$id[!au$v %in% c("ANONYMOUS", "[ANONYMOUS]", "NA")], v = au$v[!au$v %in% c("ANONYMOUS", "[ANONYMOUS]", "NA")]))
lotka <- function(g) {
  x <- as.integer(table(au$v[g[au$id]]))
  m <- poweRlaw::displ$new(x); m$setXmin(1); m$setPars(poweRlaw::estimate_pars(m))
  c(autori = length(x), alfa = m$getPars())
}
lo <- lapply(G, lotka)
conta <- function(x) sort(table(x), decreasing = TRUE)
quanti <- function(t, nome) as.numeric(t)[match(nome, toupper(names(t)))]
paesi <- lapply(G, function(g) conta(au1[g]))
sovr <- function(rif, altro, n) length(intersect(head(names(rif), n), head(names(altro), n)))
rho_top <- function(rif, altro, n) {
  v <- head(names(rif), n); b <- as.numeric(altro)[match(v, names(altro))]; b[is.na(b)] <- 0
  suppressWarnings(cor(as.numeric(rif[v]), b, method = "spearman"))
}
ist <- leggi_csv("data/institutions.csv", na.strings = "")
escl <- ist$institution[ist$excluded]
mappa <- setNames(ist$merged_into, ist$institution)[!is.na(ist$merged_into)]
tenute <- unique(ifelse(is.na(ist$merged_into), ist$institution, ist$merged_into)[ist$reviewed & !ist$excluded])
af <- dividi(A$C3)
j <- match(af$v, names(mappa)); af$v[!is.na(j)] <- mappa[j[!is.na(j)]]
ok <- !af$v %in% escl & (af$v %in% tenute | !grepl(RE_OMBRELLA, toupper(af$v)))
af <- una_volta(list(id = af$id[ok], v = af$v[ok]))
istit <- lapply(G, function(g) conta(af$v[g[af$id]]))
TC <- suppressWarnings(as.numeric(A$TC)); TC[is.na(TC)] <- 0
primi <- lapply(G, function(g) { i <- which(g); i[order(-TC[i], A$UT[i])][seq_len(min(20, length(i)))] })

louvain <- function(termini, g, seed0 = SEED) {
  termini <- intersect(termini, colnames(B))
  C <- as.matrix(Matrix::crossprod(B[g, termini, drop = FALSE])); f <- diag(C)
  Adj <- C / outer(f, f); diag(Adj) <- 0
  gr <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = TRUE, diag = FALSE)
  migliore <- NULL
  for (r in seq_len(N_AVVII)) {
    set.seed(seed0 + r - 1); cl <- igraph::cluster_louvain(gr, resolution = 1)
    if (is.null(migliore) || igraph::modularity(cl) > igraph::modularity(migliore)) migliore <- cl
  }
  list(m = setNames(as.integer(igraph::membership(migliore)), termini), q = igraph::modularity(migliore))
}
accordo <- function(p, q) { t <- intersect(names(p$m), names(q$m)); c(comuni = length(t), ARI = ari(p$m[t], q$m[t]), NMI = nmi(p$m[t], q$m[t])) }
L_soglia <- mapply(louvain, T_soglia[CORPUS], G[CORPUS], SIMPLIFY = FALSE)
L_top <- mapply(louvain, T_top, G, SIMPLIFY = FALSE)
L_rif <- louvain(T_soglia$Definitivo, G$Definitivo, seed0 = SEED + 1000)

vs <- function(f) vapply(CORPUS, f, numeric(1))
righe <- list(
  "Documents"                                         = vs(function(k) sum(G[[k]])),
  "Documents with author keywords"                    = vs(function(k) con_de[[k]]),
  "Compound growth rate 1991-2025 (%)"                = vs(function(k) 100 * cr[[k]][["composto"]]),
  "Quasi-Poisson growth rate (%)"                     = vs(function(k) 100 * cr[[k]][["qp"]]),
  "  95% CI, lower"                                   = vs(function(k) 100 * cr[[k]][["qp_inf"]]),
  "  95% CI, upper"                                   = vs(function(k) 100 * cr[[k]][["qp_sup"]]),
  "Share of documents 2021-2025 (%)"                  = vs(function(k) 100 * cr[[k]][["quota_2021_25"]]),
  "Citations per document"                            = vs(function(k) mean(TC[G[[k]]])),
  "Distinct authors"                                  = vs(function(k) lo[[k]][["autori"]]),
  "Lotka exponent"                                    = vs(function(k) lo[[k]][["alfa"]]),
  "Share of the USA (%)"                              = vs(function(k) 100 * quanti(paesi[[k]], "USA") / sum(paesi[[k]])),
  "Top 10 countries in common with the final corpus"  = vs(function(k) sovr(paesi$Definitivo, paesi[[k]], 10)),
  "Spearman correlation, top 20 countries"            = vs(function(k) rho_top(paesi$Definitivo, paesi[[k]], 20)),
  "Top 20 institutions in common"                     = vs(function(k) sovr(istit$Definitivo, istit[[k]], 20)),
  "Spearman correlation, top 20 institutions"         = vs(function(k) rho_top(istit$Definitivo, istit[[k]], 20)),
  "Top 20 documents by citations still in the final corpus" = vs(function(k) sum(G$Definitivo[primi[[k]]])),
  "Terms with at least 50 documents"                  = vs(function(k) length(T_soglia[[k]])),
  "Top 171 terms in common"                           = vs(function(k) length(intersect(T_top[[k]], T_top$Definitivo))),
  "Jaccard index, top 171 terms"                      = vs(function(k) jacc(T_top[[k]], T_top$Definitivo)),
  "Louvain communities, top 171 terms"                = vs(function(k) length(unique(L_top[[k]]$m))),
  "ARI with the final corpus, top 171 terms"          = vs(function(k) accordo(L_top[[k]], L_top$Definitivo)[["ARI"]]),
  "NMI with the final corpus, top 171 terms"          = vs(function(k) accordo(L_top[[k]], L_top$Definitivo)[["NMI"]]),
  "ARI with the final corpus, threshold 50"           = vs(function(k) accordo(L_soglia[[k]], L_soglia$Definitivo)[["ARI"]]),
  "NMI with the final corpus, threshold 50"           = vs(function(k) accordo(L_soglia[[k]], L_soglia$Definitivo)[["NMI"]]),
  "Documents without a recognizable country"          = vs(function(k) sum(is.na(au1[G[[k]]]))))
tab <- data.frame(measure = names(righe), do.call(rbind, righe), row.names = NULL, check.names = FALSE)
names(tab)[2:4] <- c("Broad", "Proximity", "Final")
scrivi_csv(tab, "results/table_S2.csv")
rif <- accordo(L_soglia$Definitivo, L_rif)
a3 <- accordo(L_top$R3, L_top$Definitivo)
scrivi_csv(data.frame(item = c("reference ARI, final corpus with other seeds", "reference NMI, final corpus with other seeds",
                               "broader rule: records", "broader rule: top 171 terms in common", "broader rule: Jaccard index", "broader rule: ARI"),
                      value = c(rif[["ARI"]], rif[["NMI"]], sum(G$R3), length(intersect(T_top$R3, T_top$Definitivo)),
                                jacc(T_top$R3, T_top$Definitivo), a3[["ARI"]])), "results/table_S2_notes.csv")
sessione("05")
