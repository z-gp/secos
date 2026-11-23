source("R/functions.R")
source("R/secos.R")
Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE", TOKENIZERS_PARALLELISM = "false")
if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) reticulate::use_condaenv("r-text", required = TRUE)
MOD <- leggi_csv("data/models.csv")
modello <- function(e) MOD$model[MOD$encoder == e]
revisione <- function(e) MOD$revision[MOD$encoder == e]
N_NEIGHBORS <- 10; MIN_DIST <- 0.1; K_MAPPA <- 10; K_FED <- 15; SEMI_UMAP <- 1:20; KS <- c(5, 10, 20)

DEF <- readRDS("wos/corpus_clean.rds")
if (!"DB" %in% names(DEF) || all(is.na(DEF$DB))) DEF$DB <- "ISI"
STOPWORD <- readLines("data/stopwords.txt", encoding = "UTF-8")
TERMINI_P06 <- readLines("results/terms_threshold_50.txt", encoding = "UTF-8")
FREQ <- leggi_csv("results/term_frequencies.csv")
cw <- cocMatrix(DEF, Field = "DE", type = "sparse", sep = ";", binary = TRUE, remove.terms = STOPWORD)
fr_cw <- Matrix::colSums(cw)
T_SEL <- setdiff(tolower(trimws(colnames(cw)[fr_cw >= 50])), c("na", ""))
rm(cw); invisible(gc())
D2 <- DEF
D2$DE <- vapply(strsplit(ifelse(is.na(DEF$DE), "", DEF$DE), ";"), function(v) { v <- trimws(v); paste(v[tolower(v) %in% T_SEL], collapse = "; ") }, character(1))
D2$DE[!nzchar(D2$DE)] <- NA
arg <- list(M = D2, field = "DE", method = "MCA", minDegree = 50, clust = "auto", k.max = 3, stemming = FALSE,
            remove.terms = STOPWORD, documents = 5, graph = FALSE)
res <- do.call(conceptualStructure, arg[names(arg) %in% names(formals(conceptualStructure))])
rm(D2); invisible(gc())
nomi <- names(res$km.res$cluster)
T_CS <- intersect(trimws(tolower(nomi)), trimws(tolower(colnames(res$net))))
TERMINI <- c(TERMINI_P06[TERMINI_P06 %in% T_CS], setdiff(T_CS, TERMINI_P06))
Y_mca <- as.matrix(res$km.res$data[nomi, 1:2, drop = FALSE])
Y_mca <- Y_mca[match(TERMINI, trimws(tolower(nomi))), ]
scrivi_csv(data.frame(term = TERMINI, documents = FREQ$documents[match(TERMINI, FREQ$term)]), "results/terms.csv")
scrivi_csv(data.frame(term = TERMINI, Dim1 = Y_mca[, 1], Dim2 = Y_mca[, 2]), "results/mca_coordinates.csv")
co <- secos_cooccurrence(DEF, TERMINI)
scrivi_csv(data.frame(documents = nrow(res$net), cooccurrence_documents = co$documents, inertia_dim1 = res$res$inertia.e[1], inertia_dim2 = res$res$inertia.e[2]), "results/mca_summary.csv")
scrivi_csv(data.frame(term = TERMINI, co$counts, check.names = FALSE), "results/cooccurrence_matrix.csv")
rm(res); invisible(gc())

E_sci <- secos_embed(TERMINI, modello("scibert"), revisione("scibert"), layer = -2L, pooling = "term", threads = 1L)
E_sb  <- secos_embed(TERMINI, modello("sbert"), revisione("sbert"), layer = -1L, pooling = "all", threads = 4L)
E_pm  <- secos_embed(TERMINI, modello("pubmed"), revisione("pubmed"), layer = -2L, pooling = "term", threads = 4L)
extra <- setdiff(FREQ$term[FREQ$documents >= 30], TERMINI)
E_x   <- secos_embed(extra, modello("sbert"), revisione("sbert"), layer = -1L, pooling = "all")

forma <- function(x) {
  x <- tolower(x); x <- gsub("\\([^)]*\\)", " ", x); x <- gsub("'", "", x, fixed = TRUE)
  x <- gsub("[^[:alnum:]]+", " ", x); trimws(gsub(" +", " ", x))
}
th <- leggi_csv("data/thesaurus.csv", na.strings = character(0), colClasses = c(author_keyword = "character", term = "character"))
th <- th[th$term %in% TERMINI, ]
FORME <- lapply(TERMINI, function(t) { f <- unique(c(forma(t), forma(th$author_keyword[th$term == t]))); f[nzchar(f) & (grepl(" ", f) | nchar(f) >= 3)] })
AB <- as.character(DEF$AB); AB <- AB[!is.na(AB) & nzchar(trimws(AB))]
reticulate::py_run_string(r"---(
import os, re
from collections import Counter
import numpy as np, torch, transformers, tokenizers
from transformers import AutoTokenizer, AutoModel

def carica(modello, revisione, n_thread):
    torch.set_num_threads(int(n_thread))
    tok = AutoTokenizer.from_pretrained(modello, revision=revisione)
    if not tok.is_fast: raise RuntimeError("a fast tokenizer is needed: " + modello)
    mod = AutoModel.from_pretrained(modello, revision=revisione, output_hidden_states=True)
    mod.eval()
    return tok, mod

VERSIONI = dict(transformers=transformers.__version__, tokenizers=tokenizers.__version__, torch=torch.__version__)
PIANO = {}

def _indice(testo):
    voc, ids, ini, fin = {}, [], [], []
    for m in re.finditer(r"[a-z0-9]+", testo):
        ids.append(voc.setdefault(m.group(0), len(voc))); ini.append(m.start()); fin.append(m.end())
    ids, ini, fin = np.array(ids), np.array(ini), np.array(fin)
    gap = np.array([testo[fin[i]:ini[i + 1]].strip(" -") == "" for i in range(len(ids) - 1)] + [False])
    ord_ = np.argsort(ids, kind="stable")
    tagli = np.cumsum(np.bincount(ids, minlength=len(voc)))[:-1]
    pos = np.split(ord_, tagli)
    return voc, ids, ini, fin, gap, pos

def _ultime(u):
    if u.endswith("y") and len(u) > 3: return [u, u[:-1] + "ies"]
    if u.endswith("s"): return [u]
    return [u, u + "s", u + "es"]

def _cerca(forma, voc, ids, gap, pos):
    w = forma.split(" "); n = len(w)
    ult = [voc[x] for x in _ultime(w[-1]) if x in voc]
    if not ult: return np.zeros(0, int)
    if n == 1: return np.concatenate([pos[i] for i in ult])
    if w[0] not in voc: return np.zeros(0, int)
    c = pos[voc[w[0]]]; c = c[c + n - 1 < len(ids)]
    for j in range(1, n):
        if j < n - 1:
            if w[j] not in voc: return np.zeros(0, int)
            ok = ids[c + j] == voc[w[j]]
        else:
            ok = np.isin(ids[c + j], ult)
        c = c[ok & gap[c + j - 1]]
    return c

def trova_occorrenze(abstract, forme, max_occ, seme):
    rnd = np.random.RandomState(int(seme))
    frasi = []
    for a in abstract:
        a = re.sub(r"\s+", " ", a.lower()).strip()
        frasi += [s for s in re.split(r"(?<=[.!?;]) ", a) if len(s) > 0]
    testo = "\n".join(frasi)
    inizi = np.cumsum([0] + [len(s) + 1 for s in frasi[:-1]])
    voc, ids, ini, fin, gap, pos = _indice(testo)
    occ, tot, varianti = [], [], []
    for t, f in enumerate(forme):
        fine = {}
        for x in sorted(set(f), key=len):
            n = len(x.split(" "))
            for c in _cerca(x, voc, ids, gap, pos).tolist(): fine[c] = max(fine.get(c, 0), c + n - 1)
        chiavi = sorted(fine)
        v = [(int(np.searchsorted(inizi, ini[c], "right") - 1), int(ini[c]), int(fin[fine[c]])) for c in chiavi]
        cont = Counter(testo[b:c] for _, b, c in v)
        varianti.append([[k, int(n)] for k, n in cont.most_common(40)])
        tot.append(len(v))
        idx = rnd.permutation(len(v))
        for g, j in enumerate(idx[:int(max_occ)].tolist()):
            s, b, c = v[j]
            occ.append((s, t, b - int(inizi[s]), c - int(inizi[s]), g % 2))
    per_frase = {}
    for s, t, b, c, g in occ: per_frase.setdefault(s, []).append((t, b, c, g))
    ordine = sorted(per_frase, key=lambda s: (len(frasi[s]), s))
    PIANO.update(frasi=frasi, per_frase=per_frase, ordine=ordine, n_termini=len(forme))
    return dict(n_frasi=len(frasi), n_frasi_usate=len(ordine), n_occ=len(occ), totali=tot, varianti=varianti)

def embed_contesto(modello, revisione, n_thread, max_token, file_parziale, lotto=32, ogni=200):
    tok, mod = carica(modello, revisione, n_thread)
    fr, pf, ordine, nT = PIANO["frasi"], PIANO["per_frase"], PIANO["ordine"], PIANO["n_termini"]
    d = int(mod.config.hidden_size)
    S = np.zeros((nT, 2, d)); N = np.zeros((nT, 2)); persi = 0; i0 = 0
    if os.path.exists(file_parziale):
        z = np.load(file_parziale)
        if z["S"].shape == S.shape and int(z["n_ordine"]) == len(ordine):
            S, N, persi, i0 = z["S"], z["N"], int(z["persi"]), int(z["i"])
    with torch.no_grad():
        for i in range(i0, len(ordine), lotto):
            ids = ordine[i:i + lotto]
            enc = tok([fr[s] for s in ids], padding=True, truncation=True, max_length=int(max_token),
                      return_tensors="pt", return_special_tokens_mask=True, return_offsets_mapping=True)
            spec = enc.pop("special_tokens_mask").bool()
            off = enc.pop("offset_mapping")
            h = mod(**enc).hidden_states[-2]
            ok = enc["attention_mask"].bool() & ~spec
            for r, s in enumerate(ids):
                for t, b, c, g in pf[s]:
                    m = ok[r] & (off[r, :, 0] < c) & (off[r, :, 1] > b)
                    if not bool(m.any()): persi += 1; continue
                    S[t, g] += h[r][m].mean(0).numpy().astype("float64")
                    N[t, g] += 1
            j = i + len(ids)
            if (i // lotto) % ogni == ogni - 1 or j >= len(ordine):
                prov = file_parziale.replace(".npz", "_prov.npz")
                np.savez(prov, S=S, N=N, persi=persi, i=j, n_ordine=len(ordine))
                os.replace(prov, file_parziale)
    return dict(S=S, N=N, persi=persi, commit=str(getattr(mod.config, "_commit_hash", None)))
)---")
pian <- reticulate::py$trova_occorrenze(as.list(AB), unname(lapply(FORME, as.list)), 300L, SEED)
r <- reticulate::py$embed_contesto(modello("scibert"), revisione("scibert"), 4L, 256L, "wos/context_partial.npz")
S <- r$S; N <- r$N; n_occ <- rowSums(N)
E_ctx <- (S[, 1, ] + S[, 2, ]) / pmax(n_occ, 1)
dimnames(E_ctx) <- list(TERMINI, NULL)
senza <- TERMINI[n_occ == 0]
if (length(senza)) E_ctx[senza, ] <- E_sci[senza, ]
unlink("wos/context_partial.npz")
v <- unlist(reticulate::py$VERSIONI)
writeLines(c(paste("python", reticulate::py_config()$version), paste(names(v), v),
             paste("commit", c("scibert", "sbert", "pubmed", "scibert in context"), c(attr(E_sci, "commit"), attr(E_sb, "commit"), attr(E_pm, "commit"), r$commit))),
           "results/python_versions.txt")

f_txt <- tempfile(fileext = ".txt")
writeLines(trimws(gsub("[^a-z0-9]+", " ", tolower(AB))), f_txt, useBytes = TRUE)
parole <- strsplit(gsub("[-_]+", " ", TERMINI), " ", fixed = TRUE)
allena <- function(seme) {
  ctrl <- fastTextR::ft_control(loss = "ns", word_vec_size = 100L, window_size = 5L, epoch = 5L, min_count = 5L, neg = 5L, nthreads = 1L, seed = as.integer(seme))
  m <- fastTextR::ft_train(f_txt, method = "skipgram", control = ctrl)
  W <- norm_l2(fastTextR::ft_word_vectors(m, unique(unlist(parole))))
  X <- t(vapply(parole, function(p) colMeans(W[p, , drop = FALSE]), numeric(100)))
  rownames(X) <- TERMINI; X
}
E_ft <- allena(SEED)
unlink(f_txt)

ENC <- list(scibert = E_sci, scibert_context = E_ctx, sbert = E_sb, pubmed = E_pm, fasttext = E_ft)
ENC <- lapply(ENC, function(E) { attributes(E) <- list(dim = dim(E), dimnames = list(TERMINI, NULL)); E })
attributes(E_x) <- list(dim = dim(E_x), dimnames = list(extra, NULL))
saveRDS(c(ENC, list(sbert_extra = E_x)), "results/embeddings.rds")
for (e in names(ENC)) scrivi_csv(data.frame(term = TERMINI, ENC[[e]], check.names = FALSE), sprintf("results/embeddings_%s.csv", e))
scrivi_csv(data.frame(term = extra, E_x, check.names = FALSE), "results/embeddings_sbert_extra_terms.csv")

D <- lapply(ENC, dist_coseno)
misure_aniso <- function(X) {
  Z <- norm_l2(X); C <- tcrossprod(Z); lt <- C[lower.tri(C)]; n <- nrow(Z)
  s2 <- svd(Z, nu = 0, nv = 0)$d^2
  contr <- sort((colSums(Z)^2 - colSums(Z^2)) / (n * (n - 1)), decreasing = TRUE)
  c(mean_cosine = mean(lt), first_direction = s2[1] / sum(s2), top_dimension = contr[1] / sum(contr))
}
COD <- mesh_alberi()
SIM <- lapply(COD, matrice_sim, termini = TERMINI)
tab_mesh <- do.call(rbind, lapply(names(COD), function(cd) do.call(rbind, lapply(KS, function(k) {
  a0 <- acc_S(D$scibert, SIM[[cd]], k)
  do.call(rbind, lapply(names(ENC), function(e) {
    a1 <- acc_S(D[[e]], SIM[[cd]], k)
    ic <- if (e == "scibert") c(differenza = 0, inf = NA, sup = NA) else ic_diff(a0, a1)
    data.frame(coding = cd, terms = nrow(SIM[[cd]]), k = k, encoder = e, agreement = mean(a1), difference = ic[[1]], ci_low = ic[[2]], ci_high = ic[[3]])
  }))
}))))
par <- lapply(parole, setdiff, c("of", "to", "and", "the", "in", "for", "on"))
P <- matrix(FALSE, length(TERMINI), length(TERMINI))
for (i in seq_along(TERMINI)) P[i, ] <- vapply(par, function(q) any(par[[i]] %in% q), NA)
diag(P) <- FALSE
quota_par <- function(Dx, k) { N <- vicini_knn(Dx, k); mean(vapply(seq_along(TERMINI), function(i) mean(P[i, N[i, ]]), 0)) }
top50 <- function(Dx) { Dd <- Dx; Dd[upper.tri(Dd, TRUE)] <- NA; o <- head(order(Dd, na.last = NA), 50); mean(P[o]) }
AS <- co$counts / outer(diag(co$counts), diag(co$counts)); diag(AS) <- 0
D_co <- max(AS) - AS; diag(D_co) <- 0
orig <- data.frame(encoder = names(ENC), t(vapply(ENC, misure_aniso, numeric(3))),
                   mesh_k10 = tab_mesh$agreement[tab_mesh$coding == "finale" & tab_mesh$k == 10],
                   cooccurrence_k10 = vapply(D, function(x) conservazione_knn(x, D_co, 10), 0),
                   shared_word_k10 = vapply(D, quota_par, 0, k = 10), shared_word_top50_pairs = vapply(D, top50, 0),
                   shared_word_all_pairs = mean(P[lower.tri(P)]))
coppie <- combn(names(ENC), 2, simplify = FALSE)
tab_coppie <- do.call(rbind, lapply(coppie, function(p) data.frame(a = p[1], b = p[2], k10 = conservazione_knn(D[[p[1]]], D[[p[2]]], 10))))

um <- function(X, seed, ...) {
  set.seed(seed)
  Y <- uwot::umap(norm_l2(X), n_neighbors = N_NEIGHBORS, min_dist = MIN_DIST, metric = "cosine", n_components = 2, ...)
  if (is.list(Y)) Y$embedding <- `rownames<-`(Y$embedding, rownames(X)) else rownames(Y) <- rownames(X)
  Y
}
mappe <- lapply(ENC, function(X) lapply(SEMI_UMAP, function(s) um(X, s)))
saveRDS(mappe, "results/umap_runs.rds")
stab <- lapply(mappe, stabilita_vicini, k = K_MAPPA)
fed <- lapply(names(mappe), function(e) do.call(rbind, lapply(mappe[[e]], function(Y) valuta_proiezione(D[[e]], Y, k_val = K_FED))))
names(fed) <- names(mappe)
piano <- lapply(mappe, function(l) rowMeans(sapply(l, function(Y) acc_S(dist_euclidea(Y), SIM$finale, K_MAPPA))))
tab_um <- do.call(rbind, lapply(names(mappe), function(e) {
  ist <- if (e == "scibert") c(0, NA, NA) else ic_diff(stab$scibert$per_termine, stab[[e]]$per_termine)
  ipi <- if (e == "scibert") c(0, NA, NA) else ic_diff(piano$scibert, piano[[e]])
  f <- fed[[e]]
  data.frame(encoder = e, stability = stab[[e]]$media, stability_median = stab[[e]]$mediana, stability_min = min(stab[[e]]$per_termine),
             trust = mean(f$trustworthiness), continuity = mean(f$continuity), neighbours = mean(f$conservazione_knn), spearman = mean(f$spearman),
             mesh_map = mean(piano[[e]]), d_stability = ist[[1]], d_stability_low = ist[[2]], d_stability_high = ist[[3]],
             d_mesh = ipi[[1]], d_mesh_low = ipi[[2]], d_mesh_high = ipi[[3]])
}))
NN <- lapply(mappe, function(l) lapply(l, function(Y) vicini_knn(dist_euclidea(Y), K_MAPPA)))
quota <- function(A, B) mean(vapply(seq_len(nrow(A)), function(i) length(intersect(A[i, ], B[i, ])) / K_MAPPA, 0))
R <- length(SEMI_UMAP)
tab_map <- outer(seq_along(NN), seq_along(NN), Vectorize(function(a, b) {
  if (a == b) mean(vapply(seq_len(R), function(s) quota(NN[[a]][[s]], NN[[a]][[s %% R + 1]]), 0))
  else mean(vapply(seq_len(R), function(s) quota(NN[[a]][[s]], NN[[b]][[s]]), 0))
}))
dimnames(tab_map) <- list(names(NN), names(NN))

tratta <- function(E, nome) {
  Ec <- sweep(E, 2, colMeans(E))
  if (nome == "centrato") return(Ec)
  v <- svd(Ec, nu = 0, nv = 1)$v
  Ec - Ec %*% v %*% t(v)
}
tr_m <- lapply(c(centrato = "centrato", abtt1 = "abtt1"), function(tr) { X <- tratta(E_sb, tr); lapply(SEMI_UMAP, function(s) um(X, s)) })
tr_st <- lapply(tr_m, stabilita_vicini, k = K_MAPPA)
tr_pi <- lapply(tr_m, function(l) rowMeans(sapply(l, function(Y) acc_S(dist_euclidea(Y), SIM$finale, K_MAPPA))))
tab_tr <- do.call(rbind, lapply(names(tr_m), function(tr) {
  a <- ic_diff(stab$sbert$per_termine, tr_st[[tr]]$per_termine); b <- ic_diff(piano$sbert, tr_pi[[tr]])
  data.frame(treatment = tr, stability = tr_st[[tr]]$media, mesh_map = mean(tr_pi[[tr]]), d_stability = a[[1]], d_stability_low = a[[2]], d_stability_high = a[[3]],
             d_mesh = b[[1]], d_mesh_low = b[[2]], d_mesh_high = b[[3]])
}))
Z <- norm_l2(E_sb); Dz <- dist_coseno(Z); Nx <- vicini_knn(Dz, N_NEIGHBORS - 1)
a <- um(E_sb, 1, ret_nn = TRUE)
annoy <- mean(vapply(seq_len(nrow(Z)), function(i) length(intersect(a$nn$cosine$idx[i, ], c(i, Nx[i, ]))) / N_NEIGHBORS, 0))

scrivi_csv(orig, "results/encoders_original_space.csv")
scrivi_csv(tab_coppie, "results/encoders_pairs.csv")
tab_mesh$coding <- unname(CODING_EN[tab_mesh$coding]); stopifnot(!anyNA(tab_mesh$coding))
scrivi_csv(tab_mesh, "results/encoders_mesh.csv")
scrivi_csv(tab_um, "results/encoders_maps.csv")
scrivi_csv(data.frame(encoder = rownames(tab_map), tab_map, check.names = FALSE), "results/encoders_maps_shared.csv")
tab_tr$treatment <- unname(c(centrato = "centered", abtt1 = "all-but-the-top, one direction")[tab_tr$treatment]); stopifnot(!anyNA(tab_tr$treatment))
scrivi_csv(rbind(tab_tr, data.frame(treatment = "none", stability = stab$sbert$media, mesh_map = mean(piano$sbert), d_stability = NA, d_stability_low = NA,
                                    d_stability_high = NA, d_mesh = NA, d_mesh_low = NA, d_mesh_high = NA)), "results/sbert_treatments.csv")
scrivi_csv(data.frame(annoy_equal_to_exact = annoy, mapped_terms = length(COD$finale), mapped_before_review = length(COD$autore),
                      mapped_automatic = length(COD$automatica)), "results/embeddings_checks.csv")
sessione("07")
