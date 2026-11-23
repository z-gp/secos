source("R/functions.R")
DEF <- readRDS("wos/corpus_clean.rds")
if (!"DB" %in% names(DEF) || all(is.na(DEF$DB))) DEF$DB <- "ISI"
ANNI <- 1989:2025; INIZIO <- 1991; ANNO_RIF <- 2026; N_TAB <- 20

ut_esci <- readLines("data/ut_esci.txt")
py <- suppressWarnings(as.integer(DEF$PY))
esci <- DEF$UT_N %in% ut_esci
conta <- function(k) as.integer(table(factor(py[k], levels = ANNI)))
tassi <- function(n) {
  d <- data.frame(anno = ANNI, documenti = n)
  m <- glm(documenti ~ anno, family = quasipoisson, data = d[d$anno >= INIZIO, ])
  b <- coef(m)[["anno"]]; se <- sqrt(vcov(m)["anno", "anno"]); q <- qt(0.975, df.residual(m))
  c(compound = (n[length(n)] / n[ANNI == INIZIO])^(1 / (max(ANNI) - INIZIO)) - 1, quasi_poisson = exp(b) - 1,
    ci_low = exp(b - q * se) - 1, ci_high = exp(b + q * se) - 1, share_2021_2025 = sum(n[ANNI >= 2021]) / sum(n))
}
n_tutti <- conta(rep(TRUE, nrow(DEF))); n_esci <- conta(esci); n_senza <- n_tutti - n_esci
scrivi_csv(data.frame(year = ANNI, all = n_tutti, esci = n_esci, without_esci = n_senza), "results/documents_per_year.csv")
scrivi_csv(data.frame(series = c("all documents", "without ESCI records"), rbind(tassi(n_tutti), tassi(n_senza)),
                      esci_records = c(sum(esci), NA)), "results/growth.csv")
d_fig <- rbind(data.frame(anno = ANNI, documenti = n_tutti, serie = "All documents"),
               data.frame(anno = ANNI, documenti = n_senza, serie = "Excluding ESCI records"))
fig <- ggplot(d_fig, aes(anno, documenti, colour = serie, linetype = serie)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.4) + geom_vline(xintercept = 2014.5, colour = "grey60", linewidth = 0.3) +
  scale_colour_manual(values = c("#2a78d6", "#b05a00")) + scale_x_continuous(breaks = seq(1990, 2025, 5)) +
  labs(x = "Publication year", y = "Documents", colour = NULL, linetype = NULL) + theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave("results/fig2_growth.pdf", fig, width = 7, height = 3.8)

res <- biblioAnalysis(DEF, sep = ";")
S <- summary(res, k = 20, pause = FALSE, verbose = FALSE)
mi <- as.data.frame(S$MainInformationDF)
dt <- sort(table(trimws(DEF$DT)), decreasing = TRUE)
scrivi_csv(data.frame(item = c(trimws(mi[, 1]), paste("DT", names(dt))), value = c(trimws(mi[, 2]), as.integer(dt))), "results/table3_main_information.csv")

au <- lapply(strsplit(ifelse(is.na(DEF$AU), "", DEF$AU), ";"), function(v) {
  v <- toupper(trimws(v)); unique(v[nzchar(v) & !v %in% c("ANONYMOUS", "[ANONYMOUS]", "NA")]) })
x <- as.integer(table(unlist(au)))
THREADS <- 7L
mA <- displ$new(x); mA$setXmin(1); mA$setPars(estimate_pars(mA)); alfaA <- mA$getPars()
bsA <- bootstrap(mA, xmins = 1, no_of_sims = 1000, threads = THREADS, seed = SEED)
icA <- quantile(bsA$bootstraps$pars, c(0.025, 0.975), na.rm = TRUE)
gofA <- bootstrap_p(mA, xmins = 1, no_of_sims = 1000, threads = THREADS, seed = SEED)
mB <- displ$new(x); mB$setXmin(estimate_xmin(mB))
bsB <- bootstrap(mB, no_of_sims = 1000, threads = THREADS, seed = SEED)
icB <- quantile(bsB$bootstraps$pars, c(0.025, 0.975), na.rm = TRUE)
gofB <- bootstrap_p(mB, no_of_sims = 1000, threads = THREADS, seed = SEED)
mL <- dislnorm$new(x); mL$setXmin(mB$getXmin()); mL$setPars(estimate_pars(mL))
vu <- compare_distributions(mB, mL)
zeta <- function(a) sum((1:1e6)^(-a)) + (1e6)^(1 - a) / (a - 1)
scrivi_csv(data.frame(authors = length(x), max_documents = max(x), exponent = alfaA, ci_low = icA[[1]], ci_high = icA[[2]], p_gof = gofA$p,
                      xmin_estimated = mB$getXmin(), exponent_xmin = mB$getPars(), ci_low_xmin = icB[[1]], ci_high_xmin = icB[[2]], p_gof_xmin = gofB$p,
                      vuong_statistic = vu$test_statistic, vuong_p_two_sided = vu$p_two_sided), "results/lotka.csv")
scrivi_csv(data.frame(documents = 1:5, observed = vapply(1:5, function(k) mean(x == k), 0),
                      expected_exponent_2 = (1:5)^-2 / zeta(2), expected_fitted = (1:5)^-alfaA / zeta(alfaA)), "results/lotka_shares.csv")

if (!"AU1_CO" %in% names(DEF)) DEF <- metaTagExtraction(DEF, Field = "AU1_CO", sep = ";")
if (!"AU_CO" %in% names(DEF)) DEF <- metaTagExtraction(DEF, Field = "AU_CO", sep = ";")
dividi <- function(x) lapply(strsplit(ifelse(is.na(x), "", x), ";"), function(v) { v <- trimws(v); v[nzchar(v)] })
co_doc <- lapply(dividi(DEF$AU_CO), unique); n_co <- lengths(co_doc)
au1 <- trimws(DEF$AU1_CO); au1[!nzchar(au1)] <- NA
paesi1 <- sort(table(au1), decreasing = TRUE)
tab4 <- data.frame(country = names(paesi1), articles = as.integer(paesi1), stringsAsFactors = FALSE)
tab4$share <- tab4$articles / sum(!is.na(au1))
tab4$SCP <- vapply(tab4$country, function(p) sum(au1 == p & n_co == 1, na.rm = TRUE), integer(1))
tab4$MCP <- vapply(tab4$country, function(p) sum(au1 == p & n_co >= 2, na.rm = TRUE), integer(1))
tab4$MCP_share <- tab4$MCP / (tab4$SCP + tab4$MCP)
grandi <- tab4[tab4$articles >= 50, ]
rho <- suppressWarnings(cor.test(grandi$articles, grandi$MCP_share, method = "spearman", exact = FALSE))
scrivi_csv(tab4, "results/table4_countries.csv")
scrivi_csv(data.frame(documents_with_country = sum(!is.na(au1)), documents_without_country = sum(is.na(au1)),
                      countries_50_documents = nrow(grandi), spearman_mcp = unname(rho$estimate), p = rho$p.value,
                      documents_with_address_country = sum(n_co >= 1), international_share = sum(n_co >= 2) / sum(n_co >= 1)), "results/countries_mcp.csv")

ist <- leggi_csv("data/institutions.csv", na.strings = "")
escl <- ist$institution[ist$excluded]
mappa <- setNames(ist$merged_into, ist$institution)[!is.na(ist$merged_into)]
tenute <- unique(ifelse(is.na(ist$merged_into), ist$institution, ist$merged_into)[ist$reviewed & !ist$excluded])
applica <- function(v) { j <- match(v, names(mappa)); v[!is.na(j)] <- mappa[j[!is.na(j)]]; setdiff(unique(v), escl) }
aff_val <- lapply(lapply(dividi(DEF$C3), unique), applica)
n_ist <- sort(table(unlist(lapply(aff_val, intersect, tenute))), decreasing = TRUE)
scrivi_csv(data.frame(rank = seq_along(n_ist), institution = names(n_ist), documents = as.integer(n_ist)), "results/table_S3_institutions.csv")

TC <- suppressWarnings(as.numeric(DEF$TC)); TC[is.na(TC)] <- 0
cit <- aggregate(list(citations = TC[!is.na(au1)]), list(country = au1[!is.na(au1)]), sum)
cit$documents <- as.integer(table(au1)[cit$country]); cit$per_document <- cit$citations / cit$documents
scrivi_csv(cit[order(-cit$citations), ], "results/table_S4_citations_by_country.csv")
coppie <- unlist(lapply(co_doc[n_co >= 2], function(v) { v <- sort(v); cm <- utils::combn(v, 2); paste(cm[1, ], cm[2, ], sep = " - ") }))
coll <- sort(table(coppie), decreasing = TRUE)
scrivi_csv(data.frame(rank = seq_along(coll), pair = names(coll), documents = as.integer(coll)), "results/table_S6_collaborations.csv")

f <- list.files("wos/cited_references", pattern = "\\.txt$", full.names = TRUE)
NUOVO <- convert2df(file = f, dbsource = "wos", format = "plaintext")
NUOVO <- NUOVO[!duplicated(toupper(trimws(NUOVO$UT))), ]
CR <- NUOVO$CR[match(toupper(trimws(DEF$UT)), toupper(trimws(NUOVO$UT)))]
CR[!is.na(CR) & toupper(trimws(CR)) %in% c("NONE", "NA", "")] <- NA
A <- readRDS("wos/corpus_broad.rds")
in_def <- A$UT %in% DEF$UT
norm <- function(x) toupper(gsub("[.]", "", trimws(x)))
re_doi <- "10\\.[0-9]{4,9}/[^ ,;\\[\\]]+"
pulisci_doi <- function(x) sub("[.)]+$", "", tolower(x))
CR <- CR[!is.na(CR) & nzchar(CR)]
refs <- strsplit(CR, ";", fixed = TRUE)
citante <- rep(seq_along(refs), lengths(refs))
ref <- trimws(unlist(refs)); ok <- nzchar(ref); citante <- citante[ok]; ref <- ref[ok]
doi <- rep(NA_character_, length(ref)); h <- regexpr(re_doi, ref, perl = TRUE) > 0
doi[h] <- pulisci_doi(regmatches(ref, regexpr(re_doi, ref, perl = TRUE)))
p <- strsplit(ref, ", ", fixed = TRUE)
g_i <- function(i) vapply(p, function(v) if (length(v) >= i) v[i] else NA_character_, "")
p1 <- g_i(1); p2 <- g_i(2); p3 <- g_i(3); p4 <- g_i(4); p5 <- g_i(5)
k5 <- paste(norm(p1), trimws(p2), norm(p3), norm(p4), norm(p5), sep = "|")
k5[is.na(p5) | !grepl("^V", p4) | !grepl("^P", p5)] <- NA
k3 <- paste(norm(p1), trimws(p2), norm(p3), sep = "|")
kd <- !is.na(doi) & !is.na(k5); dk <- !duplicated(k5[kd])
mk <- k5[kd][dk]; md <- doi[kd][dk]
id <- ifelse(!is.na(doi), doi, ifelse(!is.na(k5), ifelse(k5 %in% mk, md[match(k5, mk)], k5), k3))
id[!(!is.na(doi) | grepl("^[0-9]{4}$", trimws(p2)))] <- NA
d <- !duplicated(paste(citante, id)) & !is.na(id); conta_r <- sort(table(id[d]), decreasing = TRUE)
top_id <- names(conta_r)[seq_len(min(100, length(conta_r)))]
doi_a <- pulisci_doi(tolower(trimws(A$DI))); doi_a[is.na(doi_a) | !nzchar(doi_a)] <- NA
au1_a <- vapply(strsplit(ifelse(is.na(A$AU), "", A$AU), ";"), function(v) if (length(v)) norm(v[1]) else "", "")
campo <- function(f) if (f %in% names(A)) norm(ifelse(is.na(A[[f]]), "", A[[f]])) else rep("", nrow(A))
key_a <- paste(au1_a, A$PY, campo("J9"), paste0("V", campo("VL")), paste0("P", campo("BP")), sep = "|")
ia <- match(top_id, doi_a, incomparables = NA); ib <- match(top_id, key_a); ia[is.na(ia)] <- ib[is.na(ia)]
im <- match(top_id, id)
scrivi_csv(data.frame(rank = seq_along(top_id), citing_documents = as.integer(conta_r[top_id]), first_author = p1[im], year = p2[im], source = p3[im],
                      in_corpus = !is.na(ia) & in_def[pmax(ia, 1)], reference = substr(ref[im], 1, 90)), "results/table_S5_cited_references.csv")
scrivi_csv(data.frame(documents_with_cited_references = length(CR)), "results/cited_references.csv")

nomina <- function(x) { x <- tools::toTitleCase(tolower(x)); x <- sub("Health-a Multidisciplinary", "Health: A Multidisciplinary", x, ignore.case = TRUE); gsub("\\bOf\\b", "of", gsub("\\bAnd\\b", "and", x)) }
nome_autore <- function(x) { c <- sub(" +[A-Z]+$", "", x); paste(tools::toTitleCase(tolower(c)), substring(x, nchar(c) + 2)) }
SO <- trimws(as.character(DEF$SO)); SO[!nzchar(SO)] <- NA
n_so <- sort(table(SO), decreasing = TRUE)
hidx <- function(v) { v <- sort(v, decreasing = TRUE); sum(v >= seq_along(v)) }
h_so <- sort(tapply(TC[!is.na(SO)], SO[!is.na(SO)], hidx), decreasing = TRUE)
cum <- cumsum(as.integer(n_so)); core <- sum(cum <= sum(n_so) / 3) + 1L
scrivi_csv(data.frame(journal = names(n_so), documents = as.integer(n_so), local_h_index = as.integer(h_so[names(n_so)])), "results/journals.csv")
scrivi_csv(data.frame(journals = length(n_so), core_journals = core, core_documents = cum[core], core_share = cum[core] / sum(n_so)), "results/bradford.csv")
AU <- lapply(strsplit(ifelse(is.na(DEF$AU), "", as.character(DEF$AU)), ";", fixed = TRUE), function(v) { v <- trimws(v); unique(v[nzchar(v)]) })
n_au <- sort(table(unlist(AU)), decreasing = TRUE)
scrivi_csv(data.frame(author = names(n_au)[1:N_TAB], documents = as.integer(n_au[1:N_TAB])), "results/authors.csv")
barre <- function(n, v, xl, col = COL[1]) { d <- data.frame(n = factor(n, levels = rev(n)), v = v)
  ggplot(d, aes(v, n)) + geom_col(fill = col, width = 0.68) + geom_text(aes(label = format(v, big.mark = ",")), hjust = -0.15, size = 2.3, colour = INCHIOSTRO) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) + labs(x = xl, y = NULL) + tema + theme(panel.grid.major.y = element_blank()) }
salva(barre(nomina(names(n_so)[1:15]), as.integer(n_so[1:15]), "Documents in the corpus"), "figS1_journals", 150, 85)
salva(barre(nomina(names(h_so)[1:15]), as.integer(h_so[1:15]), "Local h-index"), "figS2_h_index", 150, 85)
d <- data.frame(r = seq_along(n_so), n = as.integer(n_so))
salva(ggplot(d, aes(r, n)) + annotate("rect", xmin = 1, xmax = core, ymin = 0, ymax = Inf, fill = CHIARO, alpha = 0.6) + geom_line(colour = COL[1], linewidth = 0.55) +
  annotate("text", x = core, y = max(d$n), label = sprintf("Core: %d journals", core), hjust = 1.05, vjust = 1, size = 2.5, colour = INCHIOSTRO) +
  scale_x_log10() + labs(x = "Journal rank (log scale)", y = "Documents") + tema, "figS3_bradford", 130, 80)
salva(barre(nome_autore(names(n_au)[1:15]), as.integer(n_au[1:15]), "Documents in the corpus"), "figS4_authors", 130, 85)
o <- as.data.frame(table(x), stringsAsFactors = FALSE); o$k <- as.numeric(o$x); o$q <- o$Freq / length(x)
kk <- seq_len(max(o$k))
tt <- rbind(data.frame(k = kk, q = kk^-2 / zeta(2), serie = "Lotka, exponent 2"),
            data.frame(k = kk, q = kk^-alfaA / zeta(alfaA), serie = sprintf("Fitted power law, exponent %.2f", alfaA)))
salva(ggplot(o, aes(k, q)) + geom_point(colour = INCHIOSTRO, size = 0.9, stroke = 0, alpha = 0.8) + geom_line(data = tt, aes(colour = serie), linewidth = 0.5) +
  scale_colour_manual(values = setNames(COL[1:2], unique(tt$serie))) + scale_x_log10() +
  scale_y_log10(labels = function(v) vapply(v, function(z) if (is.na(z)) NA_character_ else format(z, scientific = FALSE, drop0trailing = TRUE), "")) +
  labs(x = "Documents per author (log scale)", y = "Share of authors (log scale)") + tema, "figS5_lotka", 130, 85)
t4 <- tab4[1:20, ]; n <- tools::toTitleCase(tolower(t4$country)); n[n == "Usa"] <- "USA"; n[n == "Korea"] <- "South Korea"
d <- rbind(data.frame(p = n, v = t4$SCP, tipo = "Single-country documents"), data.frame(p = n, v = t4$MCP, tipo = "Multiple-country documents"))
d$p <- factor(d$p, levels = rev(n)); d$tipo <- factor(d$tipo, levels = c("Single-country documents", "Multiple-country documents"))
e <- data.frame(p = factor(n, levels = rev(n)), v = t4$SCP + t4$MCP, l = sprintf("%.1f%%", 100 * t4$MCP_share))
salva(ggplot(d, aes(v, p)) + geom_col(aes(fill = tipo), width = 0.7, position = position_stack(reverse = TRUE)) + geom_text(data = e, aes(label = l), hjust = -0.15, size = 2.2, colour = INCHIOSTRO) +
  scale_fill_manual(values = c("Single-country documents" = COL[1], "Multiple-country documents" = COL[2])) + scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Documents, by country of the corresponding author", y = NULL) + tema + theme(panel.grid.major.y = element_blank()), "figS6_countries", 150, 110)
sessione("03")
