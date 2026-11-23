source("R/functions.R")
norm_ut <- function(x) toupper(trimws(sub("^WOS:", "", x, ignore.case = TRUE)))
M <- readRDS("wos/march_2026.rds")
M$UT_N <- norm_ut(M$UT)
ut_q <- readLines("data/ut_proximity_query.txt")
tipo_ok <- function(dt) vapply(strsplit(toupper(dt), ";"), function(v) any(trimws(v) %in% c("ARTICLE", "REVIEW")), logical(1))
anno_ok <- function(py) { py <- suppressWarnings(as.numeric(py)); !is.na(py) & py >= 1989 & py <= 2025 }
R0 <- M[M$UT_N %in% ut_q, ]
R1 <- R0[anno_ok(R0$PY), ]
R2 <- R1[tipo_ok(R1$DT), ]
A1 <- M[anno_ok(M$PY), ]
A2 <- A1[tipo_ok(A1$DT), ]

tokens <- function(s) {
  if (is.na(s) || !nzchar(s)) return(character(0))
  tk <- strsplit(gsub("[^a-z]+", " ", tolower(s)), " ")[[1]]; tk[nzchar(tk)]
}
near_k <- function(testo, k) {
  vapply(testo, function(s) {
    tok <- tokens(s); if (length(tok) == 0) return(FALSE)
    i_inj <- which(tok %in% c("injury", "injuries"))
    nxt <- c(tok[-1], "")
    i_spo <- which(startsWith(tok, "sport") | startsWith(tok, "exercise") | (tok == "physical" & startsWith(nxt, "activit")))
    length(i_inj) > 0 && length(i_spo) > 0 && any(abs(outer(i_inj, i_spo, "-")) <= k)
  }, logical(1), USE.NAMES = FALSE)
}
n_inj <- function(testo) vapply(testo, function(s) sum(startsWith(tokens(s), "injur")), integer(1), USE.NAMES = FALSE)
ESCLUSIONE <- c(
  "ischemi", "ischaemi", "reperfusion", "\\bmyocard", "cardiac injur", "\\bkidneys?\\b", "\\brenal\\b",
  "\\bnephro", "\\bliver\\b", "\\bhepat", "\\bmice\\b", "\\bmouse\\b", "\\bmurine\\b", "\\brats?\\b",
  "\\brodents?\\b", "in vitro", "spinal cord injur", "\\bstrokes?\\b",
  "older adults?.{0,40}fall", "elderly.{0,40}fall", "fall prevention",
  "\\bhorses?\\b", "\\bequine\\b", "\\bthoroughbred", "\\bcanine\\b", "\\bdogs?\\b")
p_B3 <- paste(ESCLUSIONE, collapse = "|")
cat_rev <- c("SPORT SCIENCES", "ORTHOPEDICS", "REHABILITATION", "PHYSIOLOGY", "SURGERY")
ha_cat <- function(wc) vapply(strsplit(toupper(ifelse(is.na(wc), "", wc)), ";"), function(v) any(trimws(v) %in% cat_rev), logical(1))
classifica <- function(TI, AB, DE, WC) {
  forte <- near_k(TI, 10) | near_k(DE, 10)
  ni <- n_inj(paste(TI, AB))
  data.frame(A = forte, W = ha_cat(WC), B3 = !grepl(p_B3, tolower(paste(TI, AB, DE)), perl = TRUE),
             S2 = forte | ni >= 2, S3 = forte | ni >= 3, N5 = near_k(TI, 5) | near_k(AB, 5) | near_k(DE, 5))
}
C <- classifica(R2$TI, R2$AB, R2$DE, R2$WC)
DEF <- R2[C$B3, ]

f <- leggi_csv("data/coded_sample.csv", colClasses = c(UT = "character"))
i <- match(f$UT, R2$UT)
S <- classifica(R2$TI[i], R2$AB[i], R2$DE[i], R2$WC[i])
regole <- list(
  "None (NEAR/10 query only)"                  = function(d) rep(TRUE, nrow(d)),
  "A: proximity in title or author keywords"   = function(d) d$A,
  "W: at least one of the five WoS categories" = function(d) d$W,
  "B3: exclusion terms"                        = function(d) d$B3,
  "A + B3"                                     = function(d) d$A & d$B3,
  "W + B3"                                     = function(d) d$W & d$B3,
  "B3 + injury salience (>=2)"                 = function(d) d$B3 & d$S2,
  "B3 + injury salience (>=3)"                 = function(d) d$B3 & d$S3,
  "B3 + NEAR/5"                                = function(d) d$B3 & d$N5)
tot <- sum(f$consensus)
tab <- do.call(rbind, lapply(names(regole), function(nm) {
  k <- regole[[nm]](S); x <- sum(f$consensus[k]); w <- wilson(x, sum(k))
  data.frame(filter = nm, retained_sample = sum(k), relevant_retained = x,
             precision = round(w[1], 3), ci_low = round(w[2], 3), ci_high = round(w[3], 3),
             relevant_lost = tot - x, relevant_lost_pct = round(100 * (tot - x) / tot),
             corpus_size = sum(regole[[nm]](C)), stringsAsFactors = FALSE)
}))

flusso <- data.frame(step = c("broad query", "broad, years 1989-2025", "broad corpus (Article and Review)",
                              "proximity query", "proximity, matched to the March export", "proximity, years 1989-2025",
                              "proximity corpus (Article and Review)", "final corpus (exclusion terms)"),
                     records = c(nrow(M), nrow(A1), nrow(A2), length(ut_q), nrow(R0), nrow(R1), nrow(R2), nrow(DEF)))
scrivi_csv(flusso, "results/flow.csv")
scrivi_csv(tab, "results/table_S1.csv")
writeLines(ESCLUSIONE, "results/exclusion_terms.txt")
writeLines(A2$UT, "results/ut_broad_corpus.txt")
writeLines(R2$UT, "results/ut_proximity_corpus.txt")
writeLines(DEF$UT, "results/ut_final_corpus.txt")
saveRDS(A2, "wos/corpus_broad.rds")
saveRDS(R2, "wos/corpus_proximity.rds")
saveRDS(DEF, "wos/corpus_final.rds")
sessione("01")
