source("R/functions.R")
DEF <- readRDS("wos/corpus_final.rds")
th <- leggi_csv("data/thesaurus.csv", na.strings = character(0), colClasses = c(author_keyword = "character", term = "character"))
STOP <- readLines("data/stopwords.txt", encoding = "UTF-8")
kw_doc <- lapply(strsplit(ifelse(is.na(DEF$DE), "", DEF$DE), ";"), function(v) { v <- trimws(v); v[nzchar(v)] })
kw_fin <- lapply(kw_doc, function(v) { k <- unique(th$term[match(v, th$author_keyword)]); k[!is.na(k) & nzchar(k)] })
df <- table(unlist(kw_fin))
df_ns <- sort(df[!names(df) %in% STOP], decreasing = TRUE)
T50 <- df_ns[df_ns >= 50]
DEF$DE_originale <- DEF$DE
DEF$DE <- toupper(vapply(kw_fin, paste, character(1), collapse = "; "))
DEF$DE[!nzchar(DEF$DE)] <- NA
saveRDS(DEF, "wos/corpus_clean.rds")
scrivi_csv(data.frame(term = names(df_ns), documents = as.integer(df_ns)), "results/term_frequencies.csv")
writeLines(names(T50), "results/terms_threshold_50.txt")
scrivi_csv(data.frame(measure = c("documents", "documents without author keywords", "distinct terms after cleaning",
                                  "terms with at least 50 documents", "terms with at least 30, 75, 100 documents"),
                      value = c(nrow(DEF), sum(lengths(kw_fin) == 0), length(df), length(T50),
                                paste(sapply(c(30, 75, 100), function(s) sum(df_ns >= s)), collapse = ", "))),
           "results/cleaning.csv")
sessione("02")
