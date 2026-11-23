# SeCoS: Semantic Conceptual Structure — code and data

Code and derived data for the article *Mapping sports injury research using bibliometric and embedding-based analysis* (G. Zammarchi and C. Pisanu, *Computational Statistics*).

The repository contains:

- `R/secos.R`: the function `secos()`, which maps the author keywords of a bibliographic collection from the meaning of the terms. It works on the data frame produced by `bibliometrix::convert2df()`, and the language model and the weight λ of co-occurrence are arguments. The file documents each function and gives examples.
- `scripts/`: ten scripts that reproduce every number, table and figure of the article and of its Supplementary Material.
- `data/`: the inputs that can be distributed.
- `results/`: the files produced by the scripts, as obtained by the authors.

## Data and licensing of the records

The bibliographic records come from the Web of Science Core Collection and are subject to the terms of use of Clarivate, so they are not distributed. The repository gives their identifiers (UT), so that users with access to Web of Science can download them again:

| File | Content |
|---|---|
| `data/ut_march_2026.txt` | UT of the broad search exported on 23 March 2026 |
| `data/ut_proximity_query.txt` | UT of the proximity query run on 30 October 2026 |
| `data/ut_esci.txt` | UT of the records of the Emerging Sources Citation Index in the broad search |
| `results/ut_broad_corpus.txt`, `results/ut_proximity_corpus.txt`, `results/ut_final_corpus.txt` | UT of the three corpora of Table S2; the final corpus has 21,394 records |

The scripts read the records from the folder `wos/`, which is not distributed:

- `wos/march_2026.rds`: the March export as a `bibliometrix` data frame, in the order of `data/ut_march_2026.txt` (records exported in plain text and converted with `bibliometrix::convert2df(..., dbsource = "wos", format = "plaintext")`);
- `wos/cited_references/`: plain-text exports, with cited references, of the records of the proximity corpus (exported on 2 November 2026), used only for the cited references.

Results computed from the records depend on the state of Web of Science at the time of export (citations, indexing). Steps that resample documents (script 09) also depend on the order of the records.

## Requirements

R 4.6.0 with: bibliometrix 5.4.0, uwot 0.2.4, Rtsne 0.17, ca 0.71.1, ggplot2 4.0.3, poweRlaw 1.0.0, igraph 2.3.1, Matrix 1.7-5, dbscan 1.2.7, cluster 2.1.8.2, patchwork 1.3.2, ggrepel 0.9.8, fastTextR 2.1.1, reticulate 1.46.0. The authors used the reference BLAS and LAPACK libraries of R.

Scripts 07 and 10 compute embeddings in Python through reticulate (by default from the conda environment `r-text`; set `RETICULATE_PYTHON` to use another one), with transformers 4.30.2, tokenizers 0.13.3 and torch 2.0.1 (CPU). The models are fixed to the revisions in `data/models.csv`:

| Encoder | Model | Revision |
|---|---|---|
| SBERT (adopted) | `sentence-transformers/all-mpnet-base-v2` | `e8c3b32` |
| SciBERT | `allenai/scibert_scivocab_uncased` | `24f92d3` |
| PubMedBERT | `microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract` | `d673b88` |

`results/python_versions.txt` records the Python environment of the authors' run, and `results/sessionInfo_XX.txt` the R session of each script.

## How to run

From the root of the repository, each script in a new R session and in this order:

```
Rscript scripts/01_corpus.R
Rscript scripts/02_cleaning.R
...
Rscript scripts/10_figures_timing.R
```

| Script | Content | Indicative time |
|---|---|---|
| `01_corpus.R` | corpus construction, precision of the relevance check, exclusion filter | minutes |
| `02_cleaning.R` | keyword thesaurus, frequencies, the 171 terms | under a minute |
| `03_performance.R` | growth, main information, Lotka's law, countries, institutions, citations, collaborations, cited references, journals and authors | 20–30 min |
| `04_cooccurrence.R` | co-occurrence network, thematic map, stability across seeds | 20–30 min |
| `05_query_comparison.R` | broad, proximity and final corpus (Table S2), broader rule | 20–40 min |
| `06_simulation.R` | simulation study | 15–20 min |
| `07_embeddings.R` | MCA, embeddings of the five encoders, encoder comparison | 1–3 h |
| `08_benchmark.R` | SeCoS against MCA | about 10 min |
| `09_stability_hybrid.R` | sensitivity, consensus map, hybrid version, resampling | 30–60 min |
| `10_figures_timing.R` | maps, stability figures, uncertain positions, interactive map, computation times | about 40 min |

Scripts 01 to 05 and part of 07, 09 and 10 need the records in `wos/`. Scripts 08 and the parts of 09 and 10 on the semantic maps start from the files in `results/` (embeddings, co-occurrence matrix, MCA coordinates), so they can be run without the records.

### Random numbers

Every random step sets its seed (base seed 2026) before it runs. All scripts use R's default generator (Mersenne-Twister, Inversion, Rejection), except `04_cooccurrence.R` and `05_query_comparison.R`, which use L'Ecuyer-CMRG, the generator active when those analyses were run for the article. The bootstrap of Lotka's law (script 03) runs on 7 parallel workers, each with its own L'Ecuyer-CMRG stream from the seed: with another number of workers the bootstrap interval changes slightly. The second, independent set of 100 consensus replicates used for the uncertain positions (script 09) is generated with the same seeds and the L'Ecuyer-CMRG generator. Computation times (Table 7) depend on the computer and are not expected to match.

## Where each result comes from

| Article | Script | Files in `results/` |
|---|---|---|
| Fig. 1, Sect. 2.1 | 01 | `flow.csv` |
| Sect. 2.1 (keyword cleaning) | 02 | `cleaning.csv`, `term_frequencies.csv`, `terms_threshold_50.txt` |
| Table 3 | 03 | `table3_main_information.csv` (summary of `bibliometrix`), `growth.csv`; authors from `lotka.csv` (distinct author names, without anonymous entries); international co-authorships from `countries_mcp.csv` (documents with addresses in at least two countries, among those with a recognizable country) |
| Table 4 | 03 | `table4_countries.csv`, `countries_mcp.csv` |
| Table 5 | 08, 09 | `benchmark_fidelity.csv`, `benchmark_mesh.csv`, `benchmark_partitions.csv`, `mesh_reference.csv`, `consensus_summary.csv` |
| Table 6 | 10 | `table6_stability.csv` |
| Table 7 | 10 | `table7_timing.csv`, `timing_slopes.csv`, `timing_machine.txt` |
| Table 8 | 07 | `encoders_original_space.csv`, `encoders_maps.csv`, `encoders_mesh.csv` |
| Appendix table (stopwords) | — | `data/stopwords.txt` |
| Fig. 2 | 03 | `fig2_growth.pdf` |
| Figs. 3 and 4 | 04 | `fig3_network.pdf`, `fig4_thematic_map.pdf` |
| Fig. 5 | 06 | `fig5_simulation_ari.pdf`, `simulation_summary.csv` |
| Figs. 6 to 9 | 10 | `fig6_maps.pdf`, `fig7_stability.pdf`, `fig8_sensitivity.pdf`, `fig9_lambda.pdf` |
| Sect. 3.1 | 01, 03, 05 | `table_S1.csv`, `growth.csv`, `documents_per_year.csv`, `table_S2.csv` |
| Sect. 3.2 | 03 | `journals.csv`, `bradford.csv`, `authors.csv`, `lotka.csv`, `lotka_shares.csv`, `cited_references.csv`, tables S3 to S6 |
| Sect. 3.3 | 04 | `network_summary.csv`, `seed_stability.csv`, `theme_stability.csv`, `thematic_map_themes.csv`, `thematic_map_terms.csv` |
| Sect. 3.4 | 06 | `simulation_summary.csv`, `simulation_replicates.csv` |
| Sect. 3.5 | 07, 08 | `mca_summary.csv`, `benchmark_*.csv` |
| Sect. 3.6 | 09, 10 | `table_S9_sensitivity.csv`, `stability_summary.csv`, `consensus_*.csv`, `stability_terms.csv`, `resampling.csv`, `uncertain_summary.csv`, `uncertain_terms.csv` |
| Sect. 3.7 | 07, 09 | `encoders_*.csv`, `sbert_treatments.csv`, `hybrid_space.csv`, `hybrid_space_controls.csv`, `hybrid_maps.csv`, `hybrid_periods.csv`, `hybrid_periods_documents.csv` |
| Table S1 | 01 | `table_S1.csv` |
| Table S2 | 05 | `table_S2.csv`, `table_S2_notes.csv` |
| Tables S3 to S6 | 03 | `table_S3_institutions.csv`, `table_S4_citations_by_country.csv`, `table_S5_cited_references.csv`, `table_S6_collaborations.csv` |
| Tables S7 and S8 | 04 | `table_S7_network.csv`, `table_S8_thematic_map.csv` |
| Table S9 | 09 | `table_S9_sensitivity.csv` |
| Figs. S1 to S6 | 03 | `figS1_journals.pdf` … `figS6_countries.pdf` |
| Fig. S7 | 06 | `figS7_simulation_purity.pdf` |
| Figs. S8 and S9 | 10 | `figS8_zoom.pdf`, `figS9_maps_louvain.pdf` |
| Interactive map | 10 | `interactive_map.html` |

## Derived objects

| File | Content |
|---|---|
| `results/terms.csv` | the 171 terms of SeCoS and MCA, with their number of documents |
| `results/embeddings_*.csv`, `results/embeddings.rds` | embeddings of the 171 terms for the five encoders, labeled `scibert` (SciBERT on isolated terms), `scibert_context` (SciBERT in context), `sbert` (the adopted encoder, all-mpnet-base-v2), `pubmed` (PubMedBERT) and `fasttext`, and of the further terms used in the sensitivity analysis (`sbert_extra_terms`); the same labels are used in the `encoders_*.csv` files |
| `results/cooccurrence_matrix.csv` | co-occurrence counts of the 171 terms (document frequencies on the diagonal), over the 13,564 documents with at least one of them |
| `results/mca_coordinates.csv` | MCA coordinates of the terms (`bibliometrix::conceptualStructure`) |
| `results/umap_runs.rds` | the 20 single UMAP runs of each encoder (seeds 1 to 20; k = 10, min_dist = 0.1, cosine) |
| `results/louvain_runs.rds` | the 200 Louvain runs on the co-occurrence network of the 171 terms (script 08) |
| `results/secos_check.csv` | check that `secos_map()` reproduces the consensus map (maximum differences and uncertain terms) |
| `results/consensus_replicates.rds`, `results/consensus_map.csv` | the 100 replicates of the consensus map (seeds 1 to 100, exact neighbors, fixed spectral initialization), the second independent set, and the consensus map with its ellipses and stability |
| `data/mesh_mapping.csv` | mapping of the terms onto MeSH 2026, with the tree numbers of the final coding and of the two control codings (labeled `final`, `before_review` and `automatic` in `benchmark_mesh.csv` and `encoders_mesh.csv`) |
| `data/thesaurus.csv`, `data/acronyms.csv`, `data/synonyms.csv` | keyword thesaurus: each author keyword with its final term and the rule applied; validated acronyms; manual synonyms |
| `data/institutions.csv` | review of the 100 most frequent institution entries: exclusions and merges |
| `data/coded_sample.csv` | the 200 sampled records of the proximity corpus (UT) with the consensus code (1 = relevant), from which script 01 computes the precisions of Table S1; the agreement between the two coders is summarized in Section S1 of the Supplementary Material |
| `results/exclusion_terms.txt` | exclusion terms of the relevance filter (regular expressions) |

## Using SeCoS on other data

```r
source("R/secos.R")
library(bibliometrix)
M <- convert2df("savedrecs.txt", dbsource = "wos", format = "plaintext")
reticulate::use_condaenv("my-env")          # Python with torch and transformers
x <- secos(M, field = "DE", min_docs = 50, remove_terms = c("sport", "injury"))
plot(x, labels = 25)
head(x$map)                                  # coordinates, ellipses, stability, uncertain positions
```

From a matrix of embeddings computed elsewhere, `secos_map()` gives the consensus map directly, for instance `secos_map(as.matrix(read.csv("results/embeddings_sbert.csv", row.names = 1)))`, which reproduces `results/consensus_map.csv`.

## Citation

Please cite the article and the archived version of this repository: doi:10.5281/zenodo.XXXXXXX.

## License

The code (`R/` and `scripts/`) is released under the MIT License (`LICENSE`). The data files in `data/` and `results/` are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0, https://creativecommons.org/licenses/by/4.0/). Neither license applies to the Web of Science records, which are not distributed.
