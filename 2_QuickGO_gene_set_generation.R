library(httr)
library(jsonlite)
library(dplyr)
library(readr)
library(tibble)
library(stringr)



# 1. Output directory

OUTDIR <- "./2_QuickGO_gene_set_generation_results"

dir.create(
  OUTDIR,
  recursive = TRUE,
  showWarnings = FALSE
)



# 2. Selected GO terms

go_ids <- c(
  "GO:0007155",  # cell adhesion
  "GO:0005921",  # gap junction
  "GO:0033627",  # cell adhesion mediated by integrin
  "GO:0048013",  # ephrin receptor signaling pathway
  "GO:0007219",  # Notch signaling pathway
  "GO:0001764",  # neuron migration
  "GO:0048858",  # cell projection morphogenesis
  "GO:0007411",  # axon guidance
  "GO:0097485",  # neuron projection guidance
  "GO:0071526",  # semaphorin-plexin signaling pathway
  "GO:0045216",  # cell-cell junction organization
  "GO:0034332"   # adherens junction organization
)

taxon_id <- "10090"



# 3. Utility function for QuickGO output

pick_first_existing <- function(df, candidates) {

  hit <- candidates[candidates %in% colnames(df)]

  if (length(hit) == 0) {
    return(rep(NA_character_, nrow(df)))
  }

  as.character(df[[hit[1]]])
}



# 4. Retrieve mouse annotations from QuickGO
#    including descendant GO terms

fetch_quickgo_annotations <- function(
  go_id,
  taxon_id = "10090",
  limit = 100
) {

  base_url <- "https://www.ebi.ac.uk/QuickGO/services/annotation/search"

  page <- 1
  out_list <- list()
  expected_hits <- NA_integer_

  repeat {

    cat(
      "Fetching",
      go_id,
      "page",
      page,
      "...\n"
    )

    res <- tryCatch(
      GET(
        base_url,
        query = list(
          goId = go_id,
          goUsage = "descendants",
          taxonId = taxon_id,
          limit = limit,
          page = page
        ),
        accept_json()
      ),
      error = function(e) NULL
    )

    if (is.null(res)) {
      stop(
        "QuickGO request failed: ",
        go_id,
        ", page ",
        page
      )
    }

    if (status_code(res) != 200) {
      stop(
        "QuickGO returned HTTP status ",
        status_code(res),
        ": ",
        go_id,
        ", page ",
        page
      )
    }

    txt <- content(
      res,
      as = "text",
      encoding = "UTF-8"
    )

    js <- fromJSON(
      txt,
      simplifyDataFrame = TRUE
    )

    if ("numberOfHits" %in% names(js)) {

      current_hits <- as.integer(js$numberOfHits[1])

      if (is.na(expected_hits)) {
        expected_hits <- current_hits
      } else if (current_hits != expected_hits) {
        stop(
          "QuickGO numberOfHits changed during retrieval: ",
          go_id
        )
      }
    }

    if (is.null(js$results) || nrow(js$results) == 0) {
      break
    }

    df <- as_tibble(js$results)
    df$query_GO <- go_id

    out_list[[length(out_list) + 1]] <- df

    if (!is.na(expected_hits)) {

      fetched_n <- sum(
        vapply(
          out_list,
          nrow,
          integer(1)
        )
      )

      if (fetched_n >= expected_hits) {
        break
      }

    } else if (nrow(df) < limit) {

      break
    }

    page <- page + 1
    Sys.sleep(0.2)
  }

  if (length(out_list) == 0) {

    retrieved <- tibble()

  } else {

    retrieved <- bind_rows(out_list)
  }

  retrieved_n <- nrow(retrieved)

  if (!is.na(expected_hits) &&
      retrieved_n != expected_hits) {

    stop(
      "Incomplete QuickGO retrieval for ",
      go_id,
      ": expected ",
      expected_hits,
      " annotations but retrieved ",
      retrieved_n,
      "."
    )
  }

  cat(
    go_id,
    ": expected",
    ifelse(is.na(expected_hits), "NA", expected_hits),
    ", retrieved",
    retrieved_n,
    "\n"
  )

  list(
    annotations = retrieved,
    summary = tibble(
      GO_ID = go_id,
      expected_annotations = expected_hits,
      retrieved_annotations = retrieved_n,
      retrieval_complete = ifelse(
        is.na(expected_hits),
        NA,
        retrieved_n == expected_hits
      )
    )
  )
}



# 5. Retrieve annotations for all selected GO terms

quickgo_results <- lapply(
  go_ids,
  fetch_quickgo_annotations,
  taxon_id = taxon_id,
  limit = 100
)

names(quickgo_results) <- go_ids


retrieval_summary <- bind_rows(
  lapply(
    quickgo_results,
    function(x) x$summary
  )
)

write_tsv(
  retrieval_summary,
  file.path(
    OUTDIR,
    "QuickGO_retrieval_summary.tsv"
  )
)

cat("\n===== QuickGO retrieval summary =====\n")
print(retrieval_summary, n = Inf)


ann_list <- lapply(
  quickgo_results,
  function(x) x$annotations
)

names(ann_list) <- go_ids


empty_go <- names(ann_list)[
  vapply(
    ann_list,
    nrow,
    integer(1)
  ) == 0
]

if (length(empty_go) > 0) {
  stop(
    "No annotations were retrieved for: ",
    paste(empty_go, collapse = ", ")
  )
}


ann <- bind_rows(ann_list)

write_tsv(
  ann,
  file.path(
    OUTDIR,
    "QuickGO_raw_annotations.tsv"
  )
)



# 6. Extract annotation information

gene_symbol <- pick_first_existing(
  ann,
  c("geneProductSymbol", "symbol")
)

gene_id <- pick_first_existing(
  ann,
  c("geneProductId")
)

go_id_hit <- pick_first_existing(
  ann,
  c("goId")
)

go_name_hit <- pick_first_existing(
  ann,
  c("goName")
)

aspect_hit <- pick_first_existing(
  ann,
  c("aspect")
)

evidence_code <- pick_first_existing(
  ann,
  c("goEvidence", "evidenceCode")
)

assigned_by <- pick_first_existing(
  ann,
  c("assignedBy")
)

qualifier <- pick_first_existing(
  ann,
  c("qualifier")
)

reference <- pick_first_existing(
  ann,
  c("reference")
)


annotation_table <- tibble(
  query_GO = ann$query_GO,
  gene_symbol = gene_symbol,
  gene_product_id = gene_id,
  matched_GO = go_id_hit,
  matched_GO_name = go_name_hit,
  aspect = aspect_hit,
  evidence_code = evidence_code,
  assigned_by = assigned_by,
  qualifier = qualifier,
  reference = reference
) %>%
  filter(
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  distinct()



# 7. Exclude NOT annotations

is_not_annotation <- (
  !is.na(annotation_table$qualifier) &
  str_detect(
    annotation_table$qualifier,
    "(^|\\|)NOT($|\\|)"
  )
)


excluded_not_annotations <- annotation_table[
  is_not_annotation,
]

write_tsv(
  excluded_not_annotations,
  file.path(
    OUTDIR,
    "QuickGO_excluded_NOT_annotations.tsv"
  )
)


genes_by_each_go <- annotation_table[
  !is_not_annotation,
]

write_tsv(
  genes_by_each_go,
  file.path(
    OUTDIR,
    "genes_by_each_GO_mouse.tsv"
  )
)


cat(
  "\nNOT annotation rows excluded:",
  nrow(excluded_not_annotations),
  "\n"
)



# 8. Merge gene lists and remove duplicated genes

merged_union_gene_list <- genes_by_each_go %>%
  select(gene_symbol) %>%
  distinct() %>%
  arrange(gene_symbol)


write_tsv(
  merged_union_gene_list,
  file.path(
    OUTDIR,
    "merged_union_gene_list_mouse.tsv"
  )
)



# 9. Summary

cat(
  "\nSelected GO terms:",
  length(go_ids),
  "\n"
)

cat(
  "GO terms successfully retrieved:",
  nrow(retrieval_summary),
  "\n"
)

cat(
  "Total QuickGO annotation rows:",
  nrow(ann),
  "\n"
)

cat(
  "NOT annotation rows excluded:",
  nrow(excluded_not_annotations),
  "\n"
)

cat(
  "Unique genes in merged GO gene set:",
  nrow(merged_union_gene_list),
  "\n"
)

cat(
  "Output directory:",
  normalizePath(OUTDIR),
  "\n"
)

cat("\nDone.\n")