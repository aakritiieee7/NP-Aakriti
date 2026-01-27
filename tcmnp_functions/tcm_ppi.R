tcm_ppi <- function(
  targets,
  species = 9606,
  score_threshold = 700,
  degree_filter = 2
) {

  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    install.packages("STRINGdb", repos = "https://cloud.r-project.org")
  }
  library(STRINGdb)
  library(dplyr)

  message("✔ Initializing STRINGdb")

  string_db <- STRINGdb$new(
    version = "11.5",
    species = species,
    score_threshold = score_threshold,
    input_directory = ""
  )

  # Map genes
  mapped <- string_db$map(
    data.frame(gene = targets),
    "gene",
    removeUnmappedRows = TRUE
  )

  message("✔ Targets mapped to STRING:", nrow(mapped))

  # Fetch interactions
  ppi_raw <- string_db$get_interactions(mapped$STRING_id)
  message("✔ STRING interactions fetched:", nrow(ppi_raw))

  # ID → gene mapping
  id_map <- mapped %>%
    dplyr::select(STRING_id, gene)

  ppi <- ppi_raw %>%
    dplyr::left_join(
  id_map,
  by = c("from" = "STRING_id"),
  relationship = "many-to-many") %>%
    dplyr::rename(from_gene = gene) %>%
    dplyr::left_join(id_map, by = c("to" = "STRING_id")) %>%
    dplyr::rename(to_gene = gene) %>%
    dplyr::select(from_gene, to_gene, combined_score) %>%
    stats::na.omit()

  colnames(ppi) <- c("from", "to", "weight")

  # -----------------------------
  # DEGREE FILTER (anti-clutter)
  # -----------------------------
  deg <- ppi %>%
    dplyr::count(from, name = "deg_from") %>%
    dplyr::full_join(
      dplyr::count(ppi, to, name = "deg_to"),
      by = c("from" = "to")
    ) %>%
    dplyr::mutate(
      deg_from = ifelse(is.na(deg_from), 0, deg_from),
      deg_to   = ifelse(is.na(deg_to), 0, deg_to),
      degree   = deg_from + deg_to
    )

  hub_nodes <- deg %>%
    dplyr::filter(degree >= degree_filter) %>%
    dplyr::pull(from)

  ppi_filt <- ppi %>%
    dplyr::filter(from %in% hub_nodes & to %in% hub_nodes)

  message("✔ PPI after degree filter (≥ ", degree_filter, "): ", nrow(ppi_filt))

  return(ppi_filt)
}
