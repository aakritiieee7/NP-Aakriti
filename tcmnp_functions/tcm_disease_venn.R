#' Disease–Phytochemical Target Venn Analysis
#'
#' Automatically identifies overlapping proteins between
#' disease-associated genes and phytochemical targets
#'
#' @param tcm_data data.frame with column `target`
#' @param disease_name character, disease name (e.g. "Type 2 Diabetes Mellitus")
#' @param out_dir output directory
#'
#' @return character vector of overlapping genes
#' @export
#'
#' @importFrom disgenet2r disease2gene
#' @importFrom ggVennDiagram ggVennDiagram
#' @importFrom ggplot2 scale_fill_gradient theme element_text ggtitle ggsave

tcm_disease_venn <- function(
  tcm_data,
  disease_name,
  out_dir = "outputs"
) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(disgenet2r)
    library(ggVennDiagram)
    library(ggplot2)
  })

  cat("🔬 Running disease–phytochemical overlap analysis...\n")

  # -------------------------------
  # 1. Phytochemical targets
  # -------------------------------
  phyto_targets <- unique(tcm_data$target)

  cat("✔ Phytochemical targets:", length(phyto_targets), "\n")

  # -------------------------------
  # 2. Disease genes (DisGeNET)
  # -------------------------------
  disease_df <- disease2gene(
    disease = disease_name,
    database = "CURATED"
  )

  disease_genes <- unique(disease_df$gene_symbol)

  cat("✔ Disease genes retrieved:", length(disease_genes), "\n")

  # -------------------------------
  # 3. Overlap
  # -------------------------------
  common_genes <- intersect(disease_genes, phyto_targets)

  cat("✔ Common disease–phytochemical targets:", length(common_genes), "\n")

  # -------------------------------
  # 4. Save overlap genes
  # -------------------------------
  write.csv(
    data.frame(Gene = common_genes),
    file.path(out_dir, "common_disease_phyto_targets.csv"),
    row.names = FALSE
  )

  # -------------------------------
  # 5. Venn diagram
  # -------------------------------
  venn_list <- list(
    "Disease-associated proteins" = disease_genes,
    "Phytochemical targets" = phyto_targets
  )

  p <- ggVennDiagram(
    venn_list,
    label_alpha = 0
  ) +
    scale_fill_gradient(
      low = "#E8F5E9",
      high = "#1B5E20"
    ) +
    theme(
      text = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggtitle(
      paste(
        "Overlap Between",
        disease_name,
        "Proteins and Phytochemical Targets"
      )
    )

  ggsave(
    filename = file.path(out_dir, "venn_disease_phyto_targets.png"),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )

  cat("✔ Venn diagram saved\n")

  return(common_genes)
}
