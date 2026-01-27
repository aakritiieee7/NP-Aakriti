# ====== GRAPHICS ======
options(bitmapType = "cairo")

pdf.options(
  family = "Helvetica",
  useDingbats = FALSE
)

Sys.setenv(R_GSCMD = "/usr/bin/gs")

cat("Starting Network Pharmacology Analysis...\n")

# ===============================
# 1. Load required libraries
# ===============================
required_pkgs <- c(
  "dplyr", "ggplot2", "igraph", "ggraph",
  "stringr", "magrittr", "RColorBrewer",
  "clusterProfiler", "org.Hs.eg.db", "DOSE",
  "ggsankey"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

cat("✔ Core R packages loaded\n")

theme_set(
  ggplot2::theme(
    text = ggplot2::element_text(family = "sans")
  )
)

# ===============================
# 2. Load TCMNP functions
# ===============================
TCMNP_FUN_PATH <- "tcmnp_functions"

safe_files <- c(
  "tcm_net.R",
  "tcm_sankey.R",
  "degree_plot.R",
  "tcm_alluvial.R",
  "bar_plot.R",
  "tcm_ppi.R",
  "ppi_plot.R",
  "hub_gene_identification.R",
  "tcm_disease_venn.R"  
)

for (f in safe_files) {
  source(file.path(TCMNP_FUN_PATH, f))
}

cat("✔ Required TCMNP functions loaded\n")

# ===============================
# 3. Load input data
# ===============================
INPUT_FILE <- "3_tcmnp_input/tcm_input.csv"
OUTPUT_DIR <- "outputs"

dir.create(OUTPUT_DIR, showWarnings = FALSE)

tcm_data <- read.csv(INPUT_FILE, stringsAsFactors = FALSE)
colnames(tcm_data) <- tolower(colnames(tcm_data))
tcm_data <- dplyr::distinct(tcm_data)

cat("✔ Input data loaded:", nrow(tcm_data), "rows\n")

# ===============================
# 4. Herb–Compound–Target Network
# ===============================
png(
  file.path(OUTPUT_DIR, "tcm_network.png"),
  width = 4000,
  height = 3000,
  res = 300,
  type = "cairo"
)

tcm_net(
  tcm_data,
  label.degree = 0,
  rem.dis.inter = FALSE
)

dev.off()
cat("✔ TCM network plot generated\n")

# ===============================
# 5. STRING PPI NETWORK
# ===============================
cat("Starting STRING PPI analysis...\n")

ppi_data <- tcm_ppi(
  targets = unique(tcm_data$target),
  degree_filter = 2
)

write.csv(
  ppi_data,
  file.path(OUTPUT_DIR, "ppi_edges_filtered.csv"),
  row.names = FALSE
)

png(
  file.path(OUTPUT_DIR, "ppi_network.png"),
  width = 4000,
  height = 3000,
  res = 300,
  type = "cairo"
)

ppi_plot(
  data = ppi_data,
  label.degree = 5,
  node.size = c(2, 12),
  edge.color = "grey70",
  graph.layout = "kk"
)

dev.off()
cat("✔ STRING PPI network generated\n")

# ================================
# 6. Automated Hub Gene Identification
# ================================

hub_genes_auto <- hub_gene_identification(
  ppi_data,
  top_n = Inf
)

write.csv(
  hub_genes_auto,
  file.path(OUTPUT_DIR, "hub_genes_automated_all.csv"),
  row.names = FALSE
)

cat("✔ Automated hub gene identification completed\n")
cat("✔ Total hub genes:", nrow(hub_genes_auto), "\n")

# ===============================
# 6. Sankey Plot (PDF + PNG)
# ===============================
set.seed(123)
sankey_data <- tcm_data[sample(nrow(tcm_data), min(30, nrow(tcm_data))), ]

pdf(file.path(OUTPUT_DIR, "sankey_plot.pdf"), width = 10, height = 8)
tcm_sankey(sankey_data, text.size = 3, x.axis.text.size = 18)
dev.off()

png(
  file.path(OUTPUT_DIR, "sankey_plot.png"),
  width = 3000,
  height = 2400,
  res = 300,
  type = "cairo"
)
tcm_sankey(sankey_data, text.size = 3, x.axis.text.size = 18)
dev.off()

cat("✔ Sankey plot generated\n")

# ===============================
# 7. Hub Gene Analysis
# ===============================

hub_targets_all <- tcm_data %>%
  dplyr::count(target, sort = TRUE)

write.csv(
  hub_targets_all,
  file.path(OUTPUT_DIR, "hub_targets_all.csv"),
  row.names = FALSE
)

cat("✔ All hub targets saved\n")

hub_genes_all <- tcm_data %>%
  count(target, sort = TRUE)

write.csv(
  hub_genes_all,
  file.path(OUTPUT_DIR, "hub_genes_all.csv"),
  row.names = FALSE
)

pdf(file.path(OUTPUT_DIR, "degree_plot.pdf"), width = 8, height = 6)
degree_plot(tcm_data, plot.set = "horizontal")
dev.off()

png(
  file.path(OUTPUT_DIR, "degree_plot.png"),
  width = 2400,
  height = 1800,
  res = 300,
  type = "cairo"
)
degree_plot(tcm_data, plot.set = "horizontal")
dev.off()

cat("✔ Hub gene analysis completed\n")

# ===============================
# 8. Functional Enrichment
# ===============================
eg <- bitr(
  unique(tcm_data$target),
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)
eg <- eg[!is.na(eg$ENTREZID), ]

kk <- enrichKEGG(
  gene = eg$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

pdf(file.path(OUTPUT_DIR, "kegg_enrichment.pdf"), width = 9, height = 7)
bar_plot(kk, title = "KEGG Pathway Enrichment")
dev.off()

png(
  file.path(OUTPUT_DIR, "kegg_enrichment.png"),
  width = 2700,
  height = 2100,
  res = 300,
  type = "cairo"
)
bar_plot(kk, title = "KEGG Pathway Enrichment")
dev.off()
# ===============================
# 9. Pathview Visualization (RED = hubs, GREEN = others)
# ===============================

library(pathview)

cat("Starting KEGG pathway visualization with hub highlighting...\n")

# ---- Convert enrichment to data.frame ----
kegg_df <- as.data.frame(kk)

# ---- Output directory ----
PV_DIR <- file.path(OUTPUT_DIR, "kegg_pathview")
dir.create(PV_DIR, showWarnings = FALSE)

# ---- Identify PI3K–AKT pathway ----
pi3k_row <- kegg_df[grepl("PI3K-Akt", kegg_df$Description, ignore.case = TRUE), ]

if (nrow(pi3k_row) == 0) {
  stop("❌ PI3K–Akt pathway not found in KEGG enrichment")
}

pi3k_genes <- unlist(strsplit(pi3k_row$geneID[1], "/"))

# ---- Hub genes (already computed earlier) ----
hub_symbols <- hub_targets_all$target

hub_eg <- bitr(
  hub_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

hub_eg_ids <- intersect(hub_eg$ENTREZID, pi3k_genes)

# ---- Build gene vector ----
gene_fc <- rep(-1, length(pi3k_genes))   # non-hubs = green
names(gene_fc) <- pi3k_genes
gene_fc[hub_eg_ids] <- 1                 # hubs = red

# ---- Run pathview ----
old_wd <- getwd()
setwd(PV_DIR)

pathview(
  gene.data   = gene_fc,
  pathway.id  = "hsa04151",
  species     = "hsa",
  kegg.native = TRUE,
  same.layer  = FALSE,
  low  = list(gene = "#2ECC71"),  # solid green
  mid  = list(gene = "#FFFFFF"),  # white
  high = list(gene = "#E74C3C"),  # solid red
  limit = list(gene = c(-1, 1)),
  bins = 2,
  out.suffix = "PI3K_AKT_HUBS"
)

setwd(old_wd)

cat("✔ Pathview generated (red = hubs, green = others)\n")

bp <- enrichGO(
  gene = eg$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pvalueCutoff = 0.05,
  readable = TRUE
)

pdf(file.path(OUTPUT_DIR, "go_bp_enrichment.pdf"), width = 9, height = 7)
bar_plot(bp, title = "GO Biological Process Enrichment")
dev.off()

png(
  file.path(OUTPUT_DIR, "go_bp_enrichment.png"),
  width = 2700,
  height = 2100,
  res = 300,
  type = "cairo"
)
bar_plot(bp, title = "GO Biological Process Enrichment")
dev.off()

cat("✔ Functional enrichment completed\n")

# ===============================
# Disease–Phytochemical Venn
# ===============================
common_genes <- tcm_disease_venn(
  tcm_data,
  disease_name = "Type 2 Diabetes Mellitus",
  out_dir = OUTPUT_DIR
)

cat("✔ Overlap genes ready for downstream analysis\n")

# ===============================
# DONE
# ===============================
cat("\nANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Outputs saved in:", OUTPUT_DIR, "\n")
