library(ggplot2)

base_dir <- "results/enrichment/great"
outdir <- "results/enrichment/plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

datasets <- list.dirs(base_dir, recursive = FALSE)

safe_log <- function(p) {
  p[p == 0] <- 1e-20
  -log10(p)
}

results <- list()

for (d in datasets) {
  dataset_name <- basename(d)
  file <- file.path(d, "gobp.csv")

  if (!file.exists(file)) next

  df <- read.csv(file)

  if (!all(c("name", "Binom_Adjp_BH", "Binom_Fold_Enrichment") %in% colnames(df))) {
    next
  }

  df <- df[order(df$Binom_Adjp_BH), ]
  df <- head(df, 15)

  if (nrow(df) == 0) next

  df$name <- as.character(df$name)
  df$term_factor <- factor(df$name, levels = rev(df$name))
  df$neg_log_p <- safe_log(df$Binom_Adjp_BH)

  p <- ggplot(df, aes(Binom_Fold_Enrichment, term_factor, fill = neg_log_p)) +
    geom_col() +
    scale_fill_gradient(low = "peachpuff", high = "darkred", name = "-log10(adj p)") +
    labs(title = dataset_name, x = "Fold Enrichment", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 8)
    )

  ggsave(
    file.path(outdir, paste0("barplot_", dataset_name, ".png")),
    p,
    width = 9,
    height = 6,
    dpi = 150
  )

  results[[dataset_name]] <- df
}

if (length(results) == 0) {
  stop("No valid GREAT result files found.")
}

top_terms <- unique(unlist(lapply(results, function(df) head(df$name, 10))))

hm_list <- list()

for (dataset_name in names(results)) {
  df <- results[[dataset_name]]

  tmp <- data.frame(
    term = top_terms,
    dataset = dataset_name,
    value = sapply(top_terms, function(t) {
      v <- df$Binom_Fold_Enrichment[df$name == t]
      if (length(v) == 0) NA else v[1]
    }),
    stringsAsFactors = FALSE
  )

  hm_list[[dataset_name]] <- tmp
}

hm <- do.call(rbind, hm_list)

hm$dataset <- factor(hm$dataset, levels = names(results))
hm$term <- factor(hm$term, levels = rev(top_terms))

p <- ggplot(hm, aes(dataset, term, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",
    high = "darkred",
    na.value = "grey90",
    name = "Fold\nEnrichment"
  ) +
  labs(title = "GO:BP comparison", x = NULL, y = NULL) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 10),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  )

ggsave(
  file.path(outdir, "comparison_heatmap.png"),
  p,
  width = 12,
  height = 14,
  dpi = 150
)
