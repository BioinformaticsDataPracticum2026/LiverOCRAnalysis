library(ggplot2)

# input folder with GREAT results and output folder for plots
base_dir <- "results/enrichment/great"
outdir <- "results/enrichment/plots"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# get all dataset folders inside the GREAT results folder
datasets <- list.dirs(base_dir, recursive = FALSE)

# safely convert adjusted p-values to -log10 values
safe_log <- function(p) {
  p[p == 0] <- 1e-20
  -log10(p)
}

# store processed results for each dataset
results <- list()

# make one barplot for each GREAT result folder
for (d in datasets) {
  dataset_name <- basename(d)
  file <- file.path(d, "gobp.csv")

  if (!file.exists(file)) next

  # read the GO Biological Process results
  df <- read.csv(file)

  # skip files that do not have the needed columns
  if (!all(c("name", "Binom_Adjp_BH", "Binom_Fold_Enrichment") %in% colnames(df))) {
    next
  }

  # keep the top 15 most significant terms
  df <- df[order(df$Binom_Adjp_BH), ]
  df <- head(df, 15)

  if (nrow(df) == 0) next

  # prepare values for plotting
  df$name <- as.character(df$name)
  df$term_factor <- factor(df$name, levels = rev(df$name))
  df$neg_log_p <- safe_log(df$Binom_Adjp_BH)

  # create barplot for this dataset
  p <- ggplot(df, aes(Binom_Fold_Enrichment, term_factor, fill = neg_log_p)) +
    geom_col() +
    scale_fill_gradient(low = "peachpuff", high = "darkred", name = "-log10(adj p)") +
    labs(title = dataset_name, x = "Fold Enrichment", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 8)
    )

  # save the barplot as a png file
  ggsave(
    file.path(outdir, paste0("barplot_", dataset_name, ".png")),
    p,
    width = 9,
    height = 6,
    dpi = 150
  )

  # save this dataset result for the comparison heatmap
  results[[dataset_name]] <- df
}

# stop if no valid GREAT outputs were found
if (length(results) == 0) {
  stop("No valid GREAT result files found.")
}

# collect top GO terms across all datasets
top_terms <- unique(unlist(lapply(results, function(df) head(df$name, 10))))

hm_list <- list()

# build a combined table for the heatmap
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

# combine all heatmap rows into one dataframe
hm <- do.call(rbind, hm_list)

# set factor order so the plot is easier to read
hm$dataset <- factor(hm$dataset, levels = names(results))
hm$term <- factor(hm$term, levels = rev(top_terms))

# make heatmap comparing enriched GO terms across datasets
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

# save the comparison heatmap
ggsave(
  file.path(outdir, "comparison_heatmap.png"),
  p,
  width = 12,
  height = 14,
  dpi = 150
)