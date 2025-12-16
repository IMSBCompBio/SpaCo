library(dplyr)
library(ggplot2)
res_dir <- "/home/rstudio/SpaCo_paper_code/SVG/data/res/"
res_dir_p <-
  "./python/results/svg_res/"
fp_names <-
  readRDS("/home/rstudio/SpaCo_paper_code/SVG/data/twin_names.RDS")
tp_names <-
  readRDS("/home/rstudio/SpaCo_paper_code/SVG/data/totwin.RDS")
all_genes <- c(fp_names, tp_names)
n_rads <- 5
n_bm <- 10
bm_settings <- rep(1:n_rads, each = n_bm)
n_settings <- length(bm_settings)
specifs <- seq(0, 1, by = 0.01)
spaco_scores <- rep(0, n_settings)
get_ranks <- function(res_data, score_col,
                      all_genes) {
  rank_df <- data.frame(gene = rownames(res_data),
                        rank = order(res_data[, score_col])) %>%
    filter(gene %in% all_genes)
  missing_genes <- all_genes[which(!all_genes %in% rank_df$gene)]
  rank_df[missing_genes, ] <- data.frame(missing_genes,
                                         0)
  rank_df %>% mutate(rank = order(rank))
}
get_roc_data_helper <- function(res_data,
                                score_col,
                                fp = fp_names,
                                tp = tp_names) {
  all_genes <- c(fp, tp)
  n_genes <- length(all_genes)
  n_genes_sep <- length(fp)
  rank_df <- get_ranks(res_data, score_col, all_genes)
  roc_data <- t(sapply(1:n_genes, function(idx) {
    pos <- which(rank_df$rank < idx)
    fpr <-
      sum(fp %in% rank_df[pos, 1]) / n_genes_sep
    tpr <-
      sum(tp %in% rank_df[pos, 1]) / n_genes_sep
    return(c(fpr, tpr))
  }))
  # roc_data <- rbind(c(0, 0),
  #                   roc_data,
  #                   c(1, 1))
  roc_data
}
get_roc_data <- function(j,
                         rad_i,
                         type,
                         fp = fp_names,
                         tp = tp_names) {
  res <- get_res(rad_i, j, type)
  score_col <- switch(
    type,
    "SPACO" = "score",
    "SPARKX" = "combinedPval",
    "HeartSVG" = "rank",
    "SpaGFT" = "svg_rank"
  )
  get_roc_data_helper(res, score_col, fp, tp)
}
get_res <- function(rad_i, j, type) {
  if (type == "SpaGFT") {
    res <- read.csv(paste0(res_dir_p, "r_", rad_i, "_run_", j, ".csv"))
    rownames(res) <- res$X
  } else {
    res <- readRDS(paste0(res_dir, "res_r_", rad_i, "_run_", j, ".RDS"))
    res <-
      res[[switch(
        type,
        "SPACO" = 2,
        "SPARKX" = 3,
        "HeartSVG" = 4
      )]]
    if (type == "HeartSVG") {
      rownames(res) <- res$gene
    }
    if (type == "SPACO") {
      res$score <- -res$score
    }
  }
  res
}
get_plt_data <- function(i_rad, type,
                         fp = fp_names,
                         tp = tp_names) {
  perf_mat <- data.frame(matrix(nrow = 0, ncol = 2))
  colnames(perf_mat) <- c("fpr", "tpr")
  for (j in 1:n_bm) {
    new_perf <- get_roc_data(j, rad_i = i, type = type)
    colnames(new_perf) <- c("fpr", "tpr")
    new_perf_eq <- new_perf %>% as.data.frame() %>%
      group_by(fpr) %>% summarise(tpr = mean(tpr))
    perf_mat <- rbind(perf_mat,
                      new_perf_eq)
  }
  perf_mat <- cbind(perf_mat,
                    Run = rep(1:n_bm, each = length(fp) + length(tp) + 2))
  uq_perf_mat <- perf_mat %>% as.data.frame() %>%
    group_by(fpr) %>% summarise(mean_tpr = mean(tpr))
  rbind(c(0, 0),
        uq_perf_mat,
        c(1, 1))
}
type_vec <- c("SPACO", "SPARKX", "HeartSVG", "SpaGFT")
plt_data <- matrix(nrow = 0, ncol = 4)
for (i in 1:n_rads) {
  rad_perfs <- lapply(type_vec,
                      function(type)
                        get_plt_data(i, type))
  names(rad_perfs) <- type_vec
  for (type_idx in seq_along(rad_perfs)) {
    plt_data <- rbind(plt_data,
                      cbind(rad_perfs[[type_idx]],
                            type_vec[type_idx],
                            i))
  }
}
colnames(plt_data) <- c("FPR", "TPR", "Algorithm", "Radius")
plt_data <- plt_data %>%
  mutate(Radius = factor(Radius,
                         labels = paste0("r = ",
                                         c(2, 5, 10, 20, 40))))
# saveRDS(plt_data, "./python/results/svg_res/plt_data.RDS")
plt <- ggplot(plt_data %>% filter(Radius != "r = 40"),
              aes(x = FPR, y = TPR, color = Algorithm)) +
  geom_line() +
  scale_color_manual(
    values = c(
      "SPACO" = "darkolivegreen3",
      "SPARKX" = "violetred2",
      "HeartSVG" = "blue4",
      "SpaGFT" = "tan4"
    )
  ) +
  facet_wrap(~ Radius, ncol = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      color = "black",
      # border color
      fill = "white",
      # background fill
      linewidth = 0.25
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.25
    ),
    strip.text = element_text(margin = margin(0.5, 4, 2, 4),   # top, right, bottom, left (in points)
                              color = "black"),
    # plot.margin = margin(0,0,-0.025, 0, "cm"),
    legend.position = c(0.9, 0.075)
  )
# plt
ggsave("./SVG_Result.svg",
       plt,
       width = 8, height = 8)
