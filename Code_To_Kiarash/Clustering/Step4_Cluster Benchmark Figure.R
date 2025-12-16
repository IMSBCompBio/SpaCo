library(dplyr)
library(ggplot2)
load(file = "./Data/cluster_results.Rdata")

combined_df <- do.call(rbind,
                       list(
                         liver_cluster_res,
                         anterior_cluster_res,
                         posterior_cluster_res
                       ))
combined_df$metric <- sub("\\_.*", "", rownames(combined_df))
combined_df$method <- sub("\\d.*", "",
                          sub(".*_", "", rownames(combined_df)))
combined_df$tissue <- rep(c("Liver", "Anterior", "Posterior"),
                          each = nrow(anterior_cluster_res))
colnames(combined_df)[1] <- "value"
combined_df <-
  combined_df %>% filter(!startsWith(rownames(combined_df), "i"))

plt <-
  ggplot(combined_df %>% filter(method != "svg"),
         aes(x = tissue, y = value, fill = method)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 1,
    color = NA,
    median.colour = "black",
    median.linewidth = 0.2
  ) +
  scale_fill_manual(
    values = c(
      "spaco" = "darkolivegreen3",
      "pca" = "goldenrod2",
      "spatialpca" = "skyblue",
      "svg" = "orange"
    )
  ) +
  facet_wrap(~ metric, scales = "free_y") +
  coord_cartesian(ylim = c(0, 0.8)) +  # adjust range as needed
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.1, 0.15)
  ) +
  ylab("Score") +
  xlab("Tissue")
ggsave("./Cluster_Benchmark.svg",
       plt,
       width = 12, height = 6)
#
# cleveland_plot_data <-
#   combined_df %>%
#   group_by(tissue, method, metric) %>%
#   summarize(median_value = median(value), .groups = "drop")
#
# ggplot()
