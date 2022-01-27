library(readxl)
library(tidyverse)
library(ggrepel)
library(plotly)
library(vegan)
library(ape)
library(fpc)
source('R/helpers.R')

df <- read_xlsx('microbial.xlsx', sheet = "genus.Percentages") %>%
  slice_head(n = 138)
data <- df[, 7:ncol(df)]
data <- data[, -c(which(colSums(data) == 0))] / 100

# Calculate Bray-Curtis dissimilarity matrix and run PCoA
dissim <- vegdist(data, method = "bray")
pca <- pcoa(dissim)

# Get original data projections
covar_std <- project_variables(pca, df[, 7:138]) %>%
  influential_species(.)

scores <- as.data.frame(pca$vectors) %>%
  cbind(df[, "Day"], .)
var_exp <- pca$values$Relative_eig

# Heiarchial clustering
clusters <- hclust(dissim, method = "ward.D")
opt_clusters <- optimize_hclust(dissim, clusters)

scores <- scores %>%
  mutate(cluster = factor(cutree(clusters, k = opt_clusters$clusters)),
         trt = df$Trt) %>%
  filter(!is.na(trt)) %>%
  rename(PC1 = Axis.1,
         PC2 = Axis.2,
         PC3 = Axis.3) %>%
  left_join(., get_cluster_centroids(cutree(clusters, k = 2), .),
            by = "cluster")

# The caption for the plot generated below. This is mostly string magic
# to avoid hard coding the caption in the call to labs() below.
caption <- "Principal coordinate analysis of relative abundance of bacteria. PCoA
      was performed on the Bray-Curtis dissimilarity matrix./nColors represent the results of hierarchical clustering using Ward's method.
      /nCluster number was selected as that which maximized the Calinski-Harabasz
       index./nArrows represent the original data projected on the new axes." %>%
  gsub("\n", "", .) %>%
  str_split(., "\\s+") %>%
  unlist() %>%
  paste(., collapse = " ") %>%
  gsub("/n", "\n", .)

# Generate the final biplot.
plt <- ggplot(data = scores,
              aes(x = PC1,
                  y = PC2)) +
  geom_point(aes(color = cluster,
                 shape = Day),
             size = 2) +
  stat_ellipse(geom = "polygon",
               aes(color = cluster,
                   fill = cluster),
               alpha = 0.5,
               level = 0.67) +
  xlab(paste("PCo1", 
             "(", 
             round(var_exp[1] * 100, digits = 2), 
             "% of variance)", 
             sep = " ")) +
  ylab(paste("PCo2", 
             "(", 
             round(var_exp[2] * 100, digits = 2), 
             "% of variance)", 
             sep = " ")) +
  labs(caption = caption) +
  geom_vline(xintercept = 0,
             color = "firebrick",
             linetype = "dashed") +
  geom_hline(yintercept = 0,
             color = "firebrick",
             linetype = "dashed") +
  geom_segment(data = covar_std,
               x = 0,
               y = 0,
               aes(xend = PC1 / 100 ,
                   yend = PC2 / 100 ),
               arrow = arrow(length = unit(3, "mm")),
               size = 1,
               color = "steelblue") +
  geom_segment(aes(x = center_x,
                   y = center_y,
                   xend = PC1,
                   yend = PC2),
               alpha = 0.25,
               size = 0.2) +
  geom_label_repel(data = covar_std / 100,
             aes(label = rownames(covar_std)),
             nudge_x = 0.25,
             nudge_y = 0.15,
             segment.linetype = "dashed",
             segment.alpha = 0.75,
             segment.color = "steelblue",
             segment.size = 0.75) +
  guides(fill = guide_legend(title = "Cluster"),
         color = guide_legend(title = "Cluster")) +
  theme(plot.caption = element_text(hjust = 0))

plt

ggsave('pcoa_plot.pdf', 
       plot = plt,
       width = 8,
       height = 6,
       units = "in")
