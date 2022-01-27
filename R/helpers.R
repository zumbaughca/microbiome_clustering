#' Determine k in the kmeans algorithm
#' 
#' Find the optimum number of clusters for the kmeans algorithm based on 
#' total within sum of squared error. Note that this tests clusters
#' between 1 and 20.
#' 
#' @param data Data to be passed to kmeans
#' @return a data frame consisting of the number of clusters and total WSS
#' @export
optimize_kmeans <- function(data) {
  results <- data.frame(clusters = c(0),
                        wss = c(0))
  for (i in 1:20) {
    km <- kmeans(data, centers = i, nstart = 30)
    results <- results %>%
      add_row(clusters = i,
              wss = km$tot.withinss)
  }
  return(results %>% slice(-1))
}

#' Determine the number of clusters in hierarchical clustering
#' 
#' Find the optimum number of clusters for the hclust algorithm based on the
#' maximized Calinski-Harabasz index. Note that this tests clusters between
#' 2 and 20.
#' 
#' @param matrix A (dis) similarity matrix
#' @param tree A hclust object
#' @value A named list containing the number of clusters and the 
#' Calinski-Harabasz index 
#' @export
optimize_hclust <- function(matrix, tree) {
  optimum = list(clusters = 0,
                 ch_index = 0)
  for (i in 2:10) {
    ch <- fpc::calinhara(matrix, cutree(tree, k = i))
    if (ch > optimum$ch_index) {
      optimum$clusters = i
      optimum$ch_index = ch
    }
  }
  optimum
}

#' Determine the center of each cluster in a biplot
#' 
#' This calculates the x and y coordinates for the center of each cluster
#' for data that has been subjected to dimensionality reduction followed by
#' clustering.
#' 
#' @param clusters A vector containing the cluster labels for all 
#' observations
#' @param scores A data frame containing the PCA scores.
#' @value A data frame that contains the cluster, scores, center x coordinate,
#' and center y coordinate for each observation.
#' @export
get_cluster_centroids <- function(clusters, scores) {
  scores %>%
    mutate(cluster = factor(clusters)) %>%
    select(PC1, PC2, cluster) %>%
    pivot_longer(cols = c(PC1, PC2)) %>%
    group_by(cluster, name) %>%
    summarise(center = mean(value)) %>%
    pivot_wider(id_cols = "cluster",
                names_from = "name",
                values_from = "center") %>%
    rename(center_x = PC1,
           center_y = PC2)
}

#' Project the original data on a biplot
#' 
#' Project the original data on to the new Principal Coordinates. This
#' can be used to generate arrows similar to what is observed on a biplot.
#' 
#' @param pcoa An object of class pcoa
#' @param original_data The data used to generate the object passed as pcoa
#' @value A data frame containing the transformed data.
#' @export
project_variables <- function(pcoa, original_data) {
  n <- nrow(original_data)
  scaled <- scale(pcoa$vectors[, 1:2])
  covar <- cov(original_data, scaled)
  eigs <- pcoa$values$Eigenvalues[1:2]
  covar_std <- covar %*% diag((eigs / (nrow(df) - 1)) ^ -0.5) %>%
    as.data.frame() %>%
    rename(PC1 = V1,
           PC2 = V2)
  covar_std
}

#' Get the most influential species
#' 
#' Get the top 5 species ranked by the amount of variance.
#' 
#'  @param cov_table The return value of project_variables
#'  @value A covariance table
#'  @export 
influential_species <- function(cov_table) {
  # Get top 5 for PC1 and PC2
  species <- cov_table %>%
    mutate(species = rownames(.)) %>%
    pivot_longer(cols = c(PC1, PC2)) %>%
    group_by(name) %>%
    slice_max(order_by = value, n = 5) %>%
    ungroup() %>%
    select(species)
  # Return covariance table for them
  cov_table %>%
    filter(rownames(.) %in% species[[1]])
}