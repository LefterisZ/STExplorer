# ---------------------------------------------------------------------------- #
# Run GWPCA --------------------------------------------------------------------
# ---------------------------------------------------------------------------- #
## Select the sample you would like to perform a GWPCA analysis
sfe <- getSFE(msfe, "V1_2")
SpatialExperiment::imgData(sfe) <- NULL
sfe <- sfe[,1:200]
## Get the gene names that are going to be evaluated
vars = top_hvgs[["V1_2"]]
## Set a fixed bandwidth
## bw is an important parameter as it defines the neighbourhood for which the
##  PCA will be calculated. The distance is measured in ultra-high resolution
##  image pixels. The default is 6x the diameter of the Visium spot. Make sure
##  to adjust it if it is too large or too small for your setting.
bw = 3*sfe@metadata[["spotDiameter"]][["V1_2"]][["spot_diameter_fullres"]]
## Set the number of components to be retained
k = 20
## Set the kernel to be used
kernel = "gaussian"
## Set the Minkowski distance power: p = 2 --> Euclidean
p = 2
## Is the bandwidth adaptive?: No because spots are fixed
adaptive = FALSE
## Cross-Validate GWPCA?
cv = TRUE
## Calculate PCA scores?
scores = TRUE
## Run a robust GWPCA?
robust = FALSE
## Make a cluster for parallel computing (otherwise GWPCA is slow!)
my.cl <- makeClusterGWPCA(type = "FORK")

pcagw <- gwpcaSTE(sfe = sfe,
                  assay = "logcounts",
                  vars = vars,
                  p = p,
                  k = k,
                  bw = bw,
                  kernel = kernel,
                  adaptive = adaptive,
                  scores = scores,
                  robust = robust,
                  cv = cv,
                  future = FALSE,
                  strategy = "cluster",
                  workers = my.cl,
                  verbose = TRUE)

# Extract local scores ---------------------------------------------------------
# Convert the list to a 3D array
local_scores_3d <- array(unlist(pcagw$gwpca.scores), dim = c(2513, 20, 2513))

# Check the dimensions of the resulting array
dim(local_scores_3d)

# ---------------------------------------------------------------------------- #
# MEAN LOCAL SCORE -------------------------------------------------------------
# ---------------------------------------------------------------------------- #
local_scores_2d <- t(apply(local_scores_3d, 1, rowMeans))

# local_scores_2d <- as.matrix(local_scores_3d)[, 1:(dim(local_scores_3d)[2])]

umap <- uwot::umap(local_scores_2d) %>% as.data.frame()
colnames(umap) <- c("UMAP1", "UMAP2")
umap$annotation <- sfe$annotation
ggplot(umap) + geom_point(aes(x = UMAP1, y = UMAP2, colour = annotation))

dist <- dist(local_scores_2d, method = "manhattan")
dbscan_clusters <- dbscan::dbscan(dist, eps = 1, minPts = 5)
map <- data.frame(geometry = colGeometry(sfe),
                  cluster = as.factor(dbscan_clusters$cluster))
ggplot(data = map) +
  geom_sf(aes(geometry = geometry, fill = cluster))

dist <- dist(local_scores_2d, method = "manhattan")
dist <- as.matrix(dist)
colnames(dist) <- rownames(dist) <- colnames(sfe)
k = 6
knn <- dbscan::kNN(dist, k)
knn_edges <- as.vector(t(knn$id))
knn_sources <- rep(1:nrow(knn$id), each = k)
knn_edge_list <- cbind(knn_sources, knn_edges)
graph <- igraph::graph_from_edgelist(knn_edge_list, directed = FALSE)
clust <- igraph::cluster_leiden(graph,
                        objective_function = "modularity",
                        n_iterations = 4,
                        resolution_parameter = 0.5)
igraph::membership(clust)
igraph::V(graph)$leid <- as.factor(igraph::membership(clust))
col.No = length(clust)
colour.values <- getColours(col.No)
map <- data.frame(geometry = colGeometry(sfe),
                  cluster = igraph::V(graph)$leid)
ggplot(data = map) +
  geom_sf(aes(geometry = geometry, fill = cluster)) +
  scale_fill_manual(values = colour.values)

# ---------------------------------------------------------------------------- #
# MULTI GRAPH ------------------------------------------------------------------
# ---------------------------------------------------------------------------- #
# Assuming local_scores_3d is a 3D array: locations x PCs x locations
# bandwidth b is defined
b <- 11  # bandwidth value

global_graph <- igraph::make_empty_graph(directed = FALSE)

for (i in 1:dim(local_scores_3d)[3]) {
  local_matrix <- local_scores_3d[,,i]

  # Filter locations based on bandwidth
  distances <- apply(local_matrix, 1, function(row) sqrt(sum(row^2)))
  neighbor_indices <- which(distances < b)
  local_subset <- local_matrix[neighbor_indices, ]

  # Create local graph
  local_graph <- graph_from_adjacency_matrix(as.matrix(dist(local_subset)), mode = "undirected", weighted = TRUE)

  # Add local graph to global graph
  global_graph <- union(global_graph, local_graph)
}

# Perform clustering on the global graph
clusters <- igraph::cluster_leiden(graph,
                                   objective_function = "modularity",
                                   n_iterations = 4,
                                   resolution_parameter = 0.3)

igraph::membership(clusters)
# igraph::V(global_graph)$leid <- as.factor(igraph::membership(clusters))
col.No = length(clusters)
colour.values <- getColours(col.No)
map <- data.frame(geometry = colGeometry(sfe),
                  cluster = as.factor(igraph::membership(clusters)))
ggplot(data = map) +
  geom_sf(aes(geometry = geometry, fill = cluster)) +
  scale_fill_manual(values = colour.values)

# ---------------------------------------------------------------------------- #
# MULTI GRAPH: physical dist ---------------------------------------------------
# ---------------------------------------------------------------------------- #
# Assuming local_scores_3d is a 3D array: locations x PCs x locations
# bandwidth b is defined
b <- bw/1.5  # bandwidth value

global_graph <- igraph::make_empty_graph(directed = FALSE)

dMat <- sfe@metadata$dMat$euclidean

for (i in 1:dim(local_scores_3d)[3]) {
  message("Working on location: ", i, "/",dim(local_scores_3d)[3])
  local_matrix <- local_scores_3d[,,i]
  rownames(local_matrix) <- colnames(sfe)

  # Filter locations based on bandwidth
  # Fetch pairwise distances
  distances <- dMat[i,]

  # Create binary mask matrix
  mask <- distances <= b

  # Subset local matrix
  local_subset <- local_matrix[mask,]

  # Get distance matrix
  dist <- as.matrix(dist(local_subset))

  # Create local graph
  local_graph <- igraph::graph_from_adjacency_matrix(dist,
                                                     mode = "undirected",
                                                     weighted = TRUE,
                                                     add.colnames = NULL)

  # Add local graph to global graph
  global_graph <- union(global_graph, local_graph)
}

# Plot global graph
igraph::plot.igraph(global_graph,
                    vertex.label = NA,
                    vertex.size = 3,
                    vertex.color = "black",
                    edge.color = "skyblue")

# Perform clustering on the global graph
clusters <- igraph::cluster_leiden(global_graph,
                                   objective_function = "modularity",
                                   n_iterations = 4,
                                   resolution_parameter = 0.2,
                                   weights = NULL)

igraph::membership(clusters)
igraph::V(global_graph)$leid <- as.factor(igraph::membership(clusters))
col.No = length(clusters)
colour.values <- getColours(col.No)
clust_df <- data.frame(cluster = as.factor(igraph::membership(clusters))) %>%
  rownames_to_column(var = "Barcode")
map <- data.frame(geometry = colGeometry(sfe),
                  annotation = colData(sfe)$annotation) %>%
  rownames_to_column(var = "Barcode") %>%
  left_join(clust_df, by = "Barcode") %>%
  column_to_rownames(var = "Barcode")
ggplot(data = map) +
  geom_sf(aes(geometry = geometry, fill = annotation)) +
  geom_sf(aes(geometry = geometry, colour = cluster), fill = NA, linewidth = 1) +
  scale_fill_manual(values = getColours(length(unique(map$annotation)))) +
  scale_colour_manual(values = "black") +
  theme_void()

# ---------------------------------------------------------------------------- #
# MULTI GRAPH: custom union function -------------------------------------------
# ---------------------------------------------------------------------------- #
# To modify the code to sum the weights of existing edges when performing the union of the global graph and local graphs, you can create a custom union function that handles the edge weights. Here is the updated code:

# Assuming local_scores_3d is a 3D array: locations x PCs x locations
# bandwidth b is defined
b <- bw/1.5  # bandwidth value

global_graph <- igraph::make_empty_graph(directed = FALSE)

dMat <- sfe@metadata$dMat$euclidean

custom_union <- function(graph1, graph2) {
  combined_graph <- igraph::union(graph1, graph2)

  # Get the edge attributes of both graphs
  edge_attr1 <- unlist(igraph::edge_attr(graph1))
  names(edge_attr1) <- attr(E(graph1), "vnames")
  edge_attr2 <- unlist(igraph::edge_attr(graph2))
  names(edge_attr2) <- attr(E(graph2), "vnames")

  # Combine the edge attributes by summing the weights
  combined_edge_attr <- lapply(names(edge_attr1), function(attr) {
    ifelse(attr %in% names(edge_attr2),
           edge_attr1[[attr]] + edge_attr2[[attr]],
           edge_attr1[[attr]])
  })

  # Set the combined edge attributes to the combined graph
  igraph::edge_attr(combined_graph) <- combined_edge_attr

  return(combined_graph)
}

for (i in 1:dim(local_scores_3d)[3]) {
  message("Working on location: ", i, "/", dim(local_scores_3d)[3])
  local_matrix <- local_scores_3d[,,i]
  rownames(local_matrix) <- colnames(sfe)

  # Filter locations based on bandwidth
  # Fetch pairwise distances
  distances <- dMat[i,]

  # Create binary mask matrix
  mask <- distances <= b

  # Subset local matrix
  local_subset <- local_matrix[mask,]

  # Get distance matrix
  dist <- as.matrix(dist(local_subset))

  # Create local graph
  local_graph <- igraph::graph_from_adjacency_matrix(dist,
                                                     mode = "undirected",
                                                     weighted = TRUE,
                                                     add.colnames = NULL)

  # Add local graph to global graph using the custom union function
  global_graph <- custom_union(global_graph, local_graph)
}


# The main change in the code is the introduction of a custom union function called `custom_union()`. This function takes two graphs (`graph1` and `graph2`) as input and performs the following steps:
#
# 1. It combines the two graphs using `igraph::union()` to create a `combined_graph`.
# 2. It retrieves the edge attributes (weights) of both `graph1` and `graph2` using `igraph::edge_attr()`.
# 3. It combines the edge attributes by summing the weights for edges that exist in both graphs. If an edge exists only in one graph, its weight is used as is.
# 4. It sets the combined edge attributes to the `combined_graph` using `igraph::edge_attr()`.
# 5. It returns the `combined_graph`.
#
# The rest of the code remains the same, except that the `union()` function is replaced with the `custom_union()` function when adding the local graph to the global graph.
#
# With this modification, when the local graphs are added to the global graph, the weights of existing edges will be summed, resulting in a global graph where the edge weights represent the sum of weights from all the local graphs.


# Plot global graph
igraph::plot.igraph(global_graph,
                    vertex.label = NA,
                    vertex.size = 2,
                    vertex.color = "black",
                    edge.color = "skyblue")

# ---------------------------------------------------------------------------- #
# MULTI GRAPH: GNNs ------------------------------------------------------------
# ---------------------------------------------------------------------------- #
#  Here is a simple implementation of a Graph Neural Network (GNN) in R using the `igraph` and `torch` packages:

library(igraph)
library(torch)

# Create a sample graph
g <- graph(c("A" -> "B", "A" -> "C", "B" -> "D", "C" -> "D"))

# Convert the graph to an adjacency matrix
adj_matrix <- as.matrix(get.adjacency(g))

# Define the GNN model
gnn_model <- nn_module(
  "GNN",
  initialize = function() {
    self$conv1 <- nn_conv_graph(16, 16)
    self$conv2 <- nn_conv_graph(16, 16)
    self$fc <- nn_linear(16, 16)
  },
  forward = function(x, edge_index) {
    x <- torch_relu(self$conv1(x, edge_index))
    x <- torch_relu(self$conv2(x, edge_index))
    x <- self$fc(x)
    return(x)
  }
)

# Create a tensor for the node features
node_features <- torch_tensor(rnorm(4 * 16), dtype = torch_float32)

# Create a tensor for the edge index
edge_index <- torch_tensor(t(adj_matrix), dtype = torch_long)

# Move the model and tensors to the GPU (if available)
device <- torch_device("cuda" if torch_cuda_is_available() else "cpu")
gnn_model <- gnn_model$to(device)
node_features <- node_features$to(device)
edge_index <- edge_index$to(device)

# Define the loss function and optimizer
criterion <- nn_mse_loss()
optimizer <- optim_adam(gnn_model$parameters, lr = 0.01)

# Train the GNN model
for (epoch in 1:100) {
  optimizer$zero_grad()
  out <- gnn_model(node_features, edge_index)
  loss <- criterion(out, node_features)
  loss$backward()
  optimizer$step()
  cat(paste("Epoch:", epoch, "Loss:", loss$item(), "\n"))
}

# Use the trained GNN model to make predictions
predictions <- gnn_model(node_features, edge_index)

# This implementation defines a simple GNN model with two graph convolutional layers and a fully connected layer. The model takes in node features and an edge index as input, and outputs a tensor with the same shape as the input node features.
#
# The `torch` package is used to create tensors and perform computations on the GPU (if available). The `igraph` package is used to create and manipulate the graph.
#
# Note that this is a simplified example, and you may need to modify the architecture, hyperparameters, and training algorithm to suit your specific problem.
#
# Also, keep in mind that GNNs can be computationally expensive, especially for large graphs. You may need to use a more efficient implementation or optimize the model for performance.

# ---------------------------------------------------------------------------- #
# DO NOT TRY THIS AT HOME ------------------------------------------------------
# ---------------------------------------------------------------------------- #
umap <- uwot::umap(local_scores_2d) %>% as.data.frame()
colnames(umap) <- c("UMAP1", "UMAP2")
k = 6
knn <- dbscan::kNN(umap, k)
knn_edges <- as.vector(t(knn$id))
knn_sources <- rep(1:nrow(knn$id), each = k)
knn_edge_list <- cbind(knn_sources, knn_edges)
graph <- igraph::graph_from_edgelist(knn_edge_list, directed = FALSE)
clust <- igraph::cluster_leiden(graph,
                                objective_function = "modularity",
                                n_iterations = 4,
                                resolution_parameter = 0.1)
igraph::membership(clust)
igraph::V(graph)$leid <- as.factor(igraph::membership(clust))
col.No = length(clust)
colour.values <- getColours(col.No)
map <- data.frame(geometry = colGeometry(sfe),
                  cluster = igraph::V(graph)$leid)
ggplot(data = map) +
  geom_sf(aes(geometry = geometry, fill = cluster)) +
  scale_fill_manual(values = colour.values)

# ---------------------------------------------------------------------------- #
# DO NOT TRY THIS AT HOME 2-----------------------------------------------------
# ---------------------------------------------------------------------------- #
data <- local_scores_3d[,,1000]
k = 6
knn <- dbscan::kNN(data, k)
knn_edges <- as.vector(t(knn$id))
knn_sources <- rep(1:nrow(knn$id), each = k)
knn_edge_list <- cbind(knn_sources, knn_edges)
graph <- igraph::graph_from_edgelist(knn_edge_list, directed = FALSE)
clust <- igraph::cluster_leiden(graph,
                                objective_function = "modularity",
                                n_iterations = 4,
                                resolution_parameter = 0.1)
igraph::membership(clust)
igraph::V(graph)$leid <- as.factor(igraph::membership(clust))
col.No = length(clust)
colour.values <- getColours(col.No)
map <- data.frame(geometry = colGeometry(sfe),
                  cluster = igraph::V(graph)$leid)
ggplot(data = map) +
  geom_sf(aes(geometry = geometry, fill = cluster)) +
  scale_fill_manual(values = colour.values) +
  geom_sf(data = map[950,], aes(geometry = geometry), fill = "yellow")





