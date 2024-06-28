# Load data ----
samplesPr <- c("~/Documents/Projects_1D/Visium_Prostate/data/spaceranger_outs/Prostate_Canc_hs/PrsCncA1_Results",
               "~/Documents/Projects_1D/Visium_Prostate/data/spaceranger_outs/Prostate_Canc_hs/PrsCncTest_Results")
sample_idPr <- c("PrsCncA1","PrsCncTest")

sfePr <- read10xVisiumSFE(samples = samplesPr,
                          sample_id = sample_idPr,
                          type = "sparse",
                          data = "raw",
                          images = "lowres",
                          style = "W",
                          zero.policy = TRUE)
colData(sfePr)

# Prepare data ----
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfePr)$gene_name)
sfePr <- addPerLocQC(sfePr, gTruth = NULL, assay = "counts", MARGIN = 2, subsets = list(mito = is_mito))
sfePr <- addGeometries(sfePr, samples = samplesPr, sample_id = sample_idPr, res = "fullres")
sfePr <- addPerGeneQC(sfePr, assay = "counts", version = NULL, mirror = NULL)

## Add fake annotations
# For some reason if it is not a tibble then it doesn't merge right.
ground_truthPr <- as.data.frame(colData(sfePr)) %>%
  mutate(annotation = sample(c("A", "B", "C", "D", "E", "F"), ncol(sfePr), replace = TRUE)) %>%
  dplyr::filter(in_tissue == TRUE) %>%
  dplyr::select(Barcode, sample_id, annotation) %>%
  remove_rownames() %>%
  as_tibble()
merger <-  merge(colData(sfePr), DataFrame(ground_truthPr),
                 by = c("Barcode", "sample_id"),
                 all = TRUE)
merger$index2 = as.numeric(sub("spot_", "", merger$index))
merger <- merger[order(merger$index2),]
merger$index2 <- NULL

# Add into colData
colData(sfePr)$annotation <- merger$annotation
rm(merger)


## Add fake cell counts
cellCountPr <- as.data.frame(colData(sfePr)) %>%
  mutate(cellCount = sample(c(1:20), ncol(sfePr), replace = TRUE)) %>%
  dplyr::filter(in_tissue == TRUE) %>%
  dplyr::select(Barcode, sample_id, cellCount) %>%
  remove_rownames() %>%
  as_tibble()
merger <-  merge(colData(sfePr), DataFrame(cellCountPr),
                 by = c("Barcode", "sample_id"),
                 all = TRUE)
merger$index2 = as.numeric(sub("spot_", "", merger$index))
merger <- merger[order(merger$index2),]
merger$index2 <- NULL

# Add into colData
colData(sfePr)$cellCount <- merger$cellCount
rm(merger)

# Plot spatial coordinates without annotations ----
plotQC_spots(sfePr, type = "spot", sample_id = TRUE, in_tissue = FALSE, colours = NULL)
plotQC_spots(sfePr, type = "spot", sample_id = NULL, in_tissue = TRUE)
plotQC_spots(sfePr, type = "hex", sample_id = TRUE, in_tissue = FALSE, colours = c("#3C5338", "#FF9999"))
plotQC_spots(sfePr, type = "hex", sample_id = "PrsCncA1", in_tissue = TRUE)

# Plot spatial coordinates with annotations ----
plotQC_spotsAnnotation(sfe = sfePr, type = "spot", sample_id = TRUE)
plotQC_spotsAnnotation(sfe = sfePr, type = "spot", sample_id = NULL)
plotQC_spotsAnnotation(sfe = sfePr, type = "hex", sample_id = "PrsCncA1")

# Plot manual annotation or spots with tissue image ----
plotQC_tissueImg(sfePr, res = "lowres", type = "spot", sample_id = TRUE, annotate = TRUE, alpha = 0.3)
plotQC_tissueImg(sfePr, res = "lowres", type = "spot", sample_id = TRUE, annotate = FALSE, alpha = 0.3)
plotQC_tissueImg(sfePr, res = "lowres", type = "hex", sample_id = NULL, annotate = TRUE, alpha = 0.3)
plotQC_tissueImg(sfePr, res = "lowres", type = "hex", sample_id = NULL, annotate = FALSE, alpha = 0.3)
plotQC_tissueImg(sfePr, res = "lowres", type = "none", sample_id = "PrsCncA1", annotate = FALSE)

# Keep in-tissue locations ----
sfePr <- filterInTissue(sfePr, sample_id = TRUE)

# Library sizes ----
## Density and histogram
plotQC_hist(sfePr, metric = "libsize")
plotQC_hist(sfePr, metric = "libsize", limits = c(3000, 36500))
plotQC_hist(sfePr, metric = "libsize", limits = c(3000, 36500),
            hist_args = list(bins = 100),
            dens_args = list(alpha = 0.5,
                             adjust = 0.5,
                             fill = "#F0AAA8"),
            vline_args = list(colour = "blue",
                             linetype = "dashed"))
## Scatter plot library sizes vs number of cells
plotQC_scat(sfePr, metric = "libsize")
## Select library size threshold
sfePr <- setQCthresh_LibSize(sfePr, sample_id = TRUE, min_t = 3000, max_t = 37000)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfePr, metric = "libsize", sample_id = TRUE)

# Genes expressed ----
## Density and histogram
plotQC_hist(sfePr, metric = "detected")
plotQC_hist(sfePr, metric = "detected", limits = c(2000, NA))
## Scatter plot genes expressed vs number of cells
plotQC_scat(sfePr, metric = "detected")
## Select threshold
sfePr <- setQCthresh_GenesExpr(sfePr, sample_id = TRUE, min_t = 2000, max_t = NA)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfePr, metric = "detected", sample_id = TRUE)

# Mitochondrial expression ----
## Density and histogram of %
plotQC_hist(sfePr, metric = "mito")
plotQC_hist(sfePr, metric = "mito", limits = c(NA, 22))
## Scatter plot % of mitochondrial expression vs number of cells
plotQC_scat(sfePr, metric = "mito")
## Select threshold
sfePr <- setQCthresh_Mito(sfePr, sample_id = TRUE, min_t = NA, max_t = 22)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfePr, metric = "mito", sample_id = TRUE)

# Cells in each spot ----
## Density and histogram of the number
plotQC_hist(sfePr, metric = "cellCount")
## Scatter plot genes expressed vs number of cells
plotQC_scat(sfePr, metric = "detected")
## Select threshold
sfePr <- setQCthresh_CellCount(sfePr, sample_id = TRUE, min_t = NA, max_t = 25)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfePr, metric = "cellCount", sample_id = TRUE)

# NAs in annotation ----
## Select threshold
sfePr <- setQCthresh_NAs(sfePr, sample_id = TRUE)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfePr, metric = "NAs", sample_id = TRUE)

# QC Discard Locations ----
## Set the combined filtering threshold using the QC metrics
sfePr <- setQCtoDiscard_loc(sfePr, sample_id = TRUE, filters = TRUE)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfePr, metric = "discard", sample_id = TRUE)

# Apply location-level QC threshold ----
sfePr <- applyQCthresh_loc(sfePr, sample_id = TRUE)

# Compute Library Size factors ----
sfePr <- computeLibSizeFactors(sfePr)
## Density and histogram of library sizes
plotQC_sizeFactors(sfePr)

# Generate a metaSFE object ----
meta_sfePr <- asMetaSFE(sfePr)

# Normalise counts ----
meta_sfePr <- normaliseCounts(msfe = meta_sfePr)

# Calculating extra gene QC metrics ----
meta_sfePr <- perGeneLogMean(meta_sfePr)

# Zero expression genes ----
meta_sfePr <- setQCthresh_ZeroExpr(meta_sfePr)

# Lowly expressed (noise?!) genes ----
meta_sfePr <- setQCthresh_LowLogMean(meta_sfePr)

# QC discard Features ----
## Set the combined filtering threshold using the QC metrics
meta_sfePr <- setQCtoDiscard_feat(meta_sfePr, filters = TRUE)

# Apply gene-level QC threshold ----
meta_sfePr <- applyQCthresh_feat(meta_sfePr)

# Model Gene Variance ----
dec_Pr <- modelGeneVariance(meta_sfePr, method = "Var")
dec_Pr[[1]]

# Get Top HVGs ----
top_hvgs_Pr <- getTopHighVarGenes(dec_Pr,
                               var.field = "bio",
                               prop = 0.5,
                               var.threshold = 0,
                               fdr.threshold = 0.1)
plotGeneVariance(dec = dec_Pr, hvgs = top_hvgs_Pr)

# Add neighbour graphs ----
meta_sfePr <- addSpatialNeighGraphs(meta_sfePr,
                                    sample_id = TRUE,
                                    type = "knearneigh",
                                    style = "W",
                                    distMod = "raw",
                                    k = 6)
plotNeighbourGraph(meta_sfePr, sample_id = TRUE,
                   res = "lowres", plotImage = TRUE)
plotNeighbourGraph(meta_sfePr, sample_id = "PrsCncA1",
                   res = "lowres", plotImage = TRUE)
plotNeighbourGraph(meta_sfePr, sample_id = "PrsCncA1",
                   res = "lowres", plotImage = FALSE)

# Calculate a simple distance matrix ----
meta_sfePr <- addDistMat(meta_sfePr, p = 2)

# Functional Clustering ----
msigdb <- getMSigDBData("Homo sapiens")
viewCollections()
t2g <- getTerm2Gene(msig_data = msigdb, cat = "C2", subcat = "CP")
gsea_map <- gwpca_FunctionalClustering(gwpca = pcagw,
                                       pc = 1,
                                       genes_no = 1,
                                       NES = 1.5,
                                       minGSSize = 5,
                                       pvalueCutoff = 0.25,
                                       TERM2GENE = t2g,
                                       pAdjustMethod = "fdr",
                                       scoreType = "std",
                                       nPermSimple = 10000,
                                       mc.cores = 4)
plotGWPCA_FuncCLust(gsea_map, count = 5, legend = "right")

# FGWC ----
sfe_nmf <- fgwc_nmf(sfe, sample_id = "JBO019", top_hvgs = top_hvgs)
fgwc_param <- fgwc_params(algorithm = "classic", ncluster = 5)
fgwc <- fgwcSTE(sfe, "JBO019", data = sfe_nmf, dMetric = "euclidean", fgwc_param = fgwc_param)
plotFGWC_single(fgwc = fgwc, m_sfe = sfe, sample_id = "JBO019")
plotFGWC_multi(fgwc = fgwc, m_sfe = sfe, sample_id = "JBO019")
heatmap <- plotFGWC_heatmap(fgwc = fgwc, m_sfe = sfe, sample_id = "JBO019", markers = markers, cluster_no = 3)
plotFGWC_subClust(heatmap = heatmap, k = 5, clust = 4, m_sfe = sfe, sample_id = "JBO019")
plotFGWC_subHeatmap(heatmap = heatmap, k = 5, markers = markers, m_sfe = sfe, sample_id = "JBO019", cluster_no = 1)
