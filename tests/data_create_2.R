# Load data ----
sampleDir <- c("./data/test_data/Visium_Human_Liver/Healthy/JBO022_Results",
               "./data/test_data/Visium_Human_Liver/Steatotic/JBO019_Results")
sampleNames <- c("JBO022", "JBO019")

sfe_2samples <- read10xVisiumSFE(samples = sampleDir,
                          sample_id = sampleNames,
                          type = "HDF5",
                          data = "filtered",
                          images = "lowres",
                          style = "W",
                          zero.policy = TRUE)
usethis::use_data(sfe_2samples, overwrite = TRUE)
saveIn <- c("./data/test_data/Visium_Human_Liver")
saveRDS(sfe_2samples, file = file.path(saveIn, "sfe_2samples.rds"))

sfeLivH <- read10xVisiumSFE(samples = sampleDir[2],
                           sample_id = sampleNames[2],
                           type = "HDF5",
                           data = "filtered",
                           images = "lowres",
                           style = "W",
                           zero.policy = TRUE)
usethis::use_data(sfeLivH, overwrite = TRUE)

sfeLivS <- read10xVisiumSFE(samples = sampleDir[1],
                            sample_id = sampleNames[1],
                            type = "HDF5",
                            data = "filtered",
                            images = "lowres",
                            style = "W",
                            zero.policy = TRUE)
usethis::use_data(sfeLivS, overwrite = TRUE)

sfe_raw <- sfeLiv[,order(as.data.frame(spatialCoords(sfeLiv))$pxl_col_in_fullres)[2:300]]
usethis::use_data(sfe_raw, overwrite = TRUE)

# Create MSFE objects ----
msfe <- MetaSpatialFeatureExperiment()
msfe <- addSFE(msfe, sfe = sfeLivH)
msfe <- addSFE(msfe, sfe = sfeLivS)
usethis::use_data(msfe, overwrite = TRUE)

msfe_2sample <- MetaSpatialFeatureExperiment()
msfe_2sample <- addMultipleSFE(msfe_2sample, sfe = sfe_2samples)
usethis::use_data(msfe_2sample, overwrite = TRUE)

colData(sfe)

# Prepare data ----
is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$gene_name)
sfe <- addPerLocQC(sfe, gTruth = NULL, assay = "counts", MARGIN = 2, subsets = list(mito = is_mito))
sfe <- addGeometries(sfe, samples = samplesPr, sample_id = sample_idPr, res = "fullres")
sfe <- addPerGeneQC(sfe, assay = "counts", version = NULL, mirror = NULL)

## Add fake annotations
# For some reason if it is not a tibble then it doesn't merge right.
ground_truthPr <- as.data.frame(colData(sfe)) %>%
  mutate(annotation = sample(c("A", "B", "C", "D", "E", "F"), ncol(sfe), replace = TRUE)) %>%
  dplyr::filter(in_tissue == TRUE) %>%
  dplyr::select(Barcode, sample_id, annotation) %>%
  remove_rownames() %>%
  as_tibble()
merger <-  merge(colData(sfe), DataFrame(ground_truthPr),
                 by = c("Barcode", "sample_id"),
                 all = TRUE)
merger$index2 = as.numeric(sub("spot_", "", merger$index))
merger <- merger[order(merger$index2),]
merger$index2 <- NULL

# Add into colData
colData(sfe)$annotation <- merger$annotation
rm(merger)


## Add fake cell counts
cellCountPr <- as.data.frame(colData(sfe)) %>%
  mutate(cellCount = sample(c(1:20), ncol(sfe), replace = TRUE)) %>%
  dplyr::filter(in_tissue == TRUE) %>%
  dplyr::select(Barcode, sample_id, cellCount) %>%
  remove_rownames() %>%
  as_tibble()
merger <-  merge(colData(sfe), DataFrame(cellCountPr),
                 by = c("Barcode", "sample_id"),
                 all = TRUE)
merger$index2 = as.numeric(sub("spot_", "", merger$index))
merger <- merger[order(merger$index2),]
merger$index2 <- NULL

# Add into colData
colData(sfe)$cellCount <- merger$cellCount
rm(merger)

# Plot spatial coordinates without annotations ----
plotQC_spots(sfe, type = "spot", sample_id = TRUE, in_tissue = FALSE, colours = NULL)
plotQC_spots(sfe, type = "spot", sample_id = NULL, in_tissue = TRUE)
plotQC_spots(sfe, type = "hex", sample_id = TRUE, in_tissue = FALSE, colours = c("#3C5338", "#FF9999"))
plotQC_spots(sfe, type = "hex", sample_id = "PrsCncA1", in_tissue = TRUE)

# Plot spatial coordinates with annotations ----
plotQC_spotsAnnotation(sfe = sfe, type = "spot", sample_id = TRUE)
plotQC_spotsAnnotation(sfe = sfe, type = "spot", sample_id = NULL)
plotQC_spotsAnnotation(sfe = sfe, type = "hex", sample_id = "PrsCncA1")

# Plot manual annotation or spots with tissue image ----
plotQC_tissueImg(sfe, res = "lowres", type = "spot", sample_id = TRUE, annotate = TRUE, alpha = 0.3)
plotQC_tissueImg(sfe, res = "lowres", type = "spot", sample_id = TRUE, annotate = FALSE, alpha = 0.3)
plotQC_tissueImg(sfe, res = "lowres", type = "hex", sample_id = NULL, annotate = TRUE, alpha = 0.3)
plotQC_tissueImg(sfe, res = "lowres", type = "hex", sample_id = NULL, annotate = FALSE, alpha = 0.3)
plotQC_tissueImg(sfe, res = "lowres", type = "none", sample_id = "PrsCncA1", annotate = FALSE)

# Keep in-tissue locations ----
sfe <- filterInTissue(sfe, sample_id = TRUE)

# Library sizes ----
## Density and histogram
plotQC_hist(sfe, metric = "libsize")
plotQC_hist(sfe, metric = "libsize", limits = c(3000, 36500))
plotQC_hist(sfe, metric = "libsize", limits = c(3000, 36500),
            hist_args = list(bins = 100),
            dens_args = list(alpha = 0.5,
                             adjust = 0.5,
                             fill = "#F0AAA8"),
            vline_args = list(colour = "blue",
                              linetype = "dashed"))
## Scatter plot library sizes vs number of cells
plotQC_scat(sfe, metric = "libsize")
## Select library size threshold
sfe <- setQCthresh_LibSize(sfe, sample_id = TRUE, min_t = 3000, max_t = 37000)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "libsize", sample_id = TRUE)

# Genes expressed ----
## Density and histogram
plotQC_hist(sfe, metric = "detected")
plotQC_hist(sfe, metric = "detected", limits = c(2000, NA))
## Scatter plot genes expressed vs number of cells
plotQC_scat(sfe, metric = "detected")
## Select threshold
sfe <- setQCthresh_GenesExpr(sfe, sample_id = TRUE, min_t = 2000, max_t = NA)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "detected", sample_id = TRUE)

# Mitochondrial expression ----
## Density and histogram of %
plotQC_hist(sfe, metric = "mito")
plotQC_hist(sfe, metric = "mito", limits = c(NA, 22))
## Scatter plot % of mitochondrial expression vs number of cells
plotQC_scat(sfe, metric = "mito")
## Select threshold
sfe <- setQCthresh_Mito(sfe, sample_id = TRUE, min_t = NA, max_t = 22)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "mito", sample_id = TRUE)

# Cells in each spot ----
## Density and histogram of the number
plotQC_hist(sfe, metric = "cellCount")
## Scatter plot genes expressed vs number of cells
plotQC_scat(sfe, metric = "detected")
## Select threshold
sfe <- setQCthresh_CellCount(sfe, sample_id = TRUE, min_t = NA, max_t = 25)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "cellCount", sample_id = TRUE)

# NAs in annotation ----
## Select threshold
sfe <- setQCthresh_NAs(sfe, sample_id = TRUE)
## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "NAs", sample_id = TRUE)

# QC Discard Locations ----
## Set the combined filtering threshold using the QC metrics
sfe <- setQCtoDiscard_loc(sfe, sample_id = TRUE, filters = TRUE)

## Check putative spatial patterns of removed spots
plotQC_filtered(sfe, metric = "discard", sample_id = TRUE)

# Apply location-level QC threshold ----
sfe <- applyQCthresh_loc(sfe, sample_id = TRUE)

# Compute Library Size factors ----
sfe <- computeLibSizeFactors(sfe)
## Density and histogram of library sizes
plotQC_sizeFactors(sfe)

# Generate a metaSFE object ----
meta_sfePr <- asMetaSFE(sfe)

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

# SA ----

