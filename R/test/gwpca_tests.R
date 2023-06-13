#' A script to test different gwpca bits
#' 
# rm(inputPCAgw, vars, bw, k, x, wt, data, col.nm, var.idx, kernel,dMat,wpcaGPT,wpca,.x,.wtox,wox,posn,wt.median.1,medians,centered,mids,result,obj,dp.locat,dp.n,assay,scores,temp, use,w,robust,pcafun,gw_weights,i,j)
# rm(a,b,d,d1,elocat,local.PV,pca.res,res.df,scores.i,scores.ii,weighted_centered,adaptive,CV,CV2,DM.given,ep.given,ep.n,longlat,len.var,optimized_lower,optimized_upper,original_lower,original_upper,ox,p,score,theta,
#    var.n,var.names,var.nms,win.var.pc1,gwpca.cv.contr.ste,n,missing_names_list1,missing_names_list2.list1,list2)
# rm(.obj, .assay, .elocat, .vars, .p, .k, .bw, .kernel,.adaptive, .scores, .robust, .cv, i, ..pcafun, ..x,..dMat,
#    ..w,..d,..wt,..gwpca.scores,..scores.i,..k,..bw,..adaptive,..kernel,..scores,future,strategy,workers,verbose,.verbose, common_names,
#    cv, cvContrib, missing_names_list2, name, score.contrib, v, temp2)

# sfex <- sfe[,(colData(sfe)$layer %in% c("WM")) & (spatialCoords(sfe)[,"pxl_col_in_fullres"] < 4000)]
# sfex <- sfe[,spatialCoords(sfe)[,"pxl_col_in_fullres"] < 6000]

# plot(spatialCoords(sfe))
# ncol(sfe)

# inputPCAgw <- SpatialPointsDataFrame(spatialCoords(sfe), as.data.frame(t(as.matrix(assay(sfe, "counts")))), match.ID = TRUE)
# inputPCAgw <- inputPCAgw[colSums(inputPCAgw@data) > 2000]
# # data <- inputPCAgw
# n <- ncol(inputPCAgw@data)
# vars = colnames(inputPCAgw@data)[1:n]
# bw = 6*sfex@metadata[["spotDiameter"]][["151673"]][["spot_diameter_fullres"]]
# k = 5
# kernel = "gaussian"
# p = 1
# adaptive = FALSE
# cv = TRUE
# scores = FALSE
# robust = FALSE
# dMat = gw.dist(coordinates(inputPCAgw))
# colnames(dMat) <- rownames(inputPCAgw@data)
# rownames(dMat) <- rownames(inputPCAgw@data)
# 
# pcagw_original <- gwpca(data = data, 
#                         vars = vars, 
#                         p = p, 
#                         k = k, 
#                         bw = bw, 
#                         kernel = kernel,
#                         adaptive = adaptive, 
#                         scores = scores, 
#                         robust = robust,
#                         cv = cv)
# 
# my.cl <- parallel::makeCluster(availableCores() - 1, type = 'FORK')
# 
# pcagw_stex <- gwpca.ste(obj = sfex, 
#                         assay = "counts",
#                         vars = vars, 
#                         p = p, 
#                         k = k, 
#                         bw = bw, 
#                         kernel = kernel,
#                         adaptive = adaptive, 
#                         scores = scores, 
#                         robust = robust,
#                         cv = cv,
#                         future = FALSE,
#                         strategy = "cluster",
#                         workers = my.cl,
#                         verbose = TRUE)
# 
# pcagw_ste <- gwpca.ste(obj = sfe, 
#                        assay = "counts",
#                        vars = vars, 
#                        p = p, 
#                        k = k, 
#                        bw = bw, 
#                        kernel = kernel,
#                        adaptive = adaptive, 
#                        scores = scores, 
#                        robust = robust,
#                        cv = cv,
#                        future = FALSE,
#                        strategy = "cluster",
#                        workers = my.cl,
#                        verbose = TRUE)
# 
#  ## Check results
# all.equal(pcagw_original$pca$sdev, pcagw_ste$pca$sdev) # pca slot
# check.equal.list(pcagw_original, pcagw_ste)

## Internal gwpca data structures
# data <- as(data, "data.frame")
# data.x <- as(data, "data.frame")
# col.nm <- colnames(data)
# var.idx <- match(vars, col.nm)[!is.na(match(vars, col.nm))]
# x <- as.matrix(data.frame(data[,var.idx]))
# xx <- as.matrix(data.frame(data.x[,var.idx]))
# wt <- gw.weight(dMat[,1],bw,kernel,adaptive = FALSE)
# 
# ## gwpca.ste internal data structures
# obj = sfex
# assay = "logcounts"
# scores = FALSE
# i = 1
# j = 1
# .obj = sfex
# .assay = assay
# .k = k
# .bw = bw
# .vars = vars
# .p = p
# .kernel = kernel
# .adaptive = adaptive 
# .scores = scores
# .robust = robust
# .cv = cv
################################################################################
################################################################################
################################################################################

# # Original code
# original_upper <- NULL
# original_lower <- NULL
# 
# # Optimized code
# optimized_upper <- NULL
# optimized_lower <- NULL
# 
# # Test parameters
# adaptive <- FALSE
# DM.given <- TRUE
# dMat <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
# p <- 3
# theta <- 1
# longlat <- FALSE
# dp.locat <- matrix(c(0, 0, 1, 1), ncol = 2)
# dp.n <- nrow(dp.locat)
# 
# # Original code
# if (adaptive) {
#   original_upper <- dp.n
#   original_lower <- 2
# } else {
#   if (DM.given) {
#     original_upper <- range(dMat)[2]
#     original_lower <- original_upper / 5000
#   } else {
#     dMat <- NULL
#     if (p == 2) {
#       b.box <- bbox(dp.locat)
#       original_upper <- sqrt((b.box[1, 2] - b.box[1, 1])^2 + (b.box[2, 2] - b.box[2, 1])^2)
#       original_lower <- original_upper / 5000
#     } else {
#       original_upper <- 0
#       for (i in 1:dp.n) {
#         dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, p = p, theta = theta, longlat = longlat)
#         original_upper <- max(original_upper, range(dist.vi)[2])
#       }
#       original_lower <- original_upper / 5000
#     }
#   }
# }
# 
# # Optimized code
# if (adaptive) {
#   optimized_upper <- dp.n
#   optimized_lower <- 2
# } else {
#   if (DM.given) {
#     optimized_upper <- range(dMat)[2]
#     optimized_lower <- optimized_upper / 5000
#   } else {
#     dMat <- NULL
#     if (p == 2) {
#       b.box <- bbox(dp.locat)
#       optimized_upper <- sqrt((b.box[1, 2] - b.box[1, 1])^2 + (b.box[2, 2] - b.box[2, 1])^2)
#       optimized_lower <- optimized_upper / 5000
#     } else {
#       optimized_upper <- sapply(1:dp.n, function(i) {
#         dist.vi <- gw.dist(dp.locat = dp.locat, focus = i, p = p, theta = theta, longlat = longlat)
#         range(dist.vi)[2]
#       })
#       optimized_upper <- max(optimized_upper)
#       optimized_lower <- optimized_upper / 5000
#     }
#   }
# }
# 
# # Compare results
# identical(original_upper, optimized_upper) # Check if upper values are identical
# identical(original_lower, optimized_lower) # Check if lower values are identical

################################################################################
#############################Test GWPCA plots###################################
################################################################################
# plotGWPCA_global(gwpca = gwpca,
#                  comps = 1:10,
#                  type = "scree",
#                  point_args = list(size = 5, colour = "red"),
#                  line_args = list(linewidth = 2, colour = "dodgerblue"))
# 
# plotGWPCA_leadingG(gwpca = gwpca,
#                    comps = 1:4,
#                    type = "single",
#                    arrange = TRUE)
# 
# plotGWPCA_leadingG(gwpca = gwpca,
#                         comps = 1:2,
#                         type = "multi",
#                         arrange = TRUE)

################################################################################
################################################################################
################################################################################
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(tidyverse)
library(scran)
library(scater)
library(ggspavis)
library(sf)
library(spdep)
library(GWmodel)
library(future)
library(doFuture)
library(foreach)
library(progressr)
library(parallel)
library(cols4all)
library(rlist)
library(pheatmap)
# List and source all scripts in the core folder.
file.sources <- list.files(path = "../STExplorer/R/core",
                           pattern = "*.R$",
                           all.files = TRUE,
                           recursive = TRUE,
                           full.names = TRUE,
                           ignore.case = TRUE)
sapply(file.sources, source, local = .GlobalEnv)

sampleDir <- "/Users/b9047753/Documents/Projects_1D/Visium_Liver/data/spaceranger_outs/Human_Liver_Steatotic/Human_Liver_Steatotic_JBO019_Results"
sampleNames <- "JBO019"
sfe <- read10xVisiumSFE(samples = sampleDir, 
                        sample_id = sampleNames, 
                        type = "sparse", 
                        data = "filtered", 
                        images = "lowres", 
                        style = "W", 
                        zero.policy = TRUE)

ground_truth <- read_table("/Users/b9047753/Documents/Projects_1D/Visium_Liver/data/spaceranger_outs/Human_Liver_Steatotic/Human_Liver_Steatotic_JBO019_Results/outs/spatial/spotzonationGroup.txt")

is_mito <- grepl("(^MT-)|(^mt-)", rowData(sfe)$symbol)
sfe <- addPerLocQC(sfe, gTruth = ground_truth, assay = "counts", 2, subsets = list(mito = is_mito))
sfe <- addGeometries(sfe, samples = sampleDir, sample_id = sampleNames, res = "fullres")
sfe <- addPerGeneQC(sfe, assay = "counts", version = 77)

colData(sfe)
rowData(sfe)
colGeometries(sfe)

## Select library size threshold
qc_lib_size <- colData(sfe)$sum < 700
## Check how many spots are filtered out
table(qc_lib_size)
## Add threshold in colData
colData(sfe)$qc_lib_size <- qc_lib_size

## Select expressed genes threshold
qc_detected <- colData(sfe)$detected < 500
## Check how many spots are filtered out
table(qc_detected)
## Add threshold in colData
colData(sfe)$qc_detected <- qc_detected

## Select expressed genes threshold
qc_mito <- colData(sfe)$subsets_mito_percent > 25
## Check how many spots are filtered out
table(qc_mito)
## Add threshold in colData
colData(sfe)$qc_mito <- qc_mito

## Check the number of discarded spots for each metric
apply(cbind(qc_lib_size, qc_detected, qc_mito), 2, sum)
## Combine together the set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito
## Store the set in the object
colData(sfe)$discard <- discard

## Check the spatial pattern of combined set of discarded spots
plotQC(sfe, type = "spots", 
       discard = "discard")

## remove combined set of low-quality spots
sfe <- sfe[, !colData(sfe)$discard]

## add neighbour graph and distance matrix
sfe <- get.spatialNeighGraphs(sfe, sampleNames, type = "knearneigh", style = "W", distMod = "raw", k = 6)
sfe <- addDistMat(sfe, p = 2)
colGraphs(sfe)

## Calculate library size factors
sfe <- computeLibraryFactors(sfe)
## Have a look at the size factors
summary(sizeFactors(sfe))

# calculate logcounts using library size factors
sfe <- logNormCounts(sfe)

# FEATURE SELECTION
is_zero <- rowData(sfe)$total == 0
is_low <- (rowData(sfe)$JBO019.s_mean <= 1.5 & 
             rowData(sfe)$JBO019.sparsity < 0.99)
discard_gs <- is_zero | is_mito | is_low

rowData(sfe)$discard <- discard_gs
# remove mitochondrial genes
sfe <- sfe[!rowData(sfe)$discard, ]
# fit mean-variance relationship
dec <- modelGeneVar(sfe,
                    assay.type = "logcounts")
dec_CV2 <- modelGeneCV2(sfe, 
                        size.factors = colData(sfe)$sizeFactor, 
                        assay.type = "logcounts")
# select top HVGs
top_hvgs <- getTopHVGs(dec, 
                       var.field = "bio", 
                       prop = 0.5,
                       var.threshold = 0,
                       fdr.threshold = 0.05)

top_hvgs_CV2 <- getTopHVGs(dec_CV2, 
                           var.field = "ratio", 
                           prop = 0.1,
                           var.threshold = 0.5,
                           fdr.threshold = 0.05)

sfe

vars = top_hvgs
bw = 6*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]
k = 20
kernel = "gaussian"
p = 2
adaptive = FALSE
cv = TRUE
scores = FALSE
robust = FALSE

my.cl <- parallel::makeCluster(availableCores() - 1, type = 'FORK')

# >>> it returns an error when inside the markdown. Maybe run with verbose = FALSE
pcagw_ste <- gwpcaSTE(obj = sfe, 
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


plotGWPCA_global(gwpca = pcagw_ste,
                 comps = 1:10,
                 type = "scree",
                 point_args = list(size = 3, colour = "red"),
                 line_args = list(linewidth = 1, colour = "dodgerblue"))

pcagw_ste <- gwpca_LeadingGene(gwpca = pcagw_ste, 
                               sfe = sfe, 
                               pc_nos = 1:4, 
                               type = "single", 
                               names = "gene_names")

pcagw_ste <- gwpca_LeadingGene(gwpca = pcagw_ste, 
                               sfe = sfe, 
                               pc_nos = 1:4, 
                               genes_n = 4, 
                               type = "multi", 
                               method = "membership", 
                               names = "gene_names")


plotGWPCA_leadingG(gwpca = pcagw_ste,
                   comps = 1:4,
                   type = "single",
                   arrange = TRUE)

plotGWPCA_leadingG(gwpca = pcagw_ste,
                   comps = 1,
                   type = "multi",
                   arrange = TRUE)

pcagw_ste <- gwpca_PropVar(gwpca = pcagw_ste, n_comp = 2:10, sfe = sfe)

plotGWPCA_ptv(gwpca = pcagw_ste,
              comps = 1:10,
              type = "violin")

plotGWPCA_ptv(gwpca = pcagw_ste,
              comps = 1:5,
              type = "map")
# calculate discrepancies outside of GWPCA
discrepancy <- gwpca_cvSTE(bw = bw,
                            x = as.matrix(t(assay(sfe, "logcounts")[vars,])),
                            dMat = gw.dist(spatialCoords(sfe), p = 1),
                            k = 20,
                            loc = spatialCoords(sfe),
                            verbose = TRUE,
                            cvContrib = TRUE)

plotGWPCA_discr(pcagw_ste, type = "box")
plotGWPCA_discr(pcagw_ste, type = "map")

discrepancy_loc_dt <- getDiscrepancyLocData(sfe = sfe, 
                                            gwpca = pcagw_ste, 
                                            sample_id = "JBO019")

head(discrepancy_loc_dt)
focus <- discrepancy_loc_dt$barcodes[1:2]
bw = 3*sfe@metadata[["spotDiameter"]][["JBO019"]][["spot_diameter_fullres"]]

plotGWPCA_discrHeatmap(sfe = sfe,
                       assay = "logcounts",
                       vars = NULL,
                       focus = focus,
                       dMetric = "euclidean", 
                       sample_id = "JBO019",
                       bw = bw, 
                       mean.diff = 1, 
                       show.vars = "top", 
                       scale = "row", 
                       gene.names = TRUE)

discrepancy_gene_dt <- getDiscrepancyGeneData(sfe = sfe,
                                              assay = "logcounts",
                                              vars = NULL,
                                              focus = focus[2],
                                              dMetric = "euclidean", 
                                              sample_id = "JBO019",
                                              bw = bw, 
                                              mean.diff = 1, 
                                              show.vars = "top",
                                              exportExpression = TRUE)


