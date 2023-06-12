BiocManager::install("clusterProfiler")
BiocManager::install("fgsea")
BiocManager::install("msigdbr")
library(clusterProfiler)
library(fgsea)
library(msigdbr)

pc <- 1
focus <- 1
x <- 1

# Generate term2gene object
msigdbr_show_species() # check available species
msigdbr_collections() # check available collections
## get the collection of interest
msigdbr_all <- msigdbr(species = "Homo sapiens")

gs_cats <- c("H", "C2", "C4", "C5", "C6")
gs_subcats <- c("GO:BP", "CM")

msigdbr_H_df <- msigdbr_all %>%
    filter(gs_cat == gs_cats[1])

msigdbr_C2_df <- msigdbr_all %>%
    filter(gs_cat == gs_cats[2])

msigdbr_C4_df <- msigdbr_all %>%
    filter(gs_cat == gs_cats[3] & gs_subcat == gs_subcats[2])

msigdbr_C5_df <- msigdbr_all %>%
    filter(gs_cat == gs_cats[4] & gs_subcat == gs_subcats[1])

msigdbr_C6_df <- msigdbr_all %>%
    filter(gs_cat == gs_cats[5])

## make t2g data.frame
msigdbr_t2g_H <- msigdbr_H_df %>%
    select(gs_name, human_ensembl_gene) %>%
    as.data.frame()

msigdbr_t2g_C2 <- msigdbr_C2_df %>%
    select(gs_name, human_ensembl_gene) %>%
    as.data.frame()

msigdbr_t2g_CM <- msigdbr_C4_df %>%
    select(gs_name, human_ensembl_gene) %>%
    as.data.frame()

msigdbr_t2g_BP <- msigdbr_C5_df %>%
    select(gs_name, human_ensembl_gene) %>%
    as.data.frame()

msigdbr_t2g_C6 <- msigdbr_C6_df %>%
    select(gs_name, human_ensembl_gene) %>%
    as.data.frame()


msigdbr_t2g <- msigdbr_t2g_C6

# get input dataset from gwpca
inputGSEA <- pcagw_ste$loadings[,,pc]

# Run GSEA for all locations
gsea.list.H <- lapply(1:dim(inputGSEA)[1], function(x){
    message(x, " out of ", dim(inputGSEA)[1])
    # Generate geneList
    geneList <- inputGSEA[x,] %>%
        .[order(inputGSEA[x,], decreasing = TRUE)]
    # geneList #check they are descending order
    
    # Run GSEA
    gsea <- GSEA(geneList, 
                 minGSSize = 5,
                 pvalueCutoff = 0.05,
                 TERM2GENE = msigdbr_t2g,
                 pAdjustMethod = "none",
                 verbose = FALSE)
    
    if (dim(gsea@result)[1] == 0) {
        na <- c(rep(NA, dim(gsea@result)[2])) %>%
            t() %>%
            as.data.frame()
        colnames(na) <- colnames(gsea@result)
        gsea@result <- rbind(gsea@result, na)
    }
    # Add location barcode
    gsea@result$Location <- rownames(inputGSEA)[x]
    gsea@result
})

# Make the dataframe ----
gsea_H <- bind_rows(gsea.list.H)

gsea_H <-  gsea_H %>%
    dplyr::group_by(Location) %>%
    dplyr::arrange(NES, .by_group = TRUE) %>% 
    slice_max(NES, n = 1, with_ties = FALSE) %>%
    column_to_rownames(var = "Location")

gsea_H.map <- data.frame(gsea_H,
                         "geometry" = gwpca$geometry) %>%
    mutate(Pathway = gsub("HALLMARK_", "", ID))


# Plot functional clsutering
gsea.map <- gsea_H.map
col.No <- length(unique(gsea.map$Pathway))
colour.values <- get.colours(col.No)

ggplot() + 
    geom_sf(data = gsea.map$geometry, 
            aes(fill = gsea.map$Pathway),
            colour = "grey30", 
            show.legend = TRUE) + 
    scale_fill_manual(values = colour.values) +
    labs(title = NULL,
         fill = "Hallmark\nGene Sets") + 
    theme(legend.position = "right")
