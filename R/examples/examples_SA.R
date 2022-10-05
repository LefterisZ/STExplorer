## Transform counts to vst ----
vst <- varianceStabilizingTransformation(dds)
vst_df <- as.data.frame(t(assay(vst))) # transpose and transform to df

## Select most variable genes ----
row_vars <- rowVars(assay(vst))
select <- order(row_vars, decreasing = TRUE)[seq_len(500)]

## Prepare for SA calculations ----
# Get spot names
nb_names <- polygons$Barcode

# Get normalised counts from the top 500 most variable genes
counts.SA <- counts(dds, normalized = TRUE) %>% # export normalised counts
    t() %>% # transpose for downstream SA calculation
    as.data.frame() %>% # make it a df
    #.[nb_names,] %>% # Re-order rows to match the polygon object row order
    .[select] # select the 500 most variable genes

# Get VST transformed counts from the top 500 most variable genes
counts.SA.vst <- vst_df[select] %>%
    .[nb_names,]

# Find neighbours 
neighbours_knn <- knn2nb(knearneigh(polygons$geom_cntd, k = 6), 
                         row.names = nb_names)

# Calculate neighbour weights in different ways
neighbours_w_exp <- nb2listwdist(neighbours_knn, polygons$geom_cntd,
                                 type = "idw", style = "raw", alpha = 2,
                                 zero.policy = TRUE)

neighbours_w <- nb2listw(neighbours)

## Run SA: Moran's I (global/local) ----
# Global
moran(x = counts.SA$ENSMUSG00000000740, 
      listw = neighbours_w_exp, 
      n = length(neighbours_w_exp$neighbours),
      S0 = Szero(neighbours_w_exp))


# Local
test.gene.exp <- counts.SA$ENSMUSG00000000740 %>%
    setNames(rownames(counts.SA))

moran.plot(x = test.gene.exp, 
           listw = neighbours_w_exp,
           labels = FALSE,
           main = "ENSMUSG00000000740",
           xlab = "Gene expression",
           ylab = "Spatially lagged gene expression",
           zero.policy = TRUE,
           spChk = TRUE)

moran.local <- localmoran(x = test.gene.exp, 
                          listw = neighbours_w_exp, 
                          spChk = TRUE)

moran.local.perm <- localmoran_perm(x = test.gene.exp, 
                          listw = neighbours_w_exp, 
                          nsim = 999,
                          spChk = TRUE)

# Run Moran's I tests with MC perm. and Z-score ----
# Z-score
moran.test(x = counts.SA$ENSMUSG00000000740, 
           listw = neighbours_w_exp)

# Monte Carlo permutations
moran.mc(x = counts.SA$ENSMUSG00000000740, 
         listw = neighbours_w_exp, 
         nsim = 999)

#### RUN Multiple Moran's ####
moran.n <- length(neighbours_w_exp$neighbours)
moran.S0 <- Szero(neighbours_w_exp)

moran.multi <- apply(counts.SA, 2, function(x) moran(x, 
                                                     neighbours_w_exp,
                                                     n = moran.n,
                                                     S0 = moran.S0))
moran.multi.df <- moran.multi %>% 
    reduce(bind_rows) %>% 
    as.data.frame() # reduce list of lists to a df

rownames(moran.multi.df) <- colnames(counts.SA) # give row names

moran.P.SA <- moran.multi.df[order(moran.multi.df$I, decreasing = TRUE),]

# Create a test counts df to map the gene expression
test.counts <- data.frame(gene.expr.norm = counts.SA$ENSMUSG00000000740,
                          pixel_x = polygons$pixel_x,
                          pixel_y = polygons$pixel_y)

# Map gene expression
ggplot(test.counts) + 
    geom_point(aes(x = pixel_x, y = pixel_y, colour = gene.expr.norm),
               size = 3) + 
    scale_colour_continuous(limits = c(0 , 150)) + 
    labs(title = "ENSMUSG00000000740",
         subtitle = NULL,
         colour = "Normalised counts\ngene expression") +
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme

ggsave(file.path(graphDir, "gene.expression_map_SA.NegativeCtrl.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
ggsave(file.path("~/Downloads", "gene.expression_map_SA.NegativeCtrl.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Map Local Moran statistic
moran.loc.plot.Ii <- data.frame(geometry = polygons$geom_pol,
                               Ii = moran.local[,"Ii"])

ggplot() +
    geom_sf(data = moran.loc.plot.Ii$geometry,
            colour = "grey30",
            aes(fill = moran.loc.plot.Ii$Ii)) + 
    labs(title = "ENSMUSG00000000740",
         subtitle = NULL,
         fill = "Local Moran\nIi") +
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme
ggsave(file.path(graphDir, "moran.local_ENSMUSG00000000740.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
ggsave(file.path("~/Downloads", "moran.local_ENSMUSG00000000740.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Map LISA clusters
quads <- attr(moran.local.perm, "quadr") # get the quadrant attributes
moran.loc.plot.Quads <- data.frame(geometry = polygons$geom_pol,
                                   p.value = moran.local.perm[,"Pr(folded) Sim"],
                                   Quadrant = quads$pysal) %>%
    mutate(colours = Quadrant) %>%
    mutate(colours = case_when(p.value > 0.05 ~ "Not signif.",
                               Quadrant == "Low-Low" ~ "Low-Low",
                               Quadrant == "High-Low" ~ "High-Low",
                               Quadrant == "Low-High" ~ "Low-High",
                               Quadrant == "High-High" ~ "High-High"))

ggplot() +
    geom_sf(data = moran.loc.plot.Quads$geometry,
            colour = "grey30",
            aes(fill = moran.loc.plot.Quads$colours)) +
    scale_fill_manual(values = c("#ff0000", "#ffb3b3", "#0066ff", 
                                 "#ffffff", "#99c2ff")) +
    labs(title = "ENSMUSG00000000740",
         subtitle = NULL,
         fill = "LISA\nclusters") +
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme
ggsave(file.path(graphDir, "moran.local.LISAClust_ENSMUSG00000000740.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
ggsave(file.path("~/Downloads", "moran.local.LISAClust_ENSMUSG00000000740.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
