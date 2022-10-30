library(ggbeeswarm)

# non-package functions needed
source("./R/side_extras/gwpca/gwpca.prop.var.R")
source("./R/side_extras/gwpca/gwpca.plot.prop.vars.multi.R")
source("./R/side_extras/gwpca/gwpca.plot.prop.vars.single.R")
source("./R/side_extras/gwpca/prop.var.R")

dt_top.Gs <- gwpca.topLeading.Gs(gwpca = pca_gw.list$pca_gw.500.20.gau,
                                 pc.no = 1,
                                 genes.n = 3,
                                 sf.geom = polygons$geom_pol,
                                 method = "membership",
                                 gene.names = TRUE,
                                 biomart = biomart.mouse.98,
                                 check.names = FALSE)


pc = 1
genes = 3
method = "membership"
groups = count(unique(as.data.frame(dt_top.Gs$Top_lead_Gs)))

#map the leading genes
ggplot() +
    geom_sf(data = dt_top.Gs$geometry, 
            aes(fill = as.factor(dt_top.Gs$Top_lead_Gs))) +
    geom_point(size = 3) + 
    #scale_fill_manual(values = c4a("wright25", 18)) +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    labs(title = paste0("Leading Genes on PC", pc),
         subtitle = paste0("Top ", genes, " Genes"),
         caption = paste0("Grouping method: ",
                          method,
                          "\nNumber of groups: ",
                          groups),
         fill = "Top Leading Genes\nGroups") +
    my_theme +
    theme(legend.position="right")

ggsave(file.path(graphDir, "gwpca_.500.20.gau_leadingGenes3_PC1-membership.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

#calculate the PTV for multiple Components
props <- gwpca.prop.var(gwpca.obj = pca_gw.list$pca_gw.500.20.gau,
                        n.comp = c(5, 10, 20, 30, 40, 50))

#plot them all together
gwpca.plot.prop.vars.multi(select(props, -c(pixel_x, pixel_y)), theme = my_theme)

ggsave(file.path(graphDir, "gwpca_.500.20.gau_PTV_all.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

#plot PTVs one by one and in a panel
for (n in column_names) {
    gwpca.plot.prop.vars.single(data = select(props, n))
    
    name = sub("Comps_", "", n)
    
    ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_PTV_", name,".pdf")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
}

plot_list <- lapply(colnames(select(props, -c(pixel_x, pixel_y))),
                    gwpca.plot.prop.vars.single, 
                    data = select(props, -c(pixel_x, pixel_y)), ylab = NULL)

plot_list <- setNames(plot_list, colnames(select(props, -c(pixel_x, pixel_y))))

egg::ggarrange(plots = plot_list, nrow = 2, ncol = 3)

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_PTV_panel.pdf")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
# ggpubr::annotate_figure(fig,
#                         left = text_grob("Percantage of Total Variation (PTV)", 
#                                          color = "green", rot = 90))


# Map the PTV for a specific selection of components
for(i in c(5, 10, 20, 30)) {
    print(i)
    comps <- sprintf("Comps_%02d", i)
    ptv.map <- dplyr::select(props, all_of(c(comps, "geometry")))
    
    ggplot() + 
        geom_sf(data = ptv.map$geometry,
                aes(fill = ptv.map[,1])) +
        scale_fill_viridis_c(option = "magma", limits = c(0, 100)) +
        xlab("X coordinates (pixels)") +
        ylab("Y coordinates (pixels)") +
        labs(title = "Percantage of Total Variation\n(PTV)",
             fill = paste0("PTV of ", i, "\n components")) +
        my_theme
}


# calculate the discrepancies
data.mat <- as.matrix(inputPCAgw@data)
discrepancy <- gwpca.cv.contrib(data.mat, coordinates(inputPCAgw), 
                                bw = 6*spot_diameter(spatialDir), 
                                adaptive = TRUE, dMat = dist.Mat)

# plot the discrepancies in a box plot
discrepancy_df <- data.frame(disc = discrepancy)

ggplot(pivot_longer(discrepancy_df, col = "disc"),
       aes(x = name, y = value)) + 
    geom_boxplot(fill = "#D1E5F0", colour = "#2166AC", 
                 outlier.colour = "red", outlier.size = 2) + 
    geom_jitter(col = "#EF8A62", size = 2, width = 0.3, alpha = 0.8) +
    coord_flip() +
    ggtitle("Local PC Discrepancy") +
    xlab(NULL) + 
    my_theme

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_discreps.box.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# map the discrepancies
dt <- inputPCAgw@data %>%
    mutate(disc = discrepancy,
           geometry = polygons$geom_pol)

disc.map <- dplyr::select(dt, all_of(c("disc", "geometry")))

ggplot() + 
    geom_sf(data = disc.map$geometry, 
            aes(fill = disc.map$disc)) + 
    scale_fill_viridis_c(option = "inferno") +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    labs(title = "Local PC Discrepancy",
         fill = "Discrepancy\nscore") +
    my_theme

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_discreps.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# create the input data table for the outlier plot
inputPCAgw.outlier <- vst_df[,select] %>% # select top 500 variable genes
    as.data.frame() %>%                      # make it a df
    .[nb_names,]                             # order rows

# create the biomart for mouse
biomart.mouse <- create_biomart("mouse")

# plot the heatmap to visualise the genes that make this location an outlier
tiff(file.path(graphDir, paste0("gwpca_.500.20.gau_discrep.max.outlierPlot.tiff")),
     width = grDevices::dev.size(units = "in")[1],
     height = grDevices::dev.size(units = "in")[2],
     units = "in",
     res = 400)

gwpca.plot.outlier(inputPCAgw.outlier,
                   bw = 3*spot_diameter(spatialDir),
                   focus = which(discrepancy_df == max(discrepancy_df)),
                   dMat = dist.Mat,
                   show.vars = "top",
                   mean.diff = 1,
                   gene.names = TRUE, 
                   biomart = biomart.mouse.98,
                   show.data = FALSE,
                   check.names = FALSE,
                   scale = "row",
                   cutree_cols = 5,
                   cutree_rows = 6,
                   color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000)))

dev.off()
