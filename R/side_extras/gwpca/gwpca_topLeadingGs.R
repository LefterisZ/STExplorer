library(ggbeeswarm)
# library(egg)
# library(ggpubr)

# non-package functions needed
source("./R/side_extras/gwpca/gwpca.prop.var.R")
source("./R/side_extras/gwpca/gwpca.plot.prop.vars.multi.R")
source("./R/side_extras/gwpca/gwpca.plot.prop.vars.single.R")
source("./R/side_extras/gwpca/prop.var.R")

local.loadings <- pca_gw.list$pca_gw.500.20.gau$loadings[,,1]
dt <- apply(abs(local.loadings), 1, order, decreasing = TRUE)
dt[1:10, 1:5]
dt_top10 <- dt[1:10,]
dt_top10 <- apply(dt_top10, 2, function(x) colnames(local.loadings)[x])
dt_top10 <- t(dt_top10) %>%
    as.data.frame() %>% 
    unite("Top10_lead_Gs", 1:10, sep = ";", remove = TRUE)
dt_top10[1:5,]

dt_top10 <- data.frame(dt_top10) %>%
    mutate(pixel_x = polygons$pixel_x,
           pixel_y = polygons$pixel_y)

#map the leading genes
ggplot(dt_top10, 
       aes(x = pixel_x, y = pixel_y, colour = as.factor(Top10_lead_Gs))) + 
    geom_point(size = 3) + 
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    ggtitle("Leading Genes on PC1") +
    my_theme +
    theme(legend.position="none")

ggsave(file.path(graphDir, "gwpca_.500.20.gau_leadingGenes10_PC1.tiff"),
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


# map the PTV for a specific selection of components
ggplot(data = select(props, c(Comps_30, pixel_x, pixel_y)), 
       aes(x = pixel_x, y = pixel_y, colour = Comps_30)) + 
    geom_point(size = 3) +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    ggtitle("Percantage of Total Variation (PTV)") +
    my_theme


# calculate the discrepancies
data.mat <- as.matrix(inputPCAgw@data)
discrepancy <- gwpca.cv.contrib(data.mat, coordinates(inputPCAgw), 
                                bw =6*spot_diameter(spatialDir), 
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
    mutate(pixel_x = polygons$pixel_x,
           pixel_y = polygons$pixel_y,
           disc = discrepancy)

ggplot() + 
    geom_point(data = select(dt, c(disc, pixel_x, pixel_y)), 
               aes(x = pixel_x, y = pixel_y, colour = disc),
               size = 3.5) + 
    geom_sf(data = outlier.coords, colour = "darkorange", size = 3) + 
    geom_point(data = filter(dt, disc > 200000),
               aes(x = pixel_x, y = pixel_y), colour = "red", size = 3) + 
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    ggtitle("Local PC Discrepancy") +
    my_theme

ggsave(file.path(graphDir, paste0("gwpca_.500.20.gau_discreps.tiff")),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)
