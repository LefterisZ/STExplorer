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
    geom_point(size = 1) +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    ggtitle("Leading Gene on PC1") +
    my_theme +
    theme(legend.position="none")

ggsave(file.path(graphDir, "gwpca_.500.20.gau_leadingGenes10_PC1.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

#calculate the PTV for multiple Components
props <- gwpca.prop.var(gwpca.obj = pca_gw.list$pca_gw.500.20.gau,
                        n.comp = c(5, 10, 20, 30, 40, 50))
#plot them all together
gwpca.plot.prop.vars.multi(select(props, -c(pixel_x, pixel_y)), theme = my_theme)

ggsave(file.path(graphDir, "gwpca_.500.20.gau_PTV_all.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# plot PTVs one by one
gwpca.plot.prop.single(data = select(props, Comps_05), 
                       violcol = "skyblue", theme = my_theme)

# map the PTV for a specific selection of components
ggplot(props, 
       aes(x = pixel_x, y = pixel_y, colour = props)) + 
    geom_point(size = 1) +
    xlab("X coordinates (pixels)") +
    ylab("Y coordinates (pixels)") +
    ggtitle("Percantage of Total Variation (PTV)") +
    my_theme
