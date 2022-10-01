library(ggpubr)
library(RColorBrewer)

# Prepare data for MAUP examples
cor.input <- vst_df
row_vars <- rowVars(assay(vst))
select <- order(row_vars, decreasing = TRUE)[seq_len(500)]
cor.input <- cor.input[select]

# Calculate spearman correlations
cor.test(cor.input$ENSMUSG00000015090, cor.input$ENSMUSG00000052305,
         method = "spearman")

ggscatter(cor.input, x = "ENSMUSG00000015090", y = colnames(cor.input[,2:13]),
          combine = TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ENSMUSG00000015090",
          ylab = "")

ggsave(file.path(graphDir, "MAUP_cor_ENSMUSG00000015090_panel.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Calculate multiple correlations and save in df
cor.multi <- cor(cor.input, cor.input, method = "spearman")
diag(cor.multi) <- 0 #set diagonal to zero

cor.multi.highest <- data.frame(gene_X = rownames(cor.multi),
                                gene_Y = colnames(cor.multi)[max.col(cor.multi)]) #get the pairs with the highest correlation in two columns

for (i in 1:10){
    p <- ggscatter(cor.input, 
                   x = cor.multi.highest[i, 1], 
                   y = cor.multi.highest[i, 2], 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman", 
                   xlab = cor.multi.highest[i, 1], 
                   ylab = cor.multi.highest[i, 2], 
                   ggtheme = my_theme)
    
    ggsave(file.path(graphDir, paste0("MAUP_cor_",cor.multi.highest[i, 1], "_", cor.multi.highest[i, 2], ".pdf")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
    
    print(p)
}


# Add spots in Vertical and Horizontal zones
polygons$MAUP.groupsV <- c(rep(1:6, each = 200))[1:nrow(polygons)]

polygons <- polygons %>%
    mutate(MAUP.groupsH = if_else(pixel_y >= 400, 1,
                                  if_else(pixel_y < 400 & pixel_y >= 350, 2,
                                          if_else(pixel_y < 350 & pixel_y >= 300, 3,
                                                  if_else(pixel_y < 300 & pixel_y >= 250, 4,
                                                          if_else(pixel_y < 250 & pixel_y >= 200, 5, 6))))))

# Plot Vertical zone (zone 1)
ggplot() +
    geom_sf(data = polygons$geom_pol, 
            colour = "black", aes(fill = factor(polygons$MAUP.groupsV))) +
    labs(title = paste("MAUP Zonation Vertical"),
         subtitle = "",
         fill = "Vertical\nzones") + 
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme

ggsave(file.path(graphDir, "MAUP_cor_zone-v.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Plot Horizontal zone (zone 2)
ggplot() +
    geom_sf(data = polygons$geom_pol, 
            colour = "grey30", aes(fill = factor(polygons$MAUP.groupsH))) +
    labs(title = paste("MAUP Zonation Horizontal"),
         subtitle = "",
         fill = "Horizontal\nzones") + 
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme

ggsave(file.path(graphDir, "MAUP_cor_zone-h.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Add zones to the correlation input df
cor.input <- cbind(cor.input, 
                   MAUP.groupsV = polygons$MAUP.groupsV, 
                   MAUP.groupsH = polygons$MAUP.groupsH)

# Group vertically and then horizontally. Each time get a set of gene names
cor.input.groupedV <- cor.input %>% 
    group_by(MAUP.groupsV) %>%
    select(-c("MAUP.groupsH"))

genesV <- colnames(cor.input.groupedV)[1:500]

cor.input.groupedH <- cor.input %>% 
    group_by(MAUP.groupsH) %>%
    select(-c("MAUP.groupsV"))

genesH <- colnames(cor.input.groupedH)[1:500]

# Get gene expression sums or means for each zone
maup_summarise_at <- function(data, gene, method){
    if (method == "mean"){
        summarise_at(data, vars(gene), list(gene = mean))
    } else {
        summarise_at(data, vars(gene), list(gene = sum))
    }
}

cor.input.meansV <- lapply(genesV, maup_summarise_at, data = cor.input.groupedV, method = "mean")
cor.input.meansV <- reduce(cor.input.meansV, full_join, by = "MAUP.groupsV")
colnames(cor.input.meansV)[2:501] <- genesV

cor.input.meansH <- lapply(genesH, maup_summarise_at, data = cor.input.groupedH, method = "mean")
cor.input.meansH <- reduce(cor.input.meansH, full_join, by = "MAUP.groupsH")
colnames(cor.input.meansH)[2:501] <- genesH

# Calculate correlations for the same genes as aggregated in zones
for (i in 1:10){
    p <- ggscatter(cor.input.sumsH, 
                   x = cor.multi.highest[i, 1], 
                   y = cor.multi.highest[i, 2], 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   title = "Zones Horizontal",
                   xlab = cor.multi.highest[i, 1], 
                   ylab = cor.multi.highest[i, 2],
                   ggtheme = my_theme)
    
    ggsave(file.path(graphDir, paste0("MAUP_cor_",cor.multi.highest[i, 1], "_", cor.multi.highest[i, 2], "_zoneH.pdf")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
    
    print(p)
}

# Add spots in a grid made from the vertical and horizontal zones
polygons <- polygons %>%
    unite(MAUP.groups, MAUP.groupsV, MAUP.groupsH)

n <- length(unique(polygons$MAUP.groups))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Plot Grid
ggplot() +
    geom_sf(data = polygons$geom_pol, 
            colour = "grey30", aes(fill = factor(polygons$MAUP.groups))) +
    scale_fill_manual(values = sample(col_vector, 35)) +
    labs(title = paste("MAUP Zonation Grid"),
         subtitle = "",
         fill = "Groups") + 
    xlab("X coordinates (pixels)") + 
    ylab("Y coordinates (pixels)") + 
    my_theme

ggsave(file.path(graphDir, "MAUP_cor_zone-grid.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)


# Add zones to the correlation input df
cor.input <- cbind(cor.input, 
                   MAUP.groups = polygons$MAUP.groups)

# Group by grid box and get a set of gene names
cor.input.grouped <- cor.input %>% 
    group_by(MAUP.groups) %>%
    select(-c("MAUP.groupsH", "MAUP.groupsV"))

genesG <- colnames(cor.input.grouped)[1:500]

# Get gene expression sums for each grid box
cor.input.means <- lapply(genesG, maup_summarise_at, data = cor.input.grouped)
cor.input.means <- reduce(cor.input.means, full_join, by = "MAUP.groups")
colnames(cor.input.means)[2:501] <- genesG

# Calculate correlations for the same genes as aggregated in zones
for (i in 1:10){
    p <- ggscatter(cor.input.means, 
                   x = cor.multi.highest[i, 1], 
                   y = cor.multi.highest[i, 2], 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   title = "Zone Grid",
                   xlab = cor.multi.highest[i, 1], 
                   ylab = cor.multi.highest[i, 2],
                   ggtheme = my_theme)
    
    ggsave(file.path(graphDir, paste0("MAUP_cor_",cor.multi.highest[i, 1], "_", cor.multi.highest[i, 2], "_zoneG-mean.pdf")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
    
    print(p)
}

# Correlations for Grid zonation as panel
ggscatter(cor.input.means, x = "ENSMUSG00000015090", y = colnames(cor.input.means[,3:14]),
          combine = TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ENSMUSG00000015090",
          ylab = "")
ggsave(file.path(graphDir, "MAUP_cor_ENSMUSG00000015090_panel_zoneG.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Correlations for Vertical zonation as panel
ggscatter(cor.input.meansV, x = "ENSMUSG00000015090", y = colnames(cor.input.meansV[,3:14]),
          combine = TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ENSMUSG00000015090",
          ylab = "")
ggsave(file.path(graphDir, "MAUP_cor_ENSMUSG00000015090_panel_zoneV.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Correlations for Horizontal zonation as panel
ggscatter(cor.input.meansH, x = "ENSMUSG00000015090", y = colnames(cor.input.meansH[,3:14]),
          combine = TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ENSMUSG00000015090",
          ylab = "")
ggsave(file.path(graphDir, "MAUP_cor_ENSMUSG00000015090_panel_zoneH.pdf"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)

# Create manually p0, p1, p2 plots. ggarrange them for each gene as p3 and p4
# Finally ggarrange p3 and p4 to get 6 plots. For p3 add titles, but remove 
# x-axis labels and add y-axis label only to the leftmost plot. For p4 remove
# titles, add x-axis label to the middle plot and y-axis label to the leftmost.
# Arrange with ggarrange p3 on top and p4 at the bottom.
p0 <- ggscatter(cor.input, 
          x = "ENSMUSG00000015090", 
          y = c("ENSMUSG00000069917"), 
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          title = "",
          xlab = "", 
          ylab = "ENSMUSG00000069917",
          ggtheme = my_theme)

p4 <- ggarrange(p0, p1, p2, ncol = 3)
ggarrange(p3, p4, nrow = 2)

ggsave(file.path(graphDir, "MAUP_cor_ENSMUSG00000015090_panel_zoneGV.tiff"),
       width = grDevices::dev.size(units = "in")[1],
       height = grDevices::dev.size(units = "in")[2],
       units = "in",
       dpi = 400)


