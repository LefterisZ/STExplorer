library(ggpubr)

# Prepare data for MAUP examples
cor.input <- vst_df
row_vars <- rowVars(assay(vst))
select <- order(row_vars, decreasing = TRUE)[seq_len(500)]
cor.input <- cor.input[select]

# Calculate spearman correlations
cor.test(cor.input$ENSMUSG00000015090, cor.input$ENSMUSG00000052305,
         method = "spearman")

ggscatter(cor.input, x = "ENSMUSG00000015090", y = colnames(cor.input[,2:11]),
          combine = TRUE,
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "ENSMUSG00000015090")

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
    my_theme +
    theme(legend.position = "right")

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

genesH <- colnames(cor.input.groupedV)[1:500]

# Get gene expression means for each zone
maup_summarise_at <- function(data, gene){
    summarise_at(data, vars(gene), list(gene = sum))
}

cor.input.sumsV <- lapply(genesV, maup_summarise_at, data = cor.input.groupedV)
cor.input.sumsV <- reduce(cor.input.sumsV, full_join, by = "MAUP.groupsV")
colnames(cor.input.sumsV)[2:501] <- genesV

cor.input.sumsH <- lapply(genesH, maup_summarise_at, data = cor.input.groupedH)
cor.input.sumsH <- reduce(cor.input.sumsH, full_join, by = "MAUP.groupsH")
colnames(cor.input.sumsH)[2:501] <- genesH

# Calculate correlations for the same genes as aggregated in zones
for (i in 1:10){
    p <- ggscatter(cor.input.sumsV, 
                   x = cor.multi.highest[i, 1], 
                   y = cor.multi.highest[i, 2], 
                   add = "reg.line", conf.int = TRUE, 
                   cor.coef = TRUE, cor.method = "spearman",
                   title = "Zones Vertical",
                   xlab = cor.multi.highest[i, 1], 
                   ylab = cor.multi.highest[i, 2],
                   ggtheme = my_theme)
    
    ggsave(file.path(graphDir, paste0("MAUP_cor_",cor.multi.highest[i, 1], "_", cor.multi.highest[i, 2], "_zoneV.pdf")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
    
    print(p)
}


