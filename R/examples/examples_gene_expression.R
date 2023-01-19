#' In this script we will plot the expression of certain genes.
#' 
goiExprDir <- file.path(graphDir, "gene_expression_maps")

counts_annot <- annotate_data.frame(counts, add = "gene_name", biomart = biomart.mouse.98)

goi <- c("Dcx", "Omp", "Calb2", "Sox11", "Calb1", "Gad2")

for(name in goi){
    g.expr <- filter(counts_annot, gene_name == name) %>% #get goi expression
        select(-c(1)) %>% #remove ENSGIDs
        t() %>% #transpose
        as.data.frame() %>% 
        rownames_to_column() #move barcodes to column
    
    colnames(g.expr) <- g.expr[1,] #add colnames
    g.expr <- g.expr[-1,] %>% #remove first row since it is now as colnames
        rename("Barcode" = "gene_name", "Norm_expr" = name) %>% #rename the column with the barcodes
        left_join(polygons[,c("Barcode", "geom_pol")]) %>% 
        mutate(Norm_expr = as.numeric(Norm_expr))
    
    ggplot() + 
        geom_sf(data = g.expr$geom_pol,
                aes(fill = g.expr$Norm_expr)) + 
        scale_fill_distiller(palette = "YlGnBu") +
        labs(title = paste0(name, " expression"),
             fill = "Norm. Expr.") + 
        my_theme
    ggsave(file.path(goiExprDir, paste0("goi.expr.map_", name, ".tiff")),
           width = grDevices::dev.size(units = "in")[1],
           height = grDevices::dev.size(units = "in")[2],
           units = "in",
           dpi = 400)
}


