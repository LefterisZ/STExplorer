
gwpca.prop.var <- function(gwpca.obj, n.comp){
    # apply prop.var() over the n.comp vector
    ptvs.list <- lapply(n.comp, prop.var, gwpca.obj.var = gwpca.obj$var)
    
    # set names to the list objects
    ptvs.list <- setNames(ptvs.list, sprintf("Comps_%02d", n.comp))
    
    # combine into one data frame
    ptvs <- bind_rows(ptvs.list, .id = "column_labels")
    
    # print a summary per column
    print(summary(ptvs))
    
    # make the values percentages and add coordinates
    ptvs <- ptvs * 100
    ptvs <- ptvs %>% 
        mutate(pixel_x = polygons$pixel_x,
               pixel_y = polygons$pixel_y)
    
    return(ptvs)
}
