gwpca.leading.G.single <- function(input, pc.no, pc.name){
    #get local loadings for pc.no
    local.loadings <- round(input$loadings[,,pc.no], 4)
    
    #find the leading item (colname) on each location (rows)
    lead.item <- colnames(local.loadings)[max.col(abs(local.loadings))]
    
    #print how many items are leading in this PC
    print(paste0("The leading items in ", pc.name, " are: ", unique(lead.item)))
    
    lead.item <- data.frame(lead.item)
    colnames(lead.item) <- pc.name
    
    return(lead.item)
}
