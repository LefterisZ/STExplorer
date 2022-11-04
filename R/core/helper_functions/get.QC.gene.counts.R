# folders <- list.files(inputDir, 
#                       recursive = TRUE, 
#                       pattern = "Results", 
#                       include.dirs = TRUE)

get.QC.gene.counts <- function(inDir, folders){
    metrics <- list()
    for(f in folders){
        sampleDir <- f
        countsDir <- file.path(inDir, sampleDir, "filtered_feature_bc_matrix")
        message("Working on folder: ", f)
        message("Importing count data ...")
        inputD <- readSpacerangerD(countsDir)
        message("Filtering count data ...")
        inputD_filt <- inputD[!rowSums(inputD) < 10,]
        
        message("Calculating metrics ...")
        Total = sum(inputD_filt)
        perSpot = sum(inputD_filt)/dim(inputD_filt)[2]
        perGenePerSpot = sum(inputD_filt)/(dim(inputD_filt)[1]*dim(inputD_filt)[2])
        Spots.No = dim(inputD_filt)[2]
        filt.Genes.No = dim(inputD_filt)[1]
        
        df <- data.frame("Total" = Total,
                          "perSpot" = perSpot,
                          "perGenePerSpot" = perGenePerSpot,
                          "Spots.No" = Spots.No,
                          "filt.Genes.No" = filt.Genes.No)
        message("Appending ...")
        metrics <- list.append(metrics, df)
    }
    
    names <- gsub("^.*s_", "", folders) %>%
        gsub("_Results", "", .)
    names(metrics) <- names
    
    out <- bind_rows(metrics)
    rownames(out) <- names(metrics)
    message("READY!!")
    return(out)

}
