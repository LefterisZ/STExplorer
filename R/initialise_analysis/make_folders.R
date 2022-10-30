#----------------------------------------------------#
## CREATE/SET FOLDERS/PATHS ----
#----------------------------------------------------#

## Set file paths
projDir <- file.path(getwd(), "data")
inputDir <- file.path(projDir, "spaceranger_outs")
outputDir <- file.path(projDir, "rObjects")
csvDir <- file.path(projDir, "csvFiles")
graphDir <- file.path(projDir, "graphics_out")
countsDir <- file.path(inputDir, "Olfactory_Bulb/Olfactory_Bulb_A1_Results/filtered_feature_bc_matrix")

## Check if inputDir/ outputDir/ csvDir exist and create them if not.
dirs <- c(inputDir, outputDir, csvDir, graphDir, countsDir)

dirCheck <- lapply(
    dirs,
    FUN = function(x){
        if(!dir.exists(x)) {
            dir.create(x, recursive = TRUE)
            print(paste0("Folder: ", x, " CREATED!"))
        } else {
            print(paste0("Folder: ", x, " EXISTS!"))
        }
    }
)
