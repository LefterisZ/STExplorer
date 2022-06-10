#----------------------------------------------------#
## CREATE/SET FOLDERS/PATHS ----
#----------------------------------------------------#

## Set file paths
projDir <- file.path(getwd(), "data")
inputDir <- file.path(projDir, "spaceranger_outs/")
outputDir <- file.path(projDir, "rObjects/")
csvDir <- file.path(projDir, "csvFiles/")
graphDir <- file.path(projDir, "graphics_out/")

## Check if inputDir/ outputDir/ csvDir exist and create them if not.
dirs <- c(inputDir, outputDir, csvDir, graphDir)

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
