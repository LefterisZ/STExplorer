library(rhdf5)
filtered_hs <- H5Fopen("./data/rawData/H37_Steatotic/GSM5764427_filtered_feature_bc_matrix_JBO019.h5")
h5ls(filtered_hs)


## Export barcodes.tsv
write_tsv(as.data.frame(filtered_hs$matrix$barcodes), 
          file = "./data/rawData/H37_Steatotic/barcodes.tsv",
          col_names = FALSE)

## Export features.tsv
write_tsv(data.frame(id = filtered_hs$matrix$features$id,
                     name = filtered_hs$matrix$features$name,
                     feature_type = filtered_hs$matrix$features$feature_type), 
          file = "./data/rawData/H37_Steatotic/features.tsv",
          col_names = FALSE)

## Export matrix.mtx
Matrix::writeMM(obj = Matrix::sparseMatrix(i = filtered_hs$matrix$indices,          # Zero-based row index of corresponding element in data
                                           p = filtered_hs$matrix$indptr,           # IndexPointer: Zero-based index into data / indices of the start of each column, that is the data corresponding to each barcode sequence
                                           x = as.numeric(filtered_hs$matrix$data), # Nonzero UMI counts in column-major order
                                           dims = filtered_hs$matrix$shape,         # Tuple of (# rows, # columns) indicating the matrix dimensions
                                           repr = "C",
                                           index1 = FALSE),
                file = "./data/rawData/H37_Steatotic/matrix.mtx")
