#' This script contains benchmarking code for different functions.


benchmark(
    "apply_on_the_fly" <- {
        apply(bin_0,
              1,
              spot_neighbours,
              bin.0 = bin_0,
              bin.1 = bin_1)
    },
    # "apply_test_apply" <- {
    #     cbind(test_tissue_positions,
    #           "new_bins" = apply(test_tissue_positions, 1, test_apply))
    # },
    replications = 1,
    columns = c(
        "test",
        "replications",
        "elapsed",
        "relative",
        "user.self",
        "sys.self"
    )
)
