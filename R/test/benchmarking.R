#' This script contains benchmarking code for different functions.


benchmark(
    "apply_on_the_fly" <- {
        apply(bin_0,
              1,
              spot_neighbours,
              bin.1 = bin_1)
    },
    "apply_test_apply" <- {
        cbind(test_tissue_positions,
              "new_bins" = apply(test_tissue_positions, 1, test_apply))
    },
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


benchmark(
    "nbs[[]]" <- {
        nbs[["Section"]]
    },
    "nbs[]" <- {
        nbs["Section"]
    },
    replications = 1000,
    columns = c(
        "test",
        "replications",
        "elapsed",
        "relative",
        "user.self",
        "sys.self"
    )
)


benchmark(
    "test_perimeter" <- {
        test_add_perimeter(test_tissue_positions)
    },
    "perimeter" <- {
        add_perimeter(test_tissue_positions)
    },
    replications = 10,
    columns = c(
        "test",
        "replications",
        "elapsed",
        "relative",
        "user.self",
        "sys.self"
    )
)

