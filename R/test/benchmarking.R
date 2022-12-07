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






benchmark(
    "direct_call" <- {
        d <- sapply(focus.n, function(X, method, p){
            dist(data[,,X], method, p) %>% as.matrix()
        }, method = "euclidian", p = 2, simplify = "array")
    },
    "indirect_call" <- {
        get.dist.array(obs.W, focus.n = 1:2, method = "euclidean", p = 2)
    },
    replications = 100,
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
    "apply" <- {
        edge.comb <- expand.grid(nodes, nodes, stringsAsFactors = FALSE) %>% # get all possible edges
            mutate(temp = apply(., 1, sort.N.paste)) %>%
            .[!duplicated(.$temp),]
    },
    "mclapply" <- {
        edge.comb <- expand.grid(nodes, nodes, stringsAsFactors = FALSE) %>% # get all possible edges
            t() %>%
            data.frame()
        temp = mclapply(edge.comb, sort.N.paste)
        edge.comb <- edge.comb %>% 
            mutate(temp = temp) %>%
            .[!duplicated(.$temp),]
    },
    replications = 2,
    columns = c(
        "test",
        "replications",
        "elapsed",
        "relative",
        "user.self",
        "sys.self"
    )
)
