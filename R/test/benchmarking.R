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


benchmark(
    "in" <- {## Rfast::transpose needs a SQUARE numeric matrix
        if(!dim(dists.X)[1] == dim(dists.X)[2]){
            stop("The distance matrix must be a square. Check your input!")
        }
    
        ## Find closest neighbours
        ## a. get their distances
        distances <- Rfast::transpose(apply(dists.X, 1,
                                            function(x, part){
                                                Rfast::Sort(x, partial = part)[1:part]
                                            },
                                            part = k))
        ## b. get their indexes
        neighbour_idx <- t(apply(dists.X, 1,
                                 function(x, partial) {
                                     colnames(dists.X)[Rfast::Match(Rfast::Sort(x, partial = partial)[1:k], x)]
                                 },
                                 partial = k))
        
        ## Prepare output
        out = list(neighbour_idx, distances)
        names(out) <- c("indexes", "distances")
    },
    "out" <- {
        ## Rfast::transpose needs a SQUARE numeric matrix
        if (!dim(dists.X)[1] == dim(dists.X)[2]) {
            stop("The distance matrix must be a square. Check your input!")
        }
        
        ## Find closest neighbours
        ## a. get their distances
        distances <- Rfast::transpose(apply(dists.X, 1,
                                            .Sort.part,
                                            part = k))
        ## b. get their indexes
        neighbour_idx <- t(apply(dists.X, 1,
                                 function(x, partial) {
                                     colnames(dists.X)[Rfast::Match(Rfast::Sort(x, partial = partial)[1:k], x)]
                                 },
                                 partial = k))
        
        ## Prepare output
        out = list(neighbour_idx, distances)
        names(out) <- c("indexes", "distances")
    }, 
    replications = 100)

profvis::profvis({
    g <- .get.gw.counts(focus = focus, obs = obs, wdmat = wdmat) %>% # weight gene expression
        .dist.wrap(gwCounts = ., method = method, p = p) %>% # Calculate distances
        .get.knn(dists.X = ., k = k) %>% # Get nearest neighbours
        .get.gwGraph(kList = ., combine = FALSE) # Get the graphs
    
    # Store it in a list
    updt(x = graphs_list, f) <- g
    
    # free up some memory
    rm(g)
})

benchmark("future" <- {
    graphs_list <- future_lapply(focus.n[1, drop = FALSE], 
                                 get.gwGraphs, 
                                 obs = obs,
                                 wdmat = wdmat,
                                 method = method,
                                 p = p)
},
"foreach" <- {
    
}, replications = 10)

benchmark("multicore" <- {get.gwGraph.list(.focus.n = focus.n[1:9],
                                          .obs = obs,
                                          .wdmat = wdmat,
                                          .method = "euclidean",
                                          .p = 0,
                                          .k = 7,
                                          .strategy = "multicore",
                                          .workers = 9,
                                          verbose = FALSE,
                                          .global = TRUE, 
                                          .handlers = c("progress", "beepr"))
},
"sequential" <-  {get.gwGraph.list(.focus.n = focus.n[1:9],
                                   .obs = obs,
                                   .wdmat = wdmat,
                                   .method = "euclidean",
                                   .p = 0,
                                   .k = 7,
                                   .strategy = "sequential",
                                   .workers = 9,
                                   verbose = TRUE,
                                   .global = TRUE, 
                                   .handlers = c("progress", "beepr"))
    
},
replications = 1,
columns = c(
    "test",
    "replications",
    "elapsed",
    "relative",
    "user.self",
    "sys.self"
))


profvis::profvis({
    graphs_list <- get.gwGraph.list(.focus.n = focus.n[1:9],
                                       .obs = obs,
                                       .wdmat = wdmat,
                                       .method = "euclidean",
                                       .p = 0,
                                       .k = 7,
                                       .strategy = "multicore",
                                       .workers = 9,
                                       verbose = FALSE,
                                       .global = TRUE, 
                                       .handlers = c("progress", "beepr"))
})

profvis::profvis({
    get.gwGraph.table(..focus = focus, 
                      ..obs = obs,
                      ..wdmat = wdmat,
                      ..method = method,
                      ..p = p,
                      ..k = k)
})

profvis::profvis({
    get.gwGraph.table(..focus = "1",
                      ..obs = obs,
                      ..wdmat = wdmat,
                      ..method = method,
                      ..k = k,
                      threads = 9)
})

system.time(
    get.gwGraph.table(..focus = "1",
                      ..obs = obs,
                      ..wdmat = wdmat,
                      ..method = method,
                      ..k = k)
)

s <- Sys.time()
graphs_list <- get.gwGraph.list(.focus.n = focus.n,
                                .obs = obs,
                                .wdmat = wdmat,
                                .method = "euclidean",
                                .k = 7,
                                verbose = TRUE)
t <- Sys.time()


benchmark(
    "Dist" <- {
        Rfast::Dist(x = obs.W, method = "euclidean")
        },
    "parDist" <- {
        parallelDist::parDist(x = obs.W, method = "euclidean", threads = 1)
    },
    replications = 10,
    columns = c(
        "test",
        "replications",
        "elapsed",
        "relative",
        "user.self",
        "sys.self")
)


profvis::profvis({
    ## Take out a matrix from the list of graphs
    wGraph_s <- as.data.frame(graph.list["1"],
                              col.names = dimnames(graph.list["1"])[[2]])
    
    ## Transform counts and wDist columns to numerics
    wGraph_s$count <- as.numeric(wGraph_s$count)
    wGraph_s$wDist <- as.numeric(wGraph_s$wDist)
    
    ## Make a matching vector for the flags in question
    x <- match(wGraph_s[,"flag"], eges_mtx$flag, nomatch = 0) 
    
    ## Subset and replace the count
    eges_mtx[x, "count"] <- eges_mtx[x, "count"] + wGraph_s$count
    
    ## Subset and replace the weights
    eges_mtx[x, "wDistSum"] <- eges_mtx[x, "wDistSum"] + wGraph_s$wDist
})
