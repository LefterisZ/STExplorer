samples <- names(msfe@sfe_data)

for (s in samples) {
  result <- fgwc_nmfFactorNumber(m_sfe = msfe,
                                 sample_id = s,
                                 assay = "logcounts",
                                 top_hvgs = top_hvgs[[s]],
                                 k.range = seq(2, 10, 1),
                                 n.cores = 1,
                                 do.plot = FALSE,
                                 seed = 1,
                                 loss = "mse",
                                 max.iter = 250)

  print(plotFGWC_factorSelection(result))
}


# ---------------------------------------------------------------------------- #


