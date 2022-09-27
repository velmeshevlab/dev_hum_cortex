function (cds, neighbor_graph = c("knn", "principal_graph"),
    reduction_method = "UMAP", k = 25, method = c("Moran_I"),
    alternative = "greater", expression_family = "quasipoisson",
    cores = 1, verbose = FALSE) 
{
    neighbor_graph <- match.arg(neighbor_graph)
    lw <- monocle3:::calculateLW(cds, k = k, verbose = verbose, neighbor_graph = neighbor_graph, 
        reduction_method = reduction_method)
    if (verbose) {
        message("Performing Moran's I test: ...")
    }
    exprs_mat <- SingleCellExperiment::counts(cds)[, attr(lw, 
        "region.id"), drop = FALSE]
    sz <- size_factors(cds)[attr(lw, "region.id")]
    wc <- spdep::spweights.constants(lw, zero.policy = TRUE, 
        adjust.n = TRUE)
    test_res <- pbmcapply::pbmclapply(row.names(exprs_mat), FUN = function(x, 
        sz, alternative, method, expression_family) {
        exprs_val <- exprs_mat[x, ]
        if (expression_family %in% c("uninormal", "binomialff")) {
            exprs_val <- exprs_val
        }
        else {
            exprs_val <- log10(exprs_val/sz + 0.1)
        }
        df = cbind(as.data.frame(exprs_val), colData(cds)$sex, colData(cds)$region_broad)
        colnames(df) <- c("exp", "sex", "region")
        test_res <- tryCatch({
            if (method == "Moran_I") {
                mt <- suppressWarnings(monocle3:::my.moran.test(df, lw, wc, alternative = alternative))
                data.frame(status = "OK", p_value = mt$p.value, 
                  morans_test_statistic = mt$statistic, morans_I = mt$estimate[["Moran I statistic"]])
            }
            else if (method == "Geary_C") {
                gt <- suppressWarnings(my.geary.test(exprs_val, 
                  lw, wc, alternative = alternative))
                data.frame(status = "OK", p_value = gt$p.value, 
                  geary_test_statistic = gt$statistic, geary_C = gt$estimate[["Geary C statistic"]])
            }
        }, error = function(e) {
            data.frame(status = "FAIL", p_value = NA, morans_test_statistic = NA, 
                morans_I = NA)
        })
    }, sz = sz, alternative = alternative, method = method, expression_family = expression_family, 
        mc.cores = cores, ignore.interactive = TRUE)
    if (verbose) {
        message("returning results: ...")
    }
    test_res <- do.call(rbind.data.frame, test_res)
    row.names(test_res) <- row.names(cds)
    test_res <- merge(test_res, rowData(cds), by = "row.names")
    row.names(test_res) <- test_res[, 1]
    test_res[, 1] <- NULL
    test_res$q_value <- 1
    test_res$q_value[which(test_res$status == "OK")] <- stats::p.adjust(subset(test_res, 
        status == "OK")[, "p_value"], method = "BH")
    test_res$status = as.character(test_res$status)
    test_res[row.names(cds), ]
}
