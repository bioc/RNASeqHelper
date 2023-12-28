
test_that("rnaseqhelper", {
    n <- 1000
    tab <- matrix(as.integer(rnorm(n**2, 1000, 500)), ncol = n)
    tab <- tab - min(tab)
    rownames(tab) <- paste0("G", 1:n)
    colnames(tab) <- paste0("S", 1:n)
    phenotype_data <- data.frame(
        sample = paste0("S", 1:n),
        condition = c(rep("red", n / 2), rep("green", n / 2)),
        time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
    )
    rnaseqhelper(tab, phenotype_data,
                out_dir = tempdir(), "green", "red",
                heat_params = list(
                    score_thresh = c(0.2, 0.5),
                    kmeans_list = 2
                ))
})

test_that("run_deseq", {
    n <- 1000
    tab <- matrix(as.integer(rnorm(n**2, 1000, 500)), ncol = n)
    tab <- tab - min(tab)
    rownames(tab) <- paste0("G", 1:n)
    colnames(tab) <- paste0("S", 1:n)
    phenotype_data <- data.frame(
        sample = paste0("S", 1:n),
        condition = c(rep("red", n / 2), rep("green", n / 2)),
        time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
    )
    keep_genes <- paste0("G", 1:n)
    res <- run_deseq(tab, keep_genes, phenotype_data, tempdir())
})

test_that("top_n_genes", {
    dsqres <- data.frame(mLog10Padj = 1:10 / 10, gene = paste0("G", 1:10))
    RNASeqHelper:::top_n_genes(dsqres, 3, out_dir = tempdir(), prefix = "test")
})

test_that("heatmap_with_geneplots", {
    n <- 100
    norm_counts <- matrix(rnorm(n**2, mean = 5), nrow = n)
    pheno <- data.frame(
        sample = paste0("S", 1:n),
        condition = c(rep("red", n / 2), rep("green", n / 2)),
        time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
    )
    rownames(norm_counts) <- paste0("G", 1:n)
    colnames(norm_counts) <- paste0("S", 1:n)
    genes <- list(
        highlight = paste0("G", 1:10),
        interest = paste0("G", c(2, 3)),
        score_thresh = c(0.2, 0.3)
    )
    res <- heatmap_with_geneplots(norm_counts, pheno, 2,
                                genes,
                                out_dir = tempdir())
})

test_that("calculate_cluster_corr", {
    n <- 100
    scale_mat <- matrix(rnorm(n**2, mean = 0), nrow = n)
    rownames(scale_mat) <- paste0("G", 1:n)
    colnames(scale_mat) <- paste0("S", 1:n)
    clust_assign <- data.frame(
        gene = rownames(scale_mat),
        cluster = c(rep(1, n / 50), rep(2, n / 50)),
        value = rnorm(n, 100, 10)
    )
    ca <- calculate_cluster_corr(clust_assign, scale_mat, tempdir(), "red")
})

test_that("do_gene_plots", {
    n <- 100
    norm_counts <- matrix(rnorm(n**2, mean = 5), nrow = n)
    rownames(norm_counts) <- paste0("G", 1:n)
    colnames(norm_counts) <- paste0("S", 1:n)
    scale_mat <- matrix(rnorm(n**2, mean = 0), nrow = n)
    rownames(scale_mat) <- paste0("G", 1:n)
    colnames(scale_mat) <- paste0("S", 1:n)
    pheno <- data.frame(
        sample = paste0("S", 1:n),
        condition = c(rep("red", n / 2), rep("green", n / 2)),
        time = as.integer(rnorm(n, 2, 0.5) + 1) * 5
    )
    gene_cluster_scores <- data.frame(
        gene = rownames(scale_mat),
        cluster = c(rep(1, n / 50), rep(2, n / 50)),
        score = rnorm(n, 0.5, 0.25)
    )
    score_thresh <- c(0.2, 0.3)
    genes_of_interest <- paste0("G", 5:50)
    res <- do_gene_plots(norm_counts, scale_mat, pheno,
                        gene_cluster_scores, score_thresh,
                        genes_of_interest,
                        out_dir = tempdir())
})

test_that("single_kmeans_heatmap", {
    n <- 100
    scale_mat <- matrix(rnorm(n**2, mean = 5), nrow = n)
    rownames(scale_mat) <- paste0("G", 1:n)
    colnames(scale_mat) <- paste0("S", 1:n)
    top_genes <- paste0("G", n - 50:n - 20)
    res <- single_kmeans_heatmap(scale_mat, 2, top_genes, "test", "test",
                                out_dir = tempdir(), "test", 7, 7)
})

test_that("high_quality_genes", {
    sam_mat <- matrix(1:100, nrow = 10,
                    dimnames = list(
                        paste0("G", 1:10),
                        paste0("C", 1:10)
                    ))
    keep <- high_quality_genes(sam_mat, 10, 7, tempdir())
    names(keep) == paste0("G", 7:10)
})

test_that("do_volcanos", {
    n <- 100
    dsqres <- data.frame(
        gene = paste0("G", 1:n),
        data = rnorm(n, 5, 2),
        log2FoldChange = 1e-5,
        mLog10Padj = rnorm(n, 1000, 2),
        padj = 1/rnorm(n, 1000, 2)
    )
    top_genes_tohighlight <- paste0("G", (n-50):(n-20))
    res <- RNASeqHelper:::do_volcanos(dsqres, top_genes_tohighlight,
                                    "my plot", tempdir())
})

test_that("volcano_plot", {
    n <- 100
    dsqres <- data.frame(gene = paste0("G", 1:n),
                        data = rnorm(n, 5, 2),
                        log2FoldChange = 1e-5,
                        padj = 1 / rnorm(n, 5, 2))
    degenes <- paste0("G", (n-50):(n-20))
    res <- volcano_plot(dsqres, degenes, "my plot")
})

test_that("gene_clusters_by_score", {
    n <- 100
    tab <- data.frame(
        gene = paste0("G", 1:n),
        cluster = as.integer(rnorm(n, 5, 2) + 1),
        condition = c(rep("red", n/2), rep("green", n/2)),
        value = rnorm(n, 10, 2) + 1,
        time = as.integer(rnorm(n, 2, 0.5) + 1) * 5,
        score = as.integer(rnorm(n, 2, 0.5) + 1) / 4
    )
    plot <- gene_clusters_by_score(tab, c(0.2, 0.7),
                                    out_dir = tempdir())
})

test_that("generate_labelling_table", {
    tabn <- data.frame(
        cluster = c(1, 1, 2, 2),
        gene = paste0("G", 1:4),
        value = 1:4
    )
    res <- generate_labelling_table(tabn)
})

test_that("num_clusters", {
    expect_equal(
        num_clusters(data.frame(
            cluster = c(1, 1, 2, 3, 2)
        )),
        3
    )
})