#' Aggregates or summarizes omics data into higher-level functional representations that can be interpreted as functional activity scores or measures.
#'
#' Given an omics matrix and a list of functional molecular sets, this function aggregates or summarizes
#' the omics data into higher-level functional representations such as GO terms gene sets or KEGG metabolic pathways, facilitating the analysis
#' of functional molecular sets that allow reducing dimensionality and providing
#' easier and faster biological interpretations. Coordinated functional activity scores
#' can be as informative as single molecules.
#'
#' Notes:
#' - Different aggregation operators can be used, including summary statistics such as median (default), mean, sd, min, max,
#' dimensionality reduction scores such as pca, mds, pathifier, or nmf, and statistical tests such as ttest, wilcoxon test, kolmogorov test.
#' - The minimum size per molecular set is by default 10 molecules (e.g. genes or metabolites) and can be changed with the parameter minsize.
#' - If "pathifier" is chosen as pooling type, the `aggby_pathifier` function internally generates a log file named 'pathifierlog.txt' during its execution.
#' This log file may contain additional information that could be useful for troubleshooting or advanced analysis.
#' Users typically do not need to interact with this file directly, but it is mentioned here for informational purposes. For more details, this
#' function utilizes the \href{https://bioconductor.org/packages/pathifier/}{Pathifier} package.
#'
#'
#' @param omicsmat A matrix or data frame representing omics data. Rows correspond to molecular identifiers, and columns correspond to samples.
#' @param sets A list of functional sets. Each element in the list should represent a molecular set, and the elements of the set should match the row names of the omics matrix.
#' @param type The type of pooling operator to be applied for each set. Possible values include "mean" (default), "median", "sd", "min", "max", "pca", "mds", "pathifier", "nmf", "ttest", "wilcox", "kolmogorov".
#' @param minsize The minimum size per molecular set (default is 10).
#'
#' @return A matrix-like table with the activity measures for each group or set of molecules, i.e., sxn matrix,
#' for s molecular sets and n samples.
#' @export
#'
#' @examples
#' # Example usage:
#' g <- 10000
#' s <- 20
#' X <- matrix(abs(rnorm(g * s)), nrow = g, dimnames = list(paste0("g", 1:g), paste0("s", 1:s)))
#' pathways <- as.list(sample(10:100, size = 100, replace = TRUE))
#' pathways <- lapply(pathways, function(s, g) paste0("g", sample(1:g, size = s, replace = FALSE)), g)
#' names(pathways) <- paste0("pathway", seq_along(pathways))
#' pathway_activity <- summarize_pathway_level(X, pathways, type = "mean", minsize = 12)

#' @importFrom NMF nmf
#' @importFrom NMF coef
#' @importFrom pathifier quantify_pathways_deregulation
#' @importFrom stats t.test wilcox.test ks.test
#' @importFrom stats prcomp

#' @author Elisa Gomez de Lope
#'
#' @keywords omics aggregation summary functional representation pathway activity-scores

summarize_pathway_level <- function(omicsmat, sets = NULL, type = "mean", minsize = 10) {
    # Parameter checking
    check_parameter(omicsmat, "omicsmat", c("matrix", "data.frame"))
    check_parameter(sets, "sets", "list")

    omicsmat <- as.matrix(omicsmat)
    molnames <- rownames(omicsmat)

    # functional level matrix
    funmat <- matrix(0, nrow = length(sets), ncol = ncol(omicsmat))
    rownames(funmat) <- rep("", nrow(funmat))

    count <- 0
    neg_count <- 0
    message("\n", length(sets), " functional sets read.\n")
    
    if (type == "pathifier") {
        # identify functional sets < minsize
        sets_lengths <- vapply(sets, function(pathway) length(pathway), FUN.VALUE = integer(1))
        short_sets <- names(sets_lengths[sets_lengths < minsize])
        mols2remove <- setdiff(
            unique(unlist(sets[short_sets])),
            unlist(sets[setdiff(names(sets), short_sets)])
        )
        mols <- rownames(omicsmat)[!rownames(omicsmat) %in% mols2remove]
        omicsmat <- omicsmat[match(mols, rownames(omicsmat)), ]
        sets <- sets[setdiff(names(sets), short_sets)]

        funmat <- aggby_pathifier(omicsmat, sets)
        neg_count <- length(short_sets)
        count <- length(sets)
    } else {
        # Loop through functional sets
        for (i in seq_along(sets)) {
            if (i %% 100 == 0) {
                message(sprintf("iteration %i", i))
            }

            gset <- sets[[i]]
            mapid <- match(gset, molnames)
            notna <- which(!is.na(mapid))

            # Check minsize threshold
            if (length(notna) < minsize) {
                neg_count <- neg_count + 1
                next
            }

            # Aggregation based on the specified type
            ifunmat <- omicsmat[mapid[notna], ]
            if (type %in% c("mean", "median", "sd", "min", "max")) {
                rep_vec <- aggby_stat(ifunmat, type)
            } else if (type %in% c("pca", "mds", "nmf")) {
                rep_vec <- aggby_dimred(ifunmat, type)
            } else if (type %in% c("ttest", "wilcox", "kolmogorov")) {
                rep_vec <- aggby_test(ifunmat, type, mapid, notna)
            } else {
              stop(sprintf("Aggregation type %s is not supported.\nPlease check the list of supported aggregation operators.", type))
            }

            count <- count + 1
            funmat[count, ] <- rep_vec
            rownames(funmat)[count] <- names(sets)[i]
        }
    }
    # Print summary & return
    message(paste0(count, " successful functional aggregations over minsize"))
    message(paste0(neg_count, " failed functional aggregations under minsize"))
    if (count >= 1) {
        funmat <- funmat[seq_len(count), ]
        colnames(funmat) <- colnames(omicsmat)
        return(funmat)
    } else {
        message("No functional molecular sets met the criteria.")
        return(NULL)
    }
}


#' @keywords internal
check_parameter <- function(param, name, expected_type = NULL) {
    if (is.null(param)) {
        stop(sprintf("The '%s' parameter must be specified and not NULL.", name))
    }
    if (!is.null(expected_type) && !inherits(param, expected_type)) {
        stop(sprintf("The '%s' parameter must be of type '%s'.", name, expected_type))
    }
}

#' @keywords internal
aggby_stat <- function(X, aggtype) {
    switch(aggtype,
        mean = apply(X, 2, mean),
        median = apply(X, 2, median),
        min = apply(X, 2, min),
        max = apply(X, 2, max),
        sd = apply(X, 2, sd),
    )
}

#' @keywords internal
aggby_dimred <- function(X, aggtype) {
    switch(aggtype,
        mds = as.vector(cmdscale(dist(t(X)), k = 1)),
        pca = {
            rem <- which(apply(X, 1, var) == 0)
            curfunmatfilt <- X
            if (length(rem)) {
                curfunmatfilt <- X[-rem, ]
            }
            if (length(curfunmatfilt)) {
                pca <- stats::prcomp(t(curfunmatfilt), retx = TRUE, scale = TRUE)
                pca$x[, 1]
            } else {
                rep(0, ncol(X))
            }
        },
        nmf = {
            nmf_res <- NMF::nmf(X, rank = 1)
            NMF::coef(nmf_res)
        }
    )
}


#' @keywords internal
aggby_test <- function(X, aggtype, mapid, notna) {
    switch(aggtype,
        ttest = {
            path_outmat <- X[-mapid[notna], ]
            path_ttestres <- vapply(seq_len(ncol(X)), function(x) {
                dat <- t.test(X[, x], path_outmat[, x], alternative = "greater")
                list(dat$stat, dat$p.value)
            },
            FUN.VALUE = list(stat = numeric(1), p.value = numeric(1))
            )
            path_ttest <- as.numeric(path_ttestres[1, ])
            # path_ttestpval = as.numeric(path_ttestres[2,])
            path_ttest
        },
        wilcox = {
            path_outmat <- X[-mapid[notna], ]
            path_wxtestres <- vapply(seq_len(ncol(X)), function(x) {
                dat <- wilcox.test(X[, x], path_outmat[, x], alternative = "greater")
                list(dat$stat, dat$p.value)
            },
            FUN.VALUE = list(stat = numeric(1), p.value = numeric(1))
            )
            path_wxtest <- as.numeric(path_wxtestres[1, ])
            # path_wxtestpval = as.numeric(path_wxtestres[2,])
            path_wxtest
        },
        kolmogorov = {
            path_outmat <- X[-mapid[notna], ]
            path_kstestres <- vapply(seq_len(ncol(X)), function(x) {
                dat <- ks.test(X[, x], path_outmat[, x], alternative = "greater")
                list(dat$stat, dat$p.value)
            },
            FUN.VALUE = list(stat = numeric(1), p.value = numeric(1))
            )
            path_kstest <- as.numeric(path_kstestres[1, ])
            # path_kstestpval = as.numeric(path_kstestres[2,])
            path_kstest
        }
    )
}


#' @keywords internal
aggby_pathifier <- function(X, gs) {
    pathifier_agg <- quantify_pathways_deregulation(
        data = as.matrix(X),
        allgenes = rownames(X),
        syms = gs,
        pathwaynames = names(gs),
        normals = NULL,
        logfile = "pathifierlog.txt",
        attempts = 5,
        min_exp = 0
    ) # = remove effect if min_exp
    pathifier_scores <- data.frame(Reduce(rbind, pathifier_agg$scores))
    colnames(pathifier_scores) <- colnames(X)
    rownames(pathifier_scores) <- names(pathifier_agg$scores)
    return(pathifier_scores)
}
