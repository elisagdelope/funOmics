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


#' Retrieves KEGG pathway gene sets for a specified organism and gene ID type.
#'
#' This function retrieves KEGG pathway gene sets for a specified organism. 
#' It fetches all pathways available for the specified organism from the KEGG database and maps the genes involved in each pathway. 
#' Currently, the function only supports choice of gene identifiers (entrez IDs, gene symbols or Ensembl IDs) for Homo sapiens (organism = "hsa") using the org.Hs.eg.db package.
#' 
#' @param organism The organism abbreviation for which KEGG pathway gene sets are to be retrieved (e.g., "ecj" for E. coli). Default is "hsa" (Homo sapiens).
#' @param geneid_type The type of gene IDs to provide. Default is "entrez"; options are "entrez", "symbol", or "ensembl". 
#'                   This parameter is only used when the organism is "hsa" (Homo sapiens).
#'
#' @return A list where each element represents a KEGG pathway gene set. The names of the list correspond to the pathway names.
#'
#' @export
#'
#' @examples
#' # Retrieve KEGG pathway gene sets for Homo sapiens with entrez IDs (default)
#' hsa_kegg_sets_entrez <- get_kegg_sets()
#' 
#' # Retrieve KEGG molecular sets using gene symbols
#' hsa_kegg_sets_symbol <- get_kegg_sets(geneid_type = "symbol")
#'
#' # Retrieve KEGG molecular sets using Ensembl IDs
#' hsa_kegg_sets_ensembl <- get_kegg_sets(geneid_type = "ensembl")
#'
#' # Retrieve KEGG pathway gene sets for another organism (e.g., Escherichia coli)
#' ecoli_kegg_sets <- get_kegg_sets(organism = "ecj")
#' 
#' @importFrom dplyr select left_join %>% tibble
#' @importFrom AnnotationDbi mapIds
#' @import org.Hs.eg.db
#' @importFrom KEGGREST keggLink keggList
#' @importFrom stringr str_match str_extract
#'
#' @seealso \code{\link{summarize_pathway_level}}
#' @seealso \code{\link{keggLink}}, \code{\link{keggList}}
#' @seealso \code{\link{mapIds}}
#' 
get_kegg_sets <- function(organism="hsa", geneid_type="entrez") {
  # checks
  check_parameter(organism, "organism", "character")
  check_parameter(geneid_type, "geneid_type", "character")
  if (!(geneid_type %in% c("entrez" , "symbol", "ensembl"))) {
    stop("Invalid gene ID type. Please use one of: ", paste(c("entrez" , "symbol", "ensembl"), collapse = ", "))
  }

  # get all pathways and their entrez gene ids
  path_entrez <- tryCatch({ keggLink("pathway", organism) %>% 
    tibble(pathway = gsub("path:", "", .), geneID = sub(paste0(organism, ":"), "", names(.))) %>%
    dplyr::select(-.)
  }, error = function(e) {
    stop("Invalid organism abbreviation or unsupported organism: ", organism)
  })
  
  # get pathway names
  kegg_pathways <- keggList("pathway", organism) %>% 
    tibble(pathway = names(.), description = .) %>%
    mutate(description = unname(description))
  org_substring <- str_match(kegg_pathways$description, ".* - (.*)")[, 2]
  if (length(unique(org_substring)) == 1) {
    kegg_pathways$description <- str_extract(kegg_pathways$description, "^.*(?= - )")
  }
  if (organism=="hsa") {
    # get gene symbols and ensembl ids using the gene ids (entrez IDs) retrieved from kegg
    path_entrez <- path_entrez %>%
      mutate(
        symbol = suppressMessages(unname(mapIds(org.Hs.eg.db, geneID, "SYMBOL", "ENTREZID"))),
        ensembl = suppressMessages(unname(mapIds(org.Hs.eg.db, geneID, "ENSEMBL", "ENTREZID")))
      ) 
  }
  
  # merge
  KEGG_pathways <- left_join(kegg_pathways, path_entrez, by = "pathway")
  
  # Split by gene id type if hsa
  if (organism=="hsa") {
    kegg_sets <- switch(geneid_type,
                        "entrez" = split(KEGG_pathways$geneID, KEGG_pathways$description),
                        "symbol" = split(KEGG_pathways$symbol, KEGG_pathways$description),
                        "ensembl" = split(KEGG_pathways$ensembl, KEGG_pathways$description),
                        stop("Invalid geneid_type. Must be 'entrez', 'symbol', or 'ensembl'.")
    )
  }
  else {
    kegg_sets <- split(KEGG_pathways$geneID, KEGG_pathways$description)
  }
  return(kegg_sets)
}

