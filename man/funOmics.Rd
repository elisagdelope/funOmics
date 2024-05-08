
\name{funOmics}

\title{Aggregating Omics Data into Higher-Level Functional Representations}

\description{
The `funOmics` package aggregates or summarizes omics data into higher-level functional representations such as GO terms, gene sets, or KEGG metabolic pathways. The aggregated data matrix represents functional activity scores that facilitate the analysis of functional molecular sets while allowing dimensionality reduction and providing easier and faster biological interpretations. Coordinated functional activity scores can be as informative as single molecules!
}

\usage{
get_kegg_sets(organism = "hsa", geneid_type = "entrez")
summarize_pathway_level(omicsmat, sets = NULL, type = "mean", minsize = 10)
}

\details{
Pre-release and latest development version on https://github.com/elisagdelope/funomics
}

\author{
Elisa Gomez de Lope, Enrico Glaab
}
