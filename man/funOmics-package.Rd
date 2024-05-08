\name{funOmics-package}
\alias{funOmics-package}
\alias{funOmics}
\docType{package}
\title{
Aggregating Omics Data into Higher-Level Functional Representations
}
\description{
The `funOmics` package aggregates or summarizes omics data into higher-level functional representations such as GO terms, gene sets, or KEGG metabolic pathways. The aggregated data matrix represents functional activity scores that facilitate the analysis of functional molecular sets while allowing dimensionality reduction and providing easier and faster biological interpretations. Coordinated functional activity scores can be as informative as single molecules!
}
\details{
%%  ~~ An overview of how to use the package, ~~

The package offers functionalities for:
\itemize{
  \item {Data aggregation into functional representations}
  \item {Dimensionality reduction of omics data sets}
  \item {Analysis of coordinated functional activity scores}
}

Use function `get_kegg_sets` to get KEGG pathway sets for a given organism and geneid type:
  
  get_kegg_sets(organism = "hsa", geneid_type = "entrez")

Use function `summarize_pathway_level` to summarize omics data into higher-level functional representations that can be interpreted as functional activity scores or measures:
  
  summarize_pathway_level(omicsmat, sets = NULL, type = "mean", minsize = 10)
}
\author{
Elisa Gomez de Lope, Enrico Glaab

Maintainer: Elisa Gomez de Lope (elisa.gomezdelope@uni.lu)
}
\references{
%%  ~~ Literature or other references for background information ~~
}
\keyword{package}

%% Uncomment below to imitate parts of library(help = funOmics)
%\section{The \file{DESCRIPTION} File}{\packageDESCRIPTION{funOmics}}
%\section{Documentation Index}{\packageIndices{funOmics}}

\seealso{
%%\code{\link[<funOmics>:<funOmics>-package]{<funOmics>}}
\code{\link{funOmics::summarize_pathway_level}}
\code{\link{funOmics::get_kegg_sets}}

}
\examples{
%%  ~~ Optional simple examples of the most important functions ~~
}
