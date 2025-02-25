% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgwc_plot.R
\name{.int_markers}
\alias{.int_markers}
\title{Process cell type Markers}
\usage{
.int_markers(markers)
}
\arguments{
\item{markers}{A data frame containing information about cell type
markers, including gene names and types.}
}
\value{
A processed data frame with validated column names, arranged by
marker types, and duplicates removed.
}
\description{
This function processes cell type markers, ensuring required column
names, arranging by marker types, and removing duplicates.
}
\details{
The function validates the column names of the input data frame,
ensuring the presence of "gene.name," "ensg.ID," "Type," and "Subtype"
columns. It arranges the markers by types and removes duplicate entries
based on ENSG.IDs, resulting in a processed data frame ready for further
analysis.
}
\seealso{
\code{\link{.int_validateColumnNames}},
\code{\link{.int_arrangeByType}}, \code{\link{.int_removeDuplicates}}
}
\author{
Eleftherios (Lefteris) Zormpas
}
\keyword{cell}
\keyword{duplicates}
\keyword{markers,}
\keyword{type,}
\keyword{validation,}
