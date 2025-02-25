% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_utils.R
\name{dot-int_spotHex_subset}
\alias{dot-int_spotHex_subset}
\alias{.int_spotHex_subset}
\title{Internal: get info for spots on the edge}
\usage{
.int_spotHex_subset(minMaxList, .cData, .data)
}
\arguments{
\item{minMaxList}{A list containing the minimum and maximum x and y
coordinate values.}

\item{.cData}{A data frame containing the barcode and array coordinate
information.}

\item{.data}{A data frame containing the spot data.}
}
\value{
A list with the following elements:
\describe{
\item{bcd_Xmax}{Data frame with the spots corresponding to the maximum
x-coordinate value.}
\item{bcd_Xmin}{Data frame with the spots corresponding to the minimum
x-coordinate value.}
\item{bcd_Ymax}{Data frame with the spots corresponding to the maximum
y-coordinate value.}
\item{bcd_Ymin}{Data frame with the spots corresponding to the minimum
y-coordinate value.}
\item{bcd_Xmax1}{Data frame with the spots corresponding to the x-coordinate
value adjacent to the maximum x-coordinate.}
\item{bcd_Xmin1}{Data frame with the spots corresponding to the x-coordinate
value adjacent to the minimum x-coordinate.}
}
}
\description{
IMPORTANT: this function is not exported to be used alone.
This function fetches data and coordinates for the required locations
based on the provided minimum and maximum values.
}
\details{
This function is used internally to fetch data and coordinates for specific
locations based on the provided minimum and maximum values. It takes the
barcode and array coordinate information from the \code{.cData} data frame
and the spot data from the \code{.data} data frame. It returns a list of data
frames, where each data frame corresponds to a specific location based on the
minimum and maximum values.
}
\author{
Eleftherios (Lefteris) Zormpas
}
