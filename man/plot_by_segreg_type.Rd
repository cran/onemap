% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_raw_data.R
\name{plot_by_segreg_type}
\alias{plot_by_segreg_type}
\title{Draw a graphic showing the number of markers of each segregation pattern.}
\usage{
plot_by_segreg_type(x, subcateg = TRUE)
}
\arguments{
\item{x}{an object of class \code{onemap}}

\item{subcateg}{a TRUE/FALSE option to indicate if results will be plotted showing
all possible categories (only for outcrossing populations)}
}
\value{
a ggplot graphic
}
\description{
The function receives an object of class \code{onemap}.
For outcrossing populations, it can show detailed information (all 18 possible categories),
or a simplified version.
}
\examples{
\donttest{
data(onemap_example_out) #Outcrossing data
plot_by_segreg_type(onemap_example_out)
plot_by_segreg_type(onemap_example_out, subcateg=FALSE)

data(onemap_example_bc)
plot_by_segreg_type(onemap_example_bc)

data(mapmaker_example_f2)
plot_by_segreg_type(mapmaker_example_f2)
}

}
