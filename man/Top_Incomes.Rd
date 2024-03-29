% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hello.R
\name{Top_Incomes}
\alias{Top_Incomes}
\title{Top Incomea plot}
\usage{
Top_Incomes(
  data,
  p = 0.01,
  thr = seq(0.85, 0.999, by = 0.001),
  tail_index = TRUE,
  ...
)
}
\arguments{
\item{data}{dataframe obtained from \code{\link{tidy_income}} function}

\item{p}{probability level (default 0.01)}

\item{thr}{numeric vector of probability levels to model a Pareto distribution (default is \code{seq(.85,.999,by=.001)} from 0.85 up to 0.999)}

\item{tail_index}{logical to plot the tail index (default \code{TRUE})}
}
\value{
one or two graphs (depending on \code{tail==TRUE})
}
\description{
Top Incomea plot
}
\examples{
url_1 <- "https://github.com/freakonometrics/TopIncome/raw/master/data_csv/dataframe_yw_1.csv"
df <- read.table(url_1,sep=";",header=TRUE)
data_1  <-  tidy_income(income = df$y, weights = df$w)
\dontrun{Top_Incomes(data_1)}
}
