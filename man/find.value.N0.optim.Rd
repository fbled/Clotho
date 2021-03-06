% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effective population size.R
\name{find.value.N0.optim}
\alias{find.value.N0.optim}
\title{Base for N0 optimization function}
\usage{
find.value.N0.optim(
  N0,
  threshold.N,
  surv.p,
  sex.ratio,
  maturity.age.F,
  maturity.age.M,
  max.age.F = NULL,
  max.age.M = NULL
)
}
\arguments{
\item{\code{N0}}{numeric, estimate for initial population size for age-class 0}

\item{\code{threshold.N}}{numeric, adult population threshold target}

\item{\code{surv.p}}{A data frame, with 2 columns "\code{surv.F}" and "\code{surv.M}" giving for each age-sex class the corresponding survival probability}

\item{\code{sex.ratio}}{numeric, sex ratio of females in \code{N0}}

\item{\code{maturity.age.F}, \code{maturity.age.M}}{numeric, sexual maturity age for females and males}

\item{\code{max.age.F}, \code{max.age.M}}{numeric, maximum age for females and males. If not included, it is inferred from the survival probability data.}
}
\value{
a numeric, absolute difference between threshold target and simulated total adult population.
}
\description{
Provide metric to optimize initial population size for age-class 0 based on adult population threshold target
}
\examples{
N0=182723
surv.p=data.frame(  surv.F=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
                              0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0) ,
                    surv.M=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
                              0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0)
                  )
sex.ratio=0.444
maturity.age.F=4
maturity.age.M=4
N.res=find.value.N0.optim(N0,surv.p,sex.ratio,maturity.age.F,maturity.age.M)
N.res$Ntot
}
