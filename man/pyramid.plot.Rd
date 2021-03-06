% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics.R
\name{pyramid.plot}
\alias{pyramid.plot}
\title{Pyramid plot}
\usage{
pyramid.plot(N.res)
}
\arguments{
\item{\code{N.res}}{output from a "\code{eff.size.fn}" computation}
}
\value{
a population pyramid plot, with \code{Ntot} (total number of adults), \code{N0} (estimated initial population size for age-class 0)
}
\description{
Plots population pyramid plot from "\code{eff.size.fn}" output
}
\examples{
surv.p=data.frame(  surv.F=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
                              0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0) ,
                    surv.M=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
                              0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
                              0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0)
                  )

maturity.age.F=4
maturity.age.M=4
sex.ratio=0.444

threshold.N.adult=199595

N.res=eff.size.fn(threshold.N=threshold.N.adult,surv.p,sex.ratio,maturity.age.F,maturity.age.M,method="BFGS")
pyramid.plot(N.res)
}
