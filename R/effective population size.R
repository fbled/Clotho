#' Testing vector for contiguous survival data starting from the first element
#'
#' @param \code{vec} vector of non-harvest mortality rate of length Y
#'
#' @return \code{TRUE} if the vector is composed either of non missing values, or start with a contiguous set of non missing values followed by NAs.
#'
#' @examples
#' test.fn(c(4,3,6,1,7,8))   # TRUE
#' test.fn(c(4,3,6,1,7,NA))  # TRUE
#' test.fn(c(4,3,6,NA,1,8))  # FALSE
#' test.fn(c(NA,4,3,6,1,8))  # FALSE
#' test.fn(c(4,3,NA,6,1,NA)) # FALSE
#' @export
test.fn=function(vec) { length(na.contiguous(vec)) == length(na.omit(vec)) &
                         all(na.contiguous(vec)==vec[1:length(na.contiguous(vec))]) }



#' Population pyramid computation
#'
#' Calculating population pyramid based on age-sex-specific survival probabilities and initial population size for age-class 0
#'
#' @param \code{N0} numeric, initial population size for age-class 0
#' @param \code{surv.p} A data frame, with 2 columns "\code{surv.F}" and "\code{surv.M}" giving for each age-sex class the corresponding survival probability
#' @param \code{sex.ratio} numeric, sex ratio of females in N0
#' @param \code{maturity.age.F},\code{maturity.age.M} numeric, sexual maturity age for females and males
#' @param \code{max.age.F},\code{max.age.M} numeric, maximum age for females and males. If not included, it is inferred from the survival probability data.
#'
#' @return a list, with \code{Ntot} (total number of adults), \code{NaF} (total number of adults females), \code{NaM} (total number of adults males), and \code{Nx} (data frame of age-sex specific abundances).
#'
#' @examples
#' N0=182723
#' surv.p=data.frame(  surv.F=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0) ,
#'                     surv.M=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0)
#'                   )
#' sex.ratio=0.444
#' maturity.age.F=4
#' maturity.age.M=4
#' N.res=N.calc.fn(N0,surv.p,sex.ratio,maturity.age.F,maturity.age.M)
#' N.res$Ntot
#' @export
N.calc.fn=function(N0,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F=NULL,max.age.M=NULL){

    if(test.fn(surv.p$surv.F)) {surv.F=na.contiguous(surv.p$surv.F) } else { error("Missing values for some age-specific survival probabilities for females")}
    if(test.fn(surv.p$surv.M)) {surv.M=na.contiguous(surv.p$surv.M) } else { error("Missing values for some age-specific survival probabilities for males")}
    
    if(is.null(max.age.F)) { max.age.F=length(surv.F)+1 }
    if(is.null(max.age.M)) { max.age.M=length(surv.M)+1 }
    
    NxF=sex.ratio*N0*cumprod(c(1,na.contiguous(surv.p$surv.F)))
    NxM=(1-sex.ratio)*N0*cumprod(c(1,na.contiguous(surv.p$surv.M)))

    NaF=sum(NxF[maturity.age.F:max.age.F])
    NaM=sum(NxM[maturity.age.M:max.age.M])
    Ntot=NaF+NaM

    max.age=  max(max.age.F, max.age.M)
    length(NxF) = max.age
    length(NxM) = max.age
    NxF[is.na(NxF)] = 0
    NxM[is.na(NxM)] = 0

    return(list(Ntot=Ntot,NaF=NaF,NaM=NaM,Nx=data.frame(Age=1:max.age,NxF=NxF,NxM=NxM)))

}



#' Incremental search for N0
#'
#' Calculating population pyramid based on age-sex-specific survival probabilities and initial population size for age-class 0
#'
#' @param \code{threshold.N} numeric, adult population threshold target
#' @param \code{surv.p} A data frame, with 2 columns "\code{surv.F}" and "\code{surv.M}" giving for each age-sex class the corresponding survival probability
#' @param \code{sex.ratio} numeric, sex ratio of females in N0
#' @param \code{maturity.age.F},\code{maturity.age.M} numeric, sexual maturity age for females and males
#' @param \code{max.age.F},\code{max.age.M} numeric, maximum age for females and males. If not included, it is inferred from the survival probability data.
#'
#' @return a numeric, \code{N0} estimate for initial population size for age-class 0 based the set of population parameters using a complete incremental search starting at \code{N0=1}.
#'                    Warning: while this approach offers an exhaustive search, it can get really slow if \code{N0} is large.
#'
#' @examples
#' threshold.N=2000
#' surv.p=data.frame(  surv.F=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0) ,
#'                     surv.M=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0)
#'                   )
#' sex.ratio=0.444
#' maturity.age.F=4
#' maturity.age.M=4
#' find.value.N0.inc(threshold.N,surv.p,sex.ratio,maturity.age.F,maturity.age.M)
#' @export
find.value.N0.inc=function(threshold.N,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F=NULL,max.age.M=NULL)
{
    guess=1
    Ntot.est=N.calc.fn(N0=guess,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F,max.age.M)$Ntot

    while (Ntot.est <= threshold.N){
       guess=guess+1
       Ntot.est=N.calc.fn(N0=guess,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F,max.age.M)$Ntot
    }

    return(N0=guess-1)
}



#' Base for N0 optimization function
#' 
#' Provide metric to optimize initial population size for age-class 0 based on adult population threshold target
#'
#' @param \code{N0} numeric, estimate for initial population size for age-class 0
#' @param \code{threshold.N} numeric, adult population threshold target
#' @param \code{surv.p} A data frame, with 2 columns "\code{surv.F}" and "\code{surv.M}" giving for each age-sex class the corresponding survival probability
#' @param \code{sex.ratio} numeric, sex ratio of females in \code{N0}
#' @param \code{maturity.age.F},\code{maturity.age.M} numeric, sexual maturity age for females and males
#' @param \code{max.age.F},\code{max.age.M} numeric, maximum age for females and males. If not included, it is inferred from the survival probability data.
#'
#' @return a numeric, absolute difference between threshold target and simulated total adult population.
#'
#' @examples
#' N0=182723
#' surv.p=data.frame(  surv.F=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0) ,
#'                     surv.M=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0)
#'                   )
#' sex.ratio=0.444
#' maturity.age.F=4
#' maturity.age.M=4
#' N.res=find.value.N0.optim(N0,surv.p,sex.ratio,maturity.age.F,maturity.age.M)
#' N.res$Ntot
#' @export
find.value.N0.optim=function(N0,threshold.N,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F=NULL,max.age.M=NULL)
{
    abs( N.calc.fn(N0,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F,max.age.M)$Ntot - threshold.N )
}


#' Effective population size estimation
#'
#' Calculating population pyramid based on age-sex-specific survival probabilities and initial population size for age-class 0 based on target adult population threshold
#'
#' \code{eff.size.fn} estimate an appropriate initial population size for age-class 0 (\code{N0}) based on a specified targeted adult population threshold.
#' Optimization method to assess correct N0 is controlled by the argument method.
#' Method "\code{BFGS}" uses a quasi-Newton approach (also known as a variable metric algorithm), specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno. See ?optim for more details.
#' Method "\code{Incremental}" performs an exhaustive search for N0 using a complete incremental search starting at \code{N0=1}. For large N0, this approach can get really slow.
#'
#' @param \code{threshold.N} numeric, adult population threshold target
#' @param \code{surv.p} A data frame, with 2 columns "\code{surv.F}" and "\code{surv.M}" giving for each age-sex class the corresponding survival probability
#' @param \code{sex.ratio} numeric, sex ratio of females in \code{N0}
#' @param \code{maturity.age.F},\code{maturity.age.M} numeric, sexual maturity age for females and males
#' @param \code{max.age.F},\code{max.age.M} numeric, maximum age for females and males. If not included, it is inferred from the survival probability data.
#' @param method "\code{BFGS}" or "\code{Incremental}", selects method for estimation of appropriate N0
#' @param round logical, if true, reported abundance are rounded to the nearest integer (\code{default = TRUE})
#'
#' @return a list, with \code{Ntot} (total number of adults), \code{NaF} (total number of adults females), \code{NaM} (total number of adults males), \code{Nx} (data frame of age-sex specific abundances), and \code{N0} (estimated initial population size for age-class 0). A list of input parameter value is also included.
#'
#' @examples
#' surv.p=data.frame(  surv.F=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0) ,
#'                     surv.M=c( 0.501, 0.624, 0.672, 0.696, 0.71 , 0.718, 0.876, 0.88 , 0.883,
#'                               0.885, 0.886, 0.887, 0.888, 0.888, 0.889, 0.889, 0.889, 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 ,
#'                               0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0.89 , 0)
#'                   )
#'
#' maturity.age.F=4
#' maturity.age.M=4
#' sex.ratio=0.444
#'
#' threshold.N.adult=199595
#' 
#' start=Sys.time()
#' N.res=eff.size.fn(threshold.N=threshold.N.adult,surv.p,sex.ratio,maturity.age.F,maturity.age.M,method="BFGS")
#' end=Sys.time()
#' end-start    #  Time difference of 0.171541 secs
#'
#' ### Not Run: ###
#' start=Sys.time()
#' N.res=eff.size.fn(threshold.N=threshold.N.adult,surv.p,sex.ratio,maturity.age.F,maturity.age.M,method="Incremental")
#' end=Sys.time()
#' end-start   #  Time difference of 1.934507 mins
#' @export
eff.size.fn=function(threshold.N,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F=NULL,max.age.M=NULL,method=c("BFGS","Incremental"),round=T){

  if(method=="BFGS"){  N0.est=optim(200000,
                                 find.value.N0.optim,
                                 threshold.N=threshold.N,surv.p=surv.p,sex.ratio=sex.ratio,
                                 maturity.age.F=maturity.age.F,maturity.age.M=maturity.age.M,
                                 max.age.F=max.age.F,max.age.M=max.age.M, method ="BFGS") $par
                     } else  { if(method=="Incremental") {  warning("Incremental optimization can be really slow")
                                                            N0.est=find.value.N0.inc(threshold.N,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F=NULL,max.age.M=NULL)
                                                         }  else  { error("Incorrect optimization method selected") } }

                    

  res=N.calc.fn(N0=floor(N0.est),surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F,max.age.M)

  if(round==T) { res=N.calc.fn(N0=floor(N0.est),surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F,max.age.M)
                 res=lapply(append(res,list("N0"=floor(N0.est))),round) }    else   {
                              res=N.calc.fn(N0=N0.est,surv.p,sex.ratio,maturity.age.F,maturity.age.M,max.age.F,max.age.M)
                              res=append(res,list("N0"=N0.est))
                              }
  return(append(res,list("par"=list(threshold.N=threshold.N,surv.p=surv.p,sex.ratio=sex.ratio,
                                    maturity.age.F=maturity.age.F,maturity.age.M=maturity.age.M,
                                    max.age.F=max.age.F,max.age.M=max.age.M,method=method))))
  }

