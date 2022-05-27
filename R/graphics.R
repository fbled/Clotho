#' Pyramid plot
#'
#' Plots population pyramid plot from "\code{eff.size.fn}" output
#'
#' @param \code{N.res} output from a "\code{eff.size.fn}" computation
#'
#' @return a population pyramid plot, with \code{Ntot} (total number of adults), \code{N0} (estimated initial population size for age-class 0)
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
#' N.res=eff.size.fn(threshold.N=threshold.N.adult,surv.p,sex.ratio,maturity.age.F,maturity.age.M,method="BFGS")
#' pyramid.plot(N.res)
#' @export

pyramid.plot=function(N.res) {

    pyramid.data=reshape(N.res$Nx, idvar = "Age", varying = list(c("NxM","NxF")),
                         v.names = "Nx", direction = "long",times=c("Male","Female"),timevar="Sex")
    
    pyramid.data$Class=NA
     pyramid.data$Class[pyramid.data$Sex=="Male"]=ifelse(pyramid.data$Age[pyramid.data$Sex=="Male"]<N.res$ par $ maturity.age.M,"Juvenile","Adult")
     pyramid.data$Class[pyramid.data$Sex=="Female"]=ifelse(pyramid.data$Age[pyramid.data$Sex=="Female"]<N.res$ par $ maturity.age.F,"Juvenile","Adult")
    
    pyramid.data$N=ifelse(pyramid.data$Sex=="Male",-pyramid.data$Nx,pyramid.data$N)
    
    pyramid.data$Sex=factor(pyramid.data$Sex,levels=c("Male","Female"))
    
    
    
    par(mfrow=c(1,2),family = "mono")
    par(mar= c(5, 2, 2.5, 0))
    barplot(N~Age,data=subset(pyramid.data,Sex=="Male"),horiz=T,xaxt="n",yaxt="n",xlab="", ylab="",
                  col="transparent",
                  border = NA)
    do.call("segments",list(x0=axTicks(1), y0=0, y1 = 1.2*max(N.res$ Nx $ Age),col="lightgrey",lwd=0.5))
    axis(1, las=1,at = axTicks(1),labels=format(abs(axTicks(1)), scientific=FALSE),tck=F,col = "transparent" ,line= -1 ,cex.axis=0.75)
    
    barplot(N~Age,data=subset(pyramid.data,Sex=="Male"),horiz=T,xaxt="n",yaxt="n",xlab="", ylab="",
                  col=c(rep("#C1E8F2",N.res$par$maturity.age.F-1),rep("#6AC7DF",(max(pyramid.data$Age)-N.res$par$maturity.age.F+1))),
                  border = NA,add=T)
    title(sub="Male",line=2.25,cex.sub=1.5)
    

    legend("topleft",legend= as.expression(c( substitute(paste(N[adult],txt,n), list(txt=" = ",n=N.res$Ntot)),
                        substitute(paste(N[1],txt,n), list(txt=" = ",n=N.res$N0)))), 
                        box.col = "white", bg = "white",cex=1.25,xjust=0, y.intersp=1.5,adj=0)
    
    par(mar= c(5, 2.125, 2.5, 2))
    barplot(N~Age,data=subset(pyramid.data,Sex=="Female"),horiz=T,xaxt="n",yaxt="n",xlab="", ylab="",
                  col="transparent", 
                  border = NA)
    do.call("segments",list(x0=axTicks(1), y0=0, y1 = 1.2*max(N.res$ Nx $ Age),col="lightgrey",lwd=0.5))
    axis(1, las=1,at = axTicks(1),labels=format(abs(axTicks(1)), scientific=FALSE),tck=F,col = "transparent" ,line= -1 ,cex.axis=0.75 )
    barplot(N~Age,data=subset(pyramid.data,Sex=="Female"),horiz=T,xaxt="n",yaxt="n",xlab="", ylab="",
                  col=c(rep("#C3E5DE",N.res$par$maturity.age.F-1),rep("#5DDFC4",(max(pyramid.data$Age)-N.res$par$maturity.age.F+1))), 
                  border = NA, add=T)
    
    title(sub="Female",line=2.25,cex.sub=1.5)
    
    axis(2, las=1,at = 1.2*c(1,seq(5,max(N.res$ Nx $ Age),by=5))-0.5,labels=c(1,seq(5,max(N.res$ Nx $ Age),by=5)),tck=F,col = "transparent" ,hadj=0.5  ,cex.axis=1.125)

}


