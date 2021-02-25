library(DAAG)
library(EnvStats)
library(factoextra)
library(geoR)
library(GGally)
library(ggplot2)
library(gridExtra)
library(htmltools)
library(htmlwidgets)
library(leaflet)
library(lubridate)
library(MASS)
library(mgcv)
library(MuMIn)
library(ncdf4)
library(plotly)
library(plyr)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(relaimpo)
library(reshape)
library(sp)

#Define Functions:
lengthNONA=function(x){return(length(na.omit(x)))}

ExpandingExtract=function(r,SpDF,Dists=c(500,1000,2000,4000,8000)){
  OutDF=data.frame(values=rep(NA,nrow(SpDF)),Dist=rep(NA,nrow(SpDF)),N=rep(NA,nrow(SpDF)))
  nDists=length(Dists)
  cnt=1
  NAi=which(is.na(OutDF$values))
  NAsLeft=length(NAi)>0
  while(cnt<=nDists&NAsLeft){
    NAi=which(is.na(OutDF$values))
    pull=extract(x=r,y=SpDF[NAi,],buffer=Dists[cnt],na.rm=TRUE)
    Nper=unlist(lapply(pull,lengthNONA))
    OutDF$values[NAi]=unlist(lapply(pull,mean,na.rm=TRUE))
    OutDF$Dist[NAi]=Dists[cnt]
    OutDF$N[NAi]=Nper
    NAi=which(is.na(OutDF$values))
    NAsLeft=length(NAi)>0
    cnt=cnt+1
  }
  return(OutDF)
}
quant99=function(x) return(quantile(x,probs=0.99,na.rm=TRUE))
quant01=function(x) return(quantile(x,probs=0.01,na.rm=TRUE))
diff.rng=function(x,na.rm=TRUE) return(abs(diff(range(x,na.rm=na.rm))))
wkly_rng=function(x){
  diff.rng=function(x) return(abs(diff(range(x,na.rm=TRUE))))
  out=rollapply(x,7,diff.rng)
  return(mean(out,na.rm=TRUE))
}

GamThemAll=function(y,xes,df,k=3){
  #y=y[1:nrow(df)]
  require(mgcv)
  str=paste("mod=gam(",y,"~s(",sep="")
  collapseVal=paste(",k=",k,")+s(",sep="")
  str=paste(str,paste(xes,collapse=collapseVal,sep=""),sep="")
  str=paste(str,",k=",k,"),data=",df,",na.action='na.fail')",sep="")
  eval(parse(text=str))
  return(mod)
}
GamThemAll2=function(y,xes,dfdata,k=3){
  #y=y[1:nrow(df)]
  require(mgcv)
  str=paste("mod=gam(",y,"~s(",sep="")
  collapseVal=paste(",k=",k,")+s(",sep="")
  str=paste(str,paste(xes,collapse=collapseVal,sep=""),sep="")
  str=paste(str,",k=",k,"),data=dfdata,na.action='na.fail')",sep="")
  eval(parse(text=str))
  return(mod)
}
size_scale=function(values,base=1,vmin=NULL,vmax=NULL,scale=1){
  if(is.null(vmin)) vmin=min(values,na.rm=TRUE)
  if(is.null(vmax)) vmax=max(values,na.rm=TRUE)
  out=base+scale*(values-vmin)/(vmax-vmin)
  return(out)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}