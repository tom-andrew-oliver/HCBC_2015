#HCBC Bleaching Exploration
library(ggplot2)
library(plotly)
library(leaflet)
library(RColorBrewer)
library(htmltools)
library(sp)
library(raster)
library(ncdf4)
library(rasterVis)
rep=data.frame(rbind(c(NaN,NA),c(Inf,NA),c(-Inf,NA)))
source("C:/Users/Thomas.Oliver/WORK/Projects/Environmental Data Summary/HelperCode/gcdist.R")
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


#read in HCBC Dataset, convert to SpatialPointsDataFrame, and set the proj4string (i.e. Coordinate Reference System) of the dataset
B=read.csv("/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_Cluster.csv")
coordinates(B)=~Longitude_mn+Latitude_mn
#Lat Long, standard WGS84 CRS
proj4string(B)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"


###############################
###############################
# SatLayers - DHW
###############################
###############################
dhw_v2=raster("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/DHW/version2/dhw_b05kmnn_20151101.nc")
dhw_v3=raster("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/DHW/version3/b5km_dhw_20151101.nc")
mur_dhw_15=raster("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/DHW_MUR/MHI_MUR_DHW_MAX2015.grd")
mur_dhw_14=raster("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/DHW_MUR/MHI_MUR_DHW_MAX2014.grd")

proj4string(dhw_v2)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj4string(dhw_v3)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj4string(mur_dhw_14)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
proj4string(mur_dhw_15)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
dhw_v2.=crop(dhw_v2,extent(B))
dhw_v3.=crop(dhw_v3,extent(B))
mur_dhw_14.=crop(mur_dhw_14,extent(B))
mur_dhw_15.=crop(mur_dhw_15,extent(B))

#Extract Nearest non-NA value from raster
Dists=seq(10,10000,by=100)
samp_DHWv2=ExpandingExtract(dhw_v2.,B,Dists = Dists)
samp_DHWv3=ExpandingExtract(dhw_v3.,B,Dists = Dists)
samp_DHWmur14=ExpandingExtract(mur_dhw_14.,B,Dists = Dists)
samp_DHWmur15=ExpandingExtract(mur_dhw_15.,B,Dists = Dists)

B$DHW_v2_20151101=samp_DHWv2$values
B$DHW_v3_20151101=samp_DHWv3$values
B$DHW_MUR_20141101=samp_DHWmur14$values
B$DHW_MUR_20151101=samp_DHWmur15$values

###############################
###############################
# ModelLayers - Proportion of Bleaching Resistant Taxa
###############################
###############################
P_brestaxa=raster("/Users/thomas.oliver/WORK/Projects/HCBC/PropBleachingResistantTaxa/resistant_coral_predicted_gcs84.tif")

proj4string(P_brestaxa)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

P_brestaxa.=crop(P_brestaxa,extent(B))

#Extract Nearest non-NA value from raster
Dists=seq(10,10000,by=100)
samp_Pbrt=ExpandingExtract(P_brestaxa.,B,Dists = Dists)

B$Prop_BleachResTaxa=samp_Pbrt$values
write.csv(B,"C:/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_DHW_pBRT.csv")

###############################
###############################
# Lecky Layers
###############################
###############################

#Deal With Lecky Layers
LeckyPath="/Users/thomas.oliver/WORK/Projects/HCBC/LeckyLayer/FromJoey/"
LeckyFileList=list.files(LeckyPath,recursive = TRUE)
LenLFL=nchar(LeckyFileList)
LeckyLayerList=LeckyFileList[grep(".tif",substr(LeckyFileList,LenLFL-3,LenLFL))]


#coordinates(B)=c(B$Longitude_mn,B$Latitude_mn)
LLval=data.frame(val=rep(NA,nrow(B)))
LLDist=data.frame(val=rep(NA,nrow(B)))
LLN=data.frame(val=rep(NA,nrow(B)))
for(i in 1:length(LeckyLayerList)){
  r=raster(paste0(LeckyPath,LeckyLayerList[i]))
  r.=projectRaster(r,crs=crs(B))
  df=ExpandingExtract(r.,B,Dists=c(10,500,1000,1500,2000,4000,8000))
  LLval[,i]=df$values
  LLDist[,i]=df$Dist
  LLN[,i]=df$N
  names(LLval)[i]=names(r)
  names(LLDist)[i]=names(r)
  names(LLN)[i]=names(r)
  print(paste0("Completed ",i," extractions of ",length(LeckyLayerList)," layers."))
}


r=raster(paste0(LeckyPath,LeckyLayerList[grep("Sediment",LeckyLayerList)]))
r.=projectRaster(r,crs=crs(B))
r..=subs(r.,data.frame("NaN","NA"),subsWithNA=FALSE)
sed=ExpandingExtract(r..,B,Dists=c(10,500,1000,seq(1250,10000,by=250)))


r=raster(paste0(LeckyPath,LeckyLayerList[grep("OTP_MHI_Fishing_Commercial_Net",LeckyLayerList)]))
r.=projectRaster(r,crs=crs(B))
r..=subs(r.,data.frame("NaN","NA"),subsWithNA=FALSE)
net=ExpandingExtract(r..,B,Dists=c(10,500,1000,seq(1250,10000,by=250)))

image(r..)
points(B[which(is.nan(sed$values)),],pch=6)

LLval$OTP_MHI_Sedimentation=sed$values
LLval$OTP_MHI_Fishing_Commercial_Net=net$values

write.csv(LLval,paste0("/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/LeckyExtractions_1km_2020.csv"))

B@data[,names(LLval)]=LLval
names(LLval)%in%names(B)


#WAVE EXPOSURE
B.=as.data.frame(B)
WE=read.csv("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/WaveExposure/15m_contours.csv")
WE=subset(WE,BAD_FLAG==0)
WE$mn1979_2012=rowMeans(WE[,4:37])
Ds=gcdist.set1vset2(B.$Latitude_mn,B.$Longitude_mn,WE$y,WE$x)
WE_i=apply(Ds,FUN = which.min,MARGIN = 1)
WE_vec=WE$mn1979_2012[WE_i]

class(B)
B$WaveEnergy_MN1979.2012=WE_vec
names(B)
write.csv(B,paste0("C:/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_OTP_WE_DHW_pBRT.csv"))

RELOAD=FALSE
if(RELOAD){
  B=read.csv("/Users/thomas.oliver/WORK/Projects/HCBC/BleachingData/Canon/HCBC_1000m_SingleCluster_DHW.LL.WE.SST.20170829.csv")
  coordinates(B)=~Longitude_mn+Latitude_mn
  #Lat Long, standard WGS84 CRS
  proj4string(B)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
}
###############################
###############################
# K490
###############################
###############################
K490.MHI=brick("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/K490/MHI.K490_LayerBrick.grd")
K490.MHI.=subs(K490.MHI,rep,subsWithNA=FALSE)


K490val=data.frame(val=rep(NA,nrow(B)))
K490Dist=data.frame(val=rep(NA,nrow(B)))
K490N=data.frame(val=rep(NA,nrow(B)))
for(i in 1:nlayers(K490.MHI.)){
  df=ExpandingExtract(K490.MHI.[[i]],B,Dists=c(10,500,1000,1500,2000,3000,4000,8000))
  K490val[,i]=df$values
  K490Dist[,i]=df$Dist
  K490N[,i]=df$N
  print(paste0("Completed ",i," extractions of ",nlayers(K490.MHI.)," layers."))
}
names(K490val)=names(K490.MHI.)
names(K490Dist)=names(K490.MHI.)
names(K490N)=names(K490.MHI.)

B@data[,names(K490val)]=K490val


###############################
###############################
# PAR
###############################
###############################
PAR.MHI=brick("/Users/thomas.oliver/WORK/Projects/HCBC/SatData/PAR/MHI.PAR_LayerBrick.grd")
#Replace Inf and NaN with NA
rep=data.frame(rbind(c(NaN,NA),c(Inf,NA),c(-Inf,NA)))
PAR.MHI.=subs(PAR.MHI,rep,subsWithNA=FALSE)

PARval=data.frame(val=rep(NA,nrow(B)))
PARDist=data.frame(val=rep(NA,nrow(B)))
PARN=data.frame(val=rep(NA,nrow(B)))
for(i in 1:nlayers(PAR.MHI.)){
  df=ExpandingExtract(PAR.MHI.[[i]],B,Dists=c(10,500,1000,1500,2000,3000,4000,8000))
  PARval[,i]=df$values
  PARDist[,i]=df$Dist
  PARN[,i]=df$N
  print(paste0("Completed ",i," extractions of ",nlayers(PAR.MHI.)," layers."))
}
names(PARval)=names(PAR.MHI)
names(PARDist)=names(PAR.MHI)
names(PARN)=names(PAR.MHI)

B@data[,names(PARval)]=PARval

head(B)
write.csv(B,paste0("C:/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_PAR_K490_OTP_WE_DHW_pBRT.csv"))

###############################
###############################
# Get MUR SST Data
###############################
###############################
rMHI.brick=brick("I:/SatelliteDataByIsland/GridData/SST_MUR_1km_Daily_2002_2016/MHILayers/MHI.MURsst.LT.2015.2014.Layers.grd")
writeRaster(rMHI.brick,filename = "/Users/thomas.oliver/WORK/Projects/SatelliteDataByIsland/GridData/SST_MUR_1km_Daily_2002_2016/MHILayers/MHI.MURsst.LT.2015.2014.Layers.tif",format="GTiff")
#Mask Land Values
rMHI.brick.mask=rMHI.brick#mask(x=rMHI.brick,mask=rMHI.mask,maskval=2)
plot(rMHI.brick)
#plot(rMHI.mask)
plot(rMHI.brick.mask)

#Replace Inf and NaN with NA
rep=data.frame(rbind(c(NaN,NA),c(Inf,NA),c(-Inf,NA)))
rMHI.brick.mask.=subs(rMHI.brick.mask,rep,subsWithNA=FALSE)

#coordinates(B)=c(B$Longitude_mn,B$Latitude_mn)
MURval=data.frame(val=rep(NA,nrow(B)))
MURDist=data.frame(val=rep(NA,nrow(B)))
MURN=data.frame(val=rep(NA,nrow(B)))
for(i in 1:nlayers(rMHI.brick.mask.)){
  df=ExpandingExtract(rMHI.brick.mask.[[i]],B,Dists=c(10,500,1000,1500,2000,3000,4000,8000))
  MURval[,i]=df$values
  MURDist[,i]=df$Dist
  MURN[,i]=df$N
  print(paste0("Completed ",i," extractions of ",nlayers(rMHI.brick.mask.)," layers."))
}
names(MURval)=names(rMHI.brick.mask.)
names(MURDist)=names(rMHI.brick.mask.)
names(MURN)=names(rMHI.brick.mask.)

B@data[,names(MURval)]=MURval
names(B)=sub("rMHI","SST",names(B))
names(B)

write.csv(B,paste0("C:/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_MUR_PAR_K490_OTP_WE_DHW_pBRT.csv"))

PARi=grep("PAR",names(B))
K490i=grep("K490",names(B))
#Now generate PARz according to relationship in Zaneveld et al 1993
#PARz=par$par*exp(-(0.0085+1.6243*(k490$k490*(0.3048*Obs_OCC$BenthicDepthMean[i]))))PARz=calc
parzFUN=function(par,k490,depth){
  return(par*exp(-(0.0085+1.6243*(k490*(0.3048*depth)))))
}


PARz=B@data[,PARi]*exp(-(0.0085+1.6243*(B@data[,K490i]*(0.3048*B$Depth_m_mn))))
dim(PARz)
names(PARz)=sub("PAR","PARz",names(PARz))
B@data=cbind(B@data,PARz)
write.csv(as.data.frame(B),"/Users/thomas.oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_PARz_MUR_PAR_K490_OTP_WE_DHW_pBRT.csv",row.names = FALSE)


########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

