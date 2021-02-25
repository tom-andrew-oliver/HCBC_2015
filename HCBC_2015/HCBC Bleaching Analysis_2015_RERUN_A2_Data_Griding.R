#HCBC Bleaching Exploration
source("./HCBC_2015/HelperCode/gcdist.R")
source("./HCBC_2015/HelperCode/HCBC_2015_Functions.R")

#rc2cv_index=function(r,c,nrow_matrix=max(r)){return((c-1)*nrow_matrix+r)}

inpath="./HCBC_2015/"
#read in HCBC Dataset
HCBC=read.csv(paste0(inpath,"HCBC_MHI_2015_10.11.csv"))

#Generate Meaningful Dates & Ordered Depths Name
HCBC$ObservationDate_R=ymd(HCBC$ObservationDate)


##########################################################################################
# Grid Data Points By 1000 m distances, maintain DepthBins
##########################################################################################

# ##########################################################################################
# # Rasterize Data Points - Ignore Depth for now...
# ##########################################################################################
coordinates(HCBC)=cbind(HCBC$Longitude_DD,HCBC$Latitude_DD)

#2 km grid
SpatialResolution_Deg_2km=2/111#X km meter bins, maintaining Depth Bins
nc_2km=round(diff(extent(HCBC)[1:2])/SpatialResolution_Deg_2km)
nr_2km=round(diff(extent(HCBC)[3:4])/SpatialResolution_Deg_2km)
Bgrid_2km = raster(extent(HCBC),
                   ncols=nc_2km,
                   nrows=nr_2km,
                   vals=1:(nc_2km*nr_2km))

#1 km grid
SpatialResolution_Deg_1km=1/111#X km meter bins, maintaining Depth Bins
nc_1km=round(diff(extent(HCBC)[1:2])/SpatialResolution_Deg_1km)
nr_1km=round(diff(extent(HCBC)[3:4])/SpatialResolution_Deg_1km)
Bgrid_1km = raster(extent(HCBC),
               ncols=nc_1km,
               nrows=nr_1km,
               vals=1:(nc_1km*nr_1km))

#0.5 km grid
SpatialResolution_Deg_500m=.5/111#X km meter bins, maintaining Depth Bins
nc_500m=round(diff(extent(HCBC)[1:2])/SpatialResolution_Deg_500m)
nr_500m=round(diff(extent(HCBC)[3:4])/SpatialResolution_Deg_500m)
Bgrid_500m = raster(extent(HCBC),
                   ncols=nc_500m,
                   nrows=nr_500m,
                   vals=1:(nc_500m*nr_500m))

#0.25 km grid
SpatialResolution_Deg_250m=.25/111#X km meter bins, maintaining Depth Bins
nc_250m=round(diff(extent(HCBC)[1:2])/SpatialResolution_Deg_250m)
nr_250m=round(diff(extent(HCBC)[3:4])/SpatialResolution_Deg_250m)
Bgrid_250m = raster(extent(HCBC),
                   ncols=nc_250m,
                   nrows=nr_250m,
                   vals=1:(nc_250m*nr_250m))

#Check Data/Grid Overlap:
#2km
RasterLookUp_2km=as.data.frame(rasterToPoints(Bgrid_2km))
HCBC$RasterAddress_2km=extract(Bgrid_2km,HCBC)
HCBC$Longitude_RAS_2km=RasterLookUp_2km[match(HCBC$RasterAddress_2km,RasterLookUp_2km$layer),"x"]
HCBC$Latitude_RAS_2km=RasterLookUp_2km[match(HCBC$RasterAddress_2km,RasterLookUp_2km$layer),"y"]
#1km
RasterLookUp_1km=as.data.frame(rasterToPoints(Bgrid_1km))
HCBC$RasterAddress_1km=extract(Bgrid_1km,HCBC)
HCBC$Longitude_RAS_1km=RasterLookUp_1km[match(HCBC$RasterAddress_1km,RasterLookUp_1km$layer),"x"]
HCBC$Latitude_RAS_1km=RasterLookUp_1km[match(HCBC$RasterAddress_1km,RasterLookUp_1km$layer),"y"]
#500m
RasterLookUp_500m=as.data.frame(rasterToPoints(Bgrid_500m))
HCBC$RasterAddress_500m=extract(Bgrid_500m,HCBC)
HCBC$Longitude_RAS_500m=RasterLookUp_500m[match(HCBC$RasterAddress_500m,RasterLookUp_500m$layer),"x"]
HCBC$Latitude_RAS_500m=RasterLookUp_500m[match(HCBC$RasterAddress_500m,RasterLookUp_500m$layer),"y"]
#250m
RasterLookUp_250m=as.data.frame(rasterToPoints(Bgrid_250m))
HCBC$RasterAddress_250m=extract(Bgrid_250m,HCBC)
HCBC$Longitude_RAS_250m=RasterLookUp_250m[match(HCBC$RasterAddress_250m,RasterLookUp_250m$layer),"x"]
HCBC$Latitude_RAS_250m=RasterLookUp_250m[match(HCBC$RasterAddress_250m,RasterLookUp_250m$layer),"y"]

paste0("2km: ",length(unique(HCBC$RasterAddress_2km))," points")
paste0("1km: ",length(unique(HCBC$RasterAddress_1km))," points")
paste0("500m: ",length(unique(HCBC$RasterAddress_500m))," points")
paste0("250m: ",length(unique(HCBC$RasterAddress_250m))," points")

r_2km <-rasterize(HCBC,Bgrid_2km,HCBC$PctAffected,fun=mean,na.rm=TRUE)
crs(r_2km)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
r_1km <-rasterize(HCBC,Bgrid_1km,HCBC$PctAffected,fun=mean,na.rm=TRUE)
crs(r_1km)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
r_500m <-rasterize(HCBC,Bgrid_500m,HCBC$PctAffected,fun=mean,na.rm=TRUE)
crs(r_500m)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
r_250m <-rasterize(HCBC,Bgrid_250m,HCBC$PctAffected,fun=mean,na.rm=TRUE)
crs(r_250m)="+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

pal <- colorNumeric(c("green", "red"), values(r_250m),
                    na.color = "transparent")

PctAffGrid_2km=leaflet() %>% 
  addProviderTiles(provider = "Esri.WorldImagery") %>%
  addRasterImage(r_2km, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(r_2km),
            title = "Percent Affected")#%>%
PctAffGrid_2km

saveWidget(PctAffGrid_2km,
           file=paste0(inpath,"PctAffectedMaps/HCBC_2015_PctAffect_2000m.htm"))

PctAffGrid_1km=leaflet() %>% 
  addProviderTiles(provider = "Esri.WorldImagery") %>%
  addRasterImage(r_1km, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(r_1km),
            title = "Percent Affected")#%>%
PctAffGrid_1km

saveWidget(PctAffGrid_1km,
           file=paste0(inpath,"PctAffectedMaps/HCBC_2015_PctAffect_1000m.htm"))

PctAffGrid_500m=leaflet() %>% 
  addProviderTiles(provider = "Esri.WorldImagery") %>%
  addRasterImage(r_500m, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(r_500m),
            title = "Percent Affected")#%>%
PctAffGrid_500m

saveWidget(PctAffGrid_500m,
           file=paste0(inpath,"PctAffectedMaps/HCBC_2015_PctAffect_500m.htm"))

PctAffGrid_250m=leaflet() %>% 
  addProviderTiles(provider = "Esri.WorldImagery") %>%
  addRasterImage(r_250m, colors = pal, opacity = 0.8) %>%
  addLegend(pal = pal, values = values(r_250m),
            title = "Percent Affected")#%>%
PctAffGrid_250m

saveWidget(PctAffGrid_250m,
           file=paste0(inpath,"PctAffectedMaps/HCBC_2015_PctAffect_250m.htm"))


#Output Data with Grid COlumns
HCBCdf=as.data.frame(HCBC)
write.csv(x = HCBCdf,file = paste0(inpath,"HCBC_MHI_2015_10.11_GridColumnsAdded.csv"))


#Output Clustered Data
HCBC_cluster_2km=ddply(HCBCdf,.(Island_Name,RasterAddress_2km,DepthBin_5m),summarize,
                       Depth_m_mn=median(Depth_m,na.rm=TRUE),
                       Longitude_mn=mean(Longitude_DD,na.rm=TRUE),
                       Latitude_mn=mean(Latitude_RAS_2km,na.rm=TRUE),
                       Date_mn=mean(ObservationDate_R,na.rm=TRUE),
                       Longitude_ras=mean(Longitude_RAS_2km,na.rm=TRUE),
                       Latitude_ras=mean(Latitude_DD,na.rm=TRUE),
                       AreaSurveyed_m2_sum=sum(AreaSurveyed_m2,na.rm=TRUE),
                       PercentLiveCoralCover_mean=mean(PercentLiveCoralCover,na.rm=TRUE),
                       PercentLiveCoralCover_sd=sd(PercentLiveCoralCover,na.rm=TRUE),
                       PercentLiveCoralCover_N=length(PercentLiveCoralCover),
                       PctCoralUnbleached_mean=mean(PctCoralUnbleached,na.rm=TRUE),
                       PctCoralUnbleached_sd=sd(PctCoralUnbleached,na.rm=TRUE),
                       PctCoralUnbleached_N=length(PctCoralUnbleached),
                       PctCoralPartialBleached_mean=mean(PctCoralPartialBleached,na.rm=TRUE),
                       PctCoralPartialBleached_sd=sd(PctCoralPartialBleached,na.rm=TRUE),
                       PctCoralPartialBleached_N=length(PctCoralPartialBleached),
                       PctCoralFullyBleached_mean=mean(PctCoralFullyBleached,na.rm=TRUE),
                       PctCoralFullyBleached_sd=sd(PctCoralFullyBleached,na.rm=TRUE),
                       PctCoralFullyBleached_N=length(PctCoralFullyBleached),
                       PctAffected_mean=mean(PctAffected,na.rm=TRUE),
                       PctAffected_sd=sd(PctAffected,na.rm=TRUE),
                       PctAffected_N=length(PctAffected))
write.csv(x = HCBC_cluster_2km,file = paste0(inpath,"ClusteredData/HCBC_2015_2km_Cluster.csv"))

#Output Clustered Data
HCBC_cluster_1km=ddply(HCBCdf,.(Island_Name,RasterAddress_1km,DepthBin_5m),summarize,
                   Depth_m_mn=median(Depth_m,na.rm=TRUE),
                   Longitude_mn=mean(Longitude_DD,na.rm=TRUE),
                   Latitude_mn=mean(Latitude_RAS_1km,na.rm=TRUE),
                   Date_mn=mean(ObservationDate_R,na.rm=TRUE),
                   Longitude_ras=mean(Longitude_RAS_1km,na.rm=TRUE),
                   Latitude_ras=mean(Latitude_DD,na.rm=TRUE),
                   AreaSurveyed_m2_sum=sum(AreaSurveyed_m2,na.rm=TRUE),
                   PercentLiveCoralCover_mean=mean(PercentLiveCoralCover,na.rm=TRUE),
                   PercentLiveCoralCover_sd=sd(PercentLiveCoralCover,na.rm=TRUE),
                   PercentLiveCoralCover_N=length(PercentLiveCoralCover),
                   PctCoralUnbleached_mean=mean(PctCoralUnbleached,na.rm=TRUE),
                   PctCoralUnbleached_sd=sd(PctCoralUnbleached,na.rm=TRUE),
                   PctCoralUnbleached_N=length(PctCoralUnbleached),
                   PctCoralPartialBleached_mean=mean(PctCoralPartialBleached,na.rm=TRUE),
                   PctCoralPartialBleached_sd=sd(PctCoralPartialBleached,na.rm=TRUE),
                   PctCoralPartialBleached_N=length(PctCoralPartialBleached),
                   PctCoralFullyBleached_mean=mean(PctCoralFullyBleached,na.rm=TRUE),
                   PctCoralFullyBleached_sd=sd(PctCoralFullyBleached,na.rm=TRUE),
                   PctCoralFullyBleached_N=length(PctCoralFullyBleached),
                   PctAffected_mean=mean(PctAffected,na.rm=TRUE),
                   PctAffected_sd=sd(PctAffected,na.rm=TRUE),
                   PctAffected_N=length(PctAffected))
write.csv(x = HCBC_cluster_1km,file = paste0(inpath,"ClusteredData/HCBC_2015_1km_Cluster.csv"))


#Output Clustered Data
HCBC_cluster_500m=ddply(HCBCdf,.(Island_Name,RasterAddress_500m,DepthBin_5m),summarize,
                       Depth_m_mn=median(Depth_m,na.rm=TRUE),
                       Longitude_mn=mean(Longitude_DD,na.rm=TRUE),
                       Latitude_mn=mean(Latitude_RAS_500m,na.rm=TRUE),
                       Date_mn=mean(ObservationDate_R,na.rm=TRUE),
                       Longitude_ras=mean(Longitude_RAS_500m,na.rm=TRUE),
                       Latitude_ras=mean(Latitude_DD,na.rm=TRUE),
                       AreaSurveyed_m2_sum=sum(AreaSurveyed_m2,na.rm=TRUE),
                       PercentLiveCoralCover_mean=mean(PercentLiveCoralCover,na.rm=TRUE),
                       PercentLiveCoralCover_sd=sd(PercentLiveCoralCover,na.rm=TRUE),
                       PercentLiveCoralCover_N=length(PercentLiveCoralCover),
                       PctCoralUnbleached_mean=mean(PctCoralUnbleached,na.rm=TRUE),
                       PctCoralUnbleached_sd=sd(PctCoralUnbleached,na.rm=TRUE),
                       PctCoralUnbleached_N=length(PctCoralUnbleached),
                       PctCoralPartialBleached_mean=mean(PctCoralPartialBleached,na.rm=TRUE),
                       PctCoralPartialBleached_sd=sd(PctCoralPartialBleached,na.rm=TRUE),
                       PctCoralPartialBleached_N=length(PctCoralPartialBleached),
                       PctCoralFullyBleached_mean=mean(PctCoralFullyBleached,na.rm=TRUE),
                       PctCoralFullyBleached_sd=sd(PctCoralFullyBleached,na.rm=TRUE),
                       PctCoralFullyBleached_N=length(PctCoralFullyBleached),
                       PctAffected_mean=mean(PctAffected,na.rm=TRUE),
                       PctAffected_sd=sd(PctAffected,na.rm=TRUE),
                       PctAffected_N=length(PctAffected))
write.csv(x = HCBC_cluster_500m,file = paste0(inpath,"ClusteredData/HCBC_2015_500m_Cluster.csv"))


#Output Clustered Data
HCBC_cluster_250m=ddply(HCBCdf,.(Island_Name,RasterAddress_250m,DepthBin_5m),summarize,
                       Depth_m_mn=median(Depth_m,na.rm=TRUE),
                       Longitude_mn=mean(Longitude_DD,na.rm=TRUE),
                       Latitude_mn=mean(Latitude_RAS_250m,na.rm=TRUE),
                       Date_mn=mean(ObservationDate_R,na.rm=TRUE),
                       Longitude_ras=mean(Longitude_RAS_250m,na.rm=TRUE),
                       Latitude_ras=mean(Latitude_DD,na.rm=TRUE),
                       AreaSurveyed_m2_sum=sum(AreaSurveyed_m2,na.rm=TRUE),
                       PercentLiveCoralCover_mean=mean(PercentLiveCoralCover,na.rm=TRUE),
                       PercentLiveCoralCover_sd=sd(PercentLiveCoralCover,na.rm=TRUE),
                       PercentLiveCoralCover_N=length(PercentLiveCoralCover),
                       PctCoralUnbleached_mean=mean(PctCoralUnbleached,na.rm=TRUE),
                       PctCoralUnbleached_sd=sd(PctCoralUnbleached,na.rm=TRUE),
                       PctCoralUnbleached_N=length(PctCoralUnbleached),
                       PctCoralPartialBleached_mean=mean(PctCoralPartialBleached,na.rm=TRUE),
                       PctCoralPartialBleached_sd=sd(PctCoralPartialBleached,na.rm=TRUE),
                       PctCoralPartialBleached_N=length(PctCoralPartialBleached),
                       PctCoralFullyBleached_mean=mean(PctCoralFullyBleached,na.rm=TRUE),
                       PctCoralFullyBleached_sd=sd(PctCoralFullyBleached,na.rm=TRUE),
                       PctCoralFullyBleached_N=length(PctCoralFullyBleached),
                       PctAffected_mean=mean(PctAffected,na.rm=TRUE),
                       PctAffected_sd=sd(PctAffected,na.rm=TRUE),
                       PctAffected_N=length(PctAffected))
write.csv(x = HCBC_cluster_250m,file = paste0(inpath,"ClusteredData/HCBC_2015_250m_Cluster.csv"))

