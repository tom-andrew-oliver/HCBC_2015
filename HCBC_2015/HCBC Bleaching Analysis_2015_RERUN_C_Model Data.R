rm(list=ls())
####PCA Regression
library(mgcv)
library(relaimpo)
library(MuMIn)

library(RColorBrewer)
library(gridExtra)
library(DAAG)
library(ggplot2)

library(sp)
library(raster)
library(ncdf4)
library(factoextra)

source("C:/Users/Thomas.Oliver/OLD_DRIVE/thomas.oliver/WORK/CRED_WORK/R Scripts/better.biplot.R")
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

###########################################################################################
###########################################################################################
###########################################################################################

B=read.csv("C:/Users/Thomas.Oliver/WORK/Projects/HCBC/RERUN2019/ClusteredData/HCBC_2015_1km_PARz_MUR_PAR_K490_OTP_WE_DHW_pBRT.csv")
B$DepthBin_5m=factor(B$DepthBin_5m,levels=rev(unique(B$DepthBin_5m)))
head(B)
#Remove shitty NA datapoint
NA_row=table(which(is.na(B),arr.ind = TRUE)[,1])
drop=as.numeric(names(NA_row[which(NA_row>10)]))
if(length(drop)>0) B=B[-drop,]
#Transform some data
B$LogTotalEff=log10(B$OTP_MHI_OSDS_TotalEffluent+1)
B$Log_NitrogenFlux=log10(B$OTP_MHI_OSDS_NitrogenFlux+1)
B$Log_PhosphorusFlux=log10(B$OTP_MHI_OSDS_PhosphorusFlux+1)
B$Sediment_Zero=B$OTP_MHI_Sedimentation  
B$Sediment_Zero[is.na(B$Sediment_Zero)]=0
B$Log_Sediment_Zero=log10(B$Sediment_Zero+1)

# Define Driver Columns
all_x=c("Depth_m_mn",names(B)[c(which(names(B)=="DHW_v3_20151101"):ncol(B))])

#Drop Columns with lots of NA
countNA=function(x){return(length(which(is.na(x))))}
HowManyNA=apply(B[,all_x],2,countNA)
dropcols=names(which(HowManyNA>0))
#Drop any columns with NA data
if(length(dropcols)>0){all_x=all_x[-which(all_x%in%dropcols)]}

#Drop all PAR specific columns, and all 2014 columns
#dropcols2=unique(c(grep("PAR.2",all_x),grep("PAR.L",all_x),grep("2014",all_x)))
#all_x=all_x[-dropcols2]

#####################################################################Best Hypothesis Driven Model
bestmodsim2=gam(PctAffected_mean~
                  s(DHW_v3_20151101,k=5)+
                  s(PARz.2015.mn,k=5)+
                  s(SST.LT.sum_wk_rng.mn,k=5)+
                  s(PARz.LT.sum_wk_rng,k=5)+
                  s(Prop_BleachResTaxa,k=3)+
                  s(LogTotalEff,k=5)+
                  s(WaveEnergy_MN1979.2012,k=5)+
                  s(Depth_m_mn,k=3)
                ,data=B,na.action = na.fail)
sbmod2=summary(bestmodsim2)
sbmod2
#####################################################################








#PCA Analysis
pcB=prcomp(B[,all_x],scale. = T)
plot(pcB)
biplot(pcB,cex=.5)
fviz_eig(pcB)
fviz_pca_ind(pcB,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(pcB,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

pcBx=data.frame(PctAffected_mean=B$PctAffected_mean,pcB$x)
pcGAM=GamThemAll(y = "PctAffected_mean",xes = names(pcBx)[2:8],df = "pcBx",k=3)
pcD=dredge(pcGAM,rank = "BIC")
subset(pcD,delta<4)
bestmod.pc=get.models(pcD,1)[[1]]
summary(bestmod.pc)
plot(bestmod.pc,pages=1,cex=4,residuals=T)




library(dendextend)
#Check for Colinearity & Univariate Correlation
XX=1-(cor(B[,all_x])^2)
hcXX=hclust(as.dist(XX))

gcorY=rep(NA,length(all_x));names(gcorY)=all_x
for(i in 1:length(all_x)){
  gam1mod=GamThemAll(y = "PctAffected_mean",xes = c(all_x[i]),df = "B",k = 3)
  gcorY[i]=BIC(gam1mod)
}

#Cluster Variables into co-linear clusters
colin.thresh=0.6
colin_clust=cutree(tree = hcXX,h = 1-colin.thresh)

dev.off()
plot(hcXX,cex=.75)
abline(h=1-colin.thresh,col="red",lty=3)
#rect.hclust(hcXX,h=1-colin.thresh)

#build model with one variable (by univariate BIC) from each cluster
K=length(unique(colin_clust))
globalmod_vars=NULL
for(i in 1:K){
  thisclustBIC=gcorY[which(colin_clust==i)]
  globalmod_vars=c(globalmod_vars,  names(thisclustBIC[which.min(thisclustBIC)]))
}

gcorY2=gcorY[which(names(gcorY)%in%globalmod_vars)]
gcorY2=sort(gcorY2)

XX2=1-(cor(B[,names(gcorY2)])^2)
hcXX2=hclust(as.dist(XX2))
dev.off()
plot(hcXX2,cex=.75)
abline(h=1-colin.thresh,col="red",lty=3)
rect.hclust(hcXX2,h=1-colin.thresh)
length(gcorY2)
gcorY2=gcorY2[-which(names(gcorY2)%in%c("SST.LT.mn","K490.2015.mi","ShipBased_Shipping"))]

#Gam Them All, In Blocks of 6
supermod=GamThemAll(y = "PctAffected_mean",xes = names(gcorY2)[1:30],df = "B",k = 3)
superD=dredge(supermod,m.lim = c(3,6))

#Gam Them All, In Blocks of 6
Nparam=length(globalmod_vars)
Nparam
Nruns=100
NdredgeParams=10
dredgetable=NULL
dredgeDF=na.omit(B[,c("PctAffected_mean",globalmod_vars)])
for(i in 1:Nruns){
  dredgemod_vars=sample(globalmod_vars,NdredgeParams)
  dredgemod=GamThemAll(y = "PctAffected_mean",xes = dredgemod_vars,df = "dredgeDF",k=3)
  
  #Dredge Down to a Select Few
  this_dredgeset=dredge(dredgemod,m.lim=c(1,3))
  if(i==1){
    dredgetable=this_dredgeset
  }else{
    dredgetable=rbind(dredgetable,this_dredgeset)
  }
  print(paste("Run",i,"of",Nruns,"-"))
}

par(mar=c(2,2,10,2))
plot(dredgetable,cex=.75)
subset(dredgetable,delta<5)
RelMod=model.avg(dredgetable)
summary(RelMod)
Import=RelMod$importance

par(mfrow=c(1,1),mar=c(16,4,4,2))
barplot(as.numeric(Import),names.arg =names(Import),las=2,
        ylab="Relative Importance",
        main="Parameter Relative Importance, \nGlobal Model Average",
        col=c(rep("red",3),rep("gray",length(Import)-3)),cex.names =.5)
title(xlab="Hypothesized Exposure Driver",line=14)
dev.off()



## Hypthothessis Driven Model Selection ######################################
#Build Hypothesis, swap within co-linear cluster
basemod=gam(PctAffected_mean~
              s(Depth_m_mn,k=3)+
              s(DHW_v3_20151101,k=3)+
              s(PARz.2015.mn,k=3)+
              s(PARz.LT.mn,k=3)+
              s(PARz.2015.rng,k=3)+
              s(PARz.LT.rng,k=3)+
              s(SST.LT.wk_rng.mn,k=3)+
              s(WaveEnergy_MN1979.2012,k=3)+
              s(Prop_BleachResTaxa,k=3)+
              s(LogTotalEff,k=3)+
              s(Sediment_Zero,k=3)+
              s(LBSP_Urban_runoff_01,k=3)+
              s(LBSP_AgGolf_runoff_01,k=3)+
              s(OTP_MHI_Fishing_NonCommercial_ShoreBased_Total,k=3)+
              s(OTP_MHI_Fishing_NonCommercial_BoatBased_Total,k=3)+
              s(OTP_MHI_Fishing_Commercial_Total,k=3)
            ,data=B,na.action = na.fail)
summary(basemod)
testmod=lm(Btest$PctAffected_mean~predict(basemod,newdata=Btest))
summary(testmod)
vif.gam(basemod)

Dbase=dredge(basemod,rank = "BIC",m.lim = c(3,6))

bestmod=get.models(Dbase,1)[[1]]
summary(bestmod)
vis.gam(bestmod)
dev.off()
plot(bestmod,pages=1,residuals=T,cex=5)


ri=calc.relimp(bestmod,rela=T)
sort(ri$lmg)*100

AIC(basemod)
length(predict(basemod))


Btrain_i=sample(1:nrow(B),round(.9*nrow(B)),replace=F)
Btrain=B[Btrain_i,]
Btest_i=setdiff(1:nrow(B),Btrain_i)
Btest=B[Btest_i,]
bestmodsim=gam(PctAffected_mean~
                 s(DHW_v3_20151101,k=3)+
                 s(LogTotalEff,k=3)+
                 s(PARz.2015.mn,k=3)+
                 s(WaveEnergy_MN1979.2012,k=3)+
                 s(LBSP_Urban_runoff_01,k=3)+
                 s(Prop_BleachResTaxa,k=3)+
                 s(Depth_m_mn,k=3)+
                 s(OTP_MHI_Fishing_Commercial_Total,k=3),data=Btrain,na.action = na.fail)
summary(bestmodsim)
testmod=lm(Btest$PctAffected_mean~predict(basemod,newdata=Btest))
summary(testmod)


vis.gam(bestmodsim,view = c("DHW_v3_20151101","PARz.2015.mn"),theta=-45)
vis.gam(bestmodsim,view = c("LogTotalEff","Prop_BleachResTaxa"))

obsd=testset$PctAffected_mean
pred=predict(basemod,newdata=testset)
plot(obsd,pred,xlab="Observed",ylab="Predicted")
lmod=lm(obsd~pred)
summary(lmod)
abline(lmod)
################################################################################


BaseModAVG=model.avg(Dbase)
summary(BaseModAVG)
Import=BaseModAVG$importance

par(mfrow=c(1,1),mar=c(16,4,4,2))
Variables=c("Degree Heating Weeks",
            "Depth",
            "PARz 2015 mean",
            "Prop. Bleaching\nResistant Corals",
            "Commericial Fishing",
            "Total Sewage Effluent (log)",
            "Urban Runoff",
            "Wave Energy 1979-2012",
            "Sediment",
            "Non-Commericial Fishing (Shore)",
            "Non-Commericial Fishing (Boat)",
            "Ag/Golf Course Runoff",
            "SST LT weekly range")

par(mfrow=c(1,1),mar=c(16,4,4,2))
barplot(as.numeric(Import),names.arg=Variables,las=2,ylab="Relative Importance",
        main="Parameter Relative Importance, \nExposure Model Average",col=c(rep("red",3),rep("gray",length(Import)-3)))
title(xlab="Hypothesized Bleaching Driver",line=14)

bestmod=Dbase





summary(globalmod)

#Dredge Down to a Select Few
dredgeset=dredge(globalmod,m.lim=c(1,5))

#Swap Away
