rm(list=ls())
# Library/Function Loads --------------------------------------------------
source("./HCBC_2015/HelperCode/gcdist.R")
source("./HCBC_2015/HelperCode/better.biplot.R")
source("./HCBC_2015/HelperCode/HCBC_2015_Functions.R")

###########################################################################################
###########################################################################################
###########################################################################################
# Read-In, Set-Up Data --------------------------------------------------

B=read.csv("./HCBC_2015/ClusteredData/HCBC_2015_1km_PARz_MUR_PAR_K490_OTP_WE_DHW_pBRT.csv")
B$DepthBin_5m=factor(B$DepthBin_5m,levels=rev(unique(B$DepthBin_5m)))
head(B)
names(B)

#Pull EDS data
EDS=read.csv("./HCBC_2015/SatData/HCBC2015_EDS_Out_2021-03-03.csv")
# DHW_v3_20151101_sc+
#   PAR.2015.mn_sc+
#   K490.2015.mn_sc+
#   SST.LT.sum_wk_rng.mn_sc+
#   Prop_BleachResTaxa_sc+
#   Sqrt_Total_Eff_sc+
#   Sqrt_AgGolf_sc+
#   Sqrt_Urban_sc+
#   Sqrt_Wave_Energy_sc+
#   Sqrt_TourRec_sc+
#   Depth_m_mn_sc)

EDSs=EDS[,c(names(EDS)[1:11],
            "DHW.MeanMax_Degree_Heating_Weeks_MO03",
            "DHW.MeanMax_Degree_Heating_Weeks_YR10YR01",
            "mean_PAR_MODIS_Daily_MO03",
            "mean_kdPAR_VIIRS_Weekly_MO03",
            "mean_weekly_range_SST_CRW_Daily_ALLB4",
            "mean_biweekly_range_SST_CRW_Daily_ALLB4",
            "mean_monthly_range_SST_CRW_Daily_ALLB4")]

Bs=merge(B,EDSs)
Bs[which(Bs==-9991,arr.ind = T)]=NA
B=Bs

#Remove shitty NA datapoint
NA_row=table(which(is.na(B),arr.ind = TRUE)[,1])
drop=as.numeric(names(NA_row[which(NA_row>10)]))
if(length(drop)>0) B=B[-drop,]

#Transform some data

#Response # More normal as simple untransformed Proportion
asinTransform <- function(p) { asin(sqrt(p)) }
logitTransform <- function(p) { log(p/(1-p)) }
B$PropAffected_mean=B$PctAffected_mean/100
B$PropAffected_mean_arcs=asinTransform(B$PropAffected_mean)
B$PropAffected_mean_logit=logitTransform(B$PropAffected_mean+.03)
B$PropAffected_mean[is.infinite(B$PropAffected_mean_logit)]
B$PropAffected_mean_sqrt=sqrt(B$PropAffected_mean)
#explore
# par(mfrow=c(2,2))
# statA=B$PropAffected_mean
# statnameA="PropAffected_mean"
# hist(statA,50)
# qqnorm(statA,main = statnameA);qqline(statA)
# shapiro.test(statA)
# 
# hist(statB,50)
# statB=B$PropAffected_mean_sqrt
# statnameB="PropAffected_mean_sqrt"
# qqnorm(statB,main = statnameB);qqline(statB)
# shapiro.test(statB)
# 

# Predictors # having trouble getting any transform of TE to be close to normal. Log10 x+1 is simple and better than un
#OSDS
off=min(B$OTP_MHI_OSDS_TotalEffluent[B$OTP_MHI_OSDS_TotalEffluent>0])/2
B$Log_Total_Eff=log10(B$OTP_MHI_OSDS_TotalEffluent+off)
B$Sqrt_Total_Eff=sqrt(B$OTP_MHI_OSDS_TotalEffluent)
off=min(B$OTP_MHI_OSDS_NitrogenFlux[B$OTP_MHI_OSDS_NitrogenFlux>0])/2
B$Log_NitrogenFlux=log10(B$OTP_MHI_OSDS_NitrogenFlux+off)
B$Sqrt_NitrogenFlux=sqrt(B$OTP_MHI_OSDS_NitrogenFlux)
off=min(B$OTP_MHI_OSDS_PhosphorusFlux[B$OTP_MHI_OSDS_PhosphorusFlux>0])/2
B$Log_PhosphorusFlux=log10(B$OTP_MHI_OSDS_PhosphorusFlux+off)
B$Sqrt_PhosphorusFlux=sqrt(B$OTP_MHI_OSDS_PhosphorusFlux)

#OTP
B$Sediment_Zero=B$OTP_MHI_Sedimentation  
B$Sediment_Zero[is.na(B$Sediment_Zero)]=0
off=min(B$Sediment_Zero[B$Sediment_Zero>0])/2
B$Log_Sediment_Zero=log10(B$Sediment_Zero+off)
off=min(B$LBSP_AgGolf_runoff_01[B$LBSP_AgGolf_runoff_01>0])/2
B$Log_AgGolf=log10(B$LBSP_AgGolf_runoff_01+off)
off=min(B$LBSP_Urban_runoff_01[B$LBSP_Urban_runoff_01>0])/2
B$Log_Urban=log10(B$LBSP_Urban_runoff_01+off)
off=min(B$TourRec_DirectHuman_10yrAvgPUD[B$TourRec_DirectHuman_10yrAvgPUD>0])/2
B$Log_TourRec=log10(B$TourRec_DirectHuman_10yrAvgPUD+off)
off=min(B$OTP_MHI_NewDevelopment[B$OTP_MHI_NewDevelopment>0])/2
B$Log_NewDevelopment=log10(B$OTP_MHI_NewDevelopment+off)

B$Sqrt_Sediment_Zero=sqrt(B$Sediment_Zero)
B$Sqrt_AgGolf=sqrt(B$LBSP_AgGolf_runoff_01)
B$Sqrt_Urban=sqrt(B$LBSP_Urban_runoff_01)
B$Sqrt_TourRec=sqrt(B$TourRec_DirectHuman_10yrAvgPUD)
B$Sqrt_NewDevelopment=sqrt(B$OTP_MHI_NewDevelopment)

#Wave Energy
off=min(B$WaveEnergy_MN1979.2012[B$WaveEnergy_MN1979.2012>0])/2
B$Log_Wave_Energy=log10(B$WaveEnergy_MN1979.2012)
B$Sqrt_Wave_Energy=sqrt(B$WaveEnergy_MN1979.2012)

#PARz
# calculate PARz = PARs * exp(-Kpar*Z) 
# PARs = Surface PAR
# Kpar = attenuation coefficient for PAR
# Z = depth in meters
# PARz = PAR at depth
# exp(x) = e^x ; exponential function in R
B$mean_PARz_MODIS_Daily_MO03 <- B$mean_PAR_MODIS_Daily_MO03*exp(-B$mean_kdPAR_VIIRS_Weekly_MO03*B$Depth_m_mn)


# Transform Test ----------------------------------------------------------
# B$Sqrt_Total_Eff=sqrt(B$OTP_MHI_OSDS_TotalEffluent)
# bcTE=boxcoxfit(B$OTP_MHI_OSDS_TotalEffluent+off)
# B$BC_Total_Eff=boxcoxTransform(B$OTP_MHI_OSDS_TotalEffluent+off,lambda = bcTE$lambda)
# B$BC2_Total_Eff=B$BC_Total_Eff^2.5
# 
# par(mfrow=c(5,2))
# hist(B$OTP_MHI_OSDS_TotalEffluent,50)
# qqnorm(B$OTP_MHI_OSDS_TotalEffluent);qqline(B$OTP_MHI_OSDS_TotalEffluent)
# shapiro.test(B$OTP_MHI_OSDS_TotalEffluent)
# hist(B$Sqrt_Total_Eff,50)
# qqnorm(B$Sqrt_Total_Eff);qqline(B$Sqrt_Total_Eff)
# shapiro.test(B$Sqrt_Total_Eff)
# hist(B$Log_Total_Eff,50)
# qqnorm(B$Log_Total_Eff);qqline(B$Log_Total_Eff)
# shapiro.test(B$Log_Total_Eff)
# hist(B$BC_Total_Eff,50)
# qqnorm(B$BC_Total_Eff);qqline(B$BC_Total_Eff)
# shapiro.test(B$BC_Total_Eff)
# hist(B$BC2_Total_Eff,50)
# qqnorm(B$BC2_Total_Eff);qqline(B$BC2_Total_Eff)
# shapiro.test(B$BC2_Total_Eff)


# Set-Up Correlation ------------------------------------------------------
# Define Driver Columns
not_x=c(names(B)[1:(which(names(B)=="DHW_v2_20151101")-1)],"rownum","Date_mn","val")
all_x=c("Depth_m_mn",setdiff(names(B),not_x))#,names(B)[c(which(names(B)=="DHW_v3_20151101"):ncol(B))])

#Drop Columns with greater than 5 NA values
countNA=function(x){return(length(which(is.na(x))))}
HowManyNA=apply(B[,all_x],2,countNA)
dropcols=names(which(HowManyNA>5))
#Drop any columns with NA data
if(length(dropcols)>0){all_x=all_x[-which(all_x%in%dropcols)]}

#Drop any NA rows
drop_row=sort(unique(which(is.na(B[,all_x]),arr.ind = T)[,1]))
if(length(drop_row)>0){B=B[-drop_row,]}



#Drop all PAR specific columns, and all 2014 columns
#dropcols2=unique(c(grep("PAR.2",all_x),grep("PAR.L",all_x),grep("2014",all_x)))
#all_x=all_x[-dropcols2]

#Drop all PAR specific columns, and all 2014 columns
#dropcols2=unique(c(grep("PAR.2",all_x),grep("PAR.L",all_x),grep("2014",all_x)))
#all_x=all_x[-dropcols2]
#InvSD=1/B$PctAffected_sd
#InvSD[B$PctAffected_N==1]=min(InvSD,na.rm=T)
#InvSD[B$PctAffected_N==2]=2*min(InvSD,na.rm=T)

#Define Correlation Weight as Inverse Standard Error of cluster size,
#Outliers set to q05% of N>=3 and q05% of N>=3
# For N >= 3; N=1 is q05% of N>=3; N=2 is 2x q05% of N>=3; 
CorW=1/(B$PctAffected_sd/sqrt(B$PctAffected_N))
CorW_05=quantile(CorW[-is.na(CorW)],.05,na.rm=T)
CorW_95=quantile(CorW[-is.na(CorW)],.95,na.rm=T)
CorW[is.na(CorW)]=CorW_05
CorW[CorW>CorW_95]=CorW_95
CorW=CorW/CorW_95
hist(CorW,50)
B$CorW=CorW

# Scale X Values  --------------------------------------------------
for(xi in 1:length(all_x)){
  eval(parse(text=paste0("B$",all_x[xi],"_sc=as.vector(scale(B$",all_x[xi],"))")))
}
all_xsc=paste0(all_x,"_sc")
#AllLSBP
B$Log_LBSP_sc=B$Log_AgGolf_sc+B$Log_Urban_sc
B$Sqrt_LBSP_sc=B$Sqrt_AgGolf_sc+B$Sqrt_Urban_sc


# Run Full Model  --------------------------------------------------
#####################################################################Full Hypothesis Driven Model - LM
hypmod_NoInt=lm(PctAffected_mean~
                  DHW.MeanMax_Degree_Heating_Weeks_MO03_sc+
                  DHW.MeanMax_Degree_Heating_Weeks_YR10YR01_sc+
                  mean_PAR_MODIS_Daily_MO03_sc+
                  mean_kdPAR_VIIRS_Weekly_MO03_sc+
                  mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
                  Prop_BleachResTaxa_sc+
                  Sqrt_Total_Eff_sc+
                  Sqrt_AgGolf_sc+
                  Sqrt_Urban_sc+
                  Sqrt_Wave_Energy_sc+
                  Sqrt_TourRec_sc+
                  Depth_m_mn_sc,
                data=B,
                na.action = na.fail,
                weights = CorW)

sort(vif(hypmod_NoInt),decreasing = T)
modcol=names(vif(hypmod_NoInt))
pairs(B[,c("PctAffected_mean",modcol)],upper.panel=panel.cor)
# 
# hypmod_NoInt_PARZ=lm(PctAffected_mean~
#                   DHW.MeanMax_Degree_Heating_Weeks_MO03_sc+
#                   DHW.MeanMax_Degree_Heating_Weeks_YR10YR01_sc+
#                   mean_PARz_MODIS_Daily_MO03_sc+
#                   mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
#                   Prop_BleachResTaxa_sc+
#                   Sqrt_Total_Eff_sc+
#                   Sqrt_AgGolf_sc+
#                   Sqrt_Urban_sc+
#                   Sqrt_Wave_Energy_sc+
#                   Sqrt_TourRec_sc,
#                 data=B,
#                 na.action = na.fail,
#                 weights = CorW)
# 
# sort(vif(hypmod_NoInt_PARZ))
# modcol=names(vif(hypmod_NoInt_PARZ))
# pairs(B[,c("PctAffected_mean",modcol)],lower.panel=panel.cor)



# Surface PAR -------------------------------------------------------------
hypmod=lm(PctAffected_mean~
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc+
            DHW.MeanMax_Degree_Heating_Weeks_YR10YR01_sc+
            mean_PAR_MODIS_Daily_MO03_sc+
            mean_kdPAR_VIIRS_Weekly_MO03_sc+
            mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
            Prop_BleachResTaxa_sc+
            Sqrt_Total_Eff_sc+
            Sqrt_AgGolf_sc+
            Sqrt_Urban_sc+
            Sqrt_Wave_Energy_sc+
            Sqrt_TourRec_sc+
            Depth_m_mn_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:DHW.MeanMax_Degree_Heating_Weeks_YR10YR01_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:mean_PAR_MODIS_Daily_MO03_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:mean_kdPAR_VIIRS_Weekly_MO03_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Prop_BleachResTaxa_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_Total_Eff_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_AgGolf_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_Urban_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_Wave_Energy_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_TourRec_sc+
            DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Depth_m_mn_sc,
          data=B,
          na.action = na.fail,
          weights = CorW)

shypmod=summary(hypmod)
shypmod

StepHyp=stepAIC(hypmod,k = log(nrow(B)))
#StepHypMod=update(StepHyp)

SumStepHyp=summary(StepHyp)
SumStepHyp

#vif check
modints=attr(SumStepHyp$terms,"term.labels")
modints=modints[grep(":",modints)]
StepHypNoInt=dropNterm(mod = StepHyp,drop_terms = modints,drop_ints = T)

vif_final=vif(StepHypNoInt)
vif_final[order(vif_final,decreasing = T)]

# On we go... -------------------------------------------------------------
CoefStepHyp=as.data.frame(SumStepHyp$coefficients)
CoefStepHyp=CoefStepHyp[-1,]
CoefStepHyp[order(CoefStepHyp$Estimate,decreasing = T),]
CoefStepHyp$Name=row.names(CoefStepHyp)
CoefStepHyp$Name=factor(CoefStepHyp$Name,levels=CoefStepHyp[order(abs(CoefStepHyp$Estimate)-abs(CoefStepHyp$`Std. Error`),decreasing = F),"Name"])
CoefStepHyp$SigClass=cut(CoefStepHyp$`Pr(>|t|)`,breaks=c(0,.0001,.001,.01,.05,Inf),labels=c("p<0.0001","p<0.001","p<0.01","p<0.05","NS"))


# PARz Model --------------------------------------------------------------
# hypmodPz=lm(PctAffected_mean~(
#   DHW.MeanMax_Degree_Heating_Weeks_MO03_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_YR10YR01_sc+
#     mean_PARz_MODIS_Daily_MO03_sc+
#     mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
#     Prop_BleachResTaxa_sc+
#     Sqrt_Total_Eff_sc+
#     Sqrt_AgGolf_sc+
#     Sqrt_Urban_sc+
#     Sqrt_Wave_Energy_sc+
#     Sqrt_TourRec_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:mean_PARz_MODIS_Daily_MO03_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Prop_BleachResTaxa_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_Total_Eff_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_AgGolf_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_Urban_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_Wave_Energy_sc+
#     DHW.MeanMax_Degree_Heating_Weeks_MO03_sc:Sqrt_TourRec_sc+
#     Prop_BleachResTaxa_sc:mean_PARz_MODIS_Daily_MO03_sc+
#     Prop_BleachResTaxa_sc:mean_weekly_range_SST_CRW_Daily_ALLB4_sc+
#     Prop_BleachResTaxa_sc:Prop_BleachResTaxa_sc+
#     Prop_BleachResTaxa_sc:Sqrt_Total_Eff_sc+
#     Prop_BleachResTaxa_sc:Sqrt_AgGolf_sc+
#     Prop_BleachResTaxa_sc:Sqrt_Urban_sc+
#     Prop_BleachResTaxa_sc:Sqrt_Wave_Energy_sc+
#     Prop_BleachResTaxa_sc:Sqrt_TourRec_sc),
#   data=B,
#   na.action = na.fail,
#   weights = CorW)
# 
# shypmodPz=summary(hypmodPz)
# shypmodPz
# 
# StepHypPz=stepAIC(hypmodPz,k = log(nrow(B)))

#SumStepHypPz=summary(StepHypPz)
#SumStepHypPz

#SumStepHyp
#SumStepHypPz

