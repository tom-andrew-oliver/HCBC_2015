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
hypmod_NoInt=lm(PctAffected_mean~(
  DHW_v3_20151101_sc+
    PAR.2015.mn_sc+
    K490.2015.mn_sc+
    SST.LT.sum_wk_rng.mn_sc+
    Prop_BleachResTaxa_sc+
    Sqrt_Total_Eff_sc+
    Sqrt_AgGolf_sc+
    Sqrt_Urban_sc+
    Sqrt_Wave_Energy_sc+
    Sqrt_TourRec_sc+
    Depth_m_mn_sc),
  data=B,
  na.action = na.fail,
  weights = CorW)

sort(vif(hypmod_NoInt))
modcol=names(vif(hypmod_NoInt))
pairs(B[,c("PctAffected_mean",modcol)],lower.panel=panel.cor)


hypmod=lm(PctAffected_mean~(
  DHW_v3_20151101_sc+
    PAR.2015.mn_sc+
    K490.2015.mn_sc+
    SST.LT.sum_wk_rng.mn_sc+
    Prop_BleachResTaxa_sc+
    Sqrt_Total_Eff_sc+
    Sqrt_AgGolf_sc+
    Sqrt_Urban_sc+
    Sqrt_Wave_Energy_sc+
    Sqrt_TourRec_sc+
    Depth_m_mn_sc+
    DHW_v3_20151101_sc:PAR.2015.mn_sc+
    DHW_v3_20151101_sc:K490.2015.mn_sc+
    DHW_v3_20151101_sc:SST.LT.sum_wk_rng.mn_sc+
    DHW_v3_20151101_sc:Prop_BleachResTaxa_sc+
    DHW_v3_20151101_sc:Sqrt_Total_Eff_sc+
    DHW_v3_20151101_sc:Sqrt_AgGolf_sc+
    DHW_v3_20151101_sc:Sqrt_Urban_sc+
    DHW_v3_20151101_sc:Sqrt_Wave_Energy_sc+
    DHW_v3_20151101_sc:Sqrt_TourRec_sc+
    DHW_v3_20151101_sc:Depth_m_mn_sc+
    Prop_BleachResTaxa_sc:PAR.2015.mn_sc+
    Prop_BleachResTaxa_sc:K490.2015.mn_sc+
    Prop_BleachResTaxa_sc:SST.LT.sum_wk_rng.mn_sc+
    Prop_BleachResTaxa_sc:Prop_BleachResTaxa_sc+
    Prop_BleachResTaxa_sc:Sqrt_Total_Eff_sc+
    Prop_BleachResTaxa_sc:Sqrt_AgGolf_sc+
    Prop_BleachResTaxa_sc:Sqrt_Urban_sc+
    Prop_BleachResTaxa_sc:Sqrt_Wave_Energy_sc+
    Prop_BleachResTaxa_sc:Sqrt_TourRec_sc+
    Prop_BleachResTaxa_sc:Depth_m_mn_sc),
  data=B,
  na.action = na.fail,
  weights = CorW)

shypmod=summary(hypmod)
shypmod

StepHyp=stepAIC(hypmod,k = log(nrow(B)))
#StepHypMod=update(StepHyp)

SumStepHyp=summary(StepHyp)
SumStepHyp

CoefStepHyp=as.data.frame(SumStepHyp$coefficients)
CoefStepHyp=CoefStepHyp[-1,]
CoefStepHyp[order(CoefStepHyp$Estimate,decreasing = T),]
CoefStepHyp$Name=row.names(CoefStepHyp)
CoefStepHyp$Name=factor(CoefStepHyp$Name,levels=CoefStepHyp[order(abs(CoefStepHyp$Estimate)-abs(CoefStepHyp$`Std. Error`),decreasing = F),"Name"])
CoefStepHyp$SigClass=cut(CoefStepHyp$`Pr(>|t|)`,breaks=c(0,.0001,.001,.01,.05,Inf),labels=c("p<0.0001","p<0.001","p<0.01","p<0.05","NS"))

SumStepHyp
ggplot(CoefStepHyp,aes(x=Estimate,xmin=Estimate-`Std. Error`,xmax=Estimate+`Std. Error`,y=Name,color=SigClass))+
  geom_point()+
  geom_errorbarh()+
  geom_vline(xintercept = 0)+theme_bw()

StepHypNoInt=update(StepHyp,".~.-DHW_v3_20151101_sc:PAR.2015.mn_sc -DHW_v3_20151101_sc:K490.2015.mn_sc   -DHW_v3_20151101_sc:Log_Total_Eff_sc  -DHW_v3_20151101_sc:Depth_m_mn_sc      -Prop_BleachResTaxa_sc:Depth_m_mn_sc ")
summary(StepHypNoInt)

sort(vif(StepHypNoInt))
pairs(B[,c("PctAffected_mean",names(vif(StepHypNoInt)))],lower.panel=panel.cor)


sensitivity_prediction_matrix_sc=function(mod,data,pperturb=.1){
  terms=attributes(mod$terms)$term.labels[attributes(mod$terms)$order==1]
  response=as.character(mod$call$formula)[2]
  org_data=data[,c(terms,response)]
  
  #Get unscaled data
  uc_terms=substr(terms,1,nchar(terms)-3)
  uc_data=data[,c(uc_terms)]

  #set up for scaling after perterbation
  mn_data=apply(uc_data,2,mean,na.rm=T)
  sd_data=apply(uc_data,2,sd,na.rm=T)
  mn_mat=matrix(rep(mn_data,nrow(uc_data)),ncol = ncol(uc_data),nrow=nrow(uc_data),byrow = T)
  sd_mat=matrix(rep(sd_data,nrow(uc_data)),ncol = ncol(uc_data),nrow=nrow(uc_data),byrow = T)
  
  out_data=data.frame(SITE_ID=rep(1:nrow(org_data),2),Y=c(org_data[,response],predict(mod,new_data=org_data)),
                      VAR=c(rep("OBS",nrow(org_data)),rep("PRED",nrow(org_data))),
                      DIRECTION="NONE")
  for(i in 1:length(uc_terms)){
    #Nudge Up
    mod_data=uc_data
    nudge_by=mean(mod_data[,uc_terms[i]],na.rm=T)*pperturb
    mod_data[,uc_terms[i]]=mod_data[,uc_terms[i]]+nudge_by
    mod_data_sc=(mod_data-mn_mat)/sd_mat
    names(mod_data_sc)=paste0(names(mod_data),"_sc")
    append_data=data.frame(SITE_ID=1:nrow(org_data),Y=predict(mod,newdata=mod_data_sc),VAR=terms[i],DIRECTION="UP")
    out_data=rbind(out_data,append_data)
    #Nudge Down
    mod_data=uc_data
    nudge_by=mean(mod_data[,uc_terms[i]],na.rm=T)*pperturb
    mod_data[,uc_terms[i]]=mod_data[,uc_terms[i]]-nudge_by
    mod_data_sc=(mod_data-mn_mat)/sd_mat
    names(mod_data_sc)=paste0(names(mod_data),"_sc")
    append_data=data.frame(SITE_ID=1:nrow(org_data),Y=predict(mod,newdata=mod_data_sc),VAR=terms[i],DIRECTION="DOWN")
    out_data=rbind(out_data,append_data)
  }
  return(out_data)
}

ss=sensitivity_prediction_matrix_sc(mod = StepHyp,data=B,pperturb = .10)
lsort=subset(ss,VAR=="OBS")
lsort=lsort[order(lsort$Y),]
ss$SITE_ID=factor(ss$SITE_ID,levels=lsort$SITE_ID)
ss$VAR=factor(ss$VAR,levels=c("OBS","PRED",terms))

ggplot(ss,aes(x=SITE_ID,y=Y,color=DIRECTION))+geom_point()+facet_wrap("VAR")+ggtitle("25% Perturbation")+scale_x_discrete(guide=NULL)

sss=ddply(ss,.(VAR,DIRECTION),summarize,mn=mean(Y),N=length(Y),sd=sd(Y),se=sd/sqrt(N))
sss

yoff=subset(sss,VAR=="PRED")$mn
ggplot(sss,aes(y=VAR,x=mn-yoff,xmin=mn-se-yoff,xmax=mn+se-yoff,color=DIRECTION))+
  geom_point()+
  geom_errorbar()+
  theme_bw()+
  geom_vline(xintercept = 0)+
  scale_x_reverse()+
  theme(axis.text.x = element_text(angle=45,hjust=1))+ggtitle("Sensitivity Plot, +/- 10%")

modminus=update(StepHyp,.~-DHW_v3_20151101_sc)


bestmodsim2=gam(PctAffected_mean~
                  s(DHW_v3_20151101_sc,k=3)+
                  s(PAR.2015.mn_sc,k=3)+
                 s(K490.2015.mn_sc,k=3)+
                  s(SST.LT.sum_wk_rng.mn_sc,k=3)+
                 # s(Prop_BleachResTaxa_sc,k=3)+
                  #s(Log_Total_Eff_sc,k=3)+
                  #s(Log_AgGolf_sc,k=3)+
                  s(Log_Urban_sc,k=3)+
                  s(Log_TourRec_sc,k=3)+
                  #s(Log_Wave_Energy_sc,k=3)+
                  s(Depth_m_mn_sc,k=3)
                ,data=B,na.action = na.fail,weights = CorW)
sbmod2=summary(bestmodsim2)
sbmod2
plot(bestmodsim2,pages=1,residuals=T)



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
vif.gam(basemod)head(b0)

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
