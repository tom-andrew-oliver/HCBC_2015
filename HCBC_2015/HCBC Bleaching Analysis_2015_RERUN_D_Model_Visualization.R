rm(list=ls())
# Library/Function Loads --------------------------------------------------
source("./HCBC_2015/HelperCode/gcdist.R")
source("./HCBC_2015/HelperCode/better.biplot.R")
source("./HCBC_2015/HelperCode/HCBC_2015_Functions.R")
source("./HCBC_2015/HCBC Bleaching Analysis_2015_RERUN_C_Model Data.R")
#Returns StepHyp as best fit model, SumStepHyp is summary, CorW are correlation weights, B is data frame
SumStepHyp


#Plotnames
fullmod.terms=attr(terms(hypmod),"term.labels")
fulltermsLU=c("Acute Thermal Stress (DHW)",
              "Historic Thermal Stress (DHW)",
              "Surface Light (PAR)",
              "Light Attenuation (KdPAR)",
              "Thermal Variability (SST range)",
              "Prop. Bleaching Resistant Taxa",
              "Sewage Effluent",
              "Ag/Golf Runoff",
              "Urban Runoff",
              "Wave Energy",
              "Tourism/Recreation",
              "Depth (m)",
              "Acute Thermal Stress (DHW) : Historic Thermal Stress (DHW)",
              "Acute Thermal Stress (DHW) : Surface Light (PAR)",
              "Acute Thermal Stress (DHW) : Light Attenuation (KdPAR)",
              "Acute Thermal Stress (DHW) : Thermal Variability (SST range)",
              "Acute Thermal Stress (DHW) : Prop. Bleaching Resistant Taxa",
              "Acute Thermal Stress (DHW) : Sewage Effluent",
              "Acute Thermal Stress (DHW) : Ag/Golf Runoff",
              "Acute Thermal Stress (DHW) : Urban Runoff",
              "Acute Thermal Stress (DHW) : Wave Energy",
              "Acute Thermal Stress (DHW) : Tourism/Recreation",
              "Acute Thermal Stress (DHW) : Depth (m)",
              "Prop. Bleaching Resistant Taxa : Historic Thermal Stress (DHW)",
              "Prop. Bleaching Resistant Taxa : Surface Light (PAR)",
              "Prop. Bleaching Resistant Taxa : Light Attenuation (KdPAR)",
              "Prop. Bleaching Resistant Taxa : Thermal Variability (SST range)",
              "Prop. Bleaching Resistant Taxa : Sewage Effluent",
              "Prop. Bleaching Resistant Taxa : Ag/Golf Runoff",
              "Prop. Bleaching Resistant Taxa : Urban Runoff",
              "Prop. Bleaching Resistant Taxa : Wave Energy",
              "Prop. Bleaching Resistant Taxa : Tourism/Recreation",
              "Prop. Bleaching Resistant Taxa : Depth (m)")
names(fulltermsLU)=fullmod.terms

mod.terms=attr(terms(StepHyp),"term.labels")
termsLU=c("Acute Thermal Stress (DHW)",
          "Historic Thermal Stress (DHW)",
          "Surface Light (PAR)",
          "Light Attenuation (KdPAR)",
          "Thermal Variability (SST range)",
          "Sewage Effluent",
          "Urban Runoff",
          "Wave Energy",
          "Depth (m)",
          "Acute Thermal Stress : Light (PAR)",
          "Acute Thermal Stress : Light Attenuation",
          "Acute Thermal Stress : Sewage Effluent",
          "Acute Thermal Stress : Depth")
names(termsLU)=mod.terms

#Transforms
exvars=row.names(attr(terms(StepHyp),"factors"))[-1]
exvars
transformsLU=c("None",
               "None",
               "None",
               "None",
               "None",
               "sqrt",
               "sqrt",
               "sqrt",
               "None")
names(transformsLU)=exvars


#TimeStampFor Saving
timestamp=Sys.Date()
#scale for plot outputs
sc=2

#residuals plots
par(mfrow=c(2,2))
plot(StepHyp)


#Plot Estimate Figure
FullModEstPlot=plot_model_estimates(mod=hypmod,termLU=fulltermsLU,split_ints=T,response_name="Percent Reef Affected (Pale or Bleached) : Full Model")
FullModEstPlot
ggsave(filename = paste0("./HCBC_2015/Output_Plots/FullModelEstimatePlot_",timestamp,".jpg"),plot = FullModEstPlot,width=sc*8.5,height=sc*11/4)

BestModEstPlot=plot_model_estimates(mod=StepHyp,termLU=termsLU,split_ints=T,response_name="Percent Reef Affected (Pale or Bleached): AIC Selected Model")
BestModEstPlot
ggsave(filename = paste0("./HCBC_2015/Output_Plots/BestModelEstimatePlot_",timestamp,".jpg"),plot = BestModEstPlot,width=sc*8.5,height=sc*11/4)

bothmodest=FullModEstPlot/BestModEstPlot+plot_annotation(tag_levels = "A")
bothmodest
ggsave(filename = paste0("./HCBC_2015/Output_Plots/BothModelEstimatePlots_",timestamp,".jpg"),plot = bothmodest,width=sc*8.5,height=sc*11)


# Plot Partial Residual vs Predictor Plots
prp=partial_residual_plots(mod=StepHyp,data=B,termLU = termsLU,transLU = transformsLU)
prp
ggsave(filename = paste0("./HCBC_2015/Output_Plots/PartialResidualPlots_EachVar_",timestamp,".jpg"),plot = prp,width=sc*8.5,height=sc*11/2)

jpg(paste0("./HCBC_2015/Output_Plots/AddedVariablePlots_EachVar_",timestamp,".jpg"),width=sc*8.5,height=sc*11/2)
singform=formula(paste0("~",paste0(exvars,collapse="+")))
avPlots(model = StepHyp,terms =singform)
dev.off()

# Plot Interactions- Partial Residual Plots
pip=partial_interaction_plots(mod=StepHyp,data=B,termLU = termsLU,transLU = transformsLU,color_scale=c(25,75))
pip
ggsave(filename = paste0("./HCBC_2015/Output_Plots/PartialInteractionPlots_2Way_",timestamp,".jpg"),plot = prp,width=sc*8.5,height=sc*11/2)

#### Perturbation
#Prep
all_terms=row.names(attr(StepHyp$terms,"factors"))
dep_terms=all_terms[-1]
res_var=all_terms[1]
dep_terms_usc=substr(dep_terms,1,nchar(dep_terms)-3)
D_mod=B[,c(all_terms,dep_terms_usc)]
D_mod$Full_Y_Pred=predict(StepHyp,newdata=D_mod)

#Basic Performance
ggplot(D_mod,aes(x=PctAffected_mean,y=Full_Y_Pred))+geom_point()+stat_smooth(method="lm")

#Run Perturbation
PerturbP=1
plotlist=vector("list", length(dep_terms))
for(p_i in 1:length(dep_terms)){
  thisvar=dep_terms[p_i]
  thisvar_usc=substr(thisvar,1,nchar(thisvar)-3)
  pertY=pertDat_fun(var = thisvar_usc, var_s = thisvar, "none", D_mod, StepHyp, PerturbP)
  Pnames=paste0(c("Yp_INC_","Yp_DEC_","Yd_INC_","Yd_DEC_"),thisvar)
  names(pertY)[3:6]=Pnames
  D_mod[,Pnames]=pertY[,Pnames]
  
  #Generate Plots
  plotlist[[p_i]]=pertPlots_fun(pertY) + ggtitle(termsLU[thisvar])+xlab(termsLU[thisvar])+ylab("Pct. Affected")
  
  #Summarize
  se=function(x,na.rm=T){return(sd(x,na.rm=na.rm)/sqrt(length(x)))}
  PertSum_row=c(thisvar,termsLU[thisvar],apply(pertY,2,mean,na.rm=T)[2:6],apply(pertY,2,sd,na.rm=T)[2:6],apply(pertY,2,se,na.rm=T)[2:6],nrow(pertY))
  names(PertSum_row)=c("var","var_name","Y_mn","Yinc_mn","Ydec_mn","delYinc_mn","delYdec_mn",
                       "Y_sd","Yinc_sd","Ydec_sd","delYinc_sd","delYdec_sd",
                       "Y_se","Yinc_se","Ydec_se","delYinc_se","delYdec_se","N")
  if(p_i==1){PertSum=PertSum_row}else{PertSum=rbind(PertSum,PertSum_row)}
}
PertSum=as.data.frame(PertSum)
PertSum[,3:18]=apply(PertSum[,3:18],2,as.numeric)

#Summary Plot
sortmet=abs(PertSum$delYdec_mn)-abs(PertSum$delYdec_se)#apply(abs(PertSum[,c("delYdec_mn","delYinc_mn")]),1,mean)
PertSum$var_name=factor(PertSum$var_name,levels=PertSum$var_name[order(sortmet,decreasing = T)])
PertSumPlot=ggplot(PertSum,aes(x=var_name))+
  geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=Y_mn-Y_se,ymax=Y_mn+Y_se),fill="gray75",alpha=.1)+
  geom_hline(yintercept = PertSum$Y_mn[1])+
  #  geom_segment(aes(xend=var_name,yend=PertSum$Y_mn[1],y=Yinc_mn),color="orange",size=z)+
  #  geom_segment(aes(xend=var_name,yend=PertSum$Y_mn[1],y=Ydec_mn),color="skyblue",size=2)+
  geom_point(aes(y=Yinc_mn),fill="orange",color="black",size=3,shape=24)+
  geom_point(aes(y=Ydec_mn),fill="skyblue",color="black",size=3,shape=25)+
  scale_fill_manual(name="UP/DOWN",values = c("orange","skyblue"),labels=c("Increase","Decrease"))+
  geom_errorbar(aes(ymin=Yinc_mn-Yinc_sd/sqrt(N),ymax=Yinc_mn+Yinc_sd/sqrt(N)),color="orange",width=.25)+
  geom_errorbar(aes(ymin=Ydec_mn-Ydec_sd/sqrt(N),ymax=Ydec_mn+Ydec_sd/sqrt(N)),color="skyblue",width=.25)+
  theme_bw()+
  theme(axis.text=element_text(angle=45,hjust=1))+
  xlab("Variable")+ylab("Predicted Pct Affected")+
  ggtitle("Model Perturbation (+/- 1 SD)")
PertSumPlot
ggsave(filename = paste0("./HCBC_2015/Output_Plots/PertubationSensitivityPlot_",timestamp,".jpg"),plot = PertSumPlot,width=sc*8.5,height=sc*11/2)

#Generate Plots
pertYplot_grid=wrap_plots(plotlist,guides="collect")+plot_annotation(title="Model Perturbation (+/- 0.5 SD)")
pertYplot_grid
ggsave(filename = paste0("./HCBC_2015/Output_Plots/PertubationSensitivityGrid_",timestamp,".jpg"),plot = pertYplot_grid,width=sc*8.5,height=sc*11/2)



png(width = 1350, height = 1450, filename = "Drivers_Pert_.5.png")
grid.arrange(
  pertPlots_fun(pertDat_fun(var = "DHW.MeanMax_Degree_Heating_Weeks_MO03", var_s = "DHW.MeanMax_Degree_Heating_Weeks_MO03_sc", "none", B, StepHyp, PerturbP)) +
    ggtitle("Acute Thermal Stress") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("DHW.MeanMax_Degree_Heating_Weeks_YR10YR01", "DHW.MeanMax_Degree_Heating_Weeks_YR10YR01_sc", "none", B, StepHyp, PerturbP)) +
    ggtitle("Historical Thermal Stress") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("mean_PAR_MODIS_Daily_MO03", "mean_PAR_MODIS_Daily_MO03_sc", "none", B, StepHyp, PerturbP)) + 
    ggtitle("Surface Light (PAR)") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("mean_kdPAR_VIIRS_Weekly_MO03", "mean_kdPAR_VIIRS_Weekly_MO03_sc", "none", B, StepHyp, PerturbP)) + 
    ggtitle("Light Attenuation") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("mean_weekly_range_SST_CRW_Daily_ALLB4", "mean_weekly_range_SST_CRW_Daily_ALLB4_sc", "none", B, StepHyp, PerturbP)) + 
    ggtitle("Thermal Variability") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("Sqrt_Total_Eff", "Sqrt_Total_Eff_sc", "sqrt", B, StepHyp, PerturbP)) + 
    ggtitle("Effluent") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("Sqrt_Urban", "Sqrt_Urban_sc", "sqrt", B, StepHyp, PerturbP)) + 
    ggtitle("Urban Run-off") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("Sqrt_Wave_Energy", "Sqrt_Wave_Energy_sc", "sqrt", B, StepHyp, PerturbP)) + 
    ggtitle("Wave Energy") + theme(legend.position = "none"),
  pertPlots_fun(pertDat_fun("Depth_m_mn", "Depth_m_mn_sc", "sqrt", B, StepHyp, PerturbP)) + 
    ggtitle("Depth (m)") + theme(legend.position = c(0.76,0.1)),
  nrow = 3, left = textGrob("Predicted Percent Affected (%)", gp=gpar(fontsize=22), rot = 90), top = textGrob("SD Perturbation\n", gp=gpar(fontsize=22)))
dev.off()



ss=sensitivity_prediction_matrix_sc(mod = StepHyp,data=B,pperturb = .25)
lsort=subset(ss,VAR=="OBS")
lsort=lsort[order(lsort$Y),]
ss$SITE_ID=factor(ss$SITE_ID,levels=lsort$SITE_ID)
ss$VAR=factor(ss$VAR,levels=c("OBS","PRED",terms))

ggplot(ss,aes(x=SITE_ID,y=Y,color=DIRECTION))+geom_point()+facet_wrap("VAR")+ggtitle("25% Perturbation")+scale_x_discrete(guide=NULL)

sss=ddply(ss,.(VAR,DIRECTION),summarize,mn=mean(Y),N=length(Y),sd=sd(Y),se=sd/sqrt(N))
sss

