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
library(ggpubr)
library(patchwork)
library(metR)
library(scales)
library(car)
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


plot_model_estimates=function(mod,termLU=NULL,split_ints=T,response_name=NULL){
  sumod=summary(mod)
  plot_df=data.frame(terms=attr(mod$terms,"term.labels"),
                     estimates=sumod$coefficients[-1,1],
                     SE=sumod$coefficients[-1,2],
                     P_val=sumod$coefficients[-1,4],stringsAsFactors = F,
                     isInt=F)
  plot_df$isInt[grep(":",plot_df$terms)]=T
  if(!is.null(termLU)){plot_df$terms=termLU[plot_df$terms]}
  plot_df$pval_bin=cut(plot_df$P_val,breaks=c(0,.001,.01,.05,1),labels=c("p < 0.001","p < 0.01","p < 0.05","NS"))
  if(split_ints){termsort=plot_df$terms[order(!plot_df$isInt,abs(plot_df$estimates))]
    }else{termsort=plot_df$terms[order(abs(plot_df$estimates))]}
  plot_df$terms=factor(plot_df$terms,levels=termsort)
  outplot=ggplot(plot_df,aes(x=estimates,y=terms,color=pval_bin))+
    geom_point()+
    geom_errorbarh(aes(xmin=estimates-SE,xmax=estimates+SE))+
    geom_vline(xintercept = 0)+
    scale_color_discrete(name="Significance",drop=FALSE)+
    ylab("Model Terms")+
    xlab("Model Slope Estimates (scaled)")+
    ggtitle(paste0("Multiple Regression Model Estimates: ",response_name))+
    theme_pubr()+theme(legend.position="bottom")
  return(outplot)
}

drop1term=function(mod,drop_term,drop_ints=T){
  terms=attr(mod$terms,"term.labels")
  if (drop_ints){drop_terms=terms[grep(drop_term,terms)]
  }else{drop_terms=drop_term}
  for(i in 1:length(drop_terms)){
    eval(parse(text=paste0("mod=update(mod,. ~ . - ",drop_terms[i],")")))
  }
  return(mod)
}


dropNterm=function(mod,drop_terms,drop_ints=T){
  terms=attr(mod$terms,"term.labels")
  #Collect Model terms to drop
  drop_all=NULL
  for(i in 1:length(drop_terms)){
    if (drop_ints){drop_all=c(drop_all,terms[grep(drop_terms[i],terms)])
    }else{drop_all=drop_terms}
  }
  #Update Model
  for(i in 1:length(drop_all)){
    eval(parse(text=paste0("mod=update(mod,. ~ . - ",drop_all[i],")")))
  }
  return(mod)
}

partial_residual_plots=function(mod,data,termLU=NULL,transLU=NULL){
  sumod=summary(mod)
  terms=rownames(attr(terms(mod),"factors"))[-1]
  Nt=length(terms)
  plotlist=list(Nt)
  meandata=data.frame(matrix(rep(apply(data[,terms],2,mean,na.rm=T),100),nrow=100,byrow=T))
  names(meandata)=terms
  
  P_val=sumod$coefficients[-1,4]
  P_val_bin=cut(P_val,breaks=c(0,.001,.01,.05,1),labels=c("p < 0.001","p < 0.01","p < 0.05","NS"))
  
  
  for(i in 1:Nt){
    thisvar=terms[i]
    if(!is.null(termLU)){thisplotvar=termLU[thisvar]}else{thisplotvar=thisvar}

    #Generate Point Data
    thisdata=meandata
    thisdata[,thisvar]=seq(min(data[,thisvar],na.rm=T),max(data[,thisvar],na.rm=T),length.out = 100)
    redmod=drop1term(mod,drop_term = terms[i])
    thisx=data[,thisvar]
    thisx_usc=data[,substr(thisvar,1,nchar(thisvar)-3)]
    scmn=mean(thisx_usc,na.rm=T)
    scsd=sd(thisx_usc,na.rm=T)
    resY=residuals(redmod)
    point_df=data.frame(x_sc=thisx,y=resY,x_usc=thisx_usc,scmn=scmn,scsd=scsd)
    if(!is.null(transLU)&transLU[thisvar]!="None"){
      if(transLU[thisvar]=="sqrt"){
        point_df$x_utr=point_df$x_usc^2
        xlabels=pretty(point_df$x_utr)
        xbreaks=(sqrt(xlabels)-scmn)/scsd
      }else
        if(transLU[thisvar]=="log"){
          point_df$x_utr=exp(point_df$x_usc)
          xlabels=pretty(point_df$x_utr)
          xbreaks=(log(xlabels)-scmn)/scsd
        }
    }else{
      point_df$x_utr=point_df$x_usc
      xlabels=pretty(point_df$x_utr)
      xbreaks=(xlabels-scmn)/scsd
    }
    
    #Line
    predY=predict(mod,newdata=thisdata,se.fit=T)
    mod_int=mod$coefficients[1]
    line_df=data.frame(x_sc=thisdata[,thisvar],
                       y=predY$fit-mod_int,
                       ymin=predY$fit-predY$se.fit-mod_int,
                       ymax=predY$fit+predY$se.fit-mod_int,
                       P_val_bin=P_val_bin[i])
    
    plotlist[[i]]=ggplot()+
      geom_point(data=point_df,aes(x=x_sc,y=y,size=mod$weights))+
      geom_ribbon(data=line_df,aes(x=x_sc,ymin=ymin,ymax=ymax),color="gray75",alpha=.25)+
      geom_line(data=line_df,aes(x=x_sc,y=y,color=P_val_bin),size=1)+
      ggtitle(thisplotvar)+
      scale_size(guide="none")+
      scale_color_discrete(name="Significance",drop=FALSE)+
      ylab("Partial Residual")+
      xlab(thisplotvar)+
      scale_x_continuous(breaks = xbreaks,labels=xlabels)+
      theme_pubr()+
      theme(legend.position = "bottom",axis.text=element_text(size=8))
      
    
  }
  outp=wrap_plots(plotlist,guides="collect")& theme(legend.position = 'bottom')
  #outlist=list(outp,plotlist)
  return(outp)
}


partial_interaction_plots=function(mod,data,termLU=NULL,transLU=NULL,color_scale=NULL){
  sumod=summary(mod)
  terms=rownames(attr(terms(mod),"factors"))[-1]
  ints_i=grep(":",attr(terms(mod),"term.labels"))
  ints=attr(terms(mod),"term.labels")[ints_i]
  Nt=length(ints)
  plotlist=list(Nt)
  meandata_row=apply(data[,terms],2,mean,na.rm=T)
  names(meandata_row)=terms
  
#  P_val=sumod$coefficients[ints_i+1,4]
#  P_val_bin=cut(P_val,breaks=c(0,.001,.01,.05,1),labels=c("p < 0.001","p < 0.01","p < 0.05","NS"))
  for(i in 1:Nt){
    thesevar=unlist(strsplit(ints[i],":"))
    if(!is.null(termLU)){theseplotvar=termLU[thesevar]}else{theseplotvar=thesevar}
    
    #Generate Data Frame For Plotting
    Nvar=length(thesevar)
    if(Nvar>2){print(paste0("Only plotting 2-Way interactions. Skipping this ",Nvar," variable interaction."));next()}
    var1=seq(min(data[,thesevar[1]],na.rm=T),max(data[,thesevar[1]],na.rm=T),length.out = 100)
    var2=seq(min(data[,thesevar[2]],na.rm=T),max(data[,thesevar[2]],na.rm=T),length.out = 100)
    vars=expand.grid(var1,var2)
    thisdata=data.frame(matrix(rep(meandata_row,nrow(vars)),nrow=nrow(vars),byrow = T))
    names(thisdata)=terms
    thisdata[,thesevar[1]]=vars[1]
    thisdata[,thesevar[2]]=vars[2]
    
    #Update Reduced Model
    #Show residuals from model without just this interaction term
    redmod=drop1term(mod,drop_term = ints[i],drop_ints = F)
    #Show residuals from model without either term, or any other interactions it's in
    #redmod=dropNterm(mod,drop_terms = thesevar,drop_ints = T)
    
    #Build out raw, trans, scaled in X
    thisx=data[,thesevar[1]]
    thisx_usc=data[,substr(thesevar[1],1,nchar(thesevar[1])-3)]
    xscmn=mean(thisx_usc,na.rm=T)
    xscsd=sd(thisx_usc,na.rm=T)
    #Build out raw, trans, scaled in Y
    thisy=data[,thesevar[2]]
    thisy_usc=data[,substr(thesevar[2],1,nchar(thesevar[2])-3)]
    yscmn=mean(thisy_usc,na.rm=T)
    yscsd=sd(thisy_usc,na.rm=T)
    
    #Get Reduced Model Residuals    
    resZ=residuals(redmod)
    
    #Start Building Point DataFrame
    point_df=data.frame(x_sc=thisx,x_usc=thisx_usc,
                        y_sc=thisy,y_usc=thisy_usc,
                        z_res=resZ,
                        xscmn=xscmn,xscsd=xscsd,
                        yscmn=yscmn,yscsd=yscsd)
    
    #Get Axis Labels -X
    if(!is.null(transLU)&transLU[thesevar[1]]!="None"){
      if(transLU[thesevar[1]]=="sqrt"){
        point_df$x_utr=point_df$x_usc^2
        xlabels=pretty(point_df$x_utr)
        xbreaks=(sqrt(xlabels)-xscmn)/xscsd
      }else
        if(transLU[thesevar[1]]=="log"){
          point_df$x_utr=exp(point_df$x_usc)
          xlabels=pretty(point_df$x_utr)
          xbreaks=(log(xlabels)-xscmn)/xscsd
        }
    }else{
      point_df$x_utr=point_df$x_usc
      xlabels=pretty(point_df$x_utr)
      xbreaks=(xlabels-xscmn)/xscsd
    }
    
    #Get Axis Labels -Y
    if(!is.null(transLU)&transLU[thesevar[2]]!="None"){
      if(transLU[thesevar[2]]=="sqrt"){
        point_df$y_utr=point_df$y_usc^2
        ylabels=pretty(point_df$y_utr)
        ybreaks=(sqrt(ylabels)-yscmn)/yscsd
      }else
        if(transLU[thesevar[2]]=="log"){
          point_df$y_utr=exp(point_df$y_usc)
          ylabels=pretty(point_df$y_utr)
          ybreaks=(log(ylabels)-yscmn)/yscsd
        }
    }else{
      point_df$y_utr=point_df$y_usc
      ylabels=pretty(point_df$y_utr)
      ybreaks=(ylabels-yscmn)/yscsd
    }    
    
    #Surface
    mod_int=mod$coefficients[1]
    thisdata$x_sc=thisdata[,thesevar[1]]
    thisdata$y_sc=thisdata[,thesevar[2]]
    thisdata$predZ=predict(mod,newdata=thisdata,se.fit=F)
#    thisdata$predZ[which(thisdata$predZ<0)]=0
#    thisdata$predZ[which(thisdata$predZ>100)]=100
    
    thisdata$predZ_minint=thisdata$predZ-mod_int
    
    if(is.null(color_scale)){color_scale=quantile(c(thisdata$predZ,point_df$resZ),c(.1,.9),na.rm=T)}
    
    plotlist[[i]]=ggplot()+
      geom_tile(data=thisdata,aes(x=x_sc,y=y_sc,fill=predZ))+
      geom_point(data=point_df,aes(x=x_sc,y=y_sc,size=mod$weights),alpha=.1)+#,color=z_res))+#,color="white",shape=21)+
      ggtitle(paste0("Interaction Surface: \n",theseplotvar[1]," : ",theseplotvar[2]))+
      scale_size(guide="none", range=c(2,6))+
      #scale_color_divergent(name="Partial Residual",limits=color_scale,oob=squish,
      #                       low="lightgreen",mid="grey75",high="orange",midpoint=0)+
      scale_fill_continuous(name="Predicted Pct. Afffected",limits=color_scale,oob=squish,
                           low="skyblue",high="red")+
      xlab(theseplotvar[1])+
      ylab(theseplotvar[2])+
      scale_x_continuous(breaks = xbreaks,labels=xlabels)+
      scale_y_continuous(breaks = ybreaks,labels=ylabels)+
      theme_pubr()+
      theme(legend.position = "bottom",axis.text=element_text(size=8))
  }
  outp=wrap_plots(plotlist,guides="collect")& theme(legend.position = 'bottom')
  #outlist=list(outp,plotlist)
  return(outp)
}


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


## iv. MODEL PERTURBATIONS #### 
#--> can use this to calculate site susceptibility

# take model preds, compare preds of model with data as is, to preds of model with one of these variables tuned down
# i.e. if effluent was zero, how would that change the prediction of bleaching?
# x: sort blusters from least to most predicted bleaching, y: predicted bleaching (cumulative density function)
# if you perturb model by making a change to one parameter to all sites @ once, what changes in distribution of values? how big of an effect?
# sensitivity analysis: for every variable, what happens when increase/decrease by 10%?
# or choose realistic things to inform management

pertDat_fun <- function(var, var_s, trans = "none", dat, mod, pert = 1){
  v = var
  v_s = var_s
  d = dat
  t = trans
  m = mod
  p = pert
  
  d$Bluster = 1:nrow(d)
  d$predict_ble <- predict(m, d)
  
  if(t == "none"){
    mn_v = mean(d[[v]])
    sd_v = sd(d[[v]])
  }
  
  if(t == "log"){
    mn_v = mean(log(d[[v]]))
    sd_v = sd(log(d[[v]]))
  }
  
  if(t == "sqrt"){
    mn_v = mean(sqrt(d[[v]]))
    sd_v = sd(sqrt(d[[v]]))
  }
  
  #d$original_var <- d[[v_s]] # save a copy of original value (might not need this though)
  #mn_perc_p1 <- mean(d[[v]])*p # calculate x% of the mean (specified by pert val, defaults to 10%)
  
  sd_pert <- sd(d[[v]])*p
  
  ## add x% of mean
  te_pert <- d[[v]] + sd_pert
  # sometimes adding this will make the value exceed the observed range, which we do not want
  te_pert[ te_pert > max(d[[v]], na.rm=T)] = max(d[[v]], na.rm=T) # anything greater than observed range gets set to obs range
  # scale the new values accordingly
  if(t == "none"){
    d[[v_s]] <- (te_pert- mn_v)/sd_v
  }
  if(t == "log"){
    d[[v_s]] <- (log(te_pert)- mn_v)/sd_v
  }
  if(t == "sqrt"){
    d[[v_s]] <- (sqrt(te_pert)- mn_v)/sd_v
  }
  # make predictions based on new values
  d$pred.pert.plus <- predict(m, d)
  
  ## subtract x% of mean
  te_pert <- d[[v]] - sd_pert
  # or sometimes subtracting this will make value fall below observed range, also bad
  te_pert[ te_pert < min(d[[v]], na.rm=T)] = min(d[[v]], na.rm=T) # anything less than observed range gets set to obs range
  if(t == "none"){
    d[[v_s]] <- (te_pert- mn_v)/sd_v
  }
  if(t == "log"){
    d[[v_s]] <- (log(te_pert)- mn_v)/sd_v
  }
  if(t == "sqrt"){
    d[[v_s]] <- (sqrt(te_pert)- mn_v)/sd_v
  }
  d$pred.pert.minus <- predict(m, d)
  
  d2 <- d[,c("Bluster", "predict_ble", "pred.pert.plus", "pred.pert.minus")]
  # d2 <- transform(d2,
  #                 pred.pert.plus=ifelse(pred.pert.plus>10,10,pred.pert.plus))
  # 
  # d2 <- transform(d2,
  #                 pred.pert.minus=ifelse(pred.pert.minus<0,0,pred.pert.minus))
  
  
  d2$ble_inc <- d2$pred.pert.plus - d2$predict_ble
  d2$ble_dec <- d2$predict_ble - d2$pred.pert.minus
  
  colnames(d2)[c(3:6)] <- paste(colnames(d2)[c(3:6)],v)
  return(d2)
}

pertPlots_fun <- function(dat){
  d = dat
  d <- d[ order(d$predict_ble),]
  bl_ord <- d$Bluster # save the order of blusters that are low to high predicted bleaching
  d$Bluster <- factor(d$Bluster, levels = bl_ord) # make the levels of bluster equal to this order
  
  cols <- c("inc" = "orange", "dec" = "cadetblue3")
  
  ggplot() + 
    geom_point(data = d, aes(x = factor(Bluster, level = bl_ord), y = d[,3], color = "inc"), size = 3) +
    geom_point(data = d, aes(x = factor(Bluster, level = bl_ord), y = d[,4], color = "dec"), size = 3) +
    geom_line(data = d, aes(x = factor(Bluster, level = bl_ord), y = predict_ble), group = 1) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.spacing = unit(0, "lines"),
          axis.line = element_line(color = "black"),
          panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = "white"),
          axis.text.y = element_text(colour="black"),
          legend.text = element_text(size = 14),
          text = element_text(size = 14),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          title = element_text(size = 14)) +
    ylab("") +
    xlab("") +
    scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+#breaks = c(0,2.5, 5, 7.5, 10), labels = c(0,5,25,56.25,100)) +
    scale_color_manual(name = "", 
                       breaks = c("inc", "dec"), 
                       values = cols,
                       labels = c("Increase", "Decrease"))
}

