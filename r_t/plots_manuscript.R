################################################################################### #
########### TBD                                                           ######### #
################################################################################### #

### Created 03/09/2021 by Pietro Coletti   

rm(list=ls(all=TRUE))
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}
library(ggplot2)
library(dplyr)
library(plyr)
library(stringr)
library(data.table)
library(tidyr)

`%notin%` <- function(a,b) ! a %in% b  

naming_function<-function(int_id_wave,int_age_breaks,int_susc_vec,int_inf_vec){
  int_age_breaks_text <- paste(int_age_breaks,collapse=',')
  int_susc_vec_text <- paste(int_susc_vec,collapse=',')
  int_inf_vec_text <- paste(inf_vec,collapse=',')
  
  literature_susc_vec=c("0.4","0.39","0.38","0.79","0.86","0.8","0.82","0.88","0.74","0.74")
  literature_inf_vec=c("0.54","0.55","0.56","0.59","0.7","0.76","0.9","0.99","0.99","0.99")
  
  fname<-"./results_draft_NGA/R0_ratios_wave_"
  fname<-paste0(fname,id_wave,"_full_AC_",str_replace_all(int_age_breaks_text,",","_"))              
                
  inf_tag<-paste0("_varied_inf_",str_replace_all(int_inf_vec_text,",","_"))
  if(all(int_inf_vec==literature_inf_vec)){inf_tag="_lit_inf"}
  if(all(int_inf_vec==1)){inf_tag="_homogeneous_inf"}
  fname<-paste0(fname,inf_tag)
  
  susc_tag<-paste0("_varied_susc_",str_replace_all(int_susc_vec_text,",","_"))
  if(all(int_susc_vec==literature_susc_vec)){susc_tag="_lit_susc"}
  if(all(int_susc_vec==1)){susc_tag="_homogeneous_susc"}
  fname<-paste0(fname,susc_tag,".RDS")
}


### Change needed
list_id_waves<-12:23
dates_waves<-c("23/12","05/01","19/01","03/02","17/02","01/03","15/03","30/03","13/04","27/04",
               "13/05","25/05","09/06","23/06","06/07","20/07","03/08"
               )
waves_levels<-c("Round 12 (23/12)","Round 13 (05/01)","Round 14 (19/01)",
                "Round 15 (03/02)","Round 16 (17/02)","Round 17 (01/03)",
                "Round 18 (15/03)","Round 19 (30/03)","Round 20 (13/04)",
                "Round 21 (27/04)","Round 22 (13/05)","Round 23 (25/05)",
                "Round 24 (09/06)","Round 25 (23/06)","Round 26 (06/07)",
                "Round 27 (20/07)","Round 28 (03/08)"
                )

age_breaks<- c(0,6,12,18,30,40,50,60,70,80)
age_breaks_text <- paste(age_breaks,collapse=',')
q<-c(0.025,0.25,0.75,0.975)

types_of_fits<-c("Homogenoeus","fitted_inf","fitted_susc")

if(exists("df_tot")){rm(df_tot)}
for(type_data in types_of_fits){
  if(exists("df_ratios")){rm(df_ratios)}
  if(type_data=="Homogenoeus"){
    susc_vec<-rep(1,length(age_breaks))
    inf_vec<-rep(1,length(age_breaks))
  }
  if(type_data=="fitted_susc"){
    susc_vec<-c("0.182","0.550","0.603","1.000","1.172","1.009","0.880","0.869","0.846","0.805")
    inf_vec<-c("0.54","0.55","0.56","0.59","0.7","0.76","0.9","0.99","0.99","0.99")
  }
  if(type_data=="fitted_inf"){
    susc_vec<-c("0.4","0.39","0.38","0.79","0.86","0.8","0.82","0.88","0.74","0.74")
    inf_vec<-c("0.346","0.892","1.310","1,","0.645","3.783","1.32","0.266","1.277","0.099")
  }
  if(type_data=="time_varying"){
    inf_vec<-c("0.54","0.55","0.56","0.59","0.7","0.76","0.9","0.99","0.99","0.99")
  }
  
  
for(id_wave in list_id_waves){
  if(type_data=="time_varying"){
    i_data_frame<-as.integer((id_wave-12)/2)+1
    susc_vec<-susinfdataframe[[i_data_frame]]$mean
  }
  
    irow<-which(list_id_waves==id_wave)
    fname<-naming_function(id_wave,age_breaks_text,susc_vec,inf_vec)
    print(fname)
    QOI<-readRDS(file =fname)
    formatted_date<-as.Date(paste0(dates_waves[irow],"/2020"),format = "%d/%m/%Y")
    if(id_wave>12){formatted_date<-as.Date(paste0(dates_waves[irow],"/2021"),format = "%d/%m/%Y")}
    print(formatted_date)
    dum_df_ratios=data.frame(
      "ratio" =QOI,"wave"=paste0("Round ",id_wave,"\n (",dates_waves[irow],")" ),
      "date"=formatted_date,
      "id_wave"=id_wave,
      "type"=type_data
    )
    #waves_levels_full[id_wave]<-paste0("Round ",id_wave,"\n (",dates_waves[irow],")" )
    if(!exists("df_ratios")){df_ratios<-dum_df_ratios}else{
      df_ratios<-rbind(df_ratios,dum_df_ratios)}
  }
  
  df2 <- df_ratios %>%
    group_by(wave) %>%
    dplyr::summarise(
      quant2_5  = quantile(ratio, probs = q[1]),
      quant25   = quantile(ratio, probs = q[2]),
      median    = median(ratio),
      quant75   = quantile(ratio, probs = q[3]),
      quant97_5 = quantile(ratio, probs = q[4]),
      date=date, 
      id_wave=id_wave,
      type=type
    )
  df2 <- df2 %>% distinct()
  if(!exists("df_tot")){df_tot<-df2}else{
    df_tot<-rbind(df_tot,df2)}
}


  
######################################################################################################################################################
############################## COMPARISON OF R  (DSI)
######################################################################################################################################################
library(rjson)

###########
result <- fromJSON(file = "./r_data_cases.json")
national_results<-list()
i_count<-1
for(i in 1:length(result)){
  a<-result[[i]]
  if(a$region=="Belgique"){
    national_results[[i_count]]<-a
    i_count<-i_count+1
  }
}
if(exists("R_DSI")){rm(R_DSI)}
for(i in 1:length(national_results)){
  oneline<-data.frame(national_results[[i]])
  if(!exists("R_DSI")){R_DSI<-oneline}else{
    R_DSI<-rbind(R_DSI,oneline)
  }
}
R_DSI$date<-as.Date(R_DSI$date, format="%Y-%m-%d")

###########

colnames(R_DSI)<-c("t_start","t_end","Mean_R","Std_R","Quantile.0.025_R","Quantile.0.05_R","Quantile.0.25_R","Median_R",
                   "Quantile.0.75_R","Quantile.0.95_R","Quantile.0.975_R","date","region")


R_DSI$rollmean     <-frollmean(R_DSI[,8], 7)
R_DSI$rollmean_97_5<-frollmean(R_DSI[,11], 7)
R_DSI$rollmean_2_5 <-frollmean(R_DSI[,5], 7)
R_DSI<-subset(R_DSI,R_DSI$date>as.Date("2020-12-20"))



shift=7
val_ref<- R_DSI[R_DSI$date==as.Date("2020-12-23")+shift,]$rollmean
wave_calibration<-12
factor_homogeneous<-val_ref/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="Homogenoeus")$median
factor_fitted_inf<-val_ref/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="fitted_inf")$median
factor_time_var<-val_ref/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="time_varying")$median
factor_fitted_suscr<-val_ref/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="fitted_susc")$median

df_tot$factor<-factor_homogeneous
df_tot$factor[df_tot$type=="fitted_inf"]<-factor_fitted_inf
df_tot$factor[df_tot$type=="time_varying"]<-factor_time_var
df_tot$factor[df_tot$type=="fitted_susc"]<-factor_fitted_suscr
dates_plot<-df_tot$date+shift
R_DSI_sub<-subset(R_DSI,R_DSI$date%in%dates_plot)

p1<-ggplot()+
  geom_errorbar(data=df_tot,aes(x=date+shift,
                                y = median*factor,ymin=quant2_5*factor, ymax=quant97_5*factor,
                                colour = type),size=0.2,position=position_dodge2(2))+
  geom_errorbar(data=R_DSI_sub,aes(x=date,y=rollmean,ymin=rollmean_2_5,ymax=rollmean_97_5,),color="black",size=0.2)+
  ylab("Rt (from DSI)")+xlab("Date")+scale_x_date(limits = as.Date(c('2020-12-20','2021-08-31')))+ theme(legend.position = "none")


dev.off()

sub_plot_labels_a<-data.frame(
   x=as.Date(2020-12-20)
  ,y=c(1)
  ,text=c("a)")
  
)


df_tot_sub<-subset(df_tot,df_tot$type%in%c("fitted_susc","Homogenoeus"))


R_DSI_attaching<-data.frame("wave"=NA
                                  ,"quant2_5"=R_DSI$Quantile.0.025_R
                                  ,"quant25"=NA
                                  ,"median"=R_DSI$Median_R
                                  ,"quant75"=NA
                                  ,"quant97_5"=R_DSI$Quantile.0.975_R
                                  ,"date"=R_DSI$date
                                  ,"id_wave"=NA
                                  ,"type"="Infections data"
                                  ,"factor"=NA
                                  ,"factor_hosp"=NA)
df_test<-rbind(df_tot_sub,R_DSI_attaching)
df_test$type_comix<-NA
df_test$type_comix[df_test$type=="fitted_susc"]<-"fitted_susc"
df_test$type_comix[df_test$type=="Homogenoeus"]<-"Homogenoeus"

df_test$is_ext_data<-NA
df_test$is_ext_data[df_test$type=="Infections data"]<-1
df_test$date_rescaled<-df_test$date
df_test$date_rescaled[df_test$type %in% c("fitted_susc","Homogenoeus")]   <-df_test$date_rescaled[df_test$type %in% c("fitted_susc","Homogenoeus")]+shift



p3_old<-ggplot()+geom_linerange(data=df_tot_sub,aes(x=date+shift,
                                                y = median*factor,ymin=quant2_5*factor, ymax=quant97_5*factor,
                                                colour = type,alpha="CoMix"),size=0.7,position=position_dodge2(1))+
  geom_point(data=df_tot_sub,aes(x=date+shift,
                                 y = median*factor,colour = type),size=2,position=position_dodge2(1))+
  geom_line(data=R_DSI,mapping=aes(x=date,y=Mean_R,alpha="Infections data"),color="black",size=1)+
  geom_ribbon(data=R_DSI,aes(x=date,ymin = Quantile.0.025_R, ymax = Quantile.0.975_R, fill = "red"), alpha = 0.3,show.legend=FALSE)+
  theme_minimal()+
  ylab(bquote(R[t]))+xlab("Date")+scale_x_date(limits = as.Date(c('2020-12-20','2021-06-30')))+
  scale_color_manual("",breaks = c("Homogenoeus", "fitted_susc"), values=c("#00ba38","#f8766d")
                     ,labels=c("CoMix-homogenous","CoMix-heterogenous")
  )+  scale_alpha_manual(name = NULL,
                     values = c(2, NA),
                     breaks = c("Infections data", "") )+ 
  labs(tag = "a)") +
  theme(plot.tag.position = c(0.01, 0.99))
  
p3<-ggplot(df_test)+
  geom_errorbar(aes(x=date_rescaled,
                    y = median*factor,ymin=quant2_5*factor, ymax=quant97_5*factor,colour = type_comix),size=0.9,position=position_dodge2(1))+
  geom_point(aes(x=date_rescaled,
                 y = median*factor,colour = type_comix),size=2.2,position=position_dodge2(1))+
  ylab(bquote(R[t]))+xlab("Date")+scale_x_date(limits = as.Date(c('2020-12-20','2021-06-30')))+
  geom_line(data=subset(df_test,df_test$type=="Infections data"),aes(x=date_rescaled,
                                                                          y = median*is_ext_data) ,color="black",size=1.2,show.legend = TRUE)+
  geom_ribbon(data=subset(df_test,df_test$type=="Infections data"),aes(x=date_rescaled,ymin = quant2_5, ymax = quant97_5, fill = type), alpha = 0.3)+
  theme_minimal()+
  scale_color_manual("",breaks = c("Homogenoeus", "fitted_susc"), values=c("#00ba38","#f8766d")
                     ,labels=c("CoMix-homogenous","CoMix-heterogenous"),na.translate = F
  )+
  scale_fill_manual("",breaks=c("Infections data"),values=c("red"),labels=c("Infections data"))+
  guides(fill=guide_legend("",override.aes=list(linewidth=10,linetype=1),labels=waiver() ))+
  guides(colour=guide_legend(override.aes=list(shape=c(16,16), linetype=c(NA,NA))))+
  labs(tag = "a)") +
  theme(plot.tag.position = c(0.01, 0.99))+theme(legend.position = "bottom" )+
  theme(text = element_text(size = 20))     



######################################################################################################################################################
############################## COMPARISON OF R  (Sciensano)
######################################################################################################################################################
R_sciensano<-read.csv("R_sciensano.csv")
R_sciensano$date<-seq(as.Date("2020-03-23"), as.Date("2021-09-10"), by="days")

shift_hosp=14
val_ref_hosp<- R_sciensano[R_sciensano$date==as.Date("2020-12-23")+shift_hosp,]$Rt
wave_calibration_hosp<-12
factor_homogeneous<-val_ref_hosp/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="Homogenoeus")$median
factor_fitted_inf<-val_ref_hosp/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="fitted_inf")$median
factor_fitted_susc<-val_ref_hosp/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="fitted_susc")$median
factor_time_var<-val_ref_hosp/subset(df_tot,df_tot$id_wave==wave_calibration&df_tot$type=="time_varying")$median


df_tot$factor_hosp<-factor_homogeneous
df_tot$factor_hosp[df_tot$type=="fitted_inf"]<-factor_fitted_inf
df_tot$factor_hosp[df_tot$type=="fitted_susc"]<-factor_fitted_susc
df_tot$factor_hosp[df_tot$type=="time_varying"]<-factor_time_var
dates_plot<-df_tot$date+shift_hosp
R_sciensano_sub<-subset(R_sciensano,R_sciensano$date%in%dates_plot)

p2<-ggplot()+geom_errorbar(data=df_tot,aes(x=date+shift_hosp,
                                      y = median*factor_hosp,ymin=quant2_5*factor_hosp, ymax=quant97_5*factor_hosp,
                                      colour = type),size=0.2,position=position_dodge2(2))+
  geom_point(data=df_tot,aes(x=date+shift_hosp,
                                y = median*factor_hosp,colour = type),size=0.2,position=position_dodge2(11))+
  geom_errorbar(data=R_sciensano_sub,aes(x=date,y=Rt,ymin=Rt_lower,ymax=Rt_upper,),color="black",size=0.2)+
  ylab("Rt (from Sciensano)")+xlab("Date")+scale_x_date(limits = as.Date(c('2020-12-20','2021-08-31')))+geom_line(data=R_DSI,mapping=aes(x=date,y=Mean_R))

library(gridExtra)




df_tot_sub<-subset(df_tot,df_tot$type%in%c("fitted_susc","Homogenoeus"))


p4_old<-ggplot()+geom_linerange(data=df_tot_sub,aes(x=date+shift_hosp,
                                                y = median*factor_hosp,ymin=quant2_5*factor_hosp, ymax=quant97_5*factor_hosp,
                                                colour = type,alpha="CoMix"),size=0.7,position=position_dodge2(1))+
  geom_point(data=df_tot_sub,aes(x=date+shift_hosp,
                                 y = median*factor_hosp,colour = type),size=2,position=position_dodge2(1))+
  geom_line(data=R_sciensano,mapping=aes(x=date,y=Rt,alpha="Hospitalizations data"),color="black",size=1)+
  geom_ribbon(data=R_sciensano,aes(x=date,ymin = Rt_lower, ymax = Rt_upper, fill = "red"), alpha = 0.3,show.legend=FALSE)+
  theme_minimal()+
  ylab(bquote(R[0]))+xlab("Date")+scale_x_date(limits = as.Date(c('2020-12-20','2021-06-30')))+
  scale_color_manual("",breaks = c("Homogenoeus", "fitted_susc"), values=c("#00ba38","#f8766d")
                     ,labels=c("CoMix-homogenous","CoMix-heterogenous")
  )+  scale_alpha_manual(name = NULL,
                         values = c(2, NA),
                         breaks = c("Hospitalizations data", "") ,
   )+ 
  labs(tag = "b)") +
  theme(plot.tag.position = c(0.01, 0.99))+
  scale_linetype_manual(values=c(1,2,NA))+
    scale_shape_manual(values=c(NA,NA,2))
R_sciensano_attaching<-data.frame("wave"=NA
                                  ,"quant2_5"=R_sciensano$Rt_lower
                                  ,"quant25"=NA
                                  ,"median"=R_sciensano$Rt
                                  ,"quant75"=NA
                                  ,"quant97_5"=R_sciensano$Rt_upper
                                  ,"date"=R_sciensano$date
                                  ,"id_wave"=NA
                                  ,"type"="Hospitalization data"
                                  ,"factor"=NA
                                  ,"factor_hosp"=NA)
df_test<-rbind(df_tot_sub,R_sciensano_attaching)
df_test$type_comix<-NA
df_test$type_comix[df_test$type=="fitted_susc"]<-"fitted_susc"
df_test$type_comix[df_test$type=="Homogenoeus"]<-"Homogenoeus"

df_test$is_ext_data<-NA
df_test$is_ext_data[df_test$type=="Hospitalization data"]<-1
df_test$date_rescaled<-df_test$date
df_test$date_rescaled[df_test$type %in% c("fitted_susc","Homogenoeus")]   <-df_test$date_rescaled+shift_hosp



p4<-ggplot(df_test)+
  geom_errorbar(aes(x=date_rescaled,
                    y = median*factor_hosp,ymin=quant2_5*factor_hosp, ymax=quant97_5*factor_hosp,colour = type_comix),size=0.9,position=position_dodge2(1))+
  geom_point(aes(x=date_rescaled,
                 y = median*factor_hosp,colour = type_comix),size=2.2,position=position_dodge2(1))+
  ylab(bquote(R[t]))+xlab("Date")+scale_x_date(limits = as.Date(c('2020-12-20','2021-06-30')))+
  geom_line(data=subset(df_test,df_test$type=="Hospitalization data"),aes(x=date_rescaled,
                                                                          y = median*is_ext_data) ,color="black",size=1.2,show.legend = TRUE)+
  geom_ribbon(data=subset(df_test,df_test$type=="Hospitalization data"),aes(x=date_rescaled,ymin = quant2_5, ymax = quant97_5, fill = type), alpha = 0.3)+
  theme_minimal()+
  scale_color_manual("",breaks = c("Homogenoeus", "fitted_susc"), values=c("#00ba38","#f8766d")
                     ,labels=c("CoMix-homogenous","CoMix-heterogenous"),na.translate = F
  )+
  scale_fill_manual("",breaks=c("Hospitalization data"),values=c("red"),labels=c("Hospitalizations data"))+
  guides(colour=guide_legend(order=1,override.aes=list(shape=c(16,16), linetype=c(NA,NA))))+
  guides(fill=guide_legend(order=2,"",override.aes=list(linewidth=10,linetype=1),labels=waiver() ))+
  labs(tag = "b)") +
  theme(plot.tag.position = c(0.01, 0.99))+theme(legend.position = "bottom" ) +
  theme(text = element_text(size = 20))



a<-grid.arrange(p3,p4,ncol=2,widths=c(0.48, 0.52))
print(a)
ggsave("Fig5.eps",grid.arrange(p3,p4,ncol=2,widths=c(0.48, 0.52)), device=cairo_ps, fallback_resolution = 1200, width = 14*1.2, height = 6*1.2)
dev.off()
