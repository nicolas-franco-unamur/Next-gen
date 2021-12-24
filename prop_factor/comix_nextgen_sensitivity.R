###############################################################################
# coMix nextgen approach - additional script for generating sensitivity study
# Estimating susceptibility/Infectiousness - all ages
# Nicolas Franco - UHasselt
# Version:  December 223 2021
###############################################################################

# cleaning if needed
rm(list = ls()); gc()   

#path to nextgen working directory
library(rstudioapi)    
path_to_wd <- paste0(dirname(rstudioapi::getActiveDocumentContext()$path),"/")
#path to socrates (absolute or relative to wd)
path_to_socrates <- "socrates_rshiny-master/"
#path to csv file with confirmed PCR test (absolute or relative to wd)
#This file must be requested from Sciensano
PCR_file <- "RDS"

#age groups to consider (last group going up to infinity)
Age_breaks <- c(0,6,12,18,30,40,50,60,70,80)
#Choose a asssumption on Suscepticility or Infectiousness (only one)
#Fixed_Susceptibility <- c(replicate(length(Age_breaks),1))
#Fixed_Susceptibility <- c(0.4,0.39,0.38,0.79,0.86,0.8,0.82,0.88,0.74,0.74)  #Davies Nature
#Fixed_Infectiousness <- c(replicate(length(Age_breaks),1))
Fixed_Infectiousness <- c(0.54,0.55,0.56,0.59,0.7,0.76,0.9,0.99,0.99,0.99)  #Abrams
#starting and ending waves
Wavemin <- 12
Wavemax <- 23
#size of wave group
Wavegroupsize <- 12
#wave group size = min size
Wavegroupsizemin <- 0
#Number of samples for bootstrap (>1) - only the first is taken here !
N_bs <- 2
#Number of samples for sensitivity
N_sens <- 200
#stop criteria: number of iterations with no change (stop criteria)
N_ite <- 100
#Delay before starting survey and PCR considered period (default = 7)
Delay <- 7
#length of PCR considered period (default = 14)
Period_length <- 14

#libraries
library(dplyr)
library(ggplot2)
library(tidyverse)
library(xtable)

#load functions
source(paste0(path_to_wd,'comix_nextgen_compute.R'))
source(paste0(path_to_wd,'comix_nextgen_appendpcr.R'))

#apply nextgen on socrates
set.seed(20200101)
comix_nextgen <- nextgen_compute(Wavemin = Wavemin,
                Wavemax = Wavemax,
                path_to_wd = path_to_wd,
                path_to_socrates = path_to_socrates,
                Age_breaks = Age_breaks,
                N_bs = N_bs);

#append pcr positive tests
comix_nextgen <- append_pcr(comix_nextgen,
                            Wavemin = Wavemin,
                            Wavemax = Wavemax,
                            path_to_wd = path_to_wd,
                            PCR_file = PCR_file,
                            Age_breaks = Age_breaks,
                            Delay = Delay,
                            Period_length = Period_length
)

#load socrates
setwd(path_to_wd)
setwd(path_to_socrates)
source('R/socrates_main.R')
source('R/load_config_base.R')

#initialisation
wavematrix_bs <- vector(mode = "list", length = Wavemax)
Wavegroups<- ((Wavemin:Wavemax-Wavemin)%/%Wavegroupsize) +1
Wavegroupsnum <- (Wavemax-Wavemin)%/%Wavegroupsize +1
if(Wavegroupsizemin ==1 && ((Wavemax-Wavemin+1)%%Wavegroupsize) > 0){
  Wavegroups[Wavegroups==Wavegroupsnum] <- Wavegroupsnum-1
  Wavegroupsnum <- Wavegroupsnum-1
} 

#load reference matrices
set.seed(20200101)
for(i in Wavemin:Wavemax){
  survey_obj <- get_survey_object(country="Belgium 2020 CoMix (Coletti 2020)",
                                  daytype="All contacts",
                                  touch="All contacts",
                                  duration="All contacts",
                                  gender="All",
                                  cnt_location=opt_location,
                                  bool_reciprocal=TRUE,
                                  bool_suppl_professional_cnt=TRUE,
                                  bool_hhmatrix_selection=FALSE,
                                  wave=i,
                                  quiet = TRUE)
  
  matrix_out_bs <- contact_matrix(survey_obj, 
                                  age.limits = Age_breaks,
                                  symmetric  = TRUE,
                                  quiet      = TRUE,
                                  weigh.dayofweek = TRUE,
                                  weigh.age = TRUE,
                                  weight.threshold = weight_threshold,
                                  estimated.contact.age = 'sample',
                                  missing.contact.age = 'sample',
                                  n=N_bs
  )
  wavematrix_bs[[i]] <- lapply(matrix_out_bs$matrices, function(x) {x$matrix[is.na(x$matrix)] <- 0; return(x)})
  wavematrix_bs[[i]] <- lapply(wavematrix_bs[[i]], function(x) {x$matrix <- t(x$matrix)   ; return(x)})
}
 

#initialise vectors
susinflist <- vector("list", length = Wavegroupsnum) 
susinfvec <- vector("list", length = Wavegroupsnum) 
newsusinfvec <- vector("list", length = Wavegroupsnum) 
for(group in 1:Wavegroupsnum){
  susinflist[[group]] <- vector("list", length = N_sens)  
  susinfvec[[group]] <- c(replicate(length(Age_breaks),1))
}

# sensitivity loop
# stay at bs=1
bs <- 1
Fixed_Infectiousness_basis <- Fixed_Infectiousness
for(sens in 1:N_sens){
  
Fixed_Infectiousness <-  Fixed_Infectiousness_basis + runif(length(Age_breaks),min=-0.1,max=+0.1)
  
  lklh <- 1000000
  ite_without_change <- 0
  ite <- 0
  wavematrixnew_bs  <- vector(mode = "list", length = Wavemax)
  comix_bs <- vector(mode = "list", length = Wavemax)
  
  while(ite_without_change < N_ite){
    
  newsusinfvec<-lapply(susinfvec,function(x) {rnorm(length(Age_breaks),x,0.005)} )
  newsusinfvec<-lapply(newsusinfvec,function(x) {x[x<0] <- 0 ; return(x)} )

  counter<-0
  newlklh<-0
  
  for(i in Wavemin:Wavemax){
    
    if(exists("Fixed_Susceptibility")){
      Susceptibility <- Fixed_Susceptibility
    } else {
      Susceptibility <- newsusinfvec[[Wavegroups[[i-Wavemin+1]]]]
    }
    if(exists("Fixed_Infectiousness")){
      Infectiousness <- Fixed_Infectiousness
    } else {
      Infectiousness <- newsusinfvec[[Wavegroups[[i-Wavemin+1]]]]
    }
    
    #compute eigenvectors
    wavematrixnew_bs[[i]] <-  diag(Susceptibility) %*% wavematrix_bs[[i]][bs][[1]]$matrix %*% diag(Infectiousness) 
    comix_bs[[i]] <- data.frame(Age_breaks,abs(eigen(wavematrixnew_bs[[i]])$vectors[,1]))
    colnames(comix_bs[[i]] ) <- c("ageclass","nextgenmean")
    comix_bs[[i]]  <- comix_bs[[i]]   %>% mutate(nextgenmean=nextgenmean *100 / sum(comix_bs[[i]]$nextgenmean))
  
    for(ages in 1:(length(Age_breaks))){
      counter <- counter + 1
      #Hellinger distance
      newlklh <- newlklh + (sqrt(comix_bs[[i]]$nextgenmean[ages]) -sqrt(comix_nextgen[[i]]$pcr[ages]))^2
    }
    
  }
  
  ite_without_change <- ite_without_change+1
  ite <- ite +1
  
  if(newlklh<lklh){
    lklh<-newlklh
    ite_without_change <- 0
    susinfvec<-newsusinfvec
    #print(c(ite,bs,lklh,susinfvec),digits = 2)
  }
  
  }
  print(c(ite,bs,sens,lklh,susinfvec),digits = 2)
  
  for(group in 1:Wavegroupsnum){
    susinflist[[group]][[sens]] <- susinfvec[[group]]
  }

}

#Collecting final results
for(sens in 1:N_sens){
  bs <- 1

  for(i in Wavemin:Wavemax){
    if(exists("Fixed_Susceptibility")){
      Susceptibility <- Fixed_Susceptibility
    } else {
      Susceptibility <- susinflist[[Wavegroups[[i-Wavemin+1]]]][[sens]]
    }
    if(exists("Fixed_Infectiousness")){
      Infectiousness <- Fixed_Infectiousness
    } else {
      Infectiousness <- susinflist[[Wavegroups[[i-Wavemin+1]]]][[sens]]
    }
    wavematrix_bs[[i]][bs][[1]]$matrix <-  diag(Susceptibility) %*% wavematrix_bs[[i]][bs][[1]]$matrix %*% diag(Infectiousness) 
  }
}

comix_nextgenfinal <- vector(mode = "list", length = Wavemax)
for(i in Wavemin:Wavemax){
  comix_nextgenfinal[[i]]$bs <- lapply(wavematrix_bs[[i]],function(x) {data.frame(Age_breaks,abs(eigen(x$matrix)$vectors[,1]))})
  comix_nextgenfinal[[i]]$bs <- lapply(comix_nextgenfinal[[i]]$bs, setNames, c("ageclass","nextgen"))
  comix_nextgenfinal[[i]]$bs <- lapply(comix_nextgenfinal[[i]]$bs,function(x) {mutate(x,nextgen=nextgen *100 / sum(x$nextgen))  })
  comix_nextgenfinal[[i]] <-  bind_rows(comix_nextgenfinal[[i]]$bs)  %>% group_by(ageclass) %>% summarise(nextgenmean = mean(nextgen), nextgenmedian = median(nextgen),  nextgensd = sd(nextgen), nextgenlower = unname(quantile(nextgen,0)), nextgenupper = unname(quantile(nextgen,1)), .groups = 'drop')
  comix_nextgenfinal[[i]]$wave <- i
  comix_nextgenfinal[[i]]$wavename <- paste("Wave",i)
  comix_nextgenfinal[[i]]$ageclassname <- factor(matrix_out_bs$demography$age.group, levels = matrix_out_bs$demography$age.group)
  comix_nextgenfinal[[i]]$ageclasslower <- factor(matrix_out_bs$participants$lower.age.limit, levels = matrix_out_bs$participants$lower.age.limit)
  }

#append pcr positive tests
comix_nextgenfinal <- append_pcr(comix_nextgenfinal,
                            Wavemin = Wavemin,
                            Wavemax = Wavemax,
                            path_to_wd = path_to_wd,
                            PCR_file = PCR_file,
                            Age_breaks = Age_breaks,
                            Delay = Delay,
                            Period_length = Period_length
)

for(i in Wavemin:Wavemax){
  comix_nextgen[[i]]$method <- "CoMix - homogeneous"
  comix_nextgenfinal[[i]]$method <- "CoMix - heterogeneous"
  comix_nextgenfinal[[i]] <- bind_rows(comix_nextgen[[i]],comix_nextgenfinal[[i]])
  comix_nextgenfinal[[i]]$method <-  factor(comix_nextgenfinal[[i]]$method,levels=c("CoMix - homogeneous","CoMix - heterogeneous"))
}

#preparing plots
ymin <- 0 ; ymax <- 0
for(i in Wavemin:Wavemax){
  ymin <- min(ymin,comix_nextgenfinal[[i]]$pcr,comix_nextgenfinal[[i]]$nextgenlower)
  ymax <- max(ymax,comix_nextgenfinal[[i]]$pcr,comix_nextgenfinal[[i]]$nextgenupper)
}
ymax <- min(ymax,30)

comix_facets <- bind_rows(comix_nextgenfinal)
comix_facets$wavename <- factor(comix_facets$wavename, levels = paste("Wave",Wavemin:Wavemax))

# backup
susinflistcopy <- susinflist  
susinflist <- susinflistcopy

# plot relative incidence
plotrelincid <- ggplot(comix_facets) +
  geom_point(aes(x = factor(ageclasslower), y = pcr,colour = "PCR positive tests"),size=1) +
  geom_point(aes(x = factor(ageclasslower), y = nextgenmean,colour = method),position =position_dodge(width = 0.8),size=1) +
  geom_linerange(aes(x = factor(ageclasslower), y = nextgenmean,colour = method,ymin=nextgenlower, ymax=nextgenupper),position =position_dodge2(width = 0.8),size=0.6) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  theme_minimal()+
  scale_shape_identity()   +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ylab("Relative incidence")  +
  xlab("Age classes") +
  facet_wrap(~wavename, ncol=4) +
  guides(color = guide_legend(reverse = TRUE)) +
 # guides(color = guide_legend(direction = "horizontal",ncol = 1,label.position = "bottom",label.theme = element_text(angle = 90))) +
  theme(legend.position = "bottom",legend.title=element_blank()) +
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position = "bottom",axis.ticks.x=element_blank()) 

#preparing plot on susceptibility/infectiousness
susinfdataframe <- vector("list", length = Wavegroupsnum) 
for(group in 1:Wavegroupsnum){
  #normalisation
  #susinflist[[group]] <-  lapply(susinflist[[group]],function(x) {x/x[4]})  #normalisation on young adults
  #susinflist[[group]] <-  lapply(susinflist[[group]],function(x) {x/mean(x[1:3])})  #normalisation on childrens
}
  maxval <- mean(unlist(lapply(susinflist[[1]],"[",4)))
for(group in 1:Wavegroupsnum){
  susinflist[[group]] <-  lapply(susinflist[[group]],function(x) {x/maxval})
  susinfdataframe[[group]] <- t(sapply(susinflist[[group]],c))
  colnames(susinfdataframe[[group]]) <- c("[0,6)","[6,12)","[12,18)","[18,30)","[30,40)","[40,50)","[50,60)","[60,70)","[70,80)","80+")
  susinfdataframe[[group]] <- stack(as.data.frame(susinfdataframe[[group]]))
  colnames(susinfdataframe[[group]]) <- c("values","ageclass")
  susinfdataframe[[group]] <- susinfdataframe[[group]]  %>% group_by(ageclass) %>% summarise(mean = mean(values), median = median(values),  sd = sd(values), lower = unname(quantile(values,0)), upper = unname(quantile(values,1)), .groups = 'drop')
  susinfdataframe[[group]]$wave <- paste("Waves",paste(which(Wavegroups %in% group)+Wavemin-1, collapse = " "))
}
commondataframe <- bind_rows(susinfdataframe)
commondataframe$wave <- factor(commondataframe$wave)

# #plot sus/inf
# plotsusinf <- ggplot(commondataframe) +
#   geom_point(aes(x = wave, y=mean,group = 1,colour = factor(wave)), size=2)   +
#   geom_line(aes(x = wave, y=mean,group = 1), colour="darkorange",linetype = "solid", size=1)   +
#   geom_linerange(aes(x = wave, y = mean,ymin=lower, ymax=upper,colour = factor(wave)), size=0.7) +
#   geom_ribbon(aes(x = wave,ymin=lower, ymax=upper,group = 1) ,fill="darkorange", alpha=0.2) +
#   facet_wrap(~ageclass, strip.position="bottom", ncol=length(Age_breaks)) +
#   geom_hline(yintercept=1,color = "grey") +
#   theme_minimal()+
#   scale_x_discrete(breaks = NULL) +
#   guides(colour=guide_legend("CoMix waves groups:")) +
#   theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position = "bottom",axis.ticks.x=element_blank()) +
#   coord_cartesian(ylim = c(0, 10)) +
#   #ylab("Relative q-susceptibility") +
#   ylab("Relative q-infectiousness") +
#   xlab("Age classes")

#plot sus/inf
plotsusinf <- ggplot(commondataframe) +
  geom_point(aes(x = wave, y=mean,group = 1,colour = factor(wave)), size=2*1.5)   +
 # geom_line(aes(x = wave, y=mean,group = 1), colour="darkorange",linetype = "solid", size=1)   +
  geom_errorbar(aes(x = wave, y = mean,ymin=lower, ymax=upper,colour = factor(wave)), size=0.7*1.5) +
  #geom_ribbon(aes(x = wave,ymin=lower, ymax=upper,group = 1) ,fill="darkorange", alpha=0.2) +
  facet_wrap(~ageclass, strip.position="bottom", ncol=length(Age_breaks)) +
  geom_hline(yintercept=1,color = "grey") +
  theme_minimal()+
  scale_x_discrete(breaks = NULL) +
  guides(colour=guide_legend("CoMix waves groups:")) +
  theme(plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust = 0.5),legend.position = "none",axis.ticks.x=element_blank()) +
  coord_cartesian(ylim = c(0, 1.5)) +
  ylab("Relative q-susceptibility") +
  #ylab("Relative q-infectiousness") +
  xlab("Age classes")
plotsusinf
 
#save figures and workspace
plotrelincid
ggsave("fig_relincid_sensivity.eps", device=cairo_ps, width = 9, height = 6*9/11)
plotsusinf
ggsave("fig_susinf_sensivity.eps", device=cairo_ps, fallback_resolution = 1200, width = 9*0.85, height = 6*0.85)
#ggsave("fig_susinf.eps", device=cairo_ps, width = 9*0.85*0.65, height = 6*0.85*0.65)


pdf("fig_relincid_sensivity.pdf", width = 9, height = 6*9/11)
print(plotrelincid);
dev.off() 
pdf("fig_susinf_sensivity.pdf", width = 9*0.85, height = 6*0.85)
print(plotsusinf);
dev.off()
save.image("nextgen_sensivity.RData") 

#latex output
#for(group in 1:Wavegroupsnum){
#print(xtable(susinfdataframe[[group]][,1:6 ], digits = 3,caption = toString(susinfdataframe[[group]]$wave[1])), include.rownames = FALSE)
#}

