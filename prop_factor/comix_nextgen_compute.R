###############################################################################
# Function for next generation approach on CoMix data
# Load socrates data and commpute next generation 
# Nicolas Franco - UHasselt
# Version: Augustus 4 2021
###############################################################################

  
#libraries
library(dplyr)

nextgen_compute <- function(Wavemin = 9,
                            Wavemax = 9,
                            path_to_wd = ".",
                            path_to_socrates = "../socrates_rshiny-master/",
                            Age_breaks = c(18,30,40,50,60),
                            Susceptibility = c(1,1,1,1,1),
                            Infectiousness = c(1,1,1,1,1),
                            N_bs=1,
                            rds="no"
                            ) {
  
  #permanent seed (to remove if several runs)
  set.seed(20200101)
  
  setwd(path_to_wd)
  setwd(path_to_socrates)
  source('R/socrates_main.R')
  source('R/load_config_base.R')
  
  if(rds!="no"){
    matrix_rds_bs <- readRDS(file = paste0(path_to_wd,rds))
  }
  
  if(length(Age_breaks) != length(Susceptibility)){
    warning("Susceptibility not provided or inaccurate -> homogeneous susceptibility considered")
    Susceptibility <- rep(1, length(Age_breaks))
  }
  if(length(Age_breaks) != length(Infectiousness)){
    warning("Infectiousness not provided or inaccurate -> homogeneous infectiousness considered")
    Infectiousness <- rep(1, length(Age_breaks))
  }
  
  if(N_bs>1){ # case with bootstrap

    wavematrix_bs <- vector(mode = "list", length = Wavemax)
    comix_bs <- vector(mode = "list", length = Wavemax)
    
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
   if(rds=="no"){
       matrix_out_bs <- contact_matrix(survey_obj, 
                                    age.limits = Age_breaks,
                                    symmetric  = TRUE,
                                    quiet      = TRUE,
                                    weigh.dayofweek = TRUE,
                                    weigh.age = TRUE,
                                    weight.threshold = weight_threshold,
                                    estimated.contact.age = 'sample',
                                    missing.contact.age =  "sample",
                                    n=N_bs
    )
     } else {
     matrix_out_bs <- matrix_rds_bs[[i]]
     }
      wavematrix_bs[[i]] <- lapply(matrix_out_bs$matrices, function(x) {x$matrix[is.na(x$matrix)] <- 0; return(x)})
    wavematrix_bs[[i]] <- lapply(wavematrix_bs[[i]], function(x) {x$matrix <- diag(Susceptibility) %*% t(x$matrix) %*% diag(Infectiousness)  ; return(x)})
    comix_bs[[i]]$bs <- lapply(wavematrix_bs[[i]],function(x) {data.frame(Age_breaks,abs(eigen(x$matrix)$vectors[,1]))})
    comix_bs[[i]]$bs <- lapply(comix_bs[[i]]$bs, setNames, c("ageclass","nextgen"))
    comix_bs[[i]]$bs <- lapply(comix_bs[[i]]$bs,function(x) {mutate(x,nextgen=nextgen *100 / sum(x$nextgen))  })
    comix_bs[[i]] <-  bind_rows(comix_bs[[i]]$bs)  %>% group_by(ageclass) %>% summarise(nextgenmean = mean(nextgen), nextgenmedian = median(nextgen),  nextgensd = sd(nextgen), nextgenlower = unname(quantile(nextgen,0.025)), nextgenupper = unname(quantile(nextgen,0.975)), .groups = 'drop')
    comix_bs[[i]]$wave <- i
    comix_bs[[i]]$wavename <- paste("Wave",i)
    comix_bs[[i]]$ageclassname <- factor(matrix_out_bs$participants$age.group, levels = matrix_out_bs$participants$age.group)
    comix_bs[[i]]$ageclasslower <- factor(matrix_out_bs$participants$lower.age.limit, levels = matrix_out_bs$participants$lower.age.limit)
  } # end for wave
    
    return(comix_bs)
    
  } else { # case without bootstrap

    wavematrix <- vector(mode = "list", length = Wavemax)
    comix <- vector(mode = "list", length = Wavemax)
    
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
    
    matrix_out <- contact_matrix(survey_obj, 
                                    age.limits = Age_breaks,
                                    symmetric  = TRUE,
                                    quiet      = TRUE,
                                    weigh.dayofweek = TRUE,
                                    weigh.age = TRUE,
                                    weight.threshold = weight_threshold,
                                    estimated.contact.age = 'sample',
                                 missing.contact.age =  "sample",
                                 n=1
    )
    wavematrix[[i]] <- t(matrix_out$matrix)
    wavematrix[[i]][is.na(wavematrix[[i]])] <- 0
    wavematrix[[i]] <-  diag(Susceptibility) %*% t(wavematrix[[i]]) %*% diag(Infectiousness) 
    comix[[i]] <- data.frame(Age_breaks,abs(eigen(wavematrix[[i]])$vectors[,1]))
    colnames(comix[[i]] ) <- c("ageclass","nextgenmean")
    comix[[i]]  <- comix[[i]]   %>% mutate(nextgenmean=nextgenmean *100 / sum(comix[[i]]$nextgenmean))
    comix[[i]]$nextgenmedian <- comix[[i]]$nextgenmean
    comix[[i]]$nextgensd <- 0
    comix[[i]]$nextgenlower <- comix[[i]]$nextgenmean
    comix[[i]]$nextgenupper <- comix[[i]]$nextgenmean
    comix[[i]]$wave <- i
    comix[[i]]$wavename <- paste("Wave",i)
    comix[[i]]$ageclassname <- factor(matrix_out$participants$age.group, levels = matrix_out$participants$age.group)
    comix[[i]]$ageclasslower <- factor(matrix_out$participants$lower.age.limit, levels = matrix_out$participants$lower.age.limit)
    comix[[i]] <- as_tibble(comix[[i]])
    } # end for wave
    
    return(comix)
    
  } # end with/without bootstrap

  
} # end function




