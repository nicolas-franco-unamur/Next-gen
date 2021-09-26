###############################################################################
# Function for next generation approach on CoMix data
# append pcr positive tests to nextgen output 
# Nicolas Franco - UHasselt
# Version: Augustus 4 2021
###############################################################################

  
#libraries
library(dplyr)

append_pcr <- function(comix_nextgen,
                       Wavemin = 9,
                       Wavemax = 9,
                       path_to_wd = ".",
                       PCR_file = "COVID19BE_AGE_POSTCODE_SEX_LAB.csv",
                       Age_breaks = c(18,30,40,50,60),
                       Delay = 7,
                       Period_length = 14,
                       wavedate,
                       Agemax = 120
                       ) {

  Agemin <- Age_breaks[1] #to exclude some ages
  
  setwd(path_to_wd)


#starting dates of coMix survey
if(missing(wavedate)){
  wavedate <- vector(mode = "list", length = 100)
  wavedate[[1]] <- "2020-04-20"
  wavedate[[2]] <- "2020-05-04"
  wavedate[[3]] <- "2020-05-18"
  wavedate[[4]] <- "2020-06-01"
  wavedate[[5]] <- "2020-06-15"
  wavedate[[6]] <- "2020-06-29"
  wavedate[[7]] <- "2020-07-13"
  wavedate[[8]] <- "2020-07-27"
  wavedate[[9]] <- "2020-11-11"
  wavedate[[10]] <- "2020-11-25"
  wavedate[[11]] <- "2020-12-09"
  wavedate[[12]] <- "2020-12-22"
  wavedate[[13]] <- "2021-01-06"
  wavedate[[14]] <- "2021-01-20"
  wavedate[[15]] <- "2021-02-03"
  wavedate[[16]] <- "2021-02-17"
  wavedate[[17]] <- "2021-03-03"
  wavedate[[18]] <- "2021-03-17"
  wavedate[[19]] <- "2021-04-01"
  wavedate[[20]] <- "2021-04-15"
  wavedate[[21]] <- "2021-04-29"
  wavedate[[22]] <- "2021-05-13"
  wavedate[[23]] <- "2021-05-27"
  wavedate[[24]] <- "2021-06-10"
  wavedate[[25]] <- "2021-06-24"}

  
  if(PCR_file != 'RDS'){   #using Sciensano non aggregate file
    testsDatabrut=read.csv(PCR_file)
    testsDatabrut$dateusedforstatistics=as.Date(testsDatabrut$dateusedforstatistics,"%d%b%Y")-Delay
    testscount <- vector(mode = "list", length = Wavemax)
    for(i in Wavemin:Wavemax){
      if(is.null(wavedate[[i]])){next}
     testscount[[i]] <- filter(testsDatabrut,dateusedforstatistics>=as.Date(wavedate[[i]])&dateusedforstatistics<as.Date(wavedate[[i]])+Period_length&Age>=Agemin&Age<=Agemax)
     testscount[[i]] <- testscount[[i]] %>% count(ageclass = cut(Age, breaks = c(Age_breaks, Inf),include.lowest = TRUE,right=FALSE))  
     testscount[[i]] <- testscount[[i]] %>% rename(percentage = n)
     testscount[[i]] <- testscount[[i]] %>% mutate(percentage=percentage *100 / sum(testscount[[i]]$percentage)) 
     comix_nextgen[[i]]$pcr <- testscount[[i]]$percentage
    }
  #saveRDS(testscount, file = "PCR_REAGGREGATE_DATA.rds")
    
  } else {   #using re-aggregate PCR_REAGGREGATE_DATA.rds file
    testscount <- readRDS(file = "PCR_REAGGREGATE_DATA.rds")
    for(i in Wavemin:Wavemax){
      if(is.null(wavedate[[i]])){next}
      testscount[[i]] <- testscount[[i]] %>% rename(percentage = n)
      testscount[[i]] <- testscount[[i]] %>% mutate(percentage=percentage *100 / sum(testscount[[i]]$percentage)) 
      comix_nextgen[[i]]$pcr <- testscount[[i]]$percentage
    }
  }
  
  return(comix_nextgen)
  
} #end function



