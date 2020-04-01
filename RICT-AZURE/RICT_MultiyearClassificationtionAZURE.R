#--- ########################################################
#--- title: "RICT IV WHPT MULTI-YEAR PREDICTION  version 2" #
#--- output: html_document                                  #
#--- Author: Dr K. M Muyeba aka Maybin                      #  
#--- ########################################################

############################################################################################################################
# This experiment is for multi-year classification during spring and autumn                                                #
# MultiYear Classification-uses Monte Carlo 10,000 simulation cycles to calculate probability of class - of the predictions#
# and MINTA - minimum taxa and classes                                                                                     #
# ##########################################################################################################################



#Part 1: This Script reads all prediction indices for classification 
 
Allpredictions <- maml.mapInputPort(1) # class: data.frame

#install.packages(c("cowplot", "plyr", "httr", "xml2", "stringr", "xts", "rjson", "ggmap", "ggplot2", "sp", "rgdal", "parallel")) # replaced plyr by dplyr
#install.packages("src/rnrfa_1.4.0.zip", lib = ".", repos = NULL, verbose = TRUE, dependencies = TRUE)

#head(Allpredictions)
library(dplyr)
library(magrittr)
library(gridExtra)
library(data.table)

source("src/ClassificationfunctionsV2.R")
 
 start_time <- proc.time()
# Ead big CSV files
GB685_Ass_score <-  readCSVFile ("src/EndGrp_AssessScores.csv") # read.csv("src/EndGrp_AssessScores.csv")
Aj <- readCSVFile ("src/adjustParams_ntaxa_aspt.csv") # read.csv("src/adjustParams_ntaxa_aspt.csv")

# Set the seed for repeatability of experiment
set.seed (1234) #(2345)

#Enter source files 
#Use the column header as site names in the final output 
allSites <- Allpredictions[,1]
#Keep YEAR, WATERBODY
year_waterBody <- Allpredictions[,c("YEAR","WATERBODY")]

#Combine allSites with more information  - e.g. YEAR, WATERBODY
allSites <- cbind(allSites, year_waterBody)

#Allpredictions <- Allpredictions[,-1]
# Create a routine that checks column names and puts them in the "Allpredictions" dataframe
# Allpredictions <- Allpredictions[,-c(1,3:15)]
 # Change all names to upper case
 
names(Allpredictions) <- toupper(names(Allpredictions))
head(colnames(Allpredictions), 28)
#tail(colnames(Allpredictions), 28)

#Remove the "_CompFarm_" columns
Allpredictions <- select(Allpredictions, -matches("_COMPFAM_") ) # use the "-" with "match" from dplyr
#Get the biological data TL2_WHPT_NTAXA_AbW_DistFam_spr
namesBiological <-    c("SPR_SEASON_ID", "SPR_TL2_WHPT_ASPT (ABW,DISTFAM)","SPR_TL2_WHPT_NTAXA (ABW,DISTFAM)", "SPR_NTAXA_BIAS", "SUM_SEASON_ID", "SUM_TL2_WHPT_ASPT (ABW,DISTFAM)","SUM_TL2_WHPT_NTAXA (ABW,DISTFAM)", "SUM_NTAXA_BIAS", "AUT_SEASON_ID", "AUT_TL2_WHPT_ASPT (ABW,DISTFAM)","AUT_TL2_WHPT_NTAXA (ABW,DISTFAM)","AUT_NTAXA_BIAS")
 
biologicalData <- Allpredictions[,namesBiological]

# Remove biologicalData from Allpredictions
Allpredictions <- Allpredictions[,!names(Allpredictions) %in% namesBiological]
 
#Store allProbabilities in one dataframe. Use p1,p2,... etc in case data column positions change in future
probNames <- c("p1","p2","p3","p4","p5","p6","p7","p8","p9","p10","p11","p12","p13","p14","p15","p16","p17","p18","p19","p20","p21","p22","p23","p24","p25","p26","p27","p28","p29","p30","p31","p32","p33","p34","p35","p36","p37","p38","p39","p40","p41","p42","p43")

allProbabilities <- Allpredictions[,toupper(probNames)] # Needs to change when not uppercase
# Input Adjustment factors for reference site quality scores (Q1, Q2, Q3, Q4, Q5)

# Extract Ubias8 from Biological data # 
UBIAS_main <- biologicalData[,"SPR_NTAXA_BIAS"][1] # Put new AZURE

# Put the UBIAS_main default value of 1.68 if the user does not enter any value or entersa -9
if(is.na(UBIAS_main) | UBIAS_main==-9) { # For NI model, the default is ZERO
  UBIAS_main <- 1.68
}

# OBSERVED ASPT
#Obs_aspt_spr <- observed_aspt[,1]
Obs_aspt_spr    <- biologicalData[,"SPR_TL2_WHPT_ASPT (ABW,DISTFAM)"]
#Obs_aspt_aut <- observed_aspt[,2]
Obs_aspt_aut    <- biologicalData[,"AUT_TL2_WHPT_ASPT (ABW,DISTFAM)"] 

# OBSERVED NTAXA
#Obs_ntaxa_spr <- observed_ntaxa[,1]
Obs_ntaxa_spr    <- biologicalData[,"SPR_TL2_WHPT_NTAXA (ABW,DISTFAM)"]
#Obs_ntaxa_aut <- observed_ntaxa[,2]
Obs_ntaxa_aut    <- biologicalData[,"AUT_TL2_WHPT_NTAXA (ABW,DISTFAM)"] # change AZURE

# Input Multiplicative Adjustment factors Aj, 1,..,5)
Aj <- as.matrix(Aj)
Qij <- computeScoreProportions(as.data.frame(GB685_Ass_score[,-1])) # Remove the first Column
## Part 2:  Calculate AdjustedExpected from all probabilities, WE4.5 of WFD72C

# Compute Rj = sum(Pi*Qij)
Rj <- as.matrix(getWeighted_proportion_Rj(allProbabilities, Qij))# We should have five of these 

#Multiply Rj by Aj, note each row of Aj is for NTAXA, ASPT, so transpose to multiply by Rj
RjAj <- compute_RjAj(Rj, Aj)
One_over_RjAj <- 1/RjAj
#tail(One_over_RjAj,10)

# Write a function that computes aspt, ntaxa adjusted (1 = "NTAXA", 2="ASPT") or select them by name as declared in the classification functions 
ntaxa_Adjusted <- select(Allpredictions, matches("_NTAXA_")) / RjAj[,"NTAXA"]
aspt_Adjusted <- select(Allpredictions, matches("_ASPT_")) / RjAj[,"ASPT"] #Compute AdjExpected as E=Allpredictions/Sum(Rj*Aj)
 
Adjusted_Expected <- cbind(ntaxa_Adjusted, aspt_Adjusted)
Adjusted_Expected_new <- cbind(as.data.frame( allSites), Adjusted_Expected) # Include site names from Allpredictions

# Part 3:  Calculation of Exp_ref from "AdjustedExpected_new" values, divide by K ( = 1.0049 for NTAXA,  = 0.9921 for ASPT)
# run simulations from here 
N_runs <- 10000

# ******* FOR ASPT ************
Exp_ref_aspt  <- aspt_Adjusted/0.9921
Ubias8 <- UBIAS_main

# find the non-bias corrected  EQR = Obs/ExpRef
nonBiasCorrected_WHPT_aspt_spr <- Obs_aspt_spr/select(Exp_ref_aspt, matches("_spr"))
nonBiasCorrected_WHPT_aspt_aut <- Obs_aspt_aut/select(Exp_ref_aspt, matches("_aut"))

# Now do the Obs_rb withONE SITE Obs_aspt_spr[1]
sdobs_aspt <- SDObs_One_year_new(0.269, 0.279, 1)

EQRAverages_aspt_spr <- data.frame() # Store average EQRs for spr in a dataframe
EQRAverages_aspt_aut <- data.frame() # Store average EQRs for spr in a dataframe

# **************  For NTAXA   *************
Exp_ref_ntaxa <- ntaxa_Adjusted/1.0049 # select(Adjusted_Expected_new, matches("_NTAXA_"))/1.0049

# Find the non-bias corrected  EQR = Obs/ExpRef, from the raw inputs, not used but useful for output checking purposes only
nonBiasCorrected_WHPT_ntaxa_spr <- Obs_ntaxa_spr/select(Exp_ref_ntaxa, matches("_spr"))
nonBiasCorrected_WHPT_ntaxa_aut <- Obs_ntaxa_aut/select(Exp_ref_ntaxa, matches("_aut"))

# Now do the Obs_rb with ONE SITE Obs_ntaxa_spr[1]
sdobs_ntaxa <- SDObs_One_year_new(0.247, 0.211, 1)
#Define sdexp
sdexp8_ntaxa <- 0.53  
sdexp9_aspt <- 0.081
#Define for ASPT 
u_9a  <- 4.35 
u_9b <- 0.271 
u_9c <- 2.5

SiteProbabilityclasses_spr_ntaxa <- data.frame() # Store site probabilities in a dataframe
SiteProbabilityclasses_aut <- data.frame() # Store site probabilities in a dataframe
SiteProbabilityclasses_aut_ntaxa <- data.frame()
SiteProbabilityclasses_spr_aut_comb_ntaxa <- data.frame()
SiteProbabilityclasses_spr_aut_comb_aspt <- data.frame()

#MINTA
SiteMINTA_whpt_spr <- data.frame()
SiteMINTA_whpt_aut <- data.frame()
SiteMINTA_whpt_spr_aut <- data.frame()

# ASPT 
SiteProbabilityclasses_spr_aspt <- data.frame() # Store site probabilities in a dataframe
SiteProbabilityclasses_aut_aspt <- data.frame() # Store site probabilities in a dataframe
SiteProbabilityclasses_spr_aut_comb_aspt <- data.frame()
classArray_siteOne_spr_aut_ntaxa <- data.frame()
classArray_siteOne_spr_aut_aspt  <- data.frame()

#Setup biases
Ubias8r_spr <-  getUbias8r_new (N_runs, Ubias8)
Ubias8r_aut <-  getUbias8r_new (N_runs, Ubias8)
 print(head(Ubias8r_spr,3))
#Store all multiYear  
EQRAverages_ntaxa_spr_aut <- data.frame() # Store average EQRs for spr in a dataframe
EQRAverages_aspt_spr_aut  <- data.frame() # Store average EQRs for spr in a dataframe

# initalise all MultiYear 
multiYear_EQRAverages_ntaxa_spr <- data.frame(n=N_runs)
multiYear_EQRAverages_ntaxa_aut <- data.frame(n=N_runs)
multiYear_EQRAverages_ntaxa_spr_aut <- data.frame(n=N_runs) # Stores averages for nyears-use to calculate, find all spring, all autumn, then average these for each index

multiYear_EQRAverages_aspt_spr <- data.frame(n=N_runs)
multiYear_EQRAverages_aspt_aut <- data.frame(n=N_runs)
multiYear_EQRAverages_aspt_spr_aut <- data.frame(n=N_runs) # Stores averages for nyears-use to calculate, find all spring, all autumn, then average these for each index

multipleSite_encoutered <- FALSE
# Store the duplicated names of sites as single names
namesOfSites <- data.frame()

# Variable that flags if last site has not been processed
lastSiteProcessed <- FALSE
  
# Collection of indices 
indicesDistinct <- data.frame()
k<- 1 

while (k <=nrow(Allpredictions) | (lastSiteProcessed==FALSE)) { 
  # initalise all MultiYear AGAIN for each site
  multiYear_EQRAverages_ntaxa_spr <- data.frame(n=N_runs)
  multiYear_EQRAverages_ntaxa_aut <- data.frame(n=N_runs)
  multiYear_EQRAverages_ntaxa_spr_aut <- data.frame(n=N_runs) # Stores averages for nyears-use to calculate, find all spring, all autumn, then average these for each index

  multiYear_EQRAverages_aspt_spr <- data.frame(n=N_runs)
  multiYear_EQRAverages_aspt_aut <- data.frame(n=N_runs)
  multiYear_EQRAverages_aspt_spr_aut <- data.frame(n=N_runs) # Stores averages for nyears-use to calculate, find all spring, all autumn, then average these for each index

   #Declare a boolean variable that indicates multiple sites encountered, then switch it back to FALSE at start of loop
  j <- k
  print(c(" j = ",j))
  indicesDistinct <- rbind(indicesDistinct, j)
  
  if(j<nrow(Allpredictions) && (Allpredictions[j,"SITE"]==Allpredictions[j+1,"SITE"])){
        multipleSite_encoutered <- TRUE
  }
  
  #Get site out
  siteToProcess <- Allpredictions[j,"SITE"]
  while ( (Allpredictions[j,"SITE"]==siteToProcess  && j<= nrow(Allpredictions))){ 
  # print(c("Processing site j= ",j, " as ", as.character(Allpredictions[j,"SITE"]), " and site j= ",j+1," as ", as.character(Allpredictions[j+1,"SITE"])))
        
  # Part 1: Deal with NTAXA: Observed and Expcted Calculations
  ObsIDX8r_spr  <- getObsIDX8r(as.double(Obs_ntaxa_spr[j]),getZObs_r_new(sdobs_ntaxa,N_runs))
  ObsIDX8r_aut  <- getObsIDX8r(as.double(Obs_ntaxa_aut[j]),getZObs_r_new(sdobs_ntaxa,N_runs))
  Obs_site1_ntaxa_spr <- ObsIDX8r_spr + Ubias8r_spr # rename "Obs_site1_ntaxa_spr" to ObsIDX8rb_spr
  Obs_site1_ntaxa_aut <- ObsIDX8r_aut + Ubias8r_aut # rename "Obs_site1_ntaxa_aut" to ObsIDX8rb_aut
  # Part 2 . Do the RefAdjExpected bias
  ExpIDX8r_ntaxa_spr <- data.frame(val = (as.double(Exp_ref_ntaxa[j,1]) + getZObs_r_new (sdexp8_ntaxa, N_runs)))
  ExpIDX8r_ntaxa_aut <- data.frame(val = (as.double(Exp_ref_ntaxa[j,2]) + getZObs_r_new (sdexp8_ntaxa, N_runs)))
  
  EQR_ntaxa_spr <- as.data.frame(Obs_site1_ntaxa_spr/ExpIDX8r_ntaxa_spr[,1])
  EQR_ntaxa_aut <- as.data.frame(Obs_site1_ntaxa_aut/ExpIDX8r_ntaxa_aut[,1] )
  
  #Store these multi sites for ntaxa here
  multiYear_EQRAverages_ntaxa_spr <- cbind(multiYear_EQRAverages_ntaxa_spr,EQR_ntaxa_spr ) # Afterwards, use:  EQR_ntaxa_spr <- rowMeans(multiYear_EQRAverages_ntaxa_spr[,-1])
  multiYear_EQRAverages_ntaxa_aut <- cbind(multiYear_EQRAverages_ntaxa_aut,EQR_ntaxa_aut ) # Afterwardds, use  EQR_ntaxa_aut <- rowMeans(multiYear_EQRAverages_ntaxa_aut[,-1])
 
  # Part 1: Deal with ASPT: Observed and Expcted Calculations
  # ****************************************
  # **** Workout FOR ASPT STARTS HERE 
  # Part 1: Deal with ASPT : Observed and Expcted Calculations

  Ubias9r_spr <- getUbias9r_new (u_9a, u_9b,u_9c, as.double(Obs_aspt_spr[j]), N_runs, Ubias8r_spr)
  Ubias9r_aut <- getUbias9r_new (u_9a, u_9b, u_9c,as.double(Obs_aspt_aut[j]), N_runs, Ubias8r_aut)
  Ubias7r_spr <- Ubias8r_spr*Ubias9r_spr
  Ubias7r_aut <- Ubias8r_aut*Ubias9r_aut
  ObsIDX9r_spr  <- getObsIDX9r (as.double(Obs_aspt_spr[j]),getZObs_r_new(sdobs_aspt,N_runs)) 
  ObsIDX9r_aut  <- getObsIDX9r(as.double(Obs_aspt_aut[j]),getZObs_r_new(sdobs_aspt,N_runs))
  ObsIDX7r_spr <-  ObsIDX8r_spr* ObsIDX9r_spr
  ObsIDX7r_aut <-  ObsIDX8r_aut* ObsIDX9r_aut
  ObsIDX7rb_spr <- ObsIDX7r_spr+Ubias7r_spr
  ObsIDX7rb_aut <- ObsIDX7r_aut+Ubias7r_aut
  ObsIDX8rb_spr <- ObsIDX8r_spr+Ubias8r_spr
  ObsIDX8rb_aut <- ObsIDX8r_aut+Ubias8r_aut
  ObsIDX9rb_spr <- ObsIDX7rb_spr/ObsIDX8rb_spr
  ObsIDX9rb_aut <- ObsIDX7rb_aut/ObsIDX8rb_aut
  # Part 2 . Do the RefAdjExpected bias
  ExpIDX9r_aspt_spr <- data.frame(val = (as.double(Exp_ref_aspt[j,1]) + getZObs_r_new (sdexp9_aspt, N_runs)))
  ExpIDX9r_aspt_aut <- data.frame(val = (as.double(Exp_ref_aspt[j,2]) + getZObs_r_new (sdexp9_aspt, N_runs)))
  # Calculating simulated EQR
  EQR_aspt_spr <- data.frame(ObsIDX9rb_spr/ExpIDX9r_aspt_spr[,1])
  EQR_aspt_aut <- data.frame(ObsIDX9rb_aut/ExpIDX9r_aspt_aut[,1] )
  
  #Store these multi sites for aspt here
  multiYear_EQRAverages_aspt_spr <- cbind(multiYear_EQRAverages_aspt_spr,EQR_aspt_spr ) # Afterwards, use:  EQR_aspt_spr <- rowMeans(multiYear_EQRAverages_aspt_spr[,-1])
  multiYear_EQRAverages_aspt_aut <- cbind(multiYear_EQRAverages_aspt_aut,EQR_aspt_aut ) # Afterwards, use:  EQR_aspt_aut <- rowMeans(multiYear_EQRAverages_aspt_aut[,-1])
      j<- j+1
  }
  
  if(multipleSite_encoutered==FALSE){ 
    if(k==nrow(Allpredictions)){
      lastSiteProcessed <- TRUE
    }
  # Part 1: Deal with NTAXA: Observed and Expcted Calculations
  ObsIDX8r_spr  <- getObsIDX8r(as.double(Obs_ntaxa_spr[k]),getZObs_r_new(sdobs_ntaxa,N_runs))
  ObsIDX8r_aut  <- getObsIDX8r(as.double(Obs_ntaxa_aut[k]),getZObs_r_new(sdobs_ntaxa,N_runs))
  Obs_site1_ntaxa_spr <- ObsIDX8r_spr + Ubias8r_spr # rename "Obs_site1_ntaxa_spr" to ObsIDX8rb_spr
  Obs_site1_ntaxa_aut <- ObsIDX8r_aut + Ubias8r_aut # rename "Obs_site1_ntaxa_aut" to ObsIDX8rb_aut
  # Part 2 . Do the RefAdjExpected bias
  ExpIDX8r_ntaxa_spr <- data.frame(val = (as.double(Exp_ref_ntaxa[k,1]) + getZObs_r_new (sdexp8_ntaxa, N_runs)))
  ExpIDX8r_ntaxa_aut <- data.frame(val = (as.double(Exp_ref_ntaxa[k,2]) + getZObs_r_new (sdexp8_ntaxa, N_runs)))
  
  EQR_ntaxa_spr <- as.data.frame(Obs_site1_ntaxa_spr/ExpIDX8r_ntaxa_spr[,1])
  EQR_ntaxa_aut <- as.data.frame(Obs_site1_ntaxa_aut/ExpIDX8r_ntaxa_aut[,1] )
  
  # Part 1: Deal with ASPT: Observed and Expcted Calculations
  # ****************************************
  # **** Workout FOR ASPT STARTS HERE 
  # Part 1: Deal with ASPT : Observed and Expcted Calculations

  Ubias9r_spr <- getUbias9r_new (u_9a, u_9b, u_9c,as.double(Obs_aspt_spr[k]), N_runs, Ubias8r_spr)
  Ubias9r_aut <- getUbias9r_new (u_9a, u_9b, u_9c,as.double(Obs_aspt_aut[k]), N_runs, Ubias8r_aut)
  Ubias7r_spr <- Ubias8r_spr*Ubias9r_spr
  Ubias7r_aut <- Ubias8r_aut*Ubias9r_aut
  ObsIDX9r_spr  <- getObsIDX9r (as.double(Obs_aspt_spr[k]),getZObs_r_new(sdobs_aspt,N_runs)) 
  ObsIDX9r_aut  <- getObsIDX9r(as.double(Obs_aspt_aut[k]),getZObs_r_new(sdobs_aspt,N_runs))
  ObsIDX7r_spr <-  ObsIDX8r_spr* ObsIDX9r_spr
  ObsIDX7r_aut <-  ObsIDX8r_aut* ObsIDX9r_aut
  ObsIDX7rb_spr <- ObsIDX7r_spr+Ubias7r_spr
  ObsIDX7rb_aut <- ObsIDX7r_aut+Ubias7r_aut
  ObsIDX8rb_spr <- ObsIDX8r_spr+Ubias8r_spr
  ObsIDX8rb_aut <- ObsIDX8r_aut+Ubias8r_aut
  ObsIDX9rb_spr <- ObsIDX7rb_spr/ObsIDX8rb_spr
  ObsIDX9rb_aut <- ObsIDX7rb_aut/ObsIDX8rb_aut
  # Part 2 . Do the RefAdjExpected bias
  ExpIDX9r_aspt_spr <- data.frame(val = (as.double(Exp_ref_aspt[k,1])+ getZObs_r_new (sdexp9_aspt, N_runs)))
  ExpIDX9r_aspt_aut <- data.frame(val = (as.double(Exp_ref_aspt[k,2])+ getZObs_r_new (sdexp9_aspt, N_runs)))
  # Calculating simulated EQR
  EQR_aspt_spr <- as.data.frame(ObsIDX9rb_spr/ExpIDX9r_aspt_spr[,1])
  EQR_aspt_aut <- as.data.frame(ObsIDX9rb_aut/ExpIDX9r_aspt_aut[,1] )
  }
  
  if(multipleSite_encoutered==TRUE) {
      if( (j==nrow(Allpredictions)) | (j-1==nrow(Allpredictions))) { # Means last site was a duplicate, and so is processed
          lastSiteProcessed<- TRUE
      }
      # Move to current record just processed
      j<- j-1 
      # ******************************************
      # Part 1.1: for "Spring" - DO FOR NTAXA 
      # Combined ntaxa spr-aut
    
      EQR_ntaxa_spr <- data.frame(EQR_ntaxa_spr = rowMeans(data.frame(multiYear_EQRAverages_ntaxa_spr[,-1])))
      EQR_ntaxa_aut <- data.frame(EQR_ntaxa_aut = rowMeans(data.frame(multiYear_EQRAverages_ntaxa_aut[,-1])))
    
      # ******************************************
      # Part 1.1: ASPT for "Spring"
      # Find the averages of both spr and autum, declare a function to compute this
      # First find all rowMeans, and store them in EQR appropriate variables
    
      EQR_aspt_spr <- data.frame(EQR_aspt_spr =  rowMeans(data.frame(multiYear_EQRAverages_aspt_spr[,-1])))
      EQR_aspt_aut <- data.frame(EQR_aspt_aut =  rowMeans(data.frame(multiYear_EQRAverages_aspt_aut[,-1])))
  }
  
  # Calculate EQRs here, i.e. rowSums if multipleTrue else just getAvgEQR() for single season
  eqr_av_spr <- data.frame(rowMeans(getAvgEQR_SprAut (EQR_ntaxa_spr,EQR_ntaxa_aut)))
  eqr_av_spr_aspt <- data.frame(rowMeans(getAvgEQR_SprAut (EQR_aspt_spr,EQR_aspt_aut)))
  
  # START TO CALCULATE probability of class
  # Part 2: Start calculating for NTAXA probability of CLASS
  
  multiYear_EQRAverages_ntaxa_spr_aut <- data.frame(rbind(cbind(EQR_ntaxa_spr,EQR_ntaxa_aut)))
  multiYear_EQRAverages_ntaxa_spr_aut <- data.frame(EQR_ntax_aspr_aut=rowMeans(multiYear_EQRAverages_ntaxa_spr_aut))
  multiYear_EQRAverages_aspt_spr_aut  <- data.frame(rbind(cbind(EQR_aspt_spr,EQR_aspt_aut)))
  multiYear_EQRAverages_aspt_spr_aut  <- data.frame(EQR_aspt_spr_aut=rowMeans(multiYear_EQRAverages_aspt_spr_aut))
    
  classArray_siteOne_spr_aut_ntaxa <- getClassarray_ntaxa(multiYear_EQRAverages_ntaxa_spr_aut) #data.frame(EQR_ntaxa_spr))
  
  # define an array to hold probability of class for each site- how much of the site belongs to each classes, adds up to 100%
  
  probClass_spr <- matrix(0, ncol = 1, nrow = 5) # 5 is the number of classes- H, G, M, B, P, ncol=1 or 2 for two seasons or ntaxa_spr, ntaxa_aut, spr_aut_av_taxa, and spt etc
   
  for(i in 1:5) {
    probClass_spr[i] <- 100*sum(classArray_siteOne_spr_aut_ntaxa[classArray_siteOne_spr_aut_ntaxa==i,]/i)/N_runs
  }
  
  # Part 2.1: for Spring_aut
  probabilityClass <- getProbClassLabelFromEQR()
  a_ntaxa_spr_aut <- t(probClass_spr) # spr, need a_ntaxa_spr
  colnames(a_ntaxa_spr_aut) <- getProbClassLabelFromEQR()[,1]
  rownames(a_ntaxa_spr_aut) <- as.character(Allpredictions[j,"SITE"]) #c(paste0("TST-",j))
  
  #Find most probable class, i.e the maximum, and add it to the site
  mostProb <- getMostProbableClass(a_ntaxa_spr_aut)
  a_ntaxa_spr_aut <- data.frame(cbind(a_ntaxa_spr_aut, mostProb)) # add the site to the dataframe
  SiteProbabilityclasses_spr_aut_comb_ntaxa<- rbind(SiteProbabilityclasses_spr_aut_comb_ntaxa,a_ntaxa_spr_aut)
  
  #Add the averages of spr,aut
  EQRAverages_ntaxa_spr_aut <- rbind(EQRAverages_ntaxa_spr_aut, eqr_av_spr)
  
  ######## **************** Now do the ASPT from HERE - using the calculations from ASPT ABOVE*********************
  
  # Part 2: Start calculating for ASPT probability of CLASS
  # Classify these for each SITE using the EQR just for spring
  classArray_siteOne_spr_aut_aspt  <- getClassarray_aspt(multiYear_EQRAverages_aspt_spr_aut) #data.frame(EQR_ntaxa_aut))
  # define an array to hold probability of class for each site- how much of the site belongs to each classes, adds up to 100%
    
  probClass_spr <- matrix(0, ncol = 1, nrow = 5) # 5 is the number of classes- H, G, M, B, P, ncol=1 or 2 for two seasons or ntaxa_spr, ntaxa_aut, spr_aut_av_taxa, and spt etc
   
  for(i in 1:5) {
    probClass_spr[i] <- 100*sum(classArray_siteOne_spr_aut_aspt[classArray_siteOne_spr_aut_aspt==i,]/i)/N_runs
  }
    
  # Work out ASPT probability of classes
  #probabilityClass <- getProbClassLabelFromEQR()
  a_aspt_spr_aut <- t(probClass_spr) # spr
  colnames(a_aspt_spr_aut) <- getProbClassLabelFromEQR()[,1]
  #print(c(" j =",j," site = ",as.character(Allpredictions[j,"SITE"]), "pated TST = ",paste0("TST-",j)))
  rownames(a_aspt_spr_aut) <- as.character(Allpredictions[j,"SITE"]) #c(paste0("TST-",j))
  
  #Find most probable class, i.e the maximum, and add it to the site
  mostProb <- getMostProbableClass(a_aspt_spr_aut)
  # add the site to the dataframe
  a_aspt_spr_aut <- data.frame(cbind(a_aspt_spr_aut, mostProb))
  SiteProbabilityclasses_spr_aut_comb_aspt<- rbind(SiteProbabilityclasses_spr_aut_comb_aspt,a_aspt_spr_aut)
  #Add the averages of spr
  EQRAverages_aspt_spr_aut <- rbind(EQRAverages_aspt_spr_aut, eqr_av_spr_aspt)
  
  ########  Calculate the MINTA -spring aut case  worse class = 1 i.e. min of class from NTAXA and ASPT ######
    
  # Do the MINTA spr_aut case 
  minta_ntaxa_aspt_spr_aut <- getMINTA_ntaxa_aspt (as.matrix(classArray_siteOne_spr_aut_ntaxa),  as.matrix(classArray_siteOne_spr_aut_aspt))
  minta_probClass_spr_aut <- matrix(0, ncol = 1, nrow = 5) # 5 is the number of classes- H, G, M, B, P, ncol=1 or 2 for two seasons or ntaxa_spr, ntaxa_aut, spr_aut_av_taxa, and spt etc
   
  for(i in 1:5) {
      minta_probClass_spr_aut[i] <- 100*sum(minta_ntaxa_aspt_spr_aut[minta_ntaxa_aspt_spr_aut==i,]/i)/N_runs
  }
    
  #probabilityClass <- getProbClassLabelFromEQR()
  aa <- t(minta_probClass_spr_aut) # spr
  colnames(aa) <- getProbClassLabelFromEQR()[,1]
  rownames(aa) <- as.character(Allpredictions[j,"SITE"]) #c(paste0("TST-",j))
  #Find most probable MINTA class, i.e the maximum, and add it to the site
  mostProb <- getMostProbableClass(aa)
  aa <- data.frame(cbind(aa, mostProb))
  # Now bind the MINTA proportion to the dataframe
  SiteMINTA_whpt_spr_aut <- rbind(SiteMINTA_whpt_spr_aut, aa) # Error in match.names(clabs, names(xi)) : 
  ##### MINTA ENDS HERE  #####
  
  # Move the pointer k to new adjusted position for j - whether multiple or not
   k <- j+1
    
}# END of FOR LOOP

# MINTA outputs
 
colnames(SiteMINTA_whpt_spr_aut) <- c(paste0("mintawhpt_spr_aut_",names(SiteMINTA_whpt_spr_aut)))
# Combine all MINTA

allMINTA_whpt <- SiteMINTA_whpt_spr_aut 

# DO for NTAXA 
colnames(EQRAverages_ntaxa_spr_aut) <- c(paste0("NTAXA_",colnames(EQRAverages_ntaxa_spr_aut)))
whpt_ntaxa_spr_aut_averages <- data.frame(NTAXA_aver_spr_aut=rowMeans(EQRAverages_ntaxa_spr_aut))
 #Rename column names so they dont conflict
colnames(SiteProbabilityclasses_spr_aut_comb_ntaxa) <- paste0(colnames(SiteProbabilityclasses_spr_aut_comb_ntaxa), "_NTAXA_spr_aut")
#Colname "mostProb" doesnt appear in SiteProbabilityclasses_spr_ntaxa, so reassign them here again - BUG in ML Studio R version

# ****** FOR ASPT outputs ********
colnames(EQRAverages_aspt_spr_aut) <- c(paste0("ASPT_",colnames(EQRAverages_aspt_spr_aut)))
whpt_aspt_spr_aut_averages <- data.frame(ASPT_aver_spr_aut=rowMeans(EQRAverages_aspt_spr_aut))
#Rename column names so they dont conflict
colnames(SiteProbabilityclasses_spr_aut_comb_aspt) <- paste0(colnames(SiteProbabilityclasses_spr_aut_comb_aspt), "_ASPT_spr_aut")

### DO FOR ALL including MINTA
# Bind the NTAXA
allResults <- cbind(SiteProbabilityclasses_spr_aut_comb_ntaxa,whpt_ntaxa_spr_aut_averages)
# Bind the ASPT
allResults <- cbind(allResults, cbind(SiteProbabilityclasses_spr_aut_comb_aspt,whpt_aspt_spr_aut_averages))
# Change names of Sites
namesOfSites <- data.frame(SITE=Allpredictions[unlist(indicesDistinct),"SITE"])
# Change waterbody to correct number of sites left after removing duplicate sites, do a collection of indices
year_waterBody <- year_waterBody[unlist(indicesDistinct),]

#Bind waterBody, and namesOfSites
allResults <- cbind(year_waterBody,allResults)
allResults <-  cbind(namesOfSites,allResults)

#Bind MINTA
# Rename columns for MINTA, so they dont conflict
colnames(SiteMINTA_whpt_spr_aut) <- paste0(colnames(SiteMINTA_whpt_spr_aut), "_MINTA_")
allResults_ntaxa_aspt_minta_combined <- cbind(allResults, SiteMINTA_whpt_spr_aut)


print("Time = ")
print(proc.time()-start_time)
maml.mapOutputPort("allResults_ntaxa_aspt_minta_combined");
