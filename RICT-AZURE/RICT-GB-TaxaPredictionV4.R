#--- ###########################################
#--- title: "RICT IV GB TAXA Prediction version 4"
#--- output: html_document
#--- Author: Dr K. M Muyeba aka Maybin 
#--- #########################################

  ####################################################################################################################################
  # This R Script is a TAXA prediction engine for FOR SITES. A user will usually choose which sites 1,2,4,8 or 2,5,9 etc to predict
  # If the user does not choose any sites, the program will run all sites in the file
  #################################################################################################################################### 
  
start_time <- proc.time()
print("Start-time ")
 
raw.input.data <- maml.mapInputPort(1) # class: data.frame
chooseSite <- as.double(maml.mapInputPort(2)) # class: data.frame

print(chooseSite)

#install.packages(c("cowplot", "plyr", "httr", "xml2", "stringr", "xts", "rjson", "ggmap", "ggplot2", "sp", "rgdal", "parallel"))
install.packages("src/rnrfa_1.4.0.zip", lib = ".", repos = NULL, verbose = TRUE, dependencies = TRUE)
library(rnrfa, lib.loc=".", verbose=TRUE)

library(dplyr) 
library(magrittr)
library(gridExtra)
library(data.table)

source("src/PredictionfunctionsV2.R") #PredictionfunctionsV1.R")
source("src/Helperfunctionsv1.R")
source("src/MeanAirTempAirTempRangeASFunction_new.R") #MeanAirTempAirTempRangeASFunction.R

model <- c("RIVPACS IV GB") # Changed in AZURE/RSTUDIO
  
DFMean_gd <- read.delim("src/DFMEAN_GB685.DAT", header = FALSE, sep="", as.is=TRUE)
DFCoeff_gb685 <- read.delim("src/DFCOEFF_GB685.DAT", header = FALSE, sep="", as.is=TRUE)

NRefg_groups <- readCSVFile("src/EndGrp_AssessScores.csv")
AirTempGrid <- readCSVFile("src/AirTempGrid.csv")
 
#Rename the column 1 from Ã¯..SITE to "SITE"
colnames(raw.input.data)[1]  <- c("SITE") # change AZURE
raw.input.data$SITE <- as.character(raw.input.data$SITE) # change AZURE
 
# check for length <5, add a "0" to get proper Easting/Northing 5 digit codes
raw.input.data$Easting  <- getCorrectCodes(raw.input.data$Easting)   # Change AZURE
raw.input.data$Northing <- getCorrectCodes(raw.input.data$Northing)  # Change AZURE

# Change all column names to uppercase
names(raw.input.data) <- toupper(names(raw.input.data)) # 

# Filter out those sites chosen by the user
raw.input.data <- if(length(chooseSite)>0) raw.input.data[chooseSite,] else (raw.input.data)
#ifelse(length(chooseSite)>0, raw.input.data[chooseSite,], raw.input.data)

print("New rows of data = ")
print(nrow(raw.input.data))

# Get all the bioligical data
namesBiological <-    c(colnames(raw.input.data)[1],colnames(raw.input.data)[2],colnames(raw.input.data)[3],"SPR_SEASON_ID", "SPR_TL2_WHPT_ASPT (ABW,DISTFAM)","SPR_TL2_WHPT_NTAXA (ABW,DISTFAM)", "SPR_NTAXA_BIAS", "SUM_SEASON_ID", "SUM_TL2_WHPT_ASPT (ABW,DISTFAM)","SUM_TL2_WHPT_NTAXA (ABW,DISTFAM)", "SUM_NTAXA_BIAS", "AUT_SEASON_ID", "AUT_TL2_WHPT_ASPT (ABW,DISTFAM)","AUT_TL2_WHPT_NTAXA (ABW,DISTFAM)","AUT_NTAXA_BIAS")
     
biologicalData <- raw.input.data[,namesBiological]

# Choose the seasons to run 
SEASONS_TO_RUN <- c(raw.input.data$SUM_SEASON_ID[1]) 

# Data validation and conversion
# 0. Validation of Various Env. variables before transformation

raw.input.data <- validateEnvData (raw.input.data) # Change in AZURE

# # 1. MEAN_WIDTH, lower_bound=0.4, upper_bound=117
 
 valid_mean_width <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
   valid_mean_width <- rbind(valid_mean_width,getValidEnvInput(raw.input.data$MEAN_WIDTH[i], 0.4, 117, "MEAN_WIDTH"))
 }

# # Change column names to suit env variable name, and cbind to original dataset
 colnames (valid_mean_width) <- paste0("mn_width_",noquote(colnames(valid_mean_width)))
 #head(valid_mean_width,5)
 raw.input.data <- cbind(raw.input.data, valid_mean_width)
 #head(raw.input.data, 3)
 
# # Data validation
# # 2. MEAN_DEPTH, lower_bound=1.7, upper_bound=999
 
 valid_mean_depth <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
   valid_mean_depth <- rbind(valid_mean_depth,getValidEnvInput(raw.input.data$MEAN_DEPTH [i], 1.7, 300, "MEAN_DEPTH"))
 }
 colnames (valid_mean_depth) <- paste0("mn_depth_",noquote(colnames(valid_mean_depth)))
 raw.input.data <- cbind(raw.input.data, valid_mean_depth)
 
# # Data validation
# # 3. SLOPE, lower_bound=0.1, upper_bound=150
 
 valid_slope <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
   valid_slope <- rbind(valid_slope,getValidEnvInput(raw.input.data$SLOPE [i], 0.1, 150, "SLOPE"))
 }
 colnames (valid_slope) <- paste0("vld_slope_",noquote(colnames(valid_slope))) # vld = valid
 raw.input.data <- cbind(raw.input.data, valid_slope)
 
# # Data validation
# # 4. DIST_FROM_SOURCE, lower_bound=0.1, upper_bound=999
# 
 valid_dist_src <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
   valid_dist_src <- rbind(valid_dist_src,getValidEnvInput(raw.input.data$DIST_FROM_SOURCE [i], 0.1, 202.8, "DIST_FROM_SOURCE"))
 }
 colnames (valid_dist_src) <- paste0("vld_dist_src_",noquote(colnames(valid_dist_src))) # vld = valid
 raw.input.data <- cbind(raw.input.data, valid_dist_src)
 
# # Data validation
# # 5. ALTITUDE, has two sets of bounds, lower_bound=1, upper_bound=590, lower_low_bound=0, upper_up_bound = 1345
# #[0,1345] are hard coded, could be parameterised QED
#  
  
 valid_altitude <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
    valid_altitude   <- rbind(valid_altitude,getEnvVariableTwoBounds(raw.input.data$ALTITUDE [i], 1, 590, 0, 1345, "ALTITUDE"))
   #valid_altitude <- rbind(valid_altitude,getAltitude(raw.input.data$ALTITUDE [i], 1, 590)) 
 }
 colnames (valid_altitude) <- paste0("vld_alt_src_",noquote(colnames(valid_altitude))) # vld = valid
 raw.input.data <- cbind(raw.input.data, valid_altitude)

# 
# # Data validation
# # 6. ALKALINITY, has bounds, lower_bound=1.2, upper_bound=999
# # getLogAlkalinity <- function (hardness, calcium, conduct, alkal, lower_b, upper_b)
# 
 valid_alkalinity <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
   valid_alkalinity <- rbind(valid_alkalinity,getLogAlkalinity(raw.input.data$HARDNESS[i], raw.input.data$CALCIUM[i], raw.input.data$CONDUCTIVITY[i],raw.input.data$ALKALINITY[i], 1.2, 366))
 }
 

# Above loop same as # thiscopy <- as.data.frame(with(raw.input.data, mapply(getLogAlkalinity, HARDNESS, CALCIUM, CONDUCTIVITY, ALKALINITY, 1.2, 999)))
#  thiscopy$V1$msg=="Succ"
colnames (valid_alkalinity) <- paste0("vld_alkal_",noquote(colnames(valid_alkalinity))) # vld = valid
raw.input.data <- cbind(raw.input.data, valid_alkalinity)
 
# # Data validation
# # 7. Validate SUBSTRATUM for sum of values "TOTSUB" in interval [97,103] exclussive,and MSUBSTR in interval [-8, 8]. Write to a file if errors found
# # Remove the site or records with such errors, and continue the prediction
# 
# # getSubstrate <- function(bould_cob, pebbles_gr, snd, silt_cl, lower_b, upper_b) 
# 

valid_substrate <- data.frame(log=as.numeric(), msg=as.character()) # Note that we don't use log for calculation of substrate
 for(i in 1:nrow(raw.input.data)){
   valid_substrate <- rbind(valid_substrate,getSubstrate (raw.input.data$BOULDER_COBBLES[i], raw.input.data$PEBBLES_GRAVEL[i], raw.input.data$SAND[i], raw.input.data$SILT_CLAY[i], 97, 103))
 }
 colnames (valid_substrate) <- paste0("vld_substr_",noquote(colnames(valid_substrate))) # vld = valid
 raw.input.data <- cbind(raw.input.data, valid_substrate)
 
#  #raw.input.data %>%    
#  #  subset(total>=97 & total<=103) %>%
# #    select(-ends_with("total")) # Remove the column "total"
#  
# # Data validation and conversion
# # 8. Discharge category, bounds [0, 10]. Discharge calculated from velocity if not provided using width, depth
# 
 valid_discharge <- data.frame(log=as.numeric(), msg=as.character())
 for(i in 1:nrow(raw.input.data)){
   valid_discharge <- rbind(valid_discharge, getLogDischarge(raw.input.data$MEAN_DEPTH[i], raw.input.data$MEAN_WIDTH[i], raw.input.data$DISCHARGE [i], raw.input.data$VELOCITY[i],0, 10, 1,9))
 }
 colnames (valid_discharge) <- paste0("disch_",noquote(colnames(valid_discharge)))
 raw.input.data <- cbind(raw.input.data, valid_discharge)
  
# 
# # Data validation and conversion
# # 9. Calculation of Lat/Long, and validation of LAT, LONG
# 
# # Calculation of Lat/Long using BNG (British National Grids)
# concatenatedNGR <- with(raw.input.data, paste(NGR, substr(Easting,1,3), substr(Northing,1,3), sep=""))
 
 
# #Use function getLatLong()
#coord_system = c("BNG", "WGS84")
#coordsystem = "WGS84"
       # a <- paste(raw.input.data$NGR, substr(raw.input.data$EASTING,1,3), substr(raw.input.data$northing,1,3), sep="") 
       # head(a, 2)
       #lat.long <- osg_parse (a, "WGS84")
 
  
lat.long <- with(raw.input.data, getLatLong_AZURE(NGR,EASTING, NORTHING, "WGS84")) #"WGS84") )
# head(lat.long, 1)

#### Calculate Longitude #####
raw.input.data$LONGITUDE <- lat.long$lon
 #print(c("lat.long = ",lat.long))
valid_longitude <- data.frame(log=as.numeric(), msg=as.character())
for(i in 1:nrow(raw.input.data)){
   valid_longitude <- rbind(valid_longitude,getLongitude(raw.input.data$LONGITUDE [i], -9, 1.5)) 
}
colnames (valid_longitude) <- paste0("vld_long_src_",noquote(colnames(valid_longitude))) # vld = valid
raw.input.data <- cbind(raw.input.data, valid_longitude)

# head(valid_longitude)
#### Calculate Latitude #####
raw.input.data$LATITUDE <- lat.long$lat
# print(c("lat.long = ",lat.long))
valid_latitude <- data.frame(log=as.numeric(), msg=as.character())
for(i in 1:nrow(raw.input.data)){
   valid_latitude <- rbind(valid_latitude,getLatitude(raw.input.data$LATITUDE [i], 48, 61)) 
}
colnames (valid_latitude) <- paste0("vld_lat_src_",noquote(colnames(valid_latitude))) # vld = valid
# head(raw.input.data,4)
raw.input.data <- cbind(raw.input.data, valid_latitude)

# # Data validation and conversion
# # 10. Calculation of mean temperature (TMEAN), range temperature (TRANGE), using function calc.temps() from package "rnfra"

# # Use function getBNG()
BNG <- with(raw.input.data, getBNG_AZURE(NGR,EASTING, NORTHING, "BNG") ) 

# head(BNG,4)

# # Lat long used for temperature lookups, using source MeanAirTempAirTempRangeASFunction.R
# Uncomment::
my.temperatures <- calc.temps(data.frame(
   Site_ID = raw.input.data$SITE,
   Easting4 = BNG$easting/100,
   Northing4 = BNG$northing/100,
   stringsAsFactors = FALSE))
# Assign to variables as appropriate
# Uncomment:: 
raw.input.data$TMEAN <- my.temperatures$TMEAN
# Uncomment::
raw.input.data$TRANGE <- my.temperatures$TRANGE


# # Data validation and conversion
# # 12. Write to file all Warnings and Failrures: SITE, MSG, iterate through the list of all variables with vld
# 
 
# # WRITE TO LOG FILES all Warnings and Errors 
# # 1. Warnings to log file :1
# 
# # Deal with all warnings, save them in a file
# #Same as above, but using pipes, and using all the variables
 msg_columns <- names(select(raw.input.data, ends_with("_msg")))
 this_warning <- raw.input.data %>%    
          filter(substr(vld_alt_src_msg,1,5)=="Warn:"    | substr(mn_width_msg,1,5)=="Warn:" 
                 | substr(mn_depth_msg,1,5)=="Warn:"     | substr(vld_alkal_msg,1,5)=="Warn:"
                 | substr(disch_msg,1,5)=="Warn:"        | substr(vld_substr_msg,1,5)=="Warn:"
                 | substr(vld_dist_src_msg,1,5)=="Warn:" | substr(vld_slope_msg,1,5)=="Warn:" ) 
        #  select("SITE","YEAR",msg_columns) # Select some columns
# # write.csv(this_warning, file = paste0(path,"/Warnings_file_data.csv"))
# # which rows are these 
# # raw.input.data[which(this_warning[1,1] %in% raw.input.data[,c("SITE")]),]

#head(this_warning, 18) # OK

# #2. Failings to log file
# 
# # Deal with all failings, save them in a file
  this_failing <- raw.input.data %>%    
          filter(substr(vld_alt_src_msg,1,5)=="Fail:"    | substr(mn_width_msg,1,5)=="Fail:" 
                 | substr(mn_depth_msg,1,5)=="Fail:"     | substr(vld_alkal_msg,1,5)=="Fail:"
                 | substr(disch_msg,1,5)=="Fail:"        | substr(vld_substr_msg,1,5)=="Fail:"
                 | substr(vld_dist_src_msg,1,5)=="Fail:" | substr(vld_slope_msg,1,5)=="Fail:" ) 
            # select("SITE","YEAR",msg_columns) # Select some columns         
# # write.csv(this_failing, file = paste0(path,"/Failings_file_data.csv"))
# Put warnings and failures in a file of warnings_failings  
Warnings_failings <- rbind(this_warning, this_failing)

head(Warnings_failings, 1)

# Add all fails and warnings
# Check if dataframe of warnings is empty, if not write to file
if(nrow(Warnings_failings)>0) {
   newdf<- data.frame(Warnings_failings)
   pdf("Fails_Warnings2.pdf", height=11, width=18.5)
   #Output the pdf file
   data.frame(grid.table(newdf))
   print(newdf)
}
# 
# # Data validation and conversion
# # 13.2 subset the instances to run in prediction by removing "this_failing", use anti-join i.e."Return all rows from x where there are no matching values in y, keeping just columns from x.
# # This is a filtering join"
# 

#head(raw.input.data, 15)
# final.predictors1 <- anti_join(raw.input.data, this_failing, by="SITE") # This works in R Studio, but not in ML AZURE
final.predictors1 <- raw.input.data[is.na(match(raw.input.data$SITE, this_failing$SITE)), ]
# DONT SORT, if you do , dont use the SORTED array for prediction. it duplicates the results ******

# # Generate data for classification 
#  # Final Data for classification e.g. Linear discriminant Analysis (LDA) classifier/predictor
# 
 final.predictors <- data.frame(
   SITE                     <-  final.predictors1$SITE,
   LATITUDE                 <-  final.predictors1$LATITUDE,
   LONGITUDE                <-  final.predictors1$LONGITUDE,
   LOG.ALTITUDE             <-  final.predictors1$vld_alt_src_log,
   LOG.DISTANCE.FROM.SOURCE <-  final.predictors1$vld_dist_src_log,
   LOG.WIDTH                <-  final.predictors1$mn_width_log,
   LOG.DEPTH                <-  final.predictors1$mn_depth_log,
   MEAN.SUBSTRATUM          <-  final.predictors1$vld_substr_log,
   DISCHARGE.CATEGORY       <-  final.predictors1$DISCHARGE,    #raw.input.data$disch_log,
   ALKALINITY               <-  final.predictors1$ALKALINITY,
   LOG.ALKALINITY           <-  final.predictors1$vld_alkal_log,
   LOG.SLOPE                <-  final.predictors1$vld_slope_log,
   MEAN.AIR.TEMP            <-  final.predictors1$TMEAN,
   AIR.TEMP.RANGE           <-  final.predictors1$TRANGE
 )
colnames(final.predictors) <- c("SITE","LATITUDE","LONGITUDE","LOG.ALTITUDE","LOG.DISTANCE.FROM.SOURCE","LOG.WIDTH","LOG.DEPTH","MEAN.SUBSTRATUM","DISCHARGE.CATEGORY","ALKALINITY","LOG.ALKALINITY", "LOG.SLOPE","MEAN.AIR.TEMP","AIR.TEMP.RANGE")

head(final.predictors,3)

# #   Prediction Settings
# #2. Find the DFScores of each row using one line of coefficients DFCoeff_gb685[1,-1] # removes the first column
# 
# # NRefg = number of reference sites in end group g , for GB = 685, for NI = 11
#NRefg <- 685
# NRefg_groups<- read.csv(paste0(path,"/EndGrp_AssessScores.csv"))
NRefg_all <- rowSums(NRefg_groups[,-1])
# 
# #DFScore_g <- DFCoef1 * Env1 + ... + DFCoefn * Envn ; remove "SITE" col=1 from final.predictors, and  remove col=1 from DFCoeff_gb685
DFScores <- as.data.frame(getDFScores(final.predictors, DFCoeff_gb685))

# 
# # Calculate the Mahanalobis disance of point x from site g for all referene sites
MahDist_g <- getMahDist(DFScores , DFMean_gd)
probs <- paste0("p",1:43)
MahDistNames <- probs
MahDistNames <- gsub("p","Mah",MahDistNames)
colnames(MahDist_g) <- MahDistNames

# # Calculate the minimum Mahanalobis disance of point x from site g
MahDist_min <- getMahDist_min(DFScores , DFMean_gd)
# #Calculate the probability distribution
PDist_g <- PDist (NRefg_all, MahDist_g)
# #Main dataframe needed:: Calculate probabilities of sites belonging to the endgroups, prob_g, l,as last column 44 contrains the total "PGdistTot
PDistTot <- as.data.frame(PDistTotal(PDist_g)) ## ALL probabilities p1..pn,  rowsums() add to 1, except when last row which it "total" is removed i.e. rowSums(PDistTot[,-ncol(PDistTot)])=1
# Rename the columns to probabilities p1,p2,...,p43

colnames(PDistTot) <- c(probs, "Total")
 
# #   final.predictors <- cbind(final.predictors, PDist_g[,-ncol(PDist_g)]) # This is the line we need
final.predictors_try1 <- cbind(final.predictors, PDistTot[,-ncol(PDistTot)]) # sum(final.predictors_try[1,-c(1:14)]) should give 1
#head(final.predictors_try1 , 1) # OK
# 
# #3.Use chisquare to find suitability codes. Start for Britain GB, # # Could use a file for these chisquare values
# # 1 = GB 21.02606 24.05393 26.21696 32.90923 
# # 2 = NI 18.30700 21.16080 23.20930 29.58830

chiSquare_vals <- data.frame(CQ1=c(21.02606, 18.30700), CQ2=c(24.05393,21.16080), 
                              CQ3=c(26.21696,23.20930), CQ4=c(32.90923,29.58830))
                              
suitCodes <- getSuitabilityCode(MahDist_min, chiSquare_vals)
#head(suitCodes,2)

# # add suitab ility codes to the final data, using cbind 
final.predictors_try2 <- cbind(final.predictors_try1, suitCodes)

# Find max class group belongs to by getting the column name: use
BelongsTo_endGrp <- colnames(final.predictors_try2[,15:57])[apply(data.frame(matrix(unlist(final.predictors_try2[,15:57]), nrow=nrow(final.predictors_try2[,15:57]), byrow=T),stringsAsFactors=FALSE), 1, which.max)]

#Relace p with EndGr
BelongsTo_endGrp <- gsub("p","EndGr",BelongsTo_endGrp)
final.predictors_try2 <- cbind(final.predictors_try2, BelongsTo_endGrp)
 
# Remove objects not in use
rm(PDist_g)
rm(PDist)
rm(PDistTotal)

#4 Prediction: WE1.5 Algorithms for prediction of expected values of any index based on probability of end group 
# membership and average values of the index amongst reference sites in each end group. Algorithms for taxonomic 
# prediction (including prediction of occurrence and expected abundances)
 
# Read TAXAPRAB taxa data

taxa.input.data <- readCSVFile("src/TAXA_AB_APR_PRAB.csv", skip = 1, stringsAsFactors = TRUE)
print("taxa.input.data = ")  
head(taxa.input.data,2)

# Isolate for particular model eg GB model
taxa_model      <<- model
taxa.input.data <- data.table(filter(taxa.input.data, Model == ifelse(taxa_model=="RIVPACS IV GB",1,2) ) ) # 1 is GB, 2 is NI

# Declare a dataframe where we append all sites
AllSites <- data.table()
 
# Use cpomplete cases removing null values
taxa.input.data <- data.table(taxa.input.data[complete.cases(taxa.input.data),])

# Use compiler function calls

#myFuncCmp     <- cmpfun(groupSitesFunction)
#EngGrpFuncCmp <- cmpfun(EndGrpProb_Replacement)

nsites <- nrow(final.predictors_try2)

taxa.input.data <- data.table::data.table(taxa.input.data)
taxa.input.data$EndGrp_Probs <- taxa.input.data$`End Group`

###### Loop for each SITE   
#foreach(i = 1:nsites, .combine=rbind)  %dopar% {
#for (i in 1:nsites){ 
for (i in chooseSite) {
      site1 <- taxa.input.data
      siteIndex <- ifelse(i<10,paste0("0",i),i)
      
      measuredColumns <- site1[,c("Average_Numerical_Abundance", "Average_Log10_Abundance", "Prob_Occurrence","Prob_Log1","Prob_Log2","Prob_Log3","Prob_Log4","Prob_Log5")]
      names_measuredColumns <- names(measuredColumns)      
      #site1$EndGrp_Probs <- EngGrpFuncCmp (site1$EndGrp_Probs,final.predictors_try2, i) # 
      site1$EndGrp_Probs <- EndGrpProb_Replacement (site1$EndGrp_Probs,final.predictors_try2, which(chooseSite==i))
      b1 <- site1 %>% dplyr::group_by(Season_Code, TL, Furse_Code) %>% 
         dplyr::mutate(count = n(), mlist_endGrps = list(EndGrp_Probs)) # `End Group`
 
      allUniqueSites <- unique(b1[,c(4,3,5,6,7,20,21)]) # gives 3,522 rows.or unique b1 Multiply by 12 sites gives us 3,532x12 = 42,384 groups in TOTAL.

      print(tail(allUniqueSites))
      for (k in 1:nrow(allUniqueSites)) { ### loop over these unique rows per SITE
            # sitex   <- myFuncCmp (allUniqueSites, k, siteIndex, b1)
            sitex    <- groupSitesFunction (allUniqueSites, k, siteIndex, b1)
            print(head(sitex))
            AllSites <- rbind(AllSites, sitex) 
    }# foreack k
     
}# for i 
 
#Remove the "End.Group" column
AllSites <- AllSites[,-c("End.Group")]
 
end_time <- proc.time()
 
print("Time = ")
print((proc.time()-start_time)/60)
maml.mapOutputPort("AllSites")
