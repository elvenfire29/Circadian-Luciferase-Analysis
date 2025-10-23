#Version 2025.October.23.
#Luminescence (Clock) Assay. SIMPLIFIED ARSER. With ANOVA and Tukey HSD.


#PLEASE READ
##############
#Need to download packages:
# MetaCycle  ggplot2  dplyr  magrittr  stringr  filesstrings  circular  AICcmodavg  broom

#Pre-Input Instructions:
#1. Download the above packages into RStudio if you have not previously done so.
#2. Prep input by deleting all cells other than the luminescence data with the Time row and accompanying column.
#3. Time should be in whole number hours as a continuous time series.
#4. Save as CSV file.


########## User Input ###########

#User Input I: Steps 1-4 below must be changed to tailor the analysis to a specific data set on a local computer.

#Step 1: Set the working directory
setwd("C:/Users/pippi/Downloads")

#Step 2: Select your input file (should be a csv file).
plate_no <- 'NO7.csv'

#Step 3: Name the treatments and/or lines used.
#Per specific experiment, please change sample names accordingly from A to H, OR from 1 to 12 as on the plate.
#The number and order of labels matters!! This is because the code counts the number (needs 8 or 12) and uses the order.
#Same treatment labels will combine into one group in the analysis.
#Please avoid using a backslash in any treatment names!! Please avoid "NA" as a label when using ANOVA!
treatment <- c(
'Sample1',
'Sample2',
'Sample3',
'Sample4',
'Sample5',
'Sample6',
'Sample7',
'Sample8'
  )

#Step 4: Indicate the relative start time for the assay.
#For example, if subjective dawn was 7 am and the luciferase assay was started at 9 am, the user would input 2 as the relative start time.
#This number must be an integer and cannot be negative!
RelativeTime <- 0


###############
#User Input II: Options 1-7 can be changed per user's specific needs.

#Option 1: Include graphs?
#If you want basic graphs of the luminescent curves and the graphs of average statistics, type TRUE. If not, type FALSE.
GRAPHS <- TRUE

################

#Option 2: Run a One-Way ANOVA test with Tukey HSD post-hoc?
#Code will run the test if TRUE is selected. If FALSE is selected, it will not.
ANOVA <- TRUE

#Option 2 addition 1: What is your control?
#You can specify the output control for the ANOVA test.
#Leave empty to use each treatment as a control in a pairwise comparison.
A_CONTROL <- ''

#Option 2 addition 2: Number ANOVA output files?
#If we are doing all treatments as controls for an ANOVA, do you want the file names numbered for convenience?
#Type TRUE or FALSE for yes or no, respectively.
A_NumberLabel <- TRUE

################

#Option 3: Run a t-test?
#If yes, type TRUE; if no, type FALSE.
#If an error message says data is too similar, write FALSE and use ANOVA or do t-test by hand.
t_test <- TRUE

#Option 3 addition 1: Should the t-test be pairwise?
#If yes, type TRUE; if no, type FALSE
pairwise<- TRUE

#Option 3 addition 2: If the test is pairwise, are the data paired?
#If yes, type TRUE; if no, type FALSE
paired<- FALSE

#Option 3 addition 3: What treatment is the control?
#Leave the list empty to use each treatment as control -OR- write the treatment name of the single control.
T_CONTROL <- ''

#Option 3 addition 4: Number t-test output files?
#If we are doing all treatments as controls for a t-test, do you want the file names numbered for convenience?
#Type TRUE or FALSE for yes or no, respectively.
T_NumberLabel <- TRUE

################

#Option 4: Round your time points?
#If yes, type TRUE. If no, type FALSE.
#ARSER analysis can only be run if the time series is consistent. If it is not, try rounding to the nearest hour.
#Not advisable if times are off by more than a couple of minutes from the hour interval!
roundTP <- TRUE

################

#Option 5: Does your data list wells by A1,A2,A3... OR by A1,B1,C1...?
#If your data lists wells by A1,A2,A3.... then type TRUE. Otherwise, type FALSE.
orientation <- TRUE

#################################


#END OF USER INPUT
########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################


########################
#Call libraries#
#######libraries#################
library(MetaCycle)

library(ggplot2)
library(dplyr)

library(stringr)
library(filesstrings)
library(magrittr)
library(circular)

library(broom)
library(AICcmodavg)

########################
# FUNCTION DEFINITIONS #
##########functions##############

GetTimepoints <- function(name) {
  #get actual timepoints
  timepoints <- c()

  #testing timepoints for minutes by splitting at h and looking at remaining lengths
  for (n in name) {
    header <- n[1]
    header <- gsub(" ", "", header, fixed = TRUE)
    timet <- 0

    #if it starts with an "X", remove it
    if (grepl("X",header,fixed=TRUE)) {
      header <- unlist(strsplit(header, "X"))
      header<-header[2]
    }

    #if it ends in ".min", remove that whole chunk, and add the minutes
    if (grepl("min",header,fixed=TRUE)) {
      #strip the ".min" part
      header <- unlist(strsplit(header, ".min"))
      header<-header[1]
      #split on the ".h."
      if (grepl("h",header,fixed=TRUE)){
        sepTime <- unlist(strsplit(header, ".h."))
        header <- sepTime[1]
        timet <- timet + (as.numeric(sepTime[2])/60)
      }
    }

    #if it ends in ".h", remove that
    if (grepl("h",header,fixed=TRUE)) {
      header <- unlist(strsplit(header, ".h"))
      header <- header[1]
    }

    #we should be left with just the hours, so let's add them
    timet <- timet + (as.numeric(header))

    #append it to the timepoints array
    timepoints <- cbind(timepoints, timet)
  }

  return(timepoints)
}

skip<-roundTP

GetZTtim <- function(timepoints, ZT, skip) {
  #to reflect true phase, we shift time points to relative time (ZT)
  ZTtim <- c()
  for (i in timepoints) {
    ZTtim = cbind(ZTtim, i + ZT)
  }

  #if skip==TRUE, round to nearest hour (use with discretion)
  if (skip) {
    ZTtim <- round(ZTtim)
  }

  return(ZTtim)
}

# Gets the PerPhaAmp for a given plate and name.
GetPerPhaAmp <- function(plate, name, ZT, method) {
  #get actual timepoints
  timepoints <- GetTimepoints(name)

  #to reflect true phase, we shift time points to relative time (ZT)
  ZTtim <- GetZTtim(timepoints, ZT, skip)

  #function that gives us period, amplitude, and phase
  #timeUnit is auto hour, but can be "day", "minute", "hour", or "second"
  cyc_1 <- meta2d(
    infile = as.character(plate),
    filestyle = 'csv',
    timepoints = ZTtim,
    minper = 18,
    maxper = 30,
    cycMethod = method,
    outIntegration = "onlyIntegration",
    adjustPhase = "predictedPer",
    outputFile = F,
    outRawData = TRUE,
    weightedPerPha = T
  )
  PerPhaAmp = cyc_1$meta

  if (shift) {
    #if adjphase > YY, subtract by XX
    adjust_phase <- function(e) {
      if (e > YY){
        return(as.numeric(e) - XX)
        }
      else{
        return(e)
      }
    }

    if ('ARS_adjphase' %in% colnames(PerPhaAmp)){
      PerPhaAmp[, 'ARS_adjphase'] <-
        sapply(PerPhaAmp[, 'ARS_adjphase'], adjust_phase)
    }
  }

  return(PerPhaAmp)
}


GetWellSeq <- function(i, treatment, orientation) {
  #well_seq<-if(length(treatment) == 8) seq((i-1)*12+1,12*i,1) else seq(i,84+i,12)
  #change collection of samples based on the plate readout
  if (orientation == TRUE) {
    seq_12 <- seq(i, 84 + i, 12)
    seq_8 <- seq((i - 1) * 12 + 1, 12 * i, 1)
  }
  #change collection of samples based on the plate readout
  if (orientation == FALSE) {
    seq_8 <- seq(i, 88 + i, 8)
    seq_12 <- seq((i - 1) * 8 + 1, 8 * i, 1)
  }
  if (length(treatment) == 8) {
    return(seq_8)
  } else {
    return(seq_12)
  }
}

################################code starts/functions end#################################################

#reassign variable name
NumberLabel<-T_NumberLabel
pvalues<-t_test

#get working directory
wd <- c(getwd())

#create new folder for output
plate_num <- strsplit(plate_no, split = ".csv")
folder <- paste(plate_num[1], 'output')
if (!dir.exists(folder)) {
  dir.create(folder)
}

#Get location of new folder
wd <- getwd()
newloco <- paste(wd, folder, sep = "/")

#set location of new sub-folders
folder2 <- paste("#6 __", plate_num[1], 't_test')
folder3 <- paste("#5 __", plate_num[1], 'ANOVA')
folder4 <- paste("#4 __", plate_num[1], 'Individual Wells')

t_testFolder <- paste(newloco,folder2,sep='/')
ANOVAFolder <- paste(newloco,folder3,sep='/')
newloco2 <- paste(newloco,folder4,sep='/')

#create more folders within the new output folder for more organization
if(pvalues){
  if (!dir.exists(t_testFolder)) {
    dir.create(t_testFolder)
  }
}
if(ANOVA){
  if (!dir.exists(ANOVAFolder)) {
    dir.create(ANOVAFolder)
  }
}

#save treatment names as they originally are
OGTreat<-treatment
CONTROL<-T_CONTROL
OGCONTROL<-CONTROL
OGA_CONTROL<-A_CONTROL
#Take out replicated treatment names
OGTreat_unique<-unique(OGTreat)

#make treatment labels ok for file names by taking out certain symbols /  :  -  .  ~
# the backslash is never ok to use here
for (i in 1:length(treatment)){
  treatment[i]<-gsub('/', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub(':', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub('-', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub('.', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub('~', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
}
for (i in 1:length(CONTROL)){
  CONTROL[i]<-gsub('/', '_', CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  CONTROL[i]<-gsub(':', '_', CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  CONTROL[i]<-gsub('-', '_', CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  CONTROL[i]<-gsub('.', '_', CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  CONTROL[i]<-gsub('~', '_', CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
}
for (i in 1:length(A_CONTROL)){
  A_CONTROL[i]<-gsub('/', '_', A_CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  A_CONTROL[i]<-gsub(':', '_', A_CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  A_CONTROL[i]<-gsub('-', '_', A_CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  A_CONTROL[i]<-gsub('.', '_', A_CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  A_CONTROL[i]<-gsub('~', '_', A_CONTROL[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
}

#gets labels (time points) from infile
time_points = read.csv(as.character(plate_no))
time_points = time_points[, 2:ncol(time_points)]
names = colnames(time_points)
names = as.matrix(names)

#labels for the calculations we are interested in
PPA_range <-
  c(
    'ARS_period',
    'ARS_adjphase',
    'ARS_amplitude'
  )

#not 8 or 12 treatments. doesn't know how to handle that
if ((length(treatment) != 8) & (length(treatment) != 12)) {
  stop("Please input labels for 8 or 12 treatments.")
}

#program decides if there are 12 samples
if (length(treatment) == 12) {
  print("There are 12 treatments.")
} else {
  print("There are 8 treatments.")
}

#base repetitions per treatment change based on number of treatments
#this does not account for the repeated same-name treatments that are later combined
rep <- if (length(treatment) == 8){12}else{8}

plate <- plate_no

#find per, phase, amp
ZT<-RelativeTime
method <- c('ARS')
PerPhaAmp <- GetPerPhaAmp(plate, names, ZT, method)

#separate well values into per, pha, amp ARS for each treatment
#this chunk is specific to 12 treatments with 8 replicates per
preARS <- match(nomatch = F, names(PerPhaAmp), x = "ARS_period")

#ARS is not computed
if (preARS <= 0) {
  stop(
    "Something went wrong with calculating the period, phase, and amplitude. Please review your input."
  )
}

#only ARS is computed
if (preARS > 0){
  print("ARS analysis is available.")
  aver_source <- c(seq(1:rep))
  PPA_range <- c('ARS_period', 'ARS_adjphase', 'ARS_amplitude')
  for (i in 1:length(treatment)) {
    well_seq <- GetWellSeq(i, treatment, orientation)
    avery <- vector()
    for (x in PPA_range) {
      aver <- unlist(PerPhaAmp[well_seq, x])
      if (x == PPA_range[1]) {
        avery <- c(aver)
      }
      if (x != PPA_range[1]) {
        avery <- cbind(avery, aver)
      }
    }
    aver_source <- cbind(aver_source, avery)
  }
}

###write dividers for eventual outfiles
Average <- c()
SEM <- c()
P_value <- c()

for (i in 1:length(PPA_range)) {
  du <- c("..")
  Average <- cbind(Average,du)
  SEM <- cbind(SEM,du)
  P_value <- cbind(P_value,du)
}

doh <- vector()
for (i in 1:rep) {
  duh <- c("..")
  doh <- cbind(doh, duh)
}
indv_wells <- as.vector(doh)

GetDataForAvererColumn <- function(averer, x, col_count) {
  temp <- vector()

  name <- treatment[floor((x - 1) / col_count) + 1]
  offset <- mod(x - 1, col_count)

  #find everyone else in treatment with the same name
  for (other_treatment_index in 1:length(treatment)) {
    other_treatment_name <- treatment[other_treatment_index]
    if (name == other_treatment_name) {
      other_x <- ((other_treatment_index - 1) * col_count) + offset + 1
      temp <- c(temp, averer[, other_x])
    }
  }

  return(na.omit(temp))
}

##############period, phase, and amplitude############
#get data
averer <- aver_source[,-1]

#find average and sem of each treatment
#find sem
averer <- aver_source[,-1]
combined_average <- "mean"
combined_sd <- "sem"
ind_wells <- indv_wells

#do this for all col in averer
count2 = 0
for (x in 1:length(averer[1, ])) {
  count2 = count2 + 1

  #get rid of NA values
  temp <- vector()
  temp2 <- vector()

  temp2 <- GetDataForAvererColumn(averer, x, length(PPA_range))

  #do calculations
  sam_treatment <- (sd(temp2) / sqrt(length(temp2)))
  average_treatment <- mean(temp2)

  temp2 <- na.omit(averer[, x])

  #save individual values
  wells <- temp2

  #collect the calculations
  combined_average <- cbind(combined_average, average_treatment)
  combined_sd <- cbind(combined_sd, sam_treatment)

  #each column is a repetition and each section is a line and each row is a value (amp,per,phase)
  if (count2 <= length(PPA_range)) {
    ind_wells <- rbind(ind_wells, wells)
  }
  if (count2 > length(PPA_range)) {
    ind_wells <- rbind(ind_wells, indv_wells, wells)
    count2 <- 1
  }
}

#prep to reorganize averages so that every seventh (not counting the first one/label) is stacked
pretty_avg <- PPA_range
pretty_sd <- PPA_range

#take out label col
com_sd <- as.matrix(combined_sd[,-1])
com_av <- as.matrix(combined_average[,-1])
com_sd <- as.numeric(com_sd)
com_av <- as.numeric(com_av)

#actually reorganize
indv_lab <- vector()
for (i in 1:(length(treatment))) {
  raverage <-
    unlist(com_av[seq((i - 1) * length(PPA_range) + 1, length(PPA_range) * i, 1)])
  pretty_avg <- rbind(pretty_avg, raverage)

  rsem <-
    c(com_sd[seq((i - 1) * length(PPA_range) + 1, length(PPA_range) * i, 1)])
  pretty_sd <- rbind(pretty_sd, rsem)

  #organize labels for indiv wells data
  if(length(treatment)==12){
    indv_lab <- rbind(indv_lab, i)
  }
  if(length(treatment)==8){
    indv_lab <- rbind(indv_lab, LETTERS[i])
  }

  for (x in PPA_range) {
    indv_lab <- rbind(indv_lab, x)
  }
}

#redo col label
colnames(pretty_avg) <- PPA_range
colnames(pretty_sd) <- PPA_range

pretty_avg <- pretty_avg[-1, ]
pretty_sd <- pretty_sd[-1, ]

#colname for indiv wells depends on how many replicates there are
if (length(treatment) == 12) {
  wellnum <- c()
  for (i in 1:8) {
    welln <- paste('well', LETTERS[i])
    wellnum <- cbind(wellnum, welln)
  }
  colnames(ind_wells) <- wellnum
}

if (length(treatment) == 8) {
  wellnum <- c()
  for (i in 1:12) {
    welln <- paste('well', i)
    wellnum <- cbind(wellnum, welln)
  }
  colnames(ind_wells) <- wellnum
}

#label rows
row.names(ind_wells) <- indv_lab
row.names(pretty_avg) <- OGTreat
row.names(pretty_sd) <- OGTreat

pretty_avg <- pretty_avg[!duplicated(rownames(pretty_avg)), ]
pretty_sd <- pretty_sd[!duplicated(rownames(pretty_sd)), ]

#Reorganize for ANOVA: Plant_Line, Period, Phase, Amplitude columns
#Get the column of treatment names for individual well values
DataForANOVA<-c()
for (x in 1:length(OGTreat)){
  temp<-c(OGTreat[x])
  temp2<-c()
  for (i in 1:rep){
    temp2<-rbind(temp2,temp)
  }
  if(i==1 & x==1){
    DataForANOVA<-c(temp)
  }
  else{
    DataForANOVA<-rbind(DataForANOVA,temp2)
  }
}

#Get columns for period, phase, amplitude ready to add to treatment labels column
for (x in 1:length(treatment)){
  if (x==1){
    tempData2<-averer[,seq(x,x+2,1)]
  }
  else{
    tempData<-averer[,seq(((x-1)*3)+1,(((x-1)*3)+1)+2,1)]
    tempData2<-rbind(tempData2,tempData)
  }
}

#combine the data and labels
DataForANOVA<-cbind(DataForANOVA,tempData2)

#find number of unique treatments
Treatments0<-unique(DataForANOVA[,1])
Treatments<-DataForANOVA[,1]

#move same treatment names so that they are near one another (doesn't affect ANOVA code, just a preference)
if(length(Treatments0) != length(treatment)){
  potkettle<-c()
  blacklist<-c()
  for(i in 1:length(Treatments)){
    #find number of replicates per treatment
    Repeats<-which(Treatments==Treatments[i])
    #get all reps for a single treatment
    if(i==1){
      count=0
      for(k in Repeats){
        if(count==0){
          potkettle<-c(DataForANOVA[k,])
          count=count+1
        }
        else{
          potkettle<-rbind(potkettle,DataForANOVA[k,])
        }
      }
    }
    else{
      #make sure we aren't repeating treatment replicates
      if((Treatments[i] %in% blacklist)==FALSE){
        for(p in Repeats){
        potkettle<-rbind(potkettle,DataForANOVA[p,])
        }
      }
    }
    blacklist<-c(blacklist,Treatments[i])
  }
  DataForANOVA<-potkettle
}

#alphabetize data
#Alpha<-DataForANOVA[order(DataForANOVA[,1]),]

#Label the column names
colnames(DataForANOVA)<-c("Plant_Line","Period","Phase","Amplitude")
row.names(DataForANOVA)<-NULL

#write file for running ANOVA code for easier re-running in case of changes or well deletions later
plate_num <- strsplit(plate_no, split = ".csv")
titleT<-paste("DataForANOVA",".csv",sep="")
title2<-paste(plate_num,titleT)
write.csv(x=DataForANOVA,title2,quote=FALSE,row.names = FALSE)

if(ANOVA==FALSE){
  file.move(title2, newloco, overwrite = TRUE)
}

#run if selected is TRUE
if (ANOVA==TRUE){

  #read in data from input file
  my_data<-read.csv(title2,header=TRUE,colClasses=c("factor","numeric","numeric","numeric"))

  #run ANOVA code
  #check data is read in properly
  print(summary(my_data))

  #run anova for the period
  #run the anova
  period <- aov(Period~Plant_Line,data = my_data)

  #print out the summary of values from the anova
  print(summary(period))

  #ad-hoc test
  per<-TukeyHSD(period)

  #run again for the phase
  #run the anova
  phase <- aov(Phase~Plant_Line,data = my_data)

  #print out the summary of values from the anova
  print(summary(phase))

  #ad-hoc test
  pha<-TukeyHSD(phase)

  #run again for the amplitude
  #run the anova
  amplitude <- aov(Amplitude~Plant_Line,data = my_data)

  #print out the summary of values from the anova
  print(summary(amplitude))

  #ad-hoc test
  amp<-TukeyHSD(amplitude)

  #put results into outfile to read-in below
  otheroutA <- c('ANOVA_Tukey All Results.csv')
  stats<-c()
  stats<-c(per,"Period",pha,"Phase",amp,"Amplitude")
  write.csv(x = stats, otheroutA)

  #read in the all results ANOVA file
  ANOVA_PVALUES<-read.csv(otheroutA)
  #edit so that it is easier to read (i.e. take out unnecessary columns and leave p-value)
  for (i in 1:length(ANOVA_PVALUES)){
    if(i==1){
      ANOVA_PVALUES<-ANOVA_PVALUES[,-seq(2,4,1)]
    }
    if(i!=1){
      ANOVA_PVALUES<-ANOVA_PVALUES[,-seq(i+1,i+4,1)]
    }
  }
  #Change column names to be more descriptive
  ANOVA_COL<-c('Treatments Compared','Period p-value','Phase p-value','Amplitude p-value')
  colnames(ANOVA_PVALUES)<-ANOVA_COL

  #Look for ANOVAs with the control treatment
  comparisons<-ANOVA_PVALUES[,1]
  #make sure that if we can tell wt from wt with chemical, etc
  #add spaces to ends of each label to attempt to do so
  for (n in 1:length(comparisons)){
    comparisons[n] <-paste(".",comparisons[n],".",sep="")
  }

  #if there is a single control
  if (length(A_CONTROL) > 1 ||
      (length(A_CONTROL) == 1 && A_CONTROL[1] != "")) {

    cup<-c()
    for(n in 1:length(comparisons)){

      #Some treatments may include the name of other treatments, but with added stuff
      chonks1 <-comparisons[n]
      #take out control and separating - so that only the treatment tested against control remains
      test1<-paste("-",OGA_CONTROL,".",sep="")
      test2<-paste(".",OGA_CONTROL,"-",sep="")

      #if TRUE, then add the row to the collection cup for the file
      temp<-c()
      if(str_detect(chonks1, test1)==TRUE){
        #make sure ID is 100% correct
        ID<-sub(test1,"",chonks1)
          #take out the "." that is left in ID at beginning or end
          #count occurrences of "." so that we don't ruin annotation of treatment
          COUNT<-sum(charToRaw(ID) == charToRaw('.'))
          if(COUNT==1){
            #remove "."
            tempID<-str_remove(ID,"[.]")
          }
          #remove only the first or last "."
          if(COUNT>1){
            #find if "." is first or last
            if(unlist(gregexpr('.', ID))[1]==1){
              #get rid of first "." only
              tempID<-sub('.', '', ID)
            }
            last<-tail(unlist(gregexpr('.', ID)), n=1)
            temptemptemp<-unlist(ID)
            if(nchar(temptemptemp)==last){
              #get rid of last "." only by knocking out the last character, which we assume is "."
              tempID<-str_sub(ID,1,-2)
            }
          }
          #compare to other names to see if same
          for (u in 1:length(OGTreat_unique)){
            if(tempID==OGTreat_unique[u]){
              if(tempID!=OGA_CONTROL){
                temp<-ANOVA_PVALUES[n,]
                if(length(cup)==0){
                  cup<-temp
                }
                else{
                  cup<-rbind(cup,temp)
                }
              }
          }
        }
      }
      if(str_detect(chonks1, test2)==TRUE){
        #make sure ID is 100% correct
        ID<-sub(test2,"",chonks1)
          #take out the "." that is left in ID at beginning or end
          #count occurrences of "." so that we don't ruin annotation of treatment
          COUNT<-sum(charToRaw(ID) == charToRaw('.'))
          if(COUNT==1){
            #remove "."
            tempID<-str_remove(ID,"[.]")
          }
          #remove only the first or last "."
          if(COUNT>1){
            #find if "." is first or last
            if(unlist(gregexpr('.', ID))[1]==1){
              #get rid of first "." only
              tempID<-sub('.', '', ID)
            }
            last<-tail(unlist(gregexpr('.', ID)), n=1)
            temptemptemp<-unlist(ID)
            if(nchar(temptemptemp)==last){
              #get rid of last "." only by knocking out the last character, which we assume is "."
              tempID<-str_sub(ID,1,-2)
            }
          }
          #compare to other names to see if same
          for (u in 1:length(OGTreat_unique)){
            if(tempID==OGTreat_unique[u]){
              if(tempID!=OGA_CONTROL){
                temp<-ANOVA_PVALUES[n,]
                if(length(cup)==0){
                  cup<-temp
                }
                else{
                  cup<-rbind(cup,temp)
                }
              }
          }
        }
      }
    }

    #alter data to fit with PER,PHA,AMP data
    #make the comparisons easier to read by listing the one that's not the control only
    rows<-cup[,1]
    for (n in 1:length(rows)){
      rows[n] <-paste(".",rows[n],".",sep="")
    }
    spy<-c()
    for(i in 1:length(rows)){
      chonks2 <-rows[i]
      #take out control and separating - so that only the treatment tested against control remains
      test1<-paste("-",OGA_CONTROL,".",sep="")
      test2<-paste(".",OGA_CONTROL,"-",sep="")

      if(str_detect(chonks2, test1)){
        chunky<-unlist(strsplit(chonks2,split=test1))
        chunky<-chunky[1]

        #collect the new row labels (treatment names tested against control)
        if(i==1){
          if(length(chunky)!=0){
          spy<-chunky
          }
          else{
            spy<-chonks2
          }
        }
        else{
          if(length(chunky)!=0){
            spy<-rbind(spy,chunky)
          }
          else{
            spy<-rbind(spy,chonks2)
          }
        }
      }
      #if it doesn't find one orientation, try the other
      else{
        if(str_detect(chonks2, test2)){
        chunky<-unlist(strsplit(chonks2,split=test2))
        chunky<-chunky[2]

        #collect the new row labels (treatment names tested against control)
        if(i==1){
          if(length(chunky)!=0){
          spy<-chunky
          }
          else{
            spy<-chonks2
          }
        }
        else{
          if(length(chunky)!=0){
          spy<-rbind(spy,chunky)
          }
          else{
            spy<-rbind(spy,chonks2)
          }
          }
        }
      }
    }
    cup<-cup[,-1]
    #see if marker is still present and remove
    spy<-na.omit(spy)
    if(length(spy)!=0){
    for(y in 1:length(spy)){
      if(str_detect(spy[y],".")){
        temp12<-spy[y]
        temp14<-unlist(strsplit(temp12,".",fixed=TRUE))
        for(w in OGTreat_unique){
        if(temp14[1]==w){
          spy[y]<-temp14[1]
        }
         if(length(temp14)>1){
           if(temp14[-1]==w){
            spy[y]<-temp14[-1]
           }
         }
        }
      }
    }
    spy<-unique(spy)
    #combined labels and data
    cup<-cbind(spy,cup)

    #need new column names (hoho is now too short)
    hoho<-PPA_range
    hoho2<-append("Treatment",hoho)
    hoho<-hoho2[-1]

    colnames(cup)<-hoho2
    rownames(cup)<-NULL

    #set fake row names of spacers by adding col
    if (length(Average)==3){
      Average<-cbind("Average",Average)
      colnames(Average)<-hoho2
    }
    if (length(SEM)==3){
      SEM<-cbind("SEM",SEM)
      colnames(SEM)<-hoho2
    }
    if (length(P_value)==3){
      P_value<-cbind("P_value",P_value)
      colnames(P_value)<-hoho2
    }

    #make a list of 1's for the control vs self
    one<-c()
    for (o in 1:length(hoho)){
      one<-cbind(one,"1")
    }
    one<-cbind(as.character(OGA_CONTROL),one)
    colnames(one)<-hoho2

    #replace rownames as columns to avoid confusing numbering of repeat rownames
    avgRN<-rownames(pretty_avg)
    rownames(pretty_avg)<-NULL
    pretty_avg<-cbind(avgRN,pretty_avg)
    colnames(pretty_avg)<-hoho2

    semRN<-rownames(pretty_sd)
    rownames(pretty_sd)<-NULL
    pretty_sd<-cbind(semRN,pretty_sd)
    colnames(pretty_sd)<-hoho2

    #make file
    file5<-paste(A_CONTROL,"ANOVA p-values.csv")
    finalA<-rbind(Average,pretty_avg,SEM,pretty_sd,P_value,one,cup)
    rownames(finalA)<-NULL
    write.csv(x = finalA, file5)
    file.move(file5, ANOVAFolder, overwrite = TRUE)
    }
  }
  ###############ALL CONTROLS FOR ANOVA###########################
  #if there are all the controls
  else{
    for(k in 1:length(OGTreat_unique)){
      cup<-c()
      #take out control and separating - so that only the treatment tested against control remains
      test1<-paste("-",OGTreat_unique[k],".",sep="")
      test2<-paste(".",OGTreat_unique[k],"-",sep="")

      #search for each treatment separately
      for(n in 1:length(comparisons)){
        #Some treatments may include the name of other treatments, but with added stuff
        chonks1 <-comparisons[n]
        temp<-c()

        #if TRUE, then add the row to the collection cup for the file
        if(str_detect(chonks1, test1)==TRUE){
          #make sure ID is 100% correct
          ID<-sub(test1,"",chonks1)
            #take out the "." that is left in ID at beginning or end
            #count occurrences of "." so that we don't ruin annotation of treatment
            COUNT<-sum(charToRaw(ID) == charToRaw('.'))
            if(COUNT==1){
            #remove "."
              tempID<-str_remove(ID,"[.]")
            }
            #remove only the first or last "."
            if(COUNT>1){
              #find if "." is first or last
              if(unlist(gregexpr('.', ID))[1]==1){
                #get rid of first "." only
                tempID<-sub('.', '', ID)
              }
              last<-tail(unlist(gregexpr('.', ID)), n=1)
              temptemptemp<-unlist(ID)
              if(nchar(temptemptemp)==last){
                #get rid of last "." only by knocking out the last character, which we assume is "."
                tempID<-str_sub(ID,1,-2)
              }
            }
            #compare to other names to see if same
            for (u in 1:length(OGTreat_unique)){
              if(tempID==OGTreat_unique[u]){
                if(tempID!=OGTreat_unique[k]){
                  temp<-ANOVA_PVALUES[n,]
                  if(length(cup)==0){
                    cup<-temp
                  }
                  else{
                    cup<-rbind(cup,temp)
                  }
                }
              }
            }
        }
        #if the other way is TRUE, add it in!
        if(str_detect(chonks1, test2)==TRUE){
          #make sure ID is 100% correct
          ID<-sub(test2,"",chonks1)
            #take out the "." that is left in ID at beginning or end
            #count occurrences of "." so that we don't ruin annotation of treatment
            COUNT<-sum(charToRaw(ID) == charToRaw('.'))
            if(COUNT==1){
              #remove "."
              tempID<-str_remove(ID,"[.]")
            }
            #remove only the first or last "."
            if(COUNT>1){
              #find if "." is first or last
              if(unlist(gregexpr('.', ID))[1]==1){
                #get rid of first "." only
                tempID<-sub('.', '', ID)
              }
              last<-tail(unlist(gregexpr('.', ID)), n=1)
              temptemptemp<-unlist(ID)
              if(nchar(temptemptemp)==last){
                #get rid of last "." only by knocking out the last character, which we assume is "."
                tempID<-str_sub(ID,1,-2)
              }
            }
            #compare to other names to see if same
            for (u in 1:length(OGTreat_unique)){
              if(tempID==OGTreat_unique[u]){
                  if(tempID!=OGTreat_unique[k]){
                    temp<-ANOVA_PVALUES[n,]
                    if(length(cup)==0){
                      cup<-temp
                    }
                    else{
                      cup<-rbind(cup,temp)
                    }
                  }
              }
            }
          }
        }

      #alter data to fit with PER,PHA,AMP data
      #make the comparisons easier to read by only listing the one that's not the control
      rows<-cup[,1]
      for (n in 1:length(rows)){
        rows[n] <-paste(".",rows[n],".",sep="")
      }
      spy<-c()
      for(i in 1:length(rows)){
        chonks2 <-rows[i]
        #take out control and separating - so that only the treatment tested against control remains
        test1<-paste("-",OGTreat_unique[k],".",sep="")
        test2<-paste(".",OGTreat_unique[k],"-",sep="")

        if(str_detect(chonks2, test1)){
          chunky<-unlist(strsplit(chonks2,split=test1))
          chunky<-chunky[1]

          #collect the new row labels (treatment names tested against control)
          if(i==1){
            if(length(chunky)!=0){
            spy<-chunky
            }
            else{
              spy<-chonks2
            }
          }
          else{
            if(length(chunky)!=0){
            spy<-rbind(spy,chunky)
            }
            else{
              spy<-rbind(spy,chonks2)
            }
          }
        }
        #if it doesn't find one orientation, try the other
        else {
          if(str_detect(chonks2, test2)){
          chunky<-unlist(strsplit(chonks2,split=test2))
          chunky<-chunky[2]

          #collect the new row labels (treatment names tested against control)
          if(i==1){
            if(length(chunky)!=0){
            spy<-chunky
            }
            else{
              spy<-chonks2
            }
          }
          else{
            if(length(chunky)!=0){
            spy<-rbind(spy,chunky)
            }
            else{
              spy<-rbind(spy,chonks2)
            }
          }
          }
          }
      }
      cup<-cup[,-1]
      #see if marker is still present and remove
      spy<-na.omit(spy)
      if(length(spy)!=0){
      for(y in 1:length(spy)){
        if(str_detect(spy[y],".")){
          temp12<-spy[y]
          temp14<-unlist(strsplit(temp12,".",fixed=TRUE))
          for(w in OGTreat_unique){
            if(temp14[1]==w){
              spy[y]<-temp14[1]
            }
            if(length(temp14)>1){
              if(temp14[-1]==w){
                spy[y]<-temp14[-1]
              }
            }
          }
        }
      }
      spy<-unique(spy)
      #combine labels and data again
      cup<-cbind(spy,cup)

      #make file
      treatmentsU<-unique(treatment)
      if(A_NumberLabel){
        numnum<-paste("#",k,sep="")
        file6<-paste(numnum,"__",treatmentsU[k],"ANOVA p_values.csv")
      }
      else{
        file6<-paste(treatmentsU[k],"ANOVA p_values.csv")
      }

      #need new column names (hoho is now too short)
      hoho<-PPA_range
      hoho2<-append("Treatment",hoho)
      hoho<-hoho2[-1]

      colnames(cup)<-hoho2
      rownames(cup)<-NULL

      #set fake row names of spacers by adding col
      if (length(Average)==3){
        Average<-cbind("Average",Average)
        colnames(Average)<-hoho2
      }
      if (length(SEM)==3){
        SEM<-cbind("SEM",SEM)
        colnames(SEM)<-hoho2
      }
      if (length(P_value)==3){
        P_value<-cbind("P_value",P_value)
        colnames(P_value)<-hoho2
      }

      #make a list of 1's for the control vs self
      one<-c()
      for (o in 1:length(hoho)){
        one<-cbind(one,"1")
      }
      one<-cbind(as.character(treatment[k]),one)
      colnames(one)<-hoho2

      #replace rownames as columns to avoid confusing numbering of repeat rownames
      avgRN<-rownames(pretty_avg)
      rownames(pretty_avg)<-NULL
      pretty_avg<-cbind(avgRN,pretty_avg)
      colnames(pretty_avg)<-hoho2

      semRN<-rownames(pretty_sd)
      rownames(pretty_sd)<-NULL
      pretty_sd<-cbind(semRN,pretty_sd)
      colnames(pretty_sd)<-hoho2

      #make file
      finalA<-rbind(Average,pretty_avg,SEM,pretty_sd,P_value,one,cup)
      #row names set-up
      rownum<-(nrow(finalA)-3)
      treatnum<-(rownum/3)
      #set up for Average portion
      finalLAB<-c(".")
      for(i in 1:treatnum){
        numlet<-paste(i,"Avg")
        finalLAB<-rbind(finalLAB,numlet)
      }
      #set up for SEM portion
      finalLAB<-rbind(finalLAB,"..")
      for(i in 1:treatnum){
        numlet<-paste(i,"Sem")
        finalLAB<-rbind(finalLAB,numlet)
      }
      #set up for p-value portion
      finalLAB<-rbind(finalLAB,"...")
      for(i in 1:treatnum){
        numlet<-paste(i,"Pval")
        finalLAB<-rbind(finalLAB,numlet)
      }
      rownames(finalA)<-finalLAB
      #write file
      write.csv(x = finalA, file6)
      file.move(file6, ANOVAFolder, overwrite = TRUE)
    }
    }
  }

  #write out new file
  file4<-paste(plate_num,"All ANOVA Results.csv")
  write.csv(x = ANOVA_PVALUES, file4)
}

#####################rearranging individual PerPhaAmp again###############
#rearrange data for easier use in Prism GraphPad or similar
pairData<-t(my_data)

#get just that type of data
labels2<-pairData[1,]
periodPair<-t(pairData[2,])
phasePair<-t(pairData[3,])
ampPair<-t(pairData[4,])

#stack lines
periodStacked<-periodPair[seq(1,rep,1)]
phaseStacked<-phasePair[seq(1,rep,1)]
ampStacked<-ampPair[seq(1,rep,1)]

for (i in 2:length(treatment)){
  periodStacked<-rbind(periodStacked,periodPair[seq((1+(rep*(i-1))),rep*i,1)])
  phaseStacked<-rbind(phaseStacked,phasePair[seq((1+(rep*(i-1))),rep*i,1)])
  ampStacked<-rbind(ampStacked,ampPair[seq((1+(rep*(i-1))),rep*i,1)])
}
#relabel with treatments
Geoff_Pe<-cbind(OGTreat,as.data.frame(periodStacked))
Geoff_Ph<-cbind(OGTreat,as.data.frame(phaseStacked))
Geoff_Am<-cbind(OGTreat,as.data.frame(ampStacked))

#make labels
label3<-c("Treatment","Replicate")
label2<-c("Treatment")
for(i in 1:rep){
  repNum<-paste("Rep",i)
  label2<-cbind(label2,repNum)
}
#unlabel Geoff and label2 columns
colnames(Geoff_Pe)<-NULL
colnames(Geoff_Ph)<-NULL
colnames(Geoff_Am)<-NULL
colnames(label2)<-NULL
#relabel with what each column is
colnames(Geoff_Pe)<-as.data.frame(label2)
colnames(Geoff_Ph)<-as.data.frame(label2)
colnames(Geoff_Am)<-as.data.frame(label2)

#Prism would want the data transposed
PrismPeriod<-t(Geoff_Pe)
PrismPhase<-t(Geoff_Ph)
PrismAmplitude<-t(Geoff_Am)

#rownames and numbers
label2<-c("Treatment")
for(i in 2:nrow(PrismPeriod)){
  repnum<-paste("Rep",(i-1))
  label2<-rbind(label2,repnum)
}

#alphabetize data
#Alpha3<-PrismPeriod[order(PrismPeriod[1,]),]
#Alpha4<-PrismPhase[order(PrismPhase[1,]),]
#Alpha5<-PrismAmplitude[order(PrismAmplitude[1,]),]

#get rid of colnames and rownames for cleaner files
colnames(PrismPeriod)<-NULL
row.names(PrismPeriod)<-label2
colnames(PrismPhase)<-NULL
row.names(PrismPhase)<-label2
colnames(PrismAmplitude)<-NULL
row.names(PrismAmplitude)<-label2

#name outfiles
ppe<-paste(plate_num,"PrismPeriod.csv")
pph<-paste(plate_num,"PrismPhase.csv")
pam<-paste(plate_num,"PrismAmplitude.csv")
#write outfiles
write.csv(x = PrismPeriod, ppe)
write.csv(x = PrismPhase, pph)
write.csv(x = PrismAmplitude, pam)
#move to a subfolder
file.move(ppe, newloco2, overwrite = TRUE)
file.move(pph, newloco2, overwrite = TRUE)
file.move(pam, newloco2, overwrite = TRUE)

#############most csv outfiles and outfile folder#############
#delete metaout folder
unlink('metaout', recursive = TRUE)

#add individual well stats
otherout <- paste(plate_num[1], '96 Well Individual PerPhaAmp.csv')
write.csv(x = ind_wells, otherout)

#move files to new folder
file.move(otherout, newloco2, overwrite = TRUE)

#move any ANOVA files to ANOVA folder
if(ANOVA){
  file.move(title2, ANOVAFolder, overwrite = TRUE)
  file.move(file4, ANOVAFolder, overwrite = TRUE)

#delete the old file for the complicated ANOVA results otheroutA
  file.remove(otheroutA)
}

#make the Period, Phase, and Amplitude averages outfile
final <- rbind(Average, pretty_avg, SEM, pretty_sd)
#row names set-up
rownum<-(nrow(final)-2)
treatnum<-(rownum/2)
#set up for Average portion
finalLAB<-c(".")
for(i in 1:treatnum){
  numlet<-paste(i,"Avg")
  finalLAB<-rbind(finalLAB,numlet)
}
#set up for SEM portion
finalLAB<-rbind(finalLAB,"..")
for(i in 1:treatnum){
  numlet<-paste(i,"Sem")
  finalLAB<-rbind(finalLAB,numlet)
}
rownames(final)<-finalLAB
#set up name of outfile
out_file <- paste('#1 __', plate_num[1], 'Avg Per_Pha_Amp.csv')
#write outfile
write.csv(x = final, file= out_file)
file.move(out_file, newloco, overwrite = TRUE)

################loop controls for t-test#######################
#loop for all controls on the plate (shown: run each treatment as a control)
pcollect <- vector()
averer <- aver_source[,-1]

if (pvalues) {
  for (z in 1:length(treatment)) {
    #if we're in pick-a-control mode, and this isn't it, `next`
    if (treatment[z] %in% CONTROL) {
      if (duplicated(treatment)[z]) {
        next
      }
    } else if (length(CONTROL) > 1 ||
               (length(CONTROL) == 1 && CONTROL[1] != "")) {
      next
    }

    GetNumbr <- function(z) {
      #each treatment as control for the plate
      if (length(treatment) == 8) {
        ###8 treatments, 12 rep
        sample_control <- seq((z - 1) * 12 + 1, 12 * z, 1)
      } else {
        ###12 treatments, 8 rep
        sample_control <- seq(z, 84 + z, 12)
      }

      #set number for control based on treatment number
      if (length(treatment) == 8) {
        ###8 treatments, 12 rep
        len <- length(sample_control)
        numbr <- ((sample_control[len]) / len)
      } else {
        ###12 treatments, 8 rep
        numbr <- sample_control[1]
      }
      return(numbr)
    }

    #find duplicates, make numbrs an array of numbr
    numbrs <- c()
    for (cur_z in which(treatment %in% treatment[z])) {
      numbrs <- append(numbrs, GetNumbr(cur_z))
    }

    #######find pvalue here. then organize
    combined_p_values <- "p-value"

    #do this for all col in averer
    count1 = 0
    for (x in 1:length(averer[1, ])) {
      if (count1 == length(PPA_range)) {
        count1 = 0
      }
      count1 = count1 + 1
      temp2 <- GetDataForAvererColumn(averer, x, length(PPA_range))

      #need to convert control sample to new format
      control <- c()
      for (numbr in numbrs) {
        sample1 <- (((numbr - 1) * length(PPA_range)) + (count1))
        control <- append(control, na.omit(averer[, sample1]))
      }
      #does t test control vs treatment and extracts p-value
      if(pairwise==FALSE){
        #unequal variance, alt=means are not the same (two sided)
        #Welch two sample t-test
        #p-value adj method: bonferroni
        p_value <- t.test(control, temp2, var.equal = F, p.adjust.method="bonferroni")$p.value
      }
      if(pairwise==TRUE){
        #every treatment must be same-length
        data5<-c()

        if(paired){
          #get data into the same space and differentiate it
          data5<-rbind(control,temp2)
          data5<-t(data5)
          colnames(data5)<-c("Sample1","Sample2")

          #pairwise comparisons using t-tests with paired groups (no pooled SD)
          #p-value adj method: bonferroni
          p_value <- pairwise.t.test(data5, c("Sample1","Sample2"), p.adjust.method="bonferroni", var.equal = F, paired=TRUE, alternative=c("two.sided"))$p.value
          }
        else{
          #get data for the test
          data5<-cbind(t(control),t(temp2))
          data5<-t(data5)

          colin<-c()
          for (i in 1:length(control)){
            colin<-rbind(colin,"Sample1")
          }
          for (i in 1:length(temp2)){
            colin<-rbind(colin,"Sample2")
          }
          data5<-cbind(colin,data5)
          #retrieve colin column
          dataF<-data5[,1]

          #pairwise comparisons using t-tests with pooled SD
          #p-value adj method: bonferroni
          p_value <- pairwise.t.test(as.numeric(data5[,2]), data5[,1], p.adjust.method="bonferroni", var.equal = F, alternative=c("two.sided"))$p.value
        }
      }

      #collect the calculations
      #gathers p-values for all values
      combined_p_values <- cbind(combined_p_values, p_value)
    }

    #prep to reorganize averages so that every seventh (not counting the first one/label) is stacked
    pretty_pv <- PPA_range

    #take out label col
    com_pv <- as.matrix(combined_p_values[,-1])
    com_pv <- as.numeric(com_pv)

    #actually reorganize
    for (i in 1:(length(treatment))) {
      rpval <-
        c(com_pv[seq((i - 1) * length(PPA_range) + 1, length(PPA_range) * i, 1)])
      pretty_pv <- rbind(pretty_pv, rpval)
    }
    #redo col label
    colnames(pretty_pv) <- PPA_range
    pretty_pv <- pretty_pv[-1, ]
    #replace rownames as columns to avoid confusing numbering of repeat rownames
    rownames(pretty_pv)<-NULL
    colnames(pretty_pv)<-NULL

    rownames(pretty_pv)<-treatment
    pretty_pv <- pretty_pv[!duplicated(rownames(pretty_pv)), ]
    pvrow<-rownames(pretty_pv)
    rownames(pretty_pv)<-NULL
    pretty_pv<-cbind(pvrow,pretty_pv)

    colnames(pretty_pv)<-hoho2

    #collect iterations
    pcollect <- rbind(P_value, pretty_pv)

    #write outfiles with p-values if multiple CONTROLS
    if (length(pcollect) > 0) {
      #see if CONTROL is specified or all treatments
      #is it all treatments (unspecified)?
      if (length(CONTROL) == 1 && CONTROL[1] == "") {

        #making outfile
        #set up name of outfile
        out_file<-vector()
        #label p-value files with numbers if the user wants
        if(NumberLabel==TRUE){
          number<-paste("#",cur_z,sep="")
          out_file <- paste(number,treatment[cur_z],'p_value.csv')
        }
        else{
          out_file <- paste(treatment[cur_z],'p_value.csv')
        }
        out_file<-out_file[1]

        #write outfiles
        final<-rbind(Average,as.data.frame(pretty_avg),SEM,as.data.frame(pretty_sd),P_value,pretty_pv)
        #row names set-up
        rownum<-(nrow(final)-3)
        treatnum<-(rownum/3)
        #set up for Average portion
        finalLAB<-c(".")
        for(i in 1:treatnum){
          numlet<-paste(i,"Avg")
          finalLAB<-rbind(finalLAB,numlet)
        }
        #set up for SEM portion
        finalLAB<-rbind(finalLAB,"..")
        for(i in 1:treatnum){
          numlet<-paste(i,"Sem")
          finalLAB<-rbind(finalLAB,numlet)
        }
        #set up for p-value portion
        finalLAB<-rbind(finalLAB,"...")
        for(i in 1:treatnum){
          numlet<-paste(i,"Pval")
          finalLAB<-rbind(finalLAB,numlet)
        }
        rownames(final)<-finalLAB
        #write outfile
        write.csv(x=final,out_file)
        #move to a subfolder for t-test results
        file.move(file=out_file, t_testFolder, overwrite=TRUE)
      }
    }
  }
}

#write outfile
#are the p-values wanted?
if (length(pcollect) > 0) {
  #is the CONTROL specified?
  if (length(CONTROL) > 1 ||
      (length(CONTROL) == 1 && CONTROL[1] != "")) {
    final <- rbind(Average, pretty_avg, SEM, pretty_sd, pcollect)
    #set up name of outfile
    out_file <- paste(CONTROL,'t_tests.csv')
    #write outfile
    write.csv(x = final, out_file)
    #move to a subfolder
    file.move(out_file, t_testFolder, overwrite = TRUE)
  }
}

###################GRAPHS#####################################################################
#outputs graphs as PDFs

#function that rounds up to nearest multiple of base
mround <- function(x, base) {
  base * ceiling(x / base)
}

#function checks for integers. Found online.
is.wholenumber <-
  function(x, tol = .Machine$double.eps ^ 0.5)
    abs(x - round(x)) < tol

#get actual timepoints
timepoints <- GetTimepoints(names)

####cyc_1 cannot handle skipping an hour (only takes time series). until we fix this, we will have to make sure to not use data that skip#####
#function to check for sequential-ness. found online
is.sequential <- function(x) {
  all(diff(x) == diff(x)[1])
}

#if not sequential, stop
if (is.sequential(timepoints) != TRUE) {
  stop(
    "Time points are not consecutive. Please make sure there are no time skips between hours."
  )
}

#to reflect true phase, we shift time points to relative time (ZT)
skip<-roundTP
ZTtim <- GetZTtim(timepoints, ZT, skip)

###For Graph of Per,Pha,Amp##################################
#make file name
graph <- paste("#2 __", plate_num[1], 'Graphs.pdf')

###For Graph of Luminescence Curves############################
#average values for each time point then make graphs
#dataset for graphing
dat <- read.csv(plate)
xnum <- dat[, 1]
#taking out "sample x" column
dat <- dat[, -1]

#change into indiv time points for each treatment (this is specific to treatment number)
breo <- c(seq(1:rep))
for (i in 1:length(treatment)) {
  bery <- vector()
  for (x in 1:length(names)) {
    well_seq <-
      if (length(treatment) == 8)
        seq((i - 1) * 12 + 1, 12 * i, 1)
    else
      seq(i, 84 + i, 12)
    bver <- unlist(dat[well_seq, x])
    if (x == 1) {
      bery <- c(bver)
    }
    if (x != 1) {
      bery <- cbind(bery, bver)
    }
  }
  breo <- cbind(breo, bery)
}
breo <- breo[, -1]

#average each column (data for each time point)
averageTP <- vector()

for (x in 1:length(colnames(breo))) {
  #get rid of NA values
  temp3 <- vector()
  temp4 <- vector()
  temp3 <- breo[, x]
  temp4 <- na.omit(temp3)

  #do calculations
  av_time <- mean(temp4)
  averageTP <- cbind(averageTP, av_time)
}

#actually reorganize for graph
graphData <- vector()

for (i in 1:length(treatment)) {
  reorgan <-
    unlist(averageTP[seq((i - 1) * length(names) + 1,
                         (i - 1) * length(names) + length(names),
                         1)])
  graphData <- rbind(graphData, reorgan)
}

graphData <- graphData[!duplicated(treatment),]

error_bars <- c()
for (i in 1:length(row.names(graphData)))
{
  #get the error bars
  well_seq <- GetWellSeq(i, treatment, orientation)
  error_bar_points <- time_points[well_seq,]
  error_bars1<-c()
  for (j in 1:ncol(error_bar_points)) {
    error_bar_cur_column <- error_bar_points[,j]
    error_bar_value <-
      sd(error_bar_cur_column) / sqrt(length(error_bar_cur_column))
    error_bars1 <- append(error_bars1, error_bar_value)
  }
  error_bars<-rbind(error_bars,error_bars1)
}

#Print averages for each timepoint so you can make/alter graphs more easily outside of R
#name for file
DataToGraph <- paste("#3 __", plate_num[1], 'Averaged LUC Data.csv')
#labels for data
treatment_dedup <- OGTreat[!duplicated(OGTreat)]

rownames(graphData) <- treatment_dedup
colnames(graphData) <- ZTtim

rownames(error_bars) <- treatment_dedup

LumSEM<-c()
for(i in ZTtim){
  LumSEM <- c(LumSEM,'.')
}
#makes file
write.csv(x = rbind(graphData,LumSEM,error_bars), DataToGraph)
#move file
file.move(DataToGraph, newloco, overwrite = TRUE)

#get individual reps per time point in an easy to read manner for future graphing
repCollect<-c()
if (length(treatment) == 8){
  for(k in 1:length(treatment)){
    wellseq<-seq(((k-1)*12)+1,12*k,1)
    for(i in 1:length(wellseq)){
      repCollect<-rbind(repCollect,dat[wellseq[i],])
    }
  }
}
if (length(treatment) == 12){
  for(k in 1:length(treatment)){
    wellseq<-seq(k, 84 + k, 12)
    for(i in 1:length(wellseq)){
      repCollect<-rbind(repCollect,dat[wellseq[i],])
    }
  }
}

#Print reps for each timepoint so you can make/alter graphs more easily outside of R
colnames(repCollect)<-ZTtim
repCollect<-cbind(Treatments,repCollect)

#find number of unique treatments
Treatments00<-unique(repCollect[,1])
Treatments2<-repCollect[,1]

#move same treatment names so that they are near one another (just a preference)
if(length(Treatments00) != length(treatment)){
  potkettle2<-c()
  blacklist2<-c()
  for(i in 1:length(Treatments2)){
    #find number of replicates per treatment and locations
    Repeats<-which(Treatments2==Treatments2[i])
    #get all reps for a single treatment
    if(i==1){
      count=0
      for(k in Repeats){
        if(count==0){
          potkettle2<-c(repCollect[k,])
          count=count+1
        }
        else{
          potkettle2<-rbind(potkettle2,repCollect[k,])
        }
      }
    }
    else{
      #make sure we aren't repeating treatment replicates
      if((Treatments2[i] %in% blacklist2)==FALSE){
        for(p in Repeats){
          potkettle2<-rbind(potkettle2,repCollect[p,])
        }
      }
    }
    blacklist2<-c(blacklist2,Treatments2[i])
  }
  repCollector<-potkettle2
}

#if the treatments are not repeated
if(length(Treatments00) == length(treatment)){
  repCollector<-repCollect
}

#alphabetize data
#Alpha2<-repCollector[order(repCollector[,1]),]

indie<-paste(plate_num,"LUC Replicates.csv")
write.csv(x=repCollector,indie)
file.move(indie, newloco2, overwrite = TRUE)

#########################################start pdf#######################################
#take out labels from the data to graph
pretty_avg<-pretty_avg[,-1]
pretty_sd<-pretty_sd[,-1]

if (GRAPHS) {

  #Open a pdf file
  pdf(graph, width = 9, height = 35)

  #figure out inner margins based on length of longest treatment name
  #find longest treatment name
  long<-0
  for (i in 1:length(treatment)){
    #if the treatment name is longest, keep it
    if(nchar(treatment[i])>long){
      long<-nchar(treatment[i])
    }
  }
  #once we have the longest, use it to calculate the inner margin we want below the graph
  #this will allow for enough room for the treatment labels to not run into the next graph
  below_margin<-(long/2)
  if(below_margin<6){
    below_margin<-6
  }

  #set parameters for the PDF of graphs
  par(
    #make 2 columns and TBD rows of graphs
    mfrow = c(ceiling(length(colnames(pretty_avg)) / 2) + ceiling(length(row.names(graphData)) / 2), 2),
    #axis label style: perpendicular to axis
    las = 2,
    #set outer margins (all graphs vs edge of pdf): below, left, top, right
    oma = c(2, 1, 2, 1),
    #set inner margins (graph vs graph): below, left, top, right
    mar=c(below_margin,5,6,4)
  )

  #########################Graph of Per,Pha,Amp###########
  #get units for y-axis
  if (preARS > 0) {
    yAx <-
      c("Period Length (Hours)",
        "Phase (Hours)",
        "Amplitude (Luminescence)")
  }

  for (i in 1:length(colnames(pretty_avg)))
  {
    #get upper y limit for axis
    maxP <- max(as.numeric(pretty_avg[, i]))
    lim <- c(0)
    lim <- maxP + (2 * (as.numeric(max(pretty_sd[, i]))))
    if (is.wholenumber(i / 30) == FALSE) {
      lim <- mround(lim, 5)
    }
    if (is.wholenumber(i / 30) == TRUE) {
      lim <- mround(lim, 200)
    }

      #get title for graph
      label <- hoho[i]
      title <- paste(plate_num[1], label)

    #graph
    center <-
      barplot(
        as.numeric(pretty_avg[, i]),
        ylab = yAx[i],
        ylim = c(0, lim),
        names.arg = treatment_dedup,
        col = "dark green",
        main = title
      )

    #error bars hack: we draw arrows but with very special "arrowheads"
    arrows(
      center,
      as.numeric(pretty_avg[,i]) - as.numeric(pretty_sd[,i]),
      center,
      as.numeric(pretty_avg[,i]) + as.numeric(pretty_sd[,i]),
      length = 0.05,
      angle = 90,
      code = 3,
      col = 'black'
    )
  }

  #############################Box Plots################
  # Set up graphing constants - this lets us look up the Y axis label and graph title based on what PPA range we're looking at
  BOX_PLOT_YLAB = list(
    "Hours",
    "Hours",
    "Luminescence"
  )
  names(BOX_PLOT_YLAB) = c(
    "ARS_period",
    "ARS_adjphase",
    "ARS_amplitude"
  )
  BOX_PLOT_MAIN = list(
    "ARS Period",
    "ARS AdjPhase",
    "ARS Amplitude"
  )
  names(BOX_PLOT_MAIN) = c(
    "ARS_period",
    "ARS_adjphase",
    "ARS_amplitude"
  )

  # Loop over PPA_range, look up what we need to graph, and graph it
  for (i in 1:length(PPA_range)) {
    curGraphKey <- PPA_range[i]
    curGroup <-
      averer[, seq(i,
                   length(treatment_dedup) * length(PPA_range),
                   length(PPA_range))]
    boxplot(
      curGroup,
      names = treatment_dedup,
      ylab = toString(BOX_PLOT_YLAB[curGraphKey]),
      main = toString(BOX_PLOT_MAIN[curGraphKey])
    )
  }

  #############################Circular Phase Plots#####################
  #loop for all treatments
  for (i in 1:length(treatment_dedup)){

    #use shifted values or at least shift to none over 24 hr
    if ('ARS_adjphase' %in% colnames(PerPhaAmp)){
      #if adjphase > 24, subtract by 24
      phase1 = ind_wells[rownames(ind_wells)=='ARS_adjphase',]
      phaseA = phase1[i,]
      phaseA = na.omit(phaseA)
      for (x in phaseA){
        #if a phase value is over 24, subtract 24 from it
        if (as.numeric(x)>24){
          y=as.numeric(x)-24
          #replace old value with new value
          phaseA[phaseA==x] <- y
        }
      }
    }

    #ARS Graph
    if ('ARS_adjphase' %in% colnames(PerPhaAmp)){
      #set the phase values as a circular object
      y <- as.circular(as.numeric(phaseA), zero=0,units='hours',rotation='clock',template='clock24',type='angles',modulo='asis')
      #Make a circular plot of the phase of each replicate per treatment
      title <- paste(treatment_dedup[i], "ARS Phase Distribution")
      plot.circular(y,main= title,ylab="Hours",units="hours",rotation="clock",template = 'clock24')
      #plot mean as arrow
      arrows.circular(mean(y), col = "red")
    }
  }

  #############################Graph of Curves################
  for (i in 1:length(row.names(graphData)))
  {
    #get upper y limit for axis
    maxP <- max(graphData[i,])
    lim <- c(0)
    if (maxP < 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }

    #get title for graph
    title <- paste(treatment_dedup[i], "Expression")

    #graph
    plot(
      x = ZTtim,
      y = graphData[i,],
      ylim = c(0, lim),
      xlim = c(0, length(names) + 24),
      xlab = "Time (hours)",
      ylab = "Luminescence",
      main = title,
      type = 'p',
      col = "Blue",
      xaxt = 'n'
    )
    lines(ZTtim, graphData[i,], type = 'o', col = 'blue')
    #x axis
    if (length(names) < 100) {
      axis(1, c(seq(12, length(names) + 24, 12)))
    }
    if (length(names) >= 100 &
        length(names) < 200) {
      axis(1, c(seq(12, length(names) + 24, 12)))
    }
    if (length(names) >= 200) {
      axis(1, c(seq(24, length(names) + 24, 24)))
    }
    #y axis
    if (2500 > lim &
        lim >= 1500) {
      axis(2, c(seq(250, mround(lim, 250), 250)))
    }
    if (lim >= 2500 &
        lim < 5000) {
      axis(2, c(seq(500, mround(lim, 500), 500)))
    }
    if (lim >= 5000 &
        lim < 10000) {
      axis(2, c(seq(1000, mround(lim, 1000), 1000)))
    }
    if (lim >= 10000 &
        lim < 20000) {
      axis(2, c(seq(2000, mround(lim, 2000), 2000)))
    }
    if (lim >= 20000 &
        lim < 45000) {
      axis(2, c(seq(5000, mround(lim, 5000), 5000)))
    }
    if (lim >= 45000) {
      axis(2, c(seq(10000, mround(lim, 10000), 10000)))
    }
  }

  #Do it again, but now with error bars
  for (i in 1:length(row.names(graphData)))
  {
    #get upper y limit for axis
    maxP <- max(graphData[i,])
    lim <- c(0)
    if (maxP < 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }

    #get title for graph
    title <- paste(treatment_dedup[i], "Expression")

    #graph
    plot(
      x = ZTtim,
      y = graphData[i,],
      xlim = c(0, length(names) + 24),
      ylim = c(0, lim),
      xlab = "Time (hours)",
      ylab = "Luminescence",
      main = title,
      type = 'p',
      col = "Blue",
      xaxt = 'n'
    )

    #add arrow bars [sic]
    arrows(
      ZTtim,
      graphData[i,] - error_bars[i,],
      ZTtim,
      graphData[i,] + error_bars[i,],
      length = 0.05,
      angle = 90,
      code = 3
    )

    lines(ZTtim, graphData[i,], type = 'o', col = 'blue')
    #x axis
    if (length(names) < 100) {
      axis(1, c(seq(12, length(names) + 24, 12)))
    }
    if (length(names) >= 100 &
        length(names) < 200) {
      axis(1, c(seq(12, length(names) + 24, 12)))
    }
    if (length(names) >= 200) {
      axis(1, c(seq(24, length(names) + 24, 24)))
    }
    #y axis
    if (2500 > lim &
        lim >= 1500) {
      axis(2, c(seq(250, mround(lim, 250), 250)))
    }
    if (lim >= 2500 &
        lim < 5000) {
      axis(2, c(seq(500, mround(lim, 500), 500)))
    }
    if (lim >= 5000 &
        lim < 10000) {
      axis(2, c(seq(1000, mround(lim, 1000), 1000)))
    }
    if (lim >= 10000 &
        lim < 20000) {
      axis(2, c(seq(2000, mround(lim, 2000), 2000)))
    }
    if (lim >= 20000 &
        lim < 45000) {
      axis(2, c(seq(5000, mround(lim, 5000), 5000)))
    }
    if (lim >= 45000) {
      axis(2, c(seq(10000, mround(lim, 10000), 10000)))
    }
  }

  ###now make a graph with all the lines
  maxP <- max(graphData)
  lim <- c(0)
  if (maxP < 500) {
    lim <- mround(maxP, 100)
  }
  if (maxP > 500) {
    lim <- mround(maxP, 500)
  }
  title <- paste("Expression of", plate_num[1])
  plot(
    x = ZTtim,
    y = graphData[1,],
    ylim = c(0, lim),
    xlim = c(0, length(names) + 24),
    xlab = "Time (hours)",
    ylab = "Luminescence",
    main = title,
    type = 'l',
    col = "Blue",
    xaxt = 'n'
  )
  lines(ZTtim, graphData[1,], type = 'l', col = 'blue')

  #set axis scale
  #x axis
  if (length(names) < 100) {
    axis(1, c(seq(12, length(names) + 24, 12)))
  }
  if (length(names) >= 100 &
      length(names) < 200) {
    axis(1, c(seq(12, length(names) + 24, 12)))
  }
  if (length(names) >= 200) {
    axis(1, c(seq(24, length(names) + 24, 24)))
  }
  #y axis
  if (2500 > lim &
      lim >= 1500) {
    axis(2, c(seq(250, mround(lim, 250), 250)))
  }
  if (lim >= 2500 &
      lim < 5000) {
    axis(2, c(seq(500, mround(lim, 500), 500)))
  }
  if (lim >= 5000 &
      lim < 10000) {
    axis(2, c(seq(1000, mround(lim, 1000), 1000)))
  }
  if (lim >= 10000 &
      lim < 20000) {
    axis(2, c(seq(2000, mround(lim, 2000), 2000)))
  }
  if (lim >= 20000 &
      lim < 45000) {
    axis(2, c(seq(5000, mround(lim, 5000), 5000)))
  }
  if (lim >= 45000) {
    axis(2, c(seq(10000, mround(lim, 10000), 10000)))
  }

  rain <-
    c(
      'black',
      'cyan',
      'green2',
      'dark green',
      'red',
      'orange',
      'darkorchid',
      'yellow1',
      'mediumpurple',
      'red4',
      'thistle',
      'sienna3',
      'turquoise4',
      'gray33',
      'slateblue3',
      'darkgoldenrod',
      'deeppink'
    )
  colors <- c('blue')
  for (i in 2:length(row.names(graphData)))
  {
    #get upper y limit for axis
    maxP <- max(graphData[i,])
    lim <- c(0)
    if (maxP < 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }

    #graph
    lines(
      x = ZTtim,
      y = graphData[i,],
      ylim = c(0, lim),
      xlim = c(0, length(names) + 20),
      xlab = "Time (hours)",
      ylab = "Luminescence",
      main = title,
      type = 'l'
    )
    lines(ZTtim, graphData[i,], type = 'l', col = rain[i])
    color <- c(rain[i])
    colors <- cbind(colors, color)
  }

  #make legend to side
    plot(
      0,
      0,
      type = "n",
      bty = "n",
      xaxt = "n",
      yaxt = "n",
      xlab = "",
      ylab = ""
    )
    legend(
      'topleft',
      legend = c(treatment_dedup),
      col = c(colors),
      pch = 10,
      bty = 'n',
      ncol = 1
    )


  # Close the pdf file
  dev.off()
}

#############move graph file into new folder#############
if(GRAPHS){
  file.move(graph, newloco, overwrite = TRUE)
}

###################################################
#########End of code###############################
###################################################
