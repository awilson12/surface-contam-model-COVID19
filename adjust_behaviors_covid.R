
#RUN THIS PIECE FIRST BEFORE RUNNING FUNCTION CODE

# using Marco's Markov chains to simulate different contact behaviors
if("readr" %in% rownames(installed.packages())==FALSE){install.packages("readr"); require(readr)}else{require(readr)}
if("plyr" %in% rownames(installed.packages())==FALSE){install.packages("plyr"); require(plyr)}else{require(plyr)}
if("markovchain" %in% rownames(installed.packages())==FALSE){install.packages("markovchain"); require(markovchain)}else{require(markovchain)}

#library(readr)
movsdf.rbind <- read_csv("movsdf.rbind_orientationcorrected.csv")

for(i in 1:length(table(movsdf.rbind$ActivityID))){
  if (i==1){
    movsdf.rbind.new<-movsdf.rbind[movsdf.rbind$ActivityID==i & movsdf.rbind$Surface!="AlcOutside" &
                                     movsdf.rbind$Surface!="ApronOn" & movsdf.rbind$Surface!="ApronOff" &
                                     movsdf.rbind$Surface!="GlovesOn" & movsdf.rbind$Surface!="GlovesOff"&
                                     movsdf.rbind$Surface!="Alc" ,]
  }else{
    movsdf.rbindtemp<-movsdf.rbind[movsdf.rbind$ActivityID==i & movsdf.rbind$Surface!="AlcOutside" &
                                     movsdf.rbind$Surface!="ApronOn" & movsdf.rbind$Surface!="ApronOff" &
                                     movsdf.rbind$Surface!="GlovesOn" & movsdf.rbind$Surface!="GlovesOff"&
                                     movsdf.rbind$Surface!="Alc",]
    movsdf.rbind.new<-rbind(movsdf.rbind.new,movsdf.rbindtemp)
  }
}


#####
# 2.3 Aggregating surfaces into categories for Transition Matrices
#detach("package:dplyr", unload = TRUE)
#library(plyr)
movsdf.rbind.new$SurfaceCategories<-revalue(movsdf.rbind.new$Surface,c("In"="In",
                                                                       "Door"="FarPatient",
                                                                       "Other"="FarPatient",
                                                                       "Table"="NearPatient",
                                                                       "Tray"="Equipment" ,
                                                                       "Patient"="Patient",
                                                                       "Sharps"="Equipment"  ,
                                                                       "Waste"="HygieneInside" ,
                                                                       "Sink"="HygieneInside",
                                                                       "Soap"="HygieneInside" ,
                                                                       "PaperTowel"="HygieneInside" ,
                                                                       "IV"="Equipment" ,
                                                                       "Out"="Out" ,
                                                                       #"Alc"="Alcohol" ,
                                                                       "ObsTrolley"="Equipment",
                                                                       "Wipes"="HygieneInside",
                                                                       "Bed"="NearPatient",
                                                                       "Chair"="NearPatient",
                                                                       "Stethoscope"="Equipment"))

# Creating lists so we don't count transitions between "out" and first state of next observed episode
# It looksl ike Activity ID is unique per observation

TObs.left.list<-list()
TObs.right.list<-list()
TIV.left.list<-list()
TIV.right.list<-list()
TRounds.left.list<-list()
TRounds.right.list<-list()

for (i in 1:max(movsdf.rbind.new$ActivityID)){
  TObs.left.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="Obs" & movsdf.rbind.new$Orientation=="leftFacing" & movsdf.rbind.new$ActivityID==i]
  TObs.right.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="Obs" & movsdf.rbind.new$Orientation=="rightFacing" & movsdf.rbind.new$ActivityID==i]
  
  TIV.left.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="IV" & movsdf.rbind.new$Orientation=="leftFacing" & movsdf.rbind.new$ActivityID==i]
  TIV.right.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="IV" & movsdf.rbind.new$Orientation=="rightFacing" & movsdf.rbind.new$ActivityID==i]
  
  TRounds.left.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="Rounds" & movsdf.rbind.new$Orientation=="leftFacing" & movsdf.rbind.new$ActivityID==i]
  TRounds.right.list[[i]]<-movsdf.rbind.new$SurfaceCategories[movsdf.rbind.new$CareType=="Rounds" & movsdf.rbind.new$Orientation=="rightFacing" & movsdf.rbind.new$ActivityID==i]
  
}

require(markovchain)
TObs.left<-markovchainFit(TObs.left.list)
TObs.right<-markovchainFit(TObs.right.list)

TIV.left<-markovchainFit(TIV.left.list)
TIV.right<-markovchainFit(TIV.right.list)

TRounds.left<-markovchainFit(TRounds.left.list)
TRounds.right<-markovchainFit(TRounds.right.list)
