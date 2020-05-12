# Amanda Wilson, Marco-Felipe King, Mark H. Weir

# SET PARENT DIRECTORY & LOAD/INSTALL PACKAGES -----

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    

if("truncdist" %in% rownames(installed.packages())==FALSE){install.packages("truncdist"); require(truncdist)}else{require(truncdist)}
if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}

# IMPORT DATA NEEDED FOR ALL OPERATIONS DURING AUTOMATION -----

  #load in duration data
  durations<-read.csv('ContactDuration.csv')
  #read in bootstrapped values for dose-response
  exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')

# CONSTANTS TO BE USED 

SIM.iter <- 100     # Making a master Monte Carlo iteration value
suppressMessages(source("adjust_behaviors_covid.R"))


# AUTOMATION ----
# Next line will be the names of the simulations that will be run. DCLP is Doffing correct and low patient numbers

#DP=doffing poor; DO=doffing ok; DG=doffing great
#LP=low prob COVID-19; HP=high prob COVID-19
#NS=non-surge (low # patient visits); S=surge (high # patient visits)
#control=1 patient visit, 100% prob of COVID-19 patient

SIM <- c("DPLPNS", "DPLPS","DPHPNS","DPHPS",
         "DOLPNS","DOLPS","DOHPNS","DOHPS",
         "DGLPNS","DGLPS","DGHPNS","DGHPS","control")   

NUM.SIM <- length(SIM)     # Count the number of iterations for the automated simulations

for(j in 1:NUM.SIM)
{
  sim.num <- j; sim.name <- SIM[j]
  
  if(sim.name=="DPLPNS"){if(dir.exists("DPLPNS")==FALSE){dir.create("DPLPNS"); setwd("DPLPNS")}else{setwd("DPLPNS")}}
  if(sim.name=="DPLPS"){if(dir.exists("DPLPS")==FALSE){dir.create("DPLPS"); setwd("DPLPS")}else{setwd("DPLPS")}}
  if(sim.name=="DPHPNS"){if(dir.exists("DPHPNS")==FALSE){dir.create("DPHPNS"); setwd("DPHPNS")}else{setwd("DPHPNS")}}
  if(sim.name=="DPHPS"){if(dir.exists("DPHPS")==FALSE){dir.create("DPHPS"); setwd("DPHPS")}else{setwd("DPHPS")}}

  if(sim.name=="DOLPNS"){if(dir.exists("DOLPNS")==FALSE){dir.create("DOLPNS"); setwd("DOLPNS")}else{setwd("DOLPNS")}}
  if(sim.name=="DOLPS"){if(dir.exists("DOLPS")==FALSE){dir.create("DOLPS"); setwd("DOLPS")}else{setwd("DOLPS")}}
  if(sim.name=="DOHPNS"){if(dir.exists("DOHPNS")==FALSE){dir.create("DOHPNS"); setwd("DOHPNS")}else{setwd("DOHPNS")}}
  if(sim.name=="DOHPS"){if(dir.exists("DOHPS")==FALSE){dir.create("DOHPS"); setwd("DOHPS")}else{setwd("DOHPS")}}
  
  if(sim.name=="DGLPNS"){if(dir.exists("DGLPNS")==FALSE){dir.create("DGLPNS"); setwd("DGLPNS")}else{setwd("DGLPNS")}}
  if(sim.name=="DGLPS"){if(dir.exists("DGLPS")==FALSE){dir.create("DGLPS"); setwd("DGLPS")}else{setwd("DGLPS")}}
  if(sim.name=="DGHPNS"){if(dir.exists("DGHPNS")==FALSE){dir.create("DGHPNS"); setwd("DGHPNS")}else{setwd("DGHPNS")}}
  if(sim.name=="DGHPS"){if(dir.exists("DGHPS")==FALSE){dir.create("DGHPS"); setwd("DGHPS")}else{setwd("DGHPS")}}
  
  if(sim.name=="control"){if(dir.exists("control")==FALSE){dir.create("control"); setwd("control")}else{setwd("control")}}
  
  
# MASTER FUNCTION ----
behavior.sim<-function(caretype=c("IV","Obs","Rounds"),numsequence,prob.patient.infect,numvisit,prob.contam.between){
  
  #numvisit is function of shift length and number of patients, # of care episodes
  
   set.seed(34)
   #require(truncdist)
  

  #Remove "S" at end of entries and convert from character to numeric format
  durations$Duration<-as.numeric(gsub("S","",durations$Duration))
 
  #Behavior model section

  sample.space<-c("Equipment","FarPatient","HygieneInside","In","NearPatient","Out","Patient")
  
  behavior.total<-list() #creating a list to store behaviors
  exposure.frame<-list() #creating a list to store exposure.frame data.frames
  
  #set up matrix for type of care

    if (caretype=="IV"){
      prob.mat<-TIV$estimate #IV
    }else if (caretype=="Obs"){
      prob.mat<-TObs$estimate #Obs
    }else{
      prob.mat<-TRounds$estimate #Rounds
    }
 
  

  for (j in 1:numsequence){
    
    for (m in 1:numvisit){
      
      behavior<-"In" #first behavior is "In"
      k<-2 #setting up k for while loop
      
      #while the previous behavior is not "Out," make new behaviors
      while (behavior[k-1]!="Out"){
        
        #If previous step was not out, create new behavior
        #where the prob of selecting a behavior from sample.space is determined by the row
        #of probabilities corresponding to the previous behavior contact (prob.mat[behavior[k-1],])
        
        behavior[k]<-sample(sample.space,1,replace=TRUE,prob=prob.mat[behavior[k-1],])
        
        #now we have to handle issues where someone is taking GlovesOff having never put them on
        
        #if the current selected behavior is GlovesOff and they've never put GlovesOn or the maximum position of a previous GlovesOff is greater than
        # the maximum position of a previous GlovesOn
        behaviorcount<-1:length(behavior)
        
        #advance contact number by 1
        k<-k+1
      }
      
      
      ###  exposure simulation  ###
      
      #--------- TRANSFER EFFICIENCY -------------------------------------------------------
      
      #initialize transfer efficiencies
      #lambda = hand--> surf
      #beta = surf --> hand
      
      lambda<-rep(NA,length(behavior))
      beta<-rep(NA,length(behavior))
      
      #UPDATE
      
      #transfer efficiencies for non-patient surf
      
      #lambda hand--> surface and beta surface --> hand
      lambda[behavior!="Patient"]<-runif(length(lambda[behavior!="Patient"]),0.0003,.217)
      beta[behavior!="Patient"]<-runif(length(lambda[behavior!="Patient"]),0.0003,.217)
      
      #transfer efficiency patient skin contacts
      lambda[behavior=="Patient"]<-rtrunc(length(lambda[behavior=="Patient"]),spec="norm",a=0,b=1,mean=0.056,sd=0.032)
      beta[behavior=="Patient"]<-rtrunc(length(lambda[behavior=="Patient"]),spec="norm",a=0,b=1,mean=0.056,sd=0.032)
      
      #-------------- SURFACE CONCENTRATIONS -----------------------------------------------------------------------------
      
      #Initialize surface concentration (PLACE HOLDER VALUE RIGHT NOW.. WILL DO LIT REVIEW TO INFORM THIS PARAMETER)
      surfconc<-rep(NA,length(behavior))
      
      #concentrations for hand-to-surface contact moments
      numtest<-runif(1,0,1)
      if (numtest<=prob.patient.infect){
       
         #concentrations assuming patient is infected and asymptomatic
        
        #Assume Out conc ~ In conc
        surfconc[behavior!="patient"]<-10^runif(length(surfconc[behavior!="patient"]),0,3) #place holder conc
        
        #Patient surfaces
        surfconc[behavior=="Patient"]<-10^runif(length(surfconc[behavior=="Patient"]),3,6) #place holder conc
          
      }else{
        
        #concentrations assuming patient is not infected
        
        #Assume Out conc ~ In conc
        surfconc[behavior!="patient"]<-10^runif(length(surfconc[behavior!="patient"]),-2,0)
                 
        #Patient surfaces
        surfconc[behavior=="Patient"]<-rep(0,length(surfconc[behavior=="Patient"]))
      }
        
      
      #-------------- FRACTIONAL HAND SURFACE AREA ------------------------------------------------------------------
      
      #Initialize fraction of hand used for contact
      SH<-rep(NA,length(behavior))
      
      #fractional surface area for open hand grip; from AuYeung
      SH[behavior=="In" | behavior=="Out"]<-runif(length(SH[behavior=="In" | behavior=="Out"]),0.10,0.17) #min and max of left and right hands closed hand grip in AuYeung et al. (2008)
      
      #fractional surface area for patient contact (front partial fingers to full front palm with fingers)
      SH[behavior=="Patient"]<-runif(length(SH[behavior=="Patient"]),0.04,0.25) 
      
      #fractional surface area for variety of grip types (non "in"/"out" contacts)
      SH[behavior=="Equipment"|behavior=="FarPatient"|behavior=="NearPatient"|behavior=="HygieneInside"]<-runif(length(SH[behavior=="Equipment"|behavior=="FarPatient"|behavior=="NearPatient"|behavior=="HygieneInside"]),0.04/5,0.25)
      #min and max of left and right hands in AuYeung et al. (2008) for various hand grip and hand press contacts (hand immersion contacts not included)
      #from single fingertip up to full palm
      
      #------------- RIGHT HAND VS. LEFT HAND ---------------------------------------------------------------------------
      
      #initialize R or L hand
      hand<-rep(NA,length(behavior))
      hand[sample(1:length(behavior),length(behavior)/2)]<-"right"
      hand[is.na(hand)]<-"left"
      
      #------------- INACTIVATION CONSTANTS -----------------------------------------------------------------------------
      
      #surface
      t99.s<-runif(length(behavior),3,5*24)*3600 #T99 in seconds
      k.s<-(-log(1/(10^2))/t99.s)
      
      t99.h<-runif(length(behavior),1,6)*3600 #T99 in seconds
      k.h<-(-log(1/(10^2))/t99.h)
      
      #-------------- EXPOSURE SIMULATION ------------------------------------------------------------------------------
      
      #initialize conccentration on R and L hand estimates
      handR<-rep(0,length(behavior))
      handL<-rep(0,length(behavior))
      
      #placeholder for now in case we want to adjust duration
      duration<-sample(durations$Duration,length(behavior),replace=TRUE)
      
      if (m==1){
        if(hand[1]=="right"){
          handR[1]<-beta[1]*SH[1]*(surfconc[1]*exp(-k.s[1]*duration[1]))
          handL[1]<-0  
        }else{ #if hand[1] == left...
          handL[1]<-beta[1]*SH[1]*(surfconc[1]*exp(-k.s[1]*duration[1]))
          handR[1]<-0 
        }
      }else{
        handR[1]<-handRnext
        handL[1]<-handLnext
      }
      
      for (a in 2:(length(behavior))){
        
        if(hand[a]=="right"){
          
            handR[a]<-(handR[a-1]-(lambda[a]*SH[a]*handR[a-1]*exp(-k.h[a]*duration[a]))+(beta[a]*SH[a]*surfconc[a]*exp(-k.s[a]*duration[a])))
            handL[a]<-handL[a-1]

        }else{ #if hand[a] == left...
          
            handL[a]<-(handL[a-1]-(lambda[a]*SH[a]*handL[a-1]*exp(-k.h[a]*duration[a]))+(beta[a]*SH[a]*surfconc[a]*exp(-k.s[a]*duration[a])))
            handR[a]<-handR[a-1]
        }
      
      
      }
      
      
      #---------------------------------- CALCULATION OF INFECTION RISK FOR CARE EPISODE --------------------------------------------------------------------------------------------------


      
      # require('gsl')
        
        #final concentration on hands a episode of care
        handRnext<- handR[a]
        handLnext<- handL[a]
        
        #inactivation concentration on hands
        khtemp<- k.h[a]
        
        #doffing....
        
        #1/10 participants had <0.00003% of original inoculum transfered to hands during doffing
        
        randomnum<-runif(1,0,1)
        
        if(randomnum<=prob.contam.between){
          
          handRnext<-handRnext*(3*10^-7) #assume % of viral titer in study represents transfer potential when doffing
          handLnext<-handLnext*(3*10^-7) #assume % of viral titer in study represents transfer potential when doffing
          
          #ethanol hand rub (described by Kampf et al. 2020 as recommended step after doffing)
          washnum<-runif(1,0,1)
          
          if (washnum<=.5){
            #then they wash their hands
            handRnext<-handRnext*(1/10^rtrunc(1,spec="norm",a=0,b=6,mean=1.62,sd=0.12)) #based on dist from Marco
            handLnext<-handLnext*(1/10^rtrunc(1,spec="norm",a=0,b=6,mean=1.62,sd=0.12)) #based on dist from Marco
          }else{
            #otherwise they use hand sanitizer
            handRnext<-handRnext*(1/10^runif(1,2,4)) #based on Kampf et al. (2020)
            handLnext<-handLnext*(1/10^runif(1,2,4)) #based on Kampf et al. (2020)
          }
          
          #Determines which hand will be used for the hand-to-face contact
          whichhand<-runif(1,0,1)
          
          #transfer efficiency to mouth
          TEtemp<-.3390 #Rusin et al. (2002)
          
          #fraction of hand used in mouth contact
          SMtemp<-runif(1,0.04/5, 0.06/5) #single fingertip surface area fraction = front partial fingers fractional surface area / 5
          
          #hand surface area
          AHtemp<-runif(1,445,535) #Beamer et al. (2015) office model and from Exposure Factors Handbook
          
          if (m==1){
            #if this is the first episode of care, this is our first dose
            if (whichhand<=.5){
              dose<-handLnext*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handLnext<-(1-TEtemp*SMtemp)*handLnext*exp(-khtemp)
            }else{
              dose<-handRnext*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handRnext<-(1-TEtemp*SMtemp)*handRnext*exp(-khtemp)
            }
          }else{
            #otherwise, we're adding onto our cumulative dose
            if (whichhand<=.5){
              dose<-dose+handLnext*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handLnext<-(1-TEtemp*SMtemp)*handLnext*exp(-khtemp)
            }else{
              dose<-dose+handRnext*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handRnext<-(1-TEtemp*SMtemp)*handRnext*exp(-khtemp)
            }
          }
         
          
        }else{ #(if they do not self contaminate...)
          
          #if no contam while doffing, then handwashing wouldn't change concentration
          #and therefore would not change dose
          if (m==1){
            #if this is then the first episode of care, dose would be zero
            dose<-0
          }else{
            #otherwise, it's equal to the previous dose (no change in cumulative dose)
            dose<-dose
          }
          handRnext<-0 #no contamination for next care episode
          handLnext<-0 #no contamination for next care episode
          
        }
        pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
        infecttemp<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose, give=FALSE, strict=TRUE)
        
       if(infecttemp==0){
         infecttemp<-1*10^-15 #cannot have zero infection risk, so replace with small risk
       }
      
      doserep<-rep(dose,length(behavior))
      infectrep<-rep(infecttemp,length(behavior))
      
      #saving concentrations for all rooms
      if (m==1){
        doseall<-doserep
        infectall<-infectrep
        handRall<-handR
        handLall<-handL
        handall<-hand
        behaviorall<-behavior
        SHall<-SH
        lambdaall<-lambda
        betaall<-beta
        surfconcall<-surfconc
        patientnum<-rep(m,length(handR))
        k.sall<-k.s
        k.hall<-k.h
        durationall<-duration
      }else{
        doseall<-c(doseall,doserep)
        infectall<-c(infectall,infectrep)
        handRall<-c(handRall,handR)
        handLall<-c(handLall,handL)
        handall<-c(handall,hand)
        behaviorall<-c(behaviorall,behavior)
        SHall<-c(SHall,SH)
        lambdaall<-c(lambdaall,lambda)
        betaall<-c(betaall,beta)
        surfconcall<-c(surfconcall,surfconc)
        patientnum<-c(patientnum,rep(m,length(handR)))
        k.sall<-c(k.sall,k.s)
        k.hall<-c(k.hall,k.h)
        durationall<-c(durationall,duration)
      }
      
    }
    
  
    # -------------------------------- SAVE OUTPUT FOR SIMULATION FOR SINGLE PERSON ----------------------------------------------------------------------------------
    exposure.frame.temp<-data.frame(dose=doseall,infect=infectall,patientnum=patientnum,handR=handRall,handL=handLall,hand=handall,behavior=behaviorall,duration=durationall,SH=SHall,lambda=lambdaall,beta=betaall,surfconc=surfconcall,k.sall=k.sall,k.hall=k.hall)
    behavior.total[[j]]<-behavior
    exposure.frame[[j]]<-exposure.frame.temp

    rm(behavior)
  }
  # --------------------------------- SAVE ALL OUTPUTS TO GLOBAL ENV --------------------------------------------------------------------------------------------------
  behavior.total<<-behavior.total
  exposure.frame<<-exposure.frame
  
  
} # Master function end


#dawn poor, dawn ok, dawn good
probcontam<-c(0.8,0.5,0.1)

#low patient load, high patient load
patload<-c(15,30)

#non-surge, high surge
contam<-c(0.05,0.5)

if (j<5){
  prob.contam.between=probcontam[1]
}else if (j>4 & j<9){
  prob.contam.between=probcontam[2]
}else{
  prob.contam.between=probcontam[3]
}

if(j %% 2 ==0){
  prob.patient.infect=contam[1]
}else{
  prob.patient.infect=contam[2]
}

if (sim.name=="DPLPNS" | sim.name=="DPLPS" | sim.name=="DOLPNS" | sim.name=="DOLPS"|
    sim.name=="DGLPNS" | sim.name=="DGLPS"){
  numvisit=patload[1]
}else{
  numvisit=patload[2]
}

#IV scenario
behavior.sim(caretype="IV",numsequence=500,prob.patient.infect=prob.patient.infect,
             numvisit=numvisit,prob.contam.between=prob.contam.between)
IV<-exposure.frame
      
#Rounds scenario
behavior.sim(caretype="Rounds",numsequence=500,prob.patient.infect=prob.patient.infect,
             numvisit=numvisit,prob.contam.between=prob.contam.between)
Rounds<-exposure.frame
      
#Observational care scenario
behavior.sim(caretype="Rounds",numsequence=500,prob.patient.infect=prob.patient.infect,
             numvisit=numvisit,prob.contam.between=prob.contam.between)
Obs<-exposure.frame
      
#save output
saveRDS(IV, file=sprintf("%s.IV.exposure.frame.rds",sim.name))
saveRDS(Rounds, file=sprintf("%s.IV.exposure.frame.rds",sim.name))
saveRDS(Obs, file=sprintf("%s.IV.exposure.frame.rds",sim.name))

  #reset directory to parent folder so we can go to correct subfolder within parent folder for next sim run
  setwd(this.dir)
  
  } # Automation loop end