# Amanda Wilson, Marco-Felipe King, Mark H. Weir

# SET PARENT DIRECTORY & LOAD/INSTALL PACKAGES -----

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    

if("truncdist" %in% rownames(installed.packages())==FALSE){install.packages("truncdist"); require(truncdist)}else{require(truncdist)}
if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}
if("triangle" %in% rownames(installed.packages())==FALSE){install.packages("triangle"); require(triangle)}else{require(triangle)}
if("dplyr" %in% rownames(installed.packages())==FALSE){install.packages("dplyr"); require(dplyr)}else{require(dplyr)}
if("magrittr" %in% rownames(installed.packages())==FALSE){install.packages("magrittr"); require(magrittr)}else{require(magrittr)}

# IMPORT DATA NEEDED FOR ALL OPERATIONS DURING AUTOMATION -----

  #load in duration data
  durations<-read.csv('ContactDuration.csv')
  #read in bootstrapped values for dose-response
  exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')

# CONSTANTS TO BE USED 

SIM.iter <- 1000   # Making a master Monte Carlo iteration value
suppressMessages(source("adjust_behaviors_covid.R"))


# AUTOMATION ----
# Next line will be the names of the simulations that will be run. DCLP is Doffing correct and low patient numbers

#DP=doffing poor; DO=doffing ok; DG=doffing great
#LP=low prob COVID-19; HP=high prob COVID-19
#NS=non-surge (low # patient visits); S=surge (high # patient visits)
#control=1 patient visit, 100% prob of COVID-19 patient

SIM <- c("DPLPNS", "DPLPS","DPHPNS","DPHPS", "DPallNS", "DPallS",
         "DOLPNS","DOLPS","DOHPNS","DOHPS", "DOallNS", "DOallS",
         "DGLPNS","DGLPS","DGHPNS","DGHPS", "DGallNS", "DGallS",
         "control")   

NUM.SIM <- length(SIM)     # Count the number of iterations for the automated simulations

for(j in 1:NUM.SIM)
{
  sim.num <- j; sim.name <- SIM[j]
  
  if(sim.name=="DPLPNS"){if(dir.exists("DPLPNS")==FALSE){dir.create("DPLPNS"); setwd("DPLPNS")}else{setwd("DPLPNS")}}
  if(sim.name=="DPLPS"){if(dir.exists("DPLPS")==FALSE){dir.create("DPLPS"); setwd("DPLPS")}else{setwd("DPLPS")}}
  if(sim.name=="DPHPNS"){if(dir.exists("DPHPNS")==FALSE){dir.create("DPHPNS"); setwd("DPHPNS")}else{setwd("DPHPNS")}}
  if(sim.name=="DPHPS"){if(dir.exists("DPHPS")==FALSE){dir.create("DPHPS"); setwd("DPHPS")}else{setwd("DPHPS")}}
  if(sim.name=="DPallNS"){if(dir.exists("DPallNS")==FALSE){dir.create("DPallNS"); setwd("DPallNS")}else{setwd("DPallNS")}}
  if(sim.name=="DPallS"){if(dir.exists("DPallS")==FALSE){dir.create("DPallS"); setwd("DPallS")}else{setwd("DPallS")}}

  if(sim.name=="DOLPNS"){if(dir.exists("DOLPNS")==FALSE){dir.create("DOLPNS"); setwd("DOLPNS")}else{setwd("DOLPNS")}}
  if(sim.name=="DOLPS"){if(dir.exists("DOLPS")==FALSE){dir.create("DOLPS"); setwd("DOLPS")}else{setwd("DOLPS")}}
  if(sim.name=="DOHPNS"){if(dir.exists("DOHPNS")==FALSE){dir.create("DOHPNS"); setwd("DOHPNS")}else{setwd("DOHPNS")}}
  if(sim.name=="DOHPS"){if(dir.exists("DOHPS")==FALSE){dir.create("DOHPS"); setwd("DOHPS")}else{setwd("DOHPS")}}
  if(sim.name=="DOallNS"){if(dir.exists("DOallNS")==FALSE){dir.create("DOallNS"); setwd("DOallNS")}else{setwd("DOallNS")}}
  if(sim.name=="DOallS"){if(dir.exists("DOallS")==FALSE){dir.create("DOallS"); setwd("DOallS")}else{setwd("DOallS")}}
  
  if(sim.name=="DGLPNS"){if(dir.exists("DGLPNS")==FALSE){dir.create("DGLPNS"); setwd("DGLPNS")}else{setwd("DGLPNS")}}
  if(sim.name=="DGLPS"){if(dir.exists("DGLPS")==FALSE){dir.create("DGLPS"); setwd("DGLPS")}else{setwd("DGLPS")}}
  if(sim.name=="DGHPNS"){if(dir.exists("DGHPNS")==FALSE){dir.create("DGHPNS"); setwd("DGHPNS")}else{setwd("DGHPNS")}}
  if(sim.name=="DGHPS"){if(dir.exists("DGHPS")==FALSE){dir.create("DGHPS"); setwd("DGHPS")}else{setwd("DGHPS")}}
  if(sim.name=="DGallNS"){if(dir.exists("DGallNS")==FALSE){dir.create("DGallNS"); setwd("DGallNS")}else{setwd("DGallNS")}}
  if(sim.name=="DGallS"){if(dir.exists("DGallS")==FALSE){dir.create("DGallS"); setwd("DGallS")}else{setwd("DGallS")}}
  
  if(sim.name=="control"){if(dir.exists("control")==FALSE){dir.create("control"); setwd("control")}else{setwd("control")}}
  
  
# MASTER FUNCTION ----
behavior.sim<-function(caretype=c("IV","Obs","Rounds"),numsequence,prob.patient.infect,numvisit,prob.contam.between){
  
  #numvisit is function of shift length and number of patients, # of care episodes
  
  set.seed(34)
  #require(truncdist)
  
  #initializing vector to save final infection risk and dose per simulation
  finalinfectionrisk<-rep(NA,numsequence)
  finaldose<-rep(NA,numsequence)

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
      lambda[behavior!="Patient"]<-rtrunc(length(lambda[behavior!="Patient"]),spec="norm",a=0,b=1,mean=0.123,sd=0.068)#runif(length(lambda[behavior!="Patient"]),0.0003,.217)
      beta[behavior!="Patient"]<-rtrunc(length(lambda[behavior!="Patient"]),spec="norm",a=0,b=1,mean=0.118,sd=0.088)#runif(length(lambda[behavior!="Patient"]),0.0003,.217)
      #Mean and SD from: https://bmcinfectdis.biomedcentral.com/articles/10.1186/s12879-018-3425-x/tables/1
      #transfer efficiency patient skin contacts
      lambda[behavior=="Patient"]<-rtrunc(length(lambda[behavior=="Patient"]),spec="norm",a=0,b=1,mean=0.123,sd=0.068)#mean=.3,sd=.1)
      beta[behavior=="Patient"]<-rtrunc(length(lambda[behavior=="Patient"]),spec="norm",a=0,b=1,mean=0.118,sd=0.088)
      
      #-------------- SURFACE CONCENTRATIONS -----------------------------------------------------------------------------
      
      #Initialize surface concentration (PLACE HOLDER VALUE RIGHT NOW.. WILL DO LIT REVIEW TO INFORM THIS PARAMETER)
      surfconc<-rep(NA,length(behavior))
      
      #concentrations for hand-to-surface contact moments
      numtest<-runif(1,0,1)
      
      RNAinfectious<-rep(NA,length(behavior))
      
      
      if (numtest<=prob.patient.infect){
        
        #surface areas of swabbing (these are not provided by Guo et al., 2020, so we're estimating)
        
        surfacearea<-rtriangle(length(surfconc[behavior!="Patient"]),5,195,100) #cm^2
       
        #concentrations assuming patient is infected and asymptomatic
        
        #initialize vector for saving RNAinfectious values
        
        #Assume Out conc ~ In conc
        #min, max, and median of concentrations from Table 1 (up until pharmacy values) Guo et al. (2020)
        RNAinfectious[behavior!="Patient"]<-runif(length(surfconc[behavior!="Patient"]),0.001,0.1)
        surfconc[behavior!="Patient"]<-(rtriangle(length(surfconc[behavior!="Patient"]),a=3.3E3,b=6.6E4,c=2.8E4)/surfacearea)*RNAinfectious[behavior!="Patient"]
        
        #Patient surfaces
        RNAinfectious[behavior=="Patient"]<-runif(length(surfconc[behavior=="Patient"]),0.001,0.1)
        surfconc[behavior=="Patient"]<-(3.3*10^3/sample(surfacearea,1))*RNAinfectious[behavior=="Patient"] #point value, mask concentration Table 1 Guo et al. (2020)
          
      }else{
        
        #concentrations assuming patient is not infected
        
        #Assume Out conc ~ In conc
        surfconc[behavior!="Patient"]<-rep(0,length(surfconc[behavior!="Patient"]))
                 
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
      SH[behavior=="Equipment"|behavior=="FarPatient"|behavior=="NearPatient"|behavior=="HygieneInside"]<-runif(length(SH[behavior=="Equipment"|behavior=="FarPatient"|behavior=="NearPatient"|behavior=="HygieneInside"]),0.04/5,0.25) # change 06-07-2020 0.04/5
      #min and max of left and right hands in AuYeung et al. (2008) for various hand grip and hand press contacts (hand immersion contacts not included)
      #from single fingertip up to full palm
      
      #------------- RIGHT HAND VS. LEFT HAND ---------------------------------------------------------------------------
      
      #initialize R or L hand
      hand<-rep(NA,length(behavior))
      hand[sample(1:length(behavior),length(behavior)/2)]<-"right"
      hand[is.na(hand)]<-"left"
      
      #------------- INACTIVATION CONSTANTS -----------------------------------------------------------------------------
      
      #surface #This isn't t99 but t50  from Doremalen et al. 2020
      t50.s<-runif(length(behavior),4.59,8.17)*3600# 11-07-2020 change for widest CI runif(length(behavior),3,5*24)*3600 #T99 in seconds
      k.s<-log(2)/t50.s#(-log(1/(10^2))/t50.s) #corrected error with negative sign
      
      t99.h<-runif(length(behavior),1,24)*3600 #T99 in seconds
      k.h<-(-log(1/(10^2))/t99.h)
      
      #-------------- EXPOSURE SIMULATION ------------------------------------------------------------------------------
      
      #initialize conccentration on R and L hand estimates
      handR<-rep(0,length(behavior))
      handL<-rep(0,length(behavior))
      
      #placeholder for now in case we want to adjust duration
      duration<-sample(durations$Duration,length(behavior),replace=TRUE)
      timeSinceClean<-runif(length(behavior),0,24) #Since last cleaning
      
        if(hand[1]=="right"){
          handR[1]<-beta[1]*SH[1]*(surfconc[1]*exp(-k.s[1]*timeSinceClean[1]))#beta[1]*SH[1]*(surfconc[1]*exp(-k.s[1]*duration[1]))
          handL[1]<-0  
        }else{ #if hand[1] == left...
          handL[1]<-beta[1]*SH[1]*(surfconc[1]*exp(-k.s[1]*timeSinceClean[1]))
          handR[1]<-0 
        }

      for (a in 2:(length(behavior))){
        
        if(hand[a]=="right"){
          
            handR[a]<-(handR[a-1]-(lambda[a]*SH[a]*handR[a-1]*exp(-k.h[a]*duration[a]))+(beta[a]*SH[a]*surfconc[a]*exp(-k.s[a]*timeSinceClean[a])))#(handR[a-1]-(lambda[a]*SH[a]*handR[a-1]*exp(-k.h[a]*duration[a]))+(beta[a]*SH[a]*surfconc[a]*exp(-k.s[a]*duration[a])))
            handL[a]<-handL[a-1]

        }else{ #if hand[a] == left...
          
            handL[a]<-(handR[a-1]-(lambda[a]*SH[a]*handR[a-1]*exp(-k.h[a]*duration[a]))+(beta[a]*SH[a]*surfconc[a]*exp(-k.s[a]*timeSinceClean[a])))#(handL[a-1]-(lambda[a]*SH[a]*handL[a-1]*exp(-k.h[a]*duration[a]))+(beta[a]*SH[a]*surfconc[a]*exp(-k.s[a]*duration[a])))
            handR[a]<-handR[a-1]
        }
      
      
      }#end of care episode
      
      
      #---------------------------------- CALCULATION OF INFECTION RISK FOR CARE EPISODE --------------------------------------------------------------------------------------------------


      
      # require('gsl')
        
        #final concentration on hands a episode of care
      
        if (m==1){
          #if this is the first episode of care, handRbefore and handLbefore at this point are zero
          handRbefore<-0
          handLbefore<-0
        }

        #inactivation concentration on hands
        khtemp<- k.h[a]
        
        #doffing....
        
        #1/10 participants had <0.00003% of original inoculum transfered to hands during doffing
        
        #saving these var if no contact made
        TEtemp<-0
        SMtemp<-0
        AHtemp<-0
        
        randomnum<-runif(1,0,1)
        
        if(randomnum<=prob.contam.between){
          
          handRbefore<-(handR[a]*runif(1,3*10^-7,.1))+handRbefore  #06-07-2020 change (3*10^-7) #assume % of viral titer in study represents transfer potential when doffing
          handLbefore<-(handL[a]*runif(1,3*10^-7,.1))+handLbefore   #assume % of viral titer in study represents transfer potential when doffing
          
          #ethanol hand rub (described by Kampf et al. 2020 as recommended step after doffing)
          washnum<-runif(1,0,1)
          
          if (washnum<=.5){
            #then they wash their hands
            fracleft<-(1-(rtrunc(1,spec="norm",a=0,b=(1-10^(-4)),mean=abs(1-10^(-1.62)),sd=abs(1-10^0.12))))#based on change to assumed max reduction
          }else{
            fracleft<-runif(1,1E-4,1E-2) #6-07-2020 change (1/10^runif(1,2,4))  #based on Kampf et al. (2020)
          }
          reduce<-(1-fracleft)
          
          handRbefore<-handRbefore*fracleft
          handLbefore<-handLbefore*fracleft
          #Determines which hand will be used for the hand-to-face contact
          whichhand<-runif(1,0,1) #Nan's paper says that the non-dominant hand does most face touching - Could 
          
          #transfer efficiency to mouth
          TEtemp<-rtrunc(1,spec="norm",a=0,b=1,mean=.3390,sd=.1318) #mean: Rusin et al. (2002); sd: Lopez et al. (2012) 0.15. New Rusin data SD= 13.18%
          
          #fraction of hand used in mouth contact
          SMtemp<-runif(1,0.04/5, 0.06/5) #single fingertip surface area fraction = front partial fingers fractional surface area / 5
          
          #hand surface area
          AHtemp<-runif(1,445,535) #Beamer et al. (2015) office model and from Exposure Factors Handbook
          
          if (m==1){
            #if this is the first episode of care, this is our first dose
            if (whichhand<=.5){
              dose<-handLbefore*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handLbefore<-(1-TEtemp*SMtemp)*handLbefore*exp(-khtemp)
            }else{
              dose<-handRbefore*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handRbefore<-(1-TEtemp*SMtemp)*handRbefore*exp(-khtemp)
            }
          }else{
            #otherwise, we're adding onto our cumulative dose
            if (whichhand<=.5){
              dose<-dose+handLbefore*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handLbefore<-(1-TEtemp*SMtemp)*handLbefore*exp(-khtemp)
            }else{
              dose<-dose+handRbefore*TEtemp*SMtemp*AHtemp*exp(-khtemp) #assuming duration of 1 second
              handRbefore<-(1-TEtemp*SMtemp)*handRbefore*exp(-khtemp)
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
          handRbefore<-handRbefore #no contamination for next care episode
          handLbefore<-handLbefore #no contamination for next care episode
          
        }
        pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
        infecttemp<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose, give=FALSE, strict=TRUE)
        
       if(infecttemp==0){
         infecttemp<-1*10^-15 #cannot have zero infection risk, so replace with small risk
       }
        
      # end of potential self inoculation moment
      
      SMsave<-rep(SMtemp,length(behavior))
      AH.dose<-rep(AHtemp,length(behavior))
      TE.temp<-rep(TEtemp,length(behavior))
      
      doserep<-rep(dose,length(behavior))
      infectrep<-rep(infecttemp,length(behavior))
      alpha<-rep(exactbp$alpha[pair],length(behavior))
      beta.doseresp<-rep(exactbp$Beta[pair],length(behavior))
      reducesave<-rep(reduce,length(behavior))
      
      #saving concentrations for all rooms
      if (m==1){
        reduceall<-reducesave
        SMall<-SMsave
        AH.doseall<-AH.dose
        TE.tempall<-TE.temp
        alphaall<-alpha
        beta.dose.all<-beta.doseresp
        handRnoglove<-rep(handRbefore,length(behavior))
        handLnoglove<-rep(handLbefore,length(behavior))
        RNAinfectiousall<-RNAinfectious
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
        reduceall<-c(reduceall,reducesave)
        SMall<-c(SMall,SMsave)
        AH.doseall<-c(AH.doseall,AH.dose)
        TE.tempall<-c(TE.tempall,TE.temp)
        alphaall<-c(alphaall,alpha)
        beta.dose.all<-c(beta.dose.all,beta.doseresp)
        RNAinfectiousall<-c(RNAinfectiousall,RNAinfectious)
        handRnoglovetemp<-rep(handRbefore,length(behavior))
        handLnoglovetemp<-rep(handLbefore,length(behavior))
        handRnoglove<-c(handRnoglove,handRnoglovetemp)
        handLnoglove<-c(handLnoglove,handLnoglovetemp)
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
      
      
    } #end of visit m
  
    # -------------------------------- SAVE OUTPUT FOR SIMULATION FOR SINGLE PERSON ----------------------------------------------------------------------------------
    exposure.frame.temp<-data.frame(dose=doseall,infect=infectall,patientnum=patientnum,handR=handRall,handL=handLall,hand=handall,handRnoglove=handRnoglove,handLnoglove=handLnoglove,
                                    behavior=behaviorall,duration=durationall,SH=SHall,lambda=lambdaall,beta=betaall,surfconc=surfconcall,k.sall=k.sall,k.hall=k.hall,
                                    RNAinfectiousall=RNAinfectiousall,alphaall=alphaall,beta.dose.all=beta.dose.all,SMall=SMall,
                                    AH.doseall=AH.doseall,TE.tempall=TE.tempall,reduceall=reduceall)
    
    behavior.total[[j]]<-behavior
    exposure.frame[[j]]<-exposure.frame.temp
    finalinfectionrisk[j]<-max(infectall)
    finaldose[j]<-max(doseall)

    rm(behavior)
    
  }#end of iteration loop
  # --------------------------------- SAVE ALL OUTPUTS TO GLOBAL ENV --------------------------------------------------------------------------------------------------
  behavior.total<<-behavior.total
  exposure.frame<<-exposure.frame
  finalinfectionrisk<<-finalinfectionrisk
  finaldose<<-finaldose
  
  
} # Master function end


#dawn poor, dawn ok, dawn good
probcontam<-c(0.8,0.5,0.1)

#low patient load, high patient load
patload<-c(14,7) #7 based on estimates from Ian; consider double that for high load situations

#non-surge, high surge, or all COVID ward (probability of 1 for any given patient being COVID-19 positive)
contam<-c(0.05,0.5,1)

if (j<7 | j==19){
  prob.contam.between=probcontam[1]
}else if (j>6 & j<13){
  prob.contam.between=probcontam[2]
}else{
  prob.contam.between=probcontam[3]
}


if(j==1 | j==2 | j==7 | j==8 | j==13 | j==14){
  prob.patient.infect=contam[1]
}else if (j==3 | j==4 | j==9 | j==10 | j==15 | j==16){
  prob.patient.infect=contam[2]
}else{
  prob.patient.infect=contam[3]
}

if (j %% 2 ==0 & j!=19){
  numvisit=patload[1]
}else if (j!=19){
  numvisit=patload[2]
}else{
  numvisit=1
}

#so "control" scenario is probcontam (0.8), 1 probability of COVID watch patient,
#single patient visit

#-------------------------------------------------------

#IV scenario
behavior.sim(caretype="IV",numsequence=SIM.iter,prob.patient.infect=prob.patient.infect,
             numvisit=numvisit,prob.contam.between=prob.contam.between)
IV<-exposure.frame
write.csv(finalinfectionrisk,file=sprintf("%s.IV.risks.cvs",sim.name))
IV.risk<-finalinfectionrisk
IV.dose<-finaldose
      
#Rounds scenario
behavior.sim(caretype="Rounds",numsequence=SIM.iter,prob.patient.infect=prob.patient.infect,
             numvisit=numvisit,prob.contam.between=prob.contam.between)
Rounds<-exposure.frame
write.csv(finalinfectionrisk,file=sprintf("%s.Rounds.risks.cvs",sim.name))
Rounds.risk<-finalinfectionrisk
Rounds.dose<-finaldose
      
#Observational care scenario
behavior.sim(caretype="Obs",numsequence=SIM.iter,prob.patient.infect=prob.patient.infect,
             numvisit=numvisit,prob.contam.between=prob.contam.between)
Obs<-exposure.frame
Obs.risk<-finalinfectionrisk
Obs.dose<-finaldose
write.csv(finalinfectionrisk,file=sprintf("%s.Obs.risks.cvs",sim.name))
      
#save output
saveRDS(IV, file=sprintf("%s.IV.exposure.frame.rds",sim.name))
saveRDS(Rounds, file=sprintf("%s.Rounds.exposure.frame.rds",sim.name))
saveRDS(Obs, file=sprintf("%s.Obs.exposure.frame.rds",sim.name))

  if (j==1){
    caretype<-c(rep("IV",SIM.iter),rep("Rounds",SIM.iter),rep("Observations",SIM.iter))
    risk=c(IV.risk,Rounds.risk,Obs.risk)
    dose=c(IV.dose,Rounds.dose,Obs.dose)
    prob.contam.between.all=rep(prob.contam.between,length(risk))
    prob.patient.infect.all=rep(prob.patient.infect,length(risk))
    numvisit.all=rep(numvisit,length(risk))
  }else{
    caretypetemp<-c(rep("IV",SIM.iter),rep("Rounds",SIM.iter),rep("Observations",SIM.iter))
    risktemp=c(IV.risk,Rounds.risk,Obs.risk)
    dosetemp=c(IV.dose,Rounds.dose,Obs.dose)
    dose<-c(dose,dosetemp)
    risk<-c(risk,risktemp)
    caretype<-c(caretype,caretypetemp)
    prob.contam.between.temp=rep(prob.contam.between,length(risktemp))
    prob.patient.infect.temp=rep(prob.patient.infect,length(risktemp))
    numvisit.temp=rep(numvisit,length(risktemp))
    
    prob.contam.between.all<-c(prob.contam.between.all,prob.contam.between.temp)
    prob.patient.infect.all<-c(prob.patient.infect.all,prob.patient.infect.temp)
    numvisit.all<-c(numvisit.all,numvisit.temp)
  }

  #reset directory to parent folder so we can go to correct subfolder within parent folder for next sim run
  setwd(this.dir)
  
} # Automation loop end


frameall<-data.frame(risk=risk,probcontambetween=as.character(prob.contam.between.all),
                     probpatientinfect=as.character(prob.patient.infect.all),
                     numvisit=as.character(numvisit.all),caretype=caretype,dose=dose)

#summary of risks (all care types, all scenarios) referenced in results section
summary(risk*100)
summary(frameall$risk[frameall$caretype=="IV"]*100)
summary(frameall$risk[frameall$caretype=="Rounds"]*100)
summary(frameall$risk[frameall$caretype=="Observations"]*100)
aov.results<-aov(risk ~ caretype,data=frameall)
summary(aov.results)

summary(frameall$risk[frameall$numvisit==7]*100)
summary(frameall$risk[frameall$numvisit==14]*100)

summary(frameall$risk[frameall$probpatientinfect==1]*100)
summary(frameall$risk[frameall$probpatientinfect==0.5]*100)
summary(frameall$risk[frameall$probpatientinfect==0.05]*100)
sd(frameall$risk[frameall$probpatientinfect==0.05]*100)

require(ggplot2)
require(ggpubr)

# Stair plot
df<-behavior.total[[1]] %>%as.data.frame()%>%set_colnames("df")
ggplot(df) +
  geom_step(aes(x=seq_along(df), y=df,group=1)) +
  geom_point(aes(x=seq_along(df), y=df),color="red")+
  xlab("Surface contact number") +
  ylab("Surface type") + theme_pubr()

# Risk plots
# A<-ggplot(frameall[frameall$numvisit!=1,])+geom_boxplot(aes(x=probcontambetween,
#                                     fill=numvisit,y=risk))+
#   scale_y_continuous(trans="log10",name="Infection Risk")+
#   scale_x_discrete(name="Probability of Contamination Between Care Episodes")+
#   scale_fill_brewer(palette = "Paired")+#scale_fill_discrete(name="Number of Patient Visits")+
#   facet_grid(numvisit~probpatientinfect,scales="free")+
#   theme_pubr()
# 
# 
# ggplot(frameall[frameall$numvisit!=1,])+
#   geom_density(aes(fill=numvisit,x=risk),width=0.3)+
#   scale_x_continuous(trans="log10",name="Infection Risk")+
#   #scale_x_discrete(name="Probability of Contamination Between Care Episodes")+
#   scale_fill_brewer(palette = "Set1",name="Number of Patient Visits")+#scale_fill_discrete(name="Number of Patient Visits")+
#   facet_grid(numvisit~probpatientinfect,scales="free")+
#   theme_pubr()

require(data.table)
df<-setDT(frameall)[,.(Mean = mean(risk*100), SD = sd(risk*100), Min = min(risk*100), Max = max(risk*100)), by = c('caretype','probcontambetween', 'probpatientinfect', 'numvisit')]
df$numvisit <- factor(df$numvisit, levels = c(1,7,14))

# stat.test <- frameall %>%
#   group_by('caretype','probcontambetween', 'probpatientinfect', 'numvisit') %>%
#   t_test(len ~ risk, ref.group = "0.5") 

#Barplot with standard deviations
windows()
ggplot(subset(df)) +
  geom_bar( aes(x=probcontambetween,y=Mean,fill=caretype), stat="identity",  position='dodge') + #, fill="skyblue"
  geom_pointrange( aes(x=probcontambetween, y=Mean, ymin=pmax(Mean - SD, 0), ymax=Mean+SD, position=caretype), colour="black", alpha=0.9, size=0.2, position = position_dodge(width =1))+
  facet_grid(numvisit~probpatientinfect,scales="free",
             labeller = labeller(
               numvisit = c(`1` = "1 care episode",`7` = "7 care episodes", `14` = "14 care episodes"),
               probpatientinfect = c(`0.05` = "5% COVID load", `0.5` = "50% COVID load",  `1` = "100% COVID load")
             ))+
  scale_fill_brewer(palette = "Set1",name="Care Type")+
  scale_x_discrete(name="Probability of self contamination (out of 1)")+
  scale_y_continuous(name="Infection risk %")+
  theme_pubr()+theme(axis.text=element_text(size=18),axis.title=element_text(size=18),
                     strip.text = element_text(size=15),legend.text = element_text(size=18),
                     legend.title = element_text(size=18))

#windowss() https://github.com/eliocamp/ggnewscale
#A

# A.2<-ggplot(frameall[frameall$numvisit!=1 & frameall$risk>1e-15,])+geom_boxplot(aes(x=probcontambetween,
#                                                             fill=caretype,y=risk))+
#   scale_y_continuous(trans="log10",name="Infection Risk")+
#   scale_x_discrete(name="Probability of Contamination Between Care Episodes")+
#   scale_fill_grey(name="Care Type",start = 0.4,end=0.8)+
#   facet_grid(numvisit~probpatientinfect,scales="free")+
#   theme_pubr()+ggtitle("A")
# 
# #windowss()
# A.2
# 
# A.3<-ggplot(frameall[frameall$numvisit!=1,])+geom_boxplot(aes(x=probcontambetween,
#                                                                                     fill=caretype,y=risk))+
#   scale_y_continuous(trans="log10",name="Infection Risk")+
#   scale_x_discrete(name="Probability of Contamination Between Care Episodes")+
#   scale_fill_grey(name="Care Type",start = 0.4,end=0.8)+
#   facet_grid(numvisit~probpatientinfect,scales="free")+
#   theme_pubr()+ggtitle("B")
# 
# #windowss()
# A.3
# 
# ggarrange(A.2,A.3,nrow=2,ncol=1,common.legend = TRUE)
# 
# 
# B<-ggplot(frameall[frameall$numvisit!=1,])+geom_boxplot(aes(x=probcontambetween,
#                                                             fill=numvisit,y=dose))+
#   scale_y_continuous(trans="log10",name="Dose")+
#   scale_x_discrete(name="Probability of Contamination Between Care Episodes")+
#   scale_fill_discrete(name="Number of Patient Visits")+
#   facet_grid(numvisit~probpatientinfect,scales="free")+
#   theme_pubr()
# 
# #windowss()
# B



#----------- summary statistics

#summary.stats<-function(probcontambetween,numvisit,probpatientinfect){
  
#  print(signif(summary(frameall$risk[frameall$probcontambetween==probcontambetween & 
#                                 frameall$numvisit==numvisit &
#                                 frameall$probpatientinfect==probpatientinfect]),2))
#  print(signif(min(frameall$risk[frameall$probcontambetween==probcontambetween & 
#                                       frameall$numvisit==numvisit &
#                                       frameall$probpatientinfect==probpatientinfect]),2))
  
  
  
  
#  print(signif(sd(frameall$risk[frameall$probcontambetween==probcontambetween & 
#                            frameall$numvisit==numvisit &
#                            frameall$probpatientinfect==probpatientinfect]),2))
  
#}

#summary.stats(probcontambetween=0.8,numvisit=1,probpatientinfect=1)


#summary.stats(probcontambetween=0.1,numvisit=7,probpatientinfect=0.05)
#summary.stats(probcontambetween=0.1,numvisit=7,probpatientinfect=0.5)
#summary.stats(probcontambetween=0.1,numvisit=7,probpatientinfect=1)

#summary.stats(probcontambetween=0.5,numvisit=7,probpatientinfect=0.05)
#summary.stats(probcontambetween=0.5,numvisit=7,probpatientinfect=0.5)
#summary.stats(probcontambetween=0.5,numvisit=7,probpatientinfect=1)

#summary.stats(probcontambetween=0.8,numvisit=7,probpatientinfect=0.05)
#summary.stats(probcontambetween=0.8,numvisit=7,probpatientinfect=0.5)
#summary.stats(probcontambetween=0.8,numvisit=7,probpatientinfect=1)


#summary.stats(probcontambetween=0.1,numvisit=14,probpatientinfect=0.05)
#summary.stats(probcontambetween=0.1,numvisit=14,probpatientinfect=0.5)
#summary.stats(probcontambetween=0.1,numvisit=14,probpatientinfect=1)

#summary.stats(probcontambetween=0.5,numvisit=14,probpatientinfect=0.05)
#summary.stats(probcontambetween=0.5,numvisit=14,probpatientinfect=0.5)
#summary.stats(probcontambetween=0.5,numvisit=14,probpatientinfect=1)

#summary.stats(probcontambetween=0.8,numvisit=14,probpatientinfect=0.05)
#summary.stats(probcontambetween=0.8,numvisit=14,probpatientinfect=0.5)
#summary.stats(probcontambetween=0.8,numvisit=14,probpatientinfect=1)


#----------- summary statistics per caretype

#summary.stats.caretype<-function(probcontambetween,numvisit,probpatientinfect,caretype=c("IV")){
  
#  print(signif(summary(frameall$risk[frameall$probcontambetween==probcontambetween & 
#                                       frameall$numvisit==numvisit &
#                                       frameall$probpatientinfect==probpatientinfect & frameall$caretype==caretype]),2))
#  print(signif(min(frameall$risk[frameall$probcontambetween==probcontambetween & 
#                                   frameall$numvisit==numvisit &
#                                   frameall$probpatientinfect==probpatientinfect & frameall$caretype==caretype]),2))
  
  
  
  
#  print(signif(sd(frameall$risk[frameall$probcontambetween==probcontambetween & 
#                                  frameall$numvisit==numvisit &
#                                  frameall$probpatientinfect==probpatientinfect & frameall$caretype==caretype]),2))
  
#}

#caretype="IV"

#summary.stats.caretype(probcontambetween=0.8,numvisit=1,probpatientinfect=1,caretype=caretype)

#summary.stats.caretype(probcontambetween=0.1,numvisit=7,probpatientinfect=0.05,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.1,numvisit=7,probpatientinfect=0.5,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.1,numvisit=7,probpatientinfect=1,caretype=caretype)

#summary.stats.caretype(probcontambetween=0.5,numvisit=7,probpatientinfect=0.05,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.5,numvisit=7,probpatientinfect=0.5,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.5,numvisit=7,probpatientinfect=1,caretype=caretype)

#summary.stats.caretype(probcontambetween=0.8,numvisit=7,probpatientinfect=0.05,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.8,numvisit=7,probpatientinfect=0.5,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.8,numvisit=7,probpatientinfect=1,caretype=caretype)


#summary.stats.caretype(probcontambetween=0.1,numvisit=14,probpatientinfect=0.05,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.1,numvisit=14,probpatientinfect=0.5,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.1,numvisit=14,probpatientinfect=1,caretype=caretype)

#summary.stats.caretype(probcontambetween=0.5,numvisit=14,probpatientinfect=0.05,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.5,numvisit=14,probpatientinfect=0.5,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.5,numvisit=14,probpatientinfect=1,caretype=caretype)

#summary.stats.caretype(probcontambetween=0.8,numvisit=14,probpatientinfect=0.05,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.8,numvisit=14,probpatientinfect=0.5,caretype=caretype)
#summary.stats.caretype(probcontambetween=0.8,numvisit=14,probpatientinfect=1,caretype=caretype)



#---------------- dose summary stats


#summary.stats.dose<-function(probcontambetween,numvisit,probpatientinfect){
  
#  print(signif(summary(frameall$dose[frameall$probcontambetween==probcontambetween & 
#                                       frameall$numvisit==numvisit &
#                                       frameall$probpatientinfect==probpatientinfect]),2))
#  print(signif(min(frameall$dose[frameall$probcontambetween==probcontambetween & 
#                                   frameall$numvisit==numvisit &
#                                   frameall$probpatientinfect==probpatientinfect]),2))
  
  
  
  
#  print(signif(sd(frameall$dose[frameall$probcontambetween==probcontambetween & 
#                                  frameall$numvisit==numvisit &
#                                  frameall$probpatientinfect==probpatientinfect]),2))
  
#}


#summary.stats.dose(probcontambetween=0.8,numvisit=1,probpatientinfect=1)

#summary.stats.dose(probcontambetween=0.1,numvisit=7,probpatientinfect=0.05)
#summary.stats.dose(probcontambetween=0.1,numvisit=7,probpatientinfect=0.5)
#summary.stats.dose(probcontambetween=0.1,numvisit=7,probpatientinfect=1)

#summary.stats.dose(probcontambetween=0.5,numvisit=7,probpatientinfect=0.05)
#summary.stats.dose(probcontambetween=0.5,numvisit=7,probpatientinfect=0.5)
#summary.stats.dose(probcontambetween=0.5,numvisit=7,probpatientinfect=1)

#summary.stats.dose(probcontambetween=0.8,numvisit=7,probpatientinfect=0.05)
#summary.stats.dose(probcontambetween=0.8,numvisit=7,probpatientinfect=0.5)
#summary.stats.dose(probcontambetween=0.8,numvisit=7,probpatientinfect=1)


#summary.stats.dose(probcontambetween=0.1,numvisit=14,probpatientinfect=0.05)
#summary.stats.dose(probcontambetween=0.1,numvisit=14,probpatientinfect=0.5)
#summary.stats.dose(probcontambetween=0.1,numvisit=14,probpatientinfect=1)

#summary.stats.dose(probcontambetween=0.5,numvisit=14,probpatientinfect=0.05)
#summary.stats.dose(probcontambetween=0.5,numvisit=14,probpatientinfect=0.5)
#summary.stats.dose(probcontambetween=0.5,numvisit=14,probpatientinfect=1)

#summary.stats.dose(probcontambetween=0.8,numvisit=14,probpatientinfect=0.05)
#summary.stats.dose(probcontambetween=0.8,numvisit=14,probpatientinfect=0.5)
#summary.stats.dose(probcontambetween=0.8,numvisit=14,probpatientinfect=1)


#---------------- dose summary stats, care type specific


#summary.stats.dose.caretype<-function(probcontambetween,numvisit,probpatientinfect,caretype=c("IV")){
  
#  print(signif(summary(frameall$dose[frameall$probcontambetween==probcontambetween & 
#                                       frameall$numvisit==numvisit &
#                                       frameall$probpatientinfect==probpatientinfect & frameall$caretype==caretype]),2))
#  print(signif(min(frameall$dose[frameall$probcontambetween==probcontambetween & 
#                                   frameall$numvisit==numvisit &
#                                   frameall$probpatientinfect==probpatientinfect & frameall$caretype==caretype]),2))
  
  
  
  
#  print(signif(sd(frameall$dose[frameall$probcontambetween==probcontambetween & 
#                                  frameall$numvisit==numvisit &
#                                  frameall$probpatientinfect==probpatientinfect & frameall$caretype==caretype]),2))
  
#}

#caretype="IV"

#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=1,probpatientinfect=1,caretype=caretype)

#summary.stats.dose.caretype(probcontambetween=0.1,numvisit=7,probpatientinfect=0.05,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.1,numvisit=7,probpatientinfect=0.5,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.1,numvisit=7,probpatientinfect=1,caretype=caretype)

#summary.stats.dose.caretype(probcontambetween=0.5,numvisit=7,probpatientinfect=0.05,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.5,numvisit=7,probpatientinfect=0.5,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.5,numvisit=7,probpatientinfect=1,caretype=caretype)

#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=7,probpatientinfect=0.05,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=7,probpatientinfect=0.5,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=7,probpatientinfect=1,caretype=caretype)


#summary.stats.dose.caretype(probcontambetween=0.1,numvisit=14,probpatientinfect=0.05,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.1,numvisit=14,probpatientinfect=0.5,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.1,numvisit=14,probpatientinfect=1,caretype=caretype)

#summary.stats.dose.caretype(probcontambetween=0.5,numvisit=14,probpatientinfect=0.05,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.5,numvisit=14,probpatientinfect=0.5,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.5,numvisit=14,probpatientinfect=1,caretype=caretype)

#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=14,probpatientinfect=0.05,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=14,probpatientinfect=0.5,caretype=caretype)
#summary.stats.dose.caretype(probcontambetween=0.8,numvisit=14,probpatientinfect=1,caretype=caretype)




