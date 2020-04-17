# Amanda Wilson, Marco-Felipe King, Mark H. Weir



behavior.sim<-function(room.orientation=c("left","right"),caretype=c("IV","Obs","Rounds"),numsequence,reducebetween=c("glove change", "ethanol","no change"),
                       prob.patient.infect){
   set.seed(34)
   require(truncdist)
  
  #load in duration data
  durations<-read.csv('ContactDuration.csv')
  #Remove "S" at end of entries and convert from character to numeric format
  durations$Duration<-as.numeric(gsub("S","",durations$Duration))
 
  #Behavior model section

  sample.space<-c("Equipment","FarPatient","HygieneInside","In","NearPatient","Out","Patient")
  
  behavior.total<-list() #creating a list to store behaviors
  exposure.frame<-list() #creating a list to store exposure.frame data.frames
  
  #set up matrix for type of care
  if (room.orientation=="left"){
    if (caretype=="IV"){
      prob.mat<-TIV.left$estimate #left-facing, IV
    }else if (caretype=="Obs"){
      prob.mat<-TObs.left$estimate #left-facing, Obs
    }else{
      prob.mat<-TRounds.left$estimate #left-facing, Rounds
    }
  }else{
    if (caretype=="IV"){
      prob.mat<-TIV.right$estimate #right-facing IV
    }else if (caretype=="Obs"){
      prob.mat<-TObs.right$estimate #right-facing Obs
    }else{
      prob.mat<-TRounds.right$estimate #right-facing Rounds
    }
  }
  

  for (j in 1:numsequence){
    
    for (m in 1:4){
      
      #determine reduction between patients...
      if (reducebetween=="glove change"){
        reducebtwneff<-0 #final conc from previous room will be multiplied by zero - no contamination
      }else if (reducebetween=="ethanol"){
        reducebtwneff<-(1/10^runif(1,2,4)) #divide final conc by efficacy (Kampf et al., 2020 min and max efficacies from carriage tests for ethanol)
      }else{
        reducebtwneff<-1 #final conc from previous room will be multiplied by 1 - no change in contamination
      }
      
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
      
      #fractional surface area for open hand grip
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
          handR[1]<-0-beta[1]*SH[1]*(0-surfconc[1]*exp(-k.s[1]*duration[1]))
          handL[1]<-0  
        }else{ #if hand[1] == left...
          handL[1]<-0-beta[1]*SH[1]*(0-surfconc[1]*exp(-k.s[1]*duration[1]))
          handR[1]<-0 
        }
      }else{
        handR[1]<-handRnext*reducebtwneff
        handL[1]<-handLnext*reducebtwneff
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
      
      
      #saving final concentration hands to serve as initial concentration on hands for next room
      handRnext<-handR[a]
      handLnext<-handL[a]
      
      #saving concentrations for all 4 rooms
      if (m==1){
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
    exposure.frame.temp<-data.frame(patientnum=patientnum,handR=handRall,handL=handLall,hand=handall,behavior=behaviorall,duration=durationall,SH=SHall,lambda=lambdaall,beta=betaall,surfconc=surfconcall,k.sall=k.sall,k.hall=k.hall)
    behavior.total[[j]]<-behavior
    exposure.frame[[j]]<-exposure.frame.temp

    rm(behavior)
  }
  # --------------------------------- SAVE ALL OUTPUTS TO GLOBAL ENV --------------------------------------------------------------------------------------------------
  behavior.total<<-behavior.total
  exposure.frame<<-exposure.frame
  
  #---------------------------------- CALCULATION OF INFECTION RISK --------------------------------------------------------------------------------------------------
  
  dose<-rep(NA,numsequence)
  infect<-rep(NA,numsequence)
  
  #read in bootstrapped values for dose-response
  exactbp<-read.csv('Exact_BetaPoisson_Bootstrap.csv')
  
  require('gsl')
  
  for (b in 1:numsequence){
    
    #final concentration on hands after 4 episodes of cre
    righthand<- exposure.frame[[b]]$handR[length(exposure.frame[[b]]$patientnum)]
    lefthand<- exposure.frame[[b]]$handL[length(exposure.frame[[b]]$patientnum)]
    
    #inactivation concentration on hands
    kh<- exposure.frame[[b]]$k.hall[length(exposure.frame[[b]]$patientnum)]
    
    #doffing....
    
    #1/10 participants had <0.00003% of original inoculum transfered to hands during doffing
    
    randomnum<-runif(1,0,1)
    
    if(randomnum<=0.1){
      
    righthand<-righthand*(3*10^-7) #assume % of viral titer in study represents transfer potential when doffing
    lefthand<-lefthand*(3*10^-7) #assume % of viral titer in study represents transfer potential when doffing
    
    #ethanol hand rub (described by Kampf et al. 2020 as recommended step after doffing)
    washnum<-runif(1,0,1)
    
    if (washnum<=.5){
      #then they wash their hands
      righthand<-righthand*(1/10^rtrunc(1,spec="norm",a=0,b=6,mean=1.62,sd=0.12)) #based on dist from Marco
      lefthand<-lefthand*(1/10^rtrunc(1,spec="norm",a=0,b=6,mean=1.62,sd=0.12)) #based on dist from Marco
    }else{
      #otherwise they use hand sanitizer
      righthand<-righthand*(1/10^runif(1,2,4)) #based on Kampf et al. (2020)
      lefthand<-lefthand*(1/10^runif(1,2,4)) #based on Kampf et al. (2020)
    }
    
    #Determines which hand will be used for the hand-to-face contact
    whichhand<-runif(1,0,1)
    
    #transfer efficiency to mouth
    TE<-.3390 #Rusin et al. (2002)
      
    #fraction of hand used in mouth contact
    SM<-runif(1,0.04/5, 0.06/5) #single fingertip surface area fraction = front partial fingers fractional surface area / 5
      
    #hand surface area
    AH<-runif(1,445,535) #Beamer et al. (2015) office model and from Exposure Factors Handbook
    
    if (whichhand<=.5){
      dose[b]<-lefthand*TE*SM*AH*exp(-kh) #assuming duration of 1 second
    }else{
      dose[b]<-righthand*TE*SM*AH*exp(-kh) #assuming duration of 1 second
    }
      
    }else{
      
    #if no contam while doffing, then handwashing wouldn't change concentration
    #and therefore would not change dose
    dose[b]<-0
    dose[b]<-lefthand<-0
    
    }
    pair<-sample(c(1:length(exactbp$ln.alpha.)),1)
    infect[b]<-1-hyperg_1F1(exactbp$alpha[pair], exactbp$alpha[pair]+exactbp$Beta[pair], -dose[b], give=FALSE, strict=TRUE)
    
      
  }
  
  dose<<-dose
  infect<<-infect
  
  #print("Order up! :)")
  
}

#------------------------------------------ IV care---------------------------------------------------------------------

behavior.sim(room.orientation="left",caretype="IV",numsequence=500,reducebetween="no change",prob.patient.infect=.25)
infect.left.IV.nochange<-infect
behavior.sim(room.orientation="right",caretype="IV",numsequence=500,reducebetween="no change",prob.patient.infect=.25)
infect.right.IV.nochange<-infect

behavior.sim(room.orientation="left",caretype="IV",numsequence=500,reducebetween="glove change",prob.patient.infect=.25)
infect.left.IV.glovechange<-infect
behavior.sim(room.orientation="right",caretype="IV",numsequence=500,reducebetween="glove change",prob.patient.infect=.25)
infect.right.IV.glovechange<-infect

behavior.sim(room.orientation="left",caretype="IV",numsequence=500,reducebetween="ethanol",prob.patient.infect=.25)
infect.left.IV.ethanol<-infect
behavior.sim(room.orientation="right",caretype="IV",numsequence=500,reducebetween="ethanol",prob.patient.infect=.25)
infect.right.IV.ethanol<-infect


#---------------------------------------Rounds Care ------------------------------------------------------------------

behavior.sim(room.orientation="left",caretype="Rounds",numsequence=500,reducebetween="no change",prob.patient.infect=.25)
infect.left.rounds.nochange<-infect
behavior.sim(room.orientation="right",caretype="Rounds",numsequence=500,reducebetween="no change",prob.patient.infect=.25)
infect.right.rounds.nochange<-infect

behavior.sim(room.orientation="left",caretype="Rounds",numsequence=500,reducebetween="glove change",prob.patient.infect=.25)
infect.left.rounds.glovechange<-infect
behavior.sim(room.orientation="right",caretype="Rounds",numsequence=500,reducebetween="glove change",prob.patient.infect=.25)
infect.right.rounds.glovechange<-infect

behavior.sim(room.orientation="left",caretype="Rounds",numsequence=500,reducebetween="ethanol",prob.patient.infect=.25)
infect.left.rounds.ethanol<-infect
behavior.sim(room.orientation="right",caretype="Rounds",numsequence=500,reducebetween="ethanol",prob.patient.infect=.25)
infect.right.rounds.ethanol<-infect

#---------------------------------------Observational care ----------------------------------------------------------------

behavior.sim(room.orientation="left",caretype="Obs",numsequence=500,reducebetween="no change",prob.patient.infect=.25)
infect.left.obs.nochange<-infect
behavior.sim(room.orientation="right",caretype="Obs",numsequence=500,reducebetween="no change",prob.patient.infect=.25)
infect.right.obs.nochange<-infect

behavior.sim(room.orientation="left",caretype="Obs",numsequence=500,reducebetween="glove change",prob.patient.infect=.25)
infect.left.obs.glovechange<-infect
behavior.sim(room.orientation="right",caretype="Obs",numsequence=500,reducebetween="glove change",prob.patient.infect=.25)
infect.right.obs.glovechange<-infect

behavior.sim(room.orientation="left",caretype="Obs",numsequence=500,reducebetween="ethanol",prob.patient.infect=.25)
infect.left.obs.ethanol<-infect
behavior.sim(room.orientation="right",caretype="Obs",numsequence=500,reducebetween="ethanol",prob.patient.infect=.25)
infect.right.obs.ethanol<-infect

#--------------------------------------------------------------------------------------------------------------------------

infectall<-c(infect.left.IV.nochange,infect.right.IV.nochange,
             infect.left.IV.ethanol,infect.right.IV.ethanol,
             infect.left.IV.glovechange,infect.right.IV.glovechange,
             
             infect.left.rounds.nochange,infect.right.rounds.nochange,
             infect.left.rounds.ethanol,infect.right.rounds.ethanol,
             infect.left.rounds.glovechange,infect.right.rounds.glovechange,
             
             infect.left.obs.nochange,infect.right.obs.nochange,
             infect.left.obs.ethanol,infect.right.obs.ethanol,
             infect.left.obs.glovechange,infect.right.obs.glovechange)

caretype<-c(rep("IV",length(c(infect.left.IV.nochange,infect.right.IV.nochange,
                              infect.left.IV.ethanol,infect.right.IV.ethanol,
                              infect.left.IV.glovechange,infect.right.IV.glovechange))),
            
            rep("Rounds",length(c(infect.left.rounds.nochange,infect.right.rounds.nochange,
                                  infect.left.rounds.ethanol,infect.right.rounds.ethanol,
                                  infect.left.rounds.glovechange,infect.right.rounds.glovechange))),
            
            rep("Obs",length(c(infect.left.obs.nochange,infect.right.obs.nochange,
                               infect.left.obs.ethanol,infect.right.obs.ethanol,
                               infect.left.obs.glovechange,infect.right.obs.glovechange))))

intervention<-rep(c(rep("No Change",1000),rep("Ethanol",1000),rep("Glove Change",1000)),3)

frameall<-data.frame(infect=infectall,care=caretype,intervention=intervention)

print("Order up! :)")