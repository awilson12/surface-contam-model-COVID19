# sensitivity analysis


this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    


SIM <- c("DPLPNS", "DPLPS","DPHPNS","DPHPS",
         "DOLPNS","DOLPS","DOHPNS","DOHPS",
         "DGLPNS","DGLPS","DGHPNS","DGHPS","control")   

NUM.SIM <- length(SIM)     # Count the number of iterations for the automated simulations

for(j in 1:NUM.SIM){
  
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

  
  numiterations<-500 
  
  IV<-readRDS(file=sprintf("%s.IV.exposure.frame.rds",sim.name))
  Rounds<-readRDS(file=sprintf("%s.Rounds.exposure.frame.rds",sim.name))
  Obs<-readRDS(file=sprintf("%s.Obs.exposure.frame.rds",sim.name))
  

  for(i in 1:numiterations){
    IVframe<-IV[[i]]
    Roundsframe<-Rounds[[i]]
    Obsframe<-Obs[[i]]
   
    if (i ==1){
      lambda<-c(mean(IVframe$lambda),mean(Roundsframe$lambda),mean(Obsframe$lambda))
      beta<-c(mean(IVframe$beta),mean(Roundsframe$beta),mean(Obsframe$beta))
      duration<-c(mean(IVframe$duration),mean(Roundsframe$duration),mean(Obsframe$duration))
      SH<-c(mean(IVframe$SH),mean(Roundsframe$SH),mean(Obsframe$SH))
      surfconc<-c(mean(IVframe$surfconc),mean(Roundsframe$surfconc),mean(Obsframe$surfconc))
      k.sall<-c(mean(IVframe$k.sall),mean(Roundsframe$k.sall),mean(Obsframe$k.sall))
      k.hall<-c(mean(IVframe$k.hall),mean(Roundsframe$k.hall),mean(Obsframe$k.hall))
      infect<-c(max(IVframe$infect),max(Roundsframe$infect),max(Obsframe$infect))
      care<-c("IV","Rounds","Obs")
      
    }else{
      lambdatemp<-c(mean(IVframe$lambda),mean(Roundsframe$lambda),mean(Obsframe$lambda))
      betatemp<-c(mean(IVframe$beta),mean(Roundsframe$beta),mean(Obsframe$beta))
      durationtemp<-c(mean(IVframe$duration),mean(Roundsframe$duration),mean(Obsframe$duration))
      SHtemp<-c(mean(IVframe$SH),mean(Roundsframe$SH),mean(Obsframe$SH))
      surfconctemp<-c(mean(IVframe$surfconc),mean(Roundsframe$surfconc),mean(Obsframe$surfconc))
      k.salltemp<-c(mean(IVframe$k.sall),mean(Roundsframe$k.sall),mean(Obsframe$k.sall))
      k.halltemp<-c(mean(IVframe$k.hall),mean(Roundsframe$k.hall),mean(Obsframe$k.hall))
      infecttemp<-c(max(IVframe$infect),max(Roundsframe$infect),max(Obsframe$infect))
      caretemp<-c("IV","Rounds","Obs")
      
      lambda<-c(lambda,lambdatemp)
      beta<-c(beta,betatemp)
      duration<-c(duration,durationtemp)
      SH<-c(SH,SHtemp)
      surfconc<-c(surfconc,surfconctemp)
      k.sall<-c(k.sall,k.salltemp)
      k.hall<-c(k.hall,k.halltemp)
      infect<-c(infect,infecttemp)
      care<-c(care,caretemp)
    }
    
 
  }
  
  if (j==1){
    scenario<-rep(sprintf("%s",sim.name),length(lambda))
    frame<-data.frame(lambda=lambda,beta=beta,duration=duration,SH=SH,surfconc=surfconc,k.sall=k.sall,k.hall=k.hall,
                      infect=infect,care=care,scenario=scenario)
  }else{
    scenariotemp<-rep(sprintf("%s",sim.name),length(lambda))
    frametemp<-data.frame(lambda=lambda,beta=beta,duration=duration,SH=SH,surfconc=surfconc,k.sall=k.sall,k.hall=k.hall,
                          infect=infect,care=care,scenario=scenariotemp)
    frame<-rbind(frame,frametemp)
    print(length(frame$lambda)/(j*1500))
          print(j)
  }
  
  #reset directory to parent folder so we can go to correct subfolder within parent folder for next sim run
  setwd(this.dir)
}

require(reshape2)
require(ggplot2)
require(ggpubr)

framecor = subset(frame, select = -c(care,scenario) )

cormat<-round(cor(framecor,method=c("spearman")),2)
melted_cormat<-melt(cormat)
ggplot(data=melted_cormat,aes(x=Var1,y=Var2,fill=value))+geom_tile()+
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white", high = "blue") 


setwd(this.dir)

write.csv(frame,'frame_sensitivity.csv')

#--------------working on tile plots to investigate relationships between------------------------------------------
#---------------input parameters that result in highest infection risks--------------------------------------------------------
A<-ggplot(frame)+geom_tile(aes(x=round(lambda,3),y=round(beta,3),fill=log10(infect)))+
  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")+
  scale_y_continuous(name="Surface-to-Hand Transfer Efficiency")+
  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
  theme_pubr()

B<-ggplot(frame)+geom_tile(aes(x=round(lambda,3),y=round(SH,3),fill=log10(infect)))+
  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")+
  scale_y_continuous(name="Fraction of Hand Contact Area")+
  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
  theme_pubr()

C<-ggplot(frame)+geom_tile(aes(x=round(lambda,3),y=round(log10(surfconc),1),fill=log10(infect)))+
  scale_x_continuous(name="Hand-to-Surface Transfer Efficiency")+
  scale_y_continuous(name=expression("log"[10]*phantom(x)*"Surface Concentration"))+
  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
  theme_pubr()

D<-ggplot(frame)+geom_tile(aes(x=round(beta,3),y=round(log10(surfconc),1),fill=log10(infect)))+
  scale_x_continuous(name="Surface-to-Hand Transfer Efficiency")+
  scale_y_continuous(name=expression("log"[10]*phantom(x)*"Surface Concentration"))+
  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
  theme_pubr()

E<-ggplot(frame)+geom_tile(aes(x=round(SH,3),y=round(log10(surfconc),1),fill=log10(infect)))+
  scale_x_continuous(name="Fraction of Hand Contact Area")+
  scale_y_continuous(name=expression("log"[10]*phantom(x)*"Surface Concentration"))+
  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
  theme_pubr()

G<-ggplot(frame)+geom_tile(aes(x=round(SH,3),y=round(beta,3),fill=log10(infect)))+
  scale_x_continuous(name="Fraction of Hand Contact Area")+
  scale_y_continuous(name="Surface-to-Hand Transfer Efficiency")+
  scale_fill_continuous(name=expression("log"[10]*phantom(x)*"Infection Risk"))+
  theme_pubr()

windows()
ggarrange(A,B,D,C,E,G,common.legend = TRUE)

A<-ggplot(frame)+geom_point(aes(x=lambda,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()
B<-ggplot(frame)+geom_point(aes(x=beta,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()
C<-ggplot(frame)+geom_point(aes(x=surfconc,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()
D<-ggplot(frame)+geom_point(aes(x=k.sall,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()
E<-ggplot(frame)+geom_point(aes(x=k.hall,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()
G<-ggplot(frame)+geom_point(aes(x=SH,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()
H<-ggplot(frame)+geom_point(aes(x=duration,y=infect,colour=scenario))+
  scale_y_continuous(trans="log10")+theme_pubr()

windows()
ggarrange(A,B,C,D,E,G,H,common.legend = TRUE)

  