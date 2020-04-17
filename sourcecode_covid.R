
#Mark, for some reason this wasn't working for me, but the line right below does. Hopefully these are doing the same thing.
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)  

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#------------------------------------------ RUN BEHAVIOR PIECE ------------------------------------------------------------------------------------

suppressMessages(suppressWarnings(source("adjust_behaviors_covid.R")))

#---------------------------------- NOW RUN FUNCTION FILE & SIMULATIONS ---------------------------------------------------------------------------

suppressMessages(suppressWarnings(source("markov_chain_stjames_function_covid.R")))

#---------------------------------- PLOTTING -----------------------------------------------------------------------------

require(ggplot2)
require(ggpubr)

#dataframe name is "frameall" with columns for infection risks ("infect"),
#caretype ("care"), and intervention ("intervention")

#all non-zero risks
windows()
ggplot(data=frameall)+geom_violin(aes(y=infect,x=intervention,group=intervention,fill=intervention),colour="black",draw_quantiles = c(0.25, 0.5, 0.75))+
  geom_jitter(aes(y=infect,x=intervention,group=intervention),width=0.1,alpha=0.2)+
  scale_y_continuous(trans="log10",name="Infection Risk")+
  scale_x_discrete(name="Glove Scenario")+
  theme_pubr()+
  theme(legend.position = "none")+
  facet_wrap(~care)


#summary statistics

#------------- IV

#no glove change
summary(frameall$infect[frameall$care=="IV" & frameall$intervention=="No Change"])
sd(frameall$infect[frameall$care=="IV" & frameall$intervention=="No Change"])
quantile(frameall$infect[frameall$care=="IV" & frameall$intervention=="No Change"],probs=0.99)

#glove change
summary(frameall$infect[frameall$care=="IV" & frameall$intervention=="Glove Change"])
sd(frameall$infect[frameall$care=="IV" & frameall$intervention=="Glove Change"])
quantile(frameall$infect[frameall$care=="IV" & frameall$intervention=="Glove Change"],probs=0.99)

#ethanol application
summary(frameall$infect[frameall$care=="IV" & frameall$intervention=="Ethanol"])
sd(frameall$infect[frameall$care=="IV" & frameall$intervention=="Ethanol"])
quantile(frameall$infect[frameall$care=="IV" & frameall$intervention=="Ethanol"],probs=0.99)

#----------- Doctors' Rounds


#no glove change
summary(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="No Change"])
sd(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="No Change"])
quantile(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="No Change"],probs=0.99)

#glove change
summary(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="Glove Change"])
sd(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="Glove Change"])
quantile(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="Glove Change"],probs=0.99)

#ethanol application
summary(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="Ethanol"])
sd(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="Ethanol"])
quantile(frameall$infect[frameall$care=="Rounds" & frameall$intervention=="Ethanol"],probs=0.99)

#----------- Observational Care


#no glove change
summary(frameall$infect[frameall$care=="Obs" & frameall$intervention=="No Change"])
sd(frameall$infect[frameall$care=="Obs" & frameall$intervention=="No Change"])
quantile(frameall$infect[frameall$care=="Obs" & frameall$intervention=="No Change"],probs=0.99)

#glove change
summary(frameall$infect[frameall$care=="Obs" & frameall$intervention=="Glove Change"])
sd(frameall$infect[frameall$care=="Obs" & frameall$intervention=="Glove Change"])
quantile(frameall$infect[frameall$care=="Obs" & frameall$intervention=="Glove Change"],probs=0.99)

#ethanol application
summary(frameall$infect[frameall$care=="Obs" & frameall$intervention=="Ethanol"])
sd(frameall$infect[frameall$care=="Obs"& frameall$intervention=="Ethanol"])
quantile(frameall$infect[frameall$care=="Obs" & frameall$intervention=="Ethanol"],probs=0.99)
