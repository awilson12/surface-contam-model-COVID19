
#Mark, for some reason this wasn't working for me, but the line right below does. Hopefully these are doing the same thing.
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)  

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#------------------------------------------ RUN BEHAVIOR PIECE ------------------------------------------------------------------------------------

suppressMessages(suppressWarnings(source("adjust_behaviors_covid.R")))

#---------------------------------- NOW RUN FUNCTION FILE & SIMULATIONS ---------------------------------------------------------------------------

suppressMessages(suppressWarnings(source("markov_chain_stjames_function__cumulativerisk_covid_automated.R")))


