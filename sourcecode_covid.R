
#Mark, for some reason this wasn't working for me, but the line right below does. Hopefully these are doing the same thing.
#this.dir <- dirname(parent.frame(2)$ofile)
#setwd(this.dir)  

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#------------------------------------------ Run behavior code, simulation function code, and automation ---------------------------------------------------------------------------

suppressMessages(suppressWarnings(source("markov_chain_stjames_function__cumulativerisk_covid_automated.R")))


