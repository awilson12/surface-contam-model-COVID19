
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)    

if("truncdist" %in% rownames(installed.packages())==FALSE){install.packages("truncdist"); require(truncdist)}else{require(truncdist)}
if("gsl" %in% rownames(installed.packages())==FALSE){install.packages("gsl"); require(gsl)}else{require(gsl)}

# CONSTANTS TO BE USED 

SIM.iter <- 100     # Making a master Monte Carlo iteration value

# AUTOMATION ----
# Next line will be the names of the simulations that will be run. DCLP is Doffing correct and low patient numbers
# DCHP is Doffing correct high patient numbers ... for the combinations of actions that can ocurr.
SIM <- c("Doffing Correct","Doffing Incorrect","Low No Patients","High No Patients", "DCLP", "DCHP","DILP","DIHP")   
NUM.SIM <- length(SIM)     # Count the number of iterations for the automated simulations

for(j in 1:NUM.SIM)
{
  sim.num <- j; sim.name <- SIM[j]
  
  if(sim.name=="Doffing Correct"){if(dir.exists("DC")==FALSE){dir.create("DC"); setwd("DC")}else{setwd("DC")}}
  if(sim.name=="Doffing Incorrect"){if(dir.exists("DI")==FALSE){dir.create("DI"); setwd("DI")}else{setwd("DI")}}
  if(sim.name=="Low No Patients"){if(dir.exists("LP")==FALSE){dir.create("LP"); setwd("LP")}else{setwd("LP")}}
  if(sim.name=="High No Patients"){if(dir.exists("HP")==FALSE){dir.create("HP"); setwd("HP")}else{setwd("HP")}}
  if(sim.name=="DCLP"){if(dir.exists("DCLP")==FALSE){dir.create("DCLP"); setwd("DCLP")}else{setwd("DCLP")}}
  if(sim.name=="DCHP"){if(dir.exists("DCHP")==FALSE){dir.create("DCHP"); setwd("DCHP")}else{setwd("DCHP")}}
  if(sim.name=="DILP"){if(dir.exists("DILP")==FALSE){dir.create("DILP"); setwd("DILP")}else{setwd("DILP")}}
  if(sim.name=="DIHP"){if(dir.exists("DIHP")==FALSE){dir.create("DIHP"); setwd("DIHP")}else{setwd("DIHP")}}

  setwd(this.dir)
  
}