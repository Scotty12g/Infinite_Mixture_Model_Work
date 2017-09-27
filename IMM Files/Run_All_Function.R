##################################################
##### PREPARE THE MODEL SPECIFICATIONS ###########
##################################################

## Read in your data. Add the file path in quotations. Make sure data is 
## structured with the following columns:
## - A column titled "Type" for whether a sample is "train" or "test"
## - A column titled "Site" for the source of a trained individual or
## the collection location of a test individual
## - One column for EACH element with elemental concentrations. Elemental
## concentration data should already be log (x+1) transformed prior to import

data.file =read.csv(file="/Users/scottmorello/Dropbox/Archives/Academic/UMass Boston/Etter Lab/Projects/Trace Element Mussel Connectivity/R Work/IMM Display/IMM Presentation/ICPMS_dat_fortest.csv")
#Check the data fiel before hand. Imports sometimes add a column at begining called "X", so this script removes it.
data.file<-data.file[,-which(colnames(data.file)=="X")] 

# specify parameters
extra.pops = 2
number.of.chains = 2
addaptive.iterations = 10
burn_in.iterations = 10
sample.iterations = 10

##############################
###### RUN THE MODEL #########
##############################

## This runs the model with the above specifications
mod.samples1<-Run_IMM(data.file,extra.pops,number.of.chains,addaptive.iterations,burn_in.iterations,sample.iterations)

## This gets the background info on the data, like how site labels correspond to site numbers ("background.info$source.numbers"),
## and what sites each test individual was colelcted in ("background.info$test.collections")
background.info<-Match_Site_Labels(data.file)

##################################
###### LOOK AT THE RESULTS #######
##################################
##########################################################
# VARIABLES EXPORTED FROM THE MODEL:
# alpha = Concentration parameter controling the probabilities of assignment to an extra group
# samp_assign = Source assignment for sampled individual
# samp_Kprob = Probability that an individual was asssigned to a source
# K_sourced = Whether a source was assigned to
# elem_k = Elemental composition for each source
##########################################################

# the default summary with the MCMC array, only shows the means
mod.samples1

## Save the raw R file in the event you want to access it later
saveRDS(mod.samples1,"Raw_IMM_File.rds")

## To load the file in at another point
# mod.samples1<-readRDS("__FILEPATH___")

## Create a Table of Each Test Individuals Collection Location, and it's Best Assignment
Test_Best_Assign<-background.info$test.collections
Test_Best_Assign$Assignment_Source_Number<-output.summarize(raw.output=mod.samples1,varname="samp_assign",sumtype="mode",thin.rate=1)
Source.List<-c(as.character(background.info$source.numbers$Source_Site_Label),c(nrow(background.info$source.numbers)+c(1:extra.pops)))
Test_Best_Assign$Assignment_Source_Label<-Source.List[as.numeric(Test_Best_Assign$Assignment_Source_Number)]
write.csv(Test_Best_Assign,"Test_Best_Assign.csv")

## Create a Table of Each Test Individuals Collection Location, and it's Assignment Probability to Each Source
Test_Prob_Assign<-background.info$test.collections
Test_Prob_Assign<-cbind(Test_Prob_Assign,output.summarize(raw.output=mod.samples1,varname="samp_Kprob",sumtype="mean",thin.rate=1))
Source.List<-c(as.character(background.info$source.numbers$Source_Site_Label),c(nrow(background.info$source.numbers)+c(1:extra.pops)))
colnames(Test_Prob_Assign)[-c(1:3)]<-Source.List
write.csv(Test_Prob_Assign,"Test_Prob_Assign.csv")

