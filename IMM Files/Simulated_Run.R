#### THIS IS HOW YOU READ IN A CSV FILE #########
# Read in the CSV file. You'll need to change the path I have there
sim.dat<-read.csv("/Users/scottmorello/Desktop/simdat.csv")
#Check the data file before hand. Imports sometimes add a column at begining called "X", so this script removes it.
sim.dat<-sim.dat[,-which(colnames(sim.dat)=="X")] 

### NOW PRE-PROCESS THE DATA ############
sim.dat.proc<-Create_Dataframe(sim.dat,4,"id")

#### GET BACKGROUND INFO FOR EACH INDIVIDUAL AND HOW SITE NUMBERS RELATE TO ACTUAL SITE LABELS ###
sim.dat.proc.info<-Match_Site_Labels(sim.dat)


mod.model1 <- jags.model(textConnection(Morello_InMixMod_wBaseline.model), # call the model
                         data = sim.dat.proc,
                         n.chains = 2, #how many parallel chains to run
                         n.adapt = 1) # How many samples should be in the adaptive sampling period of the chain

update(mod.model1, 1) # run the model with 1000 itterations of burn in


mod.samples1<-jags.samples(mod.model1, # now take this simulation run
                           c('samp_Kprob','samp_assign','K_sourced','elem_k','alpha'), # sample these variables
                           2) # take 2000 many samples/itteratons


### Now lets look at the results

# the default summary with the MCMC array
mod.samples1

####### SUMMARIZING THE FULL DATASET IN YOUR OWN #########
## use the "output.summarize()" Function
##########################################################
## This function will take your raw model output, which is a mcmc array and difficult to
## manipulate, and summarize a specific variable based on whatever thinning you would
## like for the data. Note: a thin.rate of 1 (the default if you don't type anything), will
## give you all itterations from the model output
##########################################################
## Set the "varname" to any of the variables exported from the model:
## "K_sourced","alpha","elem_k","samp_Kprob","samp_assign"
##########################################################
## set the "sumtype" to "mean" or "mode"
##########################################################

# Example:
# Here, we look at the mode assignment
predicted.sources<-output.summarize(raw.output=mod.samples1,varname="samp_assign",sumtype="mode",thin.rate=1)

# now lets look at it as a confusion matrix based on our actual data generation
sim.act.pred.table<-table(sim.dat.proc.info$test.collections$Collection_Site_Number,predicted.sources)

sim.act.pred.table.percent<-as.matrix(sim.act.pred.table/rowSums(sim.act.pred.table))

# plot a heatmap
heatmap(sim.act.pred.table.percent, Rowv=NA, Colv=NA, col = rev(heat.colors(256)),xlab="Predicted",ylab="Actual")


###### GET THE RAW VALUES FROM AN ELEMENT OF THE OUTPUT ################
##########################################################
## These functions will take your raw model output, which is a mcmc array and difficult to
## manipulate, and extract specified elements from specificied variables. Note that each
## function starts with a call for the raw model output file "raw.output", but may also
## ask for additional info (e.g. "elem_k raw" also asks which element number you would like 
## to extract (the column from the original dataframe), as well as which site number).
##########################################################


## Extract raw K_sourced data (binary) based on site.
## use the "get.k_sourced.raw()" function
##########################################################
## K_sourced = Whether a source had any individuals assigned to it
## over one itteration (binary).
## In the "output.summarize" output, coumns are the sources from the
## model (see "Match_Site_Labels" function to find which known
## sources correspond to which site numbers in the model).
##########################################################
## set "raw.output" to the model output
##########################################################
## set "site" to the site you want to sumamrize
##########################################################

## Extract the raw alpha values generated in the model
## use the "get.alpha.raw" function
##########################################################
## alpha = Concentration parameter controling the probabilities of 
## assignment to an extra group
##########################################################
## set "raw.output" to the model output
##########################################################


## Extract the raw elem_k data (a concentration) based on site number and element number
# use the "get.elem_k.raw" function
##########################################################
## elem_k = Elemental composition for each source,
## In the "output.summarize" output, columns are elements ordered the same as
## the input data, and the rows are the sources from the
## model (see "Match_Site_Labels" function to find which known
## sources correspond to which site numbers in the model)
##########################################################
## set "raw.output" to the model output
##########################################################
## set "site.n" to the site number you want to sumamrize
##########################################################
## set "element.n" to the element number you want to sumamrize
##########################################################


## Extract the raw samp_Kprob data (a probability) based on individual and the site number
# use the "get.samp_Kprob.raw" function
##########################################################
## samp_Kprob = Probability that an individual was asssigned to a source
## In the "output.summarize" output, columns are the sources from the
## model (see "Match_Site_Labels" function to find which known
## sources correspond to which site numbers in the model), and rows
## are the individuals from the test set (generate a summary with the
## "Match_Site_Labels" function, and then subset only the "test"
## individuals. That should give the correct order that corresponds to
## samp_Kprob data).
##########################################################
## set "raw.output" to the model output
##########################################################
## set "individual.n" to the test individual number you want to sumamrize
##########################################################
## set "site.n" to the site number you want to sumamrize
##########################################################



## Extract the raw samp_assign data (site number) based on individual and the site number
## use the "get.samp_assign.raw" function
##########################################################
## samp_assign = Source assignment for sampled individual
## In the "output.summarize" output, values are the sources from the
## model (see "Match_Site_Labels" function to find which known
## sources correspond to which site numbers in the model)
## with the highest probabilty of assignment for that itteraction, 
## and the columns are the individuals 
## from the test set (generate a summary with the
## "Match_Site_Labels" function, and then subset only the "test"
## individuals. That should give the correct order that corresponds to
## samp_Kprob data).
##########################################################
## set "raw.output" to the model output
##########################################################
## set "individual.n" to the test individual number you want to sumamrize
##########################################################


# Example:
# Maybe we want the raw values from 3rd element concretraiton from the 2nd site
site2_element3.raw<-get.elem_k.raw(raw.output=mod.samples1,site.n = 2,element.n = 3)

# Or maybe the same with the assignment of individual 25 to extra site 5
ind25to5.raw<-get.samp_Kprob.raw(raw.output=mod.samples1,individual.n = 25,site.n = 5)


## plot the raw MCMC chian
## use the "mcmc.raw.plot()" function  
# mcmc.raw.plot<-function(raw.values,chains.run=1,thin.rate=1)
##########################################################
## Insert the output from one of the "Individual variable extraction" functions
## above (e.g., "get.samp_Kprob.raw"), the number of chains the were run with the model,
## and the desired thinning rate, and this will plot the MCMC as value by itteration.
## Note: a thin.rate of 1 (the default if you don't type anything), will
## give you all itterations from the model output
##########################################################
## set "raw.values" to one of the raw outputs from above (e.g. output from "get.samp_assign.raw")
##########################################################
## set "chains.run" to the number of chains you ran in your model
##########################################################
  
# Example:
# we can plot the 2 chains of MCMC itterations for the raw values from the above examples

# the "get.elem_k.raw" ourput for the 3rd element concretraiton from the 2nd site
mcmc.raw.plot(raw.values = site2_element3.raw, chains.run = 2, thin.rate = 1)
# then test the autocorrelation for a maximum lag of 10
site2_element3.raw<-acf(x=site2_element3.raw,lag.max=10,type="correlation",plot=TRUE)


# the "get.samp_Kprob.raw" ourput for the assignment of individual 25 to site 5
mcmc.raw.plot(raw.values = ind25to5.raw, chains.run = 2, thin.rate = 1)
# then test the autocorrelation for a maximum lag of 100
ind25to5.ac<-acf(x=ind25to5.raw,lag.max=100,type="correlation",plot=TRUE)

