## Create a function to get the Mode from a vector
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]}

## Create_Dataframe() function to prep the data for the model
##########################################################
## This function will take your data, which should have a column titled "Type" for whether a
## sample is "train" or "test", a column titles "Site" for the source of a trained individual or
## the collection location of a test individual, and then a column for each element with elemental
## concentrations. Elemental concentration data should already be transformed prior to import.
## the "cov.choice" setting corresponds to the covaraince matrix on the prior, but should
## be set to "id" since that's the most stable setting, and the most appropriate based on the
## litteracture. The "group_samps" setting lets the function know whether you'll want the probability of
## assignment to be influced by test individuals gathered from the same site. If so, set it to TRUE,
## but if not, the defgault setting is FALSE and individual test samples are dealt with completely independently.
##########################################################

Create_Dataframe <- function(DataFrame,xtra,cov.choice=c("id","source-univ","source-id","univ"),group_samps=FALSE) {
  
  #### Find out which packages mnight be missing and install them, then load the required packages
  list.of.packages <- c("rjags","plyr","Matrix","abind","reshape","ggplot2","sampling")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  library("rjags")
  library("plyr")
  library("Matrix")
  library("abind")
  library("reshape")
  library("ggplot2")
  library("sampling")
  
  # create a function to get the mode of a vector
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]}
  
  
  ### Separate the data into the "env" data [type (train, test) and site (training incubation
  ### locaiton, test collection location)] and the actual elemental compostion values
  DataFrame.temp<-DataFrame
  
  DataFrame.temp$Type<-factor(DataFrame.temp$Type,levels=c("train","test"))
  source.vec<-as.character(unique(DataFrame.temp$Site[which(DataFrame.temp$Type=="train")]))
  sample.vec<-as.character(unique(DataFrame.temp$Site[which(DataFrame.temp$Type=="test")]))
  sample.only<-which(is.element(sample.vec,source.vec)==FALSE)
  DataFrame.temp$Site<-factor(DataFrame.temp$Site,levels=c(source.vec,sample.only)) # order levels so sources are the first ones
  
  env.cols<-which(colnames(DataFrame.temp)=="Type" | colnames(DataFrame.temp)=="Site")
  DataFrame.temp.env<-DataFrame.temp[,env.cols]
  DataFrame.temp.val<-DataFrame.temp[,-env.cols]
  
  
  ### Separate data beteween training and test sets
  train.n<-which(DataFrame.temp.env$Type=="train")
  test.n<-which(DataFrame.temp.env$Type=="test")
  
  sampledata<-DataFrame.temp.val[test.n,]
  refdata<-DataFrame.temp.val[train.n,]
  
  ### Extract the collection locations for each test individual, and then 
  ### extract the incubation source location for trianing samples
  if (group_samps==FALSE){sampledata.collect<-c(1:length(DataFrame.temp.env$Site[test.n]))}else{sampledata.collect<-as.numeric(DataFrame.temp.env$Site[test.n])}
  refdata.source<-as.numeric(DataFrame.temp.env$Site[train.n])
  extrasite.n<-xtra
  
  
  
  ### Recombine training and test data, and Center and Scale the data, then separate again
  scaleddata<-ddply(as.data.frame(rbind(sampledata,refdata)),.(),colwise(scale))[,-1]
  sampledata<-scaleddata[c(1:nrow(sampledata)),]
  refdata<-scaleddata[-c(1:nrow(sampledata)),]
  
  sampledata<-as.data.frame(sampledata)
  refdata<-as.data.frame(refdata)
  
  
  ### Calculate the mean and sd for each element in each source
  dat.mean.raw<-cbind(data.frame(Source=refdata.source),refdata)
  dat.mean<-ddply(dat.mean.raw,.(Source),colwise(mean))
  dat.mean<-as.matrix(dat.mean[,-1])
  dat.var<-ddply(dat.mean.raw,.(Source),colwise(var))
  dat.var<-as.matrix(dat.var[,-1])
  dat.sd<-sqrt(dat.var)
  
  ### Calculate the number of unique sources we've measured
  refsite.n<-length(unique(refdata.source))
  
  ##### Generate different types of Covariance Priors. One will be chosen when calling the function
  
  # "id": Create an identity matrix array for the prior precision matrix for each source, and extra source
  id.mat.array<-array(NA,c(ncol(sampledata),ncol(sampledata),refsite.n+extrasite.n))
  id.mat<-array(0,c(ncol(sampledata),ncol(sampledata)))
  diag(id.mat)<-1
  for (i in 1:dim(id.mat.array)[3]){id.mat.array[,,i]<-id.mat}
  
  
  # "source-univ": informative source specific precision matrix, and the universal precision matrix for unknown sources
  prec.mat.array<-array(NA,c(ncol(refdata),ncol(refdata),refsite.n+extrasite.n))
  for (i in 1:refsite.n){prec.mat.array[,,i]<-cov(refdata[which(refdata.source==i),])}
  prec.mat.array[,,c((refsite.n+1):(refsite.n+extrasite.n))]<-cov(refdata)
  for (i in 1:(refsite.n+extrasite.n)){prec.mat.array[,,i]<-solve(prec.mat.array[,,i])}
  
  # "source-id" a combination where known sources are the actual precision matrix, and unknown are an identity matrix
  precid.mat.array<-prec.mat.array
  precid.mat.array[,,c((refsite.n+1):(refsite.n+extrasite.n))]<-diag(ncol(refdata))
  
  # "univ": universal covariance for all prior precision matrices
  precall.mat.array<-array(NA,c(ncol(refdata),ncol(refdata),refsite.n+extrasite.n))
  precall.mat.array[,,c(1:(refsite.n+extrasite.n))]<-cov(refdata)
  for (i in 1:(refsite.n+extrasite.n)){precall.mat.array[,,i]<-solve(precall.mat.array[,,i])}
  
  # now choose the correct prior based on the user's selection
  if (cov.choice=="id"){covprior<-id.mat.array}else if (cov.choice=="source-univ"){covprior<-prec.mat.array}else if (cov.choice=="source-id"){covprior<-precid.mat.array}else if (cov.choice=="univ"){covprior<-precall.mat.array}else{covprior<-id.mat.array}
  
  
  #### Put all the data into the correct variables and return a list from the function
  moddat.test<-sampledata
  moddat.base<-refdata
  moddat.nsites<-refsite.n
  moddat.extrasites<-extrasite.n
  moddat.base.mean<-dat.mean
  moddat.base.sd<-dat.sd
  moddat.base.cov<-covprior
  
  
  data = list('elem_n' = ncol(moddat.test),
              'k_max' = moddat.nsites+ moddat.extrasites,
              'samp_n' = nrow(moddat.test),
              'base_n' = nrow(moddat.base),
              'samp_labs' = sampledata.collect,
              'base_labs' =  refdata.source,
              'k_samp' = length(unique(sampledata.collect)),
              'k_base' = length(unique(refdata.source)),
              'samp_dat' = as.matrix(moddat.test),
              'base_dat' = as.matrix(moddat.base),
              'base_means' = moddat.base.mean,
              'base_vars' = moddat.base.sd,
              'base_covs' = moddat.base.cov)
  
  return (data)
}


## Match_Site_Labels() function to extract important info to use later
##########################################################
## This function will take the same data you entered into the "Create_Dataframe" function
## and return a dataframe outlining information about each individual, including
## "Individual" - the row number from the dataset, "Site - the actual site name/code from the
## imported file, "Site_Number" - the site number that corresponds to the site name, since
## the model will not output numbers instead od actual site names, and "Type" - whether the
## individual was a "test" or "train" sample.
##########################################################
Match_Site_Labels<- function(DataFrame){
  
  
  #### Find out which packages mnight be missing and install them, then load the required packages
  list.of.packages <- c("rjags","plyr","Matrix","abind","reshape","ggplot2","sampling")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  library("rjags")
  library("plyr")
  library("Matrix")
  library("abind")
  library("reshape")
  library("ggplot2")
  library("sampling")
  
  # create a function to get the mode of a vector
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]}
  
  
  ### Separate the data into the "env" data [type (train, test) and site (training incubation
  ### locaiton, test collection location)] and the actual elemental compostion values
  DataFrame.temp<-DataFrame
  
  DataFrame.temp$Type<-factor(DataFrame.temp$Type,levels=c("train","test"))
  source.vec<-as.character(unique(DataFrame.temp$Site[which(DataFrame.temp$Type=="train")]))
  sample.vec<-as.character(unique(DataFrame.temp$Site[which(DataFrame.temp$Type=="test")]))
  sample.only<-which(is.element(sample.vec,source.vec)==FALSE)
  DataFrame.temp$Site<-factor(DataFrame.temp$Site,levels=c(source.vec,sample.only))
  
  env.cols<-which(colnames(DataFrame.temp)=="Type" | colnames(DataFrame.temp)=="Site")
  DataFrame.temp.env<-DataFrame.temp[,env.cols]
  
  ## Create a dataframe with all of the information for the function to return
  matched.labels<-data.frame(Individual=c(1:nrow(DataFrame)),Site=DataFrame.temp.env$Site,Site_Number=as.numeric(DataFrame.temp.env$Site),Type=DataFrame.temp.env$Type)
  test.collections<-subset(matched.labels,Type=="test")[,-c(4)]
  test.collections$Individual<-c(1:nrow(test.collections))
  colnames(test.collections)<-c("Individual_Number","Collection_Site_Label","Collection_Site_Number")
  source.numbers<-unique(subset(matched.labels,Type=="train")[,-c(1,4)])
  colnames(source.numbers)<-c("Source_Site_Label","Source_Site_Number")
  meta.data<-list(source.numbers=source.numbers,test.collections=test.collections)
  return(meta.data)
}



## The actual Mixture Model to run in JAGS
## "Morello_InMixMod_wBaseline.model"
##########################################################
## This model, written in JAGS format, which is the same as WinBugs,
## takes inputs from the "Create_Dataframe" function to fit
## an infinite mixture model with baseline elemental data for train
## individuals, and assign test individuals to either a trained
## source (using element specific means and SDs as priors, and whichever
## covariance prior was specified when calling Create_Dataframe), or
## an untrained source (the max number of untrained sources possible is
## also specified in Create_Dataframe)
##########################################################
# EXTRA INFO ON VARIABLE NAMES:
# the following are variables used in then model, and fed 
#in from the Create_Dataframe output:
# elem_n = number of elements (i.e., columns) to be analyzed (value)
# k_max = Maximum number of sources to consider in the model (value: baseline + extra sources)
# samp_n = Number of individuals of unknown origin in the sample data to allocate to sources (value; number of rows in samp_dat)
# base_n = Number of individuals of known origin in the baseline data (value; number of rows in base_dat)
# samp_labs = The collection location for each sample individual (vector; length samp_n)
# base_labs = The known source of each baseline individual (vector; length base_n)
# k_samp = The number of collection locations for the sample data (value)
# k_base = The number of known sources included in the baseline data (value)
# elem_ident = idenity matrix of dimension elem_n x elem_n (matrix; dimensions elem_n x elem_n)
# samp_dat = sample elemental composition data (matrix; dimensions samp_n x elem_n)
# base_dat = baseline elemental composition data (matrix; dimensions base_n x elem_n)
# base_means = mean composition of each element by source from the baseline data (matrix; dimensions k_base x elem_n)
# base_vars = standard deviaiton for each element's composition by source from the baseline data (matrix; dimensions k_base x elem_n)
# base_covs = positive definte covariance matrix for elemental composition by source from the baseline data (array; dimensions elem_n x elem_n x k_base)
##########################################################
{
  
  Morello_InMixMod_wBaseline.model<-"model {
  #----------------------------------------------------------------------------------------------------------------
  # SOURCE PROPERTIES
  #----------------------------------------------------------------------------------------------------------------
  # The mean elemental composition value for each element for each source (K_mu).
  # Start with the known sources, and base their value on a normal distribution, and using the baseline mean and 
  # sd as priors for mu (mean) and percision (1/variance).
  # Then, for the unknown sources, approximate the composition for each element with an uniform prior
  # between -2 and 2, since that's what centered and scaled elemental concentration data (which is 
  # how concentrations are fed into the model) should range between.
  #----------------------------------------------------------------------------------------------------------------
  
  # Known sources
  for (Ki in 1:k_base){
  for (Ei in 1:(elem_n)){
  K_mu[Ki,Ei]~dnorm(base_means[Ki, Ei],(1/(base_vars[Ki,Ei]*base_vars[Ki,Ei])))
  
  
  }
  }
  
  # Unknown sources
  for (Ki in (k_base+1):k_max){
  for (Ei in 1:(elem_n)){
  K_mu[Ki,Ei]~dunif(-2,2) 
  }
  }
  
  #----------------------------------------------------------------------------------------------------------------
  # Prior on the variance-covariance matrix for each source where tau is approximated by a Wishart multivariate
  # distribution where Scale is the degrees of freedom (elem_n + 1) and matrix R (base_covs) for each source is the 
  # square-symetric positive definite identity matrix, or actual covariance matrix, depending on 
  # what is called. Scaled inverse Wishart prior
  #----------------------------------------------------------------------------------------------------------------
  for (i in 1:k_max){
  K_tau[1:elem_n,1:elem_n,i] ~ dwish(base_covs[,,i],(elem_n+1))
  }
  
  
  #----------------------------------------------------------------------------------------------------------------
  # PROBABILITY OF SOURCE ASSIGNMENT
  #----------------------------------------------------------------------------------------------------------------
  # Known Sources:
  #
  # Start by indicating an equal probability of assignment among the k_base known sources
  #----------------------------------------------------------------------------------------------------------------
  for (Ki in 1:k_base){
  prob_known[Ki]<-1/k_base
  }
  
  #----------------------------------------------------------------------------------------------------------------
  # Unknown Sources:
  #
  # Now for the remaining, possioble unknown, sources, we specify prior for the concentration parameter (alpha) 
  # as comong from an inverse gamma distributuon. We start by approximating alpha_0 from a gamma density function 
  # of shape (r) 1 and rate (mu) 1, and take the inverse (to get inverse-gamma distribution). We also 
  # truncate ('T' funciton) alpha_0 approximation to lie between 0.1 and 3.4 to constrain the alpha from
  # getting too high
  #----------------------------------------------------------------------------------------------------------------
  
  alpha_0 ~ dgamma(1,1) T(0.01,3.4)
  alpha<-1/(alpha_0)
  
  #----------------------------------------------------------------------------------------------------------------
  # Then build priors for the probability of belonging to additional, unknown, sources through an itterative 
  # stick breaking process, where the prior propbability of belonging to each new source is a function of a beta
  # distribution with shape of 1 and the concentration parameter alpha. Then use the first break-point the model
  # generated as a starting point in the stick-breaking process. We itterate the stick breaking so that
  # each possible new source depends on the size of the pervious source, creating smaller and smaller pieces/probabilites.
  #----------------------------------------------------------------------------------------------------------------
  
  # the first stick breaking point for is approximated by a beta distribution of shape parameters 1 and concentration parameter alpha
  
  for (Ki in (k_base+1):(k_max-1)){
  break_point[Ki-k_base]~dbeta(1,alpha)
  }
  
  
  # use the inital break_point as the first stick break, and then itteratively break the stick based on the prior break location
  stick_weights[1]<-break_point[1] 
  for (j in (k_base+2):(k_max-1)){
  stick_weights[j-k_base]<-break_point[j-k_base]*(1-break_point[j-k_base-1])*stick_weights[j-k_base-1]/break_point[j-k_base-1]
  }
  
  # let the last new source be the leftover stick piece, and sort the sticks by size
  stick_weights[k_max-k_base]<-1-sum(stick_weights[1:(k_max-k_base-1)])
  stick_sorted<-sort(stick_weights[1:(k_max-k_base)])
  
  
  #----------------------------------------------------------------------------------------------------------------
  # Now combine the known and unknown prior weights. During the stick breaking process, we sorted the weights
  # so that the first k_base weights (the ones for the known sources) were the highest. We want to replace those
  # so that known sources are chosen based on a Dirichlet distribution and equal probabilites for each known source. 
  # After replacing the weights for known sources, we then go through and standardize the full vector of probabilites
  # so that known sources are much more likely than the unknown sources, and we base how much more likely on the
  # alpha value chosen earlier in the model. At the end of this process, total probability of assignment to any source
  # sums to 1.
  #----------------------------------------------------------------------------------------------------------------
  
  
  
  for (Ki in 1:k_samp){
  # Assuming a dirichlet process model, and flat prior, the k_base first Ki would have this total weights
  base_weights[Ki,1:k_base]~ddirch(prob_known[1:k_base])
  for (j in 1:k_base){
  K_prob[Ki,j]<-base_weights[Ki,j]*(1-(1-(1/(1+alpha)))^k_base)
  # Now we adjust the k_base first sources based on concentration parameter alpha
  }
  # Now for the remaining sources, we use the weights from the stick breaking process before, and adjust based on alpha
  for (j in (k_base+1):k_max){
  K_prob[Ki,k_max-j+1+k_base]<-stick_sorted[j-k_base]*((1-1/(1+alpha))^k_base)
  # at the end, the total probabilioty sums to 1, with known sources having higher probs based on alpha
  }
  }
  
  
  #----------------------------------------------------------------------------------------------------------------
  # Reallocation Process
  # where samp_n is the categorical variable of reallocation of sampled individuals, samp_assign (the assignment probabilities for an individual) is described by a catagorical discrete univariate distribution based on the Dirichlet distributed weights, the stick breaking process, and alpha
  # Elemental compositions of samples (samp_dat) are described by a multivariate normal distribution related to the mean elemetnal concetrations for a source and the variance covvariance matrix.
  # This multivariate normal distribution is linked to the baseline data (base_dat).
  #----------------------------------------------------------------------------------------------------------------
  
  # assigning individuals based on the probabilities outlined before, and their location within the multivariate normal distributions of elemental concentration data
  for (ind in 1:samp_n){
  samp_assign[ind] ~ dcat(K_prob[samp_labs[ind],1:k_max])
  samp_dat[ind,1:elem_n] ~ dmnorm(K_mu[samp_assign[ind],1:elem_n], K_tau[1:elem_n,1:elem_n, samp_assign[ind]])
  }
  
  # now do the same for the reference sample, but we know the sources (base_labs) so we don't need to use dcat to describe it. Instead, this links the distributions to the actual baseline data
  for (ind in 1:base_n){
  base_dat[ind,1:elem_n] ~ dmnorm(K_mu[base_labs[ind], 1:elem_n], K_tau[1:elem_n,1:elem_n,base_labs[ind]])
  }
  
  #----------------------------------------------------------------------------------------------------------------
  # Summary statistics
  #----------------------------------------------------------------------------------------------------------------
  
  # this finds the source with the highest probability of assignment for each individual for each itteration
  for (i in 1:samp_n){
  for (j in 1:k_max){
  samp_Kprob[i,j]<-equals(samp_assign[i],j) # test for equality
  }
  }
  for (j in 1:k_max){
  K_sourced[j]<-step(sum(samp_Kprob[,j])-1) # where 'step' is a test for x>=0
  # this summarizes if any individuals were assigned to a source Ki during an itteration (binary - 1=yes, 0=no)
  
  #----------------------------------------------------------------------------------------------------------------
  # Prediction of elemental composition for each element and each source
  #----------------------------------------------------------------------------------------------------------------
  
  elem_k[j,1:elem_n]~dmnorm(K_mu[j,1:elem_n], K_tau[1:elem_n,1:elem_n,j])
  }
  
}"

  }


## RUNNING THE MODEL
##########################################################
# The following code runs JAGS through rjags. The model is called,
# and the data are fed in, along with the number of chains you would
# like to run the model. For this first run, the model 
# is run through a run-up/adaptation set of itterations which
# is used by JAGS to set the optimal paramters for the MCMC (e.g., 
# choosing the correct step-sizes). Following the adaptation phase,
# the model is updated with the number of burn-in itterations.
# finally, the model is run for the number of desired itterations to sample over,
# and the variables to monitor during sampling (i.e., what the model
# will output at the end) are set (e.g., samp_assign, samp_Kprob, K_sourced, elem_k)
##########################################################
# EXTRA INFO ON VARIABLE NAMES EXPORTED FROM MODEL:
# alpha = Concentration parameter controling the probabilities of assignment to an extra group
# samp_assign = Source assignment for sampled individual
# samp_Kprob = Probability that an individual was asssigned to a source
# K_sourced = Whether a source was assigned to
# elem_k = Elemental composition for each source
##########################################################

# mod.model1 <- jags.model(textConnection(____MODEL_____), # call the model
#           data = ___OUTPUT_FROM_Create_Dataframe____,
#           n.chains = ___, #how many parallel chains to run
#           n.adapt = ___) # How many samples should be in the adaptive sampling period of the chain

# update(mod.model1, ____) # run the model with _____ itterations of burn in


# mod.samples1<-jags.samples(mod.model1, # now take this simulation run
#                           c('samp_Kprob','samp_assign','K_sourced','elem_k','alpha'), # sample these variables
#                           _______) # take _______ many samples/itteratons

## If you run this sctript (i.e., just type "Run_IMM()" and it will run), it will take care of the running the model as seen above
Run_IMM <- function (data.file=data.file, extra.pops=extra.pops, number.of.chains=number.of.chains,addaptive.iterations=addaptive.iterations,burn_in.iterations=burn_in.iterations,sample.iterations=sample.iterations){  
  #### THIS IS HOW YOU READ IN A CSV FILE #########
  sim.dat<-data.file
  
  ### NOW PRE-PROCESS THE DATA ############
  sim.dat.proc<-Create_Dataframe(sim.dat,extra.pops,"id")
  
  #### GET BACKGROUND INFO FOR EACH INDIVIDUAL AND HOW SITE NUMBERS RELATE TO ACTUAL SITE LABELS ###
  mod.model1 <- jags.model(textConnection(Morello_InMixMod_wBaseline.model), # call the model
                           data = sim.dat.proc,
                           n.chains = number.of.chains, #how many parallel chains to run
                           n.adapt = addaptive.iterations) # How many samples should be in the adaptive sampling period of the chain
  
  update(mod.model1, burn_in.iterations) # run the model itterations of burn in
  
  
  mod.samples1<-jags.samples(mod.model1, # now take this simulation run
                             c('samp_Kprob','samp_assign','K_sourced','elem_k','alpha'), # sample these variables
                             sample.iterations) # take  samples/itteratons
  
  return(mod.samples1)}

############ To save a raw R file ##############
# saveRDS(mod.samples1,"___FILENAME____")
# raw.output<-readRDS("__FILEPATH___")

## Final Output Summary
##########################################################
## Just calling the output will give you a summary (means) of all variables monitored across
## all chains and interations run.
##########################################################
# raw.output


## output.summarize() Function
##########################################################
## This function will take your raw model output, which is a mcmc array and difficult to
## manipulate, and summarize a specific variable based on whatever thinning you would
## like for the data. Note: a thin.rate of 1 (the default if you don't type anything), will
## give you all itterations from the model output
##########################################################
output.summarize<-function(raw.output,varname=c("K_sourced","alpha","elem_k","samp_Kprob","samp_assign"),sumtype=c("mean","mode"),thin.rate=1){
  dim.n<-length(dim(raw.output[[varname]]))
  thin.mean<-function(dat,thin.rate1=thin.rate){mean(dat[seq(1,length(dat),1/thin.rate1)])}
  thin.mode<-function(dat,thin.rate1=thin.rate){getmode(dat[seq(1,length(dat),1/thin.rate1)])}
  if(sumtype=="mean"){outputsum<-apply(raw.output[[varname]],c(1:(dim.n-2)),thin.mean)}else{
    if (sumtype=="mode"){outputsum<-apply(raw.output[[varname]],c(1:(dim.n-2)),thin.mode)}else{outputsum<-apply(raw.output[[varname]],c(1:(dim.n-2)),thin.mean)}}
  return(outputsum)}


## Individual variable extraction
##########################################################
## These functions will take your raw model output, which is a mcmc array and difficult to
## manipulate, and extract specified elements from specificied variables. Note that each
## function starts with a call for the raw model output file "raw.output", but may also
## ask for additional info (e.g. "elem_k raw" also asks which element number you would like 
## to extract (the column from the original dataframe), as well as which site number).
##########################################################

## Extract raw K_sourced data (binary) based on site.
##########################################################
## K_sourced = Whether a source had any individuals assigned to it
## over one itteration (binary).
## In the "output.summarize" output, coumns are the sources from the
## model (see "Match_Site_Labels" function to find which known
## sources correspond to which site numbers in the model).
##########################################################

get.k_sourced.raw<-function(raw.output,site){
  itteration.n<-dim(raw.output$K_sourced)[2]
  chian.n<-dim(raw.output$K_sourced)[3]
  K_sourced.raw<-as.vector(raw.output$K_sourced[site,c(1:itteration.n),c(1:chian.n)])
  return(K_sourced.raw)}

## Extract the raw alpha values generated in the model
##########################################################
## alpha = Concentration parameter controling the probabilities of 
## assignment to an extra group
##########################################################

get.alpha.raw<-function(raw.output){
  itteration.n<-dim(raw.output$alpha)[2]
  chian.n<-dim(raw.output$alpha)[3]
  alpha.raw<-as.vector(raw.output$alpha[1,c(1:itteration.n),c(1:chian.n)])
  return(alpha.raw)}

## Extract the raw elem_k data (a concentration) based on site number and element number
##########################################################
## elem_k = Elemental composition for each source,
## In the "output.summarize" output, columns are elements ordered the same as
## the input data, and the rows are the sources from the
## model (see "Match_Site_Labels" function to find which known
## sources correspond to which site numbers in the model)
##########################################################
get.elem_k.raw<-function(raw.output, site.n, element.n){
  itteration.n<-dim(raw.output$elem_k)[3]
  chian.n<-dim(raw.output$elem_k)[4]
  element.raw<-as.vector(raw.output$elem_k[site.n,element.n,c(1:itteration.n),c(1:chian.n)])
  return(element.raw)}


## Extract the raw samp_Kprob data (a probability) based on individual and the site number
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
get.samp_Kprob.raw<-function(raw.output, individual.n, site.n){
  itteration.n<-dim(raw.output$samp_Kprob)[3]
  chian.n<-dim(raw.output$samp_Kprob)[4]
  samp_Kprob.raw<-as.vector(raw.output$samp_Kprob[individual.n,site.n,c(1:itteration.n),c(1:chian.n)])
  return(samp_Kprob.raw)}


## Extract the raw samp_assign data (site number) based on individual and the site number
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
get.samp_assign.raw<-function(raw.output, individual.n){
  itteration.n<-dim(raw.output$samp_assign)[2]
  chian.n<-dim(raw.output$samp_assign)[3]
  samp_assign.raw<-as.vector(raw.output$samp_assign[individual.n,c(1:itteration.n),c(1:chian.n)])
  return(samp_assign.raw)}



## mcmc.raw.plot() Function  
##########################################################
## Insert the output from one of the "Individual variable extraction" functions
## above (e.g., "get.samp_Kprob.raw"), the number of chains the were run with the model,
## and the desired thinning rate, and this will plot the MCMC as value by itteration.
## Note: a thin.rate of 1 (the default if you don't type anything), will
## give you all itterations from the model output
##########################################################
mcmc.raw.plot<-function(raw.values,chains.run=1,thin.rate=1){
  itteration.n<-length(raw.values)/chains.run
  par(mfrow=c(chains.run,1))
  for(i in 1:chains.run){
    vals.temp<-raw.values[c(((itteration.n*i)-(itteration.n-1)):(itteration.n*i))]
    itters.temp<-c(1:itteration.n)
    
    thin.raw.length<-length(vals.temp)
    thin.dat<-vals.temp[seq(1,thin.raw.length,1/thin.rate)]
    thin.length<-itters.temp[seq(1,thin.raw.length,1/thin.rate)]
    
    plot(thin.dat~thin.length,xlab="itteration",ylab="value",type="l",main=paste("Chain",i))}
  
  par(mfrow=c(1,1))}

