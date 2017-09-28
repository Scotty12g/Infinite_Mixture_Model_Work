HOW TO USE THE FILES:

The Datasets: “ICPMS_dat_fortest.csv” is a sample ICPMS data to load in and run. “simdat.csv” is the simulated dataset I walked us through.

The R Files: “IMM_Functions.R” is the file you open in R, and just select all the text and run it. This will load in all of the functions you need. BE SURE TO EXECUTE THIS ENTIRE FILE (“IMM_Functions.R”) BEFORE RUNNING ANY OF THE NEXT TWO SCRIPTS. The next file, “Simulated_Run.R”, is a step by step process of running the model and looking at the summarized and raw outputs as you need (very detailed running of the model). Note that you can load any dataset (e.g., “simdat.csv” file, or the “ICPMS_dat_fortest.csv” file) into this and have it work. To work with the model very basically, and just export and save all of the assignment probabilities for all the test individuals, I’ve included a file titled  “Run_All_Function.R”. All you’ll need to do in this file is import the file with the file path (replace my file path to do that), and specify some values on how you want the model to run (e.g., # of iterations, # of extra sources, etc.). Then, if you execute all the code, it will automatically save the following files to your working directory:

- Raw_IMM_File.rds: The raw R file output from the model… you can only read this inside R using the readRDS("___filepath__") function, but it’s a good way to save the run for later use

- Test_Best_Assign.csv: A table showing the best source assignment for each test individual, and where each test individual was gathered (i.e., where it settled)

- Test_Prob_Assign.csv: A table showing where each test individual was gathered, and the probability of assignment to all sources averaged across all model iterations

The Infinite Mixture Model in the "IMM_Functions.R" file is implement in JAGS through R (via the "rjags" package). The JAGS model has been altered and repurposed from original code by:

- Angers, C. F. R., & Rennes, C. F. R. Using otolith microchemistry within Bayesian reallocation models to explore the Allis shad (Alosa alosa) metapopulation functioning.
