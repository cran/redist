## ----eval = FALSE-------------------------------------------------------------
#  ## Load data
#  library(redist)
#  data(algdat.pfull)
#  
#  ## Run the simulations
#  mcmc.out <- redist.mcmc(adjobj = algdat.pfull$adjlist,
#                          popvec = algdat.pfull$precinct.data$pop,
#                          nsims = 10000,
#                          ndists = 3)

## ----eval = FALSE-------------------------------------------------------------
#  ## Load data
#  data(algdat.p20)
#  
#  ## -----------------------------------------------
#  ## Run mcmc algorithm - hard population constraint
#  ## (reject any sample where population parity > 20%)
#  ## -----------------------------------------------
#  mcmc.out <- redist.mcmc(adjobj = algdat.p20$adjlist,
#                          popvec = algdat.p20$precinct.data$pop,
#                          nsims = 10000,
#                          popcons = .2,
#                          ndists = 3)
#  
#  ## ---------------------------------------------------------------------
#  ## Run mcmc algorithm - draws from gibbs defined by distance from parity
#  ## (Run with no tempering)
#  ## ---------------------------------------------------------------------
#  mcmc.out.gb <- redist.mcmc(adjobj = algdat.p20$adjlist,
#                             popvec = algdat.p20$precinct.data$pop,
#                             ndists = 3,
#                             nsims = 10000,
#                             constraint = "population",
#                             constraintweights = 5.4)
#  ## Reweight draws back to the uniform distribution
#  mcmc.out.gb <- redist.ipw(mcmc.out.gb, targetpop = .2)
#  
#  ## ---------------------------------------------------------------------
#  ## Run mcmc algorithm - draws from gibbs defined by distance from parity
#  ## (Run with simulated tempering, betas power law sequence of length 10
#  ## from 0 to 1)
#  ## Also including optional beta weights, to upweight prob of transitions
#  ## to colder temperatures
#  ## ---------------------------------------------------------------------
#  betaweights <- rep(NA, 10); for(i in 1:10){betaweights[i] <- 2^i}
#  mcmc.out.st <- redist.mcmc(adjobj = algdat.p20$adjlist,
#                             popvec = algdat.p20$precinct.data$pop,
#                             ndists = 3,
#                             nsims = 10000,
#                             constraint = "population",
#                             constraintweights = 5.4,
#                             temper = TRUE,
#                             betaweights = betaweights)
#  mcmc.out.st <- redist.ipw(mcmc.out.st, targetpop = .2)

## ----eval = FALSE-------------------------------------------------------------
#  ## ----------------------------------------------------------
#  ## Constrain on population and compactness with tempering,
#  ## weight on population = 5.4 while weight on compactness = 3
#  ## Also specifying argument for ssdmat, the distance matrix
#  ## ----------------------------------------------------------
#  mcmc.out.st.multiple <- redist.mcmc(adjobj = algdat.p20$adjlist,
#                                      popvec = algdat.p20$precinct.data$pop,
#                                      ndists = 3,
#                                      nsims = 10000,
#                                      constraint = c("population", "compact"),
#                                      constraintweights = c(5.4, 3),
#                                      ssdmat = algdat.p20$distancemat,
#                                      temper = TRUE,
#                                      betaweights = betaweights)
#  mcmc.out.st <- redist.ipw(mcmc.out.st.multiple, targetpop = .2)

## ----eval = FALSE-------------------------------------------------------------
#  ## ----------------------------------------------------------
#  ## Constrain on population and compactness with parallel tempering,
#  ## weight on population = 5.4 while weight on compactness = 3
#  ## Also specifying argument for ssdmat, the distance matrix.
#  ## In addition, specifying a 20-beta tempering ladder.
#  ## Save file as "redist_mpi.RData"
#  ## ----------------------------------------------------------
#  redist.mcmc.mpi(
#    adjobj = algdat.p20$adjlist,
#    popvec = algdat.p20$precinct.data$pop,
#    ndists = 3,
#    nsims = 10000,
#    constraint = c("population", "compact"),
#    constraintweights = c(5.4, 3),
#    ssdmat = algdat.p20$distancemat,
#    betaseqlength = 20,
#    savename = "redist_mpi",
#    verbose = TRUE
#  )
#  
#  ## Note that reweighting using redist.ipw() currently
#  ## has to happen in a separate analysis file after
#  ## redist.mcmc.mpi() is run and the output file is saved.
#  ## We will change this in a future release to run the
#  ## inverse probability reweighting inside of
#  ## redist.mcmc.mpi().

## ----eval = FALSE-------------------------------------------------------------
#  ## ----------------------------------------------------------
#  ## Constrain on population and compactness with annealing,
#  ## weight on population = 5.4 while weight on compactness = 3
#  ## Also specifying argument for ssdmat, the distance matrix.
#  ## In addition, specifying a 20-beta tempering ladder.
#  ## Save file as "redist_anneal_[t].RData", where t is the
#  ## SLURM array number.
#  ##
#  ## In addition, we no longer specify a "sims" argument.
#  ## This temperature schedule involves simulating at the "hot"
#  ## temperature (beta = 0) for 500 steps ("num_hot_steps"),
#  ## then taking a linear annealing temperature schedule for
#  ## 2000 steps ("num_annealing_steps"), and finally simulating
#  ## at the cold temperature (beta = 1) for 500 steps ("num_cold_steps").
#  ## Finally, we take the last draw of the algorithm and store it as output.
#  ## ----------------------------------------------------------
#  
#  t <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#  
#  redist.mcmc.anneal(
#    adjobj = algdat.p20$adjlist,
#    popvec = algdat.p20$precinct.data$pop,
#    ndists = 3,
#    num_hot_steps = 500, num_annealing_steps = 2000, num_cold_steps = 500,
#    constraint = c("population", "compact"),
#    constraintweights = c(5.4, 3),
#    ssdmat = algdat.p20$distancemat,
#    savename = paste0("redist_anneal_", t)
#  )
#  
#  ## ----------------------------------------------------
#  ## Note that reweighting using redist.ipw() currently
#  ## has to happen in a separate analysis file after
#  ## redist.mcmc.anneal() is run and the output file is saved.
#  ## Below is sample code for combining multiple runs of
#  ## redist.mcmc.anneal() back together, which should be run
#  ## in a separate script. We can then run redist.ipw()
#  ## on the combined simulations.
#  ## ----------------------------------------------------
#  setwd("[directory with files]")
#  algout <- redist.combine.anneal("redist_anneal_")
#  algout_rw <- redist.ipw(algout, targetpop = .2)

