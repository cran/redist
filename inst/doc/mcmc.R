## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(redist)
library(ggplot2)
library(dplyr)
library(patchwork)
# set seed for reproducibility
set.seed(1)

## -----------------------------------------------------------------------------
data(fl25)
data(fl25_enum)
plan <- fl25_enum$plans[, 7241]
adj <- redist.adjacency(fl25)
set.seed(1)
sims <- redist.flip(adj = adj, init_plan = plan, total_pop = fl25$pop, 
                    nsims = 6, pop_tol = 0.10, 
                    constraint = 'compact', constraintweights = 0.02,
                    compactness_metric = 'edges-removed')

## ---- echo = FALSE------------------------------------------------------------
redist.plot.map(shp = fl25, plan = plan) + redist.plot.adj(shp = fl25, plan = plan, plot_shp = FALSE) 

## ---- echo = FALSE------------------------------------------------------------
redist.plot.adj(shp = fl25, plan = plan, plot_shp = FALSE, title = 'Initialization') +
  redist.plot.adj(shp = fl25, plan = sims$plans[, 1], plot_shp = FALSE, title = 'First Iteration') +
  redist.plot.adj(shp = fl25, plan = sims$plans[, 2], plot_shp = FALSE, title = 'Second Iteration') +
  redist.plot.adj(shp = fl25, plan = sims$plans[, 3], plot_shp = FALSE, title = 'Third Iteration') +
  redist.plot.adj(shp = fl25, plan = sims$plans[, 4], plot_shp = FALSE, title = 'Fourth Iteration') +
  redist.plot.adj(shp = fl25, plan = sims$plans[, 5], plot_shp = FALSE, title = 'Fifth Iteration')

## -----------------------------------------------------------------------------
data(iowa)
redist.plot.map(iowa, plan = cd_2010)

## -----------------------------------------------------------------------------
adj <- redist.adjacency(shp = iowa, plan = iowa$cd_2010)

## ----seed, include = FALSE----------------------------------------------------
# set the seed so that within the vignette the results don't change
set.seed(1)

## -----------------------------------------------------------------------------
sims <- redist.flip(adj = adj, total_pop = iowa$pop, init_plan = iowa$cd_2010,
                    nsims = 100, pop_tol = 0.05)

## -----------------------------------------------------------------------------
class(sims)

## -----------------------------------------------------------------------------
dim(sims$plans)

## -----------------------------------------------------------------------------
redist.plot.map(shp = iowa, plan = sims$plans[, 100])

## -----------------------------------------------------------------------------
sims_comp <- redist.flip(adj = adj, total_pop = iowa$pop, init_plan = iowa$cd_2010,
                         nsims = 100, pop_tol = 0.05,
                         constraint = 'compact', constraintweights = 0.4,
                         compactness_metric = 'edges-removed')

## -----------------------------------------------------------------------------
redist.plot.map(shp = iowa, plan = sims_comp$plans[, 100])

## -----------------------------------------------------------------------------
set.seed(1)
nchains <- 4
nsims <- 100

## -----------------------------------------------------------------------------
flip_chains <- lapply(1:nchains, function(x){
  redist.flip(adj = adj,  total_pop = iowa$pop, pop_tol = 0.05,
              nsims = nsims, init_plan = 'rsg', ndists = 4,
              constraint = 'compact', constraintweights = 0.4,
              compactness_metric = 'edges-removed', verbose = FALSE)
})

## ---- eval=F------------------------------------------------------------------
#  mcmc_chains <- parallel::mclapply(1:nchains, function(x){
#    redist.flip(adjobj = adjlist,  popvec = fl25$pop,
#                nsims = nsims, ndists = 3)
#  }, mc.set.seed = 1, mc.cores = parallel::detectCores())

## -----------------------------------------------------------------------------
data(iowa)

## -----------------------------------------------------------------------------
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol=0.01)

## ---- include = FALSE---------------------------------------------------------
set.seed(2)

## ----sims---------------------------------------------------------------------
tidy_sims <- iowa_map %>% redist_flip(nsims = 100)

## -----------------------------------------------------------------------------
cons <- flip_constraints_helper(iowa_map, constraint = NULL)

## -----------------------------------------------------------------------------
tidy_sims_no_comp <- iowa_map %>% redist_flip(nsims = 100, constraints = cons)

## -----------------------------------------------------------------------------
cons <- flip_constraints_helper(map = iowa_map, constraint = c('compact', 'similarity'),
                                constraintweight = c(0.6, 1), init_plan = cd_2010)

## -----------------------------------------------------------------------------
tidy_sims_sq_comp <- iowa_map %>% redist_flip(nsims = 100, constraints = cons)

## -----------------------------------------------------------------------------
class(tidy_sims)

## -----------------------------------------------------------------------------
plans <- get_plans_matrix(tidy_sims)

## -----------------------------------------------------------------------------
tidy_sims <- tidy_sims %>% 
  mutate(competitiveness = competitiveness(map = iowa_map, rvote = rep_08, dvote = dem_08))

## -----------------------------------------------------------------------------
tidy_sims %>% 
  ggplot(aes(x = competitiveness)) +
  geom_density() + 
  theme_bw()

## -----------------------------------------------------------------------------
seg <- redist.segcalc(plans = get_plans_matrix(tidy_sims), 
                      group_pop = iowa_map$rep_08,
                      total_pop = iowa_map$pop)

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "autocorr")

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "densplot")

## -----------------------------------------------------------------------------
redist.diagplot(seg, plot = "mean")

## -----------------------------------------------------------------------------
seg_chains <- lapply(1:nchains, 
                     function(i){redist.segcalc(plans = flip_chains[[i]], 
                                                group_pop = iowa_map$rep_08,
                                                total_pop = iowa_map$pop)})

## -----------------------------------------------------------------------------
redist.diagplot(sumstat = seg_chains, plot = "trace")

## -----------------------------------------------------------------------------
redist.diagplot(sumstat = seg_chains, plot = 'gelmanrubin')

## -----------------------------------------------------------------------------
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.02, total_pop = pop)

cons <- flip_constraints_helper(map = iowa_map, constraint = c('compact', 'population'),
                                constraintweight = c(0.5, 100))
sims <- redist_flip(map = iowa_map,  nsims = 100)

## -----------------------------------------------------------------------------
mean(sims$mhdecisions, na.rm = TRUE)

## -----------------------------------------------------------------------------
sims_new <- redist_flip(map = iowa_map, nsims = 100, pop_tol = 0.02,
                        constraints = cons, eprob = 0.10, lambda = 2, verbose = FALSE)
mean(sims_new$mhdecisions, na.rm = TRUE)

## -----------------------------------------------------------------------------
dists <- redist.distances(plans = get_plans_matrix(sims))$Hamming
dists_new <- redist.distances(plans = get_plans_matrix(sims_new))$Hamming
adj_dists <- rep(NA_integer_, 100)
adj_dists_new <- rep(NA_integer_, 100)
for(i in 1:100){
  adj_dists[i] <- dists[i, i + 1]
  adj_dists_new[i] <- dists_new[i, i + 1]
}
tibble(Hamming = c(adj_dists, adj_dists_new), 
       `lambda/eprob` = c(rep('0/0.05', 100), rep('2/0.10', 100))) %>% 
  ggplot() + 
  geom_density(aes(x = Hamming, color = `lambda/eprob`)) + 
  theme_bw()

## -----------------------------------------------------------------------------
sims <- sims %>% mutate(par = plan_parity(map = iowa_map))

## -----------------------------------------------------------------------------
sims <- sims %>% filter(par <= 0.01)

## -----------------------------------------------------------------------------
cons <- flip_constraints_helper(map = iowa_map,
                                constraint = c('compact', 'population', 'partisan', 'countsplit'),
                                constraintweight = c(0.25, 50, 10, 10),
                                counties = name, 
                                rvote = iowa$rep_08, dvote = iowa$dem_08)

## -----------------------------------------------------------------------------
sims <- redist_flip(iowa_map, 100, constraints = cons)

## -----------------------------------------------------------------------------
cons <- flip_constraints_helper(map = iowa_map,
                                constraint = c('compact', 'population', 'partisan', 'countsplit'),
                                constraintweight = c(1.5, 100, 40, 20),
                                counties = name, 
                                rvote = iowa$rep_08, dvote = iowa$dem_08)
sims <- redist_flip(iowa_map, 100, constraints = cons)

## -----------------------------------------------------------------------------
summary(sims$constraint_compact, na.rm = TRUE)

## -----------------------------------------------------------------------------
summary(sims$constraint_population, na.rm = TRUE)

