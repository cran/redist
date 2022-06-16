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
fl25$plan <- plan
fl_map <- redist_map(fl25, existing_plan = plan, pop_tol = 0.2, total_pop = pop)
constr <- redist_constr(fl_map) %>% 
    add_constr_edges_rem(0.02)
set.seed(1)
sims <- redist_flip(map = fl_map, nsims = 6, constraints = constr)

## ---- echo = FALSE------------------------------------------------------------
redist.plot.map(shp = fl25, plan = plan) + redist.plot.adj(shp = fl25, plan = plan, plot_shp = FALSE) 

## ---- echo = FALSE------------------------------------------------------------
redist.plot.adj(shp = fl25, plan = plan, plot_shp = FALSE, title = 'Initialization') +
  redist.plot.adj(shp = fl25, plan = get_plans_matrix(sims)[, 1], plot_shp = FALSE, title = 'First Iteration') +
  redist.plot.adj(shp = fl25, plan = get_plans_matrix(sims)[, 2], plot_shp = FALSE, title = 'Second Iteration') +
  redist.plot.adj(shp = fl25, plan = get_plans_matrix(sims)[, 3], plot_shp = FALSE, title = 'Third Iteration') +
  redist.plot.adj(shp = fl25, plan = get_plans_matrix(sims)[, 4], plot_shp = FALSE, title = 'Fourth Iteration') +
  redist.plot.adj(shp = fl25, plan = get_plans_matrix(sims)[, 5], plot_shp = FALSE, title = 'Fifth Iteration')

## -----------------------------------------------------------------------------
data(iowa)
redist.plot.map(iowa, plan = cd_2010)
map_ia <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.05)

## ----seed, include = FALSE----------------------------------------------------
# set the seed so that within the vignette the results don't change
set.seed(1)

## -----------------------------------------------------------------------------
sims <- redist_flip(map_ia, nsims = 100)

## -----------------------------------------------------------------------------
class(sims)

## -----------------------------------------------------------------------------
dim(get_plans_matrix(sims))

## -----------------------------------------------------------------------------
redist.plot.map(shp = iowa, plan = get_plans_matrix(sims)[, 100])

## -----------------------------------------------------------------------------
constr <- redist_constr(map_ia) %>% add_constr_edges_rem(0.4)

sims_comp <- redist_flip(map_ia, nsims = 100, constraints = constr)

## -----------------------------------------------------------------------------
redist.plot.map(shp = iowa, plan = get_plans_matrix(sims_comp)[, 100])

## -----------------------------------------------------------------------------
set.seed(1)
nchains <- 4
nsims <- 100

## -----------------------------------------------------------------------------
constr <- redist_constr(map_ia) %>% add_constr_edges_rem(0.4)
map_ia <- redist_map(iowa, ndists = 4, pop_tol = 0.05)
flip_chains <- lapply(1:nchains, function(x){
  redist_flip(map_ia, nsims = nsims,
              constraints = constr, verbose = FALSE)
})

## ---- eval=F------------------------------------------------------------------
#  mcmc_chains <- parallel::mclapply(1:nchains, function(x){
#    redist.flip(fl_map, nsims)
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
cons <- redist_constr(iowa_map)

## -----------------------------------------------------------------------------
tidy_sims_no_comp <- iowa_map %>% redist_flip(nsims = 100, constraints = cons)

## -----------------------------------------------------------------------------
class(tidy_sims)

## -----------------------------------------------------------------------------
plans <- get_plans_matrix(tidy_sims)

## -----------------------------------------------------------------------------
tidy_sims <- tidy_sims %>% 
  mutate(competitiveness = rep(competitiveness(map = iowa_map, rvote = rep_08, dvote = dem_08), each = 4))

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
                     function(i){redist.segcalc(plans = get_plans_matrix(flip_chains[[i]]), 
                                                group_pop = iowa_map$rep_08,
                                                total_pop = iowa_map$pop)})

## -----------------------------------------------------------------------------
redist.diagplot(sumstat = seg_chains, plot = "trace")

## -----------------------------------------------------------------------------
redist.diagplot(sumstat = seg_chains, plot = 'gelmanrubin')

## -----------------------------------------------------------------------------
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.02, total_pop = pop)

cons <- redist_constr(iowa_map) %>% 
    add_constr_edges_rem(0.5) %>% 
    add_constr_pop_dev(100)

sims <- redist_flip(map = iowa_map,  nsims = 100)

## -----------------------------------------------------------------------------
mean(sims$mhdecisions, na.rm = TRUE)

## -----------------------------------------------------------------------------
sims_new <- redist_flip(map = iowa_map, nsims = 100, constraints = cons, 
                        eprob = 0.10, lambda = 2, verbose = FALSE)
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
cons <- redist_constr(iowa_map) %>% 
    add_constr_edges_rem(0.25) %>% 
    add_constr_pop_dev(50) %>% 
    add_constr_compet(10, rvote = rep_08, dvote = dem_08) %>% 
    add_constr_splits(10, admin = region)

## -----------------------------------------------------------------------------
sims <- redist_flip(iowa_map, 100, constraints = cons)

## -----------------------------------------------------------------------------
cons <- cons <- redist_constr(iowa_map) %>% 
    add_constr_edges_rem(1.5) %>% 
    add_constr_pop_dev(100) %>% 
    add_constr_compet(40, rvote = rep_08, dvote = dem_08) %>% 
    add_constr_splits(20, admin = region)

sims <- redist_flip(iowa_map, 100, constraints = cons)

## -----------------------------------------------------------------------------
summary(sims$constraint_edges_removed, na.rm = TRUE)

## -----------------------------------------------------------------------------
summary(sims$constraint_pop_dev, na.rm = TRUE)

