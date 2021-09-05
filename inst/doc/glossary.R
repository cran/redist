## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(redist)
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop)

## -----------------------------------------------------------------------------
# Standard eval - 
adj <- redist.adjacency(shp = iowa)
# tidy eval - 
adj <- get_adj(iowa_map)

head(adj)

## -----------------------------------------------------------------------------
# Plot it!
redist.plot.adj(shp = iowa_map)

## -----------------------------------------------------------------------------
# Standard eval - 
ndists <- 4
# tidy eval - stored within redist_map object
attr(iowa_map, 'ndists')

## -----------------------------------------------------------------------------
nsims <- 100

## -----------------------------------------------------------------------------
# standard eval - 
pop_tol <- 0.01

# tidy eval - stored within redist_map object
# - getting
get_pop_tol(iowa_map)
# - setting
iowa_map <- set_pop_tol(iowa_map, pop_tol = 0.01)

## -----------------------------------------------------------------------------
sim <- redist.rsg(adj = adj, total_pop = iowa$pop, ndists = 4, pop_tol = 0.01)

head(sim$plan)

## -----------------------------------------------------------------------------
# standard eval -
sims <- redist.smc(adj = adj, total_pop = iowa$pop, nsims = 10, ndists = 4, silent = TRUE)
plans <- sims$plans

# tidy eval - 
sims <- redist_smc(map = iowa_map, nsims = 10, silent = TRUE)
plans <- get_plans_matrix(sims)

## -----------------------------------------------------------------------------
# standard eval - 
init_plan <- iowa$cd_2010

# tidy eval - stored within redist_map object
get_existing(iowa_map)

## -----------------------------------------------------------------------------
# standard eval
total_pop <- iowa$pop

# tidy eval - a column within the redist_map object tracked by attributes
iowa_map[[attr(iowa_map, 'pop_col')]]

## -----------------------------------------------------------------------------
iowa$white

## -----------------------------------------------------------------------------
# tidy eval - stored in redist_map object
attr(iowa_map, 'pop_bounds')

## -----------------------------------------------------------------------------


