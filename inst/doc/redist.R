## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
set.seed(5118)

## ----setup, message=FALSE-----------------------------------------------------
library(redist)
library(dplyr)
library(ggplot2)

## -----------------------------------------------------------------------------
data(iowa)
print(iowa)

## -----------------------------------------------------------------------------
iowa_map = redist_map(iowa, existing_plan=cd_2010, pop_tol=0.01, total_pop = pop)
print(iowa_map)

## ----iowa-adj, fig.width=8----------------------------------------------------
plot(iowa_map, adj=T) + plot(iowa_map)

## ----iowa-chloro, message=FALSE-----------------------------------------------
areas = as.numeric(units::set_units(sf::st_area(iowa_map$geometry), mi^2))
plot(iowa_map, fill = pop / areas) + 
    scale_fill_viridis_c(name="Population density (people / sq. mi)", 
                         trans="sqrt")
plot(iowa_map, fill = dem_08 / tot_08) +
    scale_fill_gradient2(name="Pct. Democratic '08",  midpoint=0.5)
plot(iowa_map, fill = wvap / vap, by_distr = TRUE)

## ----message=FALSE------------------------------------------------------------
iowa_plans = redist_smc(iowa_map, 500, compactness=1, runs=2)

## -----------------------------------------------------------------------------
print(iowa_plans)

## ----ia-sim-plans-------------------------------------------------------------
redist.plot.plans(iowa_plans, draws=1:6, shp=iowa_map)

## -----------------------------------------------------------------------------
iowa_plans = match_numbers(iowa_plans, iowa_map$cd_2010)
print(iowa_plans)

## ----message=F----------------------------------------------------------------
county_perims = prep_perims(iowa_map, iowa_map$adj)

iowa_plans = iowa_plans %>%
    mutate(pop_dev = abs(total_pop / get_target(iowa_map) - 1),
           comp = comp_polsby(pl(), iowa_map, perim_df=county_perims),
           pct_min = group_frac(iowa_map, vap - wvap, vap),
           pct_dem = group_frac(iowa_map, dem_08, dem_08 + rep_08))
print(iowa_plans)

## -----------------------------------------------------------------------------
summary(iowa_plans)

## -----------------------------------------------------------------------------
plan_sum = group_by(iowa_plans, draw) %>%
    summarize(max_dev = max(pop_dev),
              avg_comp = mean(comp),
              max_pct_min = max(pct_min),
              dem_distr = sum(pct_dem > 0.5))
print(plan_sum)

## ----dev-comp-plot------------------------------------------------------------
library(patchwork)

hist(plan_sum, max_dev) + hist(iowa_plans, comp) +
    plot_layout(guides="collect")

## ----signature----------------------------------------------------------------
plot(iowa_plans, pct_dem, sort=FALSE, size=0.5)

## ----scatter------------------------------------------------------------------
pal = scales::viridis_pal()(5)[-1]
redist.plot.scatter(iowa_plans, pct_min, pct_dem, 
                    color=pal[subset_sampled(iowa_plans)$district]) +
    scale_color_manual(values="black")

