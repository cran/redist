## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(5118)

## ----setup, message=F---------------------------------------------------------
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

library(redist)

data_url = "https://github.com/alarm-redist/redist-data/raw/main/data/king_county.rds"
data_path = tempfile()
download.file(data_url, data_path)
king_shp = readRDS(data_path)
print(king_shp)

## ----include=F----------------------------------------------------------------
n_city = length(unique(king_shp$city)) - 1
tot_area = as.numeric(sum(sf::st_area(king_shp)))
tot_pop = sum(king_shp$pop)
unincorp_area = as.numeric(sum(sf::st_area(king_shp[king_shp$city == "UNINCORP", ])))
unincorp_pop = sum(king_shp$pop[king_shp$city == "UNINCORP"])

## ----water-plot, fig.width=8, echo=F------------------------------------------
areas = as.numeric(units::set_units(sf::st_area(king_shp), mi^2))
pop_plot = ggplot(king_shp, aes(fill=pop/areas)) +
    geom_sf(size=0) +
    scale_fill_viridis_c(trans="sqrt", labels=comma, limits=c(0, 25e3), oob=squish) +
    labs(title="Population Density") +
    theme_void() + theme(legend.position="bottom")
water_plot = ggplot(king_shp, aes(fill=pct_water)) +
    geom_sf(size=0) +
    scale_fill_viridis_c(labels=percent) +
    labs(title="Water") +
    theme_void() + theme(legend.position="bottom")
water_plot + pop_plot

## ----city-distr-plot, echo=F--------------------------------------------------
cities = summarize(group_by(king_shp, city, distr), .groups="drop")
districts = summarize(group_by(king_shp, distr))
ggplot(cities) +
    geom_sf(aes(fill=city, alpha=city!="UNINCORP"), color="#444444") +
    geom_sf(color="black", size=1, data=districts, fill="#00000000") +
    geom_sf_text(aes(label=distr), size=10, fontface="bold", color="#000000aa",
                 data=districts) +
    scale_alpha_manual(values=c(0, 0.8)) +
    guides(fill="none", alpha="none") +
    labs(title="Cities with Council District Overlays") + 
    theme_void()

## -----------------------------------------------------------------------------
existing_parity = redist.parity(king_shp$distr, king_shp$pop)
king = redist_map(king_shp, existing_plan=distr, pop_tol=existing_parity)
print(king)

## ----king-adj, fig.height=6---------------------------------------------------
plot(king, adj=T, centroids=F, zoom_to=(city == "SEA"))

## -----------------------------------------------------------------------------
filter(king, distr %in% c(2, 4, 8))

## ----king-water---------------------------------------------------------------
plot(filter(king, pct_water >= 0.99, pop == 0)) + geom_sf_text(aes(label=id))

## ----king-land----------------------------------------------------------------
water_prec = filter(king, pct_water >= 0.99, pop == 0) %>% pull(id)
water_prec = setdiff(water_prec, c("WVPS34", "WVSP34"))
king_land = filter(king, !(id %in% water_prec))
plot(king_land)

## ----seattle-land-adj, fig.height=6-------------------------------------------
plot(king_land, adj=T, centroids=F, zoom_to=(city == "SEA"))

## -----------------------------------------------------------------------------
merge_by(king_land, city)

## ----merged-city--------------------------------------------------------------
king_merged = merge_by(king_land, city, drop_geom=FALSE)
plot(king_merged, adj=T)

## -----------------------------------------------------------------------------
cat(redist.splits(king_land$distr, king_land$city), "split cities\n")

king_land %>%
    mutate(is_unsplit = !is_county_split(distr, city)) %>%
plot(is_unsplit)

king_unsplit = king_land %>%
    mutate(unsplit_id = freeze(!is_county_split(distr, city))) %>%
    merge_by(unsplit_id, city, collapse_chr=FALSE)
print(king_unsplit)

## ----unsplit-plan, fig.width=8------------------------------------------------
plans = redist_smc(king_unsplit, 50, silent=T)
print(plans)
print(pullback(plans))
redist.plot.plans(pullback(plans), draws=1:4, geom=king_land)

## -----------------------------------------------------------------------------
pop_inside_cores = function(boundary) {
    king_land %>%
        mutate(core = make_cores(boundary=boundary)) %>%
        group_by(core) %>%
        filter(n() > 2) %>% # filter to cores only
        pull(pop) %>% sum
}
pop_inside_cores(1) / sum(king_land$pop)
pop_inside_cores(2) / sum(king_land$pop)

## ----cores--------------------------------------------------------------------
king_cores = king_land %>%
    mutate(core = make_cores(boundary=1)) %>%
    merge_by(core, drop_geom=F)
plot(king_cores)

## ----core-plans---------------------------------------------------------------
plans = redist_smc(king_cores, 50, silent=T)
redist.plot.plans(plans, draws=1:4, geom=king_cores)

