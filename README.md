
true

<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynAnthrop <img src="./Figure/SynAnthrop_logo.png" align="right" alt="" width="200" />

SynAnthrop is a R package developed to assess the sensibility of species
and communities to anthropisation from occurrence data.

## Workflow <img src="./Figure/Synanthrop_workflow.png" />

Principle: the analysis is based on the cross-referencing of
georeferenced occurrence data and a map describing the gradient to be
studied (from which values are extracted at the desired resolution).
After a selection of data by species (2), a random selection of sites in
the range and according to the sampling effort is made (n = observed
occurrences). The effect (in this case of anthropisation) is measured by
evaluating the effect size between these 2 distributions (3). The
operation is repeated 500 times (4). Synanthropy scores (from 1 to 10)
are assigned to each species from the average of the differences between
the null and observed distributions (5).

## Dependencies

### Loading packages and function

``` r
Packages <- c("tidyverse", "raster","scales","sf","ks","tabularaster", "terra", "tidyterra", "CGPfunctions")
# install.packages(Packages) # if needed

lapply(Packages, library, character.only = TRUE) # to load

source("./R/Species_Synanthropisation_Index_function.R")
# devtools::install_github("/lomorel/SynAnthrop")
```

### Data

Two types of data are required to run the Species Synanthropisation
Index (SSI): (i) a database of species occurrences, with XY coordinates
and corresponding sampling dates, and (ii) a raster describing the
spatial gradient of anthropisation of the region to be analyse.

``` r
sp_by_occ_raw <- read.table("./Data/amphibian_all_2154_to_R.csv", sep=";", h=T)
head(sp_by_occ_raw)
```

    ##                 Species Year Abundance      X       Y
    ## 1         Bufo spinosus 2015         1 259519 6807029
    ## 2       Pelophylax spp. 2015         1 259519 6807029
    ## 3 Salamandra salamandra 2015         1 259519 6807029
    ## 4         Bufo spinosus 2015         1 259519 6807029
    ## 5 Salamandra salamandra 2015         1 259519 6807029
    ## 6       Pelophylax spp. 2015         1 261357 6806526

``` r
ras_raw <- raster("./Data/CartNat_Bzh.tif")#raster
```

The map used here is the French Naturalness Map, developed by Guetté et
al.(2021)
<https://uicn.fr/CartNat/CartNat_Donnees/Note_technique_m%C3%A9thodologique/Projet%20CARTNAT_note%20technique_2021.pdf>).

``` r
# create SpatRaster (can also do `x <- rast(f)`
rast_to_plot <- rast(ras_raw)

sp_map <- st_as_sf(sp_by_occ_raw, coords = c("X", "Y"), crs = 2154)
sp_map$Year <- as.integer(sp_map$Year)

ggplot() + 
  geom_spatraster(data = rast_to_plot) +
  geom_sf(data = sp_map, aes(color = Year)) +
  scale_fill_whitebox_c(palette = "muted")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Running the SSI function

### Usage

The function was designed to calculate synanthropy scores automatically
from a data set (for a taxon, a territory and a period). Several
resolutions can be evaluated simultaneously.

### Arguments

- `- resolution`: his argument allows to test several resolutions. To do
  this, specify for each resolution you want to test (the value is the
  number of raster cells aggregated to compile species occurrences).

- `- sim`: this argument corresponds to the number of simulations to run
  to model the null distribution of occurrences under the assumption
  that the variables used for the map (raster) do not influence the
  species distribution.

- `- threshold`: this argument allows to fix the threshold of species
  occurrences under it species will be not account.

### An example with amphibian populations in western France

``` r
ssi_results <- ssi(r = ras_raw, x = sp_by_occ_raw, resolution = 200 , sim = 2, threshold = 30)
```

The SSI function produce three main data.frame :

- `[[1]]` the first data.frame is a short summary table (one line per
  species) with the mean of all effect size and the corresponding index
  (range between 1 to 10) for each resolution. The number of runs used
  to calculate the mean effect size are also specified.

``` r
head(ssi_results[[1]])
```

- `[[2]]` the second data.frame compile all the raw results, i.e all the
  effect size assessed per run, with corresponding information provided
  by the function cohen_d (rstatix package) ; n1 and n2 correspond to
  the number of occurrence compare (n1 for the null distribution and n2
  for observed data).

``` r
head(ssi_results[[2]])
```

- `[[3]]` the third data.frame archives all the occurrence randomly
  drawn to assess scores.

``` r
head(ssi_results[[3]])
```

## Visualising SSI results

### Score distribution within the studied taxa

``` r
sp_ssi <- ssi_results [[1]]
sp_ssi$Index <- as.integer(sp_ssi$Index)

mean_ssi_by_resolution <-   data.frame(sp_ssi %>%
                                      group_by(Species, Resolution) %>%
                                      summarise(Index = mean(Index), n = n()))

effsize_res <- ssi_results [[2]]

# then scale by scale
sub_effsize_res <- subset(effsize_res, Resolution == "200")
mean_ssi_by_resolution <- subset(mean_ssi_by_resolution, Resolution == "200")

sub_effsize_res <- merge(sub_effsize_res, mean_ssi_by_resolution, by = "Species")
sub_effsize_res$Index <- as.factor(sub_effsize_res$Index)

ggplot(sub_effsize_res, aes(x = reorder(Species, -effsize), y = -effsize, fill = Index)) +
  geom_hline(yintercept = 0.0, color = "darkgrey", size=0.8, linetype="dashed") +
  geom_boxplot() + 
  coord_flip() +
  scale_fill_brewer(name="Score", palette = "RdYlGn") +
  ylab("Effect size") +
  theme(axis.title.y = element_blank()) +
  theme(legend.position = c(0.9, 0.2)) +
  theme_bw()
```

![](README_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

### Resolution comparison

``` r
mean_index_by_reso <- data.frame(sp_ssi %>%
                                      group_by(Species, Resolution) %>%
                                      summarise(Index = mean(Index), n = n()))

mean_index_by_reso$Resolution <- as.factor(mean_index_by_reso$Resolution)

newggslopegraph(dataframe = mean_index_by_reso,
                Resolution,
                Index,
                Grouping = Species,
                Title = "Amphibian",
                SubTitle = NULL,
                Caption = NULL)
```

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## Credits

Package and tutorial written by Loïs Morel, Lab. DECOD, Institut Agro,
Rennes, France.

Citations:

    Morel L. 2023. SynAnthrop: Species distribution and sensitivity to anthropisation, R package version 0.1.1,         https://github.com/lomorel/SynAnthrop