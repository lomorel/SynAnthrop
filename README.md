
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SynAnthrop : a R package to analyse species distribution along anthropisation gradient <img src="./Figure/SynAnthrop_logo.png" align="right" alt="" width="200" />

SynAnthrop is a R package developed to assess the sensibility of species
and communities to anthropisation from occurrence data.

The method allows to describe the ecological affinities of species and
groups in a simple, reproducible, multi-scale and less subjective way
than with expert assessments.

The principle of the analysis is based on the comparison of the observed
distribution of species along a gradient of anthropisation within a
territory, with a distribution that would be expected if anthropisation
had no effect on this distribution (null distribution). This null
distribution is defined by species, taking into account its distribution
(points are randomly generated within its range) and taking into account
the survey intensity per site.

## Workflow <img src="./Figure/Synanthrop_workflow.png" />

Two types of data are needed as input (1): georeferenced occurrence data
and a map describing the gradient to be studied (from which we extract
values at the desired resolution). After a selection of data by species
(2), a random selection of sites in the range and according to the
sampling effort is made (n = observed occurrences). The gradient effect
(of anthropisation in our case) is measured by evaluating the effect
size between these 2 distributions (3). The operation is repeated 500
times (4). Synanthropy scores (from 1 to 10) are assigned to each
species from the average of the differences between the null and
observed distributions (5).

These scores can then be used to identify species, communities, and
related habitats that require significant conservation and restoration
effort. They are complementary to other existing prioritisation indices
such as red lists, rarity or specialization degrees.

## Dependencies

### Loading packages and function

``` r
Packages <- c("tidyverse", "raster","scales","sf","ks","tabularaster", "terra", "tidyterra",
              "CGPfunctions", "RColorBrewer")
# install.packages(Packages) # if needed

lapply(Packages, library, character.only = TRUE) # to load

source("./R/Species_Synanthropy_Index_function.R")
# devtools::install_github("/lomorel/SynAnthrop")
```

### Data

Two types of data are required to run the Species Synanthropy Index
(SSI): (i) a database of species occurrences, with XY coordinates and
corresponding sampling dates, and (ii) a raster describing the spatial
gradient of anthropisation of the region to be analyse.

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

Note that data were filter to containt only recent years (2010-2021) in
order to be coherent with information compile in used maps.

``` r
ras_raw <- raster("./Data/CartNat_Bzh.tif")#raster
```

The map used here is the French Naturalness Map, developed by [Guetté et
al.(2021)](https://uicn.fr/CartNat/CartNat_Donnees/Note_technique_m%C3%A9thodologique/Projet%20CARTNAT_note%20technique_2021.pdf).

``` r
rast_to_plot <- rast(ras_raw)

sp_map <- st_as_sf(sp_by_occ_raw, coords = c("X", "Y"), crs = 2154)
sp_map$Year <- as.integer(sp_map$Year)

ggplot() + 
  geom_spatraster(data = rast_to_plot) +
  geom_sf(data = sp_map, aes(color = Year)) +
  scale_colour_gradient(low = "orange", high = "red") +
  scale_fill_whitebox_c(palette = "muted", direction=-1)
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## The Species Synanthropy Index function

### An example with amphibian populations in western France

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

### Example and results

``` r
ssi_results <- ssi(r = ras_raw, x = sp_by_occ_raw, resolution = c(100, 200) , sim = 50, threshold = 30)
```

The SSI function produce three main data.frame :

- `[[1]]` the first data.frame is a short summary table (one line per
  species) with the mean of all effect size and the corresponding index
  (range between 1 to 10) for each resolution. The number of runs used
  to calculate the mean effect size are also specified.

``` r
head(ssi_results[[1]])
```

    ##                  Species         mean nRun Index Resolution
    ## 1    Alytes obstetricans  0.007796773   50     2        100
    ## 2          Bufo spinosus  0.057358688   50     1        100
    ## 3      Epidalea calamita -0.502315252   50    10        100
    ## 4           Hyla arborea  0.012370339   50     2        100
    ## 5 Ichthyosaura alpestris -0.034269308   50     3        100
    ## 6 Lissotriton helveticus -0.098954661   50     4        100

- `[[2]]` the second data.frame compile all the raw results, i.e all the
  effect size assessed per run, with corresponding information provided
  by the function cohen_d (rstatix package) ; n1 and n2 correspond to
  the number of occurrence compare (n1 for the null distribution and n2
  for observed data).

``` r
head(ssi_results[[2]])
```

    ##                  Species   .y. group1 group2      effsize   n1   n2  magnitude
    ## 1    Alytes obstetricans value   Null    Obs -0.053743289  220  220 negligible
    ## 2          Bufo spinosus value   Null    Obs  0.057136020 2012 2011 negligible
    ## 3      Epidalea calamita value   Null    Obs -0.584828096   95   95   moderate
    ## 4           Hyla arborea value   Null    Obs -0.001765089  689  689 negligible
    ## 5 Ichthyosaura alpestris value   Null    Obs -0.048397667  237  237 negligible
    ## 6 Lissotriton helveticus value   Null    Obs -0.102106841 1197 1197 negligible
    ##   Run Resolution
    ## 1   1        100
    ## 2   1        100
    ## 3   1        100
    ## 4   1        100
    ## 5   1        100
    ## 6   1        100

- `[[3]]` the third data.frame archives all the occurrence randomly
  drawn to assess scores.

``` r
head(ssi_results[[3]])
```

    ##    Cell             Species        x       y variable Resolution
    ## 1 10220 Triturus marmoratus 304047.6 6750984     Null        100
    ## 2  7451 Triturus marmoratus 202047.6 6786984     Null        100
    ## 3  6716 Triturus marmoratus 242047.6 6796984     Null        100
    ## 4  5058 Triturus marmoratus 248047.6 6818984     Null        100
    ## 5 11568 Triturus marmoratus 282047.6 6732984     Null        100
    ## 6  3675 Triturus marmoratus 200047.6 6836984     Null        100

## Visualise and analyse the SSI results

### Score distribution within the studied taxa

``` r
sp_ssi <- ssi_results [[1]]
sp_ssi$Index <- as.integer(sp_ssi$Index)

mean_ssi_by_resolution <-   data.frame(sp_ssi %>%
                                      group_by(Species, Resolution) %>%
                                      summarise(Index = mean(Index), n = n()))

effsize_res <- ssi_results [[2]]
```

Here we selected only data for the resolution 200.

``` r
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

![](README_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

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
                Title = "Synanthropy scores for amphibian species in western France",
                SubTitle = NULL,
                Caption = NULL)
```

![](README_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

### Distribution map

Plotting these maps allows to visualise the locations of observations
and of randomly selected sites.

``` r
sub_distri <- subset(ssi_results[[3]], Resolution == "200")
sub_distri_obs <- subset(sub_distri, variable  == "Obs")
distri <- rbind(sub_distri_obs, sample_n(subset(sub_distri, variable  == "Null" ), nrow(sub_distri_obs)))
distri <- sf_transform_xy (distri, 4326, 2154)

ggplot() + 
  geom_point(data = subset(distri, variable == "Null"), aes(x = x, y = y, color = variable), 
             size = 1.5, colour = "#ff6600",alpha = 1/2) +
  geom_point(data = subset(distri, variable == "Obs"), aes(x = x, y = y, color = variable), 
             size = 1.5, colour = "#660099",alpha = 1/2) +
  theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
  facet_wrap(. ~ Species) +
  theme_void() +
  theme(legend.position="none")
```

![](README_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

## Application : an example with amphibian communities survey in ponds

The scores may then be used at the community scale to estimate the
overall sensitivity of assemblages to land use change. These average
scores are complementary to anthropisation maps whose information must
be considered as potential (with cases where the maps present values \>
to the scores when sensitive species are not detected and cases where
the maps present values \< to the scores when in reality some
populations of sensitive species have maintained themselves locally).

Here we test the scores with data of amphibian communities surveyed in
western France (Rennes) between 2000 and 2010

### Data cooking

``` r
ex_com <- read.table("./Data/Amphibian_commnity_MNIE.csv", sep=";", h=T)

#Subset to the finest resolution
sp_ssi100 <- subset(sp_ssi, Resolution == "100")

#Subset to the finest resolution
ex_com_df <- merge(ex_com, sp_ssi100[,c("Species","Index")], by = "Species")

ex_com_df1 <- data.frame(ex_com_df %>%
                            group_by(Site, Years, X, Y) %>%
                            summarise(Richness = n_distinct(Species)))

ex_com_df2 <-   data.frame(ex_com_df %>%
                                      group_by(Site, Years) %>%
                                      summarise(Mean_Index = round(mean(Index), 2)))

ex_com_df3 <-   data.frame(ex_com_df %>%
                                      group_by(Site, Years) %>%
                                      summarise(CWM = round(weighted.mean(Index, Abundance), 2)))

ex_com_df <- merge(ex_com_df1, ex_com_df2, by = c("Site", "Years")) 
ex_com_df <- merge(ex_com_df, ex_com_df3, by = c("Site", "Years"))
```

Here, we see the summary table with species richness, unweighted
community index and the community (abundance-)weighted mean (CWM).

``` r
head(ex_com_df)
```

    ##       Site Years      X       Y Richness Mean_Index  CWM
    ## 1   Site 1  2010 335076 6788565        3       3.33 3.20
    ## 2  Site 10  2010 340610 6790480        5       3.40 3.40
    ## 3 Site 100  2010 351549 6778613        5       3.40 3.67
    ## 4 Site 101  2010 351550 6778573        3       3.33 3.33
    ## 5 Site 102  2010 351577 6780908        3       2.67 3.54
    ## 6 Site 103  2010 351641 6778430        4       2.25 2.84

### Map the results

``` r
ex_com <- read.table("./Data/Amphibian_commnity_MNIE.csv", sep=";", h=T)
ex_com_df_map <- st_as_sf(ex_com_df, coords = c("X", "Y"), crs = 2154)


#Plot the map
ex_com_df_map2 <- filter(ex_com_df_map, Richness > 1)

rast_to_plot2 <- crop(rast_to_plot, ex_com_df_map2)

ggplot() + 
  geom_spatraster(data = rast_to_plot2, alpha = 0.4) +
  geom_sf(data = ex_com_df_map2, aes(color = CWM, size = Richness)) +
  scale_colour_gradient(low = "yellow", high = "darkgreen") +
  scale_fill_whitebox_c(palette = "muted", direction=-1) +
  ggtitle("Amphibian communities in ponds ecosystems") 
```

![](README_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

Here we identify the ponds that have less synanthropic assemblages and
are therefore potentially more sensitive to urbanisation and intensified
land-use. We also distinguish the ponds presenting rather synanthropic
assemblages. It is interesting to note that these results contrast in
part with the map, illustrating the complementarity of the two
approaches.

## Credits

Package and tutorial written by Loïs Morel, Lab. DECOD, Institut Agro,
Rennes, France.

Citations:

- Morel L. 2023. SynAnthrop: Species distribution and sensitivity to
  anthropisation, R package version 0.1.1,
  <https://github.com/lomorel/SynAnthrop>
