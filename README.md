# SynAnthrop: an R package to analyse species distribution along anthropisation gradient 
<img src="man/figures/SynAnthrop_logo.png" align="right" alt="SynAnthrop logo" width="200" />


## Description
`SynAnthrop` is a R package developed to assess the sensitivity of species and communities to anthropogenization based on occurrence data, via the calculation of the Species Synanthropy Index (SSI).

`Synanthrop` allows to describe the ecological affinities of species and groups in a simple, reproducible, multi-scale and less subjective way than with expert assessments. The SSI scores can then be used to identify species, communities, and related habitats that require significant conservation and restoration effort. They are complementary to other existing priorization indices such as red lists, rarity or specialization degrees.

The analysis is based on the comparison of the observed distribution of species along a gradient of anthropogenization within a territory, and a distribution that would be expected if anthropisation had no effect on this species distribution (null distribution). Null distributions are generated per species based on their range and survey intensity per site.


## Workflow 
<img src="man/figures/Synanthrop_workflow.png" />


## Installation

### If you wish to visualize the SSI tutorial only

2 options:  
* Access the "Articles" tab of <a href= "https://lomorel.github.io/SynAnthrop">the SynAnthrop GitHub site</a>  
* Download the PDF tutorial in the vignettes folder of <a href= "https://github.com/lomorel/SynAnthrop/">the SynAnthrop GitHub repository</a>  


### If you wish to install the package and execute the tutorial
In the R terminal tab, `cd` to the folder where you want to store the project, and type:

```
git clone https://github.com/lomorel/SynAnthrop.git
```

Access the content of the package by opening the `SynAnthrop.Rproj` file in Rstudio. Install `SynAnthrop` by typing in the R console:
```
devtools::install()
```
Finally, open the tutorial file (.Rmd) in the `vignettes` folder, and execute the code chunks.


### If you wish to use the package functions only, without using the tutorial
Install the package by typing the following line of code in the R console:

``` 
devtools::install_github("lomorel/SynAnthrop", dependencies = TRUE)
library(SynAnthrop)
``` 

The `ssi` function is now ready to use, and the help file is accessible by typing `?ssi`




## Citation
Package and tutorial written by Loïs Morel, Lab. DECOD, Institut Agro, Rennes, France.

Citations:  
- Morel L. 2023. SynAnthrop: Species distribution and sensitivity to anthropisation, R package version 0.1.1, https://github.com/lomorel/SynAnthrop  
  