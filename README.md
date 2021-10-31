# sc2marker
sc2marker is a toolkit for identifying specific markers for interested cell groups. 

## Install

Install scMarkerDetect with the simple comands below:

```{r}
install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/sc2marker", build_vignettes = TRUE)
require(scMarkerDetect)
```

Identify single markers for specific cell groups

```{r}
data(pbmc_small)
Idents(pbmc_small) <- "groups"
markers.results <- Detect_single_marker(pbmc_small, "g2", do.fast = F)
```


