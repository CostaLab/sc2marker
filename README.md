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
t.markers <- Detect_single_marker(mca.spleen, id = "T cell(Spleen)", pseudo.count = 0.01)
```

Display Markers with antibody information
```{r}
get_antibody(t.markers)
```

Gene ridge plot

```{r}
plot_ridge(mca.spleen, id = "T cell(Spleen)", genes = t.markers[1:9,]$gene, ncol = 3, assay = "RNA")
```
