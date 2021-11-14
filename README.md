# sc2marker

sc2marker is a toolkit for identification of cell specific markers for clusters of cells in a single cell RNA sequencing data set. In particular, sc2marker focus on markers that have valid antibodies to be used for either imaging or flow analyis.

## Install

You can install sc2marker with the following commands: 

```{r}
install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/sc2marker", build_vignettes = TRUE)
require(sc2marker)
require(Seurat)
```

To run sc2marker you need to execute the following command, providing a clustered single  cell data sets (as Seurat object), the cell type of interest and the antibody databased (IHC, ICC or Flow). 

```{r}
t.markers <- Detect_single_marker(mca.spleen, id = "T cell(Spleen)", category = "Flow")
```

You can display the results as a table with the command:

```{r}
get_antibody(t.markers)
```

and you can generate ridge plot with the following command:

```{r}
plot_ridge(mca.spleen, id = "T cell(Spleen)", genes = t.markers[1:9,]$gene, ncol = 3, assay = "RNA")
```

To calculate markers for all cell clusters, you can do by following command:

```{r}
all.markers <- Detect_single_marker_all(mca.spleen, category = "Flow")
```

To Check T cell markers from results of all clusters, and get the antibody information, you can do following.

```{r}
t.markers <- all.markers[["T cell(Spleen)"]]
get_antibody(t.markers)
```

To automatically generate the sc2marker report of all cell clusters, you can run following command:

```{r}
generate_report(mca.spleen, all.markers, fpath = ".")
```



