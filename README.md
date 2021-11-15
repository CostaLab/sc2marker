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

## Detect markers for one cell cluster

To run sc2marker you need to execute the following command, providing a clustered single  cell data sets (as Seurat object), the cell type of interest and the antibody databased (IHC, ICC or Flow). 

```{r}
nk.markers <- Detect_single_marker(mca.spleen, id = "NK cell(Spleen)", category = "Flow")
```

You can display the results as a table with the command:

```{r}
get_antibody(nk.markers)
```

and you can generate ridge plot with the following command:

```{r}
plot_ridge(mca.spleen, id = "NK cell(Spleen)", genes = nk.markers[1:9,]$gene, ncol = 3, assay = "RNA", aggr.other = F)
```

## calculate markers for all cell clusters

To calculate markers for all cell clusters, you can do by following command:

```{r}
all.markers <- Detect_single_marker_all(mca.spleen, category = "Flow")
```

To Check T cell markers from results of all clusters, and get the antibody information, you can do following.

```{r}
t.markers <- all.markers[["T cell(Spleen)"]]
get_antibody(t.markers)
```

## Generate report

To automatically generate the sc2marker report of all cell clusters, you can run following command:

```{r}
generate_report(mca.spleen, all.markers, fpath = ".", fname = "mca.spleen")
```

Or you can only analysis subset of cell clusters (B cell and NK cell) and generate the report as following:

```{r}
all.markers <- Detect_single_marker_all(mca.spleen, category = "Flow", clusters_to_detect = c("Marginal zone B cell(Spleen)", "NK cell(Spleen)"))
generate_report(mca.spleen, all.markers, fpath = ".", fname = "mca.spleen")
```

# Reproduction

All code and R object needed to generate the results of sc2marker manuscript can be found [here](https://doi.org/10.5281/zenodo.5703652).
Also, you can find the original data of MCA, HCA and Stromal data as following:
[MCA](https://figshare.com/articles/dataset/MCA_DGE_Data/5435866),
[HCA](https://atlas.fredhutch.org/nygc/multimodal-pbmc/),
[niche stromal](https://nicheview.shiny.embl.de/)



