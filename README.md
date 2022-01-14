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
nk.markers <- Detect_single_marker(mca.spleen, id = "NK cell(Spleen)", category = "Flow", org = "mouse")
```

You can display the results as a table with the command:

```{r}
get_antibody(nk.markers, org = "mouse")
```

and you can generate ridge plot with the following command:

```{r}
plot_ridge(mca.spleen, id = "NK cell(Spleen)", genes = nk.markers[1:9,]$gene, ncol = 3, assay = "RNA", aggr.other = F)
```

To use customized gene set, you can run the following command. (Relax, sc2marker will recognise genes and ignore the cases.)

```
nk.markers <- Detect_single_marker(mca.spleen, id = "NK cell(Spleen)", category = "Flow", geneset = c("CD19", "GeneA", "welcome2022"), org = "mouse")
get_antibody(nk.markers, rm.noab = F)
```


## calculate markers for all cell clusters

To calculate markers for all cell clusters, you can do by following command:

```{r}
all.markers <- Detect_single_marker_all(mca.spleen, category = "Flow", org = "mouse")
```

To Check T cell markers from results of all clusters, and get the antibody information, you can do following.

```{r}
t.markers <- all.markers[["T cell(Spleen)"]]
get_antibody(t.markers, org = "mouse")
```

## Generate report

To automatically generate the sc2marker report of all cell clusters, you can run following command:

```{r}
generate_report(mca.spleen, all.markers, fpath = ".", fname = "mca.spleen")
```

Or you can only analysis subset of cell clusters (B cell and NK cell) and generate the report as following:

```{r}
all.markers <- Detect_single_marker_all(mca.spleen, category = "Flow", clusters_to_detect = c("Marginal zone B cell(Spleen)", "NK cell(Spleen)"), org = "mouse")
generate_report(mca.spleen, all.markers, fpath = ".", fname = "mca.spleen")
```

# Reproduction

### vignette of MCA Spleen (Mouse)
[vignette_MCA_Spleen.Rmd](vignette_MCA_Spleen.Rmd)


### vignette of HCA Bone marrow (Human)
[vignette_human_BM.Rmd](vignette_HCA_BM.Rmd)

All code and R object needed to generate the results can be found [here](https://doi.org/10.5281/zenodo.5703652).
Also, you can find the original data of MCA, human-PBMC and human-PBMC&lung and Stromal data as following:
[MCA](https://figshare.com/articles/dataset/MCA_DGE_Data/5435866),
[human-PBMC&lung](https://archive.softwareheritage.org/browse/revision/1c7fcabb18a1971dc4d6e29bc3ed4f6f36b2361f/),
[human-PBMC](https://atlas.fredhutch.org/nygc/multimodal-pbmc/),
[human-BM](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128639),
[niche stromal](https://nicheview.shiny.embl.de/)



