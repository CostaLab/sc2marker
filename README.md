# sc2marker

sc2marker is a toolkit for identification of cell specific markers for clusters of cells in a single cell RNA sequencing data set. In particular, sc2marker focus on markers that have valid antibodies to be used for either imaging or flow analyis.

## Install

You can install sc2marker with the following commands: 

```{r}
install.packages("devtools")
devtools::install_github("https://github.com/CostaLab/sc2marker", build_vignettes = TRUE)
```

## Detect markers for given cell cluster

To run sc2marker you need to execute the following command, providing a clustered single  cell data sets (as Seurat object), the cell type of interest and the antibody databased (IHC, ICC or Flow). 

```{r}
require(sc2marker)
require(Seurat)
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

To use customized gene set, you can run the following command. sc2marker will recognise genes and ignore the cases.

```
nk.markers <- Detect_single_marker(mca.spleen, id = "NK cell(Spleen)", category = "Flow", geneset = c("CD19", "GeneA", "welcome2022"), org = "mouse")
get_antibody(nk.markers, rm.noab = F)
```

## To include your own antibody database

Example of self made antibody database can be found in [Self-made database](Self_DB.csv). Be sure to save the database within your R path and update the code above to reflect the file location in your system. 

```
nk.markers <- Detect_single_marker(mca.spleen, id = "NK cell(Spleen)", category = "Flow", org = "mouse", self.db = "/Path to Self_DB.csv/")
get_antibody(nk.markers, self.db = "/Path to Self_DB.csv/", org = "mouse")
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

All code and R object needed to generate the results of sc2marker manuscript can be found in [zenodo](https://zenodo.org/record/5854754#.YeknMf7MKUm).

The scripts below reproduce individual analysis. Please note that some packages might require manual installation ( 'seurat-data' and 'ggsci'). Also, configure your R paths before running these scripts.  

### Code for MCA Spleen (Mouse)
Example for MCA-Spleen in [Rmarkdown](Example_MCA_Spleen.Rmd) format and in [html](Example_MCA_Spleen.html) format.


### Code for Human Bone marrow (human-BM)
[SeuratData required] (https://github.com/satijalab/seurat-data) [Example_human_BM.Rmd](Example_human_BM.Rmd)

### Code for Human PBMC&lung (human-PBMC&lung)
[Example_human_PBMC&lung.Rmd](Example_human-PBMC&lung.Rmd)

### Code for Human Bone marrow (human-PBMC)
[Example_human_PBMC.Rmd](Example_human_PBMC.Rmd)

Also, you can find the original data of MCA, human-PBMC and human-PBMC&lung and Stromal data as following:
[MCA](https://figshare.com/articles/dataset/MCA_DGE_Data/5435866),
[human-PBMC&lung](https://archive.softwareheritage.org/browse/revision/1c7fcabb18a1971dc4d6e29bc3ed4f6f36b2361f/),
[human-PBMC](https://atlas.fredhutch.org/nygc/multimodal-pbmc/),
[human-BM](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128639),
[niche stromal](https://nicheview.shiny.embl.de/)



