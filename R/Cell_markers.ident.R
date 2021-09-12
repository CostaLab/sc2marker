
is_contigous_true_df <- function(is_sigs){
  ret_df <- data.frame(keep=FALSE, avgIdx=-1)
  if(any(is_sigs) & table(is_sigs)['TRUE'] == 1){
    ret_df$keep=TRUE
    ret_df$avgIdx = which(is_sigs == TRUE)
    return(ret_df)
  }
  return(ret_df)
}

#' Reset cell group idents
#' @param scrna seurat obj to be used
#' @param id interested cell group
#' @param panel Ident panel in seurat obj meta.data
#' @return seurat obj with re-set Idents
#'
makeid <- function(scrna, id, panel = "seurat_clusters"){
  levels(scrna@active.ident) <- c(levels(scrna@active.ident), "other")
  scrna@active.ident[scrna@active.ident != id] <- 'other'
  return(scrna)
}

get_exprs_frac <- function(x.df, step = 100){
  x.min <-  quantile(x.df$exp, c(.01, .99))[1]
  x.max <-  quantile(x.df$exp, c(.01, .99))[2]
  seq = seq(x.min, x.max, by = (x.max - x.min)/step)
  # seq = quantile(x.seq, seq(0.01, 0.99, step))
  return(seq)
}

#' Compute PRAUC for positive markers
#' @param scrna seurat obj to be used
#' @param gene gene to test
#' @param id interested cell group
#' @param step quantile steps
#' @return PRAUC of input gene
#'

get_gene_PRAUC_pos <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  # x.max <- max(data.id$exp)
  # x.min <- min(data.id$exp)

  # for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
  for (x.val in quantile(data.mat.surf$exp, seq(0, 1, step))) {

    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    # tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]

  PRAUC <- Area_Under_Curve(gene.prauc$x.rec,  gene.prauc$x.pre,
                            method = "trapezoid", na.rm = TRUE)
  return(PRAUC)
}

get_gene_PRAUC_pos_matrix <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  # x.max <- max(data.id$exp)
  # x.min <- min(data.id$exp)

  # for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
  for (x.val in quantile(data.mat.surf$exp, seq(0, 1, step))) {

    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    # tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]
  return(gene.prauc)
}

#' Compute PRAUC for negetive markers
#' @param scrna seurat obj to be used
#' @param gene gene to test
#' @param id interested cell group
#' @param step quantile steps
#' @return PRAUC of input gene
#'
get_gene_PRAUC_neg <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  # x.max <- max(data.id$exp)
  # x.min <- min(data.id$exp)

  # for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
  for (x.val in quantile(data.mat.surf$exp, seq(0, 0.85, step))) {

    tp <- sum(data.id$exp <= x.val)
    fp <- sum(data.other$exp <= x.val)
    # tn <- sum(data.other$exp > x.val)
    fn <- sum(data.id$exp > x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]

  PRAUC <- Area_Under_Curve(gene.prauc$x.rec,  gene.prauc$x.pre,
                            method = "trapezoid", na.rm = TRUE)
  return(PRAUC)
}

#' Compute Precision and Recall matrix for input gene
#' @param scrna seurat obj to be used
#' @param gene gene to test
#' @param id interested cell group
#' @param step quantile steps
#' @return PRAUC of input gene
#'

get_gene_PRAUC_matrix <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
    # for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {

    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]
  # gene.prauc <- gene.prauc[order(gene.prauc$x.rec, )]
  return(gene.prauc)
}

get_gene_PRAUC_matrix_neg <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  # x.max <- max(data.id$exp)
  # x.min <- min(data.id$exp)

  # for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
  for (x.val in quantile(data.mat.surf$exp, seq(0, 1, step))) {

    tp <- sum(data.id$exp < x.val)
    fp <- sum(data.other$exp < x.val)
    # tn <- sum(data.other$exp > x.val)
    fn <- sum(data.id$exp >= x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]

  return(gene.prauc)
}

get_gene_PRAUC_matrix.1 <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  data.mat.surf$id <- ifelse(data.mat.surf$id == id, id, "other")
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
  # for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {

    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }
  return(gene.prauc)
}

#' Compute PRAUC for each gene in input geneset and do the ranking
#' @param scrna seurat obj to be used
#' @param gene gene to test
#' @param cellgroup interested cell group
#' @param step quantile steps
#' @param geneset input geneset
#' @return PRAUC of input gene
#'

identify_single_marker <- function(scrna, cellgroup, geneset, step = 0.01){
  geneset <- intersect(rownames(scrna[["RNA"]]), geneset)
  markers <- FindMarkers(object = scrna, ident.1 = cellgroup, features = geneset)#, min.pct = 0.3)

  markers <- subset(markers, p_val_adj < 0.05)

  markers.pos <- subset(markers, avg_log2FC > 0)
  markers.neg <- subset(markers, avg_log2FC < 0)

  geneset <- rownames(markers.pos)

  gene.prauc <- data.frame(gene <- c(),
                           prauc <- c(),
                           direction <- c())

  for (genes in geneset) {
    de <- data.frame(genes, get_gene_PRAUC_pos(scrna, genes, cellgroup, step), "+")
    names(de)<-c("gene","prauc", "direction")
    gene.prauc <- rbind(gene.prauc, de)
  }

  geneset <- rownames(markers.neg)

  for (genes in geneset) {
    de <- data.frame(genes, get_gene_PRAUC_neg(scrna, genes, cellgroup, step), "-")
    names(de)<-c("gene","prauc", "direction")
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[order(gene.prauc$prauc, decreasing = T),]

  return(gene.prauc)
}



plot_grid_2 <- function(scrna, gene1, gene2, gene1.direc, gene2.direc, id, step = 0.01, return.obj = F){
  df <- data.frame(exp1 = scrna@assays$RNA@data[gene1,],
                   exp2 = scrna@assays$RNA@data[gene2,],
                   id = as.character(makeid(scrna, id)@active.ident))
  if (gene1.direc == "+") {
    x.split <- get_split_pos(scrna, gene1, id, step = step)$x.val
  }else{
    x.split <- get_split_neg(scrna, gene1, id, step = step)$x.val
  }

  left.cells <- GetCellNames(scrna, gene1, x.split, gene1.direc)
  scrna.left <- subset(scrna, cells = left.cells)

  if (gene2.direc == "+") {
    y.split <- get_split_pos(scrna.left, gene2, id, step = step)$x.val
  }else{
    y.split <- get_split_neg(scrna.left, gene2, id, step = step)$x.val
  }

  left.cells <- GetCellNames(scrna.left, gene2, y.split, gene2.direc)
  scrna.left <- subset(scrna.left, cells = left.cells)

  tp <- sum(scrna.left@active.ident == id)
  fp <- sum(scrna.left@active.ident != id)
  x.pre <- tp/(tp+fp)
  x.pre <- round(x.pre, 3)
  x.rec <- sum(scrna.left@active.ident == id)/sum(scrna@active.ident == id)
  x.rec <- round(x.rec, 3)

  g <- ggplot(df)+
    geom_point(aes(exp1, exp2, color = id))

  if (gene2.direc == "+") {
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = max(exp2)), linetype = 2)
  }else{
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = min(exp2)), linetype = 2)
  }

  if(gene1.direc == "+") {
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = max(exp1), yend = y.split), linetype = 2)
  }else{
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = min(exp1), yend = y.split), linetype = 2)
  }

  g <- g + xlab(paste(gene1, "--Recall: ", x.rec,"; Precision: ", x.pre, sep = "")) + ylab(paste(gene2))
  if (return.obj) {
    return(g)
  }
  print(g)
}

plot_marker_c_umap <- function(scrna, gene1, gene2, gene1.direc, gene2.direc, id, step = 0.01){
  df <- data.frame(exp1 = scrna@assays$RNA@data[gene1,],
                   exp2 = scrna@assays$RNA@data[gene2,],
                   id = as.character(makeid(scrna, id)@active.ident))
  if (gene1.direc == "+") {
    x.split <- get_split_pos(scrna, gene1, id, step = step)$x.val
  }else{
    x.split <- get_split_neg(scrna, gene1, id, step = step)$x.val
  }

  left.cells <- GetCellNames(scrna, gene1, x.split, gene1.direc)
  scrna.left <- subset(scrna, cells = left.cells)

  if (gene2.direc == "+") {
    y.split <- get_split_pos(scrna.left, gene2, id, step = step)$x.val
  }else{
    y.split <- get_split_neg(scrna.left, gene2, id, step = step)$x.val
  }

  cells.1 <- GetCellNames(scrna, gene1, x.split, gene1.direc)
  cells.2 <- GetCellNames(scrna, gene2, y.split, gene2.direc)
  cells.1.2 <- intersect(cells.1, cells.2)

  Idents(scrna) <- "other"
  Idents(scrna, cells = cells.1) <- paste(gene1, "positive")
  p1 <- UMAPPlot(scrna, cols = c('red', "grey"))

  Idents(scrna) <- "other"
  Idents(scrna, cells = cells.2) <- paste(gene2, "positive")
  p2 <- UMAPPlot(scrna, cols = c('blue', "grey"))

  Idents(scrna) <- "other"
  Idents(scrna, cells = cells.1.2) <- paste(gene1, "+", gene2)
  p3 <- UMAPPlot(scrna, cols = c('black',  "grey"))

  scrna <- makeid(scrna, id)
  p4<- UMAPPlot(scrna)

  plot_grid(p4, p1, p2, p3)
}

#' UmapPlot based on results from Markers Greedy Search
#' @param scrna seurat objects
#' @param df.split results from Greedy Search
#' @param id cell groups to identify
#' @return
#'
Combine_marker_umap <- function(scrna, df.split, id){
  idents.c <- Idents(scrna)
  gene1 <- df.split$gene[1]
  gene2 <- df.split$gene[2]
  gene1.split <- df.split$split.value[1]
  gene2.split <- df.split$split.value[2]

  cells.1 <- GetCellNames(cbmc, gene1, gene1.split, df.split$direction[1])
  cells.2 <- GetCellNames(cbmc, gene2, gene2.split, df.split$direction[2])
  cells.1.2 <- intersect(cells.1, cells.2)

  Idents(scrna) <- "other"
  Idents(scrna, cells = cells.1) <- paste(gene1, "filter")
  p1 <- UMAPPlot(scrna, cols = c('red', "grey"))

  Idents(scrna) <- "other"
  Idents(scrna, cells = cells.2) <- paste(gene2, "filter")
  p2 <- UMAPPlot(scrna, cols = c('blue', "grey"))

  Idents(scrna) <- "other"
  Idents(scrna, cells = cells.1.2) <- paste(gene1, "+", gene2)
  p3 <- UMAPPlot(scrna, cols = c('black',  "grey"))

  Idents(scrna) <- idents.c
  scrna <- makeid(scrna, id)
  p4 <- UMAPPlot(scrna)

  plot_grid(p1, p2, p3, p4)
}


GetCellNames <- function(scrna, gene, value, direction){
  df <- data.frame(exp = scrna@assays$RNA@data[gene,])
  if (direction == "+") {
    df <- subset(df, exp >= value)
  }else{
    df <- subset(df, exp <= value)
  }
  return(rownames(df))
}

get_split <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {

    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    # tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    x.margin <- tp - fn - fp
    # x.margin <- sum(data.id[,1] - x.val) + sum(x.val - data.other[,1])
    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,1]
  # x.split <- max(gene.prauc[,1])
  return(x.split)
}

plot_grid_1 <- function(scrna, gene, id, step = 0.01){
  df <- data.frame(exp1 = scrna@assays$RNA@data[gene,],
                   exp2 = rnorm(length(scrna@assays$RNA@data[gene,]), sd = 0.1),
                   id = as.character(makeid(scrna, id)@active.ident))

  x.split <- get_split(scrna, gene, id, step = step)

  ggplot(df)+
    geom_point(aes(exp2, exp1, color = id)) +
    # geom_vline(xintercept=x.split)+
    geom_hline(yintercept=x.split, linetype = 2)+
    # geom_segment(aes(x = x.split, y = y.split, xend = max(exp1), yend = y.split), linetype = 2)+
    # geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = max(exp2)), linetype = 2)+
    xlab(paste(gene)) + ylab("Expression") +
    facet_grid(~id) +
    theme(#axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}



#' Detect greedy step by step filters for cell group identification
#' @param scrna seurat object
#' @param cellgroup, interested cell group
#' @param depth how much steps to filter
#' @param geneset Gene set to be used
#' @return

plot_filter_stepbystep_first2 <- function(scrna, df.split, id, step = 0.01){
  gene1 <- df.split$gene[1]
  gene2 <- df.split$gene[2]
  df <- data.frame(exp1 = scrna@assays$RNA@data[gene1,],
                   exp2 = scrna@assays$RNA@data[gene2,],
                   id = as.character(makeid(scrna, id)@active.ident))

  x.split <- df.split$split.value[1]
  y.split <- df.split$split.value[2]

  ggplot(df)+
    geom_point(aes(exp1, exp2, color = id)) +
    geom_segment(aes(x = x.split, y = y.split, xend = max(exp1), yend = y.split), linetype = 2)+
    geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = max(exp2)), linetype = 2)+
    xlab(paste(gene1)) + ylab(paste(gene2))
}

plot_filter_combination_first2 <- function(scrna, df.split, id, step = 0.01){
  gene1 <- df.split$gene[1]
  gene2 <- df.split$gene[2]
  df <- data.frame(exp1 = scrna@assays$RNA@data[gene1,],
                   exp2 = scrna@assays$RNA@data[gene2,],
                   ID = as.character(makeid(scrna, id)@active.ident))

  x.split <- df.split$split.value[1]
  y.split <- df.split$split.value[2]

  g <- ggplot(df)+
    geom_point(aes(exp1, exp2, color = ID))

  if (df.split$direction[2] == "+") {
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = max(exp2)), linetype = 2)
  }else{
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = min(exp2)), linetype = 2)
  }

  if(df.split$direction[1] == "+") {
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = max(exp1), yend = y.split), linetype = 2)
  }else{
    g <- g + geom_segment(aes(x = x.split, y = y.split, xend = min(exp1), yend = y.split), linetype = 2)
  }

  g <- g + xlab(paste(gene1)) + ylab(paste(gene2))
  print(g)
}


grid_plot <- function(scrna, df.split, cellgroup){
  x.len <- length(df.split$gene)
  plots = list()
  # scrna.id <- scrna@active.ident
  scrna <- makeid(scrna, cellgroup)
  scrna.left <- scrna
  for (gene in df.split$gene) {
    plots[[gene]] <- VlnPlot(scrna.left, features = gene) + geom_hline(yintercept=df.split[df.split$gene == gene, ]$split.value, linetype = 1)
    left.cells <- GetCellNames(scrna.left, gene, df.split[df.split$gene == gene, ]$split.value, df.split[df.split$gene == gene, ]$direction)
    scrna.left <- subset(scrna.left, cells = left.cells)
    tp <- sum(scrna.left@active.ident == cellgroup)
    fp <- sum(scrna.left@active.ident != cellgroup)
    x.pre <- round(tp/(tp+fp), 3)
    x.rec <- sum(scrna.left@active.ident == cellgroup)/sum(scrna@active.ident == cellgroup)
    x.rec <- round(x.rec, 3)
    plots[[gene]] <- plots[[gene]] + xlab(paste("Precision:", x.pre, " Recall:", x.rec, sep = ""))
    plots[[gene]] <- plots[[gene]] + ggtitle(paste(gene, df.split[df.split$gene == gene, ]$direction))
  }
  plot_grid(plotlist=plots)
}



marker_stepbystep <- function(scrna, cellgroup, depth = 2, geneset = NULL, step = 0.01){
  if (is.null(geneset)) {
    geneset <- rownames(scrna)
  }
  geneset <- intersect(rownames(scrna[["RNA"]]), geneset)
  len.id <- sum(scrna@active.ident == cellgroup)
  len.other <- sum(scrna@active.ident != cellgroup)
  df.split <- data.frame()

  for (variable in seq(depth)) {
    markers <- FindMarkers(scrna, ident.1 = cellgroup, features = geneset)
    markers <- subset(markers, p_val_adj < 0.05)

    markers.pos <- subset(markers, avg_log2FC > 0)
    markers.neg <- subset(markers, avg_log2FC < 0)

    gene.prauc <- data.frame()

    for (genes in rownames(markers.pos)) {
      de <- data.frame(genes, get_split_pos(scrna, genes, cellgroup, step))
      names(de)<-c("gene", "split.value","filtered.adj", "direction")
      gene.prauc <- rbind(gene.prauc, de)
    }
    for (genes in rownames(markers.neg)) {
      de <- data.frame(genes, get_split_neg(scrna, genes, cellgroup, step))
      names(de)<-c("gene", "split.value","filtered.adj", "direction")
      gene.prauc <- rbind(gene.prauc, de)
    }
    gene.prauc <- gene.prauc[order(gene.prauc$filtered.adj, decreasing = T),]

    df.split <- rbind(df.split, gene.prauc[1,])

    left.cells <- GetCellNames(scrna, gene.prauc[1,1], gene.prauc[1,2], gene.prauc[1,4])
    scrna <- subset(scrna, cells = left.cells)
  }
  return(df.split)
}


plot_combine_PRAUC <- function(scrna, gene1, gene2, id, step = 0.01, return.obj = F){
  x.df <- get_PRAUC_matrix_combine_markers(scrna, gene1, gene2, id = id, step = step)
  g <- ggplot() + geom_line(x.df, mapping = aes(x.rec, x.pre)) + xlim(c(0,1)) + ylim(c(0,1))
  if (!return.obj) {
    print(g)
  }else{
    return(g)
  }
}


Detect_combine_markers <- function(scrna, id, geneset, step = 0.1, max.gene = NULL){
  x.split <- identify_single_marker(scrna, id, markers.allset, step = step)
  x.split <- subset(x.split, prauc > 0.2)
  if (max.gene & length(x.split$gene) > max.gene) {
    x.split <- x.split[1:max.gene, ]
  }

  gene.combn <- t(combn(x.split$gene, 2))
  gene.combn.prauc <- data.frame()

  for (gene.iter in 1:length(gene.combn[,1])) {
    genes <- gene.combn[gene.iter,]
    gene1 <- genes[1]
    gene2 <- genes[2]
    if (x.split[x.split$gene == gene1, ]$direction == "+") {
      if (x.split[x.split$gene == gene2, ]$direction == "+") {
        gene.prauc <- get_PRAUC_matrix_combine_markers(scrna = scrna, gene1 = gene1, gene2 = gene2, id = id, step = step)

      }else{
        gene.prauc <- get_PRAUC_matrix_combine_markers_neg_pos(scrna = scrna, gene1 = gene1, gene2 = gene2, id = id, step = step)
      }
    }else{
      if (x.split[x.split$gene == gene2, ]$direction == "+") {
        gene.prauc <- get_PRAUC_matrix_combine_markers_neg_pos(scrna = scrna, gene1 = gene2, gene2 = gene1, id = id, step = step)
      }else{
        gene.prauc <- get_PRAUC_matrix_combine_markers_neg(scrna = scrna, gene1 = gene1, gene2 = gene2, id = id, step = step)
      }
    }
    PRAUC <- Area_Under_Curve(gene.prauc$x.rec,  gene.prauc$x.pre,
                              method = "trapezoid", na.rm = TRUE)
    de <- data.frame(gene1, gene2, x.split[x.split$gene == gene1, ]$direction, x.split[x.split$gene == gene2, ]$direction, PRAUC)
    colnames(de) <- c("Gene1", "Gene2", "Direction1", "Direction2", "PRAUC")
    gene.combn.prauc <- rbind(gene.combn.prauc, de)
  }
  gene.combn.prauc <- gene.combn.prauc[order(gene.combn.prauc$PRAUC, decreasing = T),]
  return(gene.combn.prauc)
}


get_split_pos_tp <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  x.factor <- length(data.other$id)/length(data.id$id)

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)

  for (x.val in seq(x.min + (x.max - x.min)/100, x.max, (x.max - x.min)/100)) {
    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tp*x.factor - fp
    }
    else{
      x.margin <- tp - fp
    }
    de <- data.frame(x.val, x.margin, "+")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}

get_split_neg_tp <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  x.factor <- length(data.other$id)/length(data.id$id)

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {

    tp <- sum(data.id$exp < x.val)
    fp <- sum(data.other$exp < x.val)
    tn <- sum(data.other$exp >= x.val)
    fn <- sum(data.id$exp >= x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tp*x.factor - fp
    }
    else{
      x.margin <- tp - fp
    }
    de <- data.frame(x.val, x.margin, "-")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}


get_PRAUC_matrix_combine_markers_neg_pos <- function(scrna, gene1,  gene2, id, step = 0.01){
  data.mat.surf <- data.frame(gene1 = scrna@assays$RNA@data[gene1,],
                              gene2 = scrna@assays$RNA@data[gene2,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  gene1.min <- min(data.id$gene1)
  gene2.max <- max(data.id$gene2)
  gene1.range <- max(data.id$gene1) - gene1.min
  gene2.range <- gene2.max - min(data.id$gene2)

  for (x.seq in seq(0.05, 0.95, step)) {

    x.val <- gene1.min + x.seq*gene1.range
    y.val <- gene2.max - x.seq*gene2.range

    # x.val = quantile(data.id$gene1, x.seq)
    # y.val = quantile(data.id$gene2, x.seq)

    tp <- sum(data.id$gene1 >= x.val & data.id$gene2 <= y.val)
    fp <- sum(data.other$gene1 >= x.val & data.other$gene2 <= y.val)
    # tn <- sum(data.other$gene1 < x.val | data.other$gene2 < y.val)
    fn <- sum(data.id$gene1 < x.val | data.id$gene2 > y.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]
  return(gene.prauc)
}

get_PRAUC_matrix_combine_markers_neg <- function(scrna, gene1,  gene2, id, step = 0.01){
  data.mat.surf <- data.frame(gene1 = scrna@assays$RNA@data[gene1,],
                              gene2 = scrna@assays$RNA@data[gene2,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  for (x.seq in seq(0.05, 1, step)) {

    # x.val <- gene1.min + x.seq*gene1.range
    # y.val <- gene2.min + x.seq*gene2.range

    x.val = quantile(data.id$gene1, x.seq)
    y.val = quantile(data.id$gene2, x.seq)

    tp <- sum(data.id$gene1 < x.val & data.id$gene2 < y.val)
    fp <- sum(data.other$gene1 < x.val & data.other$gene2 < y.val)
    # tn <- sum(data.other$gene1 < x.val | data.other$gene2 < y.val)
    fn <- sum(data.id$gene1 >= x.val | data.id$gene2 >= y.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]
  return(gene.prauc)
}

get_PRAUC_matrix_combine_markers <- function(scrna, gene1,  gene2, id, step = 0.01){
  data.mat.surf <- data.frame(gene1 = scrna@assays$RNA@data[gene1,],
                              gene2 = scrna@assays$RNA@data[gene2,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  for (x.seq in seq(0, 0.95, step)) {

    # x.val <- gene1.min + x.seq*gene1.range
    # y.val <- gene2.min + x.seq*gene2.range

    x.val = quantile(data.id$gene1, x.seq)
    y.val = quantile(data.id$gene2, x.seq)

    tp <- sum(data.id$gene1 >= x.val & data.id$gene2 >= y.val)
    fp <- sum(data.other$gene1 >= x.val & data.other$gene2 >= y.val)
    # tn <- sum(data.other$gene1 < x.val | data.other$gene2 < y.val)
    fn <- sum(data.id$gene1 < x.val | data.id$gene2 < y.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]
  return(gene.prauc)
}

get_split_pos.1 <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)
  x.num <- 1/step

  x.factor <- (mean(data.id$exp) + 0.001)/(mean(data.other$exp) + 0.001)

  x.size.factor <- 10*(length(data.id$exp)/length(data.other$exp))

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    tp <- sum(data.id$exp > x.val)
    fp <- sum(data.other$exp > x.val)
    tn <- sum(data.other$exp <= x.val)
    fn <- sum(data.id$exp <= x.val)

    x.margin <- tp - fp#*x.size.factor
    # x.margin <- tp - fp*x.size.factor
    # x.margin <- (x.factor**2)*x.margin

    # x.margin.adj <- x.margin*x.factor
    x.margin.adj <- x.margin*(tp/length(data.id$exp))*(tn/length(data.other$exp))*(x.factor**2)

    de <- data.frame(x.val, x.margin, x.margin.adj, tp, fp, "+")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}


get_split_pos.1 <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)
  x.num <- 1/step

  x.factor <- (mean(data.id$exp) + 0.001)/(mean(data.other$exp) + 0.001)

  x.size.factor <- 10*(length(data.id$exp)/length(data.other$exp))

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    tp <- sum(data.id$exp > x.val)
    fp <- sum(data.other$exp > x.val)
    tn <- sum(data.other$exp <= x.val)
    fn <- sum(data.id$exp <= x.val)

    x.margin <- tp - fp#*x.size.factor
    # x.margin <- tp - fp*x.size.factor
    # x.margin <- (x.factor**2)*x.margin

    # x.margin.adj <- x.margin*x.factor
    x.margin.adj <- x.margin*(tp/length(data.id$exp))*(tn/length(data.other$exp))*(x.factor**2)

    de <- data.frame(x.val, x.margin, x.margin.adj, tp, fp, "+")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}


get_split_neg.1 <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)
  x.num <- 1/step

  x.factor <- (mean(data.other$exp) + 0.001)/(mean(data.id$exp) + 0.001)

  x.size.factor <- 10*(length(data.id$exp)/length(data.other$exp))

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    tp <- sum(data.id$exp < x.val)
    fp <- sum(data.other$exp < x.val)
    tn <- sum(data.other$exp >= x.val)
    fn <- sum(data.id$exp >= x.val)

    x.margin <- tp - fp#*x.size.factor
    # x.margin <- tp - fp*x.size.factor
    # x.margin <- (x.factor**2)*x.margin

    # x.margin.adj <- x.margin*x.factor
    x.margin.adj <- x.margin*(tp/length(data.id$exp))*(tn/length(data.other$exp))*(x.factor**2)

    de <- data.frame(x.val, x.margin, x.margin.adj, tp, fp, "-")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}


identify_single_marker.1 <- function(scrna, id, geneset = NULL, step = 0.01){
  if (!is.null(geneset)) {
    if (length(intersect(rownames(scrna[["RNA"]]), geneset)) == 0) {
      print("Genes not find.")
    }else{
      markers.de <- FindMarkers(scrna, ident.1 = id, features = intersect(rownames(scrna[["RNA"]]), geneset), test.use = "t")
    }
  }else{
    markers.de <- FindMarkers(scrna, ident.1 = id, test.use = "t")
  }
  markers.de <- markers.de[markers.de$p_val < 0.05, ]
  markers.pos <- markers.de[markers.de$avg_log2FC > 0, ]
  markers.neg <- markers.de[markers.de$avg_log2FC < 0, ]
  df <- data.frame()
  for (gene in rownames(markers.pos)) {
    df.s <- get_split_pos.1(scrna, gene, id, step = step)
    df <- rbind(df, df.s)
  }
  for (gene in rownames(markers.neg)) {
    df.s <- get_split_neg.1(scrna, gene, id, step = step)
    df <- rbind(df, df.s)
  }
  df <- df[order(df$x.margin.adj, decreasing = T), ]
  return(df)
}

get_split_pos_neg <- function(scrna, gene1, gene2, id, step = 0.01) {
  data.mat.surf <- data.frame(exp1 = normalize(scrna@assays$RNA@data[gene1,]),
                              exp2 = normalize(scrna@assays$RNA@data[gene2,]),
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.max <- 0
  x.min <- 1
  x.num <- 1/step

  x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))

  for (x.val in seq(0, 1, step)) {
    for (y.val in seq(0, 1, step)) {
      data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 < y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 < y.val, ]
      fp <- nrow(data.o.a)
      data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 >= y.val, ]
      fn <- nrow(data.i.b)
      data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 >= y.val, ]
      tn <- nrow(data.o.b)

      x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((y.val - data.i.a$exp2))*nrow(data.o.b) +
        sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((data.o.b$exp2 - y.val))*nrow(data.i.a) -
        sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((y.val - data.o.a$exp2))*nrow(data.o.b)
      x.margin <- tp - fp
      x.margin <- x.margin/data.all.l
      x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a
      de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
      gene.prauc <- rbind(gene.prauc, de)
    }
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,3]
  y.val <- x.split[,4]
  x.margin <- x.split[,5]
  data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 < y.val, ]
  tp <- nrow(data.i.a)
  data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 < y.val, ]
  fp <- nrow(data.o.a)
  data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 >= y.val, ]
  fn <- nrow(data.i.b)
  data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 >= y.val, ]
  tn <- nrow(data.o.b)

  x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((y.val - data.i.a$exp2))*nrow(data.o.b) +
    sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((data.o.b$exp2 - y.val))*nrow(data.i.a) -
    sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((y.val - data.o.a$exp2))*nrow(data.o.b)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a

  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "-")
  return(x.split)
}


get_split_neg_neg <- function(scrna, gene1, gene2, id, step = 0.01) {
  data.mat.surf <- data.frame(exp1 = normalize(scrna@assays$RNA@data[gene1,]),
                              exp2 = normalize(scrna@assays$RNA@data[gene2,]),
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.max <- 0
  x.min <- 1
  x.num <- 1/step
  # x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))

  for (x.val in seq(0, 1, step)) {
    for (y.val in seq(0, 1, step)) {
      data.i.a <- data.id[data.id$exp1 < x.val & data.id$exp2 < y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 < x.val & data.other$exp2 < y.val, ]
      fp <- nrow(data.o.a)
      data.i.b <- data.id[data.id$exp1 >= x.val | data.id$exp2 >= y.val, ]
      fn <- nrow(data.i.b)
      data.o.b <- data.other[data.other$exp1 >= x.val | data.other$exp2 >= y.val, ]
      tn <- nrow(data.o.b)

      x.margin <- tp - fp
      x.margin <- x.margin/data.all.l

      de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
      gene.prauc <- rbind(gene.prauc, de)
    }
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,3]
  y.val <- x.split[,4]
  x.margin <- x.split[,5]
  data.i.a <- data.id[data.id$exp1 < x.val & data.id$exp2 < y.val, ]
  tp <- nrow(data.i.a)
  data.o.a <- data.other[data.other$exp1 < x.val & data.other$exp2 < y.val, ]
  fp <- nrow(data.o.a)
  data.i.b <- data.id[data.id$exp1 >= x.val | data.id$exp2 >= y.val, ]
  fn <- nrow(data.i.b)
  data.o.b <- data.other[data.other$exp1 >= x.val | data.other$exp2 >= y.val, ]
  tn <- nrow(data.o.b)

  x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
    sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
    sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)

  x.margin.a <- -x.margin.a
  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a

  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "-", "-")
  return(x.split)
}

get_split_pos_pos <- function(scrna, gene1, gene2, id, step = 0.01) {
  data.mat.surf <- data.frame(exp1 = normalize(scrna@assays$RNA@data[gene1,]),
                              exp2 = normalize(scrna@assays$RNA@data[gene2,]),
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.max <- 0
  x.min <- 1
  x.num <- 1/step

  x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))

  for (x.val in seq(0, 1, step)) {
    for (y.val in seq(0, 1, step)) {
      data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
      fp <- nrow(data.o.a)
      data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
      fn <- nrow(data.i.b)
      data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
      tn <- nrow(data.o.b)

      # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
      #   sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
      #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
      x.margin <- tp - fp
      x.margin <- x.margin/data.all.l

      # x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)

      # de <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+")
      de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
      gene.prauc <- rbind(gene.prauc, de)
    }
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,3]
  y.val <- x.split[,4]
  x.margin <- x.split[,5]
  data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
  tp <- nrow(data.i.a)
  data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
  fp <- nrow(data.o.a)
  data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
  fn <- nrow(data.i.b)
  data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
  tn <- nrow(data.o.b)

  x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
    sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
    sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)
  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+")
  return(x.split)
}


get_gene_score_pos <- function(scrna, gene, id, step = 0.01, do.magic = F) {
  if (do.magic) {
    data.mat.surf <- data.frame(exp = normalize(scrna@assays$MAGIC_RNA@data[gene,]),
                                id = as.character(scrna@active.ident))
  }else{
    data.mat.surf <- data.frame(exp = normalize(scrna@assays$RNA@data[gene,]),
                                id = as.character(scrna@active.ident))
  }
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.max <- 0
  x.min <- 1
  x.num <- 1/step

  x.factor <- (mean(data.id$exp) + 0.001)/(mean(data.other$exp) + 0.001)
  # x.size.factor <- 10*(length(data.id$exp)/length(data.other$exp))

  for (x.val in seq(0, 1, step)) {
    # tp <- sum(data.id$exp > x.val)
    # fp <- sum(data.other$exp > x.val)
    # tn <- sum(data.other$exp <= x.val)
    # fn <- sum(data.id$exp <= x.val)
    data.i.a <- data.id[data.id$exp > x.val, ]
    data.o.a <- data.other[data.other$exp > x.val, ]
    data.i.b <- data.id[data.id$exp <= x.val, ]
    data.o.b <- data.other[data.other$exp <= x.val, ]

    x.margin <- sum((data.i.a$exp - x.val))*nrow(data.o.b) +
      sum((x.val -data.o.b$exp))*nrow(data.i.a) -
      sum((data.o.a$exp - x.val))*nrow(data.o.b)
    x.margin <- x.margin/data.all.l

    # x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor**2)

    # de <- data.frame(x.val, x.margin, x.margin.adj, tp, fp, "+")
    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]
  # data.i.a <- data.id[data.id$exp > x.val, ]
  # tp <- nrow(data.i.a)
  tp <- sum(data.id$exp > x.val)
  # data.o.a <- data.other[data.other$exp > x.val, ]
  # fp <- nrow(data.o.a)
  fp <- sum(data.other$exp > x.val)

  # data.i.b <- data.id[data.id$exp <= x.val, ]
  # fn <- nrow(data.i.b)
  # fn <- sum(data.id$exp <= x.val)

  # data.o.b <- data.other[data.other$exp <= x.val, ]
  # tn <- nrow(data.o.b)
  tn <- sum(data.other$exp <= x.val)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor**2)
  x.split <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, "+")

  # rownames(x.split) <- paste(gene)
  return(x.split)
}

get_gene_score_neg <- function(scrna, gene, id, step = 0.01, do.magic = F) {
  if (do.magic) {
    data.mat.surf <- data.frame(exp = normalize(scrna@assays$MAGIC_RNA@data[gene,]),
                                id = as.character(scrna@active.ident))
  }else{
    data.mat.surf <- data.frame(exp = normalize(scrna@assays$RNA@data[gene,]),
                                id = as.character(scrna@active.ident))
  }
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.max <- 0
  x.min <- 1
  x.num <- 1/step

  x.factor <- (mean(data.id$exp) + 0.001)/(mean(data.other$exp) + 0.001)

  for (x.val in seq(0, 1, step)) {
    data.i.a <- data.id[data.id$exp < x.val, ]
    data.o.a <- data.other[data.other$exp < x.val, ]
    data.i.b <- data.id[data.id$exp >= x.val, ]
    data.o.b <- data.other[data.other$exp >= x.val, ]

    x.margin <- sum((data.i.a$exp - x.val))*nrow(data.o.b) +
      sum((x.val -data.o.b$exp))*nrow(data.i.a) -
      sum((data.o.a$exp - x.val))*nrow(data.o.b)
    x.margin <- x.margin/data.all.l

    x.margin <- 0 - x.margin

    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]
  tp <- sum(data.id$exp > x.val)
  fp <- sum(data.other$exp > x.val)
  tn <- sum(data.other$exp <= x.val)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor**2)
  x.split <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, "-")

  return(x.split)
}



marker_stepbystep.2 <- function(scrna, cellgroup, depth = 2, geneset = NULL, step = 0.01, do.magic = F){
  if (is.null(geneset)) {
    geneset <- rownames(scrna[["RNA"]])
  }
  geneset <- intersect(rownames(scrna[["RNA"]]), geneset)
  len.id <- sum(scrna@active.ident == cellgroup)
  len.other <- sum(scrna@active.ident != cellgroup)
  df.split <- data.frame()

  suppressWarnings(
    for (variable in seq(depth)) {
      df.fc <- FoldChange(scrna, ident.1 = cellgroup, features = geneset)
      # markers <- get_gene_score_pos(scrna, ident.1 = cellgroup, features = geneset)
      # markers <- subset(markers, p_val < 0.05)

      markers.pos <- subset(df.fc, avg_log2FC > 0)
      markers.neg <- subset(df.fc, avg_log2FC < 0)

      gene.prauc <- data.frame()

      for (genes in rownames(markers.pos)) {
        de <- data.frame(get_gene_score_pos(scrna, genes, cellgroup, step, do.magic = do.magic))
        names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
        gene.prauc <- rbind(gene.prauc, de)
      }
      for (genes in rownames(markers.neg)) {
        de <- data.frame(get_gene_score_neg(scrna, genes, cellgroup, step, do.magic = do.magic))
        names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
        gene.prauc <- rbind(gene.prauc, de)
      }
      gene.prauc <- gene.prauc[order(gene.prauc$x.margin.adj, decreasing = T),]

      df <- gene.prauc[1,]
      if (df$direction == "+") {
        split.step <- get_split_pos(scrna, df$gene, cellgroup, do.magic = do.magic)
      }else{
        split.step <- get_split_neg(scrna, df$gene, cellgroup, do.magic = do.magic)
      }

      names(split.step) <- c("Gene", "split.value", "direction")

      df.split <- rbind(df.split, split.step)
      left.cells <- GetCellNames(scrna = scrna, gene = df$gene, value = df.split$split.value, direction = df.split$direction, do.magic = do.magic)
      # left.cells <- GetCellNames(scrna = scrna, gene = df$gene, value = 0, direction = df.split$direction)
      scrna <- subset(scrna, cells = left.cells)
      geneset <- setdiff(geneset, df$gene)
    }
  )
  return(df.split)
}


GetCellNames <- function(scrna, gene, value, direction, do.magic = F){
  if (do.magic) {
    df <- data.frame(exp = scrna@assays$MAGIC_RNA@data[gene,])
  }else{
    df <- data.frame(exp = scrna@assays$RNA@data[gene,])
  }
  if (direction == "+") {
    df <- subset(df, exp > value)
  }else{
    df <- subset(df, exp < value)
  }
  return(rownames(df))
}

get_split_pos <- function(scrna, gene, id, step = 0.01, do.magic = F) {
  if (do.magic) {
    data.mat.surf <- data.frame(exp = scrna@assays$MAGIC_RNA@data[gene,],
                                id = as.character(scrna@active.ident))
  }else{
    data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                                id = as.character(scrna@active.ident))
  }
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)
  x.min <- 0
  x.num <- 1/step
  x.factor <- as.numeric(length(data.mat.surf$id)/length(data.id$id))
  x.factor <- ifelse(x.factor > 5, x.factor, 5)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    # tp <- sum(data.id$exp >= x.val)
    # fp <- sum(data.other$exp >= x.val)
    tn <- sum(data.other$exp <= x.val)
    fn <- sum(data.id$exp <= x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tn - fn*x.factor
    }
    else{
      x.margin <- tn - fn*x.factor
    }
    de <- data.frame(gene, x.val, x.margin, "+")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  colnames(gene.prauc) <- c("Gene", "split.value", "split.margin", "direction")
  x.split <- gene.prauc[1,c(1,2,4)]
  rownames(x.split) <- paste(gene)
  return(x.split)
}

get_split_neg <- function(scrna, gene, id, step = 0.01, do.magic = F) {
  if (do.magic) {
    data.mat.surf <- data.frame(exp = scrna@assays$MAGIC_RNA@data[gene,],
                                id = as.character(scrna@active.ident))
  }else{
    data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                                id = as.character(scrna@active.ident))
  }
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.mat.surf[data.mat.surf$exp != 0, ]$exp)
  x.min <- 0
  x.num <- 1/step

  x.factor <- as.numeric(length(data.mat.surf$id)/length(data.id$id))
  x.factor <- ifelse(x.factor > 5, x.factor, 5)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    # tp <- sum(data.id$exp <= x.val)
    # fp <- sum(data.other$exp <= x.val)
    tn <- sum(data.other$exp > x.val)
    fn <- sum(data.id$exp > x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tn - fn*x.factor
    }
    else{
      x.margin <- tn - fn*x.factor
    }
    de <- data.frame(gene, x.val, x.margin, "-")
    gene.prauc <- rbind(gene.prauc, de)
  }
  # colnames(gene.prauc) <- c("Gene", "split.value", "direction")
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  colnames(gene.prauc) <- c("Gene", "split.value", "split.margin", "direction")
  x.split <- gene.prauc[1,c(1,2,4)]
  rownames(x.split) <- paste(gene)
  return(x.split)
}

get_split_pos_pos <- function(scrna, gene1, gene2, id, step = 0.01) {
  data.mat.surf <- data.frame(exp1 = normalize(scrna@assays$RNA@data[gene1,]),
                              exp2 = normalize(scrna@assays$RNA@data[gene2,]),
                              id = as.character(scrna@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.max <- 0
  x.min <- 1
  x.num <- 1/step

  x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))

  for (x.val in seq(0, 1, step)) {
    for (y.val in seq(0, 1, step)) {
      data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
      fp <- nrow(data.o.a)
      data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
      fn <- nrow(data.i.b)
      data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
      tn <- nrow(data.o.b)

      # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
      #   sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
      #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
      x.margin <- tp - fp
      x.margin <- x.margin/data.all.l

      # x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)

      # de <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+")
      de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
      gene.prauc <- rbind(gene.prauc, de)
    }
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,3]
  y.val <- x.split[,4]
  x.margin <- x.split[,5]
  data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
  tp <- nrow(data.i.a)
  data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
  fp <- nrow(data.o.a)
  data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
  fn <- nrow(data.i.b)
  data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
  tn <- nrow(data.o.b)
  x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
    sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
    sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)
  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+")
  return(x.split)
}


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

makeid <- function(scrna, id, panel = "seurat_clusters"){
  levels(scrna@active.ident) <- c(levels(scrna@active.ident), "other")
  scrna@active.ident[scrna@active.ident != id] <- 'other'
  return(scrna)
}

normalize <- function(x, ...) {
  return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}


grid_vln_plot <- function(scrna, df.split, cellgroup, ncol = 2, do.magic = F){
  x.len <- length(df.split$gene)
  plots = list()
  # scrna.id <- scrna@active.ident
  scrna <- makeid(scrna, cellgroup)
  if (do.magic) {
    scrna@active.assay = "MAGIC_RNA"
  }
  scrna.left <- scrna
  n = 1
  for (gene in df.split$Gene) {
    plots[[gene]] <- VlnPlot(scrna.left, features = gene) + geom_hline(yintercept=df.split[df.split$Gene == gene, ]$split.value, linetype = 1)
    # plot(VlnPlot(scrna.left, features = gene) + geom_hline(yintercept=df.split[df.split$Gene == gene, ]$split.value, linetype = 1))
    left.cells <- GetCellNames(scrna.left, gene, df.split[df.split$Gene == gene, ]$split.value, df.split[df.split$Gene == gene, ]$direction)
    scrna.left <- subset(scrna.left, cells = left.cells)
    tp <- sum(scrna.left@active.ident == cellgroup)
    fp <- sum(scrna.left@active.ident != cellgroup)
    x.pre <- round(tp/(tp+fp), 3)
    x.rec <- sum(scrna.left@active.ident == cellgroup)/sum(scrna@active.ident == cellgroup)
    x.rec <- round(x.rec, 3)
    plots[[gene]] <- plots[[gene]] + xlab(paste("Precision:", x.pre, " Recall:", x.rec, sep = ""))
    plots[[gene]] <- plots[[gene]] + ggtitle(paste("Filter", n, " ", gene, df.split[df.split$Gene == gene, ]$direction, sep = ""))
    n <- n + 1
    # plot(plots[[gene]])
    # plot(VlnPlot(scrna.left, features = gene) + geom_hline(yintercept=df.split[df.split$Gene == gene, ]$split.value, linetype = 1))
  }
  plot_grid(plotlist=plots, ncol = ncol)
}

grid_pie_plot <- function(scrna, df.split, cellgroup, ncol = 2){
  x.len <- length(df.split$gene)
  plots = list()
  # scrna.id <- scrna@active.ident
  # scrna <- makeid(scrna, cellgroup)
  scrna.left <- scrna
  n = 1
  for (gene in df.split$Gene) {
    df <- data.frame(idents = c(as.character(cellgroup), "other"),
                     value = c(sum(scrna.left@active.ident == cellgroup), sum(scrna.left@active.ident != cellgroup)))
    df <- df %>%
      arrange(desc(idents)) %>%
      mutate(prop = value / sum(df$value) *100) %>%
      mutate(ypos = cumsum(prop)- 0.5*prop )

    g <- ggplot(df, aes(x="", y=prop, fill=idents)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() +
      theme(legend.position="none") +
      geom_text(aes(y = ypos, label = idents), color = "white", size=4)

    plots[[gene]] <- g
    # print(g)

    left.cells <- GetCellNames(scrna.left, gene, df.split[df.split$Gene == gene, ]$split.value, df.split[df.split$Gene == gene, ]$direction)
    scrna.left <- subset(scrna.left, cells = left.cells)
    tp <- sum(scrna.left@active.ident == cellgroup)
    fp <- sum(scrna.left@active.ident != cellgroup)
    x.pre <- round(tp/(tp+fp), 2)
    x.rec <- sum(scrna.left@active.ident == cellgroup)/sum(scrna@active.ident == cellgroup)
    x.rec <- round(x.rec, 2)
    # plots[[gene]] <- plots[[gene]] + xlab(paste("Precision:", x.pre, " Recall:", x.rec, sep = ""))
    plots[[gene]] <- plots[[gene]] + ggtitle(paste("Filter", n, " ", gene, df.split[df.split$Gene == gene, ]$direction, "  ",
                                                   "PR:", x.pre, " RC:", x.rec, sep = ""))
    n <- n + 1
  }
  plot(plot_grid(plotlist=plots[1:length(plots)], ncol = ncol))
}


detect_marker <- function(scrna, id, step = 0.1, do.magic = F, category = NULL,
                          geneset = NULL, assay = "RNA", do.fast = F, min.pct = 0.1, use.all = F){
  if (is.null(geneset)) {
    geneset <- rownames(scrna[[assay]])
  }else{
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
  }
  if (use.all) {
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
    de <- FoldChange(scrna, ident.1 = id, features = geneset)
  }
  if (!use.all) {
    if (do.fast) {
      de <- FindMarkers(scrna, ident.1 = id, features = geneset)
      de <- subset(de, p_val_adj < 0.05)
    }else{
      de <- FoldChange(scrna, ident.1 = id, features = geneset)
      de <- de[de$pct.1 > min.pct, ]
    }
  }

  markers.pos <- subset(de, avg_log2FC > 0)
  markers.neg <- subset(de, avg_log2FC < 0)

  gene.rank.list <- data.frame()

  for (genes in rownames(markers.pos)) {
    de <- data.frame(get_gene_score_pos(scrna, genes, id, step, do.magic = do.magic))
    names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
    gene.rank.list <- rbind(gene.rank.list, de)
  }
  for (genes in rownames(markers.neg)) {
    de <- data.frame(get_gene_score_neg(scrna, genes, id, step, do.magic = do.magic))
    names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
    gene.rank.list <- rbind(gene.rank.list, de)
  }
  gene.rank.list <- gene.rank.list[order(gene.rank.list$x.margin.adj, decreasing = T),]
  return(gene.rank.list)
}
