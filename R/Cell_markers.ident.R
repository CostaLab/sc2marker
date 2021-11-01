
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


get_split_pos_neg <- function(scrna, gene1, gene2, id, step = 0.01, assay = "RNA", slot = "data") {
  scrna@active.assay <- assay
  # data.mat.surf <- FetchData(scrna, vars = c(gene1, gene2, "ident"), slot = slot)
  # colnames(data.mat.surf) <- c("exp1", "exp2", "id")
  # data.mat.surf$exp1 <- normalize(data.mat.surf$exp1)
  # data.mat.surf$exp2 <- normalize(data.mat.surf$exp2)
  #
  # gene.prauc <- data.frame(x.val <- c(),
  #                          margin <- c())
  #
  # data.id <- data.mat.surf[data.mat.surf$id == id,]
  # data.other <- data.mat.surf[data.mat.surf$id != id,]
  # data.id.l <- nrow(data.id)
  # data.o.l <- nrow(data.other)
  # data.all.l <- nrow(data.mat.surf)
  #
  # x.max <- 0
  # x.min <- 1
  # x.num <- 1/step
  #
  # x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))
  #
  # for (x.val in seq(0, 1, step)) {
  #   for (y.val in seq(0, 1, step)) {
  #     data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 < y.val, ]
  #     tp <- nrow(data.i.a)
  #     data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 < y.val, ]
  #     fp <- nrow(data.o.a)
  #     data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 >= y.val, ]
  #     fn <- nrow(data.i.b)
  #     data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 >= y.val, ]
  #     tn <- nrow(data.o.b)
  #
  #     x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((y.val - data.i.a$exp2))*nrow(data.o.b) +
  #       sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((data.o.b$exp2 - y.val))*nrow(data.i.a) -
  #       sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((y.val - data.o.a$exp2))*nrow(data.o.b)
  #     x.margin <- tp - fp
  #     x.margin <- x.margin/data.all.l
  #     x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a
  #     de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
  #     gene.prauc <- rbind(gene.prauc, de)
  #   }
  # }
  # gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  # x.split <- gene.prauc[1,]
  # x.val <- x.split[,3]
  # y.val <- x.split[,4]
  # x.margin <- x.split[,5]
  # data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 < y.val, ]
  # tp <- nrow(data.i.a)
  # data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 < y.val, ]
  # fp <- nrow(data.o.a)
  # data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 >= y.val, ]
  # fn <- nrow(data.i.b)
  # data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 >= y.val, ]
  # tn <- nrow(data.o.b)
  #
  # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((y.val - data.i.a$exp2))*nrow(data.o.b) +
  #   sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((data.o.b$exp2 - y.val))*nrow(data.i.a) -
  #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((y.val - data.o.a$exp2))*nrow(data.o.b)
  #
  # x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a

  data.mat.surf <- Seurat::FetchData(scrna, vars = c(gene1, gene2, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp1", "exp2", "id")
  data.mat.surf$exp1 <- normalize(data.mat.surf$exp1)
  data.mat.surf$exp2 <- normalize(data.mat.surf$exp2)

  sp1 = boxplot.stats(data.mat.surf$exp1)$out
  sp2 = boxplot.stats(data.mat.surf$exp2)$out

  `%notin%` <- Negate(`%in%`)
  data.mat.surf <- data.mat.surf[data.mat.surf$exp2 %notin% sp2 & data.mat.surf$exp1 %notin% sp1, ]

  if(max(data.mat.surf$exp2) == 0 | max(data.mat.surf$exp1) == 0){
    x.split <- data.frame(gene1, gene2, 0, 0, 0, 0, 0, 0, "+", "-")
    return(x.split)
  }

  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  size.factor <- nrow(data.mat.surf)/nrow(data.other)
  # size.factor <- 3*nrow(data.id)/nrow(data.other)
  # size.factor <- nrow(data.id)/nrow(data.other)
  # size.factor <- ifelse(size.factor < 3, 3, size.factor)

  for (x.val in seq(0.01, 1, 0.1)) {
    for (y.val in seq(0.01, 1, 0.1)) {
      data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 < y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 < y.val, ]
      fp <- nrow(data.o.a)
      # x.margin <- tp*size.factor - fp
      x.margin <- tp - fp*size.factor
      # x.margin <- tp*nrow(data.other) - fp*nrow(data.id)
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
  data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 < y.val, ]
  tp <- nrow(data.i.a)
  data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 < y.val, ]
  fp <- nrow(data.o.a)
  data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 >= y.val, ]
  tn <- nrow(data.o.b)
  # dist.o.b = 0
  #
  # if (nrow(data.o.b) >= 1) {
  #   for (i in 1:nrow(data.o.b)) {
  #     if (data.o.b[i,]$exp1 < x.val) {
  #       if (data.o.b[i,]$exp2 > y.val) {
  #         dist.o.b = dist.o.b + sqrt((data.o.b[i,]$exp1 - x.val)**2 + (data.o.b[i,]$exp2 - y.val)**2)
  #       }else{
  #         dist.o.b = dist.o.b + (x.val - data.o.b[i,]$exp1)
  #       }
  #     }else{
  #       dist.o.b = dist.o.b + (data.o.b[i,]$exp2 - y.val)
  #     }
  #   }
  # }
  # # for (i in 1:nrow(data.o.b)) {
  # #   if (data.o.b[i,]$exp1 < x.val) {
  # #     if (data.o.b[i,]$exp2 > y.val) {
  # #       dist.o.b = dist.o.b + sqrt((data.o.b[i,]$exp1 - x.val)**2 + (data.o.b[i,]$exp2 - y.val)**2)
  # #     }else{
  # #       dist.o.b = dist.o.b + (x.val - data.o.b[i,]$exp1)
  # #     }
  # #   }else{
  # #     dist.o.b = dist.o.b + (data.o.b[i,]$exp2 - y.val)
  # #   }
  # # }
  # dist.i.a <- sum(apply(data.i.a[,c(1,2)], 1, dist_p,  a2 = c(x.val, y.val)))
  # dist.o.a <- sum(apply(data.o.a[,c(1,2)], 1, dist_p,  a2 = c(x.val, y.val)))
  # x.margin.a <- dist.i.a - dist.o.a + dist.o.b

  # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((y.val - data.i.a$exp2))*nrow(data.o.b) -#+
  #   # sum((x.val - data.o.b$exp1))*nrow(data.i.a) + sum((data.o.b$exp2 - y.val))*nrow(data.i.a) -
  #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((y.val - data.o.a$exp2))*nrow(data.o.b)

  x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.a) + sum((y.val - data.i.a$exp2))*nrow(data.o.a) -#+
    # sum((x.val - data.o.b$exp1))*nrow(data.id) + sum((data.o.b$exp2 - y.val))*nrow(data.id) -
    sum((data.o.a$exp1 - x.val))*nrow(data.i.a) - sum((y.val - data.o.a$exp2))*nrow(data.i.a)

  # x.margin.a <- x.margin.a/(nrow(data.i.a) + nrow(data.o.a))

  x.margin.a2 <- sum((x.val - data.o.b$exp1))*nrow(data.id) + sum((data.o.b$exp2 - y.val))*nrow(data.id)

  x.margin.a <- ifelse(x.margin.a > 1, x.margin.a, 1)
  x.margin.a2 <- ifelse(x.margin.a2 > 1, x.margin.a2, 1)
  x.margin.a <- log(x.margin.a)
  x.margin.a2 <- log(x.margin.a2)
  x.margin.a2 <- ifelse(x.margin.a2 > 1, x.margin.a2, 1)

  # x.margin.a <- dist.i.a*nrow(data.other) - dist.o.a*nrow(data.id)# + dist.o.b*nrow(data.id)
  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*x.margin.a2
  # x.margin.adj <- x.margin.a

  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "-")
  return(x.split)
}


get_split_neg_neg  <- function(scrna, gene1, gene2, id, step = 0.01, assay = "RNA", slot = "data") {
  scrna@active.assay <- assay
  data.mat.surf <- Seurat::FetchData(scrna, vars = c(gene1, gene2, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp1", "exp2", "id")
  data.mat.surf$exp1 <- normalize(data.mat.surf$exp1)
  data.mat.surf$exp2 <- normalize(data.mat.surf$exp2)

  sp1 = boxplot.stats(data.mat.surf$exp1)$out
  sp2 = boxplot.stats(data.mat.surf$exp2)$out

  `%notin%` <- Negate(`%in%`)
  data.mat.surf <- data.mat.surf[data.mat.surf$exp2 %notin% sp2 & data.mat.surf$exp1 %notin% sp1, ]

  if(max(data.mat.surf$exp2) == 0 | max(data.mat.surf$exp1) == 0){
    x.split <- data.frame(gene1, gene2, 0, 0, 0, 0, 0, 0, "-", "-")
    return(x.split)
  }

  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  # x.max <- 0
  # x.min <- 1
  # x.num <- 1/step
  # x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))

  # size.factor <- nrow(data.mat.surf)/nrow(data.id)
  size.factor <- nrow(data.mat.surf)/nrow(data.other)
  # size.factor <- 3*nrow(data.id)/nrow(data.other)
  # size.factor <- ifelse(size.factor < 3, 3, size.factor)

  for (x.val in seq(0.01, 1, step)) {
    for (y.val in seq(0.01, 1, step)) {
      data.i.a <- data.id[data.id$exp1 < x.val & data.id$exp2 < y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 < x.val & data.other$exp2 < y.val, ]
      fp <- nrow(data.o.a)
      # data.i.b <- data.id[data.id$exp1 > x.val | data.id$exp2 > y.val, ]
      # fn <- nrow(data.i.b)
      # data.o.b <- data.other[data.other$exp1 > x.val | data.other$exp2 > y.val, ]
      # tn <- nrow(data.o.b)

      # x.margin <- tp*size.factor - fp
      x.margin <- tp - fp*size.factor
      # x.margin <- tp*nrow(data.other) - fp*nrow(data.id)
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
  # data.i.b <- data.id[data.id$exp1 > x.val | data.id$exp2 > y.val, ]
  # fn <- nrow(data.i.b)
  data.o.b <- data.other[data.other$exp1 >= x.val | data.other$exp2 >= y.val, ]
  tn <- nrow(data.o.b)

  # dist.o.b = 0
  # if (nrow(data.o.b) >= 1) {
  #   for (i in 1:nrow(data.o.b)) {
  #     if (data.o.b[i,]$exp1 < x.val) {
  #       if (data.o.b[i,]$exp2 > y.val) {
  #         dist.o.b = dist.o.b + sqrt((data.o.b[i,]$exp1 - x.val)**2 + (data.o.b[i,]$exp2 - y.val)**2)
  #       }else{
  #         dist.o.b = dist.o.b + (x.val - data.o.b[i,]$exp1)
  #       }
  #     }else{
  #       dist.o.b = dist.o.b + (data.o.b[i,]$exp2 - y.val)
  #     }
  #   }
  # }
  # dist.i.a <- sum(apply(data.i.a[,c(1,2)], 1, dist_p,  a2 = c(x.val, y.val)))
  # dist.o.a <- sum(apply(data.o.a[,c(1,2)], 1, dist_p,  a2 = c(x.val, y.val)))
  # x.margin.a <- dist.i.a - dist.o.a + dist.o.b
  # x.margin.a <- dist.i.a*nrow(data.other) - dist.o.a*nrow(data.id)# + dist.o.b*nrow(data.id)
  # x.margin.a <- sum((x.val - data.i.a$exp1))*nrow(data.o.b) + sum((y.val - data.i.a$exp2))*nrow(data.o.b) -#+
    # sum((data.o.b$exp1 - x.val))*nrow(data.i.a) + sum((data.o.b$exp2 - y.val))*nrow(data.i.a) -
    # sum((x.val - data.o.a$exp1))*nrow(data.o.b) - sum((y.val - data.o.a$exp2))*nrow(data.o.b)

  x.margin.a <- sum((x.val - data.i.a$exp1))/nrow(data.o.a) + sum((y.val - data.i.a$exp2))*nrow(data.o.a) -#+
    # sum((data.o.b$exp1 - x.val))*nrow(data.id) + sum((data.o.b$exp2 - y.val))*nrow(data.id) -
    sum((x.val - data.o.a$exp1))*nrow(data.i.a) - sum((y.val - data.o.a$exp2))*nrow(data.i.a)

  # x.margin.a2 <- sum((data.o.b$exp1 - x.val))*nrow(data.id) + sum((data.o.b$exp2 - y.val))*nrow(data.id)

  x.margin.a <- ifelse(x.margin.a > 1, x.margin.a, 1)
  # x.margin.a2 <- ifelse(x.margin.a2 > 1, x.margin.a2, 1)
  x.margin.a <- log(x.margin.a)
  # x.margin.a2 <- log(x.margin.a2)
  # x.margin.a2 <- ifelse(x.margin.a2 > 1, x.margin.a2, 1)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*x.margin.a2
  # x.margin.adj <- x.margin.a
  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "-", "-")
  return(x.split)
}

get_split_pos_pos <- function(scrna, gene1, gene2, id, step = 0.01, assay = "RNA", slot = "data") {
  DefaultAssay(scrna) <- assay
  # data.mat.surf <- FetchData(scrna, vars = c(gene1, gene2, "ident"), slot = slot)
  # colnames(data.mat.surf) <- c("exp1", "exp2", "id")
  # data.mat.surf$exp1 <- normalize(data.mat.surf$exp1)
  # data.mat.surf$exp2 <- normalize(data.mat.surf$exp2)
  # gene.prauc <- data.frame(x.val <- c(),
  #                          margin <- c())
  #
  # data.id <- data.mat.surf[data.mat.surf$id == id,]
  # data.other <- data.mat.surf[data.mat.surf$id != id,]
  # data.id.l <- nrow(data.id)
  # data.o.l <- nrow(data.other)
  # data.all.l <- nrow(data.mat.surf)
  #
  # for (x.val in seq(0, 1, step)) {
  #   for (y.val in seq(0, 1, step)) {
  #     data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
  #     tp <- nrow(data.i.a)
  #     data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
  #     fp <- nrow(data.o.a)
  #     data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
  #     fn <- nrow(data.i.b)
  #     data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
  #     tn <- nrow(data.o.b)
  #
  #     x.margin <- tp - fp
  #     x.margin <- x.margin/data.all.l
  #     de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
  #     gene.prauc <- rbind(gene.prauc, de)
  #   }
  # }
  # gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  # x.split <- gene.prauc[1,]
  # x.val <- x.split[,3]
  # y.val <- x.split[,4]
  # x.margin <- x.split[,5]
  # data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
  # tp <- nrow(data.i.a)
  # data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
  # fp <- nrow(data.o.a)
  # data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
  # fn <- nrow(data.i.b)
  # data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
  # tn <- nrow(data.o.b)
  #
  # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
  #   sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
  #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
  # x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)
  #
  # for (x.val in seq(0.01, 1, 0.1)) {
  #   for (y.val in seq(0.01, 1, 0.1)) {
  #     print(x.val, y.val)
  #     data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
  #     tp <- nrow(data.i.a)
  #     data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
  #     fp <- nrow(data.o.a)
  #
  #     x.margin <- tp - fp
  #     x.margin <- x.margin/data.all.l
  #
  #     de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
  #     gene.prauc <- rbind(gene.prauc, de)
  #   }
  # }
  data.mat.surf <- FetchData(scrna, vars = c(gene1, gene2, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp1", "exp2", "id")
  data.mat.surf$exp1 <- normalize(data.mat.surf$exp1)
  data.mat.surf$exp2 <- normalize(data.mat.surf$exp2)

  sp1 = boxplot.stats(data.mat.surf$exp1)$out
  sp2 = boxplot.stats(data.mat.surf$exp2)$out

  `%notin%` <- Negate(`%in%`)
  data.mat.surf <- data.mat.surf[data.mat.surf$exp2 %notin% sp2 & data.mat.surf$exp1 %notin% sp1, ]

  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  # size.factor <- nrow(data.mat.surf)/nrow(data.id)
  size.factor <- nrow(data.mat.surf)/nrow(data.other)
  # size.factor <- ifelse(size.factor < 3, 3, size.factor)

  for (x.val in seq(0.01, 1, step)) {
    for (y.val in seq(0.01, 1, step)) {
      data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
      tp <- nrow(data.i.a)
      data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
      fp <- nrow(data.o.a)
      # x.margin <- tp*size.factor - fp#*size.factor
      x.margin <- tp - fp*size.factor
      # x.margin <- tp*nrow(data.other) - fp*nrow(data.id)
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
  data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
  tp <- nrow(data.i.a)
  data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
  fp <- nrow(data.o.a)
  data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
  tn <- nrow(data.o.b)
  # dist.o.b = 0
  # for (i in 1:nrow(data.o.b)) {
  #   if (data.o.b[i,]$exp1 < x.val) {
  #     if (data.o.b[i,]$exp2 < y.val) {
  #       dist.o.b = dist.o.b + sqrt((data.o.b[i,]$exp1 - x.val)**2 + (data.o.b[i,]$exp2 - y.val)**2)
  #     }else{
  #       dist.o.b = dist.o.b + (x.val - data.o.b[i,]$exp1)
  #     }
  #   }else{
  #     dist.o.b = dist.o.b + (y.val - data.o.b[i,]$exp2)
  #   }
  # }
  #
  # dist.i.a <- sum(apply(data.i.a[,c(1,2)], 1, dist_p,  a2 = c(x.val, y.val)))
  # dist.o.a <- sum(apply(data.o.a[,c(1,2)], 1, dist_p,  a2 = c(x.val, y.val)))
  # x.margin.a <- dist.i.a - dist.o.a + dist.o.b

  # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) -#+
  #   #sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
  #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)

  x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.a) + sum((data.i.a$exp2 - y.val))*nrow(data.o.a) -#+
    # sum((x.val -data.o.b$exp1))*nrow(data.id) + sum((y.val -data.o.b$exp2))*nrow(data.id) -
    sum((data.o.a$exp1 - x.val))*nrow(data.i.a) - sum((data.o.a$exp2 - y.val))*nrow(data.i.a)

  # x.margin.a <- x.margin.a/(nrow(data.i.a) + nrow(data.o.a))

  # x.margin.a2 <- sum((x.val -data.o.b$exp1))*nrow(data.id) + sum((y.val -data.o.b$exp2))*nrow(data.id)

  x.margin.a <- ifelse(x.margin.a > 1, x.margin.a, 1)
  # x.margin.a2 <- ifelse(x.margin.a2 > 1, x.margin.a2, 1)
  x.margin.a <- log(x.margin.a)
  # x.margin.a2 <- log(x.margin.a2)
  # x.margin.a2 <- ifelse(x.margin.a2 > 1, x.margin.a2, 1)


  # x.margin.a <- dist.i.a*nrow(data.other) - dist.o.a*nrow(data.id)# + dist.o.b*nrow(data.id)
  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*x.margin.a2
  # x.margin.adj <- x.margin.a

  x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+", x.margin.a)
  return(x.split)
}



#' Calculate gene positive specific score
#'
#' Calculate gene specific score in positive direction
#' @param scrna obj to use
#' @param id interested cell id
#' @param gene gene to use
#' @param step quantile step
#' @param assay which assay to use in obj, default RNA
#' @param slot which slot to use, default data
#' @export
#' @examples
#' @return seurat obj with re-set Idents
#'

get_gene_score_pos <- function(scrna, gene, id, step = 0.01, assay = "RNA", slot = "data", pseudo.count = 0.01) {
  scrna@active.assay <- assay
  data.mat.surf <- Seurat::FetchData(scrna, vars = c(gene, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp", "id")
  data.mat.surf$exp <- normalize(data.mat.surf$exp)
  data.mat.surf$exp <- round(data.mat.surf$exp, 3)

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

  x.factor <- (mean(data.id$exp) + pseudo.count)/(mean(data.other$exp) + pseudo.count)
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
      sum((x.val - data.o.b$exp))*nrow(data.i.a) -
      sum((data.o.a$exp - x.val))*nrow(data.i.b)
    x.margin <- x.margin/data.all.l

    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]

  tp <- sum(data.id$exp > x.val)
  fp <- sum(data.other$exp > x.val)
  fn <- sum(data.id$exp <= x.val)
  tn <- sum(data.other$exp <= x.val)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor**2)
  x.val <- get_original_alpha(scrna, gene = gene, x.alpha = x.val, assay = assay, slot = slot)
  x.split <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, tn, fn, "+")

  # rownames(x.split) <- paste(gene)
  return(x.split)
}


#' Calculate gene negative specific score
#'
#' Calculate gene specific score in negative direction
#' @param scrna obj to use
#' @param id interested cell id
#' @param gene gene to use
#' @param step quantile step
#' @param assay which assay to use in obj, default RNA
#' @param slot which slot to use, default data
#' @export
#' @examples
#' @return seurat obj with re-set Idents
#'
#'

get_gene_score_neg <- function(scrna, gene, id, step = 0.01, assay = "RNA", slot = "data", pseudo.count = 0.01) {
  scrna@active.assay <- assay
  data.mat.surf <- FetchData(scrna, vars = c(gene, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp", "id")
  data.mat.surf$exp <- normalize(data.mat.surf$exp)
  data.mat.surf$exp <- round(data.mat.surf$exp, 3)

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

  x.factor <- (mean(data.other$exp) + pseudo.count)/(mean(data.id$exp) + pseudo.count)

  for (x.val in seq(0, 1, step)) {
    data.i.a <- data.id[data.id$exp < x.val, ]
    data.o.a <- data.other[data.other$exp < x.val, ]
    data.i.b <- data.id[data.id$exp >= x.val, ]
    data.o.b <- data.other[data.other$exp >= x.val, ]

    x.margin <- sum((data.i.a$exp - x.val))*nrow(data.o.b) +
      sum((x.val -data.o.b$exp))*nrow(data.i.a) -
      sum((data.o.a$exp - x.val))*nrow(data.i.b)
    x.margin <- x.margin/data.all.l

    x.margin <- 0 - x.margin

    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]
  tp <- sum(data.id$exp < x.val)
  fp <- sum(data.other$exp < x.val)
  tn <- sum(data.other$exp >= x.val)
  fn <- sum(data.id$exp >= x.val)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor**2)
  x.val <- get_original_alpha(scrna, gene = gene, x.alpha = x.val, assay = assay)
  x.split <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, tn, fn, "-")

  return(x.split)
}


#' Calculate gene negative specific score
#'
#' Calculate gene specific score in negative direction
#' @param scrna obj to use
#' @param gene gene to use
#' @param assay which assay to use in obj, default RNA
#' @param slot which slot to use, default data
#' @export
#' @examples
#' @return alpha before normalization
#'
get_original_alpha <- function(scrna, gene, x.alpha, assay = "RNA", slot = "data"){
  Seurat::DefaultAssay(scrna) <- assay
  g.exprs <- FetchData(scrna, gene, slot = slot)
  g.max <- max(g.exprs[,1])
  g.min <- min(g.exprs[,1])
  alpha.ori <- x.alpha*g.max + g.min
  return(alpha.ori)
}



marker_stepbystep.2 <- function(scrna, cellgroup, depth = 2, geneset = NULL, step = 0.01, assay = "RNA", slot = "data", max.recall = T){
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
        de <- data.frame(get_gene_score_pos(scrna, genes, cellgroup, step, assay  = assay, slot = slot))
        names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
        gene.prauc <- rbind(gene.prauc, de)
      }
      for (genes in rownames(markers.neg)) {
        de <- data.frame(get_gene_score_neg(scrna, genes, cellgroup, step, assay  = assay, slot = slot))
        names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
        gene.prauc <- rbind(gene.prauc, de)
      }
      gene.prauc <- gene.prauc[order(gene.prauc$x.margin.adj, decreasing = T),]

      df <- gene.prauc[1,]
      if (df$direction == "+") {
        split.step <- get_split_pos(scrna, df$gene, cellgroup, assay = assay, slot = slot, max.recall = T)
      }else{
        split.step <- get_split_neg(scrna, df$gene, cellgroup,  assay = assay, slot = slot, max.recall = T)
      }

      names(split.step) <- c("Gene", "split.value", "direction")

      df.split <- rbind(df.split, split.step)
      left.cells <- GetCellNames(scrna = scrna, gene = df$gene, value = df.split$split.value, direction = df.split$direction, assay = assay, slot = slot)
      # left.cells <- GetCellNames(scrna = scrna, gene = df$gene, value = 0, direction = df.split$direction)
      scrna <- subset(scrna, cells = left.cells)
      geneset <- setdiff(geneset, df$gene)
    }
  )
  return(df.split)
}


GetCellNames <- function(scrna, gene, value, direction, assay = "RNA", slot = "data"){
  # if (do.magic) {
  #   df <- data.frame(exp = scrna@assays$MAGIC_RNA@data[gene,])
  # }else{
  #   df <- data.frame(exp = scrna@assays$RNA@data[gene,])
  # }
  scrna@active.assay <- assay
  df <- Seurat::FetchData(scrna, vars = c(gene), slot = slot)

  colnames(df) <- "exp"

  if (direction == "+") {
    df <- subset(df, exp > value)
  }else{
    df <- subset(df, exp < value)
  }
  return(rownames(df))
}

get_split_pos <- function(scrna, gene, id, step = 0.01, assay = "RNA", slot = "data", max.recall = F) {
  scrna@active.assay <- assay
  df <- FetchData(scrna, vars = c(gene, "ident"))
  colnames(df) <- c("exp", "id")
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- df[df$id == id,]
  data.other <- df[df$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)
  x.min <- 0
  x.num <- 1/step
  x.factor <- as.numeric(length(df$id)/length(data.id$id))
  x.factor <- ifelse(x.factor > 5, x.factor, 5)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    tp <- sum(data.id$exp > x.val)
    fp <- sum(data.other$exp > x.val)
    tn <- sum(data.other$exp <= x.val)
    fn <- sum(data.id$exp <= x.val)

    if (max.recall) {
      x.margin <- tn - fn*x.factor
    }else{
      x.margin <- tn - fn#*x.factor
    }
    # if (length(data.id$id) <= length(data.other$id)) {
    #   x.margin <- tn - fn*x.factor
    # }
    # else{
    #   x.margin <- tn - fn*x.factor
    # }
    de <- data.frame(gene, x.val, x.margin, "+")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  colnames(gene.prauc) <- c("Gene", "split.value", "split.margin", "direction")
  x.split <- gene.prauc[1,c(1,2,4)]
  rownames(x.split) <- paste(gene)
  return(x.split)
}

get_split_neg <- function(scrna, gene, id, step = 0.01, assay = "RNA", slot = "data", max.recall = F) {
  scrna@active.assay <- assay
  df <- FetchData(scrna, vars = c(gene, "ident"))
  colnames(df) <- c("exp", "id")
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- df[df$id == id,]
  data.other <- df[df$id != id,]

  x.max <- max(df$exp)
  # x.min <- min(df[df$exp != 0, ]$exp)
  x.min <- min(df$exp)
  x.min <- 0
  x.num <- 1/step

  x.factor <- as.numeric(length(df$id)/length(data.id$id))
  x.factor <- ifelse(x.factor > 5, x.factor, 5)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/x.num)) {
    tp <- sum(data.id$exp <= x.val)
    fp <- sum(data.other$exp <= x.val)
    tn <- sum(data.other$exp > x.val)
    fn <- sum(data.id$exp > x.val)

    if (max.recall) {
      x.margin <- tn - fn*x.factor
    }else{
      x.margin <- tn - fn#*x.factor
    }

    # if (length(data.id$id) <= length(data.other$id)) {
    #   x.margin <- tn - fn*x.factor
    # }
    # else{
    #   x.margin <- tn - fn*x.factor
    # }
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

# get_split_pos_pos <- function(scrna, gene1, gene2, id, step = 0.01, assay = "RNA", slot = "data") {
#   scrna@active.assay <- assay
#   data.mat.surf <- FetchData(scrna, vars = c(gene1, gene2, "ident"), slot = slot)
#   colnames(data.mat.surf) <- c("exp1", "exp2", "id")
#   data.mat.surf$exp1 <- normalize(data.mat.surf$exp1)
#   data.mat.surf$exp2 <- normalize(data.mat.surf$exp2)
#   gene.prauc <- data.frame(x.val <- c(),
#                            margin <- c())
#
#   data.id <- data.mat.surf[data.mat.surf$id == id,]
#   data.other <- data.mat.surf[data.mat.surf$id != id,]
#   data.id.l <- nrow(data.id)
#   data.o.l <- nrow(data.other)
#   data.all.l <- nrow(data.mat.surf)
#
#   x.max <- 0
#   x.min <- 1
#   x.num <- 1/step
#
#   x.size.factor <- 10*(length(data.id$exp1)/length(data.other$exp1))
#
#   for (x.val in seq(0, 1, step)) {
#     for (y.val in seq(0, 1, step)) {
#       data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
#       tp <- nrow(data.i.a)
#       data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
#       fp <- nrow(data.o.a)
#       data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
#       fn <- nrow(data.i.b)
#       data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
#       tn <- nrow(data.o.b)
#
#       # x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
#       #   sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
#       #   sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
#       x.margin <- tp - fp
#       x.margin <- x.margin/data.all.l
#
#       # x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)
#
#       # de <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+")
#       de <- data.frame(gene1, gene2, x.val, y.val, x.margin)
#       gene.prauc <- rbind(gene.prauc, de)
#     }
#   }
#   gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
#   x.split <- gene.prauc[1,]
#   x.val <- x.split[,3]
#   y.val <- x.split[,4]
#   x.margin <- x.split[,5]
#   data.i.a <- data.id[data.id$exp1 > x.val & data.id$exp2 > y.val, ]
#   tp <- nrow(data.i.a)
#   data.o.a <- data.other[data.other$exp1 > x.val & data.other$exp2 > y.val, ]
#   fp <- nrow(data.o.a)
#   data.i.b <- data.id[data.id$exp1 <= x.val | data.id$exp2 <= y.val, ]
#   fn <- nrow(data.i.b)
#   data.o.b <- data.other[data.other$exp1 <= x.val | data.other$exp2 <= y.val, ]
#   tn <- nrow(data.o.b)
#   x.margin.a <- sum((data.i.a$exp1 - x.val))*nrow(data.o.b) + sum((data.i.a$exp2 - y.val))*nrow(data.o.b) +
#     sum((x.val -data.o.b$exp1))*nrow(data.i.a) + sum((y.val -data.o.b$exp2))*nrow(data.i.a) -
#     sum((data.o.a$exp1 - x.val))*nrow(data.o.b) - sum((data.o.a$exp2 - y.val))*nrow(data.o.b)
#   x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*x.margin.a#*((x1.factor*x2.factor)**2)
#   x.split <- data.frame(gene1, gene2, x.val, y.val, x.margin, x.margin.adj, tp, fp, "+", "+")
#   return(x.split)
# }

#' Firstup gene names
#'
#' @param x character
#' @return up case first letter
#' @export
#'
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


#' make id for plotting
#'
#' Reset active.ident of obj for plotting.
#' @param scrna obj to use
#' @param id interested cell id
#' @param panel ID panel in metadata of obj
#' @export
#' @examples
#' @return seurat obj with re-set Idents
#'

makeid <- function(scrna, id, panel = "seurat_clusters"){
  levels(scrna@active.ident) <- c(levels(scrna@active.ident), "other")
  scrna@active.ident[scrna@active.ident != id] <- 'other'
  return(scrna)
}

#' Min-Max Normalization
#'
#' Normalize scRNA sequencing data.
#' @param x Gene expression
#' @export
#' @examples
#' @return Normalized expression data
#'

normalize <- function(x, ...) {
  return((x - min(x, ...)) /(max(x, ...) - min(x, ...)))
}


grid_vln_plot <- function(scrna, df.split, cellgroup, ncol = 2, assay = "RNA", slot = "data"){
  x.len <- length(df.split$gene)
  plots = list()
  # scrna.id <- scrna@active.ident
  scrna <- makeid(scrna, cellgroup)
  scrna@active.assay = assay
  scrna.left <- scrna
  n = 1
  for (gene in df.split$Gene) {
    plots[[gene]] <- VlnPlot(scrna.left, features = gene) + geom_hline(yintercept=df.split[df.split$Gene == gene, ]$split.value, linetype = 1)
    # plot(VlnPlot(scrna.left, features = gene) + geom_hline(yintercept=df.split[df.split$Gene == gene, ]$split.value, linetype = 1))
    left.cells <- GetCellNames(scrna.left, gene, df.split[df.split$Gene == gene, ]$split.value, df.split[df.split$Gene == gene, ]$direction, assay = assay, slot = slot)
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

#' Detect single markers of specific cell groups
#' @param scrna seurat obj to be used
#' @param geneset custom genes to test, all genes in obj to use if NULL
#' @param id interested cell group
#' @param category indicate the ICC/IHC/Flow, if NULL use all genes.
#' @param assay indicate the assay to compute
#' @param slot indicate the slot of data to compute
#' @param step quantile steps
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
#' @param use.all Don't do any filter on input geneset
#' @return list of markers performance
#' @export
#'
#'
Detect_single_marker <- function(scrna, id, step = 0.1,  slot = "data", category = NULL,
                                 geneset = NULL, assay = "RNA", do.fast = F, min.pct = 0.1,
                                 use.all = F, do.f1score = F, pseudo.count = 0.01, min.tnr = 0.8){
  if (is.null(geneset)) {
    geneset <- rownames(scrna[[assay]])
  }else{
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
  }
  if (use.all) {
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
    de <- Seurat::FoldChange(scrna, ident.1 = id, features = geneset)
  }
  if (!use.all) {
    if (do.fast) {
      de <- Seurat::FindMarkers(scrna, ident.1 = id, features = geneset)
      de <- subset(de, p_val_adj < 0.05)
    }else{
      de <- Seurat::FoldChange(scrna, ident.1 = id, features = geneset)
      de <- de[de$pct.1 > min.pct, ]
    }
  }

  de$gene <- rownames(de)
  markers <- de

  gene.rank.list <- data.frame()

  scrna@active.assay <- assay
  exprs.matrix <- Seurat::FetchData(scrna, vars = c(rownames(markers), "ident"), slot = slot)

  gene.list <- rownames(markers)
  # gene.list <- toupper(rownames(markers))

  for (genes in gene.list) {
    de <- data.frame(get_gene_score(exprs.matrix[, c(genes, "ident")], gene = genes, id = id,
                                    step = step, pseudo.count = pseudo.count, min.tnr = min.tnr))
    names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "TP", "FP", "TN", "FN", "direction")
    gene.rank.list <- rbind(gene.rank.list, de)
  }

  gene.rank.list <- gene.rank.list[order(gene.rank.list$x.margin.adj, decreasing = T),]
  rownames(gene.rank.list) <- seq(nrow(gene.rank.list))
  return(gene.rank.list)
}


#' Calculate gene specific score in both direction
#' @param exprs.matrix exprs matrix with ident panel
#' @param id interested cell id
#' @param gene gene to use
#' @param step quantile step
#' @export
#' @examples
#' @return gene score in both direction
#'
#'
get_gene_score <- function(exprs.matrix, id, gene, step = 0.01, pseudo.count = 0.01, min.tnr = 0.75){
  data.mat.surf <- exprs.matrix
  colnames(data.mat.surf) <- c("exp", "id")
  exprs.max <- max(data.mat.surf$exp)
  exprs.min <- min(data.mat.surf$exp)
  data.mat.surf$exp <- normalize(data.mat.surf$exp)
  data.mat.surf$exp <- round(data.mat.surf$exp, 3)

  gene.prauc.p <- data.frame(x.val <- c(),
                             margin <- c())
  gene.prauc.n <- data.frame(x.val <- c(),
                             margin <- c())

  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)

  x.factor.p <- (mean(data.id$exp) + pseudo.count)/(mean(data.other$exp) + pseudo.count)
  x.factor.n <- (mean(data.other$exp) + pseudo.count)/(mean(data.id$exp) + pseudo.count)

  for (x.val in seq(0.01, 0.99, step)) {
    data.i.a <- data.id[data.id$exp > x.val, ]
    data.o.a <- data.other[data.other$exp > x.val, ]
    data.i.b <- data.id[data.id$exp <= x.val, ]
    data.o.b <- data.other[data.other$exp <= x.val, ]
    x.margin.p <- sum((data.i.a$exp - x.val))*nrow(data.o.b) +
      sum((x.val - data.o.b$exp))*nrow(data.i.a) -
      sum((data.o.a$exp - x.val))*nrow(data.i.b)
    x.margin.p <- x.margin.p/data.all.l
    de.p <- data.frame(x.val, x.margin.p)
    gene.prauc.p <- rbind(gene.prauc.p, de.p)

    data.i.a <- data.id[data.id$exp < x.val, ]
    data.o.a <- data.other[data.other$exp < x.val, ]
    data.i.b <- data.id[data.id$exp >= x.val, ]
    data.o.b <- data.other[data.other$exp >= x.val, ]

    x.margin.n <- sum((data.i.a$exp - x.val))*nrow(data.o.b) +
      sum((x.val -data.o.b$exp))*nrow(data.i.a) -
      sum((data.o.a$exp - x.val))*nrow(data.i.b)
    x.margin.n <- x.margin.n/data.all.l

    x.margin.n <- 0 - x.margin.n

    de.n <- data.frame(x.val, x.margin.n)
    gene.prauc.n <- rbind(gene.prauc.n, de.n)
  }

  gene.prauc.p <- gene.prauc.p[order(gene.prauc.p$x.margin, decreasing = T),]
  gene.prauc.n <- gene.prauc.n[order(gene.prauc.n$x.margin, decreasing = T),]

  # pos
  x.split <- gene.prauc.p[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]
  tp <- sum(data.id$exp > x.val)
  fp <- sum(data.other$exp > x.val)
  fn <- sum(data.id$exp <= x.val)
  tn <- sum(data.other$exp <= x.val)
  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor.p**2)
  x.val <- x.val*exprs.max + exprs.min
  x.split.p <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, tn, fn, "+")

  x.split <- gene.prauc.n[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]
  tp <- sum(data.id$exp < x.val)
  fp <- sum(data.other$exp < x.val)
  tn <- sum(data.other$exp >= x.val)
  fn <- sum(data.id$exp >= x.val)

  x.margin.adj <- x.margin*(tp/data.id.l)*(tn/data.o.l)*(x.factor.n**2)
  x.val <- x.val*exprs.max + exprs.min
  if (tn/data.o.l < min.tnr) {
    x.margin.adj <- 0
  }
  x.split.n <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, tn, fn, "-")
  return(rbind(x.split.p, x.split.n))
}


#' Get Antibody information table
#' @param markers.list makers list from Detect single markers function
#' @return list of markers performance
#' @export
#'
get_antibody <- function(markers.list){
  markers.list$gene <- toupper(markers.list$gene)
  markers.list$antibody <- NULL
  markers.list$x.margin <- round(markers.list$x.margin, 2)
  markers.list$x.margin.adj <- round(markers.list$x.margin.adj, 2)
  markers.list$antibody <- ifelse(is.null(markers.list$antibody), IHC[match(markers.list$gene, IHC$Gene), ]$Antibody.RRID, markers.list$antibody)
  markers.list$antibody <- ifelse(is.null(markers.list$antibody), ICC[match(markers.list$gene, ICC$Gene), ]$Antibody.RRID, markers.list$antibody)


  markers.list <- markers.list[, -3]
  markers.list[,2] <- round(markers.list[,2], 3)
  colnames(markers.list)[2] <- "Alpha"
  colnames(markers.list)[3] <- "Score"
  # msc.markers <- distinct(msc.markers, gene, .keep_all = TRUE)

  # markers.list <- markers.list[!is.null(markers.list$antibody), ]
  markers.list <- markers.list[!is.null(markers.list$antibody),]

  DT::datatable(markers.list,
                extensions = 'Buttons',
                options = list(dom = 'Bfrtip',
                               buttons = c('copy',
                                           'csv',
                                           'excel',
                                           'pdf', 'print'),
                               autoWidth = FALSE),
                escape = F)
}


#' Get cell names which meet cutoff
#' @param scrna seurat obj to be used
#' @param geneset custom genes to test, all genes in obj to use if NULL
#' @param id interested cell group
#' @param category indicate the ICC/IHC/Flow, if NULL use all genes.
#' @param assay indicate the assay to compute
#' @param slot indicate the slot of data to compute
#' @param step quantile steps
#' @param min.pct only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1
#' @param use.all Don't do any filter on input geneset
#' @return list of markers performance
#'
#'
GetGenes <- function(scrna, slot = "data", category = NULL, geneset = NULL,
                     assay = "RNA", min.pct = 0.1, use.all = F, species = "Human"){
  if (category == "IHC") {
    genes.touse <- IHC$Gene
  }else if (category == "ICC") {
    genes.touse <- ICC$Gene
  }else if (category == "Flow") {
    genes.touse <- surface.genes
  }

  if (!is.null(category) & is.null(geneset)) {
    print("Use self define genes should set the category to NULL")
    return(NULL)
  }

  if (is.null(category) & is.null(geneset)) {
    genes.touse <- rownames(scrna[[assay]])
  }else if (is.null(category) & !is.null(geneset)){
    genes.touse <- intersect(rownames(scrna[[assay]]), geneset)
  }

  if (length(genes.touse) == 0) {
    print("Find none of genes in the obj")
    return(NULL)
  }

  if (use.all) {
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
    de <- Seurat::FoldChange(scrna, ident.1 = id, features = geneset)
  }
  if (!use.all) {
    if (do.fast) {
      de <- Seurat::FindMarkers(scrna, ident.1 = id, features = geneset)
      de <- subset(de, p_val_adj < 0.05)
    }else{
      de <- Seurat::FoldChange(scrna, ident.1 = id, features = geneset)
      de <- de[de$pct.1 > min.pct, ]
    }
  }


  de$gene <- rownames(de)
  # markers.pos <- subset(de, avg_log2FC > 0)
  # markers.neg <- subset(de, avg_log2FC < 0)
  markers.pos <- de
  markers.neg <- de

  gene.rank.list <- data.frame()

  for (genes in rownames(markers.pos)) {
    de <- data.frame(get_gene_score_pos(scrna, genes, id, step, assay  = assay, slot = slot))
    names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
    gene.rank.list <- rbind(gene.rank.list, de)
  }
  for (genes in rownames(markers.neg)) {
    de <- data.frame(get_gene_score_neg(scrna, genes, id, step, assay  = assay, slot = slot))
    names(de)<-c("gene", "split.value","x.margin", "x.margin.adj", "tp", "fp", "direction")
    gene.rank.list <- rbind(gene.rank.list, de)
  }
  gene.rank.list <- gene.rank.list[order(gene.rank.list$x.margin.adj, decreasing = T),]
  rownames(gene.rank.list) <- seq(nrow(gene.rank.list))
  return(gene.rank.list)
}


Detect_combine_markers <- function(scrna, id, step = 0.1, category = NULL,
                                   geneset = NULL, assay = "RNA", slot = "data",
                                   do.fast = F, min.pct = 0.1, use.all = F, pseudo.count = 0.01){
  if (is.null(geneset)) {
    geneset <- rownames(scrna[[assay]])
  }else{
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
  }

  if (use.all) {
    geneset <- intersect(rownames(scrna[[assay]]), geneset)
    de <- Seurat::FoldChange(scrna, ident.1 = id, features = geneset)
    de$gene <- rownames(de)
  }
  if (!use.all) {
    if (do.fast) {
      de <- Seurat::FindMarkers(scrna, ident.1 = id, features = geneset)
      de <- subset(de, p_val_adj < 0.05)
    }else{
      de <- Seurat::FoldChange(scrna, ident.1 = id, features = geneset)
      de$gene <- rownames(de)
      de <- de[de$pct.1 > min.pct, ]
    }
  }

  de.genes.p <- subset(de, avg_log2FC > 0)
  de.genes.n <- subset(de, avg_log2FC < 0)

  gene.combn <- t(combn(de.genes.p$gene, 2))
  gene.combn.prauc <- data.frame()

  star_time <- Sys.time()
  for (gene.iter in 1:nrow(gene.combn)) {
    genes <- gene.combn[gene.iter,]
    gene1 <- genes[1]
    gene2 <- genes[2]
    de <- get_split_pos_pos(scrna, gene1,	gene2, id, step = step, assay = assay, slot = slot)
    colnames(de) <- c("Gene1", "Gene2", "Gene1.split", "Gene2.split", "x.margin", "x.margin.adj", "TP", "FP", "Direction1", "Direction2")
    gene.combn.prauc <- rbind(gene.combn.prauc, de)
  }

  for (gene1 in de.genes.p$gene) {
    for (gene2 in de.genes.n$gene) {
      de <- get_split_pos_neg(scrna, gene1,	gene2, id, step = step, assay = assay, slot = slot)
      colnames(de) <- c("Gene1", "Gene2", "Gene1.split", "Gene2.split", "x.margin", "x.margin.adj", "TP", "FP", "Direction1", "Direction2")
      gene.combn.prauc <- rbind(gene.combn.prauc, de)
    }
  }

  gene.combn <- t(combn(de.genes.n$gene, 2))

  for (gene.iter in 1:nrow(gene.combn)) {
    genes <- gene.combn[gene.iter,]
    gene1 <- genes[1]
    gene2 <- genes[2]
    de <- get_split_neg_neg(scrna, gene1,	gene2, id, step = step, assay = assay, slot = slot)
    colnames(de) <- c("Gene1", "Gene2", "Gene1.split", "Gene2.split", "x.margin", "x.margin.adj", "TP", "FP", "Direction1", "Direction2")
    gene.combn.prauc <- rbind(gene.combn.prauc, de)
  }

  end_time <- Sys.time()
  # close(pb)
  run_time <- end_time - star_time
  print(run_time)

  gene.combn.prauc <- gene.combn.prauc[order(gene.combn.prauc$x.margin.adj, decreasing = T),]
  return(gene.combn.prauc)
}



Combine_markers_point_plot <- function(scrna, c.markers, id, ranking = 1, assay = "RNA", slot = "data", return.obj = F, step = 0.01, max.recall = F){
  scrna@active.assay = assay
  markers.c <- c.markers[ranking, ]
  scrna@active.assay = assay
  df <- Seurat::FetchData(scrna, vars = c(markers.c$Gene1, markers.c$Gene2, "ident"))
  gene1 <- markers.c$Gene1
  gene2 <- markers.c$Gene2
  colnames(df) <- c("exp1", "exp2", "id")
  gene1.direc = markers.c$Direction1
  gene2.direc = markers.c$Direction2
  if (max.recall) {
    if (gene1.direc == "+") {
      x.split <- get_split_pos(scrna, gene1, id, step = step, assay = assay, slot = slot, max.recall = max.recall)$split.value
    }else{
      x.split <- get_split_neg(scrna, gene1, id, step = step, assay = assay, slot = slot, max.recall = max.recall)$split.value
    }
    left.cells <- GetCellNames(scrna, markers.c$Gene1, x.split, markers.c$Direction1, assay = assay, slot = slot)
    scrna.left <- subset(scrna, cells = left.cells)

    if (gene2.direc == "+") {
      y.split <- get_split_pos(scrna.left, gene2, id, step = step, assay = assay, slot = slot, max.recall = max.recall)$split.value
    }else{
      y.split <- get_split_neg(scrna.left, gene2, id, step = step, assay = assay, slot = slot, max.recall = max.recall)$split.value
    }

    left.cells <- GetCellNames(scrna, markers.c$Gene2, y.split, markers.c$Direction2, assay = assay, slot = slot)
    scrna.left <- subset(scrna.left, cells = left.cells)


    left.cells <- GetCellNames(scrna, markers.c$Gene2, y.split, markers.c$Direction2, assay = assay, slot = slot)
    scrna.left <- subset(scrna.left, cells = left.cells)
  }else {
    if (gene1.direc == "+" & gene2.direc == "+") {
      df.c <-get_split_pos_pos(scrna, gene1 = gene1, gene2 = gene2, id = id, assay = assay, slot = slot, step = step)
      colnames(df.c) <- c("Gene1", "Gene2", "Gene1.split", "Gene2.split", "x.margin", "x.margin.adj", "TP", "FP", "Direction1", "Direction2")
    }
    if (gene1.direc == "+" & gene2.direc == "-") {
      df.c <-get_split_pos_neg(scrna, gene1 = gene1, gene2 = gene2, id = id, assay = assay, slot = slot, step = step)
      colnames(df.c) <- c("Gene1", "Gene2", "Gene1.split", "Gene2.split", "x.margin", "x.margin.adj", "TP", "FP", "Direction1", "Direction2")
    }
    if (gene1.direc == "-" & gene2.direc == "-") {
      df.c <-get_split_neg_neg(scrna, gene1 = gene1, gene2 = gene2, id = id, assay = assay, slot = slot, step = step)
      colnames(df.c) <- c("Gene1", "Gene2", "Gene1.split", "Gene2.split", "x.margin", "x.margin.adj", "TP", "FP", "Direction1", "Direction2")
    }

    x.split <- df.c$Gene1.split*max(df$exp1) + min(df$exp1)
    y.split <- df.c$Gene2.split*max(df$exp2) + min(df$exp2)

    left.cells <- GetCellNames(scrna, markers.c$Gene1, x.split, markers.c$Direction1, assay = assay, slot = slot)
    scrna.left <- subset(scrna, cells = left.cells)
    left.cells <- GetCellNames(scrna, markers.c$Gene2, y.split, markers.c$Direction2, assay = assay, slot = slot)
    scrna.left <- subset(scrna.left, cells = left.cells)

  }

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

  g <- g + xlab(paste(markers.c$Gene1, "--Recall: ", x.rec,"; Precision: ", x.pre, sep = "")) + ylab(paste(markers.c$Gene2))
  if (return.obj) {
    return(g)
  }
  print(g)
}

#' RidgePlot of selected genes, with split value
#' @param scrna seurat obj to be used
#' @param genes genes to plot
#' @param id interested cell group
#' @param step quantile steps
#' @param show_split whether to show the split line
#' @return Ridge plot
#' @export
#'
plot_ridge <- function(scrna, id, genes, ncol = 1, step = 0.01, show_split = T, assay = "RNA", slot = "data"){
  df.all <-data.frame()
  df.split <-data.frame()
  for (gene in genes) {
    fc <- Seurat::FoldChange(scrna, ident.1 = id, features = gene)$avg_log2FC
    if (fc > 0) {
      df.g <- get_gene_score_pos(scrna, gene, id, step = 0.01)
    }else{
      df.g <- get_gene_score_neg(scrna, gene, id, step = 0.01)
    }
    df.c <-Seurat::FetchData(scrna, vars = c(gene, "ident"))
    df.c$ident <- ifelse(df.c$ident == id, id, "Other")
    df.c$gene <- paste(gene)
    # rownames(df.c) <- NULL
    colnames(df.c) <- c("Exprs", "Ident", "Gene")
    df.all <- rbind(df.all,df.c)
    df.g <- df.g[,c(1,2)]
    colnames(df.g) <- c("Gene", "Split")
    df.split <- rbind(df.split,df.g)
  }

  df.s <- reshape2::melt(df.all)
  df.s[df.s == -Inf] <- 0
  df.s$Gene <- factor(df.s$Gene, levels = as.character(genes))

  g <- ggplot2::ggplot(df.s, aes(x=value, y=variable, color=Ident, point_color=Ident, fill=Ident)) +
    ggridges::geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
    # scale_y_discrete(expand = c(.01, 0)) +
    # scale_x_continuous(expand = c(0, 0), name = "Expression") +
    scale_fill_manual(values = c("#CB181D80", "#2171B580")) +
    # scale_color_manual(values = c("#CB181D80", "#2171B580"), guide = "none") +
    scale_discrete_manual("point_color", values = c("#CB181D80", "#2171B580"), guide = "none") +
    guides(fill = guide_legend(
      override.aes = list(
        fill = c("#CB181D80", "#2171B580"),
        color = NA, point_color = NA))
    ) +
    ggridges::theme_ridges(grid = FALSE) +
    ylab("") +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )


  g + facet_wrap(~Gene, ncol = 3) +
    geom_vline(data = df.split, aes(xintercept = Split), linetype="dotted",
               color = "red", size=1.5)

}



#' Calculate gene positive specific score with FBeta (Hypergate like approach)
#'
#' Calculate gene FBeta score with positive direction
#' @param scrna obj to use
#' @param id interested cell id
#' @param gene gene to use
#' @param step quantile step
#' @param assay which assay to use in obj, default RNA
#' @param slot which slot to use, default data
#' @param panel ID panel in metadata of obj
#' @export
#' @examples
#' @return seurat obj with re-set Idents
#'

get_gene_score_pos.fbeta <- function(scrna, gene, id, step = 0.01, assay = "RNA", slot = "data") {
  scrna@active.assay <- assay
  data.mat.surf <- Seurat::FetchData(scrna, vars = c(gene, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp", "id")
  data.mat.surf$exp <- normalize(data.mat.surf$exp)

  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)


  for (x.val in seq(0, 1, step)) {
    am <- data.id[data.id$exp > x.val, ]
    cm <- data.other[data.other$exp > x.val, ]
    bm <- data.id[data.id$exp <= x.val, ]
    dm <- data.other[data.other$exp <= x.val, ]

    x.pr <- nrow(am)/(nrow(am) + nrow(cm))
    x.re <- nrow(am)/(nrow(am) + nrow(bm))

    x.margin <- x.pr*x.re/(x.pr + x.re)
    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]

  tp <- sum(data.id$exp > x.val)
  fp <- sum(data.other$exp > x.val)
  fn <- sum(data.id$exp <= x.val)
  tn <- sum(data.other$exp <= x.val)

  x.margin.adj <- x.margin
  x.split <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, "+")
  return(x.split)
}



#' Calculate gene positive specific score with FBeta (Hypergate like approach)
#'
#' Calculate gene FBeta score with negative direction
#' @param scrna obj to use
#' @param id interested cell id
#' @param gene gene to use
#' @param step quantile step
#' @param assay which assay to use in obj, default RNA
#' @param slot which slot to use, default data
#' @param panel ID panel in metadata of obj
#' @export
#' @examples
#' @return seurat obj with re-set Idents
#'

get_gene_score_neg.fbeta <- function(scrna, gene, id, step = 0.01, assay = "RNA", slot = "data") {
  scrna@active.assay <- assay
  data.mat.surf <- Seurat::FetchData(scrna, vars = c(gene, "ident"), slot = slot)
  colnames(data.mat.surf) <- c("exp", "id")
  data.mat.surf$exp <- normalize(data.mat.surf$exp)

  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  data.id.l <- nrow(data.id)
  data.o.l <- nrow(data.other)
  data.all.l <- nrow(data.mat.surf)


  for (x.val in seq(0, 1, step)) {
    am <- data.id[data.id$exp < x.val, ]
    cm <- data.other[data.other$exp < x.val, ]
    bm <- data.id[data.id$exp >= x.val, ]
    dm <- data.other[data.other$exp >= x.val, ]

    x.pr <- nrow(am)/(nrow(am) + nrow(cm))
    x.re <- nrow(am)/(nrow(am) + nrow(bm))

    x.margin <- x.pr*x.re/(x.pr + x.re)
    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  x.val <- x.split[,1]
  x.margin <- x.split[,2]

  tp <- sum(data.id$exp < x.val)
  fp <- sum(data.other$exp < x.val)
  fn <- sum(data.id$exp >= x.val)
  tn <- sum(data.other$exp >= x.val)

  x.margin.adj <- x.margin
  x.split <- data.frame(gene, x.val, x.margin, x.margin.adj, tp, fp, "+")
  return(x.split)
}


