GeneBarPlot <- function(de.data, xlim = NULL, main = NULL) {
  if("avg_logFC" %in% names(de.data)){ ## compatible for seurat3
    de.data$avg_log2FC <- de.data$avg_logFC/log(2)
  }
  if (any(colnames(de.data) == "cluster")) {
    top10.up <- de.data %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>%filter(avg_log2FC > 0) %>% arrange(-avg_log2FC)
    top5.dn <- de.data %>% group_by(cluster) %>% top_n(10, -avg_log2FC) %>%filter(avg_log2FC < 0) %>% arrange(-avg_log2FC)
  } else {
    top10.up <- de.data  %>% top_n(10, avg_log2FC) %>%filter(avg_log2FC > 0) %>% arrange(-avg_log2FC)
    top5.dn <- de.data  %>% top_n(10, -avg_log2FC) %>%filter(avg_log2FC < 0) %>% arrange(-avg_log2FC)
  }
  # top.up.dn <- rbind(top5.up, top5.dn)
  top.up.dn <- top10.up
  top.up.dn$gene <- make.unique(top.up.dn$gene)
  top.up.dn$type = ifelse(top.up.dn$avg_log2FC > 0, "positive", "negative")
  top.up.dn$type <- factor(top.up.dn$type, levels = c("positive", "negative"))
  g <- ggplot(data = top.up.dn,
              aes(x = gene, y = avg_log2FC, fill = type)) +
    geom_bar(stat="identity") +
    scale_x_discrete(limits=rev(top.up.dn$gene)) +
    theme(legend.position="none", axis.text=element_text(size=15)) +
    # scale_fill_manual(values = c(positive = "#E41A1C", negative = "#377EB8")) +
    coord_flip()
  if (!is.null(main)) {
    g <- g + ggtitle(main)
  } else {
    g <- g + ggtitle("Average logFC for the top 5 up and top 5 down regulated genes")
  }
  if (!is.null(xlim)) {
    # Coordinates are flipped
    g <- g + ylim(xlim)
  }
  return(g)
}

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
  Idents(scrna) <- "seurat_clusters"
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
                              id = as.character(makeid(scrna, id)@active.ident))
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
    # tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    x.pre <- tp/(tp + fp)
    x.rec <- tp/(tp + fn)
    de <- data.frame(x.pre, x.rec)
    gene.prauc <- rbind(gene.prauc, de)
  }

  gene.prauc <- gene.prauc[complete.cases(gene.prauc), ]
  gene.prauc <- gene.prauc[order(gene.prauc$x.rec), ]
  # if (sum(gene.prauc$x.pre == min(gene.prauc$x.pre)) > 1) {
  #   small.second.pr <- min( gene.prauc$x.pre[gene.prauc$x.pre!=min(gene.prauc$x.pre)])
  #   small.second.recall <- gene.prauc$x.rec[gene.prauc$x.pre == small.second.pr]
  #   x.add <- data.frame(min(gene.prauc$x.pre), small.second.recall)
  #   names(x.add) <- names(gene.prauc)
  #   gene.prauc <- rbind(gene.prauc, x.add)
  # }

  PRAUC <- Area_Under_Curve(gene.prauc$x.rec,  gene.prauc$x.pre,
                            method = "trapezoid", na.rm = TRUE)
  return(PRAUC)
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
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  x.max <- max(data.id$exp)
  x.min <- min(data.id$exp)

  for (x.val in seq(x.min, x.max, (x.max - x.min)/100)) {
  # for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {

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
  # if (sum(gene.prauc$x.pre == min(gene.prauc$x.pre)) > 1) {
  #   small.second.pr <- min( gene.prauc$x.pre[gene.prauc$x.pre!=min(gene.prauc$x.pre)])
  #   small.second.recall <- gene.prauc$x.rec[gene.prauc$x.pre == small.second.pr]
  #   x.add <- data.frame(min(gene.prauc$x.pre), small.second.recall)
  #   names(x.add) <- names(gene.prauc)
  #   gene.prauc <- rbind(gene.prauc, x.add)
  # }

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
                              id = as.character(makeid(scrna, id)@active.ident))
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
  # if (sum(gene.prauc$x.pre == min(gene.prauc$x.pre)) > 1) {
  #   small.second.pr <- min( gene.prauc$x.pre[gene.prauc$x.pre!=min(gene.prauc$x.pre)])
  #   small.second.recall <- gene.prauc$x.rec[gene.prauc$x.pre == small.second.pr]
  #   x.add <- data.frame(min(gene.prauc$x.pre), small.second.recall)
  #   names(x.add) <- names(gene.prauc)
  #   gene.prauc <- rbind(gene.prauc, x.add)
  # }

  # gene.prauc <- gene.prauc[order(gene.prauc$x.rec, )]
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
  markers <- FindMarkers(object = scrna, ident.1 = cellgroup, features = geneset)

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



plot_grid_2 <- function(scrna, gene1, gene2, id, step = 0.01){
  df <- data.frame(exp1 = scrna@assays$RNA@data[gene1,],
                   exp2 = scrna@assays$RNA@data[gene2,],
                   id = as.character(makeid(scrna, id)@active.ident))

  x.split <- get_split(scrna, gene1, id, step = step)
  y.split <- get_split(scrna, gene2, id, step = step)

  ggplot(df)+
    geom_point(aes(exp1, exp2, color = id)) +
    geom_segment(aes(x = x.split, y = y.split, xend = max(exp1), yend = y.split), linetype = 2)+
    geom_segment(aes(x = x.split, y = y.split, xend = x.split, yend = max(exp2)), linetype = 2)+
    xlab(paste(gene1)) + ylab(paste(gene2))
}


plot_marker_c_umap <- function(scrna, gene1, gene2, id, step = 0.01){

  gene1.split <- get_split(scrna, gene1, id, step = step)
  gene2.split <- get_split(scrna, gene2, id, step = step)

  cells.1 <- GetCellNames(cbmc, gene1, gene1.split)
  cells.2 <- GetCellNames(cbmc, gene2, gene2.split)
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
  plot_grid(p1, p2, p3)
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


get_split.1 <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  x.factor <- length(data.other$id)/length(data.id$id)

  for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {

    tp <- sum(data.id$exp >= x.val)
    fp <- sum(data.other$exp >= x.val)
    tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tn - fn*x.factor
    }
    else{
      x.margin <- tn - fn
    }

    # x.pre <- tp/(tp + fp)
    # x.rec <- tp/(tp + fn)
    # x.margin <- sum(data.id[,1] - x.val) + sum(x.val - data.other[,1])
    de <- data.frame(x.val, x.margin)
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  # x.split <- max(gene.prauc[,1])
  return(x.split)
}




marker_stepbystep <- function(scrna, cellgroup, depth = 2, geneset){
  geneset <- intersect(rownames(scrna[["RNA"]]), geneset)
  len.id <- sum(scrna@active.ident == cellgroup)
  len.other <- sum(scrna@active.ident != cellgroup)
  df.split <- data.frame()

  for (variable in seq(depth)) {
    markers <- FindMarkers(scrna, ident.1 = cellgroup, features = geneset)
    markers <- subset(markers, avg_log2FC > 0)
    markers <- subset(markers, p_val_adj < 0.05)
    geneset <- intersect(rownames(markers), geneset)

    gene.prauc <- data.frame(
      gene <- c(),
      split.value <- c(),
      bestTP <- c(),
      pre <- c(),
      recall <- c()
    )
    for (genes in geneset) {
      de <- data.frame(genes, get_split.1(scrna, genes, cellgroup, step))
      names(de)<-c("gene", "split.value","filtered.adj")

      gene.prauc <- rbind(gene.prauc, de)
    }
    gene.prauc <- gene.prauc[order(gene.prauc$filtered.adj, decreasing = T),]

    df.split <- rbind(df.split, gene.prauc[1,])

    left.cells <- GetCellNames(scrna, gene.prauc[1,1], gene.prauc[1,2])
    scrna <- subset(scrna, cells = left.cells)
  }
  return(df.split)
}


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
  for (gene in df.split$gene) {
    plots[[gene]] <- VlnPlot(scrna, features = gene) + geom_hline(yintercept=df.split[df.split$gene == gene, ]$split.value, linetype = 1)
    left.cells <- GetCellNames(scrna, gene, df.split[df.split$gene == gene, ]$split.value, df.split[df.split$gene == gene, ]$direction)
    scrna <- subset(scrna, cells = left.cells)
  }
  # scrna@active.ident <- scrna.id
  plot_grid(plotlist=plots)
}




get_split_pos <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {
    x.factor <- length(data.other$id)/length(data.id$id)
    # tp <- sum(data.id$exp >= x.val)
    # fp <- sum(data.other$exp >= x.val)
    tn <- sum(data.other$exp < x.val)
    fn <- sum(data.id$exp < x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tn - fn*x.factor
    }
    else{
      x.margin <- tn - fn
    }
    de <- data.frame(x.val, x.margin, "+")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}

get_split_neg <- function(scrna, gene, id, step = 0.01) {
  data.mat.surf <- data.frame(exp = scrna@assays$RNA@data[gene,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.val <- c(),
                           margin <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]

  for (x.val in quantile(data.id$exp, seq(0.001, 0.999, step))) {
    x.factor <- length(data.other$id)/length(data.id$id)
    # tp <- sum(data.id$exp <= x.val)
    # fp <- sum(data.other$exp <= x.val)
    tn <- sum(data.other$exp > x.val)
    fn <- sum(data.id$exp > x.val)

    if (length(data.id$id) <= length(data.other$id)) {
      x.margin <- tn - fn*x.factor
    }
    else{
      x.margin <- tn - fn
    }
    de <- data.frame(x.val, x.margin, "-")
    gene.prauc <- rbind(gene.prauc, de)
  }
  gene.prauc <- gene.prauc[order(gene.prauc$x.margin, decreasing = T),]
  x.split <- gene.prauc[1,]
  rownames(x.split) <- paste(gene)
  return(x.split)
}


marker_stepbystep <- function(scrna, cellgroup, depth = 2, geneset, step = 0.01){
  geneset <- intersect(rownames(scrna[["RNA"]]), geneset)
  len.id <- sum(scrna@active.ident == cellgroup)
  len.other <- sum(scrna@active.ident != cellgroup)
  df.split <- data.frame()

  for (variable in seq(depth)) {
    markers <- FindMarkers(scrna, ident.1 = cellgroup, features = geneset)
    markers <- subset(markers, p_val_adj < 0.05)

    markers.pos <- subset(markers, avg_log2FC > 0)
    markers.neg <- subset(markers, avg_log2FC < 0)

    #geneset <- intersect(rownames(markers), geneset)
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


get_PRAUC_matrix_combine_markers <- function(scrna, gene1,  gene2, id, step = 0.01){
  data.mat.surf <- data.frame(gene1 = scrna@assays$RNA@data[gene1,],
                              gene2 = scrna@assays$RNA@data[gene2,],
                              id = as.character(makeid(scrna, id)@active.ident))
  gene.prauc <- data.frame(x.pre <- c(),
                           x.rec <- c())
  data.id <- data.mat.surf[data.mat.surf$id == id,]
  data.other <- data.mat.surf[data.mat.surf$id != id,]
  # gene1.min <- min(data.id$gene1)
  # gene2.min <- min(data.id$gene2)
  # gene1.max <- max(data.id$gene1)
  # gene2.max <- max(data.id$gene2)
  #
  # gene1.range <-gene1.max - gene1.min
  # gene2.range <-gene2.max - gene2.min

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
  # if (sum(gene.prauc$x.pre == min(gene.prauc$x.pre)) > 1) {
  #   small.second.pr <- min( gene.prauc$x.pre[gene.prauc$x.pre!=min(gene.prauc$x.pre)])
  #   small.second.recall <- gene.prauc$x.rec[gene.prauc$x.pre == small.second.pr]
  #   x.add <- data.frame(min(gene.prauc$x.pre), small.second.recall)
  #   names(x.add) <- names(gene.prauc)
  #   gene.prauc <- rbind(gene.prauc, x.add)
  # }
  # PRAUC <- Area_Under_Curve(gene.prauc$x.rec,  gene.prauc$x.pre,
  #                           method = "trapezoid", na.rm = TRUE)
  return(gene.prauc)
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


Detect_combine_markers <- function(scrna, cellgroup, geneset, step = 0.01){
  geneset <- intersect(rownames(scrna[["RNA"]]), geneset)
  len.id <- sum(scrna@active.ident == cellgroup)
  len.other <- sum(scrna@active.ident != cellgroup)
  df.split <- data.frame()

  for (variable in seq(2)) {
    markers <- FindMarkers(scrna, ident.1 = cellgroup, features = geneset)
    markers <- subset(markers, p_val_adj < 0.05)

    markers.pos <- subset(markers, avg_log2FC > 0)
    markers.neg <- subset(markers, avg_log2FC < 0)

    #geneset <- intersect(rownames(markers), geneset)
    gene.prauc <- data.frame()

    for (genes in rownames(markers.pos)) {
      de <- data.frame(genes, get_split_pos_tp(scrna, genes, cellgroup, step))
      names(de)<-c("gene", "split.value","filtered.adj", "direction")
      gene.prauc <- rbind(gene.prauc, de)
    }
    for (genes in rownames(markers.neg)) {
      de <- data.frame(genes, get_split_neg_tp(scrna, genes, cellgroup, step))
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
