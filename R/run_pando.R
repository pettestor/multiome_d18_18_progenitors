#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param brain
#' @return
#' @author Petter Storm
#' @export
run_pando <- function(brain) {
  

  # Get motif data
  data(motifs)
  
  brain.f <- brain[,brain$nFeature_RNA > 4000 & brain$nFeature_ATAC > 9000]
  
  brain.f <- NormalizeData(brain.f)
  
  # Select variable features
  set.seed(22)
  seurat_object <- Seurat::FindVariableFeatures(brain.f[], assay='RNA')
  
  # Initiate GRN object and select candidate regions
  seurat_object <- initiate_grn(seurat_object)
  
  # Scan candidate regions for TF binding motifs
  seurat_object <- find_motifs(
    seurat_object,
    pfm = motifs,
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  
  
  # Infer gene regulatory network
  seurat_object <- infer_grn(seurat_object,parallel = T,method = "glm")
  
  #saveRDS(seurat_object,file = "raw_data/pando.seu.rds")
  #seurat_object<-readRDS("raw_data/pando.seu.rds")
  
  # Print inferred coefficients
  coef(seurat_object)
  
  # Find gene and regulatory modules 
  seurat_object <- find_modules(seurat_object)
  
  # Print modules
  NetworkModules(seurat_object)
  
  modules <- NetworkModules(seurat_object)
  View(modules@meta)
  
  write.csv(modules@meta,"2023/progenitors/pando.modules.csv")

  
  seurat_object2 <- find_modules(seurat_object,p_thresh = 0.01)
  
  
  saveRDS(seurat_object2,file="seurat.pando.rds")
  muo_data <- get_network_graph(seurat_object2, graph_name='umap_graph')
  png<-plot_network_graph(muo_data, graph='umap_graph')
  png
  ggsave("2023/progenitors/pando.clusters.pdf")
  
  
  png$data$kmeans_clus <- kmeans(png$data[,c("x","y")],centers = 4)$cluster
  
  png<-plot_network_graph2(muo_data, graph='umap_graph',max_overlaps = 15,clusters = png$data$kmeans_clus )
  png
  ggsave("2023/progenitors/pando2.clusters.pdf")
  
  
  ggplot(png$data, aes(x=x,y=y,col=factor(kmeans)))+geom_point()
  
  
  library(enrichR)
  
  
  results <- data.frame()
  
  for(center in unique(png$data$kmeans)){
    
    genes <- png$data[png$data$kmeans==center,"name"]
    
    dbs <- c("GO_Biological_Process_2015")
    enriched <- enrichr(genes, dbs)
    enriched$GO_Biological_Process_2015$cluster <- center
    results <- rbind(results,enriched$GO_Biological_Process_2015)      
    pp<-plotEnrich(enriched$GO_Biological_Process_2015)+theme_cowplot()+scale_fill_viridis()+ggtitle("")
    
    ggsave(plot = pp,paste0("2023/progenitors/enrichment.",center,".pdf"),w=10,h=10)  
    
  }
  
  write.csv(results,"2023/progenitors/enrichment.tf_clusters.csv")
  
  significant_tfs <- unique(modules@meta[modules@meta$padj<0.001,"tf"])
  significant_targets <- unique(modules@meta[modules@meta$padj<0.001,"target"])
  
  DefaultAssay(seurat_object)<-"RNA"
  seurat_object <- FindVariableFeatures(seurat_object,nfeatures=3000)
  vf <- VariableFeatures(seurat_object)
  
  DefaultAssay(seurat_object)<-"peaks"
  seurat_object <- RegionStats(seurat_object, genome = BSgenome.Hsapiens.UCSC.hg38)
  seurat_object <- LinkPeaks2(seurat_object,peak.assay = "peaks",expression.assay = "RNA",genes.use = vf)
  
  saveRDS(seurat_object,file = "raw_data/pando.seu.rds")
  
  
  muo_data <- get_network_graph(
    muo_data, 
    graph_name = 'sub_graph', 
    umap_method = 'coef',
    features = unique(significant_tfs$tf,significant_targets$target)
  )
  
  plot_network_graph(muo_data, graph='sub_graph',  color_nodes=T, node_size=5)
  
  
  muo_data <- get_tf_network(muo_data, tf='EBF3', graph='umap_graph')
  plot_tf_network(muo_data, tf='EBF3')
  
  
  
  module_pos <- modules@meta %>% 
    filter(estimate>0) %>% 
    group_by(tf) %>% filter(n()>5) %>% 
    group_split() %>% {names(.) <- map_chr(., function(x) x$tf[[1]]); .} %>% 
    map(function(x) x$target)
  
  plot_network_graph(muo_data, graph='sub_graph',  color_nodes=T, node_size=5)
  
  
  object <- LinkPeaks2(seurat_object2,peak.assay = "peaks",expression.assay = "RNA")
  
  links <- Links(object)
  
  linksdf <- as.data.frame(sort(table(links$gene)))
  linksdf$ranked <- 1:nrow(linksdf)
  
  linksdf$labelNames <- as.character(linksdf$Var1)
  linksdf[linksdf$Freq<10,"labelNames"] <- ""    
  linksdf$ranked<-linksdf$ranked*10
  linksdf[grep("RPS",linksdf$labelNames),"labelNames"] <- ""
  
  
  ggplot(linksdf,aes(x=ranked,y=Freq,label=labelNames))+geom_point(col="darkblue") +geom_vline(xintercept = 900,linetype="dotted",col="darkgrey")+
    geom_hline(yintercept = 10,linetype="dotted",col="darkgrey")+ggrepel::geom_text_repel()+ylab("Number of correlated peaks")+xlab("Ranked gene list")
  ggsave("2023/progenitors/dorcs.pdf",h=9)
}



plot_network_graph2 <- function(
    object,
    network = DefaultNetwork(object),
    graph = 'module_graph',
    layout = 'umap',
    edge_width = 0.2,
    edge_color = c('-1'='darkgrey', '1'='orange'),
    node_color = pals::magma(100),
    node_size = c(1,5),
    text_size = 10,
    color_nodes = TRUE,
    label_nodes = TRUE,
    color_edges = TRUE,
    max_overlaps=20,
    clusters = kmeans_clus
){
  library(tidygraph)
  library(ggraph)
  gene_graph <- NetworkGraph(object, network=network, graph=graph)
  
  has_umap <- 'UMAP_1' %in% colnames(as_tibble(activate(gene_graph, 'nodes')))
  if (layout=='umap' & !has_umap){
    stop('No UMAP coordinates found, please run `get_network_graph()` first.')
  }
  
  if (layout=='umap'){
    p <- ggraph(gene_graph, x=UMAP_1, y=UMAP_2)
  } else {
    p <- ggraph(gene_graph, layout=layout)
  }
  
  if (color_edges){
    p <- p + geom_edge_diagonal(aes(color=factor(dir)), width=edge_width) +
      scale_edge_color_manual(values=edge_color)
  } else {
    p <- p + geom_edge_diagonal(width=edge_width, color=edge_color[1])
  }
  
  if (color_nodes){
    p <- p + geom_node_point(aes(fill=centrality, size=centrality), color='darkgrey', shape=21) +
      scale_fill_gradientn(colors=node_color)
  } else {
    p <- p + geom_node_point(
      color='darkgrey', shape=21, fill='lightgrey', size=node_size[1], stroke=0.5
    )
  }
  
  p$data$kmeans <- clusters
  p <- p + geom_node_point(aes(fill=kmeans, size=centrality), color='darkgrey', shape=21) +
    scale_fill_gradientn(colors=node_color)
  
  
  if (label_nodes){
    p <- p + geom_node_text(
      aes(label=name),
      repel=T, size=text_size/ggplot2::.pt, max.overlaps=max_overlaps
    )
  }
  p <- p + scale_size_continuous(range=node_size) +
    theme_void() + no_legend()
  
  return(p)
}



LinkPeaks2<-function (object, peak.assay, expression.assay, peak.slot = "counts", 
                      expression.slot = "data", method = "pearson", gene.coords = NULL, 
                      distance = 5e+05, min.distance = NULL, min.cells = 10, genes.use = NULL, 
                      n_sample = 200, pvalue_cutoff = 0.05, score_cutoff = 0.05, 
                      gene.id = FALSE, verbose = TRUE) 
{
  if (!requireNamespace(package = "qlcMatrix", quietly = TRUE)) {
    stop("Please install qlcMatrix: install.packages('qlcMatrix')")
  }
  if (!inherits(x = object[[peak.assay]], what = "ChromatinAssay")) {
    stop("The requested assay is not a ChromatinAssay")
  }
  if (!is.null(x = min.distance)) {
    if (!is.numeric(x = min.distance)) {
      stop("min.distance should be a numeric value")
    }
    if (min.distance < 0) {
      warning("Requested a negative min.distance value, setting min.distance to zero")
      min.distance <- NULL
    }
    else if (min.distance == 0) {
      min.distance <- NULL
    }
  }
  features.match <- c("GC.percent", "count", "sequence.length")
  if (method == "pearson") {
    cor_method <- qlcMatrix::corSparse
  }else if (method == "spearman") {
    cor_method <- SparseSpearmanCor
  }else {
    stop("method can be one of 'pearson' or 'spearman'.")
  }
  
  if (is.null(x = gene.coords)) {
    annot <- Annotation(object = object[[peak.assay]])
    if (is.null(x = annot)) {
      stop("Gene annotations not found")
    }
    gene.coords <- CollapseToLongestTranscript(ranges = annot)
  }
  meta.features <- GetAssayData(object = object, assay = peak.assay, 
                                slot = "meta.features")
  if (!(all(c("GC.percent", "sequence.length") %in% colnames(x = meta.features)))) {
    stop("DNA sequence information for each peak has not been computed.\n", 
         "Run RegionsStats before calling this function.")
  }
  if (!("count" %in% colnames(x = meta.features))) {
    data.use <- GetAssayData(object = object[[peak.assay]], 
                             slot = "counts")
    hvf.info <- FindTopFeatures(object = data.use, verbose = FALSE)
    hvf.info <- hvf.info[rownames(meta.features), , drop = FALSE]
    meta.features <- cbind(meta.features, hvf.info)
  }
  peak.data <- GetAssayData(object = object, assay = peak.assay, 
                            slot = peak.slot)
  expression.data <- GetAssayData(object = object, assay = expression.assay, 
                                  slot = expression.slot)
  peakcounts <- rowSums(x = peak.data > 0)
  genecounts <- rowSums(x = expression.data > 0)
  peaks.keep <- peakcounts > min.cells
  genes.keep <- genecounts > min.cells
  peak.data <- peak.data[peaks.keep, ]
  
  if (!is.null(x = genes.use)) {
    genes.keep <- intersect(x = names(x = genes.keep[genes.keep]), 
                            y = genes.use)
  }
  
  expression.data <- expression.data[genes.keep, , drop = FALSE]
  
  if (verbose) {
    message("Testing ", nrow(x = expression.data), " genes and ", 
            sum(peaks.keep), " peaks")
  }
  genes <- rownames(x = expression.data)
  if (gene.id) {
    gene.coords.use <- gene.coords[gene.coords$gene_name %in% 
                                     genes, ]
    gene.coords.use$gene_name <- gene.coords.use$gene_name
  } else {
    gene.coords.use <- gene.coords[gene.coords$gene_name %in% 
                                     genes, ]
  }
  if (length(x = gene.coords.use) == 0) {
    stop("Could not find gene coordinates for requested genes")
  }
  if (length(x = gene.coords.use) < nrow(x = expression.data)) {
    message("Found gene coordinates for ", length(x = gene.coords.use), 
            " genes")
  }
  peaks <- granges(x = object[[peak.assay]])
  peaks <- peaks[peaks.keep]
  peak_distance_matrix <- DistanceToTSS(peaks = peaks, genes = gene.coords.use, 
                                        distance = distance)
  
  if (!is.null(x = min.distance)) {
    peak_distance_matrix_min <- DistanceToTSS(peaks = peaks, 
                                              genes = gene.coords.use, distance = min.distance)
    peak_distance_matrix <- peak_distance_matrix - peak_distance_matrix_min
  }
  if (sum(peak_distance_matrix) == 0) {
    stop("No peaks fall within distance threshold\n", "Have you set the proper genome and seqlevelsStyle for ", 
         peak.assay, " assay?")
  }
  genes.use <- colnames(x = peak_distance_matrix)
  all.peaks <- rownames(x = peak.data)
  peak.data <- t(x = peak.data)
  coef.vec <- c()
  gene.vec <- c()
  zscore.vec <- c()
  mylapply <- lapply
  res <- mylapply(X = seq_along(along.with = genes.use), FUN = function(i) {
    peak.use <- as.logical(x = peak_distance_matrix[, genes.use[[i]]])
    gene.expression <- t(x = expression.data[genes.use[[i]], 
                                             , drop = FALSE])
    gene.chrom <- as.character(x = seqnames(x = gene.coords.use[i]))
    if (sum(peak.use) < 2) {
      return(list(gene = NULL, coef = NULL, zscore = NULL))
    }
    else {
      peak.access <- peak.data[, peak.use, drop = FALSE]
      coef.result <- cor_method(X = peak.access, Y = gene.expression)
      rownames(x = coef.result) <- colnames(x = peak.access)
      coef.result <- coef.result[abs(x = coef.result) > 
                                   score_cutoff, , drop = FALSE]
      if (nrow(x = coef.result) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      }
      else {
        peaks.test <- rownames(x = coef.result)
        trans.peaks <- all.peaks[!grepl(pattern = paste0("^", 
                                                         gene.chrom), x = all.peaks)]
        meta.use <- meta.features[trans.peaks, ]
        pk.use <- meta.features[peaks.test, ]
        bg.peaks <- lapply(X = seq_len(length.out = nrow(x = pk.use)), 
                           FUN = function(x) {
                             MatchRegionStats(meta.feature = meta.use, 
                                              query.feature = pk.use[x, , drop = FALSE], 
                                              features.match = features.match, n = n_sample, 
                                              verbose = FALSE)
                           })
        bg.access <- peak.data[, unlist(x = bg.peaks), 
                               drop = FALSE]
        bg.coef <- cor_method(X = bg.access, Y = gene.expression)
        rownames(bg.coef) <- colnames(bg.access)
        zscores <- vector(mode = "numeric", length = length(x = peaks.test))
        for (j in seq_along(along.with = peaks.test)) {
          coef.use <- bg.coef[(((j - 1) * n_sample) + 
                                 1):(j * n_sample), ]
          z <- (coef.result[j] - mean(x = coef.use))/sd(x = coef.use)
          zscores[[j]] <- z
        }
        names(x = coef.result) <- peaks.test
        names(x = zscores) <- peaks.test
        zscore.vec <- c(zscore.vec, zscores)
        gene.vec <- c(gene.vec, rep(i, length(x = coef.result)))
        coef.vec <- c(coef.vec, coef.result)
      }
      gc(verbose = FALSE)
      pval.vec <- pnorm(q = -abs(x = zscore.vec))
      links.keep <- pval.vec < pvalue_cutoff
      if (sum(x = links.keep) == 0) {
        return(list(gene = NULL, coef = NULL, zscore = NULL))
      }
      else {
        gene.vec <- gene.vec[links.keep]
        coef.vec <- coef.vec[links.keep]
        zscore.vec <- zscore.vec[links.keep]
        return(list(gene = gene.vec, coef = coef.vec, 
                    zscore = zscore.vec))
      }
    }
  })
  
  gene.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              1))
  coef.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                              2))
  zscore.vec <- do.call(what = c, args = lapply(X = res, FUN = `[[`, 
                                                3))
  if (length(x = coef.vec) == 0) {
    if (verbose) {
      message("No significant links found")
    }
    return(object)
  }
  peak.key <- seq_along(along.with = unique(x = names(x = coef.vec)))
  names(x = peak.key) <- unique(x = names(x = coef.vec))
  coef.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = coef.vec)], 
                              x = coef.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = coef.matrix) <- genes.use
  colnames(x = coef.matrix) <- names(x = peak.key)
  links <- LinksToGRanges2(linkmat = coef.matrix, gene.coords = gene.coords.use)
  z.matrix <- sparseMatrix(i = gene.vec, j = peak.key[names(x = zscore.vec)], 
                           x = zscore.vec, dims = c(length(x = genes.use), max(peak.key)))
  rownames(x = z.matrix) <- genes.use
  colnames(x = z.matrix) <- names(x = peak.key)
  z.lnk <- LinksToGRanges2(linkmat = z.matrix, gene.coords = gene.coords.use)
  links$zscore <- z.lnk$score
  links$pvalue <- pnorm(q = -abs(x = links$zscore))
  links <- links[links$pvalue < pvalue_cutoff]
  Links(object = object[[peak.assay]]) <- links
  return(object)
}

LinksToGRanges2 <- function(linkmat, gene.coords, sep = c("-", "-")) {
  # get TSS for each gene
  tss <- resize(gene.coords, width = 1, fix = 'start')
  gene.idx <- sapply(
    X = rownames(x = linkmat),
    FUN = function(x) {
      which(x = x == tss$gene_name)[[1]]
    }
  )
  tss <- tss[gene.idx]
  
  # get midpoint of each peak
  peak.ranges <- StringToGRanges(
    regions = colnames(x = linkmat),
    sep = sep
  )
  midpoints <- start(x = peak.ranges) + (width(x = peak.ranges) / 2)
  
  # convert to triplet form
  dgtm <- as(object = linkmat, Class = "TsparseMatrix")
  
  # create dataframe
  df <- data.frame(
    chromosome = as.character(x = seqnames(x = peak.ranges)[dgtm@j + 1]),
    tss = start(x = tss)[dgtm@i + 1],
    pk = midpoints[dgtm@j + 1],
    score = dgtm@x,
    gene = rownames(x = linkmat)[dgtm@i + 1],
    peak = colnames(x = linkmat)[dgtm@j + 1]
  )
  
  # work out start and end coords
  df$start <- ifelse(test = df$tss < df$pk, yes = df$tss, no = df$pk)
  df$end <- ifelse(test = df$tss < df$pk, yes = df$pk, no = df$tss)
  df$tss <- NULL
  df$pk <- NULL
  
  # convert to granges
  gr.use <- makeGRangesFromDataFrame(df = df, keep.extra.columns = TRUE)
  return(sort(x = gr.use))
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @import data.table
CollapseToLongestTranscript <- function(ranges) {
  range.df <- as.data.table(x = ranges)
  range.df$strand <- as.character(x = range.df$strand)
  range.df$strand <- ifelse(
    test = range.df$strand == "*",
    yes = "+",
    no = range.df$strand
  )
  collapsed <- range.df[
    , .(unique(seqnames),
        min(start),
        max(end),
        strand[[1]],
        gene_biotype[[1]],
        gene_name[[1]]),
    "gene_id"
  ]
  colnames(x = collapsed) <- c(
    "gene_id", "seqnames", "start", "end", "strand", "gene_biotype", "gene_name"
  )
  collapsed$gene_name <- make.unique(names = collapsed$gene_name)
  gene.ranges <- makeGRangesFromDataFrame(
    df = collapsed,
    keep.extra.columns = TRUE
  )
  return(gene.ranges)
}

# Find peaks near genes
#
# Find peaks that are within a given distance threshold to each gene
#
# @param peaks A GRanges object containing peak coordinates
# @param genes A GRanges object containing gene coordinates
# @param distance Distance threshold. Peaks within this distance from the gene
# will be recorded.
# @param sep Separator for peak names when creating results matrix
#
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom Matrix sparseMatrix
#' @importFrom GenomicRanges resize
#
# @return Returns a sparse matrix
DistanceToTSS <- function(
    peaks,
    genes,
    distance = 200000,
    sep = c("-", "-")
) {
  tss <- resize(x = genes, width = 1, fix = 'start')
  genes.extended <- suppressWarnings(
    expr = Extend(
      x = tss, upstream = distance, downstream = distance
    )
  )
  overlaps <- findOverlaps(
    query = peaks,
    subject = genes.extended,
    type = 'any',
    select = 'all'
  )
  hit_matrix <- sparseMatrix(
    i = queryHits(x = overlaps),
    j = subjectHits(x = overlaps),
    x = 1,
    dims = c(length(x = peaks), length(x = genes.extended))
  )
  rownames(x = hit_matrix) <- GRangesToString(grange = peaks, sep = sep)
  colnames(x = hit_matrix) <- genes.extended$gene_name
  return(hit_matrix)
}


