#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Petter Storm
#' @export
analyze_multiome <- function(progenitors.multiome) {
  
  ######################################################
  # Pre-processing workflow setup
  ######################################################
  
  
  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  progenitors.multiome$Day <- substr(progenitors.multiome$orig.ident,1,3)
  
  ######################################################
  # Annotate using Singler/ la manno
  ######################################################
  
  
  hESCs <- LaMannoprogenitors.multiomeData('human-embryo')
  hESCs <- logNormCounts(hESCs)
  
  pred.hesc <- SingleR(test = as.SingleCellExperiment(progenitors.multiome.rna),
                       ref = hESCs, labels = hESCs$Cell_type)
  
  fractiion<-as.data.frame(table(pred.hesc$labels))
  fractiion<-fractiion[fractiion$Freq<50,]
  pred.hesc$labelsPruned <- pred.hesc$labels
  pred.hesc[pred.hesc$labels %in% fractiion$Var1, "labels"]<-"Other"
  
  progenitors.multiome.rna$labels <- pred.hesc$labels
  
  
  m1<-reshape2::melt(prop.table(table(progenitors.multiome.rna$seurat_clusters,progenitors.multiome.rna$labels),1))
  
  p0 <- DimPlot(progenitors.multiome.rna)+scale_color_few()
  p1<-ggplot(m1,aes(x=Var2,y=factor(Var1),col=factor(Var1),size=value*100))+theme(legend.position = "none")+
    geom_point()+scale_size(range=c(-1,28))+
    ylab("Clusters")+xlab("La manno cell type")+theme_cowplot()+ggthemes::scale_color_tableau(palette = "Tableau 20")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p2<-DimPlot(progenitors.multiome.rna,group.by = "labels",label=T)+ggthemes::scale_color_tableau(palette = "Tableau 20")
  
  plot_grid(p0,p1,p2,ncol = 3)
  ggsave("multiome_results/annotation_lamanno.pdf",h=5.5)
  
  
  ######################################################
  # Single Cell Proportion Test
  ######################################################
  
  dp1<-DimPlot(object = progenitors.multiome, label = TRUE,split.by = "orig.ident",group.by = "seurat_clusters",ncol = 2,pt.size = .7) + NoLegend()+
    scale_color_tableau(palette = "Superfishel Stone")+ggtitle("")
  
  prop_test <- sc_utils(progenitors.multiome)
  
  prop_test <- permutation_test(
    prop_test, cluster_identity = "seurat_clusters",
    sample_1 = "D16", sample_2 = "D18",
    sample_identity = "Day"
  )
  
  pp<-permutation_plot(prop_test)+theme(legend.position = "top")
  
  plot_grid(dp1, pp,ncol = 2)
  ggsave("multiome_results/umap.pdf",h=5.5)
  
  ######################################################
  # Plot favourite genes
  ######################################################
  
  DefaultAssay(progenitors.multiome) <- 'peaks'
  progenitors.multiome <- RegionStats(progenitors.multiome, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  da.tfs <- c("FOXA2","LMX1A","EN1","OTX2", "FABP7","SOX2","WNT5A","RSPO2","MSX1","CORIN","ASCL1","NEUROG2","NEUROD2")
  
  
  progenitors.multiome <- NormalizeData(progenitors.multiome)
  
  # link peaks to genes
  progenitors.multiome <- LinkPeaks(
    object = progenitors.multiome,
    peak.assay = "peaks",
    expression.assay = "RNA",
    genes.use = da.tfs
  )

 
  # plot one by one 
  for(g in da.tfs){
    cp<-CoveragePlot(
      object = progenitors.multiome,
      region = g,
      features = g,
      expression.assay = "RNA",
      extend.upstream = 1000,
      extend.downstream = 1000,
      ncol = 1
    )
    cp & scale_fill_tableau(palette = "Superfishel Stone")
    ggsave(paste0("multiome_results/",g,".covplot3.pdf"))
    ggsave(paste0("multiome_results/",g,".covplot3.png"))
    
  }
  
  
  ######################################################
  # Differential expression Clusters
  ######################################################
  
  da_peaks <- FindAllMarkers(
    object = progenitors.multiome,
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nFeature_peaks',
    max.cells.per.ident = 500
  )
  
  head(da_peaks)

  da_peaks.cl <- cbind(da_peaks,ClosestFeature(progenitors.multiome, da_peaks$gene))
  da_peaks.cl.signi <- da_peaks.cl[da_peaks.cl$p_val_adj<0.05,]
  
  png("multiome_results/barplot.de_peaks_per_cluster.png",width = 400,h=800)
  barplot(table(da_peaks.cl.signi$cluster),horiz = T,xlab = "# de peaks",ylab="cluster")
  dev.off()
  
  m1 <- reshape2::melt(table(da_peaks.cl.signi$cluster,factor(da_peaks.cl.signi$gene_name)))
  
  genestoplot <- as.character(head(m1[order(m1$value,decreasing = T),"Var2"]))
  
  # link peaks to genes
  progenitors.multiome <- LinkPeaks(
    object = progenitors.multiome,
    peak.assay = "peaks",
    expression.assay = "RNA",
    genes.use = genestoplot
  )
  
  DefaultAssay(progenitors.multiome)<-"peaks"
 
  
  for(g in genestoplot){
    cp<-CoveragePlot(
      object = progenitors.multiome,
      region = g,
      features = g,
      expression.assay = "RNA",
      extend.upstream = 1000,
      extend.downstream = 1000,
      ncol = 1
    )
    cp & scale_fill_tableau(palette = "Superfishel Stone")
    ggsave(paste0("multiome_results/de_cluster/",g,".covplot.pdf"))
    ggsave(paste0("multiome_results/de_cluster/",g,".covplot.png"))
    
  }
  
  ######################################################
  # Finding overrepresented motifs
  ######################################################
  
  progenitors.multiome <- AddMotifs(
    object = progenitors.multiome,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
  
  enriched.motifs <- FindMotifs(
    object = progenitors.multiome,
    features = da_peaks[da_peaks$cluster=="5","gene"]
  )
  
  MotifPlot(
    object = progenitors.multiome,
    motifs = head(rownames(enriched.motifs))
  )
  ggsave("multiome_results/enriched_motifs.cluster5.pdf")
  
  
  enriched.motifs <- FindMotifs(
    object = progenitors.multiome,
    features = da_peaks[da_peaks$cluster=="3","gene"]
  )
  
  MotifPlot(
    object = progenitors.multiome,
    motifs = head(rownames(enriched.motifs))
  )
  ggsave("multiome_results/enriched_motifs.cluster5.pdf")
  
  
  # HEATMAP
  avg <- AverageExpression(progenitors.multiome,slot = "data",return.seurat = T)
  DoHeatmap(avg,features = da_peaks.cl.signi[,"gene"],group.colors = tableau_color_pal(palette = "Superfishel Stone")(10))+scale_fill_viridis()
  ggsave("multiome_results/heatmap.cluster.pdf")
  
  
  #####################################
  # Differential expression Days
  #####################################
  
  DefaultAssay(progenitors.multiome) <- 'peaks'
  Idents(progenitors.multiome)<-progenitors.multiome$Day
  
  da_peaks <- FindAllMarkers(
    object = progenitors.multiome,
    min.pct = 0.05,
    test.use = 'LR',
    latent.vars = 'nFeature_peaks',
    max.cells.per.ident = 5000,
    logfc.threshold = 0.15,
    only.pos = T
  )
  
  head(da_peaks)
  
  da_peaks.cl <- cbind(da_peaks,ClosestFeature(progenitors.multiome, da_peaks$gene))
  
  
  da_peaks.cl.signi <- da_peaks.cl[da_peaks.cl$p_val_adj<0.05,]
  
  genestoplot <- unique(da_peaks.cl.signi$gene_name)
  
  # link peaks to genes
  progenitors.multiome <- LinkPeaks(
    object = progenitors.multiome,
    peak.assay = "peaks",
    expression.assay = "RNA",
    genes.use = genestoplot
  )
  
  DefaultAssay(progenitors.multiome)<-"peaks"
  
  
  for(g in genestoplot){
    cp<-CoveragePlot(
      object = progenitors.multiome,
      region = g,
      features = g,
      expression.assay = "RNA",
      extend.upstream = 1000,
      extend.downstream = 1000,
      ncol = 1
    )
    cp & scale_fill_tableau(palette = "Superfishel Stone")
    ggsave(paste0("multiome_results/de_day/",g,".covplot.pdf"))
    ggsave(paste0("multiome_results/de_day/",g,".covplot.png"))
    
  }
  
  
  ######################################################
  # Motifs D16 vs 18
  ######################################################
  
  progenitors.multiome <- AddMotifs(
    object = progenitors.multiome,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    pfm = pfm
  )
  
  
  
  enriched.motifs <- FindMotifs(
    object = progenitors.multiome,
    features = da_peaks[da_peaks$cluster=="D16","gene"]
  )
  
  MotifPlot(
    object = progenitors.multiome,
    motifs = head(rownames(enriched.motifs))
  )
  ggsave("multiome_results/de_day/enriched_motifs.D16.pdf")
  
  
  enriched.motifs <- FindMotifs(
    object = progenitors.multiome,
    features = da_peaks[da_peaks$cluster=="D18","gene"]
  )
  
  MotifPlot(
    object = progenitors.multiome,
    motifs = head(rownames(enriched.motifs))
  )
  ggsave("multiome_results/de_day/enriched_motifs.D18.pdf")
  
  
  # HEATMAP
  avg <- AverageExpression(progenitors.multiome,slot = "data",return.seurat = T)
  progenitors.multiome <- ScaleData(progenitors.multiome,features = da_peaks.cl.signi[,"gene"])
  DoHeatmap(progenitors.multiome[,sample(1:42000,50)],features = da_peaks.cl.signi[,"gene"],group.colors = tableau_color_pal(palette = "Superfishel Stone")(10))+scale_fill_viridis()
  ggsave("multiome_results/heatmap.cluster.pdf")
  
  return()
  
  
}
