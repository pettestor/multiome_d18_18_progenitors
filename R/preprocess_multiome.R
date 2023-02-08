#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Petter Storm
#' @export
preprocess_multiome <- function() {

  pfm <- getMatrixSet(
    x = JASPAR2020,
    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
  )
  
  
  if(file.exists("raw_data/multiome.d16_d18.rds")){
    cat("Reading from file raw_data/multiome.d16_d18.rds")
    brain <-readRDS("raw_data/multiome.d16_d18.rds")
  }else{
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotation) <- "UCSC"
    genome(annotation) <- "hg38"
    
    counts <- Read10X_h5("/data/fate_state/multiome/aggr/aggr/D16_18_aggr/outs/filtered_feature_bc_matrix.h5")
    fragpath <- "/data/fate_state/multiome/aggr/aggr/D16_18_aggr/outs/atac_fragments.tsv.gz"
    
    # create a Seurat object containing the RNA adata
    brain <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA"
    )
    
    # create ATAC assay and add it to the object
    brain[["ATAC"]] <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = fragpath,
      annotation = annotation
    )
    
    
    #clean up
    rm(counts)
    gc()
    
    brain$orig.ident <- factor(c("D16_2","D16_1","D18_2","D18_1")[factor(stringr::str_split_fixed(rownames(brain@meta.data),"-",2)[,2])],levels=c("D16_1","D16_2","D18_1","D18_2"))
    
    
    
    DefaultAssay(brain) <- "ATAC"
    
    brain <- NucleosomeSignal(brain)
    brain <- TSSEnrichment(brain)
    
    Idents(brain)<-brain$orig.ident
    
    VlnPlot(
      object = brain,
      features = c("nFeature_RNA", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal"),
      ncol = 2,log = T,
      pt.size = 0,cols = tableau_color_pal(palette = "Winter")(10)
    )
    ggsave("qc_plots/violin.atac_qc.pdf")
    
    brain <- subset(
      x = brain,
      subset = nCount_ATAC < 100000 &
        nCount_RNA < 25000 &
        nCount_ATAC > 1000 &
        nCount_RNA > 1000 &
        nucleosome_signal < 2 &
        TSS.enrichment > 1
    )
    
    brain
    
    # call peaks using MACS2
    peaks <- CallPeaks(brain, assay = "ATAC",macs2.path = "/home/petter/anaconda3/bin/macs2")
    
    # remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)
    
    # quantify counts in each peak
    macs2_counts <- FeatureMatrix(
      fragments = Fragments(brain),
      features = peaks,
      cells = colnames(brain)
    )
    
    # create a new assay using the MACS2 peak set and add it to the Seurat object
    brain[["peaks"]] <- CreateChromatinAssay(
      counts = macs2_counts,
      fragments = fragpath,
      annotation = annotation
    )
    
    saveRDS(brain,file="raw_data/multiome.d16_d18.rds")
    
  }
  
  
  # DNA accessibility data processing
  DefaultAssay(brain) <- "peaks"
  brain <- FindTopFeatures(brain, min.cutoff = '5')
  brain <- RunTFIDF(brain)
  brain <- RunSVD(object = brain)
  
  
  DepthCor(brain)
  ggsave("qc_plots/depthcor.atac_qc.pdf")
  
  
  
  ######################################################
  # https://stuartlab.org/signac/0.2/articles/integration.html#integration-with-harmony-1
  ######################################################
  
  brain <- RunHarmony(
    object = brain,
    group.by.vars = 'orig.ident',
    reduction = 'lsi',
    assay.use = 'peaks',
    project.dim = FALSE
  )
  
  
  ######################################################
  #  Non-linear dimension reduction and clustering
  # https://stuartlab.org/signac/0.2/articles/integration.html#integration-with-harmony-1
  ######################################################
  
  brain <- RunUMAP(brain, dims = 2:30, reduction = 'harmony')
  
  
  brain <- FindNeighbors(
    object = brain,
    reduction = 'harmony',
    dims = 2:30
  )
  
  brain <- FindClusters(
    object = brain,
    
    # algorithm = 3,
    resolution = 0.2,
    verbose = T
  )
  
  DimPlot(brain)
  
  
  
  
  ######################################################
  # Preprocess and integrate RNA
  ######################################################
  
  DefaultAssay(brain) <- "RNA"
  brain <- NormalizeData(brain)
  brain <- FindVariableFeatures(brain)
  brain <- ScaleData(brain)
  brain <- RunPCA(brain)
  brain <- RunHarmony(brain,group.by.vars = "orig.ident",assay.use = "RNA")
  brain<-RunUMAP(brain,reduction = "harmony",dims=1:30)
  
  
  return(brain)
  
}
