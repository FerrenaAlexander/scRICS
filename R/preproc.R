## preproc_for_harmony
preproc_for_harmony <- function(
    sobjcat, 
    assay = 'RNA_NoCOVID',  #set to RNA for general use
    samplecolumn = 'Code',
    verbose=FALSE
){
  
  require(Seurat)
  
  # set.seed(54321) #seed will be set in cluster stabiltiy function loop at specific points
  
  #instead of passing dims
  # dims <- 1:maxdim
  DefaultAssay(sobjcat) <- assay
  
  #split
  # in case already split, wrap in try block
  try(sobjcat[[assay]] <- split(sobjcat[[assay]], f = sobjcat@meta.data[,samplecolumn]))
  
  
  #proc
  sobjcat <- NormalizeData(sobjcat, verbose=verbose)
  sobjcat <- FindVariableFeatures(sobjcat, verbose=verbose)
  
  # ## exclude covid genes from HVGs
  # hvgs <- VariableFeatures(sobjcat, verbose=verbose)
  # hvgs <- hvgs[!(hvgs %in% covidgenes)]
  # VariableFeatures(sobjcat) <- hvgs
  
  sobjcat <- ScaleData(sobjcat, verbose=verbose)
  sobjcat <- RunPCA(sobjcat, verbose=verbose)
  
  ## update; 
  sobjint <- IntegrateLayers(
    object = sobjcat, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = verbose
  )
  
  
  sobjint
  
}



integrate_cluster_harmony <- function(
    sobjint, 
    maxdim=30, 
    res=0.1, 
    run.umap = F,
    assay = 'RNA_NoCOVID',  #set to RNA for general use
    # samplecolumn = 'Code',
    verbose=FALSE
){
  
  
  require(Seurat)
  
  # set.seed(54321)
  
  #instead of passing dims
  dims <- 1:maxdim
  
  ## PREPROC SEPARATELY
  # 
  # #split
  # # in case already split, wrap in try block
  # try(sobjcat[[assay]] <- split(sobjcat[[assay]], f = sobjcat@meta.data[,samplecolumn]))
  # 
  # 
  # #proc
  # sobjcat <- NormalizeData(sobjcat, verbose=verbose)
  # sobjcat <- FindVariableFeatures(sobjcat, verbose=verbose)
  # 
  # ## exclude covid genes from HVGs
  # hvgs <- VariableFeatures(sobjcat, verbose=verbose)
  # hvgs <- hvgs[!(hvgs %in% covidgenes)]
  # VariableFeatures(sobjcat) <- hvgs
  # 
  # sobjcat <- ScaleData(sobjcat, verbose=verbose)
  # sobjcat <- RunPCA(sobjcat, verbose=verbose)
  # 
  # 
  # sobjint <- IntegrateLayers(
  #   object = sobjcat, method = HarmonyIntegration,
  #   orig.reduction = "pca", new.reduction = "harmony",
  #   verbose = verbose
  # )
  #### UPDATE - WE CAN ACTUALLY DO THIS IN THE PREPROC BRANCH - IT WILL BE A BIG SPEEDUP
  
  
  ## post int clustering and umap
  sobjint <- FindNeighbors(sobjint, reduction = "harmony", dims =  dims , verbose=verbose)
  
  if(run.umap){sobjint <- RunUMAP(sobjint, reduction = "harmony", dims = dims , reduction.name = "umap.harmony", verbose=verbose)}
  
  clusteringname <- paste0("harmony_clusters", "_PCs_1-", maxdim, '.res_', res)
  sobjint <- FindClusters(sobjint, resolution = res, cluster.name =  clusteringname, verbose=verbose)
  
  return(sobjint)
  
  
}
