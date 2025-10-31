## preproc_for_harmony
preproc_for_harmony <- function(
    sobjcat, 
    assay = 'RNA_NoCOVID',  #set to RNA for general use
    samplecolumn = 'Code',
    verbose=FALSE
){
  
  require(Seurat)
  
  # set.seed(54321)
  
  #instead of passing dims
  # dims <- 1:maxdim
  
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




cluster_stability_sweep <- function(
    sobjint, 
    numreps = 50, 
    propcells.perrep = 0.8,
    sweep_maxPCs = c(5,10,15,20,25,30,40,50),
    sweep_res = seq(0.1, 1.5, by=0.2),
    workernum = 1,
    assay = DefaultAssay(sobjint),
    samplecolumn = 'Code',
    verbose = T,
    outdir = './'
){
  
  
  
  #preproc and clustering
  require(Seurat)
  
  #data handling
  require(tidyr)
  require(dplyr)
  require(stringr)
  require(tibble)
  require(magrittr)
  
  #parallelization of sweep
  require(foreach)
  require(doParallel)
  require(parallel)
  
  #ARI
  require(mclust)
  
  
  
  #set seed
  set.seed(54321)
  
  
  #record time
  timestart = proc.time()
  
  
  
  #make dirs
  if(verbose){message('Will write results to:\n', outdir)}
  dir.create(outdir, recursive = T)
  logdir <- paste0(outdir, '/logs/'); dir.create(logdir)
  sweeplogdir <- paste0(logdir, '/SweepLogs/'); dir.create(sweeplogdir)
  rawoutsdir <- paste0(outdir, '/RawOuts/' ); dir.create(rawoutsdir, recursive = T)
  globrefdir <- paste0(outdir, '/RawOuts/GlobalReferenceClustering/'); dir.create(globrefdir, recursive = T)
  bootstrapresdir <- paste0(outdir, '/RawOuts/PerBootstrapClustering/'); dir.create(bootstrapresdir, recursive = T)
  
  
  
  
  
  
  ## for each combo of dims and nPCs, loop and combine...
  # paramfield <- paste0("PCs_1-", sweep_maxPCs, ".res_", sweep_res)
  paramfield <- unlist(lapply(sweep_maxPCs, function(PCi){paste0("PCs_1-", PCi,  ".res_", sweep_res)}))
  names(paramfield) <- paramfield
  
  
  
  ### write to global logfile...
  global_logfile <- paste0(logdir, '/GlobalLogFile.txt')
  
  write_lines_l <- list(date(),
                        
                        paste0('\n'),
                        paste0('KEY PACKAGE VERSIONS:'),
                        paste0('Seurat version - ', packageVersion('Seurat')),
                        paste0('Harmony version - ', packageVersion('harmony')),
                        
                        paste0('\n'),
                        paste0('RUN OPTIONS:'),
                        
                        paste0('numreps - ', numreps),
                        paste0('propcells.perrep - ', propcells.perrep),
                        paste0('sweep_maxPCs - ', paste(sweep_maxPCs, collapse = ', ')),
                        paste0('sweep_res - ', paste(sweep_res, collapse = ', ')),
                        
                        paste0('workernum - ', workernum),
                        paste0('assay - ', assay),
                        paste0('samplecolumn - ', samplecolumn),
                        paste0('verbose - ', verbose),
                        paste0('outdir - ', outdir)
                        
                        
                        
                        
  )
  
  write('Cluster Stability Log\n\n', global_logfile)
  
  lapply(write_lines_l, function(x){write(x, global_logfile, append = T)})
  
  
  
  
  
  
  ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  ### compute reference clusters on global ###
  ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  
  
  ## if done before, read in rather than rerun... this is all so slow
  saved_globrefs <- list.files(globrefdir)
  saved_globrefs <- gsub('\\.rds', '', saved_globrefs)
  names(saved_globrefs) <- saved_globrefs
  
  # exclude existing results from the global run
  globparamfield_run <- paramfield[!(paramfield %in% names(saved_globrefs))]
  
  
  
  ## if already all pre-run, we just skip global entirely
  if(length(globparamfield_run) > 0){
    
    
    if(verbose){
      message('\n\n\n', 'COMPUTING REFERENCE CLUSTERS ACROSS PARAMETER FIELD')
    }
    
    # try to minimize size as much as possible
    if(verbose){
      message('\nDownsizing object')
    }
    try(sobjint <- JoinLayers(sobjint, assay = assay, layers = 'counts'))
    sobjint <- CreateSeuratObject(
      GetAssayData(sobjint, assay = DefaultAssay(sobjint), layer = 'counts'),
      assay = assay,
      meta.data = sobjint@meta.data
    )
    Sys.sleep(3)
    gc(full = T)
    
    
    # preproc before reference computation
    if(verbose){
      message('\nPreproc object')
    }
    sobjint <- preproc_for_harmony(
      sobjint, 
      assay = assay,
      samplecolumn = 'Code',
      verbose = verbose
    )
    
    
    ## initiate parallelization
    
    
    
    
    
    
    
    
    
    #use just len paramfield, no need to use all cpus for this
    if(verbose){
      message('\n\n\nInitiating Global Reference Parameter Sweep')
    }
    
    refworkernum <- ifelse(workernum > length(globparamfield_run),
                           yes = length(globparamfield_run), no = workernum
    )
    cl <- parallel::makeCluster(refworkernum, outfile='')
    doParallel::registerDoParallel(cl)
    
    
    
    
    
    #### save then del from global and read in per-worker... ####
    optimize.mem = T
    if(optimize.mem){
      message('\nPrep to parallelize: save to disk...')
      path_sobjint_temp <- paste0(rawoutsdir, '/TEMP_SOBJINT.RDS')
      saveRDS(sobjint,  path_sobjint_temp )
      
      #remove from mem
      rm(sobjint); Sys.sleep(3); gc(full = T)
      
      ## in each worker, load the file:
      message('\nLoad in each worker...')
      parallel::clusterExport(cl, varlist = "path_sobjint_temp", envir = environment())
      
      clusterEvalQ(cl, {
        Sys.sleep(runif(1, 3, 5))   # jitter for a few seconds
        sobjint <- readRDS( path_sobjint_temp )
        NULL
      })
      
    }
    
    
    ## exclude 'sobjint' from export because we do the above processes..
    toexport <- c('sobjint','assay', 'sweeplogdir', 'integrate_cluster_harmony', 'preproc_for_harmony')
    if(optimize.mem){
      toexport <- c('assay', 'sweeplogdir', 'integrate_cluster_harmony', 'preproc_for_harmony')
    }
    
    
    # params_i = paramfield[1]
    # params_i = globparamfield_run[1] #for testing
    
    globref_l <- foreach(
      params_i = globparamfield_run, 
      .export = toexport,
      .packages = c('Seurat', 'stringr'),
      .verbose = verbose
    ) %dopar% {
      
      
      ### try setting this to force dependency packages not to use many threads...?
      RhpcBLASctl::blas_set_num_threads(1)
      Sys.setenv(OMP_NUM_THREADS = 1)
      
      set.seed(54321)
      
      
      ### prep per-worker logfile...
      # Per-worker logfile name (same file for all iterations on the same worker)
      pid <- Sys.getpid()
      logfile <- file.path(sweeplogdir, paste0("GlobalReference.",
                                               "ParamCombo.", params_i, '.',
                                               "worker-", pid, ".log"))
      
      # Open a connection and redirect both stdout and messages to the logfile
      con <- file(logfile, open = "a")
      sink(con)                   # redirects cat(), print(), etc.
      sink(con, type = "message") # redirects message(), warning(), condition messages
      
      # Ensure sinks are reset and connection closed even if errors occur
      on.exit({
        # restore message and output sinks in reverse order
        sink(type = "message")
        sink()
        close(con)
      }, add = TRUE)
      
      
      message(sprintf("[%s] start iter %s (pid=%s)", Sys.time(), params_i, pid))
      
      ### prep per-worker logfile (close)
      
      
      
      
      #error handling; return NULL if any error
      outdf <- NULL
      try({
        
        
        #extract PCs and res
        params <- stringr::str_split_fixed(params_i, '\\.', 2)
        params[1,] <- gsub('PCs_1-', '', params[1,])
        params[1,] <- gsub('res_', '', params[1,])
        maxdim <- as.numeric(params[1,1])
        res <- as.numeric(params[1,2])
        
        #run the preproc, int and clustering
        out <- integrate_cluster_harmony(
          sobjint, 
          maxdim, 
          res=res, 
          run.umap = F, 
          assay = assay,
          verbose = verbose
        )
        
        # clustcol <- out@meta.data[,grepl(pattern = params_i, colnames(out@meta.data))]
        outmd <- out@meta.data
        outdf <- data.frame(barcode = rownames(outmd),
                            cluster = as.vector(outmd[,grepl(pattern = params_i, colnames(outmd))])
        )
        rm(out, outmd)
        Sys.sleep(3)
        gc(full = T)
        
        
        
        
      })
      
      
      message('saving output...')
      saveRDS(outdf, paste0(globrefdir, '/', params_i, '.rds'))
      
      
      
      #print to per-worker logfile...
      message(sprintf("[%s] finished iter %s", Sys.time(), params_i))
      
      outdf
      
      
    }
    
    parallel::stopCluster(cl)
    names(globref_l) <- globparamfield_run
    
    
    
    ## read sobjint back in
    if(optimize.mem){
      message('\nDone with foreach; reading sobjint back in to global')
      
      path_sobjint_temp <- paste0(rawoutsdir, '/TEMP_SOBJINT.RDS')
      sobjint <- readRDS(path_sobjint_temp)
      
      #delete the temp
      unlink(path_sobjint_temp)
      
    }
    
    
  } else{
    globref_l <- list() ## if already all pre-run, we just skip global entirely
  }
  
  
  #bind with existing res
  if(length(saved_globrefs)>0){
    message('\nThe following global refs have been found and will be read in:\n',
            paste(saved_globrefs, collapse = '\n'))
    
    presaved_globref_l <- lapply(saved_globrefs, function(params_i){
      if( file.exists( paste0(globrefdir, '/', params_i, '.rds') ) ){
        readRDS(paste0(globrefdir, '/', params_i, '.rds'))
      } 
    })
    
    globref_l <- c(presaved_globref_l, globref_l)
    globref_l <- globref_l[paramfield]
    
  }
  
  ## exclude paramfield values which failed
  globref_l <- globref_l[lengths(globref_l) > 0]
  failedparamfield <- paramfield[!(paramfield) %in% names(globref_l)]
  
  if(length(failedparamfield)>0){
    warning('THE FOLLOWING PARAMFIELD CHOICES FAILED GLOBAL CLUSTERING:\n',
            paste(failedparamfield, collapse = ','))
    
    paramfield <- paramfield[!(paramfield %in% failedparamfield)]
  }
  
  
  ### write to global logfile
  #for logfile, report num clusters per param field...
  global_numclustsperparam <- sapply(1:length(globref_l), function(i){
    
    paste0(names(globref_l)[i], ' - ', length(unique(globref_l[[i]]$cluster)), ' clusters')
    
  })
  
  write_lines_l <- list(
    
    paste0('NUMBER OF CLUSTERS IN FULL DATASET PER PARAMETER CHOICE:'),
    global_numclustsperparam,
    
    
    paste0('\n\n'),
    
    paste0('FAILED PARAMETER FIELDS IN FULL DATASET (IF ANY):'),
    paste(failedparamfield, collapse = '\n')
    
    
    
  )
  
  write('\n\nGLOBAL REFERENCE PARAMETER SWEEP SUMMARY\n\n', global_logfile, append = T)
  
  lapply(write_lines_l, function(x){write(x, global_logfile, append = T)})
  
  
  
  
  # globref <- as.data.frame(dplyr::bind_rows(lapply(globref_l, function(globref_i){globref_i[,2]})))
  # rownames(globref) <- rownames(sobjint@meta.data)
  # rm(globref_l)
  # Sys.sleep(3)
  # gc(full = T)
  
  
  if(verbose){
    message('\n\n\nDone with Global Reference Parameter Sweep!!!')
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  ### compute bootstrap clustering across parameter field ###
  ###   ###   ###   ###   ###   ###   ###   ###   ###   ### 
  
  
  
  ## if done before, read in rather than rerun... this is all so slow
  saved_bootstrapres <- list.files(bootstrapresdir)
  saved_bootstrapres <- gsub('\\.rds', '', saved_bootstrapres)
  #extract bootstrap ID and param ID
  saved_bootstrapres <- data.frame(Filename = saved_bootstrapres,
                                   BootstrapID = stringr::str_split_fixed(saved_bootstrapres, pattern = '\\.', 2)[,1],
                                   ParamID = stringr::str_split_fixed(saved_bootstrapres, pattern = '\\.', 2)[,2]
  )
  
  #add numeric bootstrap id -> not sure if we ever end up using this
  # saved_bootstrapres$BootstrapIDNumeric <- as.numeric(gsub('Bootstrap-', '', saved_bootstrapres$BootstrapID ))
  
  
  # exclude existing results from the global run
  #define the entire search space: all BootstrapID + ParamID
  bootstrapreps_run_full <- unlist(lapply(1:length(paramfield), function(i){1:numreps}))
  bootstrapreps_run_full <- sort(bootstrapreps_run_full)
  names(bootstrapreps_run_full) <- unlist(lapply(1:numreps, function(repi){paste0('Bootstrap-', repi, '.', paramfield)}))
  
  #this could be useful, so keep it as a df
  bootstrapreps_run_full_df <- names(bootstrapreps_run_full)
  bootstrapreps_run_full_df <- data.frame(Filename = bootstrapreps_run_full_df,
                                          BootstrapID = stringr::str_split_fixed(bootstrapreps_run_full_df, pattern = '\\.', 2)[,1],
                                          ParamID = stringr::str_split_fixed(bootstrapreps_run_full_df, pattern = '\\.', 2)[,2]
  )
  
  #filter the full list, by the saved list; we just need to run on ones not pre-saved (ie, if interrupted)
  bootstrapreps_run <- bootstrapreps_run_full
  bootstrapreps_run <- bootstrapreps_run[!(names(bootstrapreps_run) %in% saved_bootstrapres$Filename)]
  
  
  
  
  ## if already all pre-run, we just skip the bootstrap entirely
  if( length(bootstrapreps_run) > 0 ){
    
    
    
    
    #get the bootstrap IDs to actually run on
    bootstrapreps_run_uniqueBSids <- unique(bootstrapreps_run)
    # we could parallelize by BS + Parameter ID, but then we'd need to re-run normalization; avoid
    
    
    
    # try to minimize size as much as possible
    if(verbose){
      message('\nDownsizing object')
    }
    try(sobjint <- JoinLayers(sobjint, assay = assay, layers = 'counts'))
    sobjint <- CreateSeuratObject(
      GetAssayData(sobjint, assay = DefaultAssay(sobjint), layer = 'counts'),
      assay = assay,
      meta.data = sobjint@meta.data
    )
    Sys.sleep(3)
    gc(full = T)
    
    
    ## initiate parallelization
    if(verbose){
      message('\n\n\nInitiating bootstrap clustering procedure across parameter field')
    }
    
    # if num cpus is > num remaining bootstraps to run, then lower num cpus
    bsworkernum <- ifelse(workernum > length(bootstrapreps_run_uniqueBSids),
                          yes = length(bootstrapreps_run_uniqueBSids), no = workernum
    )
    
    cl <- parallel::makeCluster(bsworkernum, outfile='')
    doParallel::registerDoParallel(cl)
    
    
    #get metadata
    globmd <- sobjint@meta.data
    
    
    #### save then del from global and read in per-worker... ####
    optimize.mem = T
    if(optimize.mem){
      message('\nPrep to parallelize: save to disk...')
      path_sobjint_temp <- paste0(rawoutsdir, '/TEMP_SOBJINT.RDS')
      saveRDS(sobjint,  path_sobjint_temp )
      
      #remove from mem
      rm(sobjint); Sys.sleep(3); gc(full = T)
      
      ## in each worker, load the file:
      message('\nLoad in each worker...')
      parallel::clusterExport(cl, varlist = "path_sobjint_temp", envir = environment())
      
      clusterEvalQ(cl, {
        Sys.sleep(runif(1, 3, 5))   # jitter for a few seconds
        sobjint <- readRDS( path_sobjint_temp )
        NULL
      })
      
    }
    
    ## exclude 'sobjint' from export because we do the above processes..
    toexport <- c('sobjint', 'globmd', 'assay', 'samplecolumn', 'integrate_cluster_harmony', 'preproc_for_harmony')
    if(optimize.mem){
      toexport <- c('globmd', 'assay', 'samplecolumn', 'integrate_cluster_harmony', 'preproc_for_harmony')
    }
    
    
    repi = 1
    
    # lapply(1:numreps, function(repi){
    repres <- foreach(
      repi = bootstrapreps_run_uniqueBSids, 
      .export = toexport,
      .noexport = c('globref_l'),
      .packages = c('Seurat', 'stringr'),
      .verbose = verbose
    ) %dopar% {
      
      
      ### try setting this to force dependency packages not to use many threads...?
      RhpcBLASctl::blas_set_num_threads(1)
      Sys.setenv(OMP_NUM_THREADS = 1)
      
      
      #set seed...
      set.seed(repi)
      
      #for this rep, get the associated param_fiels
      thisrep_paramIDs <- names(bootstrapreps_run[bootstrapreps_run == repi])
      thisrep_paramIDs <- str_split_fixed(thisrep_paramIDs, '\\.', 2)[,2]
      names(thisrep_paramIDs) <- thisrep_paramIDs
      
      ### prep per-worker logfile...
      # Per-worker logfile name (same file for all iterations on the same worker)
      pid <- Sys.getpid()
      logfile <- file.path(sweeplogdir, paste0("Bootstrap.", repi, '.',
                                               "worker-", pid, ".log"))
      
      # Open a connection and redirect both stdout and messages to the logfile
      con <- file(logfile, open = "a")
      sink(con)                   # redirects cat(), print(), etc.
      sink(con, type = "message") # redirects message(), warning(), condition messages
      
      # Ensure sinks are reset and connection closed even if errors occur
      on.exit({
        # restore message and output sinks in reverse order
        sink(type = "message")
        sink()
        close(con)
      }, add = TRUE)
      
      
      message(sprintf("[%s] start iter %s (pid=%s)", Sys.time(), repi, pid))
      
      ### prep per-worker logfile (close)
      
      
      
      
      
      #sample and subset
      subcells <- sample(x = nrow(globmd), 
                         size = round(nrow(globmd) * propcells.perrep), 
                         replace = F)
      submd <- globmd[subcells,]
      
      ## update 2025.09.30: confirmed that less than 2 cells breaks seurat harmony wrapper. splitting layers fails.
      # throw out...? force all cells?
      # just drop them. the goal of this resampling is to induce noise by perturbation.
      codetab <- table(submd[,samplecolumn])
      if(any(codetab) < 2){
        toolowsamps = codetab[codetab<2]
        submd <- submd[!(submd[,samplecolumn] %in% names(toolowsamps)),]
        
      }
      
      sobjsub <- sobjint[,rownames(submd)]
      
      
      
      #clean mem
      rm(subcells, submd)
      Sys.sleep(3)
      gc(full = T)
      
      
      ## run pre-proc: run normalization, scaling, PCA for this subset of cells,
      # BEFORE the loop of clustering over num PCs and res
      message('\nPreproc for bootstrap:', repi)
      sobjsub <- preproc_for_harmony(
        sobjsub, 
        assay = DefaultAssay(sobjsub),
        samplecolumn = 'Code',
        verbose = verbose
      )
      
      
      
      
      
      message('\n\nParamSweep for bootstrap: ', repi)
      
      #instead of looping over  paramfield, we loop over the non-previously-checked params only
      # params_i = paramfield[1]
      params_i <- thisrep_paramIDs[1]
      
      paramres_repi_l <- lapply(thisrep_paramIDs, function(params_i){
        
        #extract PCs and res
        params <- stringr::str_split_fixed(params_i, '\\.', 2)
        params[1,] <- gsub('PCs_1-', '', params[1,])
        params[1,] <- gsub('res_', '', params[1,])
        maxdim <- as.numeric(params[1,1])
        res <- as.numeric(params[1,2])
        
        message('\n\n\nINITIATE BOOTSTRAP PARAMSWEEP:\n',
                '- bootstrap: ', repi,'; ', 'param_i: ', params_i)
        
        outdf <- NULL
        try({
          #run the preproc, int and clustering
          out <- integrate_cluster_harmony(
            sobjsub, 
            maxdim, 
            res=res, 
            run.umap = F, 
            assay = DefaultAssay(sobjsub),
            verbose = verbose
          )
          
          # clustcol <- out@meta.data[,grepl(pattern = params_i, colnames(out@meta.data))]
          outmd <- out@meta.data
          outdf <- data.frame(barcode = rownames(outmd),
                              cluster = as.vector(outmd[,grepl(pattern = params_i, colnames(outmd))])
          )
          rm(out, outmd)
          Sys.sleep(3)
          gc(full = T)
          
        })
        
        message('saving output...')
        saveRDS(outdf, paste0(bootstrapresdir, '/Bootstrap-', repi, '.',params_i,  '.rds') )
        
        
        outdf
        
      }) #close paramfield lapply for this subsample
      
      # names(paramres_repi_l) <- thisrep_paramIDs
      
      
      #report fails...?
      paramres_repi_l <- paramres_repi_l[lengths(paramres_repi_l)>0]
      fails <- paramfield[!(paramfield) %in% names(paramres_repi_l)]
      
      if(length(fails)>0){
        message('\n\n',
                'For Bootstrap repi: ', repi, '\nThe following param combos failed:\n',
                paste(fails, collapse = ', '))
        
      }
      
      # names(paramres_repi_l) <- paramfield
      
      #clean mem
      rm(sobjsub)
      Sys.sleep(3)
      gc(full = T)
      
      
      message(sprintf("[%s] finished iter %s", Sys.time(), repi))
      
      
      #return the param sweep for this bootstrap subsample
      paramres_repi_l
      
    } #close bootstrap foreach
    
    parallel::stopCluster(cl)
    
    
    ## read sobjint back in --> ACTUALLY WE DON'T REALLY NEED ANYMORE FROM THIS POINT SO JUST UNLINK THE TMP
    if(optimize.mem){
      message('\nDone with foreach; reading sobjint back in to global')
      
      path_sobjint_temp <- paste0(rawoutsdir, '/TEMP_SOBJINT.RDS')
      # sobjint <- readRDS(path_sobjint_temp)
      
      #delete the temp
      unlink(path_sobjint_temp)
      
    }
    
    
    names(repres) <- paste0("Bootstrap-", bootstrapreps_run_uniqueBSids)
    
    ## flatten repres...
    # rm(repres)
    repres_flat <- unlist(repres, recursive = F)
    
    
    
  } else{
    repres_flat <- list() ## if already all pre-run, we just skip the bootstrap entirely
  }
  
  
  
  
  
  
  
  
  #bind with existing res
  if(nrow(saved_bootstrapres)>0){
    
    message('\nThe following bootstrap refs have been found and will be read in:\n',
            paste(saved_bootstrapres$Filename, collapse = '\n'))
    
    presaved_repres_flat_l <- lapply( saved_bootstrapres$Filename, function(bs_param_i){
      
      path_bs_param_i <- paste0(bootstrapresdir, '/', bs_param_i, '.rds')
      readRDS(path_bs_param_i)
      
    })
    
    names(presaved_repres_flat_l) <- saved_bootstrapres$Filename
    
    
    repres_flat <- c(repres_flat, presaved_repres_flat_l)
    repres_flat <- repres_flat[bootstrapreps_run_full_df$Filename]
    
  }
  
  
  
  ## exclude bootstrap_paramfield values which failed
  repres_flat <- repres_flat[lengths(repres_flat) > 0]
  failed_bsID_ParamIDs <- bootstrapreps_run_full_df$Filename[!(bootstrapreps_run_full_df$Filename) %in% names(repres_flat)]
  
  if(length(failed_bsID_ParamIDs)>0){
    warning('THE FOLLOWING PARAMFIELD CHOICES FAILED BOOTSTRAP CLUSTERING:\n',
            paste(failed_bsID_ParamIDs, collapse = ','))
  }
  #if no fails we can still set this outside
  bootstrapreps_run_full_df_NOFAILS <- bootstrapreps_run_full_df[bootstrapreps_run_full_df$Filename %in% names(repres_flat),]
  
  
  
  ### write to global logfile
  #for logfile, report num clusters per param field...
  # bootstrap_numclustsperparam <- sapply(1:length(repres_flat), function(i){
  #   paste0(names(repres_flat)[i], ' - ', length(unique(repres_flat[[i]]$cluster)), ' clusters')
  #   
  # })
  
  ## report MEAN NUMBER OF PER_PARAM CHOICES...
  # params_i = unique(bootstrapreps_run_full_df_NOFAILS$ParamID)[1]
  unique_param_id_meanclusts <- sapply(unique(bootstrapreps_run_full_df_NOFAILS$ParamID), function(params_i){
    
    params_i_filenames <- bootstrapreps_run_full_df_NOFAILS[bootstrapreps_run_full_df_NOFAILS$ParamID ==params_i, "Filename"]
    
    # params_i_filename_i = params_i_filenames[1]
    mean_numclust <- round(mean(unlist(lapply(params_i_filenames, function(params_i_filename_i){
      length(unique(repres_flat[[params_i_filename_i]]$cluster))
    }))))
    
    
    paste0(params_i, ' - ', mean_numclust, ' clusters')
    
    
  })
  
  
  write_lines_l <- list(
    
    paste0('FOR EACH PARAMETER CHOICE, AVERAGE NUMBER OF CLUSTERS ACROSS BOOTSTRAP SAMPLES:'),
    unique_param_id_meanclusts,
    
    
    paste0('\n\n'),
    
    paste0('FAILED BOOTSTRAP - PARAMETER COMBINATIONS (IF ANY):'),
    paste(failed_bsID_ParamIDs, collapse = '\n')
    
    
    
  )
  
  write('\n\nBOOTSTRAP REFERENCE PARAMETER SWEEP SUMMARY\n\n', global_logfile, append = T)
  
  lapply(write_lines_l, function(x){write(x, global_logfile, append = T)})
  
  
  
  
  
  
  
  
  
  
  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  ###   ###   ###   ###   ### 
  ### compute stability metrics usign full and bootstrapped clusters ###
  ###   ###   ###   ###   ###   ###   ###   ###   ###   ###  ###   ###   ###   ###   ### 
  
  
  
  ### using the bootstrap res, compute stability metrics
  
  ## globref_l -> the global full dataset reference clusters
  ## repres_flat -> the bootstrap paramter clusters
  ## bootstrapreps_run_full_df_NOFAILS -> a dataframe with bootstram id and paramid of the bootstrap res list
  # names(globref_l)
  # names(repres_flat)
  # bootstrapreps_run_full_df_NOFAILS
  
  ## rename the bootstrap index df
  bsindex <- bootstrapreps_run_full_df_NOFAILS
  
  
  ## first, try saving / reading these
  if(verbose){
    message('\n\nWriting clustering res to:\n',
            rawoutsdir)
  }
  
  
  
  ## actually just save a list of these three things
  # rawparamfieldout <- list(globref_l=globref_l,
  #                          repres_flat=repres_flat,
  #                          bsindex=bsindex
  # )
  # #update 2025/10/31 --> we define this abvoe and wrap whole pipeline
  # # filepath_rawparamfieldout <- paste0(rawoutsdir, '/rawparamfieldout.rds')
  # saveRDS(rawparamfieldout, filepath_rawparamfieldout)
  
  
  
  
  ## initiate stability metric computation
  
  if(verbose){message('\n\nInitiating Stability Metric Computation')}
  
  
  
  
  
  
  ### reset seed; it may no longer be important anyway though
  set.seed(54321)
  
  
  
  #extract "paramfield" from globef_l, instead of earlier paramref
  paramfield <- names(globref_l)
  names(paramfield) <- paramfield
  
  
  
  #loop over paramfield: get the full clusters, map the bootstrap clusters to them, and compute ARI / Jaccard
  
  
  
  refworkernum <- ifelse(workernum > length(paramfield),
                         yes = length(paramfield), no = workernum
  )
  cl <- parallel::makeCluster(refworkernum, outfile='')
  doParallel::registerDoParallel(cl)
  
  # params_i <- paramfield[1]
  # param_outlist_l <- lapply(paramfield, function(params_i){
  param_outlist_l <- foreach(
    params_i = paramfield, 
    .export = c("globref_l", "bsindex", "repres_flat"), 
    .noexport = c('sobjint'),
    .packages = c('stringr', 'dplyr', 'mclust'),
    .verbose = verbose
  ) %dopar% {
    
    
    
    
    ### try setting this to force dependency packages not to use many threads...?
    RhpcBLASctl::blas_set_num_threads(1)
    Sys.setenv(OMP_NUM_THREADS = 1)
    
    set.seed(54321)
    
    ### prep per-worker logfile...
    # Per-worker logfile name (same file for all iterations on the same worker)
    pid <- Sys.getpid()
    logfile <- file.path(sweeplogdir, paste0("StabilityComp.",
                                             "ParamCombo.", params_i, '.',
                                             "worker-", pid, ".log"))
    
    # Open a connection and redirect both stdout and messages to the logfile
    con <- file(logfile, open = "a")
    sink(con)                   # redirects cat(), print(), etc.
    sink(con, type = "message") # redirects message(), warning(), condition messages
    
    # Ensure sinks are reset and connection closed even if errors occur
    on.exit({
      # restore message and output sinks in reverse order
      sink(type = "message")
      sink()
      close(con)
    }, add = TRUE)
    
    
    message(sprintf("[%s] start iter %s (pid=%s)", Sys.time(), params_i, pid))
    
    ### prep per-worker logfile (close)
    
    if(verbose){message('\n',params_i)}
    
    
    #get fulldatares
    globref <- globref_l[[params_i]]
    
    ## loop over bootstrap res, get the bootstrap_parameter combo result
    bsindex_params_i = bsindex[bsindex$ParamID == params_i,]
    bsidx_i = bsindex_params_i$Filename[1]
    
    thisparam_bsres_l <- lapply(bsindex_params_i$Filename, function(bsidx_i){
      
      if(verbose){message(' - ',bsidx_i)}
      
      #get the repres for this bootstrap_param combo
      bsidx_i_df <- repres_flat[[bsidx_i]]
      
      #prepare to subset and match the barcodes
      globref_bsidx_i <- globref
      
      #match them
      intbarcodes <- intersect(globref_bsidx_i$barcode, bsidx_i_df$barcode)
      bsidx_i_df <- bsidx_i_df[match(intbarcodes, bsidx_i_df$barcode),]
      globref_bsidx_i <- globref_bsidx_i[(match(intbarcodes, globref_bsidx_i$barcode)),]
      
      #first, get ARI
      ari_bsidx_i <- mclust::adjustedRandIndex(
        as.vector(globref_bsidx_i$cluster),
        as.vector(bsidx_i_df$cluster)
      )
      
      
      ## next, get jaccard
      
      
      #first , match up clusters; then compute jaccard of all clusters; finally, get mean jaccard
      
      
      #we need to match up the clusters...
      ## UPDATE 2025_10_31 --> this line caused an issue, if some small full cluster was totally excluded from a bootstrap?
      # fullclusts <- stringr::str_sort(unique(globref_bsidx_i$cluster), numeric = T)
      fullclusts <- stringr::str_sort(unique(globref$cluster), numeric = T)
      bsclusts <- stringr::str_sort(unique(bsidx_i_df$cluster), numeric = T)
      
      fci = fullclusts[1]
      jaccard_l <- lapply(fullclusts, function(fci){
        
        #i think this can fail... let's wrap in a try
        jaccardscore <- NA
        try({
          
          #subset the full res 
          globref_bsidx_i_fci = globref_bsidx_i[globref_bsidx_i$cluster == fci,]
          
          
          #try see which bootstrap cluster this is
          bci = bsclusts[1]
          bstabmatch_l <- lapply(bsclusts, function(bci){
            bsidx_i_df_bci = bsidx_i_df[bsidx_i_df$cluster == bci,]
            tabmatch = table(factor(bsidx_i_df_bci$barcode %in% globref_bsidx_i_fci$barcode, levels = c('FALSE', 'TRUE')))
            data.frame(BootstrapClust = bci,
                       Match = tabmatch['TRUE'],
                       NoMatch = tabmatch['FALSE'], row.names = paste0('BSclust.', bci))
          })
          bstabmatch <- dplyr::bind_rows(bstabmatch_l)
          
          #select max as the best match...
          bsmatch <- bstabmatch[which.max(bstabmatch$Match),]
          
          #jaccard of this full cluster and the matching bootstrap cluster
          fci_barcodes <- globref_bsidx_i_fci$barcode
          bci_barcodes <- bsidx_i_df[bsidx_i_df$cluster == bsmatch$BootstrapClust, "barcode"]
          
          jaccardscore <- length(intersect(fci_barcodes, bci_barcodes)) / length(union(fci_barcodes, bci_barcodes))
          
          
          
          
        }) #close try loop
        
        
        jaccardscore
      }) # close per-full cluster lapply
      
      #unlist, remove NAs, get mean; if all NA, just return 0?
      jaccard_vec <- unlist(jaccard_l)
      jaccard_vec_NONA <- na.omit(jaccard_vec)
      mean_jaccard_score <- ifelse(length(jaccard_vec_NONA)==0,
                                   yes=0,
                                   no=mean(jaccard_vec_NONA))
      
      #also get per-cluster jaccard
      perclust_jaccard_df <- data.frame(Jaccard = jaccard_vec, row.names = fullclusts)
      colnames(perclust_jaccard_df) <- bsidx_i
      
      
      
      
      
      ### prepare an output dataframe for this parameter combo
      ari_meanjaccard_df <- data.frame(params_i = params_i,
                                       bsidx_i = bsidx_i,
                                       ARI = ari_bsidx_i,
                                       MeanJaccard_PerClust = mean_jaccard_score,
                                       nClust_Ref = length(fullclusts),
                                       nClust_Bootstrap = length(bsclusts)
      )
      bsidx_i_outlist <- list(ari_meanjaccard_df = ari_meanjaccard_df,
                              perclust_jaccard_df = perclust_jaccard_df)
      
      
      bsidx_i_outlist
      
      
    }) # close Filename loop, break out to per-param loop
    
    #bind the per-barcode res
    ari_meanjaccard_df_l <- dplyr::bind_rows(lapply(thisparam_bsres_l,function(bsidx_i_outlist){bsidx_i_outlist[[1]]}))
    rownames(ari_meanjaccard_df_l) <- NULL
    
    # compute the mean of ARI and MeanJaccard_PerClust - ie, of this param combo across bootstrap
    parammeandf <- data.frame(params_i = params_i, 
                              nClust_Ref = ari_meanjaccard_df_l$nClust_Ref[1],
                              nClust_Bootstrap_Mean = mean(ari_meanjaccard_df_l$nClust_Bootstrap),
                              nClust_Bootstrap_SD = sd(ari_meanjaccard_df_l$nClust_Bootstrap),
                              
                              ARI_mean = mean(ari_meanjaccard_df_l$ARI),
                              ARI_sd = sd(ari_meanjaccard_df_l$ARI),
                              Jaccard_mean_of_clustermeans = mean(ari_meanjaccard_df_l$MeanJaccard_PerClust),
                              Jaccard_sd_of_clustermeans = sd(ari_meanjaccard_df_l$MeanJaccard_PerClust)
                              
    )
    
    
    #also, bind the per-cluster jaccard res
    perclust_jaccard_df <- dplyr::bind_cols(lapply(thisparam_bsres_l,function(bsidx_i_outlist){bsidx_i_outlist[[2]]}))
    #make it long
    perclust_jaccard_long <- perclust_jaccard_df %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(cols = dplyr::starts_with('Bootstrap')) %>%
      as.data.frame()
    
    #add cluster, bootstrap and param info
    colnames(perclust_jaccard_long)[1] <- 'cluster'
    
    perclust_jaccard_long$Bootstrap <- stringr::str_split_fixed(perclust_jaccard_long$name, '\\.', 2)[,1]
    perclust_jaccard_long$Param <- stringr::str_split_fixed(perclust_jaccard_long$name, '\\.', 2)[,2]
    
    
    #prep a simple per-cluster mean
    perclust_jaccard_mean = data.frame(cparams_i = params_i,
                                       cluster = rownames(perclust_jaccard_df), 
                                       perclust_jaccard_mean = base::rowMeans(perclust_jaccard_df),
                                       perclust_jaccard_sd = apply(perclust_jaccard_df, 1, sd))
    
    
    
    
    param_outlist <- list(ari_meanjaccard_df_l=ari_meanjaccard_df_l, #per bootstrap ARI and mean Jaccard of this param
                          parammeandf = parammeandf, #mean (and sd) across bootstraps of ARI and jaccard of this param
                          perclust_jaccard_long = perclust_jaccard_long, #per-cluster jaccard for each bootstrap of this param
                          perclust_jaccard_mean = perclust_jaccard_mean #per-cluster mean and sd jaccard from across bootstraps
    )
    
    
    param_outlist
    
    message(sprintf("[%s] finished iter %s", Sys.time(), params_i))
    
    #})
  }
  
  
  parallel::stopCluster(cl)
  names(param_outlist_l) <- params_i
  
  
  #per bootstrap ARI and mean Jaccard of each param -> kind of raw data, can be used for plotting
  ari_meanjaccard_df_l <- lapply(param_outlist_l, function(param_outlist){
    param_outlist$ari_meanjaccard_df_l
  })
  ari_meanjaccard_df <- dplyr::bind_rows(ari_meanjaccard_df_l)
  
  
  
  
  
  #mean (and sd) across bootstraps of ARI and jaccard of this param
  parammeandf_l <- lapply(param_outlist_l, function(param_outlist){
    param_outlist$parammeandf
  })
  parammeandf <- dplyr::bind_rows(parammeandf_l)
  
  
  
  
  
  #per-cluster jaccard for each bootstrap of each param
  perclust_jaccard_long_l <- lapply(param_outlist_l, function(param_outlist){
    param_outlist$perclust_jaccard_long
  })
  perclust_jaccard_long <- dplyr::bind_rows(perclust_jaccard_long_l)
  
  
  
  
  
  
  ##per-cluster mean and sd jaccard from across bootstraps for each param
  perclust_jaccard_mean_l <- lapply(param_outlist_l, function(param_outlist){
    param_outlist$perclust_jaccard_mean
  })
  perclust_jaccard_mean <- dplyr::bind_rows(perclust_jaccard_mean_l)
  rownames(perclust_jaccard_mean) <- NULL
  
  
  
  
  # parammeandf - for each param combo, the ARI and jaccard score.
  # we can compute a weighted score of these if we want...
  parammeandf$combinedscore <- (0.5 * parammeandf$ARI_mean) + (0.5 * parammeandf$Jaccard_mean_of_clustermeans)
  
  
  #add the max
  parammeandf$MAX = ''
  parammeandf[which.max(parammeandf$combinedscore),'MAX'] <- '*'
  
  
  
  #give it a better name!
  
  #2. mean (and sd) across bootstraps of ARI and jaccard of this param
  perbootstrap_perparam_scores = ari_meanjaccard_df
  
  #1. mean (and sd) across bootstraps of ARI and jaccard of this param
  perparam_meanscores = parammeandf
  
  #4. per-cluster jaccard for each bootstrap of each param
  perclust_jaccard_per_bootstrap_long = perclust_jaccard_long
  
  #3. per-cluster mean and sd jaccard from across bootstraps for each param
  perclust_jaccard_mean_acrossbootstraps = perclust_jaccard_mean
  
  
  
  
  
  fulloutlist <- list(
    
    #1. mean (and sd) across bootstraps of ARI and jaccard of each param
    perparam_meanscores=perparam_meanscores, 
    
    #2. mean (and sd) across bootstraps of ARI and jaccard of each param
    perbootstrap_perparam_scores = perbootstrap_perparam_scores,
    
    #3. per-cluster mean and sd jaccard from across bootstraps for each param
    perclust_jaccard_mean_acrossbootstraps = perclust_jaccard_mean_acrossbootstraps,
    
    #4. per-cluster jaccard for each bootstrap of each param
    perclust_jaccard_mean_acrossbootstraps = perclust_jaccard_mean_acrossbootstraps
    
  )
  
  
  i = 1
  lapply(1:length(fulloutlist), function(i){
    
    objname <- names(fulloutlist)[i]
    dfobj <- fulloutlist[[i]]
    filepath_obj <- paste0(outdir, '/', objname, '.csv')
    
    write.csv(dfobj, filepath_obj, row.names = F)
    
    objname
    
  })
  
  
  ## write to global logfile
  
  #compute final mem
  finalmem <- gc(verbose = T, full = T)
  mb <- sum(finalmem[,ncol(finalmem)])
  gb <- mb / 1000
  
  write('\n\nCLUSTER STABILITY METRIC SUMMARY:\n\n', global_logfile, append = T)
  
  write_lines_l <- list(
    
    paste0('\n\nDataframe object "perparam_meanscores"'),
    paste0('mean (and sd) across bootstraps of ARI and jaccard of this param'),
    # perparam_meanscores,
    
    
    
    paste0('\n\nMax scoring parameter:"'),
    paste0('params_i:' , perparam_meanscores[perparam_meanscores$MAX=='*','params_i']),
    paste0('ARI_mean; ', perparam_meanscores[perparam_meanscores$MAX=='*','ARI_mean']),
    paste0('Jaccard_mean_of_clustermeans; ', perparam_meanscores[perparam_meanscores$MAX=='*','Jaccard_mean_of_clustermeans']),
    paste0('combinedscore; ', perparam_meanscores[perparam_meanscores$MAX=='*','combinedscore']),
    
    
    
    paste0('\n\n'),
    # hourspassed <- (proc.time() - timestart)[3]/60/60
    paste0('Hours Passed: ', (proc.time() - timestart)[3]/60/60 ),
    
    #see gb computed above
    paste0('GB used (approximately): ', gb),
    
    paste0('\n\n')
    
    
  )
  
  
  
  lapply(write_lines_l, function(x){write(x, global_logfile, append = T)})
  
  
  
  
  
  #return fulloutlist
  return(fulloutlist)
  
  
  
}