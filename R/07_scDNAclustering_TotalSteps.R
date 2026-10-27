#' NEW_pqArmClustering
#' 
#' Perform pqArm clustering step in scDNA-seq data clustering Workflow. This 
#' function executes the pqArm clustering step in the single-cell DNA sequencing  
#' (scDNA-seq) data clustering workflow, grouping cells based on chromosomal 
#' arm-level copy number variations.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on  
#'  Giemsa-stained chromosomes. It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'
#' @return A data frame recording the pqArm clustering results for each cell, 
#'  tracking the clustering history at this stage.
#'
NEW_pqArmClustering <- function(input, pqArm_file)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  cluster <- getOption("cluster")
  
  #-------------------------- main code starts below ---------------------------
  
  message("=== Step 02: pqArm Clustering ===")

  Clustering_output <- data.frame(cluster = cluster, cell = names(input))
  pqArm_result <- NULL
  
  for (Label in unique(Clustering_output$cluster)) {
    ptm <- startTimed("pqArm Clustering ...")
    Smooth_pqCN <- NEW_pqArm_CN(input = input, 
                                Cluster_label = Label, 
                                Clustering_output = Clustering_output, 
                                pqArm_file = pqArm_file)
    pqArm_cluster <- pqArm_clustering(matrix = Smooth_pqCN, Label = Label)
    pqArm_cluster_summary <- NEW_pqArm_clustering_summary(matrix = pqArm_cluster, 
                                                          Label = Label)
    pqArm_cluster <- merge(pqArm_cluster, pqArm_cluster_summary)
    pqArm_result <- rbind(pqArm_result, pqArm_cluster)
  }

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(pqArm_result)
}


#' clusterConsolidation
#' 
#' Perform Cluster Consolidation in scDNA-seq data clustering workflow. This
#' function executes the Cluster Consolidation step in the single-cell DNA
#' sequencing (scDNA-seq) data clustering workflow, consolidating clusters based 
#' on differences in chromosomal arm-level copy number variations.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param pqArm_output A data frame recording the pqArm clustering results for 
#'  each cell. This table tracks the clustering history up to this stage.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#'
#' @return A data frame recording both pqArm clustering and consolidating results 
#'  for each cell, maintaining the clustering history across steps.
#'
clusterConsolidation <- function(input, pqArm_output, pqArm_file)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  difratio_chr <- getOption("difratio_chr")
  
  #-------------------------- main code starts below ---------------------------
  
  message("=== Step 03: Cluster consolidation ===")
  
  Recluster_Output <- NULL
  for(Label in unique(pqArm_output$cluster)) {
    ptm <- startTimed("Consolidating clusters to make cluster bigger ...")
    pqArm_merge_target <- NULL
    New_pqArm_clustering <- pqArm_recluster(pqArm_cluster = pqArm_output, 
                                            Cluster = Label)
    
    # potential merge target must have more than one (with less10 clusters
    # (>= 2 & < 10), or more than one more10 clusters(> 2))

    if (length(New_pqArm_clustering) > 1) {
      pqArm_dif <- pqArm_reclustering_dif(
                     input = input, 
                     pqArm_recluster_sim = New_pqArm_clustering,
                     pqArm_cluster = pqArm_output, 
                     Cluster = Label, 
                     pqArm_file = pqArm_file)
      pqArm_merge_target <- pqArm_reclusterBy_ratio_target(
                              pqArm_cluster = pqArm_output,
                              Cluster = Label, 
                              pqReclsut_sim = pqArm_dif, 
                              difratio_chr = difratio_chr)
    }
    # is.null(pqArm_merge_target)
    if (is.null(pqArm_merge_target) == FALSE) {
      Recluster_CellID <- pqArm_recluster_result(pqArm_cluster = pqArm_output, 
                                                 pqReclsut_target = pqArm_merge_target)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    } else {
      Recluster_CellID <- pqArm_output %>%
        dplyr::filter(.data$cluster %in% Label) %>%
        dplyr::mutate(Recluster_pattern = .data$pqArm_pattern,
                      Recluster_cellnum = .data$pqArm_cellnum,
                      Recluster_cluster = .data$pqArm_cluster)
      Recluster_Output <- rbind(Recluster_Output, Recluster_CellID)
    }
  }
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Recluster_Output)
}


#' filterCells
#'
#' Filters consolidated clustering results, fetches cluster-specific breakpoints, 
#' transforms the data into wide format, and converts non-missing values into 
#' a binary presence/absence matrix.
#'
#' @param input A named list where each element represents a cell, containing its 
#'  corresponding genomic segments as a `GRanges` object.
#' @param Consolidating_output A data frame containing recluster labels and cell 
#'  counts (must include columns `Recluster_cellnum` and `Recluster_cluster`).
#'
#' @return A data frame with cell IDs as row names, where columns represent genomic 
#'  regions (`seqnames_start_end`) and values are binary (`1` for present, `0` 
#'  for missing).
#'
#' @importFrom dplyr filter mutate bind_rows
#' @importFrom tidyr pivot_wider
#' 
filterCells <- function(input, Consolidating_output) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  min_cell_GMM <- getOption("min_cell_GMM")
  overlap_region <- getOption("overlap_region")
  
  #-------------------------- main code starts below ---------------------------
  
  binsize <- as.data.frame(input[[1]]$bins)$width[1]
  overlap_bp <- round((overlap_region/binsize), digits = 0)

  R <- Consolidating_output %>% dplyr::filter(.data$Recluster_cellnum >= min_cell_GMM)
  
  breakpoints_wide_list <- list()
  for (Recluster_label in unique(R$Recluster_cluster)) {
    # Fetch breakpoints
    bp <- collect_cluster_bp(input = input, 
                             Clustering_output = R, 
                             Recluster_label = Recluster_label)
    # Check if bp is NULL or contains no rows; skip iteration if true
    if (is.null(bp) || nrow(bp) == 0) {
      next
    }
    if (!"copy.number" %in% names(bp)) {
      bp$copy.number <- as.numeric(gsub("-somy", "", bp$state.left))
    }
    
    # Create concatenated region column and pivot
    bp_wide <- bp %>%
               mutate(new_column_name = paste(seqnames, start, end, sep = "_")) %>%
               pivot_wider(id_cols = cellID,
                           names_from = new_column_name,
                           values_from = copy.number,
                           values_fill = NA)

    # Store result in list
    breakpoints_wide_list[[as.character(Recluster_label)]] <- bp_wide
  }
  # Combine all wide data frames into one single data frame
  final_breakpoints_wide <- bind_rows(breakpoints_wide_list)
  # Convert non-NA to 1 and NA to 0
  cellID <- final_breakpoints_wide$cellID
  final_breakpoints_wide <- final_breakpoints_wide[,-1]
  breakpoints_binary <- as.data.frame(ifelse(is.na(final_breakpoints_wide), 0, 1))
  rownames(breakpoints_binary) <- cellID
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(breakpoints_binary)
}


#' Perform PCA, Mclust Clustering, and Filter DNA Data
#'
#' Filters non-zero variance columns from a binary breakpoint matrix, runs PCA 
#' on the log1p-transformed data, performs Gaussian Mixture Modeling (Mclust) 
#' on the top 30 principal components, and filters the consolidation and DNA 
#' data based on the predominant cluster assignment.
#'
#' @param breakpoints_binary A binary matrix or data frame where rows represent 
#'  cells and columns represent genomic features.
#' @param Consolidation_result A data frame containing cell metadata including 
#'  at least `cellID` and `Recluster_cluster` columns.
#' @param input A named list where each element represents a cell, containing its 
#'  corresponding genomic segments as a `GRanges` object.
#'
#' @return A list containing:
#' \item{Consolidation_output_sub_mod}{Filtered and annotated consolidation data frame.}
#' \item{DNA_data_sub_mod}{Subset of DNA data matching the filtered cell IDs.}
#' \item{mclust_model}{The fitted `Mclust` model object.}
#'
#' @importFrom stats prcomp var
#' @importFrom mclust Mclust mclustBIC
#' @importFrom dplyr mutate case_when add_count filter select group_by ungroup
#' 
GMMcluster <- function(breakpoints_binary, Consolidation_result, input) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.4_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  num_pcs <- getOption("num_pcs")
  g_clusters <- getOption("g_clusters")
  entropy_threshold <- getOption("entropy_threshold")
  
  #-------------------------- main code starts below ---------------------------
  
  # 1. Filter out zero-variance breakpoint features & run PCA + GMM
  breakpoints_binary_1 <- breakpoints_binary[, apply(breakpoints_binary, 2, 
                                                     function(x) var(x, na.rm = TRUE) != 0)]
  breakpoints_binary_1_pca <- prcomp(log1p(breakpoints_binary_1), scale. = TRUE)
  breakpoints_binary_1_PCs <- breakpoints_binary_1_pca$x[, 1:min(30, ncol(breakpoints_binary_1_pca$x))]
  
  mod <- Mclust(breakpoints_binary_1_PCs, G = 1:g_clusters) 
  
  MclustSelect <- data.frame(mod$classification)
  select_1 <- rownames(subset(MclustSelect, mod.classification %in% c("1")))
  select_2 <- rownames(subset(MclustSelect, mod.classification %in% c("2")))
  select_3 <- rownames(subset(MclustSelect, mod.classification %in% c("3")))
  select_4 <- rownames(subset(MclustSelect, mod.classification %in% c("4")))
  select_5 <- rownames(subset(MclustSelect, mod.classification %in% c("5")))
  select_6 <- rownames(subset(MclustSelect, mod.classification %in% c("6")))
  select_7 <- rownames(subset(MclustSelect, mod.classification %in% c("7")))
  select_8 <- rownames(subset(MclustSelect, mod.classification %in% c("8")))
  select_9 <- rownames(subset(MclustSelect, mod.classification %in% c("9")))
  
  Consolidation_output <- Consolidation_result %>%
                          mutate(
                            GMM_select = case_when(
                            cellID %in% select_1 ~ "cluster1",
                            cellID %in% select_2 ~ "cluster2",
                            cellID %in% select_3 ~ "cluster3",
                            cellID %in% select_4 ~ "cluster4",
                            cellID %in% select_5 ~ "cluster5",
                            cellID %in% select_6 ~ "cluster6",
                            cellID %in% select_7 ~ "cluster7",
                            cellID %in% select_8 ~ "cluster8",
                            cellID %in% select_9 ~ "cluster9",
                            TRUE                 ~ "dopped"))
  
  Consolidation_output_sub <- subset(Consolidation_output, GMM_select != "dopped")
  DNA_data_sub <- input[Consolidation_output_sub$cellID]
  
  # 2. Compute Shannon Entropy per group to find the least noisy cluster
  calc_mean_entropy <- function(cell_ids) {
    mat <- breakpoints_binary_1[rownames(breakpoints_binary_1) %in% cell_ids, , drop = FALSE]
    if (nrow(mat) == 0) return(Inf)
    
    # Proportion of cells containing a breakpoint at each genomic location
    p <- colMeans(mat, na.rm = TRUE)
    # Filter out fixed positions (p=0 or p=1 have 0 entropy)
    p <- p[p > 0 & p < 1]
    if (length(p) == 0) return(0)
    
    # Position-wise binary Shannon entropy
    entropy_per_pos <- -p * log2(p) - (1 - p) * log2(1 - p)
    mean(entropy_per_pos)
  }
  select_list <- mget(paste0("select_", 1:9))
  cluster_summary <- data.frame(
    cluster = paste0("cluster", 1:9),
    entropy = sapply(select_list, calc_mean_entropy),
    n = sapply(select_list, length)
  )
  
  # --- Pause Pipeline & Display Values ---
  message("\n============================================")
  message("         ENTROPY CHECKPOINT PAUSE           ")
  message("============================================")
  message(sprintf("%-12s | %-15s | %s", "cluster", "entropy_threshold", "n"))
  message("--------------------------------------------")
  for (i in seq_len(nrow(cluster_summary))) {
    nm <- cluster_summary$cluster[i]
    val <- cluster_summary$entropy[i]
    val_n <- cluster_summary$n[i]
    val_str <- ifelse(is.infinite(val), "NA", sprintf("%.4f", val))
    message(sprintf("%-12s | %-15s | %d", nm, val_str, val_n))
  }
  message("--------------------------------------------")
  message("Note: If entropy_threshold is set to maximum '1', no data will be filtered out from GMM.")
  message("--------------------------------------------")
  
  # Prompt user and pause execution until input is received
  input_val <- readline(prompt = "Enter your desired entropy threshold (enter 0 for default 0.2): ")
  entropy_threshold <- suppressWarnings(as.numeric(input_val))
  
  # Handle default fallback if input is 0, empty, or invalid
  if (is.na(entropy_threshold) || entropy_threshold == 0) {
    entropy_threshold <- 0.2
    message("-> Resuming pipeline with default entropy threshold: 0.2\n")
  } else {
    message(paste("-> Resuming pipeline with user entropy threshold:", entropy_threshold, "\n"))
  }
  # --- Resume Pipeline Excution ---
  entropy_values <- cluster_summary$entropy
  dropped_group <- entropy_values[
    entropy_values > entropy_threshold | is.infinite(entropy_values)]
  
  ck <- cluster_summary[cluster_summary$entropy %in% dropped_group,]
  dropped_group <- ck$cluster
  
  if (length(dropped_group) == 0) {
    dropped_group <- NULL
  }
  Consolidation_output_sub_mod <- Consolidation_output_sub %>% 
    filter(!GMM_select %in% dropped_group) %>% 
    group_by(Recluster_cluster) %>%
    mutate(GMM_cellnum = n()) %>%
    ungroup()
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Consolidation_output_sub_mod)
}


#' NEW_SubClustering
#' 
#' Perform subclustering in scDNA-seq data clustering workflow. This function 
#' executes the subclustering step in the single-cell DNA sequencing (scDNA-seq) 
#' data clustering workflow, refining clusters based on copy number variations 
#' (CNVs) within subpopulations of cells.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Consolidation_output A data frame recording pqArm clustering and consolidating 
#'  results for each cell. This table tracks the clustering history up to this stage.
#' @param min_cell_subclone Minimum cell count required to define a valid subclone.
#' @param min_cell_DNA An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#'
#' @return A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, 
#'  consolidating, and subclustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template 
#'  and its copy number. It includes the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the segment.
#'   - `end`: End position of the segment.
#'   - `region`: The defined region index
#'   - `Subclone`: Subclone identifier.
#'   - `CN`: Copy number value of the segment.
#'
NEW_SubClustering <- function(input, Consolidation_output, min_cell_DNA,
                              min_cell_subclone)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.5_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  overlap_region <- getOption("overlap_region")
  dif_ratio <- getOption("dif_ratio")
  
  #-------------------------- main code starts below ---------------------------
  
  message("=== Step 04: Subclustering ===")
  
  Final_output <- list()
  Subclone_cluster <- NULL
  Subclone_CNr <- NULL
  cell_clustering <- NULL

  # from overlap_region to number of bins in coverage
  binsize <- as.data.frame(input[[1]]$bins)$width[1]
  overlap_bp <- round((overlap_region/binsize), digits = 0)
  
  R <- Consolidation_output %>% dplyr::filter(.data$Recluster_cellnum >= min_cell_DNA)
  
  for (Recluster_label in unique(R$Recluster_cluster)) {
    message("Subclustering for Recluster ", Recluster_label, " ... \n")
    breakpoints <- collect_cluster_bp(input = input, 
                                      Clustering_output = R, 
                                      Recluster_label = Recluster_label)
    if (is.null(breakpoints)) {
      warning("Breakpoints returned NULL for Recluster label: ", Recluster_label, 
              ". Skipping iteration.", call. = FALSE)
      next
    }
    Cell_num <- R %>% 
                dplyr::filter(.data$Recluster_cluster == Recluster_label) %>% 
                dplyr::pull(.data$Recluster_cellnum)
    bp <- output_bp_covers(Template = breakpoints, 
                           binsize = binsize, 
                           overlap = overlap_bp, 
                           # overlap_times = 0.8*cell
                           overlap_times = Cell_num[1] * 0.8)
    if (is.null(bp)) { 
      print(paste0("bp is NULL on cluser: ", R$Recluster_label))
      next 
    }  
    event_region <- bp_events(input = input, Template = bp, binsize = binsize)
    consensus_bp_template <- bp_region(event = event_region, binsize = binsize)
    cell_CNregion <- Region_CN(input = input, 
                               Reclustering_output = Consolidation_output, 
                               Recluster_label = Recluster_label, 
                               events = consensus_bp_template)
    
    if (is.null(cell_clustering) == TRUE) {
      cell_clustering <- Subclone_clustering(CN_incells_input = cell_CNregion, 
                                             event_region = consensus_bp_template,
                                             dif_ratio = dif_ratio, 
                                             Subclone_num = 0)
    } else {
      Subclone_num = max(Subclone_cluster$Subclone)
      cell_clustering <- Subclone_clustering(CN_incells_input = cell_CNregion, 
                                             event_region = consensus_bp_template,
                                             dif_ratio = dif_ratio, 
                                             Subclone_num = Subclone_num)
    }
    Subclone_cluster <- rbind(Subclone_cluster, cell_clustering)
    # Output: Copy number in each region in each subclone
    
    s_CN <- NEW_Subclone_CNregion(sep_region = consensus_bp_template,
                                  CN_region = cell_CNregion, 
                                  each_subclone = cell_clustering,
                                  min_cell_subclone = min_cell_subclone,
                                  #Two options:"SubcloneCNVRegion", "SubcloneRegionCN" 
                                  #set in function 5.7
                                  output = "SubcloneRegionCN")
    
    Subclone_CNr <- rbind(Subclone_CNr, s_CN)
    
    rm(breakpoints, Cell_num, bp, event_region, consensus_bp_template, s_CN )
  }
  colnames(Subclone_cluster) <- c("Subclone_cluster", "cellID", "Subclone_cellnum")
  Final_output <- list(final_cluster_output = dplyr::left_join(Consolidation_output,
                                                               Subclone_cluster, 
                                                               by = "cellID"),
                       Subclone_CN = Subclone_CNr)

  #---------------------------- main code ends here ----------------------------

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(Final_output)
}


#' scDNA_Output_s3
#' 
#' Outputting scDNA clustering final outputs. Generate final output in scDNA-seq 
#' data clustering workflow. This function produces the final output for the 
#' single-cell DNA sequencing (scDNA-seq) clustering workflow, integrating clustering 
#' results and generating visualizations.
#'
#' @param input A named list where each element is a `GRanges` object representing 
#'  a single cell.
#' @param Summary A list containing two data frames:
#'   - `final_cluster_output`: Records the clustering history from the pqArm, 
#'  consolidating, and subclustering steps.
#'   - `Subclone_CN`: Records each subclone's unique chromosome segment template 
#'  and its copy number. It includes the following columns:
#'   - `chr`: Chromosome name (chr1, chr2, ...).
#'   - `start`: Start position of the segment.
#'   - `end`: End position of the segment.
#'   - `region`: The defined region index
#'   - `Subclone`: Subclone identifier.
#'   - `CN`: Copy number value of the segment.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param output_dir A character string specifying the directory where the output 
#'  files will be saved.
#' @param min_cell_DNA An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#'
#' @return The function generates final clustering results and visualizations, 
#'  saving them to the specified output directory.
#'
scDNA_Output_s3 <- function(input, Summary, pqArm_file, output_dir, min_cell_DNA, 
                            cnv_ratio_DNA)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.6_s3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)

  #-------------------------- main code starts below ---------------------------

  message("=== Step 05: Output cnvTree results ===")
  timestamp <- format(Sys.time(), "%m%d_%H")
  # 1. cellID summary
  filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_grouping")
  writeOutput(data = Summary$final_cluster_output, 
              filename = filename, 
              path = output_dir)
  message("Output /cnvTree.scDNAseq_grouping.txt is done.")

  # 2. Each region copy number to subclone
  filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_SubcloneRegionCN")
  writeOutput(data = Summary$Subclone_CN, 
              filename = filename, 
              path = output_dir)
  message("Output /cnvTree.scDNAseq_SubcloneRegionCN.txt is done.")

  # 3. All CNV region in the sample
  CNV_Data <- Total_cnvRegion_s3(input = input, #input = input_smooth
                                 Template = Summary$Subclone_CN, 
                                 pqArm_file = pqArm_file)
  if (nrow(CNV_Data) != 0){
    CNV_Data <- cnvRegion.toPQarm(FILE = CNV_Data, pqArm_file = pqArm_file)
    filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_DefinedCNVregion")
    writeOutput(data = CNV_Data, 
                filename = filename, 
                path = output_dir)
    message("Output /cnvTree.scDNAseq_DefinedCNVregion.txt is done.")
  } else {
    message("Warning: Based on the length of consecutive_region, 
             no CNV regions remain in the sample. \n
             Consider lowering the consecutive_region parameter. \n
             The file cnvTree.scDNAseq_DefinedCNVregion.txt and 
             cnvTree.scDNAseq_DefinedCNVregion.txt was not generated.")
  }

  # 4. DNA superimpose
  if (nrow(CNV_Data) != 0){
    DNA_superimpose <- scDNA.superimpose_s3(Template = Summary, DefinedCNVs = CNV_Data,
                                            cnv_ratio_DNA = cnv_ratio_DNA)
    DNA_superimpose <- scDNA.clustering(Template = DNA_superimpose,
                                        min_cell_DNA = min_cell_DNA)
    filename = paste0("/", timestamp, "_", "cnvTree.scDNAseq_DNAcluster")
    writeOutput(data = DNA_superimpose$DNA_cluster, 
                filename = filename, 
                path = output_dir)
    message("Output /cnvTree.scDNAseq_DNAcluster.txt is done.")
  } else {
    message("Warning: Based on the length of consecutive_region, 
             no CNV regions remain in the sample. \n
             Consider lowering the consecutive_region parameter. \n
             The file cnvTree.scDNAseq_DNAcluster.txt and 
             cnvTree.scDNAseq_DefinedCNVregion.txt were not generated.")
  }
  # 5. Final cluster heatmap
  FILEname = paste0("/", timestamp, "_", "cnvTree.scDNAseq_Grouping_fig.pdf")
  Totalcluster_pdf(input = input, 
                   Template = Summary$final_cluster_output,
                   pqArm_file = pqArm_file,
                   min_cell_DNA = min_cell_DNA,
                   step = "Subclone",
                   FILEpath = output_dir, 
                   FILEname = FILEname)
  message("Output /cnvTree.scDNAseq_Grouping_fig.pdf is done.")

  # 6. cluster heatmap with dendrogram
  FILEname = paste0("/", timestamp, "_", "cnvTree.scDNAseq_heatmap.png")
  scDNA_CNVpattern(input = input, 
                   final_cluster = Summary$final_cluster_output, 
                   min_cell_DNA = min_cell_DNA,
                   pqArm_file = pqArm_file, 
                   FILEpath = output_dir, 
                   FILEname = FILEname)
  message("Output /cnvTree.scDNAseq_heatmap.png is done.")

  print(paste("cnvTree Results were saved at", output_dir))
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(DNA_superimpose)
}


#' cnvTree_scDNAclustering
#' 
#' One-Step scDNA-seq Cell Clustering Pipeline. This function executes the entire 
#' single-cell DNA sequencing (scDNA-seq) clusteringworkflow in one step, including 
#' NEW_pqArm_clustering, consolidating, subclustering, and final CNV-based output 
#' generation.
#'
#' @param input_dir_DNA A named list where each element is a `GRanges` object 
#'  representing a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param output_dir A character string specifying the directory where output 
#'  files will be saved.
#' @param min_cell_DNA An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#' @param min_cell_subclone Minimum cell count required to define a valid subclone.
#' @param GMM_mode Logical. If \code{TRUE}, the pipeline incorporates the Gaussian 
#'  Mixture Model (GMM) calculation.
#'
#' @return The function performs complete clustering and CNV analysis, saving final 
#'  results and visualizations in the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config_path <- system.file("extdata", "cnvTree_config.yaml", package = "cnvTree")
#' config <- read_yaml(config_path)
#' selected <- run_cnvTree_Pipeline(output_dir = config$output_dir,
#'                                  input_dir_RNA = config$input_dir_RNA) 
#' cnvTree_scDNAclustering(input_dir_DNA = config$input_dir_DNA,
#'                         cnv_ratio_DNA = config$cnv_ratio_DNA,
#'                         output_dir = config$output_dir,
#'                         pqArm_file = config$pqArm_file,
#'                         min_cell_DNA = config$min_cell_DNA,
#'                         min_cell_subclone = config$min_cell_subclone,
#'                         GMM_mode = config$GMM_mode)
#' }
#'
cnvTree_scDNAclustering <- function(input_dir_DNA, pqArm_file, output_dir,
                                    cnv_ratio_DNA, min_cell_DNA, min_cell_subclone,
                                    GMM_mode) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.7_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  difratio_chr <- getOption("difratio_chr")
  overlap_region <- getOption("overlap_region")
  dif_ratio <- getOption("dif_ratio")
  smoothheatmap <- getOption("smoothheatmap")
  sexchromosome <- getOption("sexchromosome")
  consecutive_region <- getOption("consecutive_region")
  
  #-------------------------- main code starts below ---------------------------
  
  ptm <- startTimed("Start cnvTree_scDNAclustering...")
  
  input <- ProcessHmmList(input_dir_DNA = config$input_dir_DNA)
  pqArm_result <- NEW_pqArmClustering(input = input, pqArm_file = config$pqArm_file)
  Consolidation_result <- clusterConsolidation(input = input, 
                                               pqArm_output = pqArm_result,
                                               pqArm_file = config$pqArm_file)
  if (GMM_mode == FALSE) {
    # extracting CNVs (main algorithm)
    Subclone_output <- NEW_SubClustering(input = input,
                                         Consolidation_output = Consolidation_result,
                                         min_cell_DNA = min_cell_DNA,
                                         min_cell_subclone = min_cell_subclone)
    scDNA_Output_s3(input = input,
                    Summary = Subclone_output,
                    pqArm_file = pqArm_file, 
                    output_dir = output_dir,
                    min_cell_DNA = min_cell_DNA,
                    cnv_ratio_DNA = cnv_ratio_DNA)
    
  } else if (GMM_mode == TRUE) {
    
    # GMM Removing poor quality cells
    breakpoints_binary <- filterCells(input = input, 
                                      Consolidating_output = Consolidation_result)
    Consolidation_result_GMM <- GMMcluster(input = input, 
                                           breakpoints_binary = breakpoints_binary, 
                                           Consolidation_result = Consolidation_result)
    
    ## Pass GMM-selected data to the main algorithm
    Subclone_output_GMM <- NEW_SubClustering(input = input,
                                             Consolidation_output = Consolidation_result_GMM,
                                             min_cell_DNA = min_cell_DNA,
                                             min_cell_subclone = min_cell_subclone)
    
    ## New ".scDNAseq_grouping.txt"
    df_export <- data.frame(Subclone_output_GMM$final_cluster_output)
    Consolidation_result_sub <- subset(Consolidation_result, !(cellID %in% df_export$cellID))
    cols_add <- c("GMM_select", "GMM_cellnum", "Subclone_cluster", "Subclone_cellnum")
    Consolidation_result_sub[cols_add] <- NA
    Subclone_output_GMM$final_cluster_output <- rbind(Consolidation_result_sub, df_export)
    
    scDNA_Output_s3(input = input,
                    Summary = Subclone_output_GMM,
                    pqArm_file = pqArm_file, 
                    output_dir = output_dir,
                    min_cell_DNA = min_cell_DNA,
                    cnv_ratio_DNA = cnv_ratio_DNA)
  } 
  endTimed(ptm)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}


#' cnvTree_scDNAclustering_df
#' 
#' One-Step scDNA-seq Cell Clustering Pipeline. This function executes the entire 
#' single-cell DNA sequencing (scDNA-seq) clustering workflow in one step, including 
#' pqArm clustering, consolidating, subclustering, and final CNV-based output 
#' generation.
#'
#' @param input_dir_DNA A named list where each element is a `GRanges` object 
#'  representing a single cell.
#' @param pqArm_file In-build cytoband template for selection: `hg38`, `hg19`, 
#'  `mm10`, `mm39`. Or a filepath of a table for cytoband information seen on 
#'  Giemsa-stained chromosomes. It should include the following columns:
#'   - `chrom`: Reference sequence chromosome or scaffold.
#'   - `chromStart`: Start position in genoSeq.
#'   - `chromEnd`: End position in genoSeq.
#'   - `name`: Name of cytogenetic band.
#'   - `gieStain`: Giemsa stain results.
#' @param output_dir A character string specifying the directory where output 
#'  files will be saved.
#' @param min_cell_DNA An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#' @param min_cell_subclone Minimum cell count required to define a valid subclone.
#' @param cnv_ratio_DNA Numeric threshold (between 0 and 1) representing the 
#'  minimum overlap ratio required to classify a defined CNV region as present 
#'  (1) versus absent (0).
#' @param GMM_mode Logical. If \code{TRUE}, the pipeline incorporates the Gaussian 
#'  Mixture Model (GMM) calculation
#'
#' @return The function performs complete clustering and CNV analysis, saving final 
#'  results and visualizations in the specified output directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config_path <- system.file("extdata", "cnvTree_config.yaml", package = "cnvTree")
#' config <- read_yaml(config_path)
#' selected <- run_cnvTree_Pipeline(output_dir = config$output_dir,
#'                                  input_dir_RNA = config$input_dir_RNA) 
#' cnvTree_scDNAclustering_df(input_dir_DNA = config$input_dir_DNA,
#'                            cnv_ratio_DNA = config$cnv_ratio_DNA,
#'                            output_dir = config$output_dir,
#'                            pqArm_file = config$pqArm_file,
#'                            min_cell_DNA = config$min_cell_DNA,
#'                            min_cell_subclone = config$min_cell_subclone,
#'                            GMM_mode = config$GMM_mode)
#' }                         
#'
cnvTree_scDNAclustering_df <- function(input_dir_DNA, pqArm_file, output_dir,
                                       min_cell_DNA, min_cell_subclone, GMM_mode,  
                                       cnv_ratio_DNA) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 7.8_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  message("Start data transform...")
  ptm <- startTimed("Start cnvTree_scDNAclustering...")
  
  input_df <- NEW_changeFormat(input_dir_DNA = input_dir_DNA)
  pqArm_result <- NEW_pqArmClustering(input = input_df, pqArm_file = pqArm_file) 
  Consolidation_result <- clusterConsolidation(input = input_df, 
                                               pqArm_output = pqArm_result,
                                               pqArm_file = pqArm_file)
  
  if (GMM_mode == FALSE) {
    # extracting CNVs (main algorithm)
    Subclone_output <- NEW_SubClustering(input = input_df,
                                         Consolidation_output = Consolidation_result,
                                         min_cell_DNA = min_cell_DNA,
                                         min_cell_subclone = min_cell_subclone)
    scDNA_Output_s3(input = input_df,
                    Summary = Subclone_output,
                    pqArm_file = pqArm_file, 
                    output_dir = output_dir,
                    min_cell_DNA = min_cell_DNA,
                    cnv_ratio_DNA = cnv_ratio_DNA)
    
  } else if (GMM_mode == TRUE) {
    # GMM Removing poor quality cells
    breakpoints_binary <- filterCells(input = input_df, 
                                      Consolidating_output = Consolidation_result)
    Consolidation_result_GMM <- GMMcluster(input = input_df, 
                                           breakpoints_binary = breakpoints_binary, 
                                           Consolidation_result = Consolidation_result)

    ## Pass GMM-selected data to the main algorithm
    Subclone_output_GMM <- NEW_SubClustering(input = input_df,
                                             Consolidation_output = Consolidation_result_GMM,
                                             min_cell_DNA = min_cell_DNA,
                                             min_cell_subclone = min_cell_subclone)
    
    ## New ".scDNAseq_grouping.txt"
    df_export <- data.frame(Subclone_output_GMM$final_cluster_output)
    Consolidation_result_sub <- subset(Consolidation_result, !(cellID %in% df_export$cellID))
    cols_add <- c("GMM_select", "GMM_cellnum", "Subclone_cluster", "Subclone_cellnum")
    Consolidation_result_sub[cols_add] <- NA
    Subclone_output_GMM$final_cluster_output <- rbind(Consolidation_result_sub, df_export)
    
    scDNA_Output_s3(input = input_df,
                    Summary = Subclone_output_GMM,
                    pqArm_file = pqArm_file, 
                    output_dir = output_dir,
                    min_cell_DNA = min_cell_DNA,
                    cnv_ratio_DNA = cnv_ratio_DNA)
  } 
  endTimed(ptm)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}

