#' changeFormat
#'
#' Convert copy number data to GRanges Format. This function reads a `.rds` file 
#' or `.txt` file containing copy number variation (CNV) data and converts it into 
#' a list of `GRanges` objects, where each element corresponds to a single cell.
#'
#' @param cores An integer specifying the number of worker cores to use for parallel 
#'  processing.
#' @param input_dir_DNA A `.rds` file or `.txt` file containing a data frame with 
#'  the following required columns:
#'  
#'    - `cellID`: Unique identifier for each cell.
#'    - `seqnames`: Chromosome or sequence name.
#'    - `start`: Start position of the segment.
#'    - `end`: End position of the segment.
#'    - `copy.number`: Copy number value for the segment.
#'   
#' @return A named list where each element represents a cell, containing its 
#'  corresponding genomic segments as a `GRanges` object.
#'
changeFormat <- function(input_dir_DNA, cores) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.0_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  sexchromosome <- getOption("sexchromosome")
  
  #---------------------- main code start below ---------------------
  
  if (endsWith(input_dir_DNA, ".rds") == TRUE) {
    Bin_CN <- readRDS(input_dir_DNA)
  } else if(endsWith(input_dir_DNA, ".txt") == TRUE) {
    Bin_CN <- utils::read.delim2(input_dir_DNA, sep = " ")
  } else {
    "Please input either .rds or .txt format as input."
    stop(changeFormat)
  }

  if (sexchromosome == FALSE) {
    rows_to_remove <- grepl("chrX|chrY|chrM", Bin_CN$seqnames)
    Bin_CN <- Bin_CN[!rows_to_remove, ]
  }
  data.table::setDT(Bin_CN)
  clean_chroms <- unique(str_remove(Bin_CN$seqnames, "(?i)chr"))
  is_numeric <- !is.na(suppressWarnings(as.numeric(clean_chroms)))
  sorted_nums <- sort(as.numeric(clean_chroms[is_numeric]))
  sorted_sex  <- sort(clean_chroms[!is_numeric])
  correct_level_order <- c(as.character(sorted_nums), sorted_sex)
  Bin_CN <- Bin_CN %>%
            mutate(seqnames_sort = factor(str_remove(seqnames, "(?i)chr"), 
                   levels = correct_level_order)) %>%
            arrange(seqnames_sort)
  
  # Set up parallel environment
  future::plan(future::multisession, workers = cores)
  
  # Split data by "cellID"
  Bin_CN_list <- split(Bin_CN, by = "cellID", keep.by = FALSE)
  template_data <- Bin_CN_list[[1]][,setdiff(names(Bin_CN_list[[1]]), 
                                                   "copy.number"), 
                                                   with = FALSE]
  template_gr <- GenomicRanges::makeGRangesFromDataFrame(template_data, 
                                                         keep.extra.columns = FALSE)
  NewFormat <- list()
  NewFormat <- lapply(seq_along(Bin_CN_list), function(i) {
    cell_data <- Bin_CN_list[[i]]
    Bins <- template_gr
    S4Vectors::mcols(Bins)$copy.number <- cell_data$copy.number
    # Compute breakpoints
    shifted_cn <- c(NA, cell_data$copy.number[-nrow(cell_data)]) # 向前平移
    breakpoint_rows <- cell_data$copy.number != shifted_cn & !is.na(shifted_cn)
    breakpoints <- cell_data[breakpoint_rows, 
                             c("seqnames", "start", "end", "copy.number"), 
                             with = FALSE]
    # Convert breakpoints to GRanges if not empty
    Breakpoints <- if (nrow(breakpoints) == 0) {
      NULL
    } else {
      GenomicRanges::GRanges(seqnames = breakpoints$seqnames,
                             ranges = IRanges::IRanges(start = breakpoints$start, 
                                                       end = breakpoints$end),
                             copy.number = breakpoints$copy.number)
    }
    list(ID = names(Bin_CN_list)[i], bins = Bins, breakpoints = Breakpoints)
  })
  names(NewFormat) <- names(Bin_CN_list)
  on.exit(future::plan(future::sequential), add = TRUE)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(NewFormat)
}


#' NEW_CN_seq
#'
#' Generate copy number matrix for selected cells; this function extracts copy 
#' number variations from a list of `GRanges` objects and organizes them into an 
#' integer matrix. The matrix contains selected cells as columns, with genomic 
#' regions (fixed bins) as rows.
#'
#' @param input A named list where each element represents a single cell as a 
#'  `GRanges` object.
#' @param Template A character vector containing the `cellID`s of selected cells 
#'  to be included in the matrix.
#'
#' @return An integer matrix:
#'  - Columns represent the selected `cellID`s.
#'  - Rows represent genomic regions, separated into fixed bins.
#'
NEW_CN_seq <- function(input, Template)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #---------------------- main code start below ---------------------
  
  c <- lapply(Template, function(i) {
    obj <- input[[i]]
    if (!is.null(obj) && "bins" %in% names(obj) && 
        "copy.number" %in% names(S4Vectors::mcols(obj[["bins"]]))) {
      S4Vectors::mcols(obj[["bins"]])$copy.number
    } else {
      warning(sprintf("Missing data or 'copy.number' column for template '%s'", i))
      NA_real_  # Use real if copy numbers are decimals
    }
  })
  c <- as.data.frame(c, stringsAsFactors = FALSE) #避免字符串變為因子
  names(c) <- Template
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(c)
}


#' run_cnvTree_Pipeline
#' 
#' This initializes the output environment, sets up necessary folder structures
#' for RNA data, and extracts the available cell groups/types for selection.
#' 
#' @param output_dir Character. Path to the directory where results and output 
#'  folders will be created and saved.
#' @param input_dir_RNA Character. Path to the input single-cell RNA-seq 
#'  data directory (e.g., inferCNV output folder).
#'
#' @return A character vector or list of identified cell groups/types extracted 
#'  by \code{select_groups()} for the next steps of the pipeline.
#' 
#' @export
#' 
run_cnvTree_Pipeline <- function(output_dir, input_dir_RNA) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  RNAdataSource <- getOption("RNAdataSource")

  #---------------------- main code start below ---------------------
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)
    message("Created output folder: ", output_dir)
  }
  
  if (!is.null(input_dir_RNA)) {
    selected_groups <- select_groups(input_dir_RNA = input_dir_RNA, 
                                     RNAdataSource = RNAdataSource)
  } else {
    selected_groups <- NULL
    message("No RNA directory provided. Skipping RNA data integration.")
  }
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(selected_groups)
}


#' select_groups
#'
#' Selects Cell Groups Interactively from CNV region File; this function scans a 
#' specified directory for CNV prediction files, extracts unique cell group names, 
#' and provides an interactive menu for the user to select one or more groups for 
#' downstream analysis.
#'
#' User Flow:
#' 1. Select an initial group. If '0' is pressed, the function execution stops.
#' 2. A repeat loop allows for multiple additional selections.
#'
#' @param input_dir_RNA A character string specifying the path to the directory  
#'   containing the CNV results (specifically looking for files containing 
#'   "cnv_regions" and "HMM_CNV_predictions").
#' @param RNAdataSource Integer or Character. The type of RNA data source being used:
#'   \itemize{
#'     \item \code{1}: Tabular data file / matrix structure.
#'     \item \code{2}: InferCNV output folder structure.
#'   }
#'
#' @return A character vector containing the unique names of the selected cell 
#' groups.
#' 
select_groups <- function(input_dir_RNA, RNAdataSource)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.2.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #---------------------- main code start below ---------------------

  if (RNAdataSource == "1") { # dataframe input
    print("datatype: dataframe")
    ext <- tools::file_ext(input_dir_RNA)
    if (ext == "rds") {
      cnv_regions <- readRDS(input_dir_RNA)
    } else if (ext == "txt") {
      cnv_regions <- read.table(input_dir_RNA, sep = " ", header = TRUE) 
    } else {
      stop("Unsupported format. Please provide a .rds, .txt, or .tsv file.")
    }
    
  } else if (RNAdataSource == "2") { # InferCNV input
    print("datatype: InferCNV files")
    matching_files <- list.files(path = input_dir_RNA, 
                                 pattern = "HMM_CNV_predictions.*cnv_regions", 
                                 recursive = TRUE, full.names = TRUE)
    if (length(matching_files) == 0) {
      stop("No matching CNV prediction files found in: ", input_dir_RNA)
    }
    if (length(matching_files) > 1) {
      # Count path separators (handles both '/' and '\')
      file_depths <- lengths(gregexpr("[/\\\\]", matching_files))
      RNA_file <- matching_files[which.max(file_depths)]
    } else {
      RNA_file <- matching_files[1]
    }
    cnv_regions <- read.table(RNA_file, sep = "\t", quote = "", comment.char = "", 
                              fill = TRUE, header = TRUE)
    cat("Successfully loaded file from:\n", RNA_file, "\n")

  } else {
    stop("Please set `config$RNAdataSource = `1` or `2`")
  }
  cell_group_names <- cnv_regions$cell_group_name
  extracted_cell_group_names <- unique(sub("\\..*", "", cell_group_names))
  
  cat ("\ncnvTree Selection:",
       "  1. Enter a group number to select it.",
       "  2. Enter 0 to finish.",
       "  [Constraint]: Do not select the group used as the InferCNV reference.\n", 
       sep = "\n")
  
  menu_choices <- c(extracted_cell_group_names)
  exit_index <- length(menu_choices) 
  choice_index_0 <- menu(choices = menu_choices, title = "Select a group index:")
  
  if (choice_index_0 == 0) {
    stop("No group selected, stopping.") 
  } else {
    print(paste0("First Group Selected: ", extracted_cell_group_names[choice_index_0])) 
  }
  repeat {
    n <- length(extracted_cell_group_names)
    choice_index <- menu(
      choices = menu_choices, 
      title = "cnvTree: Choose more groups or enter 0 to finish.")
    if (choice_index == 0) {
      DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
      break 
    } else if (choice_index %in% choice_index_0) {
      if (!identical(sort(unique(choice_index_0)), seq_len(n))) {
        print("cnvTree: Group already selected. Select another group, or enter 0 to finish.") 
      } else {
        print("cnvTree: All groups selected. Select 0 to finish.") 
      }
    } else {
      choice_index_0 <- c(choice_index_0, choice_index)
      print(paste0("Group selected: ", extracted_cell_group_names[choice_index])) 
    }
  }
  selected_groups <- unique(extracted_cell_group_names[choice_index_0])
  message("Final selected groups: ")
  message(paste0(selected_groups, collapse = ", "))
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(selected_groups)
}


#' ProcessHmmList
#'
#' Load and Standardize AneuFinder HMM Objects. This function imports AneuFinder 
#' Hidden Markov Model (HMM) files from a specified directory. It ensures that 
#' the list names and internal object IDs match the file basenames rather than 
#' absolute system paths.
#'
#' @param input_dir_DNA A string specifying the directory containing the \code{.RData} 
#'   or \code{.rds} HMM files (typically the 'method-edivisive' output folder).
#'
#' @return A named list of AneuFinder HMM objects with standardized names and IDs.
#'
ProcessHmmList <- function(input_dir_DNA) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #---------------------- main code start below ---------------------
  
  files <- list.files(input_dir_DNA, full.names = TRUE)
  if (length(files) == 0) stop ("No hmm files found at: ", input_dir_DNA)
  hmms <- loadFromFiles(files)
  names(hmms) <- basename(names(hmms)) # Strip Path from Names
  
  # Sync Internal IDs
  for (i in seq_along(hmms)) {
    hmms[[i]]$ID <- names(hmms)[i]
  }
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(hmms)
}


#' cnvTree: Run Complete Single-Cell DNA and RNA CNV Tree Pipeline. 
#'
#' This function orchestrates the end-to-end CNV tree pipeline for single-cell data. 
#' It runs pipeline initialization, processes single-cell DNA clustering based on 
#' the selected data source (`"1"` or `"2"`), performs single-cell RNA clustering 
#' if DNA clusters are successfully detected, constructs and evaluates CNV trees 
#' (either RNA-only or combined DNA-RNA structures), exports the final clustering 
#' and visualization results, and saves the pipeline configuration log.
#'
#' @param cores An integer specifying the number of worker cores to use for parallel 
#'  processing.
#' @param DNAdatasource A character string or numeric value specifying the DNA data 
#'  source mode (`"1"` or `"2"`).
#' @param pqArm_file A character string specifying the file path for chromosome 
#'  arm (p/q arm) definitions.
#' @param output_dir A character string specifying the directory where the output 
#'  files will be saved.
#' @param input_dir_DNA A string specifying the directory containing the \code{.RData} 
#'  or \code{.rds} HMM files (typically the 'method-edivisive' output folder).
#' @param input_dir_RNA A character string specifying the main directory path 
#'  containing the InferCNV outputs.
#' @param cnv_ratio_DNA Numeric value between 0 and 1 indicating the minimum  
#'  coverage ratio threshold required to consider a CNV present 
#'  (\code{final_cnv = 1}).
#' @param cnv_ratio_RNA A numeric value defining the parameter filtering scale or  
#'  abundance ratio required to validate a superimposed CNV region call.
#' @param cellcutoff_DNA A numeric value defining the minimum number of cells required 
#'  for a cluster to be included.
#' @param cellcutoff_RNA An integer specifying the minimum number of single cells  
#'  required to keep a specific CNV pattern.
#' @param min_cell An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#' @param min_cell_subclone Minimum cell count required to define a valid subclone.
#' @param consecutive_region A numeric value defining the minimum length threshold 
#'  for filtering copy number variation (CNV) regions.
#' @param ZeroCol A logical value (\code{TRUE} or \code{FALSE}). If \code{TRUE},  
#'  removes columns where the sum of all cell clusters is zero.
#' @param cost_calculate A character string specifying the cost calculation method 
#'  for tree edge weights.
#' @param treeType A character string specifying the tree construction type, 
#'  either `"RNA"` or `"DNARNA"`.
#' @param difratio_chr A numeric value defining the threshold for acceptable 
#'  difference ratios across different chromosomes.
#' @param dif_ratio A numeric value specifying the tolerance threshold for 
#'  differences in copy number patterns between cells during clustering.
#' @param overlap_region A numeric value determines how much genomic space is 
#'  masked when identifying frequent CNV regions.
#'
#' @return No return value.
#'
#' @export
#' 
cnvTree <- function(cores, DNAdatasource, pqArm_file, output_dir, input_dir_DNA, 
                    input_dir_RNA, cnv_ratio_DNA, cnv_ratio_RNA, cellcutoff_DNA, 
                    cellcutoff_RNA, min_cell, min_cell_subclone, consecutive_region, 
                    overlap_region, ZeroCol, difratio_chr, dif_ratio, treeType,
                    cost_calculate) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.5_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #---------------------- main code start below ---------------------

  selected <- run_cnvTree_Pipeline(output_dir = output_dir,
                                   input_dir_RNA = input_dir_RNA)

  if (DNAdatasource == "1") {
    cnvTree_scDNAclustering_df(cores = cores,
                               input_dir_DNA = input_dir_DNA,
                               cnv_ratio_DNA = cnv_ratio_DNA,
                               output_dir = output_dir,
                               pqArm_file = pqArm_file, 
                               cellcutoff = cellcutoff_DNA,
                               min_cell = min_cell,
                               min_cell_subclone = min_cell_subclone,
                               consecutive_region = consecutive_region,
                               overlap_region = overlap_region,
                               difratio_chr = difratio_chr,
                               dif_ratio = dif_ratio)
  } else if (DNAdatasource == "2") {
    cnvTree_scDNAclustering(input_dir_DNA = input_dir_DNA,
                            output_dir = output_dir,
                            pqArm_file = pqArm_file, 
                            cellcutoff_DNA = cellcutoff_DNA,
                            cnv_ratio_DNA = cnv_ratio_DNA,
                            min_cell = min_cell,
                            min_cell_subclone = min_cell_subclone,
                            consecutive_region = consecutive_region,
                            overlap_region = overlap_region,
                            difratio_chr = difratio_chr,
                            dif_ratio = dif_ratio)
  }
  ckPoint_dna <- list.files(path = output_dir, 
                            pattern = "scDNAseq_DNAcluster.*\\.txt$", 
                            recursive = TRUE, 
                            full.names = TRUE)
  if (length(ckPoint_dna) == 0) {
    cat("No CNVs were identified.\n")
    stop("cnvTree scDNA_cluster Results are saved to: ", config$output_dir)
  } else {
    cnvTree_scRNAclustering(input_dir_RNA = input_dir_RNA,
                            selected_groups = selected,
                            output_dir = output_dir,
                            cellcutoff_DNA = cellcutoff_DNA,
                            cellcutoff_RNA = cellcutoff_RNA,
                            cnv_ratio_RNA = cnv_ratio_RNA,
                            ZeroCol = ZeroCol)
    ckPoint_rna <- list.files(path = output_dir, 
                              pattern = "scRNAseq_RNAcluster.*\\.txt$", 
                              recursive = TRUE, 
                              full.names = TRUE)
    
    ck <- read.table(ckPoint_rna[1], header = TRUE, sep = "")
    if (nrow(ck) == 1L && sum(as.numeric(ck[1, 1:(ncol(ck) - 2)])) == 0) {
      cat("cnvTree scRNA_Integrage Results are saved to: ", output_dir)
      stop("CNV patterns are too simplified to generate a tree.\n")
    } else {
      generate_Tree(output_dir = output_dir,
                    treeType = treeType,
                    cost_calculate = cost_calculate,
                    selected = selected[1])
    }
  }
  config_txt_path <- file.path(output_dir, "cnvTree_config.txt")
  config_lines <- capture.output(str(config))
  writeLines(c("=== CNVTree Configuration ===", "", config_lines), config_txt_path)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}
