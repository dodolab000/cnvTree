#' run_cnvTree_Pipeline
#' 
#' Run the cnvTree Setup Pipeline; this initializes the output environment, sets 
#' up folder structures for RNA data, and extracts the available cell groups/types
#' for selection.
#' 
#' @param input_dir_RNA Path to the input single-cell RNA-seq data directory
#'  (e.g., inferCNV output folder).
#' @param output_dir Path to the directory where results and output folders will
#'  be created and saved.
#'  
#' @return A character vector or list of identified cell groups/types extracted 
#'  by \code{select_groups()} for the next steps of the pipeline.
#' 
#' @export 
#' 
run_cnvTree_Pipeline <- function(input_dir_RNA, output_dir) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.0_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  RNAdataSource <- getOption("RNAdataSource")
  
  #-------------------------- main code starts below ---------------------------
  
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    message("Created output folder: ", output_dir)
  }
  selected_groups <- select_groups(input_dir_RNA = input_dir_RNA, 
                                   RNAdataSource = RNAdataSource)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(selected_groups)
}


#' select_groups
#'
#' This function scans a single-cell RNA-seq data directory, extracts unique cell
#' group names, and provides an interactive menu for the user to select one or 
#' more groups for downstream analysis. User Flow:
#'   1. Select an initial group. If '0' is pressed, the function execution stops.
#'   2. A repeat loop allows for multiple additional selections.
#'
#' @param input_dir_RNA Path to the input single-cell RNA-seq data directory
#'  (e.g., inferCNV output folder).
#' @param RNAdataSource: The type of RNA data source being used:
#'  '1' Tabular data file / matrix structure.
#'  '2' InferCNV output folder structure.
#'
#' @return A character vector containing the unique names of the selected cell 
#'  groups.
#' 
select_groups <- function(input_dir_RNA, RNAdataSource)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.0.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  if (RNAdataSource == "1") {
    print("datatype: dataframe")
    ext <- tools::file_ext(input_dir_RNA)
    if (ext == "rds") {
      cnv_regions <- readRDS(input_dir_RNA)
    } else if (ext == "txt") {
      cnv_regions <- read.table(input_dir_RNA, sep = " ", header = TRUE) 
    } else {
      stop("Unsupported format. Please provide a .rds, .txt, or .tsv file.")
    }
    
  } else if (RNAdataSource == "2") { 
    print("datatype: InferCNV files")
    all_dirs <- list.dirs(path = input_dir_RNA, recursive = TRUE, full.names = TRUE)
    dir_depths <- lengths(gregexpr("/", all_dirs))
    max_depth <- max(dir_depths)
    folders_pth <- all_dirs[dir_depths == (max_depth-1)][1]
    File <- list.files(folders_pth, full.names = TRUE)
    File <- File[grepl("cnv_regions", File) & grepl("HMM_CNV_predictions", File)]
    cnv_regions <- read.table(File, sep = "\t", quote = "", comment.char = "", 
                              fill = TRUE, header = TRUE)
    
  } else {
    stop("Please set RNAdataSource = 1 or 2 (Note: RNAdataSource = 1 has not been built yet)")
  }
  cell_group_names <- cnv_regions$cell_group_name
  extracted_cell_group_names <- unique(sub("\\..*", "", cell_group_names))
  
  cat (
    "\ncnvTree Selection:",
    "  1. Enter a group number to select it.",
    "  2. Enter 0 to finish.",
    "  [Constraint]: Do not select the group used as the InferCNV reference.\n", 
    sep = "\n"
  )
  
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
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(selected_groups)
}


#' NEW_changeFormat
#'
#' Convert copy number data to GRanges Format. This function reads a `.rds` file 
#' or `.txt` file containing copy number variation (CNV) data and converts it into 
#' a list of `GRanges` objects, where each element corresponds to a single cell.
#'
#' @param input_dir_DNA A `.rds` file or `.txt` file containing a data frame with 
#'  the following required columns:
#'  - `cellID`: Unique identifier for each cell.
#'  - `seqnames`: Chromosome or sequence name.
#'  - `start`: Start position of the segment.
#'  - `end`: End position of the segment.
#'  - `copy.number`: Copy number value for the segment.
#'   
#' @return A named list where each element represents a cell, containing its 
#'  corresponding genomic segments as a `GRanges` object.
#'  
#' @import data.table
#'
NEW_changeFormat <- function(input_dir_DNA)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.1a_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  sexchromosome <- getOption("sexchromosome")
  
  #-------------------------- main code starts below ---------------------------
  
  message("=== Step 1: NEW_changeFormat ===")
  
  input_dir_DNA <- config$input_dir_DNA
  
  if (endsWith(input_dir_DNA, ".rds")) {
    Bin_CN <- data.table::as.data.table(readRDS(input_dir_DNA))
  } else if (endsWith(input_dir_DNA, ".txt")) {
    Bin_CN <- data.table::fread(input_dir_DNA, sep = " ")
  } else {
    stop("Please input either .rds or .txt format as input.")
  }
  
  if (!sexchromosome) {
    Bin_CN <- Bin_CN[!grepl("chrX|chrY|chrM", seqnames)]
  }
  
  # Adjust 0-based BED coordinates to 1-based
  if (any(Bin_CN$start == 0, na.rm = TRUE)) {
    Bin_CN[, start := start + 1L]
  }
  
  # Chromosome Level Sorting
  has_chr_prefix <- any(grepl("(?i)^chr", Bin_CN$seqnames))
  raw_seqs <- unique(as.character(Bin_CN$seqnames))
  clean_chroms <- unique(sub("(?i)^chr", "", trimws(raw_seqs)))
  
  is_numeric <- !is.na(suppressWarnings(as.numeric(clean_chroms)))
  sorted_nums <- sort(as.numeric(clean_chroms[is_numeric]))
  sorted_sex <- sort(clean_chroms[!is_numeric])
  clean_order <- c(as.character(sorted_nums), sorted_sex)
  
  correct_level_order <- if (has_chr_prefix) paste0("chr", clean_order) else clean_order
  
  # Select and format final input columns
  structured_Bin_CN <- Bin_CN[,.(cellID = as.character(cellID),
                                 seqnames = factor(seqnames, levels = correct_level_order),
                                 start = as.integer(start),
                                 end = as.integer(end),
                                 copy.number = as.integer(copy.number))]
  
  # Split by cellID and create base template GRanges
  structured_Bin_CN_list <- split(structured_Bin_CN, by = "cellID", keep.by = TRUE)
  
  # Base template GRanges
  template_data <- structured_Bin_CN_list[[1]]
  template_gr <- GenomicRanges::makeGRangesFromDataFrame(template_data, keep.extra.columns = TRUE)
  cols <- c("seqnames", "start", "end", "copy.number")
  
  NewFormat <- lapply(seq_along(structured_Bin_CN_list), function(i) {
    cell_data <- structured_Bin_CN_list[[i]]
    cellid <- names(structured_Bin_CN_list)[i]
    
    # Update Bins
    Bins <- template_gr
    S4Vectors::mcols(Bins)$copy.number <- cell_data$copy.number
    
    chr_list <- unique(cell_data$seqnames)
    # Process chromosomes for the current cell
    breakpoint_list <- lapply(chr_list, function(chr) {
      sub_data <- cell_data[seqnames == chr]
      lag_cp <- shift(sub_data$copy.number, type = "lag")
      is_bp <- is.na(lag_cp) | (sub_data$copy.number != lag_cp)
      
      return(sub_data[is_bp, ..cols])
    })
    breakpoints <- rbindlist(breakpoint_list)
    
    # Construct Breakpoint GRanges
    Breakpoints <- if (nrow(breakpoints) == 0) {
      NULL
    } else {
      GenomicRanges::GRanges(seqnames = breakpoints$seqnames,
                             ranges = IRanges::IRanges(start = breakpoints$start,
                                                       end = breakpoints$end),
                             copy.number = breakpoints$copy.number)
    }
    list(ID = names(structured_Bin_CN_list)[i], bins = Bins, breakpoints = Breakpoints)
  })
  names(NewFormat) <- names(structured_Bin_CN_list)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(NewFormat)
}


#' ProcessHmmList
#'
#' Load and Standardize AneuFinder HMM Objects. This function imports AneuFinder 
#' Hidden Markov Model (HMM) files from a specified directory. It ensures that 
#' the list names and internal object IDs match the file basenames rather than 
#' absolute system paths.
#'
#' @param input_dir_DNA A string specifying the directory containing the 
#'  \code{.RData} or \code{.rds} HMM files (typically the 'method-edivisive' 
#'  output folder).
#'
#' @return A named list of AneuFinder HMM objects with standardized names and IDs.
#'
ProcessHmmList <- function(input_dir_DNA) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.1b_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  message("=== Step 1: ProcessHmmList ===")
  
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


#' NEW_CN_seq
#'
#' This function extracts copy number variations from a list of `GRanges` objects 
#' and organizes them into an integer matrix. The matrix contains selected cells
#' as columns, with genomic regions (fixed bins) as rows.
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
  fucStep <- paste0(" 1.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  c <- lapply(Template, function(i) {
    obj <- input[[i]]
    if (!is.null(obj) && "bins" %in% names(obj) && 
        "copy.number" %in% names(S4Vectors::mcols(obj[["bins"]]))) 
    {
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


#' cnvTree: Run Complete Single-Cell DNA and RNA CNV Tree Pipeline. 
#'
#' This function orchestrates the end-to-end CNV tree pipeline for single-cell data. 
#' It runs pipeline initialization, processes single-cell DNA clustering based on 
#' the selected data source (`"1"` or `"2"`), performs single-cell RNA clustering 
#' if DNA clusters are successfully detected, constructs and evaluates CNV trees 
#' (either RNA-only or combined DNA-RNA structures), exports the final clustering 
#' and visualization results, and saves the pipeline configuration log.
#'
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
#' @param min_cell_RNA An integer specifying the minimum number of single cells  
#'  required to keep a specific CNV pattern.
#' @param min_cell_DNA An integer specifying the minimum cell count required for a 
#'  cluster to be included in the output.
#' @param min_cell_subclone Minimum cell count required to define a valid subclone.
#' @param ZeroCol A logical value (\code{TRUE} or \code{FALSE}). If \code{TRUE},  
#'  removes columns where the sum of all cell clusters is zero.
#' @param cost_calculate A character string specifying the cost calculation method 
#'  for tree edge weights.
#' @param treeType A character string specifying the tree construction type, 
#'  either `"RNA"` or `"DNARNA"`.
#' @param GMM_mode  Logical. If \code{TRUE}, the pipeline incorporates the Gaussian 
#'  Mixture Model (GMM) calculation.
#'
#' @return No return value.
#'
#' @examples
#' \dontrun{
#' config_path <- system.file("extdata", "cnvTree_config.yaml", package = "cnvTree")
#' file.edit(config_path)
#' 
#' config <- yaml::read_yaml(config_path)
#' cnvtree_args <- intersect(names(config), names(formals(cnvTree)))
#' do.call(cnvTree, config[cnvtree_args])
#' }
#' 
#' @export
#' 
cnvTree <- function(DNAdatasource, pqArm_file, output_dir, input_dir_DNA, 
                    input_dir_RNA, cnv_ratio_DNA, cnv_ratio_RNA, min_cell_DNA, 
                    min_cell_subclone, GMM_mode, ZeroCol, treeType, cost_calculate,
                    min_cell_RNA) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 1.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #--------------------------- main code starts below --------------------------
  
  selected <- run_cnvTree_Pipeline(output_dir = output_dir,
                                   input_dir_RNA = input_dir_RNA)
  
  if (DNAdatasource == "1") {
    cnvTree_scDNAclustering_df(input_dir_DNA = input_dir_DNA,
                               cnv_ratio_DNA = cnv_ratio_DNA,
                               output_dir = output_dir,
                               pqArm_file = pqArm_file,
                               min_cell_DNA = min_cell_DNA,
                               min_cell_subclone = min_cell_subclone,
                               GMM_mode = GMM_mode)
  } else if (DNAdatasource == "2") {
    cnvTree_scDNAclustering(input_dir_DNA = input_dir_DNA,
                            cnv_ratio_DNA = cnv_ratio_DNA,
                            output_dir = output_dir,
                            pqArm_file = pqArm_file,
                            min_cell_DNA = min_cell_DNA,
                            min_cell_subclone = min_cell_subclone,
                            GMM_mode = GMM_mode)
  }
  ckPoint_dna <- list.files(path = output_dir, 
                            pattern = "scDNAseq_DNAcluster.*\\.txt$", 
                            recursive = TRUE, 
                            full.names = TRUE)
  if (length(ckPoint_dna) == 0) {
    cat("No CNVs were identified.\n")
    stop("cnvTree scDNA_cluster Results are saved to: ", config$output_dir)
  } else {
    cnvTree_scRNAclustering_s3(input_dir_RNA = input_dir_RNA,
                               selected_groups = selected,
                               output_dir = output_dir,
                               min_cell_RNA = min_cell_RNA,
                               cnv_ratio_RNA = cnv_ratio_RNA,
                               ZeroCol = ZeroCol)
    ckPoint_rna <- list.files(path = output_dir, 
                              pattern = "scRNAseq_RNAcluster.*\\.txt$", 
                              recursive = TRUE, 
                              full.names = TRUE)
    
    ck <- read.table(ckPoint_rna[1], header = TRUE, sep = "")
    if (nrow(ck) == 1L && sum(as.numeric(ck[1, 1:(ncol(ck) - 2)])) == 0) {
      cat("cnvTree scRNA_Integrage Results are saved to: ", config$output_dir)
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
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
}

