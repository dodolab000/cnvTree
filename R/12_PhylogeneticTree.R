#' cnvTree_constructing
#'
#' This function constructs, evaluates, and iteratively refines single-cell CNV  
#' trees based on DNA and/or RNA cluster data. It performs initial data preprocessing,  
#' builds a weighted adjacency relationship graph, generates the shortest path  
#' tree structure, and evaluates the support from single cells. If the proportion  
#' of supported cells falls below a specified threshold, it iteratively removes 
#' unsupported or low-confidence CNVs until the quality criteria are met, and 
#' finally produces a tree visualization and a structured summary list.
#'
#' @param scRNA_cluster A data frame containing single-cell RNA cluster data.
#' @param scDNA_cluster A data frame containing single-cell DNA cluster data, or  
#'  `NULL` if only RNA data is available.
#' @param DefinedCNVs A vector or data frame defining the set of valid or reference 
#'  CNVs used for shortest path tree building.
#' @param cost_calculate A function or method specifying how to calculate edge  
#'  costs/weights for CNV relationships.
#' @param percentage Numeric value between 0 and 1 specifying the minimum required 
#'  proportion of cells supporting the tree structure (default is `0.5`).
#'
#' @return A `list` named `Tree_construction` containing:
#'  \item{Final_cluster}{The final processed and filtered cluster data frame.}
#'  \item{SupportedCells}{The number of cells supporting the final tree structure.}
#'  \item{TotalCells}{The total number of cells across all clusters.}
#'  \item{cnvTree_fig}{A ggplot or graphical object visualizing the constructed 
#'     CNV tree.}
#'
cnvTree_constructing <- function(scRNA_cluster, scDNA_cluster, DefinedCNVs, 
                                 cost_calculate, percentage)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.0_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  # Locked variable
  percentage <- getOption("percentage")
  
  #-------------------------- main code starts below ---------------------------
  
  # first round data preprocessing 
  if (is.null(scDNA_cluster) == TRUE) {
    Final_cluster <- cnvTree.RNApattern(clusterRNAData = scRNA_cluster) 
    Type = "RNA"
  } else {
    # merge DNA and RNA data
    Final_cluster <- cnvTree.DNARNApattern(clusterDNAData = scDNA_cluster, 
                                           clusterRNAData = scRNA_cluster)   
    Type = "DNARNA"
  }
  weighted_adj <- cnvTree.relation(inputFILE = Final_cluster, 
                                   cost_cal = cost_calculate)
  Tree <- cnvTree.shortestpath(inputFILE = Final_cluster, EditDis = weighted_adj, 
                               DefinedCNVs = DefinedCNVs)
  Total_cell <- sum(Final_cluster$Total_cellnum)
  evaluate <- cnvTree.evaluate(inputFILE = Final_cluster, tree_structure = Tree)

  # Lack of supported cells for tree structure
  Test = evaluate/Total_cell
  if (Test < percentage) {
    Data <- cnvTree.iteratively(inputFILE = Final_cluster, 
                                Type = Type, 
                                cost_calculate = cost_calculate, 
                                DefinedCNVs = DefinedCNVs)
    CNV_rmList <- Data$Iter$evaluate %>% as.data.frame(.) %>% t(.) %>%
                  as.data.frame(.) %>% setNames("cellnum") %>%
                  arrange(desc(cellnum)) %>% rownames(.)
    
    # Iteratively remove CNV for tree structure 
    count = 0
    while (Test < percentage) {
      count = count + 1
      Final_cluster <- cnvTree.removeCNV(inputFILE = Final_cluster, Type = Type,
                                         rmCNV = CNV_rmList[count])
      # new weighted_adj and Tree 
      weighted_adj <- cnvTree.relation(inputFILE = Final_cluster, cost_cal = cost_calculate)
      Tree <- cnvTree.shortestpath(inputFILE = Final_cluster, EditDis = weighted_adj, 
                                   DefinedCNVs = DefinedCNVs)
      evaluate <- cnvTree.evaluate(inputFILE = Final_cluster, tree_structure = Tree)
      Test = evaluate/Total_cell
    }
  } 
  # final output 
  cnvTree_fig <- cnvTree.plottree(inputFILE = Final_cluster, final_tree = Tree)
  Tree_construction <- list(Final_cluster = Final_cluster,
                            SupportedCells = evaluate, 
                            TotalCells = Total_cell, 
                            cnvTree_fig = cnvTree_fig)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)  
  
  return(Tree_construction)
}


#' cnvTree.DNARNApattern
#'
#' This function processes individual DNA and RNA cluster datasets containing copy  
#' number variation (CNV) profiles, merges them based on their unified CNV patterns,  
#' calculates aggregate statistics, classifies the data origin type (`"DNA_only"`,  
#' `"RNA_only"`, `"Both"`, or `"Neither"`), and prepares the data structure for
#' downstream CNV tree construction.
#'
#' @param clusterDNAData A data frame containing DNA cluster data, including 
#'  `DNA_cluster`, `DNA_Cellnum`, and specific CNV columns.
#' @param clusterRNAData A data frame containing RNA cluster data, including 
#'  `RNA_cluster`, `RNA_Cellnum`, `sum_cnv`, and corresponding CNV columns.
#'
#' @return A processed `data.frame` (`Final_cluster`) containing combined CNV 
#'  columns, merged cluster identifiers, data source classification (`Data`), 
#'  individual and total cell counts (`DNA_Cellnum`, `RNA_Cellnum`, `Total_cellnum`), 
#'  and tree ordering indices (`Tree_cluster`, `cellnum_rank`).
#'
#' @import tidyr
#' 
cnvTree.DNARNApattern <- function(clusterDNAData, clusterRNAData) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  # DNA data 
  cnvPattern <- clusterDNAData %>% 
                select(-c(DNA_cluster, DNA_Cellnum)) %>% 
                rowwise() %>%
                mutate(cnvPattern = paste(c_across(everything()), collapse = "_")) %>% 
                select(cnvPattern)
  DNA <- cbind(DNA_cluster = clusterDNAData$DNA_cluster, 
               DNA_Cellnum = clusterDNAData$DNA_Cellnum, 
               cnvPattern = cnvPattern)
  # RNA data
  cnvPattern <- clusterRNAData %>% 
                select(-c(RNA_Cellnum, sum_cnv, RNA_cluster)) %>% 
                rowwise() %>%
                mutate(cnvPattern = paste(c_across(everything()), collapse = "_")) %>% 
                select(cnvPattern)
  RNA <- cbind(RNA_Cellnum = clusterRNAData$RNA_Cellnum, 
               RNA_cluster = clusterRNAData$RNA_cluster, 
               cnvPattern = cnvPattern)
  
  # merge DNA, RNA data
  Final_cluster <- merge(DNA, RNA, by = "cnvPattern", all = TRUE)
  num_CNVs <- ncol(clusterDNAData) - 2
  Final_cluster <- separate(Final_cluster, cnvPattern, 
                            into = paste0("CNV", 1:num_CNVs), sep = "_")
  Final_cluster$sum_cnv <- Final_cluster %>%
                           select(-c(DNA_cluster, DNA_Cellnum, RNA_Cellnum, 
                                     RNA_cluster)) %>%
                           mutate_all(as.numeric) %>%
                           rowSums()
  # add normal cluster 
  if (min(Final_cluster$sum_cnv) != 0) {
    Final_cluster <- rbind(rep(0, ncol(clusterDNAData)), Final_cluster)
  }
  
  Final_cluster <- Final_cluster %>% 
    arrange(sum_cnv) %>% 
    mutate(
      row_cluster = c(1:n()),
      Data = case_when(
        is.na(DNA_cluster) == TRUE  & is.na(RNA_cluster) == FALSE ~ "RNA_only",
        is.na(DNA_cluster) == FALSE & is.na(RNA_cluster) == TRUE  ~ "DNA_only",
        is.na(DNA_cluster) == TRUE  & is.na(RNA_cluster) == TRUE  ~ "Neither",
        is.na(DNA_cluster) == FALSE & is.na(RNA_cluster) == FALSE ~ "Both"),
      DNA_Cellnum = ifelse(is.na(DNA_cluster) == TRUE, 0, DNA_Cellnum), 
      RNA_Cellnum = ifelse(is.na(RNA_Cellnum) == TRUE, 0, RNA_Cellnum),
      Tree_cluster = c(0:(n()-1)))
  Final_cluster$Total_cellnum <- rowSums(Final_cluster[, c("DNA_Cellnum", "RNA_Cellnum")])
  Final_cluster <- Final_cluster %>% 
                   arrange(Total_cellnum) %>% 
                   mutate(cellnum_rank = 1:n()) %>% 
                   arrange(Tree_cluster)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)  
  
  return(Final_cluster)
} 


#' cnvTree.RNApattern
#'
#' This function processes single-cell RNA cluster datasets containing copy number 
#' variation (CNV) profiles. It ensures the inclusion of a normal (zero-CNV) cluster 
#' if absent, ranks clusters by RNA cell numbers and CNV sums, assigns row clusters, 
#' classifies the data origin type (`"RNA_only"` or `"Neither"`), sets total cell 
#' numbers, and indexes the tree cluster positions.
#'
#' @param clusterRNAData A data frame containing RNA cluster data, including  
#'  `RNA_cluster`, `RNA_Cellnum`, `sum_cnv`, and corresponding CNV columns.
#'
#' @return A processed `data.frame` (`Final_cluster`) containing the adjusted  
#'  cluster data with additional columns for cell number ranking (`cellnum_rank`), 
#'  row clustering order (`row_cluster`),  data origin classification (`Data`), 
#'  total cell numbers (`Total_cellnum`), and tree cluster indices  (`Tree_cluster`).
#'
cnvTree.RNApattern <- function(clusterRNAData)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.2_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  if (min(clusterRNAData$sum_cnv) != 0) {
    clusterRNAData <- rbind(rep(0, ncol(clusterRNAData)), clusterRNAData)
    clusterRNAData$RNA_cluster[1] <- NA
  } 
  Final_cluster <- clusterRNAData %>% 
                   arrange(RNA_Cellnum) %>% 
                   mutate(cellnum_rank = 1:n()) %>% 
                   arrange(sum_cnv) %>% 
                   mutate(row_cluster = 1:n(),
                          Data = ifelse(is.na(RNA_cluster) == TRUE, 
                                        "Neither", "RNA_only"),
                          Total_cellnum = RNA_Cellnum, 
                          Tree_cluster = 0:(nrow(clusterRNAData) - 1))
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Final_cluster)
}


#' cnvTree.relation
#'
#' This function evaluates directed relationships and calculates edit distances 
#' between CNV clusters. It constructs a directed adjacency matrix checking 
#' directionality based on CNV profiles (where a parent cluster's CNV features 
#' are subsets of a child's), computes edge weights/costs using specified scaling 
#' ethods (e.g., cell counts or ranks), and combines them into a final weighted 
#' edit distance matrix for shortest-path tree construction.
#'
#' @param inputFILE A data frame containing cluster data, including `Total_cellnum`, 
#'  `cellnum_rank`, and CNV columns starting with `"CNV"`.
#' @param cost_cal A character string specifying the cost calculation method for 
#'  edge weights. Options include:
#'   \itemize{
#'     \item \code{"2nodeMul"}: Multiplies raw cell counts of both nodes. 
#'     Heavily penalizes merging large clusters.
#'     \item \code{"Rank_Mul"}: Multiplies the size ranks of both nodes. Reduces 
#'     the effect of size outliers while penalizing large-cluster merges.
#'     \item \code{"Rank_Add"}: Adds the size ranks of both nodes. Imposes a balanced, 
#'     linear penalty based on cluster size.
#'     \item \code{"Start_num"}: Uses only the cell count of the target node being 
#'     attached.
#'     \item \code{"NULL"}: Assigns a uniform cost of 1 to all merges, ignoring 
#'     cluster sizes.
#'   }
#'
#' @return A numeric matrix (`Final_EditDis`) representing the directed, weighted 
#'  edit distance between clusters, where unreachable paths are set to `Inf`.
#'
cnvTree.relation <- function(inputFILE, cost_cal)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.3_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  cnv_pattern <- inputFILE %>% select(starts_with("CNV"))
  
  # 確認cluster間是否相連並具有方向性: x-axis(from), y-axis(to)
  C_num <- nrow(cnv_pattern)
  DirGraph <- matrix(ncol = C_num, nrow = C_num)
  DirGraph[1, ] <- c(0, rep(1, C_num-1)) # cluster0對其他cluster 都可以相連
  for (x in 2:C_num) {
    for (y in 1:C_num) {
      if (x == y) {
        DirGraph[x, y] <- 0
      } else {
        from = which(cnv_pattern[x, ] == 1)
        to = which(cnv_pattern[y, ] == 1)
        DirGraph[x, y] <- ifelse(all(from %in% to) == TRUE, 1, -1)
      }
    }
  }
  # 確認cluster間的edit distance 
  EditDis <- matrix(ncol = C_num, nrow = C_num)
  for (x in 1:C_num) {
    for (y in 1:C_num) {
      # cost, larger makes weight smaller
      cost <- case_when(
        # 兩cluster細胞數相乘
        cost_cal == "2nodeMul" ~ (inputFILE$Total_cellnum[x]+1)*(inputFILE$Total_cellnum[y]+1), 
        # 細胞數先rank, rank值相乘,
        cost_cal == "Rank_Mul" ~ (inputFILE$cellnum_rank[x])*(inputFILE$cellnum_rank[y]), 
        # 細胞數先rank, rank值相加, 
        cost_cal == "Rank_Add" ~ (inputFILE$cellnum_rank[x])+(inputFILE$cellnum_rank[y]), 
        # 以接上cluster的細胞數
        cost_cal == "Start_num" ~ (inputFILE$Total_cellnum[y]+1), 
        cost_cal == "NULL" ~ 1)  
        
        EditDis[x, y] <- length(which(cnv_pattern[x, ] != cnv_pattern[y, ]))/cost
        # EditDis[x, y] <- cost
    }
  }
  # 將有方向性的矩陣與edit distance matrix結合 找出有最少CNV變化的路徑
  Final_EditDis <- ifelse(DirGraph == -1, Inf, EditDis)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)  

  return(Final_EditDis)
}


#' cnvTree.shortestpath
#'
#' This function calculates the shortest path tree structure for CNV clusters using 
#' a directed From-To transition matrix based on edit distances. It extracts utilized 
#' CNV changes, maps them to chromosomal cytogenetic bands using reference data 
#' (`DefinedCNVs`), formats amplification and deletion statuses, and merges multiple 
#' CNV changes per edge into a single formatted string for tree visualization.
#'
#' @param inputFILE A data frame containing cluster data with CNV profile columns 
#'  starting with `"CNV"`.
#' @param EditDis A numeric matrix of directed edit distances between clusters.
#' @param DefinedCNVs A data frame mapping CNV regions (`CNV_region`) to cytogenetic  
#'  bands (`chr`, `first_band`, `last_band`, `CN`).
#'
#' @return A processed `data.frame` (`final_edge`) containing the selected edges 
#'  for the shortest path tree, along with mapped chromosomal cytogenetic band 
#'  information and formatted change labels.
#'
cnvTree.shortestpath <- function(inputFILE, EditDis, DefinedCNVs)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.4_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  CNmatrix <- inputFILE %>% select(starts_with("CNV"))
  FromTo <- cnvTree.FromToMatrix(inputFILE = inputFILE, EditDis = EditDis, 
                                 CNmatrix = CNmatrix)
  # 從最前面的edge開始選
  edge <- NULL
  selected_from = 1
  while (nrow(FromTo) > 0) {
   From_to <- FromTo %>% filter(From == selected_from)
     if (nrow(From_to) == 0) {
       From_to <- FromTo
       available_node <- unique(c(edge$From, edge$To))
       if (!From_to$From[1] %in% available_node) { # 如果剩下的最小值並不是從1開始，代表他不是從1開始作為root出發
         break
       }
     }
     selected_edge <- FromTo %>% filter(combine %in% c(From_to[1, 6]))
     selected_from <- selected_edge[1, 2]
     edge <- rbind(edge, selected_edge)
     usedCNV <- selected_edge %>% pull(CNVchange)
     del <- FromTo %>% 
            filter(CNVchange %in% usedCNV) %>% 
            pull(combine)
     FromTo <- FromTo %>% filter(!combine %in% del)
  }
  # edge change into chr banding level
  edge$CNV_region <- as.numeric(str_extract(edge$CNVchange, "\\d+$"))
  edge <- left_join(edge, DefinedCNVs, by = "CNV_region")
  edge$band_chr <- str_extract(edge$chr, "\\d+$")
  
  edge <- edge %>% 
          mutate(status = ifelse(CN == "amp", "+", "-"), 
                Band = paste0(status, " ", band_chr, first_band, "~", last_band))
  final_edge <- NULL
  for (i in unique(edge$combine)) {
   selected_edge <- edge %>% filter(combine %in% c(i))
   if (nrow(selected_edge) > 1) {
    CNVchange <- paste(selected_edge$CNVchange, collapse ="\n")
    Band <- paste(selected_edge$Band, collapse ="\n")
    selected_edge <- selected_edge[1, ]
    selected_edge$CNVchange <- CNVchange 
    selected_edge$Band <- Band
  }
  final_edge <- rbind(final_edge, selected_edge)
 }
 #---------------------------- main code ends here -----------------------------
  
 DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
 return(final_edge)
}


#' cnvTree.FromToMatrix
#'
#' This function iterates through all pairwise combinations of clusters in an edit 
#' distance matrix, identifies valid directed edges (excluding zeros and infinite 
#' weights), computes individual CNV changes using `cnvTree.cnvchange()`, and compiles 
#' them into a structured data frame containing source nodes, target nodes, weights, 
#' change descriptions, and the number of differing CNVs. The resulting matrix is 
#' sorted by source node, number of differences, and weight.
#'
#' @param inputFILE A data frame containing cluster data with CNV profile columns.
#' @param EditDis A numeric matrix of directed edit distances between clusters.
#' @param CNmatrix A matrix or data frame subset of CNV profiles corresponding to 
#'  the clusters.
#'
#' @return A processed `data.frame` (`FromTo`) with columns:
#'  \item{From}{The source cluster index.}
#'  \item{To}{The target cluster index.}
#'  \item{weight}{The calculated edge weight/distance.}
#'  \item{CNVchange}{The specific CNV differences between the clusters.}
#'  \item{dif_CNVs}{The total count of differing CNV features.}
#'  \item{combine}{A character string combining the source and target indices 
#'  (`"From_To"`).}
#'
cnvTree.FromToMatrix <- function(inputFILE, EditDis, CNmatrix)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.4.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  C_num <- nrow(EditDis)
  FromTo <- NULL
  for (i in 1:C_num) {
    for (j in 1:C_num) {
      if (!EditDis[i, j] %in% c(0, Inf)) {
        CNVchange <- cnvTree.cnvchange(inputFILE = inputFILE, cluster1 = i, 
                                       cluster2 = j)
        r <- cbind(i, j, EditDis[i, j], CNVchange, length(CNVchange))
        FromTo <- rbind(FromTo, r)
      }
    }
  }
  FromTo <- as.data.frame(FromTo) %>% 
            setNames(c("From", "To", "weight", "CNVchange", "dif_CNVs")) %>% 
            mutate_at(vars(From, To, weight), as.numeric) %>% 
            mutate(combine = paste(From, To, sep = "_")) %>% 
            arrange(From, dif_CNVs, weight)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(FromTo)
}


#' cnvTree.cnvchange
#'
#' This function extracts the copy number variation (CNV) profiles for two specified 
#' row clusters from an input dataset and identifies which CNV features differ 
#' between them.
#'
#' @param inputFILE A data frame containing cluster data with `row_cluster` and  
#'  CNV columns starting with `"CNV"`.
#' @param cluster1 Integer or numeric value representing the first row cluster 
#'  index.
#' @param cluster2 Integer or numeric value representing the second row cluster 
#'  index.
#'
#' @return A character vector (`dif`) containing the names of the CNV columns where 
#'   the profiles of `cluster1` and `cluster2` differ.
#'   
cnvTree.cnvchange <- function(inputFILE, cluster1, cluster2)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.4.1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  pattern <- inputFILE %>% 
             filter(row_cluster %in% c(cluster1, cluster2)) %>%
             select(starts_with("CNV"))
  Name <- names(pattern)
  dif <- Name[pattern[1, ] != pattern[2, ]]
  # dif <- paste0(dif, collapse ="\n")

  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(dif)
}


#' cnvTree.plottree
#'
#' This function takes the final tree structure and input cluster data to construct 
#' an igraph object and plot a visual CNV tree layout using `ggraph`. Nodes are 
#' sized based on cell count ranges (with customized legend labels for percentages) 
#' and colored according to their data origin classification (`"DNA_only"`, 
#' `"RNA_only"`, `"Both"`, or `"Neither"`). Edges are labeled with chromosomal 
#' band changes.
#'
#' @param inputFILE A data frame containing cluster data, including `Total_cellnum`, 
#'  `Tree_cluster`, and `Data` classification columns.
#' @param final_tree A data frame representing the final tree edges, including 
#'  `From`, `To`, `weight`, `CNVchange`, and `Band`.
#'
#' @return A `ggraph`/`ggplot` object (`tree_fig`) visualizing the CNV tree.
#'
cnvTree.plottree <- function(inputFILE, final_tree)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.5_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  C_num <- nrow(inputFILE)
  
  # get edges
  edges <- final_tree %>% 
           arrange(From) %>% 
           mutate(From = From - 1, To = To - 1) %>% 
           select(From, To, weight,  CNVchange, Band)
  selected_node <- unique(c(edges$From, edges$To))
  
  # get vertices info 
  maxcell <- sum(inputFILE$Total_cellnum)
  Range <- round(c(maxcell * 0.25, maxcell * 0.5, maxcell * 0.75), digits = 0)
  Groups <- inputFILE %>% arrange(Tree_cluster) %>% select(Data)
  vertices <- data.frame(
    name = c(0:(C_num-1)),
    Data  = Groups,  #reload
    value = case_when(
      inputFILE$Total_cellnum <= Range[1] ~ 5,
      inputFILE$Total_cellnum > Range[1] & inputFILE$Total_cellnum <= Range[2] ~ 10, 
      inputFILE$Total_cellnum > Range[2] & inputFILE$Total_cellnum <= Range[3] ~ 15,
      inputFILE$Total_cellnum > Range[3] ~ 20))
  vertices <- vertices %>% filter(name %in% selected_node)
  final_graph <- igraph::graph_from_data_frame(edges, vertices = vertices)
  
  # change figure legend names
  Label <- c(paste("≤", Range[1], "(≤25%)"),
             paste(Range[1]+1, "~", Range[2], "(25%~50%)"),
             paste(Range[2]+1, "~", Range[3], "(50%~75%)"),
             paste(">", Range[3], "(>75%)"))
  
  # plot tree
  tree_fig <- ggraph(final_graph, layout = 'dendrogram') + 
    geom_edge_link(aes(label = Band), width = 1, color = "grey") +  # 折線或曲線+
    geom_node_point(aes(size=value, color = Data), alpha = 1) +
    geom_node_text(aes(label = name), vjust = 0.5, size = 5, colour = "white", 
                   fontface = "bold") +  # 加入文字標籤在每個點上方
    scale_size_continuous(labels = Label, 
                          range = c(10, 20),
                          limits = c(5, 20)) +  # 設定最小和最大的點大小
    scale_colour_manual(values = c("DNA_only" = "#36558F", "RNA_only" = "#C179B9", 
                                   "Both" = "#7C67A4", "Neither" = "#AEB7B3"), 
                        breaks = c("Neither", "DNA_only", "RNA_only", "Both")) +
    theme(legend.position = "right",
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 14),
          panel.background = element_blank()) +
    guides(size = guide_legend(title = "# cells"),
           colour = guide_legend(title = "Data", override.aes = list(size = 14)))
  plot(tree_fig)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(tree_fig)
}


#' cnvTree.evaluate
#'
#' This function calculates the total number of cells supported by a constructed 
#' CNV tree structure by summing the cell counts (`Total_cellnum`) of all unique 
#' nodes (clusters) involved in the tree edges.
#'
#' @param inputFILE A data frame containing cluster data, including `row_cluster` 
#'  and `Total_cellnum`.
#' @param tree_structure A data frame representing the tree structure, containing 
#'  at least `From` and `To` node columns corresponding to row cluster indices.
#'
#' @return A numeric value (`Totalcell`) representing the total sum of cells 
#'  supported by the tree nodes.
#'
cnvTree.evaluate <- function(inputFILE, tree_structure)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.6_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  selected_node <- unique(c(tree_structure$From, tree_structure$To))
  # number of cells
  Totalcell <- inputFILE %>% 
               filter(row_cluster %in% selected_node) %>% 
               select(Total_cellnum) %>% 
               sum(.)
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Totalcell)
}


#' cnvTree.iteratively
#'
#' Performs an iterative workflow across all detected CNV columns in the input 
#' dataset. For each CNV set, it removes the specified CNV, calculates weighted 
#' adjacency relationships, constructs shortest-path phylogenetic trees, generates 
#' tree plots, and evaluates the resulting tree structures.
#'
#' @param inputFILE A data frame containing cell profiles and columns starting with 
#'   \code{"CNV"} representing distinct copy number variations.
#' @param Type Character string specifying the type or method for CNV removal 
#'  (passed to \code{cnvTree.removeCNV}).
#' @param cost_calculate Character string specifying the cost calculation method 
#'  for the relationship matrix (e.g., \code{"2nodeMul"}).
#' \itemize{
#'     \item \code{"2nodeMul"}: Multiplies raw cell counts of both nodes. 
#'     Heavily penalizes merging large clusters.
#'     \item \code{"Rank_Mul"}: Multiplies the size ranks of both nodes. Reduces 
#'     the effect of size outliers while penalizing large-cluster merges.
#'     \item \code{"Rank_Add"}: Adds the size ranks of both nodes. Imposes a balanced, 
#'     linear penalty based on cluster size.
#'     \item \code{"Start_num"}: Uses only the cell count of the target node being 
#'     attached.
#'     \item \code{"NULL"}: Assigns a uniform cost of 1 to all merges, ignoring 
#'     cluster sizes.
#'   }
#' @param DefinedCNVs A data frame or vector of defined CNV regions used for 
#'  shortest-path tree construction.
#'
#' @return A nested list (\code{Data_iter$Iter}) containing the following elements 
#'  for each CNV iteration:
#'   \describe{
#'     \item{\code{RM}}{Modified dataset objects after CNV removal.}
#'     \item{\code{weighted_adj}}{Weighted adjacency matrices.}
#'     \item{\code{Tree}}{Inferred shortest-path phylogenetic trees.}
#'     \item{\code{cnvTree_fig}}{Generated tree plot visualization objects.}
#'     \item{\code{evaluate}}{Evaluation metrics for the tree structures.}
#'   }
#'
cnvTree.iteratively <- function(inputFILE, Type, cost_calculate, DefinedCNVs)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.7_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  # new cluster data
  cnv_vec <- colnames(inputFILE)
  CNV <- cnv_vec[grep("^CNV", cnv_vec)]
  
  Data_iter <- NULL
  for (k in 1:length(CNV)) {
    Data_iter$Iter$RM[[CNV[k]]] <- cnvTree.removeCNV(inputFILE = inputFILE, 
                                                     rmCNV = CNV[k], Type = Type)
  
    # new weighted_adj and Tree 
    Final_cluster <- as.data.frame(Data_iter$Iter$RM[[CNV[k]]])

    weighted_adj <- cnvTree.relation(inputFILE = Final_cluster, cost_cal = cost_calculate)
    Data_iter$Iter$weighted_adj[[CNV[k]]] <- weighted_adj
  
    Tree <- cnvTree.shortestpath(inputFILE = Final_cluster, EditDis = weighted_adj, 
                                 DefinedCNVs = DefinedCNVs)
    Data_iter$Iter$Tree[[CNV[k]]] <- Tree
  
    cnvTree_fig <- cnvTree.plottree(inputFILE = Final_cluster, final_tree = Tree)
    Data_iter$Iter$cnvTree_fig[[CNV[k]]] <- cnvTree_fig
  
    evaluate <- cnvTree.evaluate(inputFILE = Final_cluster, tree_structure = Tree)
    Data_iter$Iter$evaluate[[CNV[k]]] <- evaluate
  }
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)
  
  return(Data_iter)
}


#' cnvTree.removeCNV
#'
#' Removes a target CNV column from the input dataset, checks for duplicate 
#' patterns, combines redundant clusters, and recalculates cell number counts, 
#' sums, and cluster rankings.
#'
#' @param inputFILE A data frame containing CNV columns (starting with \code{"CNV"}) 
#'  and associated metadata.
#' @param rmCNV Character string specifying the name of the CNV column to remove.
#' @param Type Character string indicating the data type (e.g., \code{"RNA"}) 
#'  used to handle cell count updates appropriately.
#'
#' @return A modified data frame containing the updated dataset with the specified 
#'  CNV removed, recalculated cell numbers, and updated ranking metrics.
#'   
cnvTree.removeCNV <-function(inputFILE, rmCNV, Type)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.7.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  CNmatrix = inputFILE %>% select(starts_with("CNV")) %>% mutate_all(as.numeric)
  rmCNmatrix = CNmatrix %>% select(-rmCNV)
  
  # preprocessed data 
  rm1 <- apply(rmCNmatrix, 1, function(x) paste(x, collapse = "")) 
  new_sumcnv <- apply(rmCNmatrix, 1, function(x) sum(x))
  Input <- inputFILE %>% 
           select(-c(rmCNV, "cellnum_rank", "sum_cnv","row_cluster", 
                     "Total_cellnum", "Tree_cluster")) %>% 
           mutate(checkpattern = rm1, sum_cnv = new_sumcnv)
  check <- as.data.frame(table(rm1)) %>% filter(Freq != 1) %>% pull(rm1) %>%
           as.character()
         
  # deal with duplicated data
  if (length(check) != 0) {
    # need to combine clusters from inputFILE file
    for (k in 1:length(check)) {
      Input <- cnvTree.cnvCombine(original = Input, selected_pattern = check[k], type = Type)
    }
    # Total_cellnum
    Input$RNA_Cellnum <- as.numeric(Input$RNA_Cellnum)
    if (Type == "RNA") {
      Input$Total_cellnum <- Input$RNA_Cellnum
    } else {
      Input$DNA_Cellnum <- as.numeric(Input$DNA_Cellnum)
      Input$Total_cellnum <- rowSums(Input[, c("DNA_Cellnum", "RNA_Cellnum")])
    }
    # sum_cnv
    Input$sum_cnv <- Input %>% 
                     select(starts_with("CNV")) %>% 
                     mutate_all(as.numeric) %>% 
                     apply(., 1, function(x) sum(x)) 
    # cellnum_rank, row_cluster
    Input <- Input %>% 
             arrange(Total_cellnum) %>% 
             mutate(cellnum_rank = 1:n()) %>% 
             arrange(sum_cnv) %>% 
             mutate(row_cluster = 1:n(), Tree_cluster = (1:n()) - 1)
  } else {
    Input <- inputFILE %>% select(-rmCNV)
  }
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(Input)
}


#' cnvTree.cnvCombine
#'
#' Combine Duplicate Clusters During Iterative CNV Removal. Filters rows matching 
#' a specific check pattern, aggregates cell numbers and cluster names based on 
#' the analysis type (\code{"RNA"} or \code{"DNARNA"}), resolves data sources, 
#' and merges the combined cluster back into the dataset.
#'
#' @param original A data frame containing the current dataset and metadata.
#' @param selected_pattern Character string specifying the check pattern whose 
#'  matching rows need to be combined.
#' @param type Character string indicating the data type context, either 
#'   \code{"RNA"} or \code{"DNARNA"}. Default is \code{c("RNA", "DNARNA")}.
#'
#' @return A modified data frame where rows matching the selected pattern 
#'   have been merged into a single consolidated cluster row.
#' 
cnvTree.cnvCombine <- function(original, selected_pattern, type = c("RNA", "DNARNA"))
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.7.1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  ss <- original %>% 
        filter(checkpattern %in% selected_pattern)
  s <- ss %>% select(-grep("^CNV", names(.))) 

  if (type == "RNA") {
    newcluster <- NULL
    Names <- colnames(s)
    for (i in 1:length(Names)) {
      if (Names[i] == "RNA_Cellnum") {
        newcluster = cbind(newcluster, sum(as.numeric(s$RNA_Cellnum)))
      } else if(Names[i] == "RNA_cluster") {
        newcluster = cbind(newcluster, paste0(s$RNA_cluster, collapse = "_"))
      } else if(Names[i] == "Data") {
        pattern = unique(s$Data)
        newcluster = cbind(newcluster, cnvTree.check_vector(vec = pattern))
      } else {
        newcluster = cbind(newcluster, s[1, i])
      }
    }
  } else if (type == "DNARNA") {
    newcluster <- NULL
    Names <- colnames(s)
    for (i in 1:length(Names)) {
      if(Names[i] == "RNA_Cellnum"){
        newcluster = cbind(newcluster, sum(as.numeric(s$RNA_Cellnum)))
      } else if (Names[i] == "DNA_Cellnum"){
        newcluster = cbind(newcluster, sum(as.numeric(s$DNA_Cellnum)))
      } else if (Names[i] == "RNA_cluster"){
        newcluster = cbind(newcluster, paste0(s$RNA_cluster, collapse = "_"))
      } else if (Names[i] == "DNA_cluster"){
        newcluster = cbind(newcluster, paste0(s$DNA_cluster, collapse = "_"))
      } else if(Names[i] == "Data"){
        pattern = unique(s$Data)
        newcluster = cbind(newcluster, cnvTree.check_vector(vec = pattern))
      } else {
        newcluster = cbind(newcluster, s[1, i])
      }
    }
  }
  cnv_pattern <- ss %>% select(starts_with("CNV"))
  newcluster <- as.data.frame(newcluster) %>% setNames(Names) 
  newcluster <- cbind(cnv_pattern[1, ], newcluster)
  original <- original %>% filter(!checkpattern == selected_pattern) %>% 
              rbind(newcluster)

  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg)

  return(original)
}


#' cnvTree.check_vector
#'
#' Check and Categorize Data Source Vector. Evaluates a vector of data source 
#' classifications (e.g., \code{"RNA_only"}, \code{"DNA_only"}, \code{"Both"}, or 
#' \code{"Neither"}) and consolidates them into a single unified category string.
#'
#' @param vec A character vector of data source labels to evaluate.
#'
#' @return A character string representing the consolidated data source category: 
#'   \code{"Both"}, \code{"RNA_only"}, \code{"DNA_only"}, or \code{"Neither"}.
#'   
cnvTree.check_vector <- function(vec) 
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.7.1.1.1_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  if (any(vec == "RNA_only" & !is.na(vec)) & any(vec == "DNA_only" & !is.na(vec))) { # RNA+DNA
    return("Both")
  } else if (any(vec == "Both" & !is.na(vec))) { # both出現
    return("Both")
  } else if (any(vec == "RNA_only" & !is.na(vec)) & any(vec == "Neither" & !is.na(vec))) { #
    return("RNA_only")
  } else if (any(vec == "DNA_only" & !is.na(vec)) & any(vec == "Neither" & !is.na(vec))) {
    return("DNA_only")
  } else if (all(vec == "RNA_only" & !is.na(vec))) {
    return("RNA_only")
  } else if (all(vec == "DNA_only" & !is.na(vec))) {
    return("DNA_only")
  } else {
    return("Neither")
  }
  #---------------------------- main code ends here ----------------------------
  
  cat(DebugMsg(fucStep, "end", cnvTree_msg = cnvTree_msg))
}

  
#' cnvTree_Output
#'
#' This function saves the final clustering results of a CNV tree analysis to an 
#' output directory, customizes the tree plot with evaluation metrics and titles, 
#' and exports the resulting plot as a high-resolution TIFF image.
#'
#' @param Summary A list containing tree construction summary results, including 
#'  `Final_cluster`, `SupportedCells`, `TotalCells`, and `cnvTree_fig`.
#' @param output_dir A character string specifying the destination directory path 
#'  for the output files.
#' @param name A character string prefix used for naming output files.
#' @param plot_height Numeric value specifying the height of the saved plot 
#'  (default is `10`).
#' @param plot_width Numeric value specifying the width of the saved plot 
#'  (default is `15`).
#'
#' @return No return value, called for writing output files and saving plots.
#'
#' @importFrom ggplot2 ggsave
#' 
cnvTree_Output <- function(Summary, output_dir, name, plot_height = 10, plot_width = 15)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.8_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  message("=== Output cnvTree Tree structure results ===")
  timestamp <- format(Sys.time(), "%m%d")
  
  # 1. Final clustering output 
  FILE = paste0("/", timestamp, "_", name, "_cnvTree.Final_grouping")
  writeOutput(data = Summary$Final_cluster, filename = FILE, path = output_dir)

  message("Output /cnvTree.scDNAseq_grouping.txt is done.")
  
  # 2. plot output
  Title = paste0("cnvTree structure\n (Supported cells/Total cells: ", 
                 Summary$SupportedCells, "/",Summary$TotalCells, ")")
  Summary$cnvTree_fig <- Summary$cnvTree_fig +
                         ggtitle(Title) +
                         theme(plot.title = element_text(hjust = 0.5, size = 20, 
                                                         face = "bold"),
                               plot.margin = margin(1.5,2,1.5,2, "cm"))
  
  ggsave(Summary$cnvTree_fig, 
         filename = paste0("/", timestamp, "_", name, "_cnvTree_fig.tiff"),
         path = output_dir,
         height = plot_height, 
         width = plot_width,
         limitsize = FALSE)
  
  message("Output /cnvTree.scDNAseq_cnvTree_fig.tiff is done.")
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
}


#' generate_Tree
#'
#' Generate CNV Tree Structure. Reads scDNA, scRNA, and defined CNV data files 
#' from the output directory, constructs a copy number variation (CNV) tree based 
#' on the designated tree type, and writes the results to disk.
#'
#' @param output_dir A character string specifying the directory where the output 
#'  files will be saved.
#' @param treeType Character string. The type of tree to construct. Options are:
#'   \itemize{
#'     \item \code{"RNA"}: Constructs the tree using single-cell RNA cluster data 
#'     only.
#'     \item \code{"DNARNA"}: Constructs the tree using both single-cell DNA and 
#'     RNA cluster data.
#'   }
#' @param cost_calculate A function or method specifying how to calculate edge  
#'  costs/weights for CNV relationships.
#'  \itemize{
#'     \item \code{"2nodeMul"}: Multiplies raw cell counts of both nodes. 
#'     Heavily penalizes merging large clusters.
#'     \item \code{"Rank_Mul"}: Multiplies the size ranks of both nodes. Reduces 
#'     the effect of size outliers while penalizing large-cluster merges.
#'     \item \code{"Rank_Add"}: Adds the size ranks of both nodes. Imposes a balanced, 
#'     linear penalty based on cluster size.
#'     \item \code{"Start_num"}: Uses only the cell count of the target node being 
#'     attached.
#'     \item \code{"NULL"}: Assigns a uniform cost of 1 to all merges, ignoring 
#'     cluster sizes.
#'   }
#' @param selected Character. Sample type identifiers.
#'
#' @examples
#' \dontrun{
#'   generate_Tree(output_dir = config$output_dir,
#'                 treeType = config$treeType,
#'                 cost_calculate = config$cost_calculate,
#'                 selected = "sample_name")
#' }
#' 
#' @export
#' 
generate_Tree <- function(output_dir, treeType, cost_calculate, selected)
{
  cnvTree_v_num <- getOption("cnvTree_v_num")
  cnvTree_msg   <- getOption("cnvTree_msg")
  fucStep <- paste0(" 12.9_cnvTree_", cnvTree_v_num)
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
  
  #-------------------------- main code starts below ---------------------------
  
  scDNA_file <- list.files(path = output_dir, 
                           pattern = "DNAcluster.*\\.txt$", 
                           recursive = TRUE, 
                           full.names = TRUE)
  scRNA_file <- list.files(path = output_dir,
                           pattern = "RNAcluster.*\\.txt$", 
                           recursive = TRUE, 
                           full.names = TRUE)
  DefinedCNV_filie <- list.files(path = output_dir, 
                                 pattern = "DefinedCNVregion.*\\.txt$", 
                                 recursive = TRUE, 
                                 full.names = TRUE)
  
  ck <- read.table(scRNA_file[1], header = TRUE, sep = "")
  if (nrow(ck) == 1L && sum(as.numeric(ck[1, 1:(ncol(ck) - 2)])) == 0) {
    stop("CNV patterns are too simplified to generate a tree.\n")
  }
  scDNA_cluster <- read.table(scDNA_file, header = TRUE)
  scRNA_cluster <- read.table(scRNA_file, header = TRUE) %>%
                   mutate(sum_cnv = rowSums(pick(starts_with("CNV")), na.rm = TRUE))
  DefinedCNVs <- read.table(DefinedCNV_filie, header = TRUE)

  if (treeType == "RNA") {
    cnvTree <- cnvTree_constructing(scRNA_cluster = scRNA_cluster, 
                                    scDNA_cluster = NULL, 
                                    DefinedCNVs = DefinedCNVs,
                                    cost_calculate = cost_calculate)
  } else if (treeType == "DNARNA") {
    cnvTree <- cnvTree_constructing(scRNA_cluster = scRNA_cluster, 
                                    scDNA_cluster = scDNA_cluster, 
                                    DefinedCNVs = DefinedCNVs,
                                    cost_calculate = cost_calculate)
  }
  cnvTree_Output(Summary = cnvTree, 
                 name = paste0(selected[1], "_", treeType, "_", cost_calculate),
                 output_dir = output_dir)
  
  #---------------------------- main code ends here ----------------------------
  
  DebugMsg(fucStep, "start", cnvTree_msg = cnvTree_msg)
}

