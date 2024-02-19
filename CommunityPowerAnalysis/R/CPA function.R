
#' Process KEGG Results in a Directory
#'
#' This function processes KEGG (Kyoto Encyclopedia of Genes and Genomes) 
#' annotation results located in multiple subdirectories. It reads files, 
#' filters out lines starting with '#', and computes the count of each KEGG Orthology (KO).
#'
#'
#' @param inputDir A string specifying the path to the input directory containing subdirectories of KEGG result files.
#' @param filePattern An optional regular expression pattern to match file names within each subdirectory. If NULL, all files are considered.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A data frame with KOs as rows and each file's KEGG Orthology count as columns.
#' 
#' @examples
#' # Assuming 'result_dir' is your directory containing KEGG result files
#' result_dt <- KeggResultProcess(inputDir = result_dir, filePattern = '.fa.faa.annotation.table.txt')
#'
#' @importFrom data.table fread
#' @importFrom dplyr filter group_by summarize mutate anti_join
#' @importFrom tidyr unite
#' @export
KeggResultProcess <- function(inputDir, filePattern = NULL, ...) {
  # Validate input directory
  if (!dir.exists(inputDir)) {
    stop("The specified directory does not exist.")
  }
  
  # Get list of subdirectories
  subDirs <- list.dirs(path = inputDir, full.names = TRUE, recursive = FALSE)
  if (length(subDirs) == 0) {
    stop("No subdirectories found in the directory.")
  }
  
  result_dt <- data.frame('KO'=character())
  
  # Iterate over subdirectories
  for (subDir in subDirs) {
    # Get list of files in the subdirectory
    files <- list.files(path = subDir, pattern = filePattern,  full.names = TRUE)
    if (length(files) == 0) {
      next # Skip subdirectories with no relevant files
    }
    
    # Process each file in the subdirectory
    for (file in files) {
      # Read the file, skipping lines that start with '#'
      data <- fread(file, header = TRUE, sep = "\t")
      data <- data %>% filter(gene != '#')
      community_name <- sub(".*/(.*)", "\\1", subDir)
      sample_name <- sub(paste0("^.*/(.*?)", filePattern, "$"), "\\1",file)
      column_name <- paste0(community_name, "_", sample_name)
      # Update the list of all KOs
      data <- data  %>% group_by(KO) %>% summarize(!!column_name := n())
      result_dt <- result_dt %>% full_join(data, by='KO')
    }
  }
  # Replace NA with 0 in all columns except 'KO'
  result_dt <- result_dt %>%
    mutate(across(-KO, ~ replace(., is.na(.), 0)))
  return(result_dt)
}


#' Calculate Unique KEGG Counts for Each Sample in Communities
#'
#' This function calculates the number of unique KEGG Orthology (KO) identifiers 
#' for each sample within each community. It identifies KOs that are present 
#' in a sample but not in others within the same community.
#'
#' 
#' 
#' @param keggData A data frame or data table with KOs as rows and samples as columns.
#'
#' @return A data frame with each row representing a community-sample pair and 
#'         the corresponding count of unique KOs in that sample.
#'
#' @examples
#' # Assuming 'keggData' is your data frame or data table
#' uniqueKeggCPA <- UniqueKeggCPA(keggData)
#'
#' @export
UniqueKeggCPA <- function(keggData) {
  keggDT <- as.data.table(keggData)
  
  # Extract community names
  communityNames <- unique(sapply(strsplit(names(keggDT)[-1], "_", fixed = TRUE), `[`, 1))
  
  # Initialize an empty DataFrame to store the results
  result_df <- data.frame(Community = character(), Sample = character(), UniqueKeggCount = numeric())
  
  for (community in communityNames) {
    # Select samples for the current community
    sampleCols <- grep(paste0("^", community, "_"), names(keggDT), value = TRUE)
    communityData <- keggDT[, c("KO", sampleCols), with = FALSE]
    communityData <- communityData[rowSums(communityData[, -1, with = FALSE]) > 0, ]
    
    for (sample in sampleCols) {
      # Read the target sample KOs
      targetKOs <- communityData[, .SD, .SDcols = c('KO',sample)]
      targetKOs <- targetKOs[rowSums(targetKOs[,-1, with = FALSE]) > 0, ]
      
      # Combine KOs from other samples
      otherKOs <- communityData[, .SD, .SDcols = c('KO', setdiff(sampleCols, sample))]
      otherKOs <- otherKOs[rowSums(otherKOs[,-1, with = FALSE]) > 0, ]
    
      # Find unique KOs in the target sample
      uniqueKOs <- anti_join(targetKOs, otherKOs, by = "KO")
      
      # Count unique KOs
      uniqueKeggCount <- nrow(uniqueKOs)
      
      # Extract the sample name
      sampleName <- sub("^.*?_", "", sample)
      
      # Add to results
      result_df <- rbind(result_df, data.frame(Community = community, Sample = sampleName, UniqueKeggCount = uniqueKeggCount))
    }
  }
  
  return(result_df)
}



#' Calculate Community Power Values for Each KO in Each Sample
#'
#' This function calculates the Community Power (CP) values for each KEGG 
#' Orthology (KO) in each sample within a community. CP value is the proportion 
#' of a specific KO in a sample relative to the total count of that KO in the community.
#'
#' 
#' @param keggData A data frame or data table with KOs as rows and samples as columns.
#'
#' @return A data frame with each row representing a community-sample-KO triplet and 
#'         the corresponding CP value, along with the CP Score which is the sum of 
#'         CP values for each sample.
#'
#' @examples
#' # Assuming 'keggData' is your data frame or data table
#' CPAResult <- CommunityPowerAnalysis(keggData)
#' 
#' @export
CommunityPowerAnalysis <- function(keggData) {
  keggDT <- as.data.table(keggData)
  
  # Extract community names
  communityNames <- unique(sapply(strsplit(names(keggDT)[-1], "_", fixed = TRUE), `[`, 1))
  
  # Initialize list to store results
  resultList <- list()
  
  for (community in communityNames) {
    # Select samples for the current community
    communityCols <- grep(paste0("^", community, "_"), names(keggDT), value = TRUE)
    communityData <- keggDT[, c("KO", communityCols), with = FALSE]
    communityData <- communityData[rowSums(communityData[, -1, with = FALSE]) > 0, ]
    
    # Compute CP values
    communityDT <- melt(communityData, id.vars = "KO", variable.name = "Sample", value.name = "Count")
    communityDT[, TotalKO := sum(Count, na.rm = TRUE), by = KO]
    communityDT[, CP_value := Count / TotalKO]
    communityDT[, Sample := sub("^.*?_", "", Sample)]
    communityDT[, Community := community]
    
    # Append to the result list
    resultList[[community]] <- communityDT[, .(Community, Sample, KO, CP_value)]
  }
  
  # Combine results
  resultDT <- rbindlist(resultList)
  
  # Calculate CP Score for each sample
  resultDT[, CP_Score := sum(CP_value), by = Sample]
  
  return(as.data.frame(resultDT))
}




