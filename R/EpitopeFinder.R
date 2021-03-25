#' Epitope Mapper
#'
#' Maps epitopes contained in an IEDB export file to the query protein, by taking the output of the Hydrophobicity_predictor function as input, and selecting the desired IEDB export .csv file from the prompt.
#' @param Hydrophobicity_predictor_output The output list generated from the Hydrophobicity_predictor_output function
#' @param rm_hydrophobic_regions Logical argument to specify if the epitopes found inside hydrophobic regions should be included in the resulting dataframe
#' @return Dataframe consisting of epitopes, and their location, found in the query protein sequence as well as metadata about the original epitope 
#' @examples 
#' Protein_name_epitope_dataframe <- Epitope_finder(Protein_name_hydrophobicity_output, rm_hydrophobic_regions = FALSE)
#' @export

Epitope_finder <- function(Hydrophobicity_predictor_output, rm_hydrophobic_regions = FALSE){
  require(stringr)
  inputseq <- Hydrophobicity_predictor_output[[3]]
  inputseq <- paste(inputseq, collapse = "")
  
  # Bug occured as package was built from source
  # This fixes the union_df which in the packaged version was turned into a list of four elements
  Hydrophobicity_predictor_output[[4]][[1]] <- append(Hydrophobicity_predictor_output[[4]][[1]], Hydrophobicity_predictor_output[[4]][[3]])
  Hydrophobicity_predictor_output[[4]][[2]] <- append(Hydrophobicity_predictor_output[[4]][[2]], Hydrophobicity_predictor_output[[4]][[4]])
  Hydrophobicity_predictor_output[[4]][[3]] <- NULL
  Hydrophobicity_predictor_output[[4]][[3]] <- NULL
  
  df_union <- data_frame("Residue" = Hydrophobicity_predictor_output[[4]][[1]], "Hydrophobicity" = Hydrophobicity_predictor_output[[4]][[2]])
  
  # Prompt to select epitope file, and reading it in as a .csv as well as creating a vector of the epitope sequences
  epitope_file <- choose.files()
  epitope_df <- read_csv(epitope_file, skip = 1)
  epitope <- unlist(epitope_df[1:nrow(epitope_df),3])
  epitope <- unlist(strsplit(paste(epitope, collapse = " "), split = " "))
  
  # Creates a vector of the epitopes mapped to the query protein sequence
  epi_seqs <- str_extract_all(inputseq, epitope)
  epi_seqs <- unlist(epi_seqs[lapply(epi_seqs,length)>0])
  
  # Stops the function if there are no epitope matches 
  if(is.null(epi_seqs)){
    stop("There are no epitope matches")
  }
  
  # Function to find the start and end of epitope matches
  s_e_func <- function(inputseq, epitope){
    
    # using StringR to locate the matching epitopes
    locations <- str_locate_all(inputseq, epitope)
    locations <- locations[lapply(locations,length)>0]
    start_end <- vector()
    
    for (posish in seq(from = 1, to = length(locations))) {
      start_end <- append(start_end, locations[[posish]][1])
      start_end <- append(start_end, locations[[posish]][2])
    }
    # Vector containing the start and end residue of all mapped epitopes
    start_end
  }
  start_end <- s_e_func(inputseq, epitope)
  
  # Function to determine if the matches are located in a hydrophillic or soluble region
  is_hydro <- function(startendvector){  
    truefalsevector <- vector()
    
    for (posish2 in seq(from = 1, to = length(startendvector), by = 2)){
      
      if(df_union[startendvector[posish2], 2] < 0 & df_union[startendvector[posish2+1], 2] < 0){
        truefalsevector <- append(truefalsevector, print("Soluable_region"))
      } else {
        truefalsevector <- append(truefalsevector, print("Hydrophobic_region"))
      }
    }
    truefalsevector
  }
  is_hydro_vector <- is_hydro(start_end)
  
  # Creates a data frame of the results
  epitope_finder_df <- data_frame("Description" = epi_seqs, "Epitope_structural_location" = is_hydro_vector, "Start (residue)" = start_end[seq(from = 1, to = length(start_end), by = 2)], "End (residue)" = start_end[seq(from = 2, to = length(start_end), by = 2)])
  
  # Creates a data frame of epitope metadata
  epitope_df <- epitope_df[,c(1,3,11,12)]
  
  # Using dplyr to join the two dataframes based on the epitope sequence
  epitope_finder_df <- left_join(epitope_finder_df, epitope_df, by = "Description")
  
  # Omitts the epitopes located in hydrophobic region if this option is chosen
  if(rm_hydrophobic_regions == TRUE){
    epitope_finder_df <- filter(epitope_finder_df, Epitope_structural_location != "TM_region")
  }
  return(epitope_finder_df)
}
