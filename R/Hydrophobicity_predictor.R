#' Hydrophobic region predictor
#'
#' Predict hydrophobic region of a protein by selecting a protein sequence file from the prompt
#' @param length is the desired size of the rolling average
#' @return A list of several elements including the hydrophobicity plot, the hydrophobicity dataframe, the input sequence and a secondary dataframe used in downstream analysis
#' @examples 
#' Protein_name_hydrophobicity_output <- Hydrophobicity_predictor(20)
#' Protein_name_hydrophobicity_plot <- Hydrophobicity_predictor(20)[[1]];
#' Protein_name_hydrophobicity_df <- Hydrophobicity_predictor(20)[[2]];
#' Protein_name_input_sequence <- Hydrophobicity_predictor(20)[[3]];
#' @export

Hydrophobicity_predictor <- function(length = 20){
  require(tidyverse) 
  #AA_values <- read.csv("AA_values")
  
  # Prompt to select protein sequence file
  # Extracts protein name from fasta title and protein sequence from file
  filename <- file.choose()
  Protein_name <- read_lines(filename)
  Protein_name <- Protein_name[1] %>% strsplit(" ") 
  Protein_name <- Protein_name[[1]][2:3] %>% paste(collapse = " ")
  inputseq <- read_lines(filename, skip = 1)
  sequence <- unlist(strsplit(inputseq, ""))
  
  # Charge identifier function
  charge_df <- function(filename){
    inputseq <- read_lines(filename, skip = 1)
    sequence <- unlist(strsplit(inputseq, ""))
    
    #  For loop to create a vector of charge identifiers based on if the amino acid is positivly charged or not
    charge_vector <- vector()
    AA_seq <- seq(from = 1, to = length(sequence))
    residue_number_vector <- vector()
    
    for (residue in AA_seq) {
      if(sequence[residue] == "K" | sequence[residue] == "R" | sequence[residue] == "H"){
        charge_vector[residue] <- "positive"
      } else {
        charge_vector[residue] <- "uncharged"
      }
      residue_number_vector[residue] <- residue
    } 
    charge_vector
  }
  charge_vector <- charge_df(filename)
  
  # Function for converting AA into hydrophobicity values
  AA_converter <- function(filename){
    
    inputseq <- read_lines(filename, skip = 1)
    protseq <- unlist(strsplit(inputseq, ""))
    
    value_vector <- seq(from = 1, to = 20)  
    hydro_vector <- vector()
    seqAA_vector <- seq(from = 1, to = length(protseq))
    
    for (seqAA in seqAA_vector) {
      for (AA in value_vector) {
        
        if (protseq[seqAA] == AA_values[AA,1]) {
          value <- AA_values[AA,2]
        }
      }
      hydro_vector[seqAA] <- value
    }
    hydro_vector
  }
  seq_hydro_values <- AA_converter(filename)
  
  # Rolling average function
  rolling_average <- function(sequence, length){
    
    posvector <- seq(from=1, to=length(sequence)-(length-1))
    avgvector <- vector()
    avg_size <- seq(from = 1, to = (length-1))
    
    for (position in posvector) {
      plusone_vector <- vector()
      for (plusone in avg_size) {
        plusone_vector[plusone] <- sequence[position + plusone]
        
      }
      plusone_vector <- append(plusone_vector, sequence[position])
      average1 <- mean(plusone_vector)
      avgvector[position] <- average1
      
    }
    avgvector
  }
  
  # Creating a data frame using the rolling average function and combining it with the charge identifier and residue numbers
  df <- data_frame("Residue" = seq(from = 1, to = length(seq_hydro_values)-(length-1)),
                   "Hydrophobicity" = rolling_average(seq_hydro_values, length), 
                   "Charge" = charge_vector[1:(length(charge_vector)-(length-1))])
  
  # Plotting the data
  g <- ggplot(df) +
    geom_line(mapping = aes(Residue, Hydrophobicity)) +
    theme_bw() +
    scale_x_continuous(breaks = seq(from = 0, to = length(seq_hydro_values), by = 20)) +
    scale_y_continuous(breaks = seq(from = round(min(df$Hydrophobicity), digits = 0), to = round(max(df$Hydrophobicity), digits = 0), by = 0.5)) +
    geom_hline(yintercept=0) +
    labs(title = paste(Protein_name, "Hydrophobicity Plot with sampling width", length)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    geom_point(mapping = aes(Residue, Hydrophobicity, color = Charge, alpha = Charge)) +
    scale_alpha_discrete(limits = rev, guide = "none")
  
  # Make a second df with the addition of the unaveraged hydrovalues for later use
  hydro_df <- data_frame("Residue" = seq(from = 1, to = length(seq_hydro_values)), "Hydrophobicity" = seq_hydro_values)
  df_union <- union(df[,c(1,2)], hydro_df[(1+nrow(hydro_df)-(nrow(hydro_df)-nrow(df))):nrow(hydro_df),1:2])
  
  return(list(g, df, inputseq, df_union))
  
}

