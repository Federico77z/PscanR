#' @keywords internal
#' @importFrom TFBSTools PFMatrixList
#' 
.ps_checks <- function(x, pfms, type)
{
 # .ps_required_packages()
  
  if(type == 4 && !is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PSMatrixList class")
  
  if(type == 4)
  {
    nabg <- vapply(pfms, ps_bg_avg, FUN.VALUE = numeric(length = 1L))
    nabg2 <- vapply(pfms, ps_bg_std_dev, FUN.VALUE = numeric(length = 1L))
    
    nabg <- is.na(nabg)
    nabg2 <- is.na(nabg2)
    nabg <- nabg | nabg2
    
    if(any(nabg))
      warning(paste("\nNo Pscan background for", name(pfms)[nabg], 
                    ID(pfms)[nabg]))
  }
  
  if(!is(pfms, "PFMatrixList") && !is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PFMatrixList or PSMatrixList class")  
  
  if(is.character(x) && type == 2)
  {
    if(file.access(x, mode = 4) != 0)
      stop(paste("Cannot access file path:",x))
    else
    {
      first_line <- readLines(x, n = 1)
        if(first_line != "[SHORT TFBS MATRIX]")
          stop(paste(x, " does not look like a Pscan .short_matrix file"))
    }
  }
  else if(!is(x, "DNAStringSet") && (type == 1 || type == 4))
    stop("x is not an object of DNAStringSet class")
  
  else if(is.data.frame(x) && type == 3) {
    req_cols <- c("BG_SIZE", "BG_MEAN", "BG_STDEV")
    
    if(!all(req_cols %in% colnames(x)))
      stop("x does not contain required columns",
           "\"BG_SIZE\", \"BG_MEAN\", \"BG_STDEV\"")
    
    if(!all(is.numeric(x$BG_SIZE), is.numeric(x$BG_MEAN), 
            is.numeric(x$BG_STDEV)))
      stop("Required columns of x must be of numeric type")
  }
}

#' @keywords internal
.ps_checks2 <- function(pfms, file = NULL, ...)
{
  #.ps_required_packages()

  if(!is(pfms, "PSMatrixList"))
    stop("pfms is not an object of PSMatrixList class")  
  
  if(!is.null(file) && !is.character(file)){
    stop('file must be a character string indicating the file path')
    
    file_dir <- dirname(file)
    
    if (!file.exists(file)) {
      if (file.access(file_dir, mode = 2) != 0) {
        stop(paste("Can't write to directory:", file_dir))
        }
      } else {
        if (file.access(file, mode = 2) != 0)
          stop(paste("Can't write to existing file:", file))
      }
  }
}

#' @keywords internal
.clean_sequence <- function(x){
  
  seq_widths <- Biostrings::width(x)
  ref_width <- max(Biostrings::width(x))
  
  diff_length_seq <- x[seq_widths != ref_width]
  
  if(length(diff_length_seq) != 0){
    warning(paste(length(diff_length_seq), 'sequences found with length 
                  different from the reference. Removing the following 
                  sequences:', 
                  paste(names(diff_length_seq), collapse = ", ")))
  }
  x <- x[seq_widths == ref_width]
  
  freqs <- Biostrings::alphabetFrequency(x, as.prob = FALSE)
  n_proportions <- freqs[, "N"] / ref_width
  rem_names <- names(x[n_proportions > 0.5])
  
  if(length(rem_names)> 0){
    warning(paste('Found', length(rem_names), 'sequences with more than 50% of N. 
                  Removing the following sequences:', 
                  paste(rem_names, collapse = ", ")))
  }
  x <- x[n_proportions <= 0.5]
  return(x)
}

#' @keywords internal
.mapping_unique_names <- function(x,pfms){

  original_names <- names(x)
  
  # popolazione vettore con tutti i nomi di sequenze presenti nell'org di studio
  all_sequences_ID <- setNames(character(length(x)), original_names)
  
  # eliminazione doppioni 
  unique_x <- BiocGenerics::unique(x)
  unique_names <- names(unique_x)  
  
  # Trasformazione in vettore di caratteri nominato per fare il confronto tra i 
  # due vettori con la funzione match
  x_char <- as.character(x)
  unique_x_char <- as.character(unique_x)
  
  # Confronti: quando 2 seq sono uguali, viene memorizzato un indice numerico che verrà 
  # usato per estrarre il nome della sequenza che è stata mantenuta da unique(). 
  # questo nome sarà il VALORE assegnato al vettore con il NOME della sequenza originale. 
  all_sequences_ID <- setNames(unique_names[match(x_char, unique_x_char)], original_names)
  
  # applico la funzione per rimuovere le seq con %N > 50%. nel vettore, in corrispondenza 
  # del nome della sequenza originale, al posto del nome della sequenza mantenuta 
  # da unique() ci sarà un NA --> nella funzione principale questo NA serve come 
  # flag per rimuovere il nome di questa sequenza, nel caso sia stata passata in 
  # input, ed emettere un warning.
  x <- .clean_sequence(x)
  removed_sequences <- setdiff(original_names, names(x))
  
  all_sequences_ID[removed_sequences] <- NA
  
  # Assegno i valori solo all'array della prima PSM per questioni di spazio 
  # --> non credo vada bene, sarebbe meglio assegnarlo (una sola volta) a tutta la 
  # PSMatrixList (campo metadati?)
  pfms@transcriptIDLegend <- all_sequences_ID
  
  return(pfms)
}

.download_background <- function(file, destfile = file){
  
  URL <- paste0('https://raw.githubusercontent.com/dianabetelli/PscanR_backgrounds/refs/heads/main/BG_files/',file)
  download.file(URL, destfile, mode = "wb")
  
  return(destfile)
}

#.ps_required_packages <- function()
#{
#  require("Biostrings")
#  require("TFBSTools")
#  require("BiocParallel")
#  require("BSDA")
#}
