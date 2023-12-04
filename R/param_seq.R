#' pick_rep
#'
#' @description The function to investigate the number of repeats in each motif
#' @param seq A character of the sequence of the repeat region
pick_rep <- function(seq){
  seq_sep <- strsplit(seq, " ")[[1]]
  reps <- integer(0)
  for(i in 1:length(seq_sep)){
    motif <- seq_sep[i]
    pos_bracket <- regexpr("]", motif)[1]
    if(pos_bracket != -1){
      reps <- c(reps, as.numeric(strsplit(motif, ']')[[1]][2]))
    }
  }
  return(reps)
}

#' calc_weighted_lus
#'
#' @description The function to calculate the weighted LUS of an allele
#' @param lus_one_al A vector of the LUS values of an allele
#' @param count_one_al A vector of the number of counts of an allele
calc_weighted_lus <- function(lus_one_al, count_one_al){
  return(sum(lus_one_al * count_one_al / sum(count_one_al)))
}

#' analyze_seq
#'
#' @description The function to analyze sequence data
#' @param dt_seq_init The initial data.table of sequence data
#' @param markers The vector of markers to estimate parameters
analyze_seq <- function(dt_seq_init, par_markers){
  lus_all_mk <- list()
  reps_all_mk <- list()
  seq_al_all_mk <- list()
  seq_count_all_mk <- list()

  for(i in 1:length(par_markers)){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq_init[Marker == par_markers[i]]
    alleles <- dt_seq_one_mk[, Allele]
    al_names <- sort(unique(alleles))
    counts <- dt_seq_one_mk[, Count]
    sequences <- dt_seq_one_mk[, RepeatRegion]

    # Investigate the number of repeats in each motif
    reps_one_mk <- lapply(sequences, pick_rep)

    # Calculate weighted LUS values
    lus_all_reps <- sapply(reps_one_mk, max)
    lus_one_mk <- rep(0, length(al_names))
    for(j in 1:length(al_names)){
      index_al <- which(alleles == al_names[j])
      lus_one_mk[j] <- calc_weighted_lus(lus_all_reps[index_al], counts[index_al])
    }
    names(lus_one_mk) <- al_names

    # Assign data of one marker
    reps_all_mk[[i]] <- reps_one_mk
    lus_all_mk[[i]] <- lus_one_mk
    seq_al_all_mk[[i]] <- alleles
    seq_count_all_mk[[i]] <- counts
  }

  names(lus_all_mk) <- par_markers
  names(reps_all_mk) <- par_markers
  names(seq_al_all_mk) <- par_markers
  names(seq_count_all_mk) <- par_markers

  return(list(lus_all_mk, reps_all_mk, seq_al_all_mk, seq_count_all_mk))
}
