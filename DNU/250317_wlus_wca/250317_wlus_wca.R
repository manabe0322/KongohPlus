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

#' calc_weighted_ca
#'
#' @description The function to calculate the weighted corrected allele number
#' @param min_us Minimum uninterrupted stretch
#' @param reps_one_al Repeat motifs of one allele
#' @param count_one_al Counts for sequence data in one allele
calc_weighted_ca <- function(min_us, reps_one_al, count_one_al){
  # Calculate corrected allele numbers in each sequence
  n_seq_al <- length(reps_one_al)
  ca_one_al <- rep(0, n_seq_al)
  for(i in 1:n_seq_al){
    tmp <- reps_one_al[[i]] - min_us
    ca_one_al[i] <- sum(tmp[tmp > 0])
  }

  # Calculate the weighted corrected allele number
  weighted_ca <- sum(ca_one_al * count_one_al / sum(count_one_al))
  return(weighted_ca)
}

#' correct_al_rep
#'
#' @description The function to correct allele repeat numbers
#' @param dt_seq The data.table of sequence data
#' @param markers_auto A vector of the autosomal markers
#' @param min_us_auto A vector of the minimum uninterrupted stretches of the autosomal markers
correct_al_rep <- function(dt_seq, markers_auto, min_us_auto){
  setindex(dt_seq, Marker)
  n_mk <- length(markers_auto)
  wlus_auto <- vector(mode = "list", length = n_mk)
  wca_auto <- vector(mode = "list", length = n_mk)

  for(i in 1:n_mk){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq[Marker == markers_auto[i]]
    alleles <- as.numeric(dt_seq_one_mk$Allele)
    al_names <- sort(unique(alleles))
    n_al <- length(al_names)
    counts <- dt_seq_one_mk$Count
    sequences <- dt_seq_one_mk$RepeatRegion
    min_us <- min_us_auto[i]

    # Investigate the number of repeats in each motif
    reps_one_mk <- lapply(sequences, pick_rep)

    # Calculate weighted LUS values
    lus_all_reps <- sapply(reps_one_mk, max)
    lus_one_mk <- rep(0, n_al)
    for(j in 1:n_al){
      index_al <- which(alleles == al_names[j])
      lus_one_mk[j] <- calc_weighted_lus(lus_all_reps[index_al], counts[index_al])
    }
    names(lus_one_mk) <- al_names

    # Calculate weighted CA values
    ca_one_mk <- rep(0, n_al)
    for(j in 1:n_al){
      index_al <- which(alleles == al_names[j])
      ca_one_mk[j] <- calc_weighted_ca(min_us, reps_one_mk[index_al], counts[index_al])
    }
    names(ca_one_mk) <- al_names

    # Assign data of one marker
    wlus_auto[[i]] <- lus_one_mk
    wca_auto[[i]] <- ca_one_mk
  }

  names(wlus_auto) <- markers_auto
  names(wca_auto) <- markers_auto

  return(list("wlus_auto" = wlus_auto,
              "wca_auto" = wca_auto))
}
