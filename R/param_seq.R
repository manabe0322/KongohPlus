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
#' @param markers_auto The vector of autosomal markers
analyze_seq <- function(dt_seq_init, markers_auto){
  lus_auto <- list()
  reps_auto <- list()
  seq_al_auto <- list()
  seq_count_auto <- list()

  for(i in 1:length(markers_auto)){
    # Extract data of one marker
    dt_seq_one_mk <- dt_seq_init[Marker == markers_auto[i]]
    alleles <- as.numeric(dt_seq_one_mk[, Allele])
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
    lus_auto[[i]] <- lus_one_mk
    reps_auto[[i]] <- reps_one_mk
    seq_al_auto[[i]] <- alleles
    seq_count_auto[[i]] <- counts
  }

  names(lus_auto) <- markers_auto
  names(reps_auto) <- markers_auto
  names(seq_al_auto) <- markers_auto
  names(seq_count_auto) <- markers_auto

  return(list("lus_auto" = lus_auto,
              "reps_auto" = reps_auto,
              "seq_al_auto" = seq_al_auto,
              "seq_count_auto" = seq_count_auto))
}

#' calc_weighted_ca
#'
#' @description The function to calculate the weighted corrected allele number
#' @param min_us Minimum uninterrupted stretch
#' @param reps_one_al Repeat motifs of one allele
#' @param seq_count_one_al Counts for sequence data in one allele
calc_weighted_ca <- function(min_us, reps_one_al, seq_count_one_al){
  # Calculate corrected allele numbers in each sequence
  n_seq_al <- length(reps_one_al)
  ca_one_al <- rep(0, n_seq_al)
  for(i in 1:n_seq_al){
    tmp <- reps_one_al[[i]] - min_us
    ca_one_al[i] <- sum(tmp[tmp > 0])
  }

  # Calculate the weighted corrected allele number
  weighted_ca <- sum(ca_one_al * seq_count_one_al / sum(seq_count_one_al))
  return(weighted_ca)
}

#' make_sr_ca_one_mk
#'
#' @description The function to make data of corrected allele numbers in one marker
#' @param min_us Minimum uninterrupted stretch
#' @param reps_one_mk Repeat motifs of one marker
#' @param seq_al_one_mk Alleles for sequence data in one marker
#' @param seq_count_one_mk Counts for sequence data in one marker
#' @param sr_al_one_mk A vector of the parent alleles of each stutter ratio
make_sr_ca_one_mk <- function(min_us, reps_one_mk, seq_al_one_mk, seq_count_one_mk, sr_al_one_mk){
  sr_ca_one_mk <- rep(0, length(sr_al_one_mk))
  al_names <- sort(unique(seq_al_one_mk))
  for(i in 1:length(al_names)){
    index_seq_al <- which(seq_al_one_mk == al_names[i])
    sr_ca_one_mk[which(sr_al_one_mk == al_names[i])] <- calc_weighted_ca(min_us, reps_one_mk[index_seq_al], seq_count_one_mk[index_seq_al])
  }
  return(sr_ca_one_mk)
}
