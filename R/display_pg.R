#' aggregate_pg_contrib
#'
#' @description The function to change the format of the pg table from per combination to per contributor
#' @param dt_pg The data.table of the probabilistic genotyping in one marker
#' @param contrib The name of the selected contributor
#' @export
aggregate_pg_contrib <- function(dt_pg, contrib){
  aggregated_data <- dt_pg %>%
    group_by_at(which(names(dt_pg) == contrib)) %>%
    summarise(Probability = sum(Probability))
  aggregated_data <- aggregated_data %>%
    arrange(desc(Probability))
  return(aggregated_data)
}
