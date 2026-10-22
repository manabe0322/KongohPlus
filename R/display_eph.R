#' create_fig_eph
#'
#' @description The function to create the figure for expected peak heights in each contributor
#' @param dt_param The data.table for estimated parameters
#' @param min_size The minimum size
#' @param max_size The maximum size
#' @param base_size The base size
#' @export
create_fig_eph <- function(dt_param, min_size, max_size, base_size){
  min_size_round <- round(min_size, 0)
  max_size_round <- round(max_size, 0)
  x_val <- min_size_round:max_size_round

  temp_total <- dt_param[["Estimated value"]][dt_param$Parameter == "Template amount"]
  mp <- dt_param[["Estimated value"]][dt_param$Parameter == "Mixture proportion"]
  temp <- temp_total * mp
  deg <- dt_param[["Estimated value"]][dt_param$Parameter == "Degradation"]
  rep_ae <- dt_param[["Estimated value"]][dt_param$Parameter == "Replicate amplification efficiency"]
  hnc <- length(temp)
  n_rep <- length(rep_ae)

  colors <- c(U1 = "#00BFC4", U2 = "#7CAE00", U3 = "#E9BC00", U4 = "#F8766D")

  if(n_rep == 1){
    y_list <- list()
    for(i in 1:hnc){
      y_list[[i]] <- temp[i] * deg[i]^((x_val - base_size) / 100)
    }
    names(y_list) <- paste0("U", 1:hnc)
    max_y_plot <- max(unlist(y_list))

    base_plot <- ggplot() + xlab("Size") + ylab("Average peak height")
    for(nm in names(y_list)){
      dt_plot <- data.table(X = x_val, Y = y_list[[nm]], group = nm)
      base_plot <- base_plot +
        geom_line(
          data = dt_plot,
          aes(
            x = X, y = Y,
            color = group
          )
        )
    }

    base_plot <- base_plot +
      scale_color_manual(values = colors, name = "Contributor") +
      coord_cartesian(ylim = c(0, max_y_plot))
  }

  else{
    x_all <- numeric(0)
    y_min_all <- numeric(0)
    y_max_all <- numeric(0)
    group_all <- character(0)
    for(i in 1:hnc){
      y_min <- min(rep_ae) * temp[i] * deg[i]^((x_val - base_size) / 100)
      y_max <- max(rep_ae) * temp[i] * deg[i]^((x_val - base_size) / 100)

      x_all <- c(x_all, x_val)
      y_min_all <- c(y_min_all, y_min)
      y_max_all <- c(y_max_all, y_max)
      group_all <- c(group_all, rep(paste0("U", i), length(x_val)))
    }
    max_y_plot <- max(y_max_all)

    dt_plot <- data.table(X = x_all, Ymin = y_min_all, Ymax = y_max_all, Group = group_all)

    base_plot <- ggplot(dt_plot, aes(x = X, color = Group, fill = Group)) +
      geom_ribbon(aes(ymin = Ymin, ymax = Ymax), alpha = 0.3, linewidth = 0) +
      geom_line(aes(y = Ymin), size = 1) +
      geom_line(aes(y = Ymax), size = 1) +
      scale_fill_manual(values = colors, name = "Contributor") +
      scale_color_manual(values = colors, name = "Contributor", guide = "none") +
      coord_cartesian(ylim = c(0, max_y_plot))
  }

  return(base_plot)
}
