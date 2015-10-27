make_signal_classification_plots <- function(signal_vars, class_vars, data, sampling_freq, subsample_rate = 10, panel_minutes, main_title, class_var_labels,
    class_var_palette_types, class_var_palettes, font_size = 44, line_height = 2, panel_legend_relative_width = 8, plot_file_path, temp_plot_file_path) {
  if(!missing(temp_plot_file_path)) {
    if(identical(substr(temp_plot_file_path, nchar(temp_plot_file_path) - 2, nchar(temp_plot_file_path)), "pdf")) {
      pdf(temp_plot_file_path)
    } else if(identical(substr(temp_plot_file_path, nchar(temp_plot_file_path) - 2, nchar(temp_plot_file_path)), "png")) {
      png(temp_plot_file_path)
    }
    dev.control(displaylist="enable")
  }
  
  recorded_plots <- vector("list", length(signal_vars) * length(class_vars))
  record_ind <- 0
  for(signal_var in signal_vars) {
    for(class_var_ind in seq_along(class_vars)) {
      make_one_signal_classification_plot(signal_var, class_vars[class_var_ind], data, sampling_freq, subsample_rate, panel_minutes, main_title,
        class_var_labels[class_var_ind], class_var_palette_types[class_var_ind], class_var_palettes[[class_var_ind]], font_size, line_height, panel_legend_relative_width)
      record_ind <- record_ind + 1
      recorded_plots[[record_ind]] <- recordPlot()
    }
  }
  
  if(!missing(temp_plot_file_path)) {
    dev.off()
    file.remove(temp_plot_file_path)
  }
  
  if(identical(substr(plot_file_path, nchar(plot_file_path) - 2, nchar(plot_file_path)), "pdf")) {
    pdf(file = plot_file_path, width = 30, height = 20)
  } else if(identical(substr(plot_file_path, nchar(plot_file_path) - 2, nchar(plot_file_path)), "png")) {
    png(file = plot_file_path, width = 30, height = 20, units = "in", res = 300, type = "cairo-png")
  }

  lapply(recorded_plots, replayPlot)
  dev.off()
  
  return(recorded_plots)
}

make_one_signal_classification_plot <- function(signal_var, class_var, data, sampling_freq, subsample_rate, panel_minutes, main_title, signal_var_label = toupper(signal_var), class_var_label,
    class_var_palette_type, class_var_palette, font_size = 44, line_height = 2, panel_legend_relative_width = 8) {
  # subsample the data
  subsampled_signal_var <- data[seq(from = 1, to = nrow(data), by = subsample_rate), signal_var]
  subsampled_class_var <- data[seq(from = 1, to = nrow(data), by = subsample_rate), class_var]
  
  # number of subsampled data points in each panel and inds for subsampled observations in each panel
  panel_length <- round(sampling_freq * 60 * panel_minutes / subsample_rate)
  num_panels <- ceiling(length(subsampled_signal_var) / panel_length)
  panel_inds <- rep(seq_len(num_panels), each = panel_length)[seq_along(subsampled_signal_var)]
  panel_inds <- lapply(seq_len(num_panels), function(panel) which(panel_inds == panel))
  
  # update panel minutes to reflect actual value
  panel_minutes <- panel_length * subsample_rate / (sampling_freq * 60)
  
  # set up plot
  plot_layout <- grid.layout(nrow = 2 + num_panels + 1, ncol = 3, widths = unit(c(1, panel_legend_relative_width, 1), c("lines", "null", "null")),
    heights = unit(c(line_height, line_height, rep(10, num_panels)), c("lines", "lines", rep("null", num_panels))))
  
  grid.newpage()
  pushViewport(viewport(layout = plot_layout))
  
  # get min_y and max_y -- same for all panels
  temp <- pretty(subsampled_signal_var)
  min_y <- temp[1]
  max_y <- temp[length(temp)]
  
  # print plot legend.  steps are:
  #  (1) print a ggplot that includes the legend along with full data plot
  #  (2) get the legend grob
  #  (3) print the legend grob to the appropriate location in the plot
  #  (4) remove the unneeded plot created in step (1)
  # note: after printing plot, grid.ls() lists all viewports and grobs
  
  # step (1) -- print a ggplot that includes the legend along with full data plot
  panel <- 1L
  class_polygons <- get_class_polygons(subsampled_class_var[panel_inds[[panel]]], min_y, max_y)
  
  temp_plot_df <- data.frame(x = seq_along(panel_inds[[panel]]), y = subsampled_signal_var[panel_inds[[panel]]])
  
  p <- ggplot() +
    geom_polygon(aes(x = x, y = y, fill = value, group = id), data = class_polygons) +
    geom_line(aes(x = x, y = y), size = 0.1, data = temp_plot_df)# +
#    xlab("Time (minutes)") + ylab(toupper(signal_var))
  
  if(identical(class_var_palette_type, "manual")) {
    p <- p + scale_fill_manual(name = class_var_label, breaks = levels(data[, class_var]), drop = FALSE, values = class_var_palette, limits = levels(data[, class_var]))
  } else {
    p <- p + scale_fill_brewer(class_var_label, breaks = levels(data[, class_var]), drop = FALSE, type = class_var_palette_type, palette = class_var_palette)
  }

  p <- p + scale_x_continuous(limits = c(0, panel_length + 0.5), breaks = seq(from = 0, len = floor(panel_minutes / 5 + 1), by = sampling_freq * 60 * 5 / subsample_rate),
                              labels = seq(from = (panel - 1) * panel_minutes, len = floor(panel_minutes / 5 + 1), by = 5)) +
    scale_y_continuous(limits = c(min_y, max_y), breaks = seq(from = min_y, to = max_y, by = 1)) +
    theme_bw() +
    theme(legend.title = element_text(size = font_size), legend.text = element_text(size = font_size),
          axis.text.x = element_text(size = font_size), axis.text.y = element_text(size = font_size),
          axis.title.x = element_text(size = font_size), axis.title.y = element_text(size = font_size))

  print(p, vp = viewport(layout.pos.row = 2 + panel + 1, layout.pos.col = 1))
  
  # step (2) -- get the legend grob
  legend_grob <- grid.get("guide-box.3-5-3-5")
  
  # step (3) -- remove the unneeded plot created in step (1)
#  grid.remove(capture.output(grid.ls())[3])
  
  # step (4) -- set up a new grid page and print the legend grob to the appropriate location in the plot
  grid.newpage()
  pushViewport(viewport(layout = plot_layout))

  pushViewport(viewport(layout.pos.row = 2 + seq_len(num_panels), layout.pos.col = 3))
  grid.draw(legend_grob)
  upViewport()
  
  # print plot title
  grid.text(main_title, just = "center", gp = gpar(fontsize = 1.2 * font_size), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid.text(paste(signal_var_label, "and", class_var_label, "vs. Time"), just = "center", gp = gpar(fontsize = 1.2 * font_size), vp = viewport(layout.pos.row = 2, layout.pos.col = 2))

  # print X and Y axis labels
  grid.text("Time (minutes)", just = "center", gp = gpar(fontsize = font_size), vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  grid.text(signal_var_label, just = "center", rot = 90, gp = gpar(fontsize = font_size), vp = viewport(layout.pos.row = seq(from = 3, length = num_panels), layout.pos.col = 1))


  # print plot panels
  for(panel in seq_len(num_panels)) {
    class_polygons <- get_class_polygons(subsampled_class_var[panel_inds[[panel]]], min_y, max_y)
    
    temp_plot_df <- data.frame(x = seq_along(panel_inds[[panel]]), y = subsampled_signal_var[panel_inds[[panel]]])
    
    p <- ggplot() +
      geom_polygon(aes(x = x, y = y, fill = value, group = id), data = class_polygons) +
      geom_line(aes(x = x, y = y), data = temp_plot_df) +
	  xlab("") + ylab("")
#      xlab("Time (minutes)") + ylab(toupper(signal_var))
    
    if(identical(class_var_palette_type, "manual")) {
      p <- p + scale_fill_manual(name = class_var_label, breaks = levels(data[, class_var]), drop = FALSE, values = class_var_palette, limits = levels(data[, class_var]))
    } else {
      p <- p + scale_fill_brewer(class_var_label, breaks = levels(data[, class_var]), drop = FALSE, type = class_var_palette_type, palette = class_var_palette)
    }
    
    p <- p + scale_x_continuous(limits = c(0, panel_length + 0.5), breaks = seq(from = 0, len = panel_minutes / 5 + 1, by = sampling_freq * 60 * 5 / subsample_rate),
                                labels = seq(from = (panel - 1) * panel_minutes, len = panel_minutes / 5 + 1, by = 5)) +
      scale_y_continuous(limits = c(min_y, max_y), breaks = seq(from = min_y, to = max_y, by = 1)) +
      theme_bw() +
      theme(legend.position = "none",
        axis.text.x = element_text(size = font_size), axis.text.y = element_text(size = font_size),
        axis.title.x = element_text(size = font_size), axis.title.y = element_text(size = font_size))
    
    print(p, vp = viewport(layout.pos.row = 2 + panel, layout.pos.col = 1))
  }
}






make_signal_multi_classification_plots <- function(signal_vars, class_vars, data, sampling_freq, subsample_rate = 10, panel_minutes, main_title, class_var_labels,
    class_var_palette_types, class_var_palettes, font_size = 44, line_height = 2, panel_legend_relative_width = 8, plot_file_path, temp_plot_file_path) {
  if(!missing(temp_plot_file_path)) {
    if(identical(substr(temp_plot_file_path, nchar(temp_plot_file_path) - 2, nchar(temp_plot_file_path)), "pdf")) {
      pdf(temp_plot_file_path)
    } else if(identical(substr(temp_plot_file_path, nchar(temp_plot_file_path) - 2, nchar(temp_plot_file_path)), "png")) {
      png(temp_plot_file_path)
    }
    dev.control(displaylist="enable")
  }
  
  recorded_plots <- vector("list", length(signal_vars))
  record_ind <- 0
  for(signal_var in signal_vars) {
    make_one_signal_multi_classification_plot(signal_var, class_vars, data, sampling_freq, subsample_rate, panel_minutes, main_title,
      class_var_labels, class_var_palette_types, class_var_palettes, font_size, line_height, panel_legend_relative_width)
    record_ind <- record_ind + 1
    recorded_plots[[record_ind]] <- recordPlot()
  }
  
  if(!missing(temp_plot_file_path)) {
    dev.off()
    file.remove(temp_plot_file_path)
  }
  
  if(identical(substr(plot_file_path, nchar(plot_file_path) - 2, nchar(plot_file_path)), "pdf")) {
    pdf(file = plot_file_path, width = 30, height = 20)
  } else if(identical(substr(plot_file_path, nchar(plot_file_path) - 2, nchar(plot_file_path)), "png")) {
    png(file = plot_file_path, width = 30, height = 20, units = "in", res = 300, type = "cairo-png")
  }

  lapply(recorded_plots, replayPlot)
  dev.off()
  
  return(recorded_plots)
}

make_one_signal_multi_classification_plot <- function(signal_var, class_vars, data, sampling_freq, subsample_rate, panel_minutes, main_title, class_var_label,
    class_var_palette_type, class_var_palette, font_size = 44, line_height = 2, panel_legend_relative_width = 8) {
  # subsample the data
  subsampled_signal_var <- data[seq(from = 1, to = nrow(data), by = subsample_rate), signal_var]
  subsampled_class_var1 <- data[seq(from = 1, to = nrow(data), by = subsample_rate), class_vars[1]]
  subsampled_class_var2 <- data[seq(from = 1, to = nrow(data), by = subsample_rate), class_vars[2]]
  
  # number of subsampled data points in each panel and inds for subsampled observations in each panel
  panel_length <- round(sampling_freq * 60 * panel_minutes / subsample_rate)
  num_panels <- ceiling(length(subsampled_signal_var) / panel_length)
  panel_inds <- rep(seq_len(num_panels), each = panel_length)[seq_along(subsampled_signal_var)]
  panel_inds <- lapply(seq_len(num_panels), function(panel) which(panel_inds == panel))
  
  # update panel minutes to reflect actual value
  panel_minutes <- panel_length * subsample_rate / (sampling_freq * 60)
  
  # set up plot
  plot_layout <- grid.layout(nrow = 2 + num_panels, ncol = 2, widths = unit(c(panel_legend_relative_width, 1), c("null", "null")),
    heights = unit(c(line_height, line_height, rep(10, num_panels)), c("lines", "lines", rep("null", num_panels))))
  
  grid.newpage()
  pushViewport(viewport(layout = plot_layout))
  
  # get min_y and max_y -- same for all panels
  temp <- pretty(subsampled_signal_var)
  min_y <- temp[1]
  max_y <- temp[length(temp)]
  
  # print plot title
  grid.text(main_title, just = "center", gp = gpar(fontsize = 1.2 * font_size), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(paste(toupper(signal_var), "and", class_var_label, "vs. Time"), just = "center", gp = gpar(fontsize = 1.2 * font_size), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  
  # print plot legend.  steps are:
  #  (1) print a ggplot that includes the legend along with full data plot
  #  (2) get the legend grob
  #  (3) print the legend grob to the appropriate location in the plot
  #  (4) remove the unneeded plot created in step (1)
  # note: after printing plot, grid.ls() lists all viewports and grobs
  
  # step (1) -- print a ggplot that includes the legend along with full data plot
  panel <- 1L
  class1_polygons <- get_class_polygons(subsampled_class_var1[panel_inds[[panel]]], min_y, 1)
  class2_polygons <- get_class_polygons(subsampled_class_var2[panel_inds[[panel]]], 1, max_y)
  
  temp_plot_df <- data.frame(x = seq_along(panel_inds[[panel]]), y = subsampled_signal_var[panel_inds[[panel]]])
  
  p <- ggplot() +
    geom_polygon(aes(x = x, y = y, fill = value, group = id), data = class1_polygons) +
    geom_polygon(aes(x = x, y = y, fill = value, group = id), data = class2_polygons) +
    geom_line(aes(x = x, y = y), size = 0.1, data = temp_plot_df) +
    xlab("Time (minutes)") + ylab(toupper(signal_var))
  
  if(identical(class_var_palette_type, "manual")) {
    p <- p + scale_fill_manual(name = class_var_label, breaks = levels(data[, class_vars[1]]), drop = FALSE, values = class_var_palette, limits = levels(data[, class_vars[1]]))
  } else {
    p <- p + scale_fill_brewer(class_var_label, breaks = levels(data[, class_vars[1]]), drop = FALSE, type = class_var_palette_type, palette = class_var_palette)
  }

  p <- p + scale_x_continuous(limits = c(0, panel_length + 0.5), breaks = seq(from = 0, len = floor(panel_minutes / 5 + 1), by = sampling_freq * 60 * 5 / subsample_rate),
                              labels = seq(from = (panel - 1) * panel_minutes, len = floor(panel_minutes / 5 + 1), by = 5)) +
    scale_y_continuous(limits = c(min_y, max_y), breaks = seq(from = min_y, to = max_y, by = 1)) +
    theme_bw() +
    theme(legend.title = element_text(size = font_size), legend.text = element_text(size = font_size),
          axis.text.x = element_text(size = font_size), axis.text.y = element_text(size = font_size),
          axis.title.x = element_text(size = font_size), axis.title.y = element_text(size = font_size))

  print(p, vp = viewport(layout.pos.row = 2 + panel + 1, layout.pos.col = 1))
  
  # step (2) -- get the legend grob
  legend_grob <- grid.get("guide-box.3-5-3-5")
  
  # step (3) -- print the legend grob to the appropriate location in the plot
  pushViewport(viewport(layout.pos.row = 2 + seq_len(num_panels), layout.pos.col = 2))
  grid.draw(legend_grob)
  upViewport()
  
  # step (4) -- remove the unneeded plot created in step (1)
  grid.remove(capture.output(grid.ls())[3])
  
  # print plot panels
  for(panel in seq_len(num_panels)) {
    class1_polygons <- get_class_polygons(subsampled_class_var1[panel_inds[[panel]]], min_y, mean(c(min_y, max_y)))
	class2_polygons <- get_class_polygons(subsampled_class_var2[panel_inds[[panel]]], mean(c(min_y, max_y)), max_y)
    
    temp_plot_df <- data.frame(x = seq_along(panel_inds[[panel]]), y = subsampled_signal_var[panel_inds[[panel]]])
    
    p <- ggplot() +
      geom_polygon(aes(x = x, y = y, fill = value, group = id), data = class1_polygons) +
      geom_polygon(aes(x = x, y = y, fill = value, group = id), data = class2_polygons) +
      geom_line(aes(x = x, y = y), data = temp_plot_df) +
      xlab("Time (minutes)") + ylab(toupper(signal_var))
    
    if(identical(class_var_palette_type, "manual")) {
      p <- p + scale_fill_manual(name = class_var_label, breaks = levels(data[, class_vars[1]]), drop = FALSE, values = class_var_palette, limits = levels(data[, class_vars[1]]))
    } else {
      p <- p + scale_fill_brewer(class_var_label, breaks = levels(data[, class_vars[1]]), drop = FALSE, type = class_var_palette_type, palette = class_var_palette)
    }
    
    p <- p + scale_x_continuous(limits = c(0, panel_length + 0.5), breaks = seq(from = 0, len = panel_minutes / 5 + 1, by = sampling_freq * 60 * 5 / subsample_rate),
                                labels = seq(from = (panel - 1) * panel_minutes, len = panel_minutes / 5 + 1, by = 5)) +
      scale_y_continuous(limits = c(min_y, max_y), breaks = seq(from = min_y, to = max_y, by = 1)) +
      theme_bw() +
      theme(legend.position = "none",
        axis.text.x = element_text(size = font_size), axis.text.y = element_text(size = font_size),
        axis.title.x = element_text(size = font_size), axis.title.y = element_text(size = font_size))
    
    print(p, vp = viewport(layout.pos.row = 2 + panel, layout.pos.col = 1))
  }
}







get_class_polygons <- function(class_var, min_y, max_y) {
  n <- length(class_var)
  if(n < 2)
    stop("Invalid length for class_var: must be >= 2")
  
  class_change_points <- which(class_var[1:(n - 1)] != class_var[2:n])
  class_startend_points <- rep(class_change_points, each = 2) + rep(c(0, 1), times = length(class_change_points))
  class_startend_points <- c(1, class_startend_points, n)
	class_startend_points <- class_startend_points + rep(c(-0.5, 0.5), times = length(class_change_points) + 1)

  class_polygon_positions <- data.frame(
    id = rep(seq(length(class_change_points) + 1), each = 4),
    x = rep(class_startend_points, each = 2),
    y = rep(c(min_y, max_y, max_y, min_y), length = 2 * length(class_startend_points))
  )
  
  class_polygon_values <- data.frame(
    id = seq(length(class_change_points) + 1),
    value = class_var[c(1, class_change_points + 1)]
  )
  
  class_polygons <- merge(class_polygon_values, class_polygon_positions, by = "id")
  
  return(class_polygons)
}
