# a function for plotting a histogram from counts with unequal bin widths using ggplot2 
#' @param x a vector of counts for each age class (A, B, C, D, EF, G, HI)
#' @return a histogram with frequency density on the y-axis and age in years on the x-axis, with bins corresponding to the age classes defined by Payne (1973)
#' @examples
#'   counts <- c(10, 20, 30, 25, 15, 5, 2)
#'   mortality.histogram.gg(counts)
#' @export 
#' 

mortality.histogram.gg <- function(x) {
  # Calculate relative frequency
  rel.freq <- function(x) {
    tot <- sum(x)
    x / tot
  }
  
  # Bin widths in years based on age classes (A = 0-2, B = 2-6, C = 6-12, D = 12-24, EF = 24-48, G = 48-72, HI = 72-120 months)
  bin.widths <- c(1/6, 1/3, 1/2, 1, 2, 2, 4)
  
  # Calculate frequency density
  freq.density <- rel.freq(x) / bin.widths
  
  s = 1 - cumsum(rel.freq(x))
  # Ages in years for x-axis
  breaks <- c(0, 2, 6, 12, 24, 48, 72, 120) / 12
  
  # Create data frame for plotting
  plot_data <- data.frame(
    xmin = breaks[-length(breaks)],
    xmax = breaks[-1],
    density = freq.density,
    mid = (breaks[-length(breaks)] + breaks[-1]) / 2,
    s = s,
    n = x, 
    label = c("A", "B", "C", "D", "EF", "G", "HI")
  )
  
  lab_counts <- with(plot_data,
                     paste0(label, ": n = ", round(n, 1), " (", round(rel.freq(n) * 100), "%)"))
  lab_counts <- paste(lab_counts, collapse = "\n")
  
  
  # Create ggplot
  ggplot(plot_data, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = density)) +
    geom_rect(fill = "grey80", color = "grey60", linewidth = 0.2) +
    # geom_text(aes(x = mid, y = max(density) * 1.05, label = label), 
    #           size = 4, inherit.aes = FALSE) +
    # xlim(c(0.03, 8)) +
    
    scale_x_continuous(
      name = "Age Class",
      breaks = plot_data$mid,
      labels = plot_data$label,
      # expand = c(0, 0)
      expand = expansion(mult = c(0.01, 0))
    ) +
    scale_y_continuous(
      name = "Frequency density",
      limits = c(0, 1.2)
      # expand = expansion(mult = c(0,  0.1)),
    ) +
    # geom_line( colour = rcartocolor::carto_pal(10,"Safe")[1], lwd = 1) +
    geom_line(aes(x = mid, y = s), colour = rcartocolor::carto_pal(10,"Safe")[1], size = 0.5) +
    # include counts (n) for each label as text in rows in plot area, center-right
    # geom_text(aes(x = 6, y = 1.2, label = "Counts:"), 
    #           size = 3, inherit.aes = FALSE) +
    annotate("text",
             x = 5.7, y = 0.6,
             label = lab_counts,
             size = 2.5, hjust = 0) +
    # labs(title = "Histogram from counts with unequal bin widths") +
    theme_classic() + # minimal() +
    theme(
      # plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 0, size = 6),
      axis.title = element_text(size = 8),
      axis.line = element_line( linewidth = 0.5)
    )
}
