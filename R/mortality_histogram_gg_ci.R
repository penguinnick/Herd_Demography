# a function for plotting a histogram from counts with unequal bin widths using ggplot2 
#' @param x a data.frame containing a column with of Payne's age classes, with multiple age classes separated by a comma. 
#' @param CI logical. Whether to calculate confidence intervals for the mortality profile using bootstrapping. Default TRUE.
#' @param iter integer. Number of bootstrap iterations to perform when calculating confidence intervals. Default= 1000
#' @param ci numeric. Confidence level for intervals (e.g., 95 for 95% confidence intervals). Default = 95.
#' @param probability.correction logical. Whether to correct counts based on Brochier's (2013) probability correction rule. Default TRUE. If TRUE, combines age groups EF, and G, and HI into a single group (EF-G-HI) with a probability correction factor of 2 for EF and 4 for HI.
#' @return a histogram with frequency density on the y-axis and age in years on the x-axis, with bins corresponding to either age classes defined by Payne (1973) or separates EF, and HI (if probability correction is FALSE)
#' @examples
#'   counts <- c(10, 20, 30, 25, 15, 5, 2)
#'   mortality.histogram.gg(counts)
#' @export 
#' 

mortality.histogram.gg.CI <- function(df, CI = TRUE, iter = 1000, ci = 95, probability.correction = TRUE) {
  
  expand_payne_groups <- function(df, probability.correction) {
    df_expanded <- df %>%
      # Split Payne Group column by comma and trim whitespace
      mutate(Payne_Group_split = str_split(Payne.Group, ",")) %>%
      unnest(Payne_Group_split) %>%
      mutate(Payne_Group = str_trim(Payne_Group_split),
             Site_Period = paste0(Site, "_", Period)) %>%
      select(-Payne_Group_split) 
    
    if(probability.correction) {
      df_expanded <- df_expanded %>%
        mutate(Ageclass = case_when(
          Payne_Group == "A" ~ 1,
          Payne_Group == "B" ~ 2,
          Payne_Group == "C" ~ 3,
          Payne_Group == "D" ~ 4,
          Payne_Group == "E" ~ 5,
          Payne_Group == "F" ~ 5,
          Payne_Group == "G" ~ 6,
          Payne_Group == "H" ~ 7,
          Payne_Group == "I" ~ 7,
          .default = 0))
      } else {
          df_expanded <- df_expanded %>%
            mutate(Ageclass = case_when(
              Payne_Group == "A" ~ 1,
              Payne_Group == "B" ~ 2,
              Payne_Group == "C" ~ 3,
              Payne_Group == "D" ~ 4,
              Payne_Group == "E" ~ 5,
              Payne_Group == "F" ~ 6,
              Payne_Group == "G" ~ 7,
              Payne_Group == "H" ~ 8,
              Payne_Group == "I" ~ 9,
              .default = 0
            ))
        }
    return(df_expanded)
  }
  
  df_expanded <- expand_payne_groups(df, probability.correction )
  dat = df_expanded$Ageclass
  # df.c <- correct.counts(df_expanded$Payne_Group, probability.correction = FALSE)
  df.c <- correct.counts(df_expanded$Payne_Group, probability.correction = probability.correction)
  
  x <- df.c$n
  
  # # Calculate relative frequency
  rel.freq <- function(x) {
    tot <- sum(x)
    x / tot
  }
  
  df.c$qx <- rel.freq(x)
  # Bin widths in years based on age classes (A = 0-2, B = 2-6, C = 6-12, D = 12-24, EF = 24-48, G = 48-72, HI = 72-120 months)
  if(probability.correction) {
    bin.widths <- c(1/6, 1/3, 1/2, 1, 2, 2, 4)
  } else {
    bin.widths <- c(1/6, 1/3, 1/2, 1, 1, 1, 1, 1, 1)
  }
  # Calculate frequency density
  df.c$fd <- rel.freq(x) / bin.widths
  freq.density <- df.c$fd
  
  s = 1 - cumsum(rel.freq(x))
  # s = 1 - cumsum(df.c$qx)
  # Ages in years for x-axis
  if(probability.correction) {
     breaks <- c(0, 2, 6, 12, 24, 48, 72, 120) / 12
  } else {
     breaks <- c(0, 2, 6, 12, 24, 36, 48, 60, 72, 84) / 12 # for probability correction with separate EF and G groups
  }
  
  # plot_data$label = ifelse(probability.correction, c("A", "B", "C", "D", "EF", "G", "HI"), LETTERS[1:9])
  
  # adapted function from ZooArch Package
  mortprof <- function(d){
    N.ages = ifelse(probability.correction, 7, 9)
    vector <- rep(NA, N.ages)   
    for(i in 1:N.ages) {
      vector[i] <- sum(d==i)/length(d)
    }
    vector[is.na(vector)] <- 0
    round(vector,4) 
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    xmin = breaks[-length(breaks)],
    xmax = breaks[-1],
    density = freq.density,
    mid = (breaks[-length(breaks)] + breaks[-1]) / 2,
    s = s,
    n = x, 
    qx = rel.freq(x),
    label = df.c$v # ifelse(probability.correction, c("A", "B", "C", "D", "EF", "G", "HI"), LETTERS[1:9]) # c("A", "B", "C", "D", "EF", "G", "HI")
    # label = LETTERS[1:9]
  )
  
  #-- confidence intervals for mortality profiles 
  #-- function adapted from Zooarch Package
  if(CI) {
    mortality.matrix <- matrix(NA, ncol = length(bin.widths), nrow = iter)
    mortality.matrix[1,] <- rel.freq(x)
    for(i in 2:iter){
      bootstrap <- sample(1:length(dat), length(dat), replace = TRUE)
      mortality.matrix[i,] <- unlist(mortprof(dat[bootstrap]))
    }  
    ci<-ci/100  
    upCI<-apply(mortality.matrix[,], MARGIN = 2, FUN = quantile, 
                probs = ci+((1-ci)/2))
    loCI<-apply(mortality.matrix[,], MARGIN = 2, FUN = quantile, 
                    probs = ((1-ci)/2))
    # # Add confidence intervals to plot_data
    # ci_data <- ci_species %>%
    #   filter(group == "All") %>%
    #   select(estimate, lower, upper) %>%
    #   slice(rep(1:n(), each = nrow(plot_data))) %>%
    #   bind_cols(plot_data)
    
    plot_data$upCI <- upCI
    plot_data$loCI <- loCI
  }
  
  lab_counts <- with(plot_data,
                     paste0(label, ": n = ", round(n, 1), " (", round(rel.freq(n) * 100), "%)"))
  lab_counts <- paste(lab_counts, collapse = "\n")
  
  # Create ggplot
  ggplot(plot_data, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = qx)) +
    geom_rect(fill = "grey80", color = "grey60", linewidth = 0.2) +
    scale_x_continuous(
      name = "Age Class",
      breaks = plot_data$mid,
      labels = plot_data$label,
      expand = expansion(mult = c(0.01, 0))
    ) +
    scale_y_continuous(
      name = "Rel Frequency",
      limits = c(0, 1.2)
    ) +
    geom_line(aes(x = mid, y = s), colour = rcartocolor::carto_pal(10,"Safe")[1], size = 0.5) +
    annotate("text",
             x = 5.7, y = 0.7,
             label = lab_counts,
             size = 2.5, hjust = 0) +
    {
      if(CI) geom_errorbar(aes(x = mid, ymin = qx - loCI, ymax = qx + upCI), width = 0.2, color = "black", size = 0.3)
      } +
    theme_classic() + 
    theme(
      axis.text.x = element_text(angle = 0, size = 6),
      axis.title = element_text(size = 8),
      axis.line = element_line( linewidth = 0.5)
    )
}
