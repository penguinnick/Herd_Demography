# calculate confidence intervals for mortality profiles

# Function to take raw wear data and expand Payne Groups into separate rows
expand_payne_groups <- function(df) {
  df_expanded <- df %>%
    # Split Payne Group column by comma and trim whitespace
    mutate(Payne_Group_split = str_split(Payne.Group, ",")) %>%
    unnest(Payne_Group_split) %>%
    mutate(Payne_Group = str_trim(Payne_Group_split),
           Site_Period = paste0(Site, "_", Period) ) %>%
    select(-Payne_Group_split)
  
  return(df_expanded)
}

df_expanded <- expand_payne_groups(raw.wear)

# Function to count age classes by grouping variable
bootstrap_counts_ci <- function(df_expanded, group_var = "Site_Period", 
                                n_bootstrap = 1000, conf_level = 0.95) {
  
  # Summarize counts first
  count_summary <- df_expanded %>%
    count(!!sym(group_var), name = "count")
  
  # Bootstrap function for a single group
  boot_count <- function(data, indices) {
    sampled_data <- data[indices, ]
    sampled_counts <- sampled_data %>%
      count(!!sym(group_var), name = "count") %>%
      complete(!!sym(group_var) := unique(df_expanded[[group_var]]), 
               fill = list(count = 0)) %>%
      pull(count)
    return(sampled_counts)
  }

  # Run bootstrap
  boot_obj <- boot(data = df_expanded, statistic = boot_count, 
                   R = n_bootstrap)
  
  # Extract bootstrap results by group
  results <- tibble(
    group = unique(count_summary[[group_var]]),
    estimate = count_summary$count,
    lower = apply(boot_obj$t, 2, quantile, probs = (1-conf_level)/2, na.rm = TRUE),
    upper = apply(boot_obj$t, 2, quantile, probs = 1-(1-conf_level)/2, na.rm = TRUE)
  )
  
  return(results)
}

ci_species <- bootstrap_counts_ci(df_expanded, "Site_Period")
