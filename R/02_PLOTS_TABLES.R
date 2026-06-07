
###########################################################################
#                                                                         #
#                               FIGURES                                   #
#                                                                         #                 
###########################################################################
#
#-- create directory for storing Figures
if(!dir.exists("./Figures")){
  dir.create("./Figures")
}
figs.dir = paste0(r_path,"Figures/")
#-- Figure 1: Offtake models
ggsave(
  # plot = HerdDynamics::culling_multiplot(HerdDynamics::offtake_models),
  plot = culling_multiplot2(offtake_models),
  filename = paste0(figs.dir,"Fig1_TheoreticalCullingStrategies_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)

#-- Figure 2: Archaeological Mortality Profiles
mort.plots = lapply(ages.list, function(a){mortality.histogram.gg(a$n)})
#-- rearrange plots to match order of strategies in Fig 1
mort.plots = c(mort.plots[3], mort.plots[4], mort.plots[1], mort.plots[2], mort.plots[5])
#-- put plots together in a grid
Fig2 = cowplot::plot_grid(plotlist = mort.plots, labels = names(mort.plots), 
                          label_size = 10, 
                          vjust = 1.5, hjust = c(-1.5, -1.5, -0.75, -0.75, -0.75))
ggsave(
  plot = Fig2,
  filename = paste0(figs.dir,"Fig2_ArchaeologicalMortProfiles_v2.jpg"),
  dpi = 800,
  width = 8,
  height = 4
)

#-- Figure 3: Simple plot bootstrapped lambda results
Fig3 = repro.boot.df %>% facet_factor_fun() %>% 
  ggplot(aes(x = strategy, y = Lambda, color = Taxon)) +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = cbPalette[10]) +
  geom_point(size = 2,
             position = position_dodge(width = 0.5),
             alpha = 0.8) +
  geom_errorbar(aes(ymin = low.lambda, ymax = up.lambda),
                width = 0.5,
                position = "dodge") +
  scale_color_manual(
    breaks = c("goat", "sheep"), 
    values = taxon_pal, #c("#d95f02", "#1b9e77"),
    labels = c("Goat", "Sheep")
    # values = taxon_pal cbPalette[c(1, 2)]
    ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(
  Fig3,
  filename = paste0(figs.dir,"Fig3_lambda-boot_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)

#-- Figure 4: Herd structure at initialization of projections
Fig4 = xini.df %>% facet_factor_fun() %>% 
  filter(n > 0) %>%
  ggplot() +
  # geom_col( aes( x = AgeClass, y = n, fill = Taxon ), position = "dodge" ) +
  geom_col(aes(x = Age, y = n, fill = Taxon), position = "dodge", alpha = 0.4) +
  scale_fill_manual(
    breaks = c("goat", "sheep"), 
    values = taxon_pal, #c("#d95f02", "#1b9e77"),
    labels = c("Goat", "Sheep")
    # values = cbPalette, labels = c("Goat", "Sheep")
    ) +
  # xlab( "Age Class" ) +
  xlab("Age (years)") +
  ylab("Number of Animals") +
  # scale_y_continuous( breaks = c(0, 25, 50)) +
  facet_manual(facets = vars(Author, strategy), design = design, strip = strip_nested()) +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    line = element_line(linewidth = 0.2),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3), 
    legend.title = element_blank(),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = "white", color = "grey60"),
    panel.spacing = unit(0.2, "lines")
  ) 

ggsave(
       Fig4,
       filename = paste0(figs.dir,"Fig4_herd_structure_v2.jpg"),
       dpi = 800,
       width = 6,
       height = 4
       )

#-- Figure 5: Projections of total population size through time for each strategy
Fig5 = tot.pop.df %>%
  mutate(Taxon = factor(Taxon, levels = c("goat", "sheep"))) %>%
  facet_factor_fun() %>% 
  # mutate(strategy = factor(strategy, levels = c(names(offtake.mortality))))  %>%
  ggplot() +
  geom_col(aes(time, pop, fill = Taxon), alpha = 0.4) +
  scale_fill_manual(
    breaks = c("goat", "sheep"),
    values = taxon_pal, #cbPalette[c(1, 2)],
    labels = c("Goat", "Sheep")
  ) +
  labs(x = "Time (years)", y = "Number of animals") +
  facet_manual(
    facets = vars(Author, strategy), 
    design = design, 
    strip = strip_nested(),
    scales = "free_y") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    line = element_line(linewidth = 0.2),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3), 
    legend.title = element_blank(),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = "white", color = "grey60"),
    panel.spacing = unit(0.2, "lines")
  ) 

ggsave(
  plot =  Fig5,
  filename = paste0(figs.dir,"Fig5_projections_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)

#-- Figure 6: Projections of total population size through time for each strategy with stochastic replicates
Fig6 <- summarize.pop.df(sto.res.df)  %>% 
  facet_factor_fun() %>% 
  mutate(Taxon = taxon) %>%
  ggplot(aes(x = time, y = mean)) +
  geom_ribbon(
    data = . %>% filter(Taxon == "sheep"),
    aes(ymin = low, ymax = up, fill = Taxon),
    alpha = 0.2, 
  ) +
  geom_ribbon(
    data = . %>% filter(Taxon == "goat"),
    aes(ymin = low, ymax = up, fill = Taxon),
    alpha = 0.2, 
  ) +
  labs(x = "Year", y = "Population size") +
  geom_line(
    data = . %>% filter(Taxon == "sheep"),
    aes(x = time, y = mean, color = Taxon), 
    lwd = 0.2,
    alpha = 0.8, show.legend = FALSE
  ) +
  geom_line(
    data = . %>% filter(Taxon == "goat"),
    aes(x = time, y = mean, color = Taxon), 
    lwd = 0.2,
    alpha = 0.8, show.legend = FALSE
  ) + 
  facet_manual(
    facets = vars(Author, strategy), 
    design = design, strip = strip_nested(), scales = "free_y"
  ) +
  ylab("Number of Animals") +
  xlab("Year") +
  scale_fill_manual(
    breaks = c("goat", "sheep"),
    values = taxon_pal, #  c("#d95f02", "#1b9e77"),
    labels = c("Goat", "Sheep")) +
  scale_color_manual(
    breaks = c("goat", "sheep"),
    values = taxon_pal, # c("#d95f02", "#1b9e77"),
    labels = c("Goat", "Sheep")) + 
  theme_classic() +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    line = element_line(linewidth = 0.2),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3), 
    legend.title = element_blank(),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = "white", color = "grey60"),
    panel.spacing = unit(0.2, "lines")
  ) 

ggsave(
  plot =  Fig6,
  filename = paste0(figs.dir,"Fig6_stochastic_rep_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)

#-- Figure 7: New projections with optimized culling rates for each strategy
Fig7 = all.res.df %>%  
  group_by(time, strategy, offtake) %>%
  reframe(time, strategy, offtake, Taxon, pop, tot.pop = sum(pop)) %>% unique() %>%
  facet_factor_fun() %>%
  # mutate(strategy = factor(strategy, levels = names(offtake.mortality)),
  #        Taxon = factor(Taxon, levels = c("goat", "sheep"))) %>%
  ggplot(aes(group = Taxon)) +
  geom_col(
    data = . %>% filter(offtake == "adjusted"),
    aes(time, pop, fill = Taxon),
    alpha = 0.4,
  ) +
  theme_minimal() +
  ylab("Number of Animals") +
  xlab("Time (years)") +
  scale_fill_manual(
    breaks = c("goat", "sheep"),
    values = taxon_pal,  #cbPalette[c(1, 2)],
    labels = c("Goat", "Sheep")
  ) +
  facet_manual(
    facets = vars(Author, strategy), 
    design = design, 
    strip = strip_nested(),
    scales = "free_y") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    line = element_line(linewidth = 0.2),
    legend.position = "inside",
    legend.position.inside = c(0.8, 0.3), 
    legend.title = element_blank(),
    strip.text = element_text(size = 6),
    strip.background = element_rect(fill = "white", color = "grey60"),
    panel.spacing = unit(0.2, "lines")
  ) 

ggsave(
  plot =  Fig7,
  filename = paste0(figs.dir,"Fig7_reprojections_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)

#-- Figure 8: Mean herd size over 200 years with optimized offtake for each strategy
Fig8 =  new.vprod %>%
  facet_factor_fun() %>% 
  # new.vprod %>%
  filter(strategy != "Baseline") %>%
  ggplot(aes(x = strategy, y = xmean, fill = Taxon)) +
  geom_col(position = position_dodge(width = c(1)), alpha = 0.4) +
  geom_hline(yintercept = 150, linetype = "dashed") +
  ylab("Mean herd size over 200 years with optimized offtake") +
  xlab("Strategy") +
  theme_minimal() +
  annotate(geom = "text",
           x = 8,
           y = 135,
           label = "Initial \n Herd Size") +
  scale_fill_manual(
    breaks = c("goat", "sheep"),
    values = taxon_pal, #cbPalette,
    labels = c("Goat", "Sheep")
  ) +
  coord_flip() +
  theme(legend.position = "bottom")

ggsave(
  filename = paste0(figs.dir,"Fig8_mean_herd_size_v2.jpg"),
  plot =  Fig8,
  dpi = 800,
  width = 6,
  height = 4
)

#-- Figures 9 and 10: Herd Growth rates and Lambda for each strategy with optimized offtake for goats and sheep, respectively
ggsave(
  plot =  plot.m("goat"),
  filename = paste0(figs.dir,"Fig9_Lambda_goat_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)
# Fig 10
ggsave(
  plot =  plot.m("sheep"),
  filename = paste0(figs.dir,"Fig10_Lambda_sheep_v2.jpg"),
  dpi = 800,
  width = 6,
  height = 4
)

#-- Fig 11

#-- Fig 11 grid
fig11.design = matrix(c(
  1, 2,
  3, 4,
  5, 6 ), nrow = 3, byrow = TRUE )

Fig11 <- sobol_indices %>%
  facet_factor_fun() %>% 
  mutate(res = factor(res, levels = c("mean", "sd"))) %>% 
  ggplot(aes(group = taxon, colour = param)) +
  geom_point(
    # aes(shape = taxon, x = param, y = total_effect),
    aes(shape = taxon, x = param, y = main_effect),
    position = position_dodge(width = 0.8), 
    size = 2,
    alpha = 0.6
  ) +
  geom_errorbar(
    aes(x = param, 
        # ymin = total_effect_low, ymax = total_effect_high
        ymin = main_effect_low, ymax = main_effect_high
    ),
    position = position_dodge(width = 0.8),
    width = 0.4) +
  xlab("Sensitivity index") +
  ylab("Main effect on Herd Size") +
  scale_color_manual(
    values = scales::hue_pal()(length(param_labels)),
    breaks = names(param_labels),
    labels = param_labels 
  ) +
  theme_minimal() +
  facet_manual(
    facets = vars( strategy, res),
    design = fig11.design,
    strip = strip_nested(),
    scales = "free_y") +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_blank(), # element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    line = element_line(linewidth = 0.2),
    # legend.position = "inside",
    legend.position.inside = c(0.8, 0.3), 
    legend.title = element_blank(),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "white", color = "grey60"),
    panel.spacing = unit(0.2, "lines"),
    guides(colour = guide_legend(override.aes = list(linetype = 0)),
           shape  = guide_legend(override.aes = list(linetype = 0)))
  )
# save Fig 11
ggsave(
  plot =  Fig11 , filename = "./Figures/Fig11_main_effects.jpg",
  # filename = paste0(figs.dir,"Fig11_Sensitivity_indices_mean_pop.jpg"),
  dpi = 300,
  width = 6,
  height = 4
)

Fig12 <- sobol_indices %>%
  facet_factor_fun() %>% 
  mutate(res = factor(res, levels = c("mean", "sd"))) %>% 
  ggplot(aes(group = taxon, colour = param)) +
  geom_point(
    aes(shape = taxon, x = param, y = total_effect),
    # aes(shape = taxon, x = param, y = main_effect),
    position = position_dodge(width = 0.8), 
    size = 2,
    alpha = 0.6
  ) +
  geom_errorbar(
    aes(x = param, ymin = total_effect_low, ymax = total_effect_high),
    position = position_dodge(width = 0.8),
    width = 0.4) +
  xlab("Sensitivity index") +
  ylab("Total effect on Herd Size") +
  scale_color_manual(
    values = scales::hue_pal()(length(param_labels)),
    breaks = names(param_labels),
    labels = param_labels 
  ) +
  theme_minimal() +
  facet_manual(
    facets = vars( strategy, res),
    design = fig11.design,
    strip = strip_nested(),
    # scales = "free_y"
  ) +
  theme(
    axis.title = element_text(size = 8),
    axis.text.x = element_blank(), # element_text(hjust = 1, size = 6),
    axis.text.y = element_text(size = 8),
    line = element_line(linewidth = 0.2),
    # legend.position = "inside",
    legend.position.inside = c(0.8, 0.3), 
    legend.title = element_blank(),
    strip.text = element_text(size = 10),
    strip.background = element_rect(fill = "white", color = "grey60"),
    panel.spacing = unit(0.2, "lines"),
    guides(colour = guide_legend(override.aes = list(linetype = 0)),
           shape  = guide_legend(override.aes = list(linetype = 0)))
  )

# save Fig 12 
ggsave(
  plot =  Fig12 , filename = "./Figures/Fig12_total_effects.jpg",
  dpi = 300,
  width = 6,
  height = 4
)
###########################################################################
#                                                                         #
#                               TABLES                                    #
#                                                                         #                 
###########################################################################
#-- create Tables directory
if(!dir.exists(paste0(r_path,"tables"))){
  dir.create(paste0(r_path,"tables"))
}
table.dir = paste0(r_path,"tables/")
#-- These are default settings for tables created with flextable package. 
set_flextable_defaults(
  font.family = "Helvetica",
  font.size = 10,
  padding = 2,
  line_spacing = 1.2,
  table.layout = "autofit"
)

#-- properties for portrait table word doc output
sect_properties_portrait <- prop_section(
  page_size = page_size(
    orient = "portrait",
    width = 8.5, height = 11
  ),
  type = "continuous",
  page_margins = page_mar(left = 1)
)
#-- properties for landscape table word doc output
sect_properties_landscape <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 8.5,
    height = 11
  ),
  type = "continuous",
  page_margins = page_mar(top = 1)
)

#-- create Table 1 of parameters used in Lefkovitch matrix
t1 = as_grouped_data(
  param.dat %>%
    mutate("Parturition Rate" = ParturitionRate) %>%
    select(Taxon, Age, Sex, "Parturition Rate", Prolificacy, Mortality) %>%
    pivot_wider(
      names_from = c(Sex),
      values_from = c("Parturition Rate", Prolificacy, Mortality)
    ) %>%
    select(
      Taxon,
      Age,
      `Parturition Rate_Female`,
      Prolificacy_Female,
      contains("Mortality")
    ),
  groups = "Taxon"
) %>% 
  as_flextable() %>%
  separate_header() %>%
  bold(
    j = 1,
    i = ~ !is.na(Taxon),
    bold = TRUE,
    part = "body"
  ) %>% # make Taxon bold
  bold(part = "header", bold = TRUE) %>%  # make headers bold
  colformat_double(j = c(4, 5), digits = 3) %>% # set mortality column to 3 digits
  colformat_double(j = c(2, 3), digits = 2) %>%
  autofit() %>%
  set_table_properties(
    layout = "autofit",
    align = "center",
    # opts_word = list("split" = FALSE, "keep_with_next" = TRUE)
  ) %>%
  paginate(
    group = "Taxon",
    # keep Sample.ID group across pages
    group_def = "rle"
  ) %>%
  theme_vanilla() %>%
  align(j = 1, align = "left") %>% 
  set_caption(
  "Table 1. Parameters used to construct the Lefkovitch population projection matrix for sheep and goats. Parturition and prolificacy rates one possible outcome of the sampling procedure."
)

#-- create Table 2 of survivorship probabilities for each culling strategy 
T2 = cbind.data.frame(c("", Payne_ages$ageClasses) ,
                      do.call(cbind.data.frame, args = c(offtake_models)))
colnames(T2)[1] = "Age Class" # = ageClasses # c(0, ages)

T2ft = flextable(T2[-1, ]) %>%
  set_header_labels(
    MeatA = "Meat A",
    MeatB = "Meat B",
    MilkA = "Milk A",
    MilkB = "Milk B"
  ) %>%
  footnote(
    i = 1,
    j = 2:3,
    value = as_paragraph(c("Redding (1981)")),
    ref_symbols = c("a"),
    part = "header",
    inline = T
  ) %>%
  footnote(
    i = 1,
    j = c(4:6),
    value = as_paragraph(c("Payne (1973)")),
    ref_symbols = c("b"),
    part = "header",
    inline = T
  ) %>%
  footnote(
    i = 1,
    j = 7:11,
    value = as_paragraph(c("Vigne and Helmer (2007)")),
    ref_symbols = c("c"),
    part = "header",
    inline = T
  ) %>%
  autofit() %>%
  set_table_properties(
    layout = "autofit",
    align = "center"
    # opts_word = list("split" = FALSE, "keep_with_next" = TRUE)
  ) %>%
  paginate(
    group = "Age Class",
    # keep Sample.ID group across pages
    group_def = "rle"
  ) %>%
  theme_vanilla() %>%
  set_caption(
    "Table 2. Percentage of herd predicted to survive under ten theoretical culling strategies associated with different production strategies, as presented by Marom and Bar-Oz (2009a). Survival probabilities derived from aRedding (1981), bPayne (1973), and cVigne and Helmer (2007) correspond to Payne's age class system."
  )


#-- Table 3: survival probabilities from culling data allocated to distribution of ages used in model
T3 = cbind.data.frame(c("", Payne_ages$ageClasses) ,
                      do.call(cbind.data.frame, args = c(culling.profiles)))
colnames(T3)[1] = "Age Class"

T3ft = flextable(T3[-1, ]) %>%
  colformat_double(digits = 2) %>%
  set_header_labels(ageClasses = "Age Class", Age = "Age (years)") %>%
  autofit() %>%
  set_table_properties(
    layout = "autofit",
    align = "center",
  ) %>%
  paginate(
    group = "Age Class",
    # keep Sample.ID group across pages
    group_def = "rle"
  ) %>%
  theme_vanilla() %>% 
  set_caption("Table 3. Survival probabilities for caprines at four Neolithic Sites in Dalmatia.")



#-- Tables 4 and 5: bootstrapped Lambda estimates and sex proportions for each strategy
colnames(repro.boot.df)[c(6, 7)] = c("Proportion Female", "Proportion Male")

#-- create table 3 (theoretical harvest profiles)
t4 = repro.boot.df %>%
  filter(!strategy %in% names(culling.profiles)) %>%
  select(Taxon, strategy, Lambda, "Proportion Female", "Proportion Male") %>%
  pivot_longer(cols = c(3:5)) %>%
  pivot_wider(names_from = c(strategy))

#-- create table 4 (empirical harvest profiles)
t5 = repro.boot.df %>%
  filter(strategy %in% names(culling.profiles)) %>%
  select(Taxon, strategy, Lambda, "Proportion Female", "Proportion Male") %>%
  pivot_longer(cols = c(3:5)) %>%
  pivot_wider(names_from = c(strategy))

#-- use lambda.tables function to format tables 4 and 5 for output
t4 = lambda.tables(t4)
t5 = lambda.tables(t5)

#-- Table S4.1 - summary of mandibles used to calculate age-at-death distributions for ovicaprids at four Neolithic sites in Dalmatia.
TableS4.1 <- wear.df %>% 
  filter(Payne.Group != "") %>%
  group_by(Site, Period, Species) %>%
  reframe(N = n()) %>%
  pivot_wider(names_from = Species, values_from = N) %>% 
  group_by(Site, Period) %>%
  reframe(Site, Period, 
          N = sum(CA, OA, OC),
          "Goat:Sheep:Indeterminate" = paste0(CA, ":", OA, ":", OC)) %>% 
  as_flextable() %>%
  set_caption("Table S4. Sample information for mandibles used to calculate age-at-death distributions for ovicaprids at four Neolithic sites in Dalmatia.") %>% 
  theme_vanilla()

#-- create output word docs
save_as_docx(t1, path = "./tables/Table1_parameters.docx", pr_section = sect_properties_portrait)
save_as_docx(T2ft, path = "./tables/Table2_survival_rates.docx", pr_section = sect_properties_portrait)
save_as_docx(T3ft, path = "./tables/Table3_Site_survival_rates.docx", pr_section = sect_properties_portrait)
save_as_docx(t4, path = "./tables/Table4_lambda_theory.docx", pr_section = sect_properties_landscape)
save_as_docx(t5, path = "./tables/Table5_lambda_culling.docx", pr_section = sect_properties_landscape)
save_as_docx(TableS4.1, path = "./tables/TableS4-1_mandible_counts.docx", pr_section = sect_properties_portrait)
