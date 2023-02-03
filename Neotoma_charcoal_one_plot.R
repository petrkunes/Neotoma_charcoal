#### Charcoal workshop - neotoma2 ####


# Install and load packages ------------------------------------------------------------

library(devtools)
devtools::install_github('NeotomaDB/neotoma2')
library(neotoma2)
library(riojaPlot)

# Get just one site ------------------------------------------------

# Select Prášilské jezero (CZ) and extract samples
PRA_ds <- get_datasets(get_sites(sitename = "Prášilské jezero"), all_data = TRUE)
PRA_pol_ch <- PRA_ds %>% neotoma2::filter(datasettype %in% c("pollen", "charcoal"))
PRA_rec <- PRA_pol_ch %>% get_downloads()
PRA_samp <- PRA_rec %>% samples()

# Select only taxa identified from pollen (and only trees/shrubs) for plotting
PRA_plot_taxa <- taxa(PRA_rec) %>% dplyr::filter(ecologicalgroup %in% c("TRSH")) %>%
  dplyr::filter(elementtype == "pollen") %>%
  dplyr::arrange(desc(samples)) %>% 
  head(n = 15) 

# Select charcoal area
PRA_charcoal <- samples(PRA_rec) %>% dplyr::filter(ecologicalgroup == "CHAR" & elementtype == "area") %>%
  dplyr::select(depth, age, variablename, elementtype, value) %>%
  dplyr::mutate(variable = paste(variablename, elementtype, sep = " "), .keep = "unused") %>%
  dplyr::group_by(depth, age) %>%
  dplyr::arrange(depth)

# Select pollen taxa for calculation pollen sum and %
PRA_taxa_sum <- taxa(PRA_rec) %>% dplyr::filter(ecologicalgroup %in% c("TRSH", "UPHE")) %>%
  dplyr::filter(elementtype == "pollen") %>%
  dplyr::arrange(desc(samples))

# Clean up. Select only pollen measured using NISP.
PRA_samp_short <- samples(PRA_rec) %>% 
  dplyr::filter(variablename %in% PRA_taxa_sum$variablename) %>% 
  dplyr::filter(ecologicalgroup %in% c("TRSH","UPHE")) %>%
  dplyr::filter(elementtype == "pollen") %>%
  dplyr::filter(units == "NISP")

# Transform to proportion values.
PRA_pollen_perc <- PRA_samp_short %>%
  dplyr::group_by(depth, age) %>%
  dplyr::mutate(pollencount = sum(value, na.rm = TRUE)) %>%
  dplyr::group_by(variablename) %>% 
  dplyr::mutate(prop = value / pollencount) %>% 
  dplyr::arrange(desc(age))

# Spread the data to a "wide" table, with taxa as column headings. Both pollen and charcoal 
PRA_pollen_perc <- PRA_pollen_perc %>%
  dplyr::select(depth, age, variablename, prop) %>% 
  dplyr::mutate(prop = as.numeric(prop))

PRA_pollen_wide <- tidyr::pivot_wider(PRA_pollen_perc,
                             id_cols = c(depth, age),
                             names_from = variablename,
                             values_from = prop,
                             values_fill = 0)

PRA_charcoal_wide <- tidyr::pivot_wider(PRA_charcoal,
                                        id_cols = c(depth, age),
                                        names_from = variable,
                                        values_from = value,
                                        values_fill = 0)

# Plot stratigraphic diagram in RiojaPlot ---------------------------------

PRA_plot <- riojaPlot(PRA_pollen_wide[,-1:-2]*100, PRA_pollen_wide[,1:2],
          selVars = PRA_plot_taxa$variablename,
          scale.percent = TRUE,
          sec.yvar.name="age",
          plot.sec.axis = TRUE,
          xRight = 0.9)
riojaPlot(PRA_charcoal_wide[,-1:-2], PRA_charcoal_wide[,1:2], 
          minmax = data.frame(0,6),
          plot.poly = FALSE,
          plot.bar = TRUE,
          plot.line = FALSE,
          lwd.bar = 2,
          col.bar = "black",
          riojaPlot = PRA_plot)
