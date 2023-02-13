#### Charcoal workshop - neotoma2 ####
#### Göttingen 15.-16.2.2023.     ####
#### Petr Kuneš & Thomas Giesecke ####

# Install and load packages ------------------------------------------------------------

library(devtools)
devtools::install_github('NeotomaDB/neotoma2')
library(neotoma2)
library(ggplot2)
library(riojaPlot)

# Get all charcoal sites --------------------------------------------------

char_sites <- get_sites(taxa = 'Charcoal%', all_data = TRUE)
# char_sites <- get_sites(datasettype = "charcoal",all_data = TRUE)

neotoma2::summary(head(char_sites, 10))

char_ds <- char_sites %>% get_datasets(all_data = TRUE)
table(as.data.frame(datasets(char_ds))$datasettype)


char_mc <- char_ds %>% neotoma2::filter(datasettype %in% c("macrocharcoal", "charcoal", "microcharcoal"))
neotoma2::plotLeaflet(char_mc)

char_mc %>% neotoma2::filter(sitename == "Kettle Lake")

pol_ds <- get_datasets(datasettype = "pollen", all_data = TRUE)
pol_ds_char <- pol_ds %>%  filter(siteid %in% getids(char_mc)$siteid)
# pol_ds_char %>% neotoma2::filter(sitename == "Kettle Lake")
char_mc_pol <- char_mc %>% neotoma2::filter(siteid %in% getids(pol_ds_char)$siteid)

char_and_pollen <- c(pol_ds_char, char_mc_pol)

records <- char_and_pollen %>% get_downloads()
samples_ch_p <- records %>% samples()
head(records, 20)



