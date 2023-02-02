#### Charcoal workshop - neotoma2 ####
#### Simon Goring charcoal 

# Download charcoal ------------------------------------------------------------


library(neotoma2)

char_sites <- get_sites(taxa = 'Charcoal%', all_data = TRUE)
char_sites %>% neotoma2::filter(sitename == "Kettle Lake")

char_ds <- char_sites %>% get_datasets(all_data = TRUE)
# char_ds %>% neotoma2::filter(sitename == "Kettle Lake")
table(as.data.frame(datasets(char_ds))$datasettype)

char_mc <- char_ds %>% neotoma2::filter(datasettype %in% c("macrocharcoal", "charcoal", "microcharcoal"))
# char_mc %>% neotoma2::filter(sitename == "Kettle Lake")

char_pol <- char_ds %>% filter(datasettype == "pollen" & siteid %in% getids(char_mc)$siteid)
# char_pol <- char_ds %>% filter(datasettype %in% "pollen")
char_pol %>% neotoma2::filter(sitename == "Kettle Lake")

summary(char_pol)

char_mc_small <- char_mc %>% filter(siteid %in% getids(char_pol)$siteid)
summary(char_mc_small)

char_and_pollen <- c(char_pol, char_mc_small)
summary(char_and_pollen)

records <- char_and_pollen %>% get_downloads()

# To just get one site:
pop_samp <- records %>% filter(sitename == "Popradské  pleso") %>% samples()

library(ggplot2)
ggplot(data = pop_samp %>% dplyr::filter(variablename %in% c("Picea", "Charcoal"))) +
  geom_path(aes(x = value, y = ageyounger, color = variablename))



