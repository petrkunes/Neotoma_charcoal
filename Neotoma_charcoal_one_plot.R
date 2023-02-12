#### Charcoal workshop - neotoma2 ####


# Install and load packages ------------------------------------------------------------

library(devtools)
devtools::install_github('NeotomaDB/neotoma2')
library(neotoma2)
library(rioja)
library(riojaPlot)
library(vegan)

# Get just one site ------------------------------------------------

# Select Prášilské jezero (CZ) 
PRA_ds <- get_datasets(get_sites(sitename = "Prášilské jezero"), all_data = TRUE)
plotLeaflet(PRA_ds)

# Filter our the datasets from pollen and charcoal
PRA_pol_ch <- PRA_ds %>% neotoma2::filter(datasettype %in% c("pollen", "charcoal"))

# Pull in and extract sample data
PRA_rec <- PRA_pol_ch %>% get_downloads()
PRA_samp <- PRA_rec %>% samples()

# Select charcoal area
PRA_charcoal <- samples(PRA_rec) %>% dplyr::filter(ecologicalgroup == "CHAR" & elementtype == "area") %>%
  dplyr::select(depth, age, variablename, elementtype, value) %>%
  dplyr::mutate(variable = paste(variablename, elementtype, sep = " "), .keep = "unused") %>%
  dplyr::group_by(depth, age) %>%
  dplyr::arrange(depth)

# Select only pollen samples measured using NISP and taxa for calculating %.
PRA_samp_short <- samples(PRA_rec) %>% 
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

# Select only pollen taxa (and only trees/shrubs) for plotting
PRA_plot_taxa <- taxa(PRA_rec) %>% dplyr::filter(ecologicalgroup %in% c("TRSH")) %>%
  dplyr::filter(elementtype == "pollen") %>%
  dplyr::arrange(desc(samples)) %>% 
  head(n = 15) 

# Calculate cluster analysis of the pollen data
PRA_clust <- chclust(vegdist(PRA_pollen_wide[,-1:-2], "euclidean"), method = "coniss")
plot(PRA_clust, hang = -1)
bstick(PRA_clust)
PRA_zones<-cutree(PRA_clust,7)

# Plot
PRA_plot1 <- riojaPlot(PRA_pollen_wide[,-1:-2]*100, PRA_pollen_wide[,1:2],
          selVars = PRA_plot_taxa$variablename,
          scale.percent = TRUE,
          sec.yvar.name="age",
          plot.sec.axis = TRUE,
          srt.xlabel = 60,
          xRight = 0.8)
PRA_plot2 <- addRPClust(PRA_plot1, PRA_clust, xLeft=0.8, xRight = 0.9)
addRPClustZone(PRA_plot2, PRA_clust, nZone=7, xRight = 0.9, col = "red")

PRA_plot3 <- riojaPlot(PRA_charcoal_wide[,-1:-2], PRA_charcoal_wide[,1:2], 
          minmax = data.frame(0,6),
          plot.poly = FALSE,
          plot.bar = TRUE,
          plot.line = FALSE,
          lwd.bar = 2,
          col.bar = "black",
          xLeft = 0.9,
          riojaPlot = PRA_plot2)

# Comparing the trend in charcoal to variability in pollen composition --------

# getting the charcoal trend

plot(PRA_charcoal_wide$age, log(PRA_charcoal_wide$`Charcoal area`+1), ylim=c(0, 6), xlab = "Age cal BP", ylab = "CHAR")
lines(PRA_charcoal_wide$age, lowess(PRA_charcoal_wide$`Charcoal area`, f=0.1)$y, col="red")
lines(PRA_charcoal_wide$age, loess(PRA_charcoal_wide$`Charcoal area`~PRA_charcoal_wide$age, span = 0.2)$fitted, col="blue")
lines(PRA_charcoal_wide$age, smooth.spline(PRA_charcoal_wide$`Charcoal area`, df = 25)$y, col="green")
legend(1,6, legend = c("LOWESS", "LOESS", "smooth spline"), col = c("red", "blue", "green"), lty=1, cex=0.8)

# estimate values for given ages
PRA_char_loess <- loess(PRA_charcoal_wide$`Charcoal area`~PRA_charcoal_wide$age, span = 0.20)
PRA_char_pred <- stats::predict(PRA_char_loess, PRA_pollen_wide$age)

# making the ordination 

PRA_pca <- rda(decostand(PRA_pollen_wide[,-1:-2], method = "hellinger"))
PRA_PCA_plot <- biplot(PRA_pca)
plot(envfit(PRA_pca, data.frame(CHAR = PRA_char_pred)))
orditorp(PRA_PCA_plot, display = "species")

# Fit smooth sutface of charcoal to ordination plot

PRA_ordisurf <- ordisurf(PRA_pca, PRA_char_pred, knots = 2)

# creating a nice ordination plot

tax.sc<-scores(PRA_pca, display="species", choice=1:2)
tax<-c("Pinus", "Betula", "Picea", "Fagus", "Poaceae", "Corylus")

ordiplot(PRA_pca, type = "n")
points(PRA_pca, col = PRA_zones + 6, pch = PRA_zones, lwd = 2)
for (i in tax) {arrows(0,0,tax.sc[i,1], tax.sc[i,2], length = 0.1, angle = 20, code = 2, lwd = 2, col = "red")
  text(tax.sc[i,1]+0.05, tax.sc[i,2]+0.05, i , cex = 1.2, col = "red")
  }
plot(PRA_ordisurf, add=TRUE)
legend(-1, 0, legend=c("zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6", "zone 7"),
       pch = c(1,2,3,4,5,6,7), pt.cex = 1.2, pt.lwd=2, bty = "n", col=c(7,8,9,10,11,12,13))

# checking the gradient length with DCA
PRA_dca <- decorana(decostand(PRA_pollen_wide[,-1:-2], method = "hellinger"))
PRA_dca
plot(PRA_dca, display = "sites")

# Check correlation between different pollen taxa and charcoal?
install.packages("corrplot")
library(corrplot)
install.packages("Hmisc")
library(Hmisc)

PRA_cor_data <- merge(PRA_pollen_wide, PRA_charcoal_wide)
PRA_cor_data$`Charcoal area` <- log(PRA_cor_data$`Charcoal area` + 1)
PRA_cor <- rcorr(as.matrix(PRA_cor_data[,-1:-2]), type = "pearson")

cor.data <- as.matrix(PRA_cor$r[89,which(PRA_cor$P[89,-89]<0.05)])
colnames(cor.data) <- "Charcoal"

corrplot(t(cor.data))
