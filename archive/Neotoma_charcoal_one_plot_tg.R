#### Charcoal workshop - neotoma2 ####


# Install and load packages ------------------------------------------------------------

#library(devtools)
#devtools::install_github('NeotomaDB/neotoma2')
library(neotoma2)
library(riojaPlot)
library(rioja)
library(dplyr)
library(vegan)
# Get just one site ------------------------------------------------

# Select Prášilské jezero (CZ) and extract samples
PRA_ds <- get_datasets(get_sites(sitename = "Prášilské jezero"), all_data = TRUE)
plotLeaflet(PRA_ds)
datasets(PRA_ds)
PRA_pol_ch <- PRA_ds %>% neotoma2::filter(datasettype %in% c("pollen", "charcoal"))
PRA_rec <- PRA_pol_ch %>% get_downloads()
PRA_samp <- PRA_rec %>% samples()

shortSamples <- PRA_samp %>% 
  dplyr::filter(ecologicalgroup %in% c("TRSH", "UPHE")) %>%
  dplyr::filter(elementtype == "pollen") %>%
  dplyr::filter(units == "NISP")

eco.g<-shortSamples %>% dplyr::select(variablename, ecologicalgroup) %>% distinct()
 

counts <- tidyr::pivot_wider(shortSamples,
                               id_cols = c(depth, age),
                               names_from = variablename,
                               values_from = value,
                               values_fill = 0)




# Select charcoal area
PRA_charcoal <- samples(PRA_rec) %>% dplyr::filter(ecologicalgroup == "CHAR" & elementtype == "area") %>%
  dplyr::arrange(depth)


charcoal.a<-PRA_charcoal %>% dplyr::select(depth, age, value) %>% distinct()




# compute clusters of the pollen data ----------------------

p.n<-decostand(perc, "normalize")
p.dist<-vegdist(p.n, "euclidean") 
p.clust<-chclust(p.dist, method = "coniss")
plot(p.clust, hang=-1)
bstick(p.clust)
p.z<-cutree(p.clust,7)


# Plot stratigraphic diagram in RiojaPlot ---------------------------------

# getting some parameters right
perc<-counts[,-c(1,2)]/rowSums(counts[,-c(1,2)])*100

max5 <- sapply(perc, max) > 5 
max5_names <- eco.g[max5,1] 

age.depth<-counts[,c(1,2)]
colnames(age.depth)<-c("Depth", "Age in cal. BP")

ch.age<-tibble(charcoal.a$age) 
colnames(ch.age)<-"Age"
c.a<-tibble(charcoal.a$value)
colnames(c.a)<-"Charcoal area"

# in case you want to plot into a PDF use the below
#pdf("polle_char.pdf", width=20/2.54, height=12/2.54, useDingbats=FALSE)

PRA_plot <-riojaPlot(perc, age.depth[,2], 
                     selVars=max5_names, 
                     groups=eco.g,
                     scale.percent=TRUE, 
                     plot.groups=TRUE, 
                     wa.order="bottomleft",
                     plot.cumul=TRUE, 
                     plot.zones=p.z,
                     xRight = 0.9)
addRPClustZone(PRA_plot, p.clust, nZone=7)
riojaPlot(c.a, ch.age, 
          minmax = data.frame(0,6), 
          plot.bar = TRUE,  
          lwd.bar = 2, 
          plot.line = FALSE, 
          col.bar = "black",
          riojaPlot = PRA_plot)
#dev.off() 

# comparing the trend in charcoal to variability in pollen composition

# 1. getting the charcoal trend

char<-as.data.frame(charcoal.a)
p.age.depth<-as.data.frame(age.depth)
plot(char[,2], log(char[,3]+1), ylim=c(0, 6))
lines(char[,2], lowess(char[,3], f=0.1)$y, col="red")
lines(char[,2], loess(char[,3]~char[,2], span = 0.20)$fitted, col="blue")
lines(char[,2], smooth.spline(char[,3], df=25)$y, col="green")

ch.l<-loess(char[,3]~char[,2], span = 0.20)
char.p<-stats::predict(ch.l, p.age.depth[,2])

# 2. making the ordination 
rda.p<-rda(sqrt(perc))
ord<-ordisurf(rda.p, char.p, knots = 2)
plot(rda.p, type = "t")
sam.sc<-scores(rda.p, display="sites",scaling = 2, choice=1:2)
tax.sc<-scores(rda.p, display="species",scaling = 2, choice=1:2)

# 3. creating a nice ordination plot

tax.nam<-colnames(perc)
tax<-c(25, 26,28, 29, 30, 32, 34, 35)

plot(rda.p, type = "n", xlim=c(-3, 3), ylim=c(-2, 3) )
points(sam.sc[,1], sam.sc[,2], pch = p.z, lwd=2, col=p.z+6)
#lines(sam.sc[1:89,1], sam.sc[1:89,2], col="gray")
for (i in tax) {arrows(0,0,tax.sc[i,1], tax.sc[i,2], length = 0.1, angle = 20, code = 2, lwd = 2)
  text(tax.sc[i,1]*1.3, tax.sc[i,2]*1.1, tax.nam[i] , cex = 1.2)}
plot(ord, add=TRUE)
legend(-3.5, 3, legend=c("zone 1", "zone 2", "zone 3", "zone 4", "zone 5", "zone 6", "zone 7"),
       pch = c(1,2,3,4,5,6,7), pt.cex = 1.2, pt.lwd=2, bty = "n", col=c(7,8,9,10,11,12,13))


# checking the gradient length
dca<-decorana(sqrt(perc))
dca.sc<-scores(dca, display="sites")
range(dca.sc[,1])

# is there a correlation between different pollen taxa and charcoal?
p.ind<-which(char$depth %in% counts$depth)
poll.cor<-mat.or.vec(88,3)
poll.cor[,1]<-colnames(perc)
for (i in 1:88) {t<-cor.test(log(char.a[p.ind,3]+1), perc[p.ind,i])
poll.cor[i,2]<-t$p.value
poll.cor[i,3]<-t$estimate}

which(poll.cor[,2]<0.01)

write.csv(poll.cor, "poll_corr.csv")
