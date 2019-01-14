rm(list = ls())
gc()

library(gstat)
library(maps)
library(mapdata)
library(sp)
library(dplyr)
library(ggplot2)
library(ggmap)

## get data
ozone = read.csv("../data/oz96_utm.csv")

## plot map
labasin = get_stamenmap(bbox = c(-119, 33, -116, 35), zoom = 9, maptype = "toner")
ggmap(labasin) +
  geom_point(data = ozone, 
             aes(x = LON, y = LAT),
             color = "red",
             size = 2) +
  theme_void()
ggsave("../plots/data_map.png", width = 5, height = 3.2)

## Fit variogram model
la.var = variogram(MAXDAY ~ 1, locations = ~ LAT + LON ,data = ozone)
la.fit = fit.variogram(la.var, model=vgm(15, "Pow"))

dists = seq(0, max(la.var$dist), length.out = 100)
var.pts = data.frame(gamma = la.var$gamma, dist = la.var$dist)
var.fit = data.frame(pow = la.fit$psill * dists^la.fit$range, dist = dists)
  
ggplot() +
  geom_point(data = data.frame(gamma = la.var$gamma, dist = la.var$dist),
             aes(x = dist,
                 y = gamma)) +
  geom_line(data = var.fit,
            aes(x = dist,
                y = pow)) +
  theme_classic() +
  xlab("Distance") +
  ylab("Semivariance")
ggsave("../plots/variogram.png", width = 5, height = 3.2)
  
print(la.fit)


## Kriging Example
example.map = get_stamenmap(bbox = c(-118.2, 33.8, -117.8, 34.2), zoom = 11, maptype = "toner")
example.data = ozone %>%
  filter(between(LON, -118.2, -117.8),
         between(LAT, 33.8, 34.2))

ggmap(example.map) +
  # coord_cartesian() +
  geom_point(data = example.data, 
             aes(x = LON, y = LAT),
             color = "red",
             size = 4) +
  geom_segment(data = example.data, 
               aes(x = LON, y = LAT, xend = -118, yend = 34),
               color = "red") +
  geom_point(aes(x=-118, y=34), 
             colour="blue",
             size = 4) +
  geom_label(data = example.data,
             aes(x = LON, y = LAT, label = MAXDAY),
             hjust = c(0.5, -0.3, 1.3, 0.6, -0.2, -0.2, 0.3), 
             vjust = c(-0.4, 0.3, 0.2, -0.4, 0.4, -0, -0.4),
             size = 4,
             fontface = "bold") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Ozone measurements")
ggsave("../plots/example_map.png", width = 3, height = 3.2)

# use kriging from SspatialExtremes to get the weights
example.krig = SpatialExtremes::kriging(example.data$MAXDAY, 
                                        as.matrix(example.data[,c("LON", "LAT")]), 
                                        krig.coord = matrix(c(-118, 34), 1, 2),
                                        cov.mod = "powexp",
                                        sill = 14.08, 
                                        range = 0.255,
                                        smooth = 0.5)

ggmap(example.map) +
  # coord_cartesian() +
  geom_point(data = example.data, 
             aes(x = LON, y = LAT),
             color = "red",
             size = 4) +
  geom_segment(data = example.data, 
               aes(x = LON, y = LAT, xend = -118, yend = 34),
               color = "red") +
  geom_point(aes(x=-118, y=34), 
             colour="blue",
             size = 4) +
  geom_label(data = cbind(example.data, Weight = round(example.krig$weights, 2)),
            aes(x = LON, y = LAT, label = Weight),
            hjust = c(0.5, -0.3, 1.3, 0.6, -0.2, -0.2, 0.3), 
            vjust = c(-0.4, 0.3, 0.2, -0.4, 0.4, -0, -0.4),
            size = 4,
            fontface = "bold") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Kriging Weights")
ggsave("../plots/weights_map.png", width = 3, height = 3.2)


ggmap(example.map) +
  # coord_cartesian() +
  geom_point(data = example.data, 
             aes(x = LON, y = LAT),
             color = "red",
             size = 4) +
  geom_segment(data = example.data, 
               aes(x = LON, y = LAT, xend = -118, yend = 34),
               color = "red") +
  geom_point(aes(x=-118, y=34), 
             colour="blue",
             size = 4) +
  geom_label(aes(x=-118, y=34, label = round(example.krig$krig.est, 1)),
             hjust = -0.3, 
             vjust = 0.6,
             size = 4,
             fontface = "bold") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Prediction")
ggsave("../plots/example_pred.png", width = 3, height = 3.2)

# Does everything add up? (these should be the same)
sum(example.krig$weights * example.data$MAXDAY)
example.krig$krig.est

## Estinate interpolation points (Kriging)
lats = seq(33, 35, length.out = 200)
lons = seq(-119, -116, length.out = 200)
la.grid  = expand.grid(lats, lons)

coordinates(ozone) = ~ LAT + LON
coordinates(la.grid) = ~ Var1 + Var2

la.krig = krige(MAXDAY ~ 1, ozone, la.grid, model = la.fit)

ozone = as.data.frame(ozone)
la.krig = as.data.frame(la.krig)


## Make Ozone into levels
la.krig[["Ozone"]] = cut(la.krig[["var1.pred"]], 10, right = F)
levels(la.krig[["Ozone"]]) = sapply(levels(la.krig[["Ozone"]]), function(x) {
  gsub(",", " - ", substr(x, 2, nchar(x)-1))
})

## Make Uncertainty into levels
la.krig[["Variance"]] = cut(la.krig[["var1.var"]], 10, right = F)
levels(la.krig[["Variance"]]) = sapply(levels(la.krig[["Variance"]]), function(x) {
  gsub(",", " - ", substr(x, 2, nchar(x)-1))
})

## plot predictions
ggmap(labasin) +
  geom_raster(data = la.krig, 
              aes(Var2, Var1, fill = Ozone), 
              alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_brewer(palette = "Spectral") +
  geom_point(data = ozone, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 2) +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill="Ozone predictions")
ggsave("../plots/prediction_map.png", width = 5, height = 3.2)

## plot uncertainty
ozone = as.data.frame(ozone)

ggmap(labasin) +
  geom_raster(data = la.krig, 
              aes(Var2, Var1, fill = Variance), 
              alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_brewer(palette = "Spectral") +
  geom_point(data = ozone, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 2) +
  theme_void() +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill="Prediction variance")
ggsave("../plots/uncertainty_map.png", width = 5, height = 3.2)


