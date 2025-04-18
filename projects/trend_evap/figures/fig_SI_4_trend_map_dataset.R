# Figure to support figure 4 - trend by dataset ----
source('source/evap_trend.R')

library(rnaturalearth)

# Map preparation -----
## World and Land borders ----
PATH_SAVE_PARTITION_EVAP <- paste0(PATH_SAVE, "partition_evap/")
load(paste0(PATH_SAVE_PARTITION_EVAP, "paths.Rdata"))

earth_box <- readRDS(paste0(PATH_SAVE_PARTITION_EVAP_SPATIAL,
                            "earth_box.rds")) %>%
  st_as_sf(crs = "+proj=longlat +datum=WGS84 +no_defs")
world_sf <- ne_countries(returnclass = "sf")

## Labels ----
labs_y <- data.frame(lon = -165, lat = c(50, 25, -5, -35, -65))
labs_y_labels <- seq(60, -60, -30)
labs_y$label <- ifelse(labs_y_labels == 0, "°", ifelse(labs_y_labels > 0, "°N", "°S"))
labs_y$label <- paste0(abs(labs_y_labels), labs_y$label)
labs_y <- st_as_sf(labs_y, coords = c("lon", "lat"),
                   crs = "+proj=longlat +datum=WGS84 +no_defs")

labs_x <- data.frame(lon = seq(120, -120, -60), lat = -82)
labs_x$label <- ifelse(labs_x$lon == 0, "°", ifelse(labs_x$lon > 0, "°E", "°W"))
labs_x$label <- paste0(abs(labs_x$lon), labs_x$label)
labs_x <- st_as_sf(labs_x, coords = c("lon", "lat"),
                   crs = "+proj=longlat +datum=WGS84 +no_defs")

## read data ----
evap_trend <- readRDS(paste0(PATH_SAVE_EVAP_TREND, "global_grid_per_dataset_evap_slope_bootstrap.rds"))  
evap_trend[, dataset := toupper(dataset)]
evap_trend[dataset == "ETMONITOR", dataset := "ETMonitor"]
evap_trend[dataset == "SYNTHESIZEDET", dataset := "SynthesizedET"]
evap_trend[dataset == "ERA5-LAND", dataset := "ERA5-land"]
evap_trend[dataset == "MERRA2", dataset := "MERRA-2"]
evap_trend[dataset == "JRA55", dataset := "JRA-55"]
evap_trend[dataset == "TERRACLIMATE", dataset := "TerraClimate"]

evap_trend[p < 1e-6, p := 1e-6]

evap_trend[, pval_brk := cut(p, breaks = c(0, 0.01, 0.05, 0.1, 0.2, 1))]

cols_trend <- c("(-117,-5]" = "darkblue","(-5,-3]" = "royalblue3",  "(-3,-1]" = "royalblue1",  
                "(-1,0]" = "lightblue", "(0,1]" = "orange", "(1,3]" = "lightcoral" , "(3,5]" = "darkred", "(5,142]" = "#330000")


dataset_list <- evap_trend[, unique(dataset)]

for(dataset_sel in dataset_list){

to_plot_sf <- evap_trend[dataset == dataset_sel, .(lon, lat, slope)]
to_plot_sf <- to_plot_sf[, .(lon, lat, slope)] %>% 
  rasterFromXYZ(res = c(0.25, 0.25),
                crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_as_stars() %>% st_as_sf()

to_plot_sf$slope_brk <-  cut(to_plot_sf$slope, breaks = c(-117,-5,-3,-1, 0, 1, 3, 5, 142))


fig_slope <- ggplot(to_plot_sf) +
  geom_sf(data = world_sf, fill = "light gray", color = "light gray") +
  geom_sf(aes(color = as.factor(slope_brk), fill = as.factor(slope_brk))) +
  geom_sf(data = earth_box, fill = NA, color = "black", lwd = 0.1) +
  scale_fill_manual(values = cols_trend) +
  scale_color_manual(values = cols_trend,
                     guide = "none") +
  labs(x = NULL, y = NULL, fill = expression(paste("ET trend \n[mm year"^-~2,"]   "))) +
  coord_sf(expand = FALSE, crs = "+proj=robin") +
  scale_y_continuous(breaks = seq(-60, 60, 30)) +
  geom_sf_text(data = labs_y, aes(label = label), color = "gray20", size = 4) +
  geom_sf_text(data = labs_x, aes(label = label), color = "gray20", size = 4) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA), panel.ontop = TRUE,
        panel.border = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        panel.grid.major = element_line(colour = "gray60"),
        axis.text = element_blank(), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),        
        legend.spacing.x = unit(1, "cm"),
        legend.spacing.y = unit(1, "cm"),
        plot.title = element_text(size = 20))+
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))


to_plot_sf <- evap_trend[dataset == dataset_sel, .(lon, lat, lower)]
to_plot_sf <- to_plot_sf[, .(lon, lat, lower)] %>% 
  rasterFromXYZ(res = c(0.25, 0.25),
                crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_as_stars() %>% st_as_sf()

to_plot_sf$lower_brk <-  cut(to_plot_sf$lower, breaks = c(-117,-5,-3,-1, 0, 1, 3, 5, 142))

to_plot_sf <- evap_trend[dataset == dataset_sel, .(lon, lat, pval_brk)
][, value := as.numeric(pval_brk)]
to_plot_sf <- to_plot_sf[, .(lon, lat, value)] %>% 
  rasterFromXYZ(res = c(0.25, 0.25),
                crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
  st_as_stars() %>% st_as_sf()

cols_pval <- c("darkblue","darkorchid4","darkorchid1","lightcoral","gold")

fig_pval <- ggplot(to_plot_sf) +
  geom_sf(data = world_sf, fill = "light gray", color = "light gray") +
  geom_sf(aes(color = as.factor(value), fill = as.factor(value))) +
  geom_sf(data = earth_box, fill = NA, color = "black", lwd = 0.1) +
  scale_fill_manual(values = cols_pval, labels = levels(evap_trend$pval_brk)) +
  scale_color_manual(values = cols_pval,
                     guide = "none") +
  labs(x = NULL, y = NULL, fill = expression(paste("P-value"))) +
  coord_sf(expand = FALSE, crs = "+proj=robin") +
  scale_y_continuous(breaks = seq(-60, 60, 30)) +
  geom_sf_text(data = labs_y, aes(label = label), color = "gray20", size = 4) +
  geom_sf_text(data = labs_x, aes(label = label), color = "gray20", size = 4) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA), panel.ontop = TRUE,
        panel.border = element_blank(),
        axis.ticks.length = unit(0, "cm"),
        panel.grid.major = element_line(colour = "gray60"),
        axis.text = element_blank(), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16),        
        legend.spacing.x = unit(1, "cm"),
        legend.spacing.y = unit(1, "cm"),
        plot.title = element_text(size = 20))+
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))


fig_col <- ggarrange(fig_slope, fig_pval,  labels = c("a", "b"), ncol = 1, nrow = 2, 
                     font.label = list(size = 20))

annotate_figure(fig_col, top = text_grob(dataset_sel, 
                                      color = "black", size = 28))

ggsave(paste0(PATH_SAVE_EVAP_TREND_FIGURES_SUPP, "fig3_SI_maps_trends_",dataset_sel,".png"), 
       width = 8, height = 8)

}
