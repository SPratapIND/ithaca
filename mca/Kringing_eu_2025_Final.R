# https://adelieresources.com/2022/10/making-contour-maps-in-r/
rm(list = ls())
dev.off()

# Libraries
library(sf)
library(gstat)
library(tidyverse)
library(stars)
library(ggplot2)
library(viridis)
library(raster)
library(doParallel)
library(foreach)

setwd("~/ithaca/mca")

# -------- Parallel Setup ----------
library(doParallel)

ncores <- max(1, detectCores() - 1)

cl <- makeCluster(ncores)
registerDoParallel(cl)


# --- Load data ---
europe <- read.csv("eu_log_new_4_ts_box_new_final_spei.csv")

df <- europe %>% 
  filter(data_type == "Proxy") %>%
  dplyr::select(Longitude, Latitude, z_score_R) %>%
  drop_na()

# Convert to sf
df_sf <- st_as_sf(df, coords = c("Longitude", "Latitude"), crs = 4326)

# Project to UTM
df_xy <- st_transform(df_sf, 32615)

# Bounding box
bbox <- st_bbox(df_xy)
Dx <- (bbox$xmax - bbox$xmin) * 0.1
Dy <- (bbox$ymax - bbox$ymin) * 0.1
bbox_expanded <- st_bbox(c(xmin = bbox$xmin - Dx,
                           ymin = bbox$ymin - Dy,
                           xmax = bbox$xmax + Dx,
                           ymax = bbox$ymax + Dy),
                         crs = st_crs(df_xy))

# Grid
grid <- st_make_grid(df_xy, cellsize = 10000, what = "centers")

# Convert grid into list for parallel loop
grid_list <- split(grid, seq_along(grid))

# --- IDW interpolation setup ---
fit_IDW <- gstat(id = "Rain", formula = z_score_R ~ 1, data = df_xy, set = list(idp = 2))

# -------- TRUE PARALLEL IDW (this uses all cores) ----------
cat("Running parallel IDW interpolation...\n")

pred_list <- foreach(i = seq_along(grid_list),
                     .combine = rbind,
                     .packages = c("sf", "gstat")) %dopar% {
                       p <- predict(fit_IDW, grid_list[[i]])
                       as.data.frame(p)
                     }

# Convert predictions back to sf
interp_points <- do.call(rbind, lapply(seq_len(nrow(pred_list)), function(i) {
  st_as_sf(st_sfc(st_point(c(pred_list[i, "coords.x1"],
                             pred_list[i, "coords.x2"]))),
           crs = 32615) %>%
    mutate(Rain = pred_list[i, "Rain"])
}))

# Convert to stars raster
interp_IDW_stars <- st_rasterize(interp_points["Rain"], dx = 10000, dy = 10000)
interp_IDW_stars <- st_transform(interp_IDW_stars, 4326)

# Convert to dataframe
interp_df <- as.data.frame(interp_IDW_stars, xy = TRUE) %>% drop_na()

# --- Plot ---
ggplot() +
  geom_raster(data = interp_df, aes(x = x, y = y, fill = Rain)) +
  geom_contour(data = interp_df, aes(x = x, y = y, z = Rain),
               color = "black", size = 0.3) +
  geom_sf(data = df_sf, color = "black", size = 1) +
  geom_sf(data = df_sf, aes(color = z_score_R), size = 3) +
  scale_fill_gradientn(colors = rainbow(7), name = "Rain",
                       limits = c(min(df$z_score_R), max(df$z_score_R))) +
  scale_color_gradientn(colors = rainbow(7), name = "Rain") +
  coord_sf() +
  labs(title = "Final X-Y plot of IDW Interpolation (Parallel)",
       x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 14)

# -------- Stop cluster ----------
stopCluster(cl)
cat("Parallel cluster stopped.\n")


# rm(list=ls(all=TRUE))
# dev.off()
# 
# #dependencies
# library(sf)
# library(sp)
# library(rgdal)
# library(tidyr)
# library(akima)
# library(gstat)
# library(tidyverse)
# library(raster)
# library(ggplot2)
# library(viridis)
# library(stars)
# 
# #set working directory
# setwd("D:/work/MCA/p_t_MCA_cp/Europe")
# 
# # ### set cores ###
# library(doParallel)   # Parallel processing
# mc <- makeCluster(detectCores())
# registerDoParallel(mc)
# # #stopCluster(mc) ## stop cluster
# 
# ##### Data file ####
# # file <-read.csv("Europe_spatial_points_100_iso_log.csv", header = T)
# europe <-read.csv("eu_log_new_4_ts_box_new_final_spei.csv", header = T)
# 
# head(europe)
# 
# filtered_data <- europe %>%
#   filter(data_type == "Proxy")
# 
# df <- filtered_data %>%
#   dplyr::select(Longitude, Latitude, z_score_R) %>%
#   drop_na()
# 
# # A few handy crs codes
# googlecrs <- "EPSG:4326"
# webcrs <- "EPSG:3857"
# localcrs <- "EPSG:26915" # UTM 15N NAD83
# localUTM <- "EPSG:32615" # a WGS84 similar to EPSG:26915 - same UTM zone
# 
# # Locations and rainfall by lat long
# df_sf <- sf::st_as_sf(df, coords=c("Longitude", "Latitude"), crs=googlecrs, agr = "identity")
# 
# # Locations and rainfall by X, Y in webcrs
# df_web <- sf::st_transform(df_sf, crs=webcrs)
# 
# # Locations and rainfall by X, Y
# df_xy <- sf::st_transform(df_sf, crs=localcrs)
# 
# # Locations and rainfall by X, Y
# df_xy_2 <- sf::st_transform(df_sf, crs=localUTM)
# 
# # Let's define a bounding box (an AOI) based on the data (in Lat-Long)
# bbox <- sf::st_bbox(df_sf)
# 
# # Expand box by 20% to give a little extra room
# Dx <- (bbox[["xmax"]]-bbox[["xmin"]])*0.1
# Dy <- (bbox[["ymax"]]-bbox[["ymin"]])*0.1
# bbox["xmin"] <- bbox["xmin"] - Dx
# bbox["xmax"] <- bbox["xmax"] + Dx
# bbox["ymin"] <- bbox["ymin"] - Dy
# bbox["ymax"] <- bbox["ymax"] + Dy
# 
# bb <- c(bbox["xmin"], bbox["ymin"], bbox["xmax"], bbox["ymax"])
# 
# #countries shape file
# countries <- st_read("D:/work/MCA/p_t_MCA_cp/shp/ne_110m_admin_0_countries", layer="ne_110m_admin_0_countries")
# sub_Europe <- countries[countries$CONTINENT == "Europe",]
# 
# # IDW 
# fit_IDW <- gstat::gstat( 
#   formula = z_score_R ~ 1,
#   data = df_xy, 
#   #nmax = 10, nmin = 3, # can also limit the reach with these numbers
#   set = list(idp = 2) # inverse distance power
# )
# 
# bbox_xy = bbox %>%
#   sf::st_as_sfc() %>%
#   sf::st_transform(crs = localcrs) %>%
#   sf::st_bbox()
# 
# # sfc POINT file
# grd_sf <- sf::st_make_grid(x=bbox_xy,
#                            what="corners",
#                            cellsize=1000,
#                            crs=localcrs)
# 
# # We set debug.level to turn off annoying output
# interp_IDW <- predict(fit_IDW, grd_sf, debug.level=0)
# 
# # Convert to a stars object so we can use the contouring in stars
# interp_IDW_stars <- stars::st_rasterize(interp_IDW %>% dplyr::select(Rain=var1.pred, geometry))
# 
# 
# interp_IDW_stars2 <-interp_IDW_stars <- readRDS("interp_IDW_stars.rds")
# 
# 
# 
# ggplot() +
#   stars::geom_stars(data=interp_IDW_stars2) +  
#   geom_sf(data=df_xy, aes(color=Rain), size=5) +
#   geom_sf(data=df_xy, color="Black", size=1) +
#   scale_fill_gradientn(colors=rainbow(5), limits=c(0,3)) +
#   scale_color_gradientn(colors=rainbow(5), limits=c(0,3)) +
#   labs(title="Inverse Distance Weighting, Power = 2")