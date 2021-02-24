# plot make map for figure 1


# Libraries ---------------------------------------------------------------

library(sf)
library(ggspatial)
library(tidyverse)
library(mapview)
library(ggthemes)
library(ggrepel)


# Data --------------------------------------------------------------------

load(file = "data_output/fig_01_map_nhd2.rda")

sites <- read_csv("data_output/Table_01_samplesites_samples.csv") %>% st_as_sf(coords=c("lon", "lat"), crs=4326, remove=FALSE)

# dams
dams <- read_sf("data/shps/Selected_CA_Dams/Selected_CA_Dams.shp") %>% 
   st_transform(st_crs(h8))

dams <- st_intersection(dams, h8) %>% 
   filter(!NID %in% c("CA00351", "CA00369", "CA00347", "CA00540",
                      "CA00253", "CA00246", "CA00964", "CA10105", 
                      "CA00820", "CA00818", "CA00819", "CA10306", 
                      "CA00546"))

# add slate creek dam
slateck <- tibble(OBJECTID_1 = 1, NAME="Slate", RIVER="SLate Creek", lon=-121.05142, lat=39.61557) %>% 
   st_as_sf(coords=c("lon", "lat"), crs=4326)

# bind w dams 
dams <- bind_rows(dams, slateck)

# filter ends from river stretches:
nhdnfy <- nhdnfy %>% filter(hydroseq > 10011419)
nhdslate <- riv_clip %>% 
   filter(hydroseq %in% c(10017870,10018175,10018486,10018814, 10019178))   
nhdnfa <- nhdnfa %>% filter(hydroseq > 10009770)
nhdrub <- nhdrub %>% filter(hydroseq < 10018430)
nhdorck <- riv_clip %>% filter(hydroseq %in% c(10030745, 10029489, 10028376, 
                                               10027377, 10026457, 10025623))
nhdmfa_peaking <- nhdmfa %>% filter(hydroseq < 10011305)
nhdmfa_bypass <- nhdmfa %>% filter(hydroseq <= 10021711, hydroseq >=10011305)

# Quick Mapview -----------------------------------------------------------

# quick map
# mapview(nhdmfa_bypass, color="orange") + 
#    mapview(nhdmfa_peaking, color="red3") + mapview(nhdrub, color="orange2") + 
#    mapview(nhdnfa, color="steelblue") + mapview(dams, col.regions="black") +
#    mapview(nhdmfy, color="orange") + mapview(nhdorck, color="orange") +
#    mapview(nhdbear, color="orange") + mapview(nhdnfy, color="steelblue") + 
#    mapview(nhdslate, color="orange") + mapview(nhdsfy, color="orange")+
#    mapview(sites, col.regions="white") + mapview(riv_clip, color="steelblue", lwd=0.8)

# Main Map No baselayer --------------------------------------------------

## ggplot version
# (fig1 <- ggplot() + 
#     #annotation_map_tile(type = "osm", zoomin = -1) +
#     coord_sf(label_axes = "--EN", datum = 4326) +
#     geom_sf(data=riv_clip, col="steelblue", lwd=0.25, alpha=0.75) +
#     geom_sf(data=h8, aes(fill=HU_8_NAME), color="gray80", lwd=0.7, alpha=0.5) + 
#     scale_fill_grey("HUC8 Watershed", guide=guide_legend(order=1)) +
#     # add mainstems:
#     geom_sf(data=nhdnfa, aes(color="Unregulated"), lwd=2) +
#     geom_sf(data=nhdnfy, aes(color="Unregulated"), lwd=2) +
#     geom_sf(data=nhdmfa_peaking, aes(color="Hydropeaking"), lwd=2) +
#     geom_sf(data=nhdmfa_bypass, aes(color="Bypass"), lwd=2) +
#     geom_sf(data=nhdsfy, aes(color="Bypass"), lwd=2) +
#     geom_sf(data=nhdrub, aes(color="Bypass"), lwd=2) +
#     geom_sf(data=nhdslate, aes(color="Bypass"), lwd=2) +
#     geom_sf(data=nhdorck, aes(color="Bypass"), lwd=2) +
#     geom_sf(data=nhdmfy, aes(color="Bypass"), lwd=2) +
#     geom_sf(data=nhdbear, aes(color="Bypass"), lwd=2) +
#     
#     scale_color_manual("", values = c("Unregulated"= "steelblue", "Bypass"="orange", "Hydropeaking"="red3"), guide=guide_legend(order=2)) +
#     # add points:
#     geom_sf(data=dams, aes(shape="Dams"), fill="black", color="black", size=3)+
#     # add labels
#     geom_text_repel(data = sites,
#                     aes(x=lon, y=lat, label=SiteID),
#                     fontface="bold", family="Roboto Condensed",
#                     segment.alpha = .8,
#                     segment.color="gray10",
#                     color = "white",     # text color
#                     bg.color = "grey30", # shadow color
#                     bg.r = 0.15,
#                     size=4.5, 
#                     box.padding = 0.6,
#                     force=0.7,
#                     min.segment.length = .25,
#                     segment.inflect = FALSE, segment.square=FALSE,
#                     segment.ncp = 3,
#                     segment.curvature = -0.1) +
#     geom_sf(data=sites, aes(shape="Sample Sites"), color="black", fill="white", size=3, lwd=1.5) +
#     scale_shape_manual("", values=c("Dams"=24, "Sample Sites"=21), guide = guide_legend(override.aes = list(alpha=c(1,1), pch=c(24, 21), color=c("black", "black"), fill=c("black", "white")), order=3)) +
#     labs(y="", x="") +
#     #coord_sf(datum = NA) +
#     # theme
#     #hrbrthemes::theme_ipsum_rc() +
#     theme_minimal(base_family = "Roboto Condensed", base_size = 14) +
#     # spatial-aware automagic scale bar
#     annotation_scale(location = "br", style="bar") +
#     # spatial-aware automagic north arrow
#     annotation_north_arrow(width = unit(0.8,"cm"), 
#                            height=unit(1.2, "cm"),
#                            pad_y = unit(1.2, "cm"),
#                            location = "br", 
#                            which_north = "true") +
#     theme(legend.background = element_rect(fill = NA, color = NA),
#           # change grid lines to "transparent" or alpha them down
#           panel.grid.major = element_line(colour = alpha("gray90", .2)),
#           legend.key = element_rect(fill = NA, color = NA),
#           legend.position = c(0.92, 0.82), 
#           legend.spacing.y = unit(0.1,"cm"),
#           legend.margin = margin(0, 0, 0, 0.1, "cm"))
# )


#ggsave(filename="figs/fig_01A_study_area_map_final.pdf", width=11, height = 8, units = "in", dpi=600, device=cairo_pdf)
#ggsave(filename="figs/fig_01_study_area_map_NHD2_wide.png", width=8, height = 7, units = "in", dpi=300)

# Get GGMAP Baselayer ------------------------------------------------

library(ggmap)

# make a bounding box in lat/lon
mapRange1 <- riv_clip

(mapRange1 <- c(range(st_coordinates(mapRange1)[,1]),range(st_coordinates(mapRange1)[,2])))

# # this only works with an API KEY
# map_bw <- get_map(location=c(mapRange1[1], mapRange1[3],mapRange1[2], mapRange1[4]), 
#                   crop = FALSE, force = TRUE, 
#                   color="bw",
#                   maptype="terrain-background", # terrain/google terrain-background-stamen
#                   source="stamen", # osm #stamen
#                   zoom=10)
# #ggmap(map_bw)
# 
# map_col <- get_map(location=c(mapRange1[1], mapRange1[3],mapRange1[2], mapRange1[4]), 
#                   crop = FALSE, force = TRUE,
#                   color="color",
#                   maptype="terrain", # terrain/google terrain-background-stamen
#                   source="google", # osm #stamen
#                   zoom=10)
# 
# 
# #ggmap(map_col)
# 
# # save as an object for later
# ggmap_bw_z10 <- map_bw
# save(ggmap_bw_z10, file = "data_output/ggmap_base_layer_zoom10.rda")
# 
# ggmap_color_z10 <- map_col
# save(ggmap_color_z10, file = "data_output/ggmap_base_layer_color_zoom10.rda")

# Map with Baselayer ------------------------------------------------------

load("data_output/ggmap_base_layer_color_zoom10.rda")
load("data_output/ggmap_base_layer_zoom10.rda")

fig1_col <- ggmap(ggmap_color_z10) +
    # coordinate system
    coord_sf(label_axes = "--EN", datum = 4326) +
    # hucs
    geom_sf(data=h8, aes(fill=HU_8_NAME), color="gray80", size=0.7, alpha=0.5, inherit.aes = FALSE) +
    scale_fill_grey("HUC8 Watershed", 
                    guide=guide_legend(override.aes = list(order=0, color=NA))) +
    # rivers base
    geom_sf(data=riv_clip, col="dodgerblue", size=0.3, alpha=.9, inherit.aes = FALSE) +
    # river mainstems
    geom_sf(data=nhdnfa, aes(color="Unregulated"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdnfy, aes(color="Unregulated"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdmfa_peaking, aes(color="Hydropeaking"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdmfa_bypass, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdsfy, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdrub, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdslate, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdorck, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdmfy, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    geom_sf(data=nhdbear, aes(color="Bypass"), size=2, inherit.aes = FALSE) +
    # fix color scale for mainstems
    scale_color_manual("", values = c("Unregulated"= "steelblue", "Bypass"="orange", "Hydropeaking"="red3"), 
                       guide=guide_legend(override.aes=list(order=1))) +
    # add dams
    geom_sf(data=dams, aes(shape="Dams"), fill="black", color="black", size=3, inherit.aes = FALSE)+
    # add sample sites: labels
    geom_text_repel(data = sites,
                    aes(x=lon, y=lat, label=SiteID),
                    fontface="bold", family="Helvetica",
                    segment.alpha = .8,
                    segment.color="grey95",
                    color = "white",     # text color
                    bg.color = "black", # shadow color
                    bg.r = 0.15,
                    size=4.5,
                    box.padding = 0.5,
                    force=1.5,
                    min.segment.length = .25,
                    segment.inflect = FALSE, segment.square=FALSE,
                    segment.ncp = 3,
                    segment.curvature = -0.1, inherit.aes = FALSE) +
    # add sample sites: points
    geom_sf(data=sites, aes(shape="Sample Sites"), color="black", fill="white", size=3, inherit.aes = FALSE) +
    # fix legened
    scale_shape_manual("", values=c("Dams"=24, "Sample Sites"=21), 
                       guide = guide_legend(override.aes = list(alpha=c(1,1), 
                                                                pch=c(24, 21), 
                                                                color=c("black", "black"), 
                                                                fill=c("black", "white")), 
                                            order=99)) +
    labs(y="", x="") +
    theme_minimal(base_family = "Helvetica", base_size = 9) +
    annotation_scale(location = "br", style="bar") +
    annotation_north_arrow(width = unit(0.8,"cm"),
                           height=unit(1.2, "cm"),
                           pad_y = unit(1.2, "cm"),
                           location = "br",
                           which_north = "true") +
    theme(legend.box.background = element_rect(fill = "white", 
                                               color = "gray10"),
          legend.box.margin = margin(0.2, 0.2, 0.1, 0.1, "cm"),
          # change grid lines to "transparent" or alpha them down
          panel.grid.major = element_line(colour = alpha("gray90", .2)),
          legend.key.size = unit(0.5, 'cm'),
          legend.title = element_text(vjust=1.5),
          axis.text = element_text(color="black", size=10),
          legend.position = c(0.86, 0.78),
          legend.spacing.y = unit(-0.05,"cm"),
          legend.margin = margin(0, 0, 0, 0.1, "cm"))

# plot
fig1_col

#ggsave(filename="figs/fig_01_study_area_map_google_color_z10.pdf", width=11, height = 8, units = "in", dpi=300, device=cairo_pdf)


# ADD INSET ---------------------------------------------------------------

us <- USAboundaries::us_boundaries(type="state", resolution = "low") %>% 
    filter(!state_abbr %in% c("PR", "AK", "HI"))

# make a box around rivers (a grid with an n=1) for inset
ca_box <- st_make_grid(riv_clip, n = 1) %>% st_centroid()

# Inset map: US
p2 <- ggplot() + 
    geom_sf(data = us, colour = "grey10", fill = "tan", alpha=0.4) +
    coord_sf() +
    theme_minimal() + 
    geom_sf(data=ca_box, fill="black", color="white", pch=23, size=5) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_line(colour = "transparent"),
          plot.background = element_rect(color = "black", fill="white"),
          #plot.background = element_blank(),
          panel.border = element_blank(),
          #panel.background = element_rect(fill="white"),
          plot.margin = unit(c(0, 0, 0 ,0), "mm"))
p2

# make with rounded edges
library(grid)
g <- ggplotGrob(p2)
bg <- g$grobs[[1]]
round_bg <- roundrectGrob(x=bg$x, y=bg$y, width=bg$width, height=bg$height,
                          r=unit(0.1, "snpc"),
                          just=bg$just, name=bg$name, gp=bg$gp, vp=bg$vp)
g$grobs[[1]] <- round_bg
cowplot::plot_grid(g)


# COMBINED MAP: PATCHWORK -------------------------------------------------

library(patchwork)

fig1_col + 
    inset_element(g, left = 0.12, bottom = 0.77, right = 0.38, top = 1.05, 
                  align_to = 'plot') +
    theme(plot.background = element_blank())

ggsave(filename="figs/fig_01_study_area_map_inset_col_final_patch.pdf", 
       width=15.6, height = 16, units = "cm", dpi=600, scale = 1)

 # FINAL COMBINED MAP ------------------------------------------------------

# add diff libraries
# library(grid)
# library(gridExtra)
# 
# 
# cairo_pdf(filename = "figs/fig_01_study_area_map_inset_col.pdf", width = 11, height = 8, family = "Helvetica", pointsize = 12, fallback_resolution = 300)
# 
# #png(filename = "figs/fig_01A_study_area_map_inset_col.png", width = 11, height = 8, family = "Helvetica", res = 300, pointsize = 12, units = "in")
# 
# grid.newpage()
# mainmap <- viewport(width = 1, height = 1, x = 0.5, y = 0.5) # main map
# insetmap <- viewport(width = 0.2, height = 0.26, x = 0.31, y = 0.89) # inset
# print(fig1_col, vp = mainmap) 
# print(g, vp = insetmap)
# 
# dev.off()
