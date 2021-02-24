# COMBINE THE WAVELET AND COLWELL DATA AND PLOT

# LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(viridis)
library(ggrepel)
library(ggthemes)
#library(sf)


# GET DATA AND CLEAN ------------------------------------------------------

load("models/wavelet_colwell_combined.rda")

# filter to study rivers only and 
# fix error in "MFA_11427765_1973_2016" (actually RUBICON)

combined_df <- combined_df %>% filter(!river %in% c("SFA", "NFMFA", "ORCK")) %>% 
  mutate(site = if_else(site=="MFA_11427765_1973_2016", "RUB_11427765_1973_2016", site),
         river = if_else(site=="RUB_11427765_1973_2016", "RUB", river)) %>% 
  mutate(regtype=case_when(
    river %in% c("RUB","SFY", "MFY", "BEAR") ~ "Bypass",
    river %in% c("NFY", "NFA") ~ "Unregulated",
    river %in% c("MFA") ~ "Hydropeaking"
  )) 

combined_df$regtype <- as.factor(combined_df$regtype)

## DROP SOME SITES TO CLEAN UP PICTURE:
# gages_to_drop <- c(11422500, 11421900, 11414196, 11421790,11433212, 11433300, 11413100, 11410500)

# Keep non-duplicate gages (like Bear Sites which have less gaging)
gages_to_drop <- c(11421900, 11414196, 11433212, 11433300, 11413100)

df <- combined_df %>% filter(!gageID %in% gages_to_drop)

# PLOT BY END YEAR ---------------------------------------------------

# plot color by end year and reg/unreg

# seasonality vs. predictability
# ggplot() +
#   geom_point(data=combined_df, aes(x=MP, y=maxPowerAvg, shape=reg, fill=as.integer(year_end)), 
#              pch=ifelse(combined_df$reg=="Y", 24, 21), size=7.5, show.legend = T)+
#   geom_label_repel(data=combined_df, aes(x=MP, y=maxPowerAvg, label=paste0(river, ":",gageID)), nudge_x = 0.07, size=2)+
#   scale_fill_viridis_c("End Year (Flow)", direction = -1, option = "A")+
#   theme_classic(base_family = "Roboto Condensed", base_size = 14) + 
#   xlab("Seasonality (M/P)") + ylab("Predictability (Avg. Power at 12 Months)")+
#   theme(legend.position=c(0.2,.8),
#         legend.text = element_text(size=14)) +
#   geom_vline(xintercept = 0.75, col="gray", lty=2, lwd=1.5) +
#   annotate("text", x=0.70, y=25, label="Post-Dam") +
#   annotate("text", x=0.79, y=25, label="Pre-Dam")

#ggsave("./figs/seasonality_predictability_pre-post_dots.png", width = 10, height=8, units = "in", dpi = 300)

# FIG 2: COMBINED COWPLOT -------------------------------------------

# (p1 <- ggplot() +
#   geom_text_repel(data=df %>% filter(reg=="N"), 
#                   aes(x=MP, y=maxPowerAvg, label=glue::glue("{river}: {gageID} ({year_end})")),
#                   fontface="bold", family="Roboto Condensed",
#                   segment.alpha = .8,
#                   segment.color="grey10",
#                   color = "black",     # text color
#                   #bg.color = "grey30", # shadow color
#                   #bg.r = 0.15, # only with geom_text
#                   size=2.5,
#                   force = .5,
#                   segment.curvature = -0.1,
#                   segment.ncp = 3,
#                   segment.angle = 20,
#                   box.padding = 0.7, 
#                   point.padding = 0.7,
#                   nudge_y = 0.25,
#                   nudge_x = 0.15,
#                   min.segment.length = .2)+
#   geom_point(data=df %>% filter(reg=="N"), 
#              aes(x=MP, y=maxPowerAvg, 
#                  shape=regtype), fill="#56B4E9", size=7, show.legend = F, alpha=0.1)+
#   geom_point(data=df %>% filter(reg=="N"), 
#              aes(x=MP, y=maxPowerAvg, 
#                  shape=regtype), size=7, show.legend = F, alpha=1) +
#   geom_point(data=df %>% filter(gageID %in% c("11413000", "11427000")), 
#              aes(x=MP, y=maxPowerAvg, 
#                  shape=regtype), size=7, fill="#56B4E9", show.legend = F, alpha=1) +
#   #scale_fill_colorblind("Flow Regime")+
#   scale_shape_manual("Flow Regime", values=c("Bypass"=22, "Hydropeaking"=24, "Unregulated"=21))+
#   xlim(c(0, 1)) + 
#   scale_y_continuous(limits = c(0,18), breaks = c(seq(0,18,3)))+
#   theme_classic(base_family = "Roboto Condensed", base_size = 12) + 
#     #xlab("Seasonality (M/P)") + 
#     xlab("") +
#     ylab("Predictability")+
#   theme(legend.position=c(0.2,.75),
#         legend.background = element_rect(fill="transparent"),
#         panel.grid.major = element_line(color=alpha("gray", 0.2)),
#         legend.text = element_text(size=12)) + 
#   labs(subtitle = "Historical (before alteration)"))
# 
# (p2 <- ggplot() +
#     geom_text_repel(data=df %>% filter(reg=="Y" | gageID %in% c("11413000", "11427000")), 
#                     aes(x=MP, y=maxPowerAvg, label=glue::glue("{river}: {gageID} ({year_end})")),
#                     fontface="bold", family="Roboto Condensed",
#                     segment.alpha = .8,
#                     segment.color="grey10",
#                     color = "black",     # text color
#                     #bg.color = "grey30", # shadow color
#                     #bg.r = 0.15, # only with geom_text
#                     size=2.5,
#                     force = 1.2,
#                     segment.curvature = -0.2,
#                     segment.ncp = 3,
#                     box.padding = 0.7, 
#                     point.padding = 0.7,
#                     nudge_y = 0.35,
#                     min.segment.length = .3)+
#     geom_point(data=df %>% filter(reg=="Y" | gageID %in% c("11413000", "11427000")), 
#                aes(x=MP, y=maxPowerAvg, 
#                    shape=regtype, fill=regtype), color="gray40", size=7, show.legend = T)+
#     
#     geom_point(data=df %>% filter(gageID %in% c("11413000", "11427000")), 
#                aes(x=MP, y=maxPowerAvg, 
#                    shape=regtype), size=7, fill="#56B4E9", show.legend = F, alpha=1) +
#     scale_fill_colorblind("Flow Regime") +
#     scale_shape_manual("Flow Regime", values=c("Bypass"=22, "Hydropeaking"=24, "Unregulated"=21))+
#     xlim(c(0, 1)) + 
#     scale_y_continuous(limits = c(0,18), breaks = c(seq(0,18,3)))+
#     theme_classic(base_family = "Roboto Condensed", base_size = 12) + xlab("Seasonality (M/P)") + ylab("Predictability")+
#     theme(legend.position=c(0.90,.85),
#           legend.background = element_rect(fill="white", color="gray60", size = .5),
#           panel.grid.major = element_line(color=alpha("gray", 0.2)),
#           legend.text = element_text(size=10)) + 
#     labs(subtitle = "Contemporary (with flow alteration)"))
# 
# 
# 
# library(cowplot)
# cowplot::plot_grid(p1, p2, nrow = 2, align = "v", labels = "AUTO")
# ggsave(filename = "figs/fig_02_pre_post_predic_seas_select_sites.png", width = 8, height = 7, units = "in", dpi = 300) #device=cairo_pdf)

#ggsave(filename = "figs/fig_02_pre_post_predic_seas_all_sites.png", width = 8, height = 7, units = "in", dpi = 300) #device=cairo_pdf)

# FIG 2: SINGLE PLOT W ARROWS -----------------------------------------

# add centroid data for paired sites and field for unregulated by value/year
df <- df %>% 
  mutate(regtype2 = case_when(
    grepl("N", reg) ~ "Unregulated",
    TRUE ~ as.character(regtype)
  )) %>% 
  mutate(regtype2 = forcats::fct_relevel(regtype2, "Bypass", "Hydropeaking", "Unregulated")) %>% 
  mutate(
    arrow_from_x = case_when(
      (gageID %in% c("11417500","11433500","11433200") & regtype2=="Unregulated") ~ MP),
    arrow_from_y = case_when(
      (gageID %in% c("11417500","11433500","11433200") & regtype2=="Unregulated") ~ maxPowerAvg),
    arrow_to_x = case_when(
      (gageID %in% c("11417500","11433500","11433200") & regtype2!="Unregulated") ~ MP),
    arrow_to_y = case_when(
      (gageID %in% c("11417500","11433500","11433200") & regtype2!="Unregulated") ~ maxPowerAvg))

# fill values down
df_trim <- df %>% 
  filter(gageID %in% c("11417500","11433500","11433200")) %>%
  group_by(regtype) %>%  
  fill(arrow_from_x, .direction = "down") %>% 
  fill(arrow_from_y, .direction = "down") %>% 
  fill(arrow_to_x, .direction = "up") %>% 
  fill(arrow_to_y, .direction = "up") 
  

df_rev <- df %>% filter(!gageID %in% c("11417500","11433500","11433200")) %>%
  bind_rows(., df_trim)

# save
#write_csv(df_rev, path = "data_output/Table_S2_usgs_gage_sites.csv")

# Make the Plot -----------------------------------------------------------

# THIS IS ORIGINAL PLOT BY REG TYPE

## PLOT
(pAll <- ggplot() +
    
    geom_segment(data=df_rev %>% 
                   filter(gageID %in% c("11417500","11433500","11433200")),
                 aes(x = arrow_from_x, y = arrow_from_y, xend = arrow_to_x, yend = arrow_to_y),
                 color="gray50", lty=2, alpha=0.4, arrow.fill = "gray40", 
                 arrow = ggplot2::arrow(ends = "last", angle = 20, length = unit(0.24, "inches"), type="open")) +
    
    geom_text_repel(data=df_rev, #%>% filter(!gageID %in% c("11417500","11433500","11433200")), 
                    aes(x=MP, y=maxPowerAvg, label=glue::glue("{river}: {gageID} \n({year_end})")),
                    fontface="bold", family="Helvetica",
                    segment.alpha = .8,
                    segment.color="grey10",
                    color = "black",     # text color
                    bg.color = "gray99", # shadow color
                    bg.r = 0.15, # only with geom_text
                    size=2.7,
                    force = 1.4,
                    segment.curvature = -0.2,
                    segment.ncp = 3,
                    box.padding = 0.5, 
                    point.padding = 0.7,
                    nudge_y = 0.9,
                    nudge_x = -0.01,
                    min.segment.length = .3) +
    
    # geom_text_repel(data=df_rev %>% filter(gageID %in% c("11417500","11433500","11433200"), reg=="Y"), 
    #                 aes(x=MP, y=maxPowerAvg, label=glue::glue("{river}: {gageID} \n({year_end})")),
    #                 fontface="bold", family="Helvetica",
    #                 segment.alpha = .8,
    #                 segment.color="grey10",
    #                 color = "black",     # text color
    #                 bg.color = "gray99", # shadow color
    #                 bg.r = 0.15, # only with geom_text
    #                 size=2.7,
    #                 force = 1.2,
    #                 segment.curvature = -0.2,
    #                 segment.ncp = 3,
    #                 box.padding = 0.7, 
    #                 point.padding = 0.5,
    #                 nudge_y = 0.35,
    #                 nudge_x = -0.01,
    #                 min.segment.length = .3) +
    geom_point(data=df_rev,
               aes(x=MP, y=maxPowerAvg, 
                   shape=regtype2, fill=regtype2), color="gray40", size=7, show.legend = T) + 
    scale_fill_colorblind("Flow Regime") +
    scale_shape_manual("Flow Regime", values=c("Bypass"=22, "Hydropeaking"=24, "Unregulated"=21))+
    xlim(c(0, 1)) + 
    scale_y_continuous(limits = c(0,18), breaks = c(seq(0,18,3)))+
    theme_classic(base_family = "Helvetica", base_size = 12) + xlab("Seasonality (M/P)") + ylab("Predictability")+
    theme(legend.position=c(0.1,.8),
          legend.background = element_rect(fill="white", color="gray60", size = .5),
          panel.grid.major = element_line(color=alpha("gray", 0.2)),
          legend.text = element_text(size=10)))


ggsave(filename = "figs/fig_02_pre_post_predic_seas_combined_v2.png", width = 9, height = 6.5, units = "in", dpi = 300)


ggsave(filename = "figs/fig_02_pre_post_predic_seas_combined_v2.pdf", width = 9, height = 6.5, units = "in", dpi = 300, device=cairo_pdf)



# PLOT: ALTERED vs. UNREGULATED -------------------------------------------

# THIS IS MODIFED TO COLLAPSE ALL REG TYPES TO "ALTERED"
# orange: #E69F00
# blue: #56B4E9
# black: 


## PLOT
pAll <- ggplot() +
   
   geom_segment(data=df_rev %>% 
                  filter(gageID %in% c("11417500","11433500","11433200")),
                aes(x = arrow_from_x, y = arrow_from_y, xend = arrow_to_x, yend = arrow_to_y),
                color="gray50", lty=2, alpha=0.4, arrow.fill = "gray40", 
                arrow = ggplot2::arrow(ends = "last", angle = 20, length = unit(0.24, "inches"), type="open")) +
   
   geom_text_repel(data=df_rev, #%>% filter(!gageID %in% c("11417500","11433500","11433200")), 
                   aes(x=MP, y=maxPowerAvg, label=glue::glue("{river}: {gageID} \n({year_end})")),
                   fontface="bold", family="Helvetica",
                   segment.alpha = .8,
                   segment.color="grey10",
                   color = "black",     # text color
                   bg.color = "gray99", # shadow color
                   bg.r = 0.15, # only with geom_text
                   size=2, # was 2.7
                   force = 1.4,
                   segment.curvature = -0.2,
                   segment.ncp = 3,
                   box.padding = 0.5, 
                   point.padding = 0.7,
                   nudge_y = 0.9,
                   nudge_x = -0.01,
                   min.segment.length = .3) +
   
 geom_point(data=df_rev,
            aes(x=MP, y=maxPowerAvg, 
                shape=reg, fill=reg), color="gray40", size=5, show.legend = T) + 
   scale_fill_manual("Flow Regime", values=c("Y"="black", "N"="#56B4E9"), labels=c("Y"="Altered", "N"="Unregulated")) +
   scale_shape_manual("Flow Regime", values=c("Y"=22, "N"=21), labels=c("Y"="Altered", "N"="Unregulated"))+
   xlim(c(0, 1)) + 
   scale_y_continuous(limits = c(0,18), breaks = c(seq(0,18,3)))+
   theme_classic(base_family = "Helvetica", base_size = 11) + xlab("Seasonality (M/P)") + ylab("Predictability")+
   theme(legend.position=c(0.2,.85),
         axis.text = element_text(size=9, color="black"),
         legend.background = element_rect(fill="white", color=NA, size = .5),
         panel.grid.major = element_line(color=alpha("gray", 0.2)),
         legend.text = element_text(size=9, color="black"))


ggsave(plot = pAll, filename = "figs/final_figs/fig_04_pre_post_predic_seas_combined_v3.pdf", 
       width=15.6, height = 11, units="cm", dpi=600, device=cairo_pdf)
       #width = 9, height = 6.5, units = "in", dpi = 300)


#ggsave(filename = "figs/fig_02_pre_post_predic_seas_combined_v3.pdf", width = 9, height = 6.5, units = "in", dpi = 300, device=cairo_pdf)
