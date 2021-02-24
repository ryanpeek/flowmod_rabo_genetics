# make some FST plots

library(tidyverse)
#library(plotly)
library(ggthemes)
library(ggrepel)

# Get Raw FST & Distance Values for Plot ----------------------------------

load("data_output/fst_dist_pw.rda")

## Add Regulation TYPE using mutate case_when
fst_dist_pw <- fst_dist_pw %>% 
  mutate(
    dist_km = ifelse(is.na(dist_km), 20.0216, dist_km),
    dist_km = round(dist_km, 3),
    fst_adj = round(fst_adj, 4),
    regType = case_when(
      river == "nfmfa" ~ "Bypass",
      regtype == "unregulated" ~ "Unregulated",
      regtype == "bypass" ~ "Bypass",
      regtype == "hydropeaking" ~ "Hydropeaking"),
    REG = case_when(
      grepl("slate", siteA) ~ "REG",
      grepl("slate", siteB) ~ "REG",
      TRUE ~ REG
    )
  )

# fix levels
fst_dist_pw$regType <- factor(fst_dist_pw$regType, levels = c("Unregulated","Bypass", "Hydropeaking"))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# drop sites w low sample numbers
fst_dist_pw <- fst_dist_pw %>% 
  filter(!grepl("scotchman", site_pair),
         !grepl("slate_onion", site_pair),
         !grepl("sfy_mcki", site_pair))

# scale the FST value for plotting:
fst_dist_pw <- fst_dist_pw %>% 
  mutate(fst_adj2 = fst_adj/(1-fst_adj))

# add a unique ID for labeling pairs
fst_dist_pw$sitepairID <- 1:nrow(fst_dist_pw)

# make list of sites:
lsA <- unique(fst_dist_pw$siteA)
lsB <- unique(fst_dist_pw$siteB)

siteList <- unique(c(lsA,lsB)) %>% sort()


# Create Centroids --------------------------------------------------------

# see here: https://stackoverflow.com/questions/23463324/r-add-centroids-to-scatter-plot

reg_centroids <- aggregate(cbind(dist_km, fst_adj2)~regType,fst_dist_pw,mean)
f <- function(z)sd(z)/sqrt(length(z)) # function to calculate std.err
se <- aggregate(cbind(se.x=dist_km,se.y=fst_adj2) ~ regType, fst_dist_pw, f)
reg_centroids <- merge(reg_centroids, se, by="regType") # add SE column

save(reg_centroids, file = "data_output/fst_centroids_all_pairwise_combos.rda")


# Basic Plot --------------------------------------------------------------


# ggplot(fst_dist_pw,aes(y=fst_adj2, x=dist_km, fill=regType)) +
#   geom_point(size=3, pch=21, show.legend = F)+ 
#   # geom_errorbar(data=reg_centroids,color="gray10",
#   #                aes(ymin=fst_adj2-se.y,ymax=fst_adj2+se.y),width=0.03)+
#   #geom_point(data=reg_centroids,size=7, aes(fill=regType), pch=21, color="gray10", alpha=0.9,
#   #           show.legend = F) +
#   #geom_errorbarh(data=reg_centroids,aes(xmin=dist_km-se.x,xmax=dist_km+se.x),height=0.02)+
#   geom_point(data=reg_centroids, size=7, aes(color=regType), pch=23, stroke=2, fill=alpha("gray10",0.5),
#              show.legend = F) +
#   theme_classic(base_size = 14, base_family = "Roboto Condensed") +
#   ylab(expression(paste("Scaled F"["ST"]))) +
#   xlab("River Distance (km)") +
#   scale_color_manual("", values = c("Unregulated"= "steelblue", "Bypass"="orange", "Hydropeaking"="red3")) +
#   scale_fill_manual("", values = c("Unregulated"= "steelblue", "Bypass"="orange", "Hydropeaking"="red3")) #+
#   #scale_color_viridis_d(option="C")+
#   #scale_fill_viridis_d(option = "C") + 
#   #facet_grid(.~regType)

#ggsave(filename = "figs/fig_04a_fst_faceted_w_centroid.png", width = 8, height= 5.5, units = "in", dpi = 300)

# print centroids?
#reg_centroids

# # Plotly by RegType -------------------------------------------------------
# 
# # plotly
# # ggplotly(
# #   ggplot() + 
# #     stat_smooth(data=fst_dist_pw,
# #                 aes(y=fst_adj2, x=dist_km, color=regType, group=regType), 
# #                 method = "glm", show.legend = T, lty=2, lwd=1,
# #                 alpha=0.2, level = 0.89) +
# #     geom_point(data=fst_dist_pw, aes(y=fst_adj2, x=dist_km,
# #                                      fill=regType, group=regType, 
# #                                      shape=regType, label=site_pair),
# #                show.legend = T, size=4, alpha=0.9) +
# #     scale_fill_manual("Regulation Type", values = c("Unregulated"=cbbPalette[3], "Bypass"=cbbPalette[2], "Hydropeaking"=cbbPalette[7])) +
# #     scale_shape_manual("Regulation Type", values=c(21,22,23))+
# #     scale_color_manual("Regulation Type", values = c("Unregulated"=cbbPalette[3], "Bypass"=cbbPalette[2], "Hydropeaking"=cbbPalette[7])) +
# #     theme_bw(base_family = "Roboto Condensed")
# # )
# 
# # Plotly By REG Only ------------------------------------------------------
# 
# # ggplotly(
# #   ggplot() + 
# #     stat_smooth(data=fst_dist_pw,
# #                 aes(y=fst_adj2, x=dist_km, color=REG, group=REG), 
# #                 method = "glm", show.legend = T, lty=2, lwd=1,
# #                 alpha=0.2, level = 0.89) +
# #     geom_point(data=fst_dist_pw, aes(y=fst_adj2, x=dist_km,
# #                                      fill=REG, group=REG, shape=REG, label=site_pair),
# #                show.legend = T, size=4, alpha=0.9) +
# #     scale_fill_manual("Regulated", values = c("REG"=cbbPalette[7], "U"=cbbPalette[3])) +
# #     scale_shape_manual("Regulated", values=c(21,22))+
# #     scale_color_manual("Regulated", values = c("U"=cbbPalette[3], "REG"=cbbPalette[7])) +
# #     theme_bw(base_family = "Roboto Condensed")
# # )

# Print Plot by Reg Type --------------------------------------------------

# print version
(fig4A <- ggplot() + 
    # stat_smooth(data=fst_dist_pw,
    #             aes(y=fst_adj2, x=dist_km, color=regType, group=regType), 
    #             method = "glm", show.legend = T, lty=2, lwd=1,
    #             alpha=0.2, level = 0.89, se = TRUE) +
    geom_point(data=fst_dist_pw, aes(y=fst_adj2, x=dist_km,
                                     fill=regType, group=regType, shape=regType),
               show.legend = T, size=4, alpha=0.9) +
    scale_fill_manual("Flow Alteration", 
                      values = c("Unregulated"=cbbPalette[3], "Bypass"=cbbPalette[2], 
                                 "Hydropeaking"=cbbPalette[7])) +
    #ylim(0,0.35)+
    scale_shape_manual("Flow Alteration", values=c(21,22,23))+
    scale_color_manual("Flow Alteration", 
                       values = c("Unregulated"=cbbPalette[3], "Bypass"=cbbPalette[2], 
                                  "Hydropeaking"=cbbPalette[7])) +
    geom_point(data=reg_centroids, aes(y=fst_adj2, x=dist_km,
                                       group=regType, shape=regType, fill=regType),
               stroke=1.5, color="gray10",
               show.legend = F, size=7, alpha=0.9)+
    
    # # add ellipses
    # geom_mark_ellipse(data=fst_dist_pw, aes(y=fst_adj2, x=dist_km, group=regType, color=regType),
    #                   fill=NA, show.legend = F, alpha=0.2) +
    
    # add labels
    # geom_text_repel(data=fst_dist_pw, 
    #                 aes(y=fst_adj2, x=dist_km, label=sitepairID, group=regType), 
    #                 alpha=0.8, size=2.8, point.padding = 0.1) +
    theme_bw(base_family = "Helvetica", base_size = 12) +
    theme(legend.position = c(0.85, 0.18)) +
    labs(#title=expression(paste("Mean F" ["ST"], " vs Mean Distance (km)")),
      #y=expression(paste("Mean F" ["ST"], " / (1 - Mean F" ["ST"],")")),
      y=expression(paste("Scaled F" ["ST"])),
      x="River Distance (km)") #+
  #facet_grid(.~regType)
)

# save(fig4A, reg_centroids, fst_dist_pw, cbbPalette, file = "figs/fig_04a_fst_plot_unfaceted.rda")
# 
# ggsave(filename = "figs/fig_04a_rev_fst_vs_dist_by_flowtype_centroid.png", width = 8, height= 6, units = "in", dpi = 300)
# 
# ggsave(filename = "figs/fig_04a_rev_fst_vs_dist_by_flowtype_centroid_no_regression.png", width = 8, height= 6, units = "in", dpi = 300)
# 
# ggsave(filename = "figs/fig_04a_rev_fst_vs_dist_by_flowtype_faceted_wide_centroid.png", width = 7, height= 4, units = "in", dpi = 300)
# ggsave(filename = "figs/fig_03_rev_fst_vs_dist_by_flowtype.tiff", width = 5, height= 4, units = "in", dpi = 600)
# ggsave(filename = "figs/fig_03a_fst_vs_dist_by_flowtype_mfa_nfmfa_resize.pdf", width = 3.42, height = 2.9, units = "in", dpi = 600, scale=1.3)

# save out for table
fst_dist_pw %>% select(siteA, siteB, REG, river, regtype, dist_km, regType, fst_adj, fst_adj2, sitepairID) %T>% 
  View() %>% 
  write_csv(., path = "manuscript/supplemental/S4_fst_dist_all_pairwise_combos_revised.csv")

### OLD CODE -------
# # * Load FST Data ------------------------------------------------------------
# 
# # the FST data
# load("data_output/fst_within_all_v3.rda")
# 
# # plot Fsts
# 
# # get names for matrix
# nameVals <- sort(unique(unlist(fst[c(4:5)]))) # cols w pair names
# 
# # construct 0 matrix of correct dimensions with row and column names
# myMat <- matrix(NA, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
# 
# # fill in the matrix with matrix indexing on row and column names
# #myMat[as.matrix(fst[c("siteA", "siteB")])] <- fst[["fst_adj"]]
# 
# #plotly::ggplotly(
#   ggplot(data=fst, aes(x=siteA, y=siteB, fill=fst_adj)) + 
#     geom_tile() +
#     scale_fill_viridis_c()
# #)
# 
# # * BEAR dist_matrix -------------------------
# 
# # read in dist matrix: BEAR
# dist_bear <-read_tsv("data_output/site_dist_matrix_BEAR.txt") %>%
#   as.matrix()
# colnames(dist_bear) <- gsub("-", replacement = "_", colnames(dist_bear))
# colnames(dist_bear)<-tolower(colnames(dist_bear))
# rownames(dist_bear)<-colnames(dist_bear)
# 
# # melt: 
# dist_bear_df <- reshape2:::melt.matrix(dist_bear, varnames = c("siteA", "siteB"), value.name = "dist_km") %>% 
#   mutate(dist_km=dist_km/1000,
#          site_pair=paste0(siteA, ".", siteB)) %>% 
#   filter(dist_km!=0) %>% 
#   mutate(dist_km = round(dist_km, 3)) %>%
#   distinct(dist_km, .keep_all = TRUE)
# 
# rm(dist_bear)
# 
# # * MFA dist_matrix -------------------------
# 
# # read in dist matrix: AMER
# dist_amer <-read_tsv("data_output/site_dist_matrix_AMER.txt") %>% 
#   as.matrix()
# colnames(dist_amer) <- gsub("-", replacement = "_", colnames(dist_amer))
# colnames(dist_amer)<-tolower(colnames(dist_amer))
# rownames(dist_amer)<-colnames(dist_amer)
# 
# # get MFA drainage only
# (dist_mfa <- dist_amer[grepl("mfa|nfmfa|rub",rownames(dist_amer)),grepl("mfa|nfmfa|rub",colnames(dist_amer))])
# 
# # melt: 
# dist_mfa_df <- reshape2:::melt.matrix(dist_mfa, varnames = c("siteA", "siteB"), value.name = "dist_km") %>% 
#   mutate(dist_km=dist_km/1000,
#          site_pair=paste0(siteA, ".", siteB)) %>% 
#   filter(dist_km!=0) %>% 
#   mutate(dist_km = round(dist_km, 3)) %>%
#   distinct(dist_km, .keep_all = TRUE)
# 
# rm(dist_mfa)
# rm(dist_amer)
# 
# # * NFA dist_matrix -------------------------
# 
# # read in dist matrix: NFA
# #dist_nfa <-read_tsv("data_output/site_dist_matrix_NFA.txt") %>% 
# dist_nfa <-read_tsv("data_output/site_dist_matrix_AMER.txt") %>% 
#   as.matrix()
# colnames(dist_nfa) <- gsub("-", replacement = "_", colnames(dist_nfa))
# colnames(dist_nfa)<-tolower(colnames(dist_nfa))
# rownames(dist_nfa)<-colnames(dist_nfa)
# 
# # calc rowmeans for single dist 
# (dist_nfa <- dist_nfa[grepl("nfa",rownames(dist_nfa)),grepl("nfa",colnames(dist_nfa))])
# 
# # melt: 
# dist_nfa_df <- reshape2:::melt.matrix(dist_nfa, varnames = c("siteA", "siteB"), value.name = "dist_km") %>% 
#   mutate(dist_km=dist_km/1000,
#          site_pair=paste0(siteA, ".", siteB)) %>% 
#   filter(dist_km!=0) %>% 
#   mutate(dist_km = round(dist_km, 3)) %>%
#   distinct(dist_km, .keep_all = TRUE)
# 
# rm(dist_nfa)
# 
# # * YUBA dist_matrix -------------------------
# 
# # read in dist matrix: YUBA
# dist_yuba <-read_tsv("data_output/site_dist_matrix_YUBA.txt") %>% 
#   as.matrix()
# colnames(dist_yuba) <- gsub("-", replacement = "_", colnames(dist_yuba))
# colnames(dist_yuba)<-tolower(colnames(dist_yuba))
# rownames(dist_yuba)<-colnames(dist_yuba)
# dist_yuba <- dist_yuba[!grepl("deer_clearck", rownames(dist_yuba)),
#                        !grepl("deer_clearck", colnames(dist_yuba))]
# 
# # calc rowmeans for single dist 
# #dist_nfy <- dist_yuba[grepl("mfy|nfy",rownames(dist_yuba)),grepl("mfy|nfy",colnames(dist_yuba))]
# dist_nfy <- dist_yuba[grepl("nfy",rownames(dist_yuba)),grepl("nfy",colnames(dist_yuba))]
# 
# # filter out sites: 
# # NFY
# (dist_nfy <- dist_nfy[!grepl("nfy_slate$|nfy_slate_onion$", rownames(dist_nfy)),!grepl("nfy_slate$|nfy_slate_onion$",colnames(dist_nfy))])
# # melt: 
# dist_nfy_df <- reshape2:::melt.matrix(dist_nfy, varnames = c("siteA", "siteB"), value.name = "dist_km") %>% 
#   mutate(dist_km=dist_km/1000,
#          site_pair=paste0(siteA, ".", siteB)) %>% 
#   filter(dist_km!=0) %>% 
#   mutate(dist_km = round(dist_km, 3)) %>%
#   distinct(dist_km, .keep_all = TRUE)
# 
# rm(dist_nfy)
# 
# 
# # MFY
# (dist_mfy <- dist_yuba[grepl("mfy",rownames(dist_yuba)),grepl("mfy",colnames(dist_yuba))])
# # melt: 
# dist_mfy_df <- reshape2:::melt.matrix(dist_mfy, varnames = c("siteA", "siteB"), value.name = "dist_km") %>% 
#   mutate(dist_km=dist_km/1000,
#          site_pair=paste0(siteA, ".", siteB)) %>% 
#   filter(dist_km!=0) %>% 
#   mutate(dist_km = round(dist_km, 3)) %>%
#   distinct(dist_km, .keep_all = TRUE)
# 
# rm(dist_mfy)
# 
# # SFY
# (dist_sfy <- dist_yuba[grepl("sfy",rownames(dist_yuba)),grepl("sfy",colnames(dist_yuba))])
# 
# # melt: 
# dist_sfy_df <- reshape2:::melt.matrix(dist_sfy, varnames = c("siteA", "siteB"), value.name = "dist_km") %>% 
#   mutate(dist_km=dist_km/1000,
#          site_pair=paste0(siteA, ".", siteB)) %>% 
#   filter(dist_km!=0) %>% 
#   mutate(dist_km = round(dist_km, 3)) %>% 
#   distinct(dist_km, .keep_all = TRUE)
# 
# rm(dist_sfy)
# rm(dist_yuba)
# 
# 
# # * Now BIND DISTANCES TOGETHER ---------------------------------------------
# 
# # bind into one DF
# 
# dists_all <- bind_rows(dist_bear_df, dist_mfa_df, dist_mfy_df, dist_nfa_df, dist_nfy_df, dist_sfy_df)
# 
# #rm(list = ls(pattern = "dist\\_"))
# 
# # fix names NFA and RUB
# dists_all$siteA <- gsub("nfa_euchds$", "nfa_euch", dists_all$siteA)
# dists_all$siteB <- gsub("nfa_euchds$", "nfa_euch", dists_all$siteB)
# dists_all$siteA <- gsub("rub_lc_us$", "rub_lcus", dists_all$siteA)
# dists_all$siteB <- gsub("rub_lc_us$", "rub_lcus", dists_all$siteB)
# 
# # fix names BEAR
# fst$siteA <- gsub("^bear_sthohaw$", "bear_stho_haw", fst$siteA) # sthc is lower
# fst$siteB <- gsub("^bear_sthohaw$", "bear_stho_haw", fst$siteB) #sth2 is upper
# dists_all$siteA <- gsub("bear_stha", "bear_stho_haw", dists_all$siteA)
# dists_all$siteB <- gsub("bear_stha", "bear_stho_haw", dists_all$siteB)
# 
# #fst$siteA <- gsub("bear_sthc$", "bear_stho", fst$siteA) # sthc is lower
# #fst$siteB <- gsub("bear_sthc$", "bear_stho", fst$siteB) #sth2 is upper
# 
# fst$siteA <- gsub("^bear_stha$", "bear_stho_haw", fst$siteA) # sthc is lower
# fst$siteB <- gsub("^bear_stha$", "bear_stho_haw", fst$siteB) #sth2 is upper
# 
# # filter out:
# dists_all <- dists_all %>% 
#   filter(!siteA=="nfa_euchus", !siteB=="nfa_euchus", !siteA=="deer_clearck")
# 
# # filter out:
# fst <- fst %>% 
#   filter(!siteA=="nfa_euchus", !siteB=="nfa_euchus", !siteA=="deer_clearck") %>% 
#   filter(riverA==riverB)
# 
# 
# 
# 
# 
# 
# 
