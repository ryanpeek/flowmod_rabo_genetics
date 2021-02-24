# fst basic plot

fst_dist_ggplot <- function(data){
  ggplot() + stat_smooth(data=data,
                         aes(y=fst_adj/(1-fst_adj), x=dist_km, 
                             color=regtype, group=regtype), 
                         method = "glm", show.legend = T, lty=2, lwd=1,
                         alpha=0.2, level = 0.89) +
    geom_point(data=data, aes(y=fst_adj/(1-fst_adj), x=dist_km,
                                          fill=regtype, group=regtype, shape=regtype),
               show.legend = T, size=4, alpha=0.9) +
    scale_fill_colorblind("Regulation Type") + 
    #scale_fill_colorblind("Regulation Type")+
    scale_shape_manual("Regulation Type", values=c(21,22,23))+
    scale_color_colorblind("Regulation Type") +
    theme_bw(base_family = "Roboto Condensed") +
    labs(title=expression(paste("F" ["ST"], " vs Distance (km)")),
         y=expression(paste("F" ["ST"], " / (1 - F" ["ST"],")")),
         x="River Distance (km)")
}


# plotly version:
fst_dist_plotly <- function(data){
  
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  plotly::ggplotly(ggplot() + 
                     stat_smooth(data=data,
                                 aes(y=fst_adj/(1-fst_adj), x=dist_km, color=regtype, group=regtype), 
                                 method = "glm", show.legend = T, lty=2, lwd=1,
                                 alpha=0.2, level = 0.89) +
                     geom_point(data=data, aes(y=fst_adj/(1-fst_adj), x=dist_km,
                                                      fill=regtype, group=regtype, shape=regtype, label=site_pair),
                                show.legend = T, size=4, alpha=0.9) +
                     scale_fill_manual("Regulation Type", values = c("unregulated"=cbbPalette[3], "bypass"=cbbPalette[2], "hydropeaking"=cbbPalette[7])) +
                     scale_shape_manual("Regulation Type", values=c(21,22,23))+
                     scale_color_manual("Regulation Type", values = c("unregulated"=cbbPalette[3], "bypass"=cbbPalette[2], "hydropeaking"=cbbPalette[7])) +
                     theme_bw(base_family = "Roboto Condensed")
  )
  
}
