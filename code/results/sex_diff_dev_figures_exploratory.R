library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(gratia) 
library(knitr)
library(mgcv)
library(RColorBrewer)
library(scales)
library(stringr)
library(rjson)
library(tidyr)

######################################
# Figures 
######################################

# format mean matrix
format_mean_matrix <- function(mat) {
  # rename schaefer 200x17 to be schaefer 200x7
  rownames(mat) <- schaefer200.labels
  colnames(mat) <- schaefer200.labels 
  return(mat)
}

# map parcels to network names for heatmap
map_network7 <- function(x) {
  dplyr::case_when(
    str_detect(x, "Vis") ~ "Visual",
    str_detect(x, "SomMot") ~ "SomMot",
    str_detect(x, "DorsAttn") ~ "DorsAttn",
    str_detect(x, "SalVentAttn") ~ "SalVentAttn",
    str_detect(x, "Cont") ~ "Frontoparietal",
    str_detect(x, "Limb") ~ "Limbic",
    str_detect(x, "Default") ~ "Default",
    TRUE ~ x
  )
}

map_network17 <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "VisCent")           ~ "visualCentral",
    stringr::str_detect(x, "VisPeri")           ~ "visualPeripheral",
    stringr::str_detect(x, "SomMotA")           ~ "somatomotorA",
    stringr::str_detect(x, "SomMotB")           ~ "somatomotorBAuditory",
    stringr::str_detect(x, "DorsAttnA")         ~ "dorsalAttentionA",
    stringr::str_detect(x, "DorsAttnB")         ~ "dorsalAttentionB",
    stringr::str_detect(x, "SalVentAttnA")      ~ "salienceVentralAttentionA",
    stringr::str_detect(x, "SalVentAttnB")      ~ "salienceVentralAttentionB",
    stringr::str_detect(x, "LimbicA_TempPole")  ~ "limbicTemporopolar",
    stringr::str_detect(x, "LimbicB_OFC")       ~ "limbicOrbitofrontal",
    stringr::str_detect(x, "ContA")             ~ "frontoparietalControlA",
    stringr::str_detect(x, "ContB")             ~ "frontoparietalControlB",
    stringr::str_detect(x, "ContC")             ~ "frontoparietalControlC",
    stringr::str_detect(x, "DefaultA")          ~ "defaultA",
    stringr::str_detect(x, "DefaultB")          ~ "defaultB",
    stringr::str_detect(x, "DefaultC")          ~ "defaultC",
    stringr::str_detect(x, "TempPar")           ~ "temporoparietal",
    TRUE                                  ~ x     # pass through original label if no match
  )
}
# helper function for plot_mean_matrix - compute network block boundaries/centers for lines + labels
runs_info <- function(groups) {
  r <- rle(groups)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  centers <- starts + (r$lengths - 1) / 2
  list(starts = starts, ends = ends, centers = centers, labels = r$values)
}

# plot mean matrix heatmap
plot_mean_matrix <- function(mat, ylim1, ylim2, sex) {
  
  # save parcel names in a variable
  orig_rows <- rownames(mat)
  orig_cols <- colnames(mat)
  
  # map parcel labels to network names 
  row_net <- map_network7(orig_rows)
  col_net <- map_network7(orig_cols)    
  
  # set display order of networks  
  net_order <- c("Visual","SomMot","DorsAttn","SalVentAttn", "Frontoparietal","Default", "Limbic")
  
  # build ordering indices: first by network, then keep original within-network order
  row_order <- order(factor(row_net, levels = net_order), seq_along(row_net))
  col_order <- order(factor(col_net, levels = net_order), seq_along(col_net))
  
  # reorder matrix and corresponding network vectors
  mat_ord <- mat[row_order, col_order]
  row_net_o <- row_net[row_order]
  col_net_o <- col_net[col_order]
  
  # compute network block boundaries/centers for lines + labels
  r_rows <- runs_info(row_net_o)
  r_cols <- runs_info(col_net_o)
  y_bounds <- r_rows$ends + 0.5
  x_bounds <- r_cols$ends + 0.5
  
  # plot full 200x200 heatmap with network blocks
  df <- expand.grid(Col = seq_len(ncol(mat_ord)), Row = seq_len(nrow(mat_ord)))
  df$Value <- as.vector(mat_ord)
  
  heatmap <- ggplot(df, aes(x = Col, y = Row, fill = Value)) +
    geom_tile() +
    geom_hline(yintercept = y_bounds, linewidth = 0.4, color = "black") +
    geom_vline(xintercept = x_bounds, linewidth = 0.4, color = "black") +
    scale_y_continuous(breaks = r_rows$centers, labels = r_rows$labels, expand = c(0,0)) +
    scale_x_continuous(breaks = r_cols$centers, labels = r_cols$labels, expand = c(0,0), position = "bottom") +
    scale_fill_gradientn(colours = c("#6200D9", "white", "#FF681C"), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    coord_equal() +
    theme_classic() + 
    theme(plot.title = element_text(size = 20, hjust = 0.5),
          axis.text.x = element_text(size = 12, angle = 45, vjust = 1.1, hjust = 1, margin = margin(t = 10)),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = 'none') + ggtitle(sex)  
  
  return(heatmap)
}


# prep for network x network difference heatmap: collapse a parcel-level matrix to a network x network matrix
collapse_to_network <- function(mat, map_network7, net_order = NULL) {
  row_net <- map_network7(rownames(mat))
  col_net <- map_network7(colnames(mat))
  
  # indices per network
  row_idx <- split(seq_len(nrow(mat)), row_net)
  col_idx <- split(seq_len(ncol(mat)), col_net)
  
  # pick an order
  if (is.null(net_order)) {
    net_order <- union(names(row_idx), names(col_idx)) |> unique()
  }
  row_idx <- row_idx[net_order]
  col_idx <- col_idx[net_order]
  
  # mean across all parcel pairs for each network pair
  net_mat <- matrix(NA_real_, nrow = length(row_idx), ncol = length(col_idx),
                    dimnames = list(names(row_idx), names(col_idx)))
  for (i in seq_along(row_idx)) {
    for (j in seq_along(col_idx)) {
      net_mat[i, j] <- mean(mat[row_idx[[i]], col_idx[[j]]], na.rm = TRUE)
    }
  }
  net_mat
}

# plot network x network difference heatmap (Female – Male)
plot_network_diff <- function(female_mat, male_mat, map_network7,
                              net_order = c("Default","DorsAttn","Frontoparietal",
                                            "Limbic","SalVentAttn","SomMot","Visual"),
                              triangle = c("upper","lower","full"),
                              limits, title) {
  
  triangle <- match.arg(triangle)
  
  female_net <- collapse_to_network(female_mat, map_network7, net_order)
  male_net   <- collapse_to_network(male_mat, map_network7, net_order)
  diff_net   <- female_net - male_net
  
  # Melt to long format for ggplot
  df <- melt(diff_net, varnames = c("Row","Col"), value.name = "diff")
  
  # keep one triangle if requested
  df$Ri <- as.integer(factor(df$Row, levels = rownames(diff_net)))
  df$Ci <- as.integer(factor(df$Col, levels = colnames(diff_net)))
  if (triangle == "upper") df <- df[df$Ri <= df$Ci,]
  if (triangle == "lower") df <- df[df$Ri >= df$Ci,]
  
  ggplot(df, aes(x = Col, y = Row, fill = diff)) +
    geom_tile() +
    coord_fixed() +
    scale_fill_gradientn(
      colours = c("#02B0DB","white","#FF8D6E"),
      values  = rescale(c(limits[1], 0, limits[2])),
      limits  = limits,
      oob     = squish,
      na.value = "white") +
    scale_x_discrete(position = "bottom") +
    ggtitle(title) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(size = 20, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none")
}


# function for calculating number of significant regions (FDR corrected)
# @param df A dataframe with GAM results
# @param fc_measure, A string for FC measure (i.e. "GBC")
# @param covariate, A string for covariate of interest (i.e. "sex_diff")
sig.effects <- function(df, fc_measure, covariate, dataset){
  df$Anova.cov.pvalue.fdr <- p.adjust(df$GAM.cov.pvalue, method = c("fdr")) 
  df$significant.fdr <- df$Anova.cov.pvalue.fdr < 0.05
  df$significant.fdr[df$significant.fdr == TRUE] <- 1
  df$significant.fdr[df$significant.fdr == FALSE] <- 0
  df <- df %>% mutate(significant_status = ifelse(significant.fdr == 1, "Sig", "Nonsig"))
  df <- df %>% mutate(cov_sig_effect = ifelse(significant.fdr == 1, GAM.cov.tvalue, NA))
  sigeffect.totaln <- sum(df$significant.fdr)
  sigeffect.percent <- round(sigeffect.totaln/length(df$significant.fdr)*100, 2)
  print(sprintf("%s: %s regions (%s percent) show a significant association between %s and %s", dataset, sigeffect.totaln, sigeffect.percent, fc_measure, covariate))
  return(df)
}
 
# function for calculating number of significant regions (FDR corrected) - for analysis with aggregated datasets
sig.effects.all_datasets <- function(df, fc_measure, covariate, dataset, pvalue_colname, effect_colname){
  df$Anova.int.pvalue.fdr <- p.adjust(df[[pvalue_colname]], method = c("fdr")) 
  df$significant.fdr <- df$Anova.int.pvalue.fdr < 0.05
  df$significant.fdr[df$significant.fdr == TRUE] <- 1
  df$significant.fdr[df$significant.fdr == FALSE] <- 0
  df <- df %>% mutate(significant_status = ifelse(significant.fdr == 1, "Sig", "Nonsig"))
  df <- df %>% mutate(int_sig_effect = ifelse(significant.fdr == 1, .[[effect_colname]], NA))
  sigeffect.totaln <- sum(df$significant.fdr)
  sigeffect.percent <- round(sigeffect.totaln/length(df$significant.fdr)*100, 2)
  print(sprintf("%s: %s regions (%s percent) show a significant association between %s and %s", dataset, sigeffect.totaln, sigeffect.percent, fc_measure, covariate))
  return(df)
}


# function for calculating number of significant interactions (FDR corrected)
# @param df A dataframe with GAM results
# @param fc_measure, A string for FC measure (i.e. "GBC")
# @param covariate, A string for covariate of interest (i.e. "sex_diff")
sig.effects.interaction <- function(df, measure, covariate, dataset){
  df$Anova.int.pvalue.fdr <- p.adjust(df$GAM.int.pvalue, method = c("fdr")) 
  df$significant.fdr <- df$Anova.int.pvalue.fdr < 0.05
  df$significant.fdr[df$significant.fdr == TRUE] <- 1
  df$significant.fdr[df$significant.fdr == FALSE] <- 0
  df <- df %>% mutate(significant_status = ifelse(significant.fdr == 1, "Sig", "Nonsig"))
  df <- df %>% mutate(int_sig_effect = ifelse(significant.fdr == 1, GAM.int.Fvalue, NA))
  sigeffect.totaln <- sum(df$significant.fdr)
  sigeffect.percent <- round(sigeffect.totaln/length(df$significant.fdr)*100, 2)
  print(sprintf("%s: %s regions (%s percent) show a significant interaction between %s and %s", dataset, sigeffect.totaln, sigeffect.percent, measure, covariate))
  return(df)
}


# plot covariate effects on cortical surface
plot_cortex <- function(df, hemi, measure, ylim1, ylim2, atlas) {
  
  if(hemi == "left") {
    cortical_pos1 <- "left lateral" 
    cortical_pos2 <- "left medial"
  } else{
    cortical_pos1 <- "right lateral" 
    cortical_pos2 <- "right medial"
  }
  plot_lateral <- ggplot() + 
    geom_brain(data = df, atlas = get(atlas), 
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.7)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos1)) +
    scale_fill_gradientn(colours = c("#FF8D6E", "white", "#02B0DB"), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey50")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"),
          plot.title = element_blank()) 
  
  plot_medial <- ggplot() + 
    geom_brain(data = df, atlas= get(atlas), 
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.7)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos2)) +
    scale_fill_gradientn(colours = c("#FF8D6E", "white", "#02B0DB"), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey50")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"),
          plot.title = element_blank())
  return(plot_grid(plot_lateral, plot_medial, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr"))
}


# plot values on cortical surface (for interaction plot - temporary)
plot_cortex_int <- function(df, hemi, measure, ylim1, ylim2, atlas) {
  
  if(hemi == "left") {
    cortical_pos1 <- "left lateral" 
    cortical_pos2 <- "left medial"
  } else{
    cortical_pos1 <- "right lateral" 
    cortical_pos2 <- "right medial"
  }
  plot_lateral <- ggplot() + 
    geom_brain(data = df, atlas = get(atlas), 
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.7)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos1)) +
    scale_fill_paletteer_c("ggthemes::Blue-Green Sequential", na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey50")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"),
          plot.title = element_blank()) 
  
  plot_medial <- ggplot() + 
    geom_brain(data = df, atlas= get(atlas), 
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.7)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos2)) +
    scale_fill_paletteer_c("ggthemes::Blue-Green Sequential", na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey50")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"),
          plot.title = element_blank())
  return(plot_grid(plot_lateral, plot_medial, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr"))
}



# make df for plotting: need to switch network 1 and 2 ordering for salience-dorsalattn and limbic-dorsal for visualization purposes
# also does fdr correction on non-duplicated pairs only
prep_networkpair_df <- function(df, fc_measure, covariate, dataset) {
  df$sorted_pairs <- apply(df[c("network1", "network2")], 1, 
                           function(x) toString(sort(x)))  
  df <- df %>% group_by(sorted_pairs) %>% distinct(sorted_pairs, .keep_all = T)  
  df <- sig.effects(df, fc_measure, covariate, dataset)

  # switch network1 and 2 labeling for salience + dorsal 
  switch1 <- which(df$network1=="dorsalAttention" & df$network2=="salienceVentralAttention")
  orig_network1 <- df$network1[switch1]
  orig_network2 <- df$network2[switch1]
  df$network1[switch1] <- orig_network2
  df$network2[switch1] <- orig_network1
  
  # switch network1 and 2 labeling for limbic + dorsal 
  switch2 <- which(df$network1=="dorsalAttention" & df$network2=="limbic")
  orig_network1 <- df$network1[switch2]
  orig_network2 <- df$network2[switch2]
  df$network1[switch2] <- orig_network2
  df$network2[switch2] <- orig_network1
  
  # switch network1 and 2 labeling for frontoparietal + dorsal 
  switch3 <- which(df$network1=="dorsalAttention" & df$network2=="frontoparietalControl")
  orig_network1 <- df$network1[switch3]
  orig_network2 <- df$network2[switch3]
  df$network1[switch3] <- orig_network2
  df$network2[switch3] <- orig_network1
  return(df)
}


# make plots 
make_networkpair7_plot <- function(df, dataset, ylim1, ylim2) {
  # pretty labels
  lab_map <- c(
    Visual="Visual", SomMot="Somatomotor", DorsAttn="Dorsal Attn",
    SalVentAttn="Salience Ventral Attn", Limbic="Limbic",
    Frontoparietal="Frontoparietal Control", Default="Default"
  )
  df <- df |>
    dplyr::mutate(
      network1_label = lab_map[network1],
      network2_label = lab_map[network2]
    )
  
  # desired order
  net_levels <- c("Default","Dorsal Attn","Frontoparietal Control",
                  "Limbic","Salience Ventral Attn","Somatomotor","Visual")
  
  dfp <- df |>
    dplyr::mutate(
      # IMPORTANT: same (non-reversed) levels on both axes
      x = factor(network2_label, levels = net_levels),  # <- swap
      y = factor(network1_label, levels = net_levels)   # <- swap
    )
  
  ggplot(dfp, aes(x = x, y = y, fill = GAM.cov.tvalue)) +
    geom_tile() +
    geom_text(aes(label = dplyr::case_when(
      Anova.cov.pvalue.fdr < 1e-4 ~ "***",
      Anova.cov.pvalue.fdr < 1e-3 ~ "**",
      Anova.cov.pvalue.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )), colour = "black", size = 10) +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +  # not reversed
    coord_fixed() +
    theme_classic(base_size = 24) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 24),
      axis.text.y = element_text(size = 24),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.key.width = unit(2, 'cm'), legend.key.height = unit(1, 'cm'),
      legend.title = element_blank(), legend.text = element_text(size = 24),
      legend.direction = "horizontal", legend.position = "bottom"
    ) +
    scale_fill_gradientn(
      colours = c("#FF8D6E", "white", "#02B0DB"),
      values = scales::rescale(c(ylim1, 0, ylim2)),
      limits = c(ylim1, ylim2),
      oob = scales::squish, na.value = "white"
    ) +
    ggtitle(dataset)
}


make_networkpair17_plot <- function(df, dataset, ylim1, ylim2) {
  # pretty labels
  lab_map <- c(
    visualCentral="Visual Central", visualPeripheral="Visual Peripheral", 
    somatomotorA="Somatomotor A", somatomotorBAuditory="Somatomotor B Auditory", 
    dorsalAttentionA="Dorsal Attn A",dorsalAttentionB="Dorsal Attn B",
    salienceVentralAttentionA="Salience Ventral Attn A", salienceVentralAttentionB="Salience Ventral Attn B", 
    limbicOrbitofrontal="Limbic OFC", limbicTemporopolar="Limbic TempPolar",
    temporoparietal="Temporoparietal",
    frontoparietalControlA="Frontoparietal A",  frontoparietalControlB="Frontoparietal B",  frontoparietalControlC="Frontoparietal C", 
    defaultA="Default A", defaultB="Default B", defaultC="Default C"
  )
  df <- df |>
    dplyr::mutate(
      network1_label = lab_map[network1],
      network2_label = lab_map[network2]
    )
  
  # desired order
  net_levels <- c("Default A","Default B","Default C", "Dorsal Attn A", "Dorsal Attn B",
                  "Frontoparietal A", "Frontoparietal B","Frontoparietal C",
                  "Limbic OFC", "Limbic TempPolar", "Salience Ventral Attn A", "Salience Ventral Attn B", 
                  "Temporoparietal", "Somatomotor A", "Somatomotor B Auditory", "Visual Central", "Visual Peripheral")
  
  dfp <- df |>
    dplyr::mutate(
      # IMPORTANT: same (non-reversed) levels on both axes
      x = factor(network2_label, levels = net_levels),  # <- swap
      y = factor(network1_label, levels = net_levels)   # <- swap
    )
  
  ggplot(dfp, aes(x = x, y = y, fill = GAM.cov.tvalue)) +
    geom_tile() +
    geom_text(aes(label = dplyr::case_when(
      Anova.cov.pvalue.fdr < 1e-4 ~ "***",
      Anova.cov.pvalue.fdr < 1e-3 ~ "**",
      Anova.cov.pvalue.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )), colour = "black", size = 6) +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +  # not reversed
    coord_fixed() +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.key.width = unit(2, 'cm'), legend.key.height = unit(1, 'cm'),
      legend.title = element_blank(), legend.text = element_text(size = 12),
      legend.direction = "horizontal", legend.position = "bottom"
    ) +
    scale_fill_gradientn(
      colours = c("#FF8D6E", "white", "#02B0DB"),
      values = scales::rescale(c(ylim1, 0, ylim2)),
      limits = c(ylim1, ylim2),
      oob = scales::squish, na.value = "white"
    ) +
    ggtitle(dataset)
}

 


# edge level analysis: 
# make a double df for symmetry - parcel1 and parcel2.SA.vec repeated in opposite ordering (see lines 500-510 in Adam's Edge-level-Age.md)
make_double_df <- function(df, fc_measure, covariate, dataset) {
  age_SA.diff <- cbind(df, SA.diff)
  age_SA.diff <- sig.effects(age_SA.diff, fc_measure, covariate, dataset) 
  
  # make a double df for symmetry - parcel1 and parcel2.SA.vec repeated in opposite ordering (see lines 500-510 in Adam's Edge-level-Age.md)
  age_SA.diff2 <- age_SA.diff 
  age_SA.diff2$parcel1.SA.vec <- age_SA.diff$parcel2.SA.vec
  age_SA.diff2$parcel2.SA.vec <- age_SA.diff$parcel1.SA.vec
  
  # "stacked" df
  double_age_SA.diff <- rbind(age_SA.diff2,age_SA.diff)

  return(double_age_SA.diff)
}

 

# model the surface
make_surface_plot <- function(double_SA.diff, dataset, measure, ylim1, ylim2) {
  g2 <- gam(get(measure) ~ te(parcel2.SA.vec,parcel1.SA.vec,k=3), data = double_SA.diff)
  
  surface.plot.fig <- gg_tensor(g2) + theme_classic() +
    theme(
      axis.title.x=element_blank(),
      axis.title.y=element_blank(),
      axis.line = element_line(color = "black"),
      axis.text=element_text(size=24, color = "black")) + 
    scale_fill_gradientn(colours = c("#FF8D6E", "white", "#02B0DB"), 
                         values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    geom_vline(xintercept = mean(range(double_SA.diff$parcel1.SA.vec)),
               linetype='dashed',size=1) + 
    geom_hline(yintercept = mean(range(double_SA.diff$parcel1.SA.vec)),
               linetype='dashed',size=1) + labs(title=dataset) + 
    theme(legend.position='none',
          legend.key.height = unit(1, 'cm'),
          legend.key.width = unit(3, 'cm'),
          legend.title = element_blank(),
          legend.text = element_text(size=20),
          strip.background = element_blank(),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          plot.title = element_text(size=20, hjust=0.5))  
  return(surface.plot.fig)
}


# extract and format correlation estimate
get_corr <- function(x, y) {
  round(cor.test(x, y, complete.obs = T)$estimate, 3)
}


# plot S-A scatterplot
plot_corr_SA <- function(df, measure, annot_text, title, ylim1, ylim2, y_label) {
  SA_plot <- ggplot(df, aes(x = SA.axis_rank, y = get(measure), fill = SA.axis_rank_sig)) +
    geom_point(color = "gray", shape = 21, size=3.5) + 
    paletteer::scale_fill_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(df$SA.axis_rank), max(df$SA.axis_rank)), oob = squish) +
    paletteer::scale_color_paletteer_c("grDevices::RdYlBu", direction = -1, limits = c(min(df$SA.axis_rank), max(df$SA.axis_rank)), oob = squish) +
    
    geom_smooth(data = df, method='lm', se=TRUE, fill=alpha(c("gray70"),.9), col="black") + 
    labs(title = title) + theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20, color = "black"),
          axis.title = element_blank(),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black")) + 
    annotate(geom="text", x=100, y=y_label, label=annot_text, color="black", size=7) +
    ylim(ylim1, ylim2)
  
  return(SA_plot)
}
 

# developmental trajectories by sex_diff: for male and female  
plot_trajectories_bycov <- function(region_name, df, covariate, y.limits, color){
  
  df_plot <- df %>% filter(region %in% region_name)
  plot <- ggplot(df_plot, aes(x = age, y = fitted, group = interaction(sex, region))) + 
    geom_line(color = color, linewidth = 1, aes(linetype = sex)) +
    theme_classic() +
    scale_y_continuous(limits = y.limits) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.title = element_blank()) + ggtitle(gsub("[A-Z0-9]+Networks ", "", gsub("_", " ", region_name)))
  return(plot)  
}

# format edge GAM output for sex effect seed based maps 
prep_edge_dfs <- function(df, fc_measure, covariate, dataset) {
  df <- sig.effects(df, fc_measure, covariate, dataset)
  df$significant.fdr <- gsub("0", "Non_Sig", df$significant.fdr)
  df$significant.fdr <- gsub("1", "Sig", df$significant.fdr)
  df <- df %>% mutate(network1=sub("(_to_).*", "", region), network2=sub(".*(_to_)", "", region)) 
  df$network1 <- gsub("_to", "", df$network1)
  df$network2 <- gsub("_to_", "", df$network2)
  
  return(df)
}

# sex effect maps for a seed region
view_age_effect <- function(df, parcel, measure, legend_title, dataset){
  gam.edge <- df
  
  parcel_df <- gam.edge %>% 
    filter(network1 == parcel | network2 == parcel) %>% 
    arrange(network2)
  parcel_df$region <- ifelse(parcel_df$network1 == parcel, parcel_df$network2, parcel_df$network1) # get regions that parcel is connected to 
  
  parcel_df[nrow(parcel_df) + 1,] <- c(rep(NA, ncol(parcel_df)-4), "Sig", c(0,0), parcel)
  parcel_df$region <- ifelse(str_detect(parcel_df$region, "LH"), paste0("lh_7", parcel_df$region), paste0("rh_7", parcel_df$region)) # reformat 'region' to match ggseg schaefer 'label'
  
  
  # add lh_Background+FreeSurfer_Defined_Medial_Wall and rh_Background+FreeSurfer_Defined_Medial_Wall
  parcel_df <- rbind(parcel_df, data.frame(region = "rh_Background+FreeSurfer_Defined_Medial_Wall",
                              GAM.cov.tvalue = 0,
                              GAM.cov.pvalue = 0,
                              ANOVA.cov.pvalue = 0,
                              Effectsize.cov.AdjRsq = 0,
                              Effectsize.cov.partialRsq = 0,
                              Anova.cov.pvalue.fdr = 0,
                              significant.fdr = "Sig",
                              significant_status = "Sig",
                              cov_sig_effect = 0,
                              network1 = "rh_Background+FreeSurfer_Defined_Medial_Wall",
                              network2 = "rh_Background+FreeSurfer_Defined_Medial_Wall"))
  parcel_df <- rbind(parcel_df, data.frame(region = "lh_Background+FreeSurfer_Defined_Medial_Wall",
                                           GAM.cov.tvalue = 0,
                                           GAM.cov.pvalue = 0,
                                           ANOVA.cov.pvalue = 0,
                                           Effectsize.cov.AdjRsq = 0,
                                           Effectsize.cov.partialRsq = 0,
                                           Anova.cov.pvalue.fdr = 0,
                                           significant.fdr = "Sig",
                                           significant_status = "Sig",
                                           cov_sig_effect = 0,
                                           network1 = "lh_Background+FreeSurfer_Defined_Medial_Wall",
                                           network2 = "lh_Background+FreeSurfer_Defined_Medial_Wall"))
   
  # need to make a label column for the medial wall - a weird ggseg thing
  parcel_df$label <- parcel_df$region 
  
  parcel_df <- parcel_df %>% select(-region)
  
  parcel_df[[measure]] <- as.numeric(parcel_df[[measure]])
  
  plot <- ggplot() + geom_brain(data=parcel_df,   
                                atlas=schaefer7_200, 
                                mapping=aes(fill=get(measure), 
                                            colour=significant.fdr, 
                                            size=I(0.3)), 
                                show.legend=TRUE, 
                                position = position_brain(side ~ hemi)) + 
    theme_void() + 
    scale_colour_manual(values = c(Sig = "black", Non_Sig = "grey84"), guide = "none") + 
    scale_fill_gradientn(colors= c("#C75DAAFF", "#FFFFFF","#00A3A7FF"), na.value="red", 
                         limits = c(min(parcel_df[[measure]]), max(parcel_df[[measure]])), 
                         oob = squish,
                         values=rescale(c(min(parcel_df[[measure]], na.rm=T), 0, max(parcel_df[[measure]], na.rm=T)))) + ggtitle(dataset) +
    labs(fill = legend_title) + 
    theme(plot.title=element_text(hjust=0.5))
  return(plot)
}


# visualize schaefer 17 network on the cortical surface for reference
visualize_schaefer17_network <- function(df, network_to_plot, color1) {
  df <- df %>% mutate(network_of_interest = ifelse(network == network_to_plot, "Yes", "No"))
  
  ggplot() + geom_brain(data=df, atlas=schaefer17_200, 
                        mapping=aes(fill= network_of_interest,
                                    size=I(0.3)), position = position_brain(side ~ hemi)) + theme_void() + theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    scale_fill_manual(values = c(Yes = color1, No = "white"), guide = "none") + ggtitle(network_to_plot)
}


