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

# "#C75DAAFF", "#FFFFFF","#00A3A7FF"
######################################
# set color palette
######################################
red_teal_white_base <- c( "#A51122FF","#E14C21FF", "#F8A86BFF", "#FFFFFFFF","#ACDFD3FF", "#0B8CACFF","#324DA0FF")
red_teal_white <- colorRampPalette(red_teal_white_base)

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
      colours = red_teal_white(256),
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
plot_cortex <- function(df, hemi, measure, ylim1, ylim2, atlas, horizontal_layout = NULL) {
  
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
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey70")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"),
          plot.title = element_blank()) 
  
  plot_medial <- ggplot() + 
    geom_brain(data = df, atlas= get(atlas), 
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.8)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos2)) +
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey70")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"),
          plot.title = element_blank())
  if (is.null(horizontal_layout)) {
    return(plot_grid(plot_lateral, plot_medial, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr"))
  } else {
    return(plot_grid(plot_lateral, plot_medial, ncol = 2, rel_heights = c(1, 1), align = "v", axis = "lr"))
  }
  
}

# plot cortex overlap in datasets that are significant
plot_cortex_overlap <- function(df, hemi, atlas, horizontal_layout = NULL) {
  
  if(hemi == "left") {
    cortical_pos1 <- "left lateral" 
    cortical_pos2 <- "left medial"
  } else {
    cortical_pos1 <- "right lateral" 
    cortical_pos2 <- "right medial"
  }
  
  base_plot <- function(pos) {
    ggplot() +
      geom_brain(
        data = df,
        atlas = get(atlas),
        mapping = aes(fill = signif_count_factor, size=I(0.8)),
        show.legend = FALSE,
        hemi = hemi,
        position = position_brain(pos)
      ) +
      scale_fill_manual(values = count_colors, name = "Number of Datasets",  na.value = "white") +
      theme_void() +
      theme(legend.position = "right",
            plot.margin = unit(c(0.1, -1, 0.1, -1), "cm"))
  }
  
  plot_lateral <- base_plot(cortical_pos1)
  plot_medial  <- base_plot(cortical_pos2)
  
  if (is.null(horizontal_layout)) {
    return(plot_grid(plot_lateral, plot_medial, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr"))
  } else {
    return(plot_grid(plot_lateral, plot_medial, ncol = 2, rel_heights = c(1, 1), align = "v", axis = "lr"))
  }
  
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


# make network pair plot 
make_networkpair17_plot <- function(df, dataset, ylim1, ylim2, fontsize, fontsize_label, measure) {
  
  if (dataset == "All Datasets") {
    # pretty labels
    lab_map <- c(
      visualCentral="Visual Central", visualPeripheral="Visual Peripheral", 
      somatomotorA="Somatomotor A", somatomotorBAuditory="Somatomotor B Auditory", 
      dorsalAttentionA="Dorsal Attn A", dorsalAttentionB="Dorsal Attn B",
      salienceVentralAttentionA="Salience Ventral Attn A", salienceVentralAttentionB="Salience Ventral Attn B", 
      limbicOrbitofrontal="Limbic OFC", limbicTemporopolar="Limbic TempPolar",
      temporoparietal="Temporoparietal",
      frontoparietalControlA="Frontoparietal A", frontoparietalControlB="Frontoparietal B", frontoparietalControlC="Frontoparietal C", 
      defaultA="Default A", defaultB="Default B", defaultC="Default C"
    )
    
    df <- df |>
      dplyr::mutate(
        network1_label = lab_map[network1],
        network2_label = lab_map[network2]
      )
    
    # desired order
    net_levels <- c(
      "Default A","Default B","Default C", "Dorsal Attn A", "Dorsal Attn B",
      "Frontoparietal A", "Frontoparietal B","Frontoparietal C",
      "Limbic OFC", "Limbic TempPolar", "Salience Ventral Attn A", "Salience Ventral Attn B", 
      "Temporoparietal", "Somatomotor A", "Somatomotor B Auditory", "Visual Central", "Visual Peripheral"
    )
    
  } else {
    # pretty labels
    lab_map <- c(
      visualCentral="VisCent", visualPeripheral="VisPeriph", 
      somatomotorA="SM A", somatomotorBAuditory="SM B Aud", 
      dorsalAttentionA="DorsAttn A", dorsalAttentionB="DorsAttn B",
      salienceVentralAttentionA="SalVentAttn A", salienceVentralAttentionB="SalVentAttn B", 
      limbicOrbitofrontal="LimbOFC", limbicTemporopolar="LimbTempPolar",
      temporoparietal="TempParietal",
      frontoparietalControlA="FPN A", frontoparietalControlB="FPN B", frontoparietalControlC="FPN C", 
      defaultA="DMN A", defaultB="DMN B", defaultC="DMN C"
    )
    
    df <- df |>
      dplyr::mutate(
        network1_label = lab_map[network1],
        network2_label = lab_map[network2]
      )
    
    # desired order
    net_levels <- c(
      "DMN A","DMN B","DMN C", "DorsAttn A", "DorsAttn B",
      "FPN A", "FPN B","FPN C",
      "LimbOFC", "LimbTempPolar", "SalVentAttn A", "SalVentAttn B", 
      "TempParietal", "SM A", "SM B Aud", "VisCent", "VisPeriph"
    )
  }
  
  # ---- ensure plotting in bottom-right triangle ----
  dfp <- df |>
    dplyr::rowwise() |>
    dplyr::mutate(
      n1_index = match(network1_label, net_levels),
      n2_index = match(network2_label, net_levels),
      # flip orientation so the filled cells appear bottom-right of diagonal
      plot_x = ifelse(n1_index < n2_index, network2_label, network1_label),
      plot_y = ifelse(n1_index < n2_index, network1_label, network2_label)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      x = factor(plot_x, levels = net_levels),
      y = factor(plot_y, levels = net_levels)
    )
  
  # ---- plot ----
  ggplot(dfp, aes(x = x, y = y, fill = get(measure))) +
    geom_tile() +
    geom_text(aes(label = dplyr::case_when(
      Anova.cov.pvalue.fdr < 1e-4 ~ "***",
      Anova.cov.pvalue.fdr < 1e-3 ~ "**",
      Anova.cov.pvalue.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )), colour = "black", size = fontsize_label) +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +
    coord_fixed() +
    theme_classic(base_size = fontsize) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
      axis.text.y = element_text(size = fontsize),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 24),
      legend.position = "none",
      plot.margin = margin(0.1, -0.1, 0.5, -0.1, "cm")
    ) +
    scale_fill_gradientn(
      colours = red_teal_white(256),
      values = scales::rescale(c(ylim1, 0, ylim2)),
      limits = c(ylim1, ylim2),
      oob = scales::squish,
      na.value = "white"
    ) +
    ggtitle(dataset)
}

# make network pair overlap plot
plot_networkpair_overlap <- function(df, num_networks = 17, fontsize = 12, fontsize_label = 4) {
  
  # use your lab_map and net_levels logic from make_networkpair17_plot()
  lab_map <- c(
    visualCentral="VisCent", visualPeripheral="VisPeriph", 
    somatomotorA="SM A", somatomotorBAuditory="SM B Aud", 
    dorsalAttentionA="DorsAttn A", dorsalAttentionB="DorsAttn B",
    salienceVentralAttentionA="SalVentAttn A", salienceVentralAttentionB="SalVentAttn B", 
    limbicOrbitofrontal="LimbOFC", limbicTemporopolar="LimbTempPolar",
    temporoparietal="TempParietal",
    frontoparietalControlA="FPN A", frontoparietalControlB="FPN B", frontoparietalControlC="FPN C", 
    defaultA="DMN A", defaultB="DMN B", defaultC="DMN C"
  )
  
  net_levels <- c(
    "DMN A","DMN B","DMN C", "DorsAttn A", "DorsAttn B",
    "FPN A", "FPN B","FPN C",
    "LimbOFC", "LimbTempPolar", "SalVentAttn A", "SalVentAttn B", 
    "TempParietal", "SM A", "SM B Aud", "VisCent", "VisPeriph"
  )
  
  dfp <- df %>%
    mutate(
      network1_label = lab_map[network1],
      network2_label = lab_map[network2]
    ) %>%
    rowwise() %>%
    mutate(
      n1_index = match(network1_label, net_levels),
      n2_index = match(network2_label, net_levels),
      plot_x = ifelse(n1_index < n2_index, network2_label, network1_label),
      plot_y = ifelse(n1_index < n2_index, network1_label, network2_label)
    ) %>%
    ungroup() %>%
    mutate(
      x = factor(plot_x, levels = net_levels),
      y = factor(plot_y, levels = net_levels)
    )
  
  ggplot(dfp, aes(x = x, y = y, fill = signif_count_factor)) +
    geom_tile() +
    scale_fill_manual(values = count_colors, na.value = "white") +
    coord_fixed() +
    theme_classic(base_size = fontsize) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
      axis.text.y = element_text(size = fontsize),
      axis.title = element_blank(),
      legend.position = "none"
    )
}

# edge level analysis: 
# make a double df for symmetry - parcel1 and parcel2_SA repeated in opposite ordering (see lines 500-510 in Adam's Edge-level-Age.md)
make_double_df <- function(df, fc_measure, covariate, dataset) {
  age_SA.diff <- cbind(df, SA.diff)
  age_SA.diff <- sig.effects(age_SA.diff, fc_measure, covariate, dataset) 
  
  # make a double df for symmetry - parcel1 and parcel2_SA repeated in opposite ordering (see lines 500-510 in Adam's Edge-level-Age.md)
  age_SA.diff2 <- age_SA.diff 
  age_SA.diff2$parcel1_SA <- age_SA.diff$parcel2_SA
  age_SA.diff2$parcel2_SA <- age_SA.diff$parcel1_SA
  
  # "stacked" df
  double_age_SA.diff <- rbind(age_SA.diff2,age_SA.diff)
  double_age_SA.diff$parcel2_SA <- as.numeric(double_age_SA.diff$parcel2_SA)
  double_age_SA.diff$parcel1_SA <- as.numeric(double_age_SA.diff$parcel1_SA)
  double_age_SA.diff$SA_diff <- as.numeric(double_age_SA.diff$SA_diff)
  
  return(double_age_SA.diff)
}

 

# model the surface
make_surface_plot <- function(double_SA.diff, dataset, measure, ylim1, ylim2, p_label, legend = NULL, significant_only = FALSE) {
  # fit tensor model using ALL edges
  model_formula <- as.formula(paste0(measure, " ~ te(parcel2_SA, parcel1_SA, k = 3)"))
  g2 <- gam(model_formula, data = double_SA.diff)
  
 
  if (!significant_only) {
      surface.plot.fig <- gg_tensor(g2)
  } else {
    
    # Significant-only plot: show only the edges that survive FDR correction
    tensor_terms <- predict(g2, newdata = double_SA.diff, type = "terms")
    # Identify the tensor smooth column
    tensor_term_name <- grep("^te\\(", colnames(tensor_terms), value = TRUE)
    if (length(tensor_term_name) != 1) {
      stop("Could not uniquely identify the tensor smooth term in the GAM.")
    }
    
    # add the tensor smooth value to the original edge data
    plot_data <- double_SA.diff
    plot_data$tensor_fit <- tensor_terms[, tensor_term_name]
    
    # keep only significant edges for the colored raster
    sig_data <- plot_data[!is.na(plot_data$significant.fdr) & plot_data$significant.fdr == 1,]
    
    # build plot
    tensor_surface <- pammtools::tidy_smooth2d(g2, ci = FALSE)
    surface.plot.fig <- ggplot() + geom_raster(data = sig_data, aes(x = parcel2_SA, y = parcel1_SA, fill = tensor_fit)) +
      geom_contour(data = tensor_surface, aes(x = x, y = y, z = fit), color = "grey30", linewidth = 0.4) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  }
  
  surface.plot.fig <- surface.plot.fig + theme_classic() + 
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = scales::squish, na.value = "white") +
    geom_vline(xintercept = mean(range(double_SA.diff$parcel2_SA, na.rm = TRUE)), linetype = "dashed", size = 1) + 
    geom_hline(yintercept = mean(range(double_SA.diff$parcel1_SA, na.rm = TRUE)), linetype = "dashed", size = 1)
  
  if (is.null(legend)) {
    surface.plot.fig <- surface.plot.fig + labs(title = dataset, subtitle = p_label) +
      theme(legend.position = "none", 
            strip.background = element_blank(), 
            plot.margin = margin(0, 0, -0.1, 0.5, "cm"), 
            axis.text.y = element_text(size = 20, color = "black"), 
            axis.text.x = element_text(size = 20, color = "black"),
            plot.title = element_text(size = 20, hjust = 0.5),
        plot.subtitle = element_text(size = 18, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black"))
    } else {
    surface.plot.fig <- surface.plot.fig + labs(title = dataset, x = "S-A Axis Rank", y = "S-A Axis Rank") +
      theme(legend.position = "bottom", 
            legend.key.height = unit(0.5,"cm"), 
            legend.key.width = unit(1.5, "cm"), 
            legend.title = element_blank(),
            legend.text = element_text(size = 20),
            strip.background = element_blank(),
            plot.margin = margin(0, 0, 1, 0.5, "cm"),
            axis.title = element_text(size = 20),
            axis.text.y = element_text(size = 20, color = "black"),
            axis.text.x = element_text(size = 20, color = "black"),
            plot.title = element_text(size = 20, hjust = 0.5),
            plot.subtitle = element_text(size = 18, hjust = 0.5),
            axis.line = element_line(color = "black"))
    }
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
    scale_fill_gradientn(colours = red_teal_white(256),, na.value="red", 
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



# load pfactor data and merge with demographics
load_pfactor_data <- function(file_path) {
  pfactor <- fread(file_path)
  pfactor <- pfactor %>% rename(sub = participant_id)
  pfactor$sub <- paste0("sub-", pfactor$sub)
  pfactor <- pfactor %>% select(sub, p_factor_mcelroy_harmonized_all_samples, internalizing_mcelroy_harmonized_all_samples, externalizing_mcelroy_harmonized_all_samples, attention_mcelroy_harmonized_all_samples)
  pfactor <- pfactor[!duplicated(pfactor), ]
  pfactor <- pfactor[complete.cases(pfactor), ]
  return(pfactor)
}


# hex plot function for pfactor vs sex edges
hex_plot <- function(df, x, y, text, ylim1, ylim2, xlim1, xlim2, x_text, y_text, x_axis, y_axis, bin_size) {
  plot <- ggplot(df, aes_string(x = x, y = y)) +
    geom_hex(bins = bin_size) +
    paletteer::scale_fill_paletteer_c("grDevices::RdPu", direction = 1) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +   
    annotate("text", x = x_text, y = y_text, 
             label = text, 
             hjust = 0, vjust = 1, size = 6, color = "black") + theme_classic() +
    theme(legend.position = "bottom", 
          legend.key.width = unit(1, 'cm'), legend.key.height = unit(0.5, 'cm'),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.margin = unit(c(1, 0.5, 0.2, 1), "cm")) +
    labs(x = x_axis, y = y_axis, fill = "Count") + ylim(ylim1, ylim2) + xlim(xlim1, xlim2)
  return(plot)
}

 

# identify overlapping WNC and BNC regions for sex effect and get average RESI effect size
identify_sexeffect_overlap <- function(metric = c("WNC", "BNC"), p_thresh = 0.05) {
  metric <- match.arg(metric)
  
  pnc  <- get(paste0("sex.main.PNC.",  metric))
  hcpd <- get(paste0("sex.main.HCPD.", metric))
  nki  <- get(paste0("sex.main.NKI." , metric))
  hbn  <- get(paste0("sex.main.HBN." , metric))
  
  # pvalues
  df_list_p <- list(
    PNC  = pnc  |> dplyr::select(region, p_PNC  = Anova.cov.pvalue.fdr),
    HCPD = hcpd |> dplyr::select(region, p_HCPD = Anova.cov.pvalue.fdr),
    NKI  = nki  |> dplyr::select(region, p_NKI  = Anova.cov.pvalue.fdr),
    HBN  = hbn  |> dplyr::select(region, p_HBN  = Anova.cov.pvalue.fdr)
  )
  
  # RESI
  df_list_resi <- list(
    PNC  = pnc  |> dplyr::select(region, resi_PNC  = RESI),
    HCPD = hcpd |> dplyr::select(region, resi_HCPD = RESI),
    NKI  = nki  |> dplyr::select(region, resi_NKI  = RESI),
    HBN  = hbn  |> dplyr::select(region, resi_HBN  = RESI)
  )
  
  # combine 
  overlap_p <- Reduce(function(x, y) dplyr::full_join(x, y, by = "region"), df_list_p)
  overlap_resi <- Reduce(function(x, y) dplyr::full_join(x, y, by = "region"), df_list_resi)
  
  # merge
  overlap_df <- dplyr::left_join(overlap_p, overlap_resi, by = "region") %>%
    dplyr::mutate(
      signif_count = rowSums(dplyr::across(dplyr::starts_with("p_"), ~ . < p_thresh), na.rm = TRUE),
      signif_2plus = ifelse(signif_count >= 2, 1, 0),
      mean_resi       = rowMeans(dplyr::across(dplyr::starts_with("resi")), na.rm = TRUE),
      mean_resi_sig2  = ifelse(signif_2plus == 1, mean_resi, NA)
    ) %>%
    dplyr::mutate(
      signif_count_factor = factor(
        signif_count, 
        levels = 0:4,
        labels = as.character(0:4)
      )
    )
  
  return(overlap_df)
}

# developmental trajectories by sex_diff: for male and female  
plot_trajectories_bycov <- function(region_name, df, covariate, y.limits, colorF, colorM, margin){
  
  df_plot <- df %>% filter(region %in% region_name)
  plot <- ggplot(df_plot, aes(x = age, y = fitted, group = interaction(sex, region))) + 
    geom_line(aes(color = sex), linewidth = 2) +
    theme_classic() +
    scale_y_continuous(limits = y.limits) +
    scale_color_manual(values = c("Female" = colorF, "Male" = colorM)) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.title = element_blank(), 
          plot.margin = unit(c(0, margin, 1, margin), "cm")) + ggtitle(gsub("[A-Z0-9]+Networks ", "", gsub("_", " ", region_name)))
  return(plot)  
}
 

# gradients analysis figures
##################
# Global dispersion
##################

plot_global_dispersion <- function(dataset, title) {
  
  outputs_root <- sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/output/%s/gradient_dispersion", dataset)
  resi_outputs_dir <- paste0(outputs_root, "/resi")
  
  global_result <- read.csv(sprintf("%s/%s_resi_sex_dispersion_global.csv", resi_outputs_dir, dataset))
  gam_df <- read.csv(sprintf("%s/dispersion_global_%s.csv", outputs_root, dataset))
  
  gam_df <- gam_df %>% drop_na(sex)
  gam_df$sex <- factor(gam_df$sex, levels = c("Male", "Female"))
  
  p_val <- global_result$p
  
  p_label <- case_when(
    p_val < 0.0001 ~ "p < 0.0001",
    p_val < 0.001  ~ "p < 0.001",
    p_val < 0.01   ~ "p < 0.01",
    p_val < 0.05   ~ "p < 0.05",
    TRUE           ~ "p = N.S."
  )
  
  # adjust dispersion for covariates
  gam_df$resid_adj <- residuals(mgcv::gam(dispersion ~ s(age, k = 3) + meanFD_avgSes, data = gam_df))
  
  p_global <- ggplot(gam_df, aes(x = sex, y = resid_adj, fill = sex, color = sex)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.1, alpha = 0.15, size = 0.8) +
    scale_fill_manual(values = c(Female = "#C53121", Male = "#0B8CACFF")) +
    scale_color_manual(values = c(Female = "#C53121", Male = "#0B8CACFF")) +
    labs(title = title, y = "", x = NULL,
         subtitle = sprintf("RESI = %.3f, %s", global_result$RESI_Xinyu, p_label)) +
    theme_classic(base_size = 20) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size = 20, color = "grey30", hjust = 0.5),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black"),
          axis.title.y = element_blank()) + ylim(-600, 1500)
  
  return(p_global)
}


##################
# Within-network dispersion
##################

plot_within_dispersion <- function(dataset, title) {
  
  resi_outputs_dir <- sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/output/%s/gradient_dispersion/resi", dataset)
  wn_results <- read.csv(sprintf("%s/%s_resi_sex_dispersion_within_network.csv", resi_outputs_dir, dataset))
  
  wn_results$significant <- wn_results$p_fdr < 0.05
  
  net_label_map <- function(x) {
    case_when(
      str_detect(x, "visualCentral") ~ "VisCent",
      str_detect(x, "visualPeripheral") ~ "VisPeri",
      str_detect(x, "somatomotorA") ~ "SM A",
      str_detect(x, "somatomotorBAuditory") ~ "SM B Aud",
      str_detect(x, "dorsalAttentionA") ~ "DorsAttn A",
      str_detect(x, "dorsalAttentionB") ~ "DorsAttn B",
      str_detect(x, "salienceVentralAttentionA") ~ "SalVentAttn A",
      str_detect(x, "salienceVentralAttentionB") ~ "SalVentAttn B",
      str_detect(x, "limbicTemporopolar") ~ "LimbTempPolar",
      str_detect(x, "limbicOrbitofrontal") ~ "LimbOFC",
      str_detect(x, "frontoparietalControlA") ~ "FPN A",
      str_detect(x, "frontoparietalControlB") ~ "FPN B",
      str_detect(x, "frontoparietalControlC") ~ "FPN C",
      str_detect(x, "defaultA") ~ "DMN A",
      str_detect(x, "defaultB") ~ "DMN B",
      str_detect(x, "defaultC") ~ "DMN C",
      str_detect(x, "temporoparietal") ~ "TempParietal",
      TRUE ~ x
    )
  }
  
  wn_plot_df <- wn_results %>%
    mutate(net_label = net_label_map(gsub("^wn_", "", outcome)))
  
  p_within <- ggplot(wn_plot_df, aes(x = RESI_Xinyu, y = reorder(net_label, RESI_Xinyu))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey20") +
    geom_errorbarh(aes(xmin = LCI, xmax = UCI, color = significant), height = 0.18, linewidth = 1.2, alpha = 0.7) +
    geom_point(aes(fill = significant, color = significant), shape = 21, size = 4.5, stroke = 1.2, alpha = 0.7) +
    scale_color_manual(values = c(`TRUE` = "#e87b96", `FALSE` = "grey50")) +
    scale_fill_manual(values = c(`TRUE` = "#e87b96", `FALSE` = "grey50")) +
    labs(title = title,  x = NULL, y = NULL) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black"),
          legend.position = "none") + xlim(-0.25, 0.17)
  
  return(p_within)
}

plot_within_overlap <- function() {
  datasets <- c("PNC", "HCPD", "NKI", "HBN")
  
  net_label_map <- function(x) {
    case_when(
      str_detect(x, "visualCentral") ~ "VisCent",
      str_detect(x, "visualPeripheral") ~ "VisPeri",
      str_detect(x, "somatomotorA") ~ "SM A",
      str_detect(x, "somatomotorBAuditory") ~ "SM B Aud",
      str_detect(x, "dorsalAttentionA") ~ "DorsAttn A",
      str_detect(x, "dorsalAttentionB") ~ "DorsAttn B",
      str_detect(x, "salienceVentralAttentionA") ~ "SalVentAttn A",
      str_detect(x, "salienceVentralAttentionB") ~ "SalVentAttn B",
      str_detect(x, "limbicTemporopolar") ~ "LimbTempPolar",
      str_detect(x, "limbicOrbitofrontal") ~ "LimbOFC",
      str_detect(x, "frontoparietalControlA") ~ "FPN A",
      str_detect(x, "frontoparietalControlB") ~ "FPN B",
      str_detect(x, "frontoparietalControlC") ~ "FPN C",
      str_detect(x, "defaultA") ~ "DMN A",
      str_detect(x, "defaultB") ~ "DMN B",
      str_detect(x, "defaultC") ~ "DMN C",
      str_detect(x, "temporoparietal") ~ "TempParietal",
      TRUE ~ x
    )
  }
  
  df <- bind_rows(lapply(datasets, function(d) {
    path <- sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/output/%s/gradient_dispersion/resi/%s_resi_sex_dispersion_within_network.csv", d, d)
    read.csv(path) %>% mutate(dataset = d)
  })) %>%
    mutate(net_label = net_label_map(gsub("^wn_", "", outcome)),
           significant = p_fdr < 0.05) %>%
    group_by(net_label) %>%
    summarise(n_sig = sum(significant, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(n_sig), net_label) %>%
    mutate(net_label = factor(net_label, levels = net_label))
  
  ggplot(df, aes(x = n_sig, y = net_label)) +
    geom_segment(aes(x = 0, xend = n_sig-0.05, yend = net_label),
                 color = "#ad2a72", linewidth = 1.2) +
    geom_point(shape = 21, size = 5, stroke = 1.2, fill = "#ad2a72", color = "#ad2a72", alpha = 0.8) +
    scale_x_continuous(breaks = 0:4, limits = c(0, 4.2)) +
    labs(title = "Overlap of Significant Effects Across Datasets",
         x = "Number of Datasets", y = NULL) +
    theme_classic(base_size = 20) +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          axis.text.x = element_text(size = 20, color = "black"),
          axis.text.y = element_text(size = 20, color = "black"))
}

 
##################
# Between-network dispersion
##################
plot_between_dispersion <- function(dataset, title, fontsize = 12, fontsize_label = 4) {
  resi_outputs_dir <- sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/output/%s/gradient_dispersion/resi", dataset)
  bn_results <- read.csv(sprintf("%s/%s_resi_sex_dispersion_between_network.csv", resi_outputs_dir, dataset))
  
  net_label_map <- function(x) {
    case_when(
      str_detect(x, "visualCentral") ~ "VisCent",
      str_detect(x, "visualPeripheral") ~ "VisPeriph",
      str_detect(x, "somatomotorA") ~ "SM A",
      str_detect(x, "somatomotorBAuditory") ~ "SM B Aud",
      str_detect(x, "dorsalAttentionA") ~ "DorsAttn A",
      str_detect(x, "dorsalAttentionB") ~ "DorsAttn B",
      str_detect(x, "salienceVentralAttentionA") ~ "SalVentAttn A",
      str_detect(x, "salienceVentralAttentionB") ~ "SalVentAttn B",
      str_detect(x, "limbicTemporopolar") ~ "LimbTempPolar",
      str_detect(x, "limbicOrbitofrontal") ~ "LimbOFC",
      str_detect(x, "frontoparietalControlA") ~ "FPN A",
      str_detect(x, "frontoparietalControlB") ~ "FPN B",
      str_detect(x, "frontoparietalControlC") ~ "FPN C",
      str_detect(x, "defaultA") ~ "DMN A",
      str_detect(x, "defaultB") ~ "DMN B",
      str_detect(x, "defaultC") ~ "DMN C",
      str_detect(x, "temporoparietal") ~ "TempParietal",
      TRUE ~ x
    )
  }
  
  net_levels <- c("DMN A","DMN B","DMN C","DorsAttn A","DorsAttn B","FPN A","FPN B","FPN C",
                  "LimbOFC","LimbTempPolar","SalVentAttn A","SalVentAttn B","TempParietal",
                  "SM A","SM B Aud","VisCent","VisPeriph")
  
  bn_plot_df <- bn_results %>%
    mutate(pair = gsub("^bn_", "", outcome)) %>%
    separate(pair, into = c("net1", "net2"), sep = "__") %>%
    mutate(net1 = net_label_map(net1), net2 = net_label_map(net2),
           net1_ord = match(net1, net_levels), net2_ord = match(net2, net_levels),
           plot_x = if_else(net1_ord < net2_ord, net2, net1),
           plot_y = if_else(net1_ord < net2_ord, net1, net2),
           x = factor(plot_x, levels = net_levels),
           y = factor(plot_y, levels = net_levels),
           sig_label = case_when(p_fdr < 1e-4 ~ "***", p_fdr < 1e-3 ~ "**", p_fdr < 0.05 ~ "*", TRUE ~ ""))
  
  ggplot(bn_plot_df, aes(x = x, y = y, fill = RESI_Xinyu)) +
    geom_tile() +
    geom_text(aes(label = sig_label), colour = "black", size = fontsize_label) +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +
    coord_fixed() +
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(-0.22, 0, 0.2)),
                         limits = c(-0.2, 0.2), oob = squish, na.value = "white") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = fontsize) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
          axis.text.y = element_text(size = fontsize),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(0.1, -0.1, 0.5, -0.1, "cm"))
}


plot_between_overlap <- function(fontsize = 12) {
  datasets <- c("PNC", "HCPD", "NKI", "HBN")
  
  net_label_map <- function(x) {
    case_when(
      str_detect(x, "visualCentral") ~ "VisCent",
      str_detect(x, "visualPeripheral") ~ "VisPeriph",
      str_detect(x, "somatomotorA") ~ "SM A",
      str_detect(x, "somatomotorBAuditory") ~ "SM B Aud",
      str_detect(x, "dorsalAttentionA") ~ "DorsAttn A",
      str_detect(x, "dorsalAttentionB") ~ "DorsAttn B",
      str_detect(x, "salienceVentralAttentionA") ~ "SalVentAttn A",
      str_detect(x, "salienceVentralAttentionB") ~ "SalVentAttn B",
      str_detect(x, "limbicTemporopolar") ~ "LimbTempPolar",
      str_detect(x, "limbicOrbitofrontal") ~ "LimbOFC",
      str_detect(x, "frontoparietalControlA") ~ "FPN A",
      str_detect(x, "frontoparietalControlB") ~ "FPN B",
      str_detect(x, "frontoparietalControlC") ~ "FPN C",
      str_detect(x, "defaultA") ~ "DMN A",
      str_detect(x, "defaultB") ~ "DMN B",
      str_detect(x, "defaultC") ~ "DMN C",
      str_detect(x, "temporoparietal") ~ "TempParietal",
      TRUE ~ x
    )
  }
  
  net_levels <- c("DMN A","DMN B","DMN C","DorsAttn A","DorsAttn B","FPN A","FPN B","FPN C",
                  "LimbOFC","LimbTempPolar","SalVentAttn A","SalVentAttn B","TempParietal",
                  "SM A","SM B Aud","VisCent","VisPeriph")
  
  df <- bind_rows(lapply(datasets, function(d) {
    path <- sprintf("/cbica/projects/network_replication/covariate_analyses/sex_diff/output/%s/gradient_dispersion/resi/%s_resi_sex_dispersion_between_network.csv", d, d)
    read.csv(path) %>% mutate(dataset = d)
  })) %>%
    mutate(pair = gsub("^bn_", "", outcome)) %>%
    separate(pair, into = c("net1", "net2"), sep = "__") %>%
    mutate(net1 = net_label_map(net1), net2 = net_label_map(net2),
           n1_index = match(net1, net_levels), n2_index = match(net2, net_levels),
           plot_x = ifelse(n1_index < n2_index, net2, net1),
           plot_y = ifelse(n1_index < n2_index, net1, net2),
           significant = p_fdr < 0.05) %>%
    group_by(plot_x, plot_y) %>%
    summarise(n_sig = sum(significant, na.rm = TRUE), .groups = "drop") %>%
    mutate(sig_count = case_when(
      n_sig <= 1 ~ "0-1",
      n_sig == 2 ~ "2",
      n_sig == 3 ~ "3",
      n_sig == 4 ~ "4"
    ),
    sig_count = factor(sig_count, levels = c("0-1", "2", "3", "4")),
    x = factor(plot_x, levels = net_levels),
    y = factor(plot_y, levels = net_levels))
  
  count_colors <- c("0-1" = "white","2" = "#F6D0B5","3" = "#D95F8D","4" = "#AD2A72")
  
  ggplot(df, aes(x = x, y = y, fill = sig_count)) +
    geom_tile() +
    scale_fill_manual(values = count_colors, na.value = "white") +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +
    coord_fixed() +
    labs(title = "Overlap of Significant Effects Across Datasets", x = NULL, y = NULL) +
    theme_classic(base_size = fontsize) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
          axis.text.y = element_text(size = fontsize),
          plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(0.1, -0.1, 0.5, -0.1, "cm"))
}
