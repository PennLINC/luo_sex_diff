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
library(pammtools)
library(RColorBrewer)
library(rjson)
library(scales)
library(stringr)
library(tidyr)


# set color palette
red_teal_white_base <- rev(c( "#A51122FF","#E14C21FF", "#F8A86BFF", "#FFFFFFFF","#ACDFD3FF", "#0B8CACFF","#324DA0FF"))
red_teal_white <- colorRampPalette(red_teal_white_base)


# function for formatting network names
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
 
# function for calculating number of significant regions (FDR corrected)
# @param df A dataframe with GAM results
# @param fc_measure, A string for FC measure (i.e. "FC_strength")
# @param covariate, A string for covariate of interest (i.e. "sex_diff")
sig.effects <- function(df, fc_measure, covariate, dataset){
  df$cov.pvalue.fdr <- p.adjust(df$GAM.cov.pvalue, method = c("fdr")) 
  df$significant.fdr <- df$cov.pvalue.fdr < 0.05
  df$significant.fdr[df$significant.fdr == TRUE] <- 1
  df$significant.fdr[df$significant.fdr == FALSE] <- 0
  df <- df %>% mutate(significant_status = ifelse(significant.fdr == 1, "Sig", "Nonsig"))
  df <- df %>% mutate(cov_sig_effect = ifelse(significant.fdr == 1, GAM.cov.tvalue, NA))
  sigeffect.totaln <- sum(df$significant.fdr)
  sigeffect.percent <- round(sigeffect.totaln/length(df$significant.fdr)*100, 2)
  print(sprintf("%s: %s regions (%s percent) show a significant association between %s and %s", dataset, sigeffect.totaln, sigeffect.percent, fc_measure, covariate))
  return(df)
}
 
# function for calculating number of significant regions (FDR corrected) - for analysis with pooled datasets
sig.effects.all_datasets <- function(df, fc_measure, covariate, dataset, pvalue_colname, effect_colname){
  df$int.pvalue.fdr <- p.adjust(df[[pvalue_colname]], method = c("fdr")) 
  df$significant.fdr <- df$int.pvalue.fdr < 0.05
  df$significant.fdr[df$significant.fdr == TRUE] <- 1
  df$significant.fdr[df$significant.fdr == FALSE] <- 0
  df <- df %>% mutate(significant_status = ifelse(significant.fdr == 1, "Sig", "Nonsig"))
  df <- df %>% mutate(int_sig_effect = ifelse(significant.fdr == 1, .[[effect_colname]], NA))
  sigeffect.totaln <- sum(df$significant.fdr)
  sigeffect.percent <- round(sigeffect.totaln/length(df$significant.fdr)*100, 2)
  print(sprintf("%s: %s regions (%s percent) show a significant association between %s and %s", dataset, sigeffect.totaln, sigeffect.percent, fc_measure, covariate))
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
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.3)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos1)) +  
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey70")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -2, 0.1, -2), "cm"),
          plot.title = element_blank()) 
  
  plot_medial <- ggplot() + 
    geom_brain(data = df, atlas= get(atlas), 
               mapping=aes(fill=get(measure), colour=significant_status, size=I(0.3)), 
               show.legend=TRUE, 
               hemi = hemi,
               position = position_brain(cortical_pos2)) +
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = squish, na.value = "white") + 
    scale_colour_manual(values = c("Sig" = "black", "Nonsig" = "grey70")) + 
    theme_void() +
    theme(legend.position = "none",
          legend.title = element_blank(),
          plot.margin = unit(c(0.1, -2, 0.1, -2), "cm"),
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
        mapping = aes(fill = signif_count_factor, size=I(0.3)), 
        colour = "black",
        show.legend = FALSE,
        hemi = hemi,
        position = position_brain(pos)
      ) +
      scale_fill_manual(values = count_colors, name = "Number of Datasets",  na.value = "white") +
      theme_void() +
      theme(legend.position = "right",
            plot.margin = unit(c(0.1, -0.5, 0.1, -0.5), "cm"))
  }
  
  plot_lateral <- base_plot(cortical_pos1)
  plot_medial  <- base_plot(cortical_pos2)
  
  if (is.null(horizontal_layout)) {
    return(plot_grid(plot_lateral, plot_medial, ncol = 1, rel_heights = c(1, 1), align = "v", axis = "lr"))
  } else {
    return(plot_grid(plot_lateral, plot_medial, ncol = 2, rel_heights = c(1, 1), align = "v", axis = "lr"))
  }
  
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
      cov.pvalue.fdr < 1e-4 ~ "***",
      cov.pvalue.fdr < 1e-3 ~ "**",
      cov.pvalue.fdr < 0.05 ~ "*",
      TRUE ~ ""
    )), colour = "black", size = fontsize_label, nudge_y = -0.2) +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +
    coord_fixed(ratio = 0.7) +
    theme_classic(base_size = fontsize) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
      axis.text.y = element_text(size = fontsize),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "none",
      plot.margin = margin(0.1, -0.1, 0, -0.1, "cm")
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
plot_networkpair_overlap <- function(df, num_networks = 17, fontsize = 7, fontsize_label = 2.6) {
  
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
    coord_fixed(0.7) +
    theme_classic(base_size = fontsize) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
      axis.text.y = element_text(size = fontsize),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
      legend.position = "none",
      plot.margin = margin(0.1, -0.1, 0, -0.1, "cm")
    ) + ggtitle("Overlap of Significant Effects")
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
    tensor_surface <- pammtools::tidy_smooth2d(g2, ci = FALSE)
    
    surface.plot.fig <- ggplot(tensor_surface, aes(x = x, y = y)) +
      geom_raster(aes(fill = fit)) +
      geom_contour(aes(z = fit), color = "grey30", linewidth = 0.25) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0))
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
      geom_contour(data = tensor_surface, aes(x = x, y = y, z = fit), color = "grey30", linewidth = 0.01) +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))
  }
  
  surface.plot.fig <- surface.plot.fig + theme_classic() + 
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(ylim1, 0, ylim2)), limits = c(ylim1, ylim2), oob = scales::squish, na.value = "white") +
    geom_vline(xintercept = mean(range(double_SA.diff$parcel2_SA, na.rm = TRUE)), linetype = "dashed", size = 0.25) + 
    geom_hline(yintercept = mean(range(double_SA.diff$parcel1_SA, na.rm = TRUE)), linetype = "dashed", size = 0.25)
  
  if (is.null(legend)) {
    surface.plot.fig <- surface.plot.fig + labs(title = dataset, subtitle = p_label) +
      theme(legend.position = "none", 
            strip.background = element_blank(), 
            strip.text = element_blank(),
            plot.margin = margin(0, 0.2, -0.1, 0.15, "cm"), 
            axis.text.y = element_text(size = 9, color = "black"), 
            axis.text.x = element_text(size = 9, color = "black"),
            plot.title = element_text(size = 9, hjust = 0.5, face = "bold", margin = margin(b = -0.05)),
        plot.subtitle = element_text(size = 9, hjust = 0.5, margin = margin(t = -0.05, b = -1)),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25))
    } else {
    surface.plot.fig <- surface.plot.fig + labs(title = dataset, x = "S-A Axis Rank", y = "S-A Axis Rank") +
      theme(legend.position = "bottom", 
            legend.key.height = unit(0.5,"cm"), 
            legend.key.width = unit(1.5, "cm"), 
            legend.title = element_blank(),
            legend.text = element_text(size = 9),
            strip.background = element_blank(),
            strip.text = element_blank(),
            plot.margin = margin(0, 0, 1, 0.5, "cm"),
            axis.title = element_text(size = 9),
            axis.text.y = element_text(size = 9, color = "black"),
            axis.text.x = element_text(size = 9, color = "black"),
            plot.title = element_text(size = 9, hjust = 0.5),
            plot.subtitle = element_text(size = 9, hjust = 0.5),
            axis.line = element_line(color = "black"))
    }
  return(surface.plot.fig)
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
  


# Between-network dispersion
plot_between_dispersion <- function(dataset, title, fontsize = 7, fontsize_label = 2.6) {
  resi_outputs_dir <- sprintf("%s/covariate_analyses/sex_diff/output/%s/gradient_dispersion/resi", root, dataset)
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
    geom_text(aes(label = sig_label), colour = "black", size = fontsize_label, nudge_y = -0.2) +
    scale_x_discrete(limits = net_levels, drop = FALSE) +
    scale_y_discrete(limits = net_levels, drop = FALSE) +
    coord_fixed(ratio = 0.7) +
    scale_fill_gradientn(colours = red_teal_white(256), values = scales::rescale(c(-0.22, 0, 0.2)),
                         limits = c(-0.2, 0.2), oob = squish, na.value = "white") +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = fontsize) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
          axis.text.y = element_text(size = fontsize),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(0.1, -0.1, 0, -0.1, "cm"))
}


# plot between network dispersion overlap
plot_between_overlap <- function(fontsize = 7) {
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
    path <- sprintf("%s/covariate_analyses/sex_diff/output/%s/gradient_dispersion/resi/%s_resi_sex_dispersion_between_network.csv", root, d, d)
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
    coord_fixed(ratio = 0.7) +
    labs(title = "Overlap of Significant Effects", x = NULL, y = NULL) +
    theme_classic(base_size = fontsize) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = fontsize),
          axis.text.y = element_text(size = fontsize),
          plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
          legend.position = "none",
          panel.grid = element_blank(),
          plot.margin = margin(0.1, -0.1, 0, -0.1, "cm"))
}
