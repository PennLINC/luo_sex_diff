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

################
# Supp Fig. 1
################
# order each RESI map to match the parcel order of the rotation file (schaefer200x17.labels).
# merge() sorts rows alphabetically by region, which is NOT the atlas order, so this step matters for the spin test.
get_resi_map <- function(ds, measure) {
  df <- get(sprintf("sex.main.%s.%s", ds, measure))
  idx <- match(schaefer200x17.labels, df$region)
  if (any(is.na(idx))) idx <- match(paste0("17", schaefer200x17.labels), df$region)  # handles the "17" prefix
  stopifnot(!any(is.na(idx)))
  df$RESI[idx]
}

make_sexeffect_corr <- function(measure) {
  maps <- sapply(datasets4, get_resi_map, measure = measure)   # 200 x 4 matrix, atlas order
  corr <- round(cor(maps, method = "pearson"), 5)
  
  p.mat <- matrix(0, 4, 4, dimnames = list(names(datasets4), names(datasets4)))
  for (i in 1:3) for (j in (i + 1):4) {
    p.mat[i, j] <- p.mat[j, i] <- perm.sphere.p(maps[, i], maps[, j], perm.id.full_schaefer200x17, corr.type = "pearson")
  }
  list(corr = corr, p.mat = p.mat)
}

plot_sexeffect_corr <- function(res, show_legend = FALSE) {
  ggcorrplot(res$corr, lab = TRUE, lab_col = "black", type = "lower",
             outline.col = "white", ggtheme = ggplot2::theme_void,
             lab_size = 4, hc.order = FALSE, insig = "blank", p.mat = res$p.mat) +
    scale_fill_gradient2(low = "white", mid = "#E5B9ADFF", high = "#A50026FF",
                         midpoint = 0.5, limit = c(0, 1), breaks = c(0, 0.5, 1),
                         guide = guide_colorbar(title = "Correlation", title.position = "bottom",
                                                title.hjust = 0.5)) +
    theme(axis.text = element_text(size = 10, color = "black"),
          axis.text.y = element_text(hjust = 1),
          axis.text.x = element_text(hjust = 0.5, angle = 0),
          legend.position = if (show_legend) "bottom" else "none",
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.key.height = unit(0.4, "cm"),
          legend.key.width = unit(1, "cm"),
          plot.margin = margin(0, 0.2, 0, 0.2, "cm"))
}

# header: panel letter on the left, title centered, same baseline
panel_header <- function(letter, title) {
  ggdraw() +
    draw_label(letter, x = 0.02, y = 0.95, hjust = 0, vjust = 1, fontface = "bold", size = 10) +
    draw_label(title, x = 0.55, y = 0.95, hjust = 0.5, vjust = 1, fontface = "bold", size = 10,
               lineheight = 0.9)
}


################
# Supp Fig. 2
################

# identify overlapping WNC and BNC regions for sex effect and get average RESI effect size
identify_sexeffect_overlap <- function(metric = c("WNC", "BNC"), p_thresh = 0.05) {
  metric <- match.arg(metric)
  
  pnc  <- get(paste0("sex.main.PNC.",  metric))
  hcpd <- get(paste0("sex.main.HCPD.", metric))
  nki  <- get(paste0("sex.main.NKI." , metric))
  hbn  <- get(paste0("sex.main.HBN." , metric))
  
  # pvalues
  df_list_p <- list(
    PNC  = pnc  |> dplyr::select(region, p_PNC  = cov.pvalue.fdr),
    HCPD = hcpd |> dplyr::select(region, p_HCPD = cov.pvalue.fdr),
    NKI  = nki  |> dplyr::select(region, p_NKI  = cov.pvalue.fdr),
    HBN  = hbn  |> dplyr::select(region, p_HBN  = cov.pvalue.fdr)
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


################
# Supp Fig. 3 and 4
################
# Schaefer region label -> network key (your mapping)
region_to_network <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "VisCent")          ~ "visualCentral",
    stringr::str_detect(x, "VisPeri")          ~ "visualPeripheral",
    stringr::str_detect(x, "SomMotA")          ~ "somatomotorA",
    stringr::str_detect(x, "SomMotB")          ~ "somatomotorBAuditory",
    stringr::str_detect(x, "DorsAttnA")        ~ "dorsalAttentionA",
    stringr::str_detect(x, "DorsAttnB")        ~ "dorsalAttentionB",
    stringr::str_detect(x, "SalVentAttnA")     ~ "salienceVentralAttentionA",
    stringr::str_detect(x, "SalVentAttnB")     ~ "salienceVentralAttentionB",
    stringr::str_detect(x, "LimbicA_TempPole") ~ "limbicTemporopolar",
    stringr::str_detect(x, "LimbicB_OFC")      ~ "limbicOrbitofrontal",
    stringr::str_detect(x, "ContA")            ~ "frontoparietalControlA",
    stringr::str_detect(x, "ContB")            ~ "frontoparietalControlB",
    stringr::str_detect(x, "ContC")            ~ "frontoparietalControlC",
    stringr::str_detect(x, "DefaultA")         ~ "defaultA",
    stringr::str_detect(x, "DefaultB")         ~ "defaultB",
    stringr::str_detect(x, "DefaultC")         ~ "defaultC",
    stringr::str_detect(x, "TempPar")          ~ "temporoparietal",
    TRUE                                       ~ x
  )
}

# "17Networks_LH_DefaultA_pCunPCC_1" -> "LH default mode A\npCunPCC 1"
pretty_region <- function(region) {
  parts <- strsplit(sub("^17Networks_", "", region), "_")[[1]]
  hemi <- parts[1]
  net_label <- unname(pretty_network[region_to_network(region)])
  sub_region <- paste(parts[-(1:2)], collapse = " ")
  paste0(hemi, " ", net_label, "\n", sub_region)
}

# panel order: association -> sensorimotor, using the same network keys
network_order <- c("defaultA", "defaultB", "defaultC",
                   "frontoparietalControlA", "frontoparietalControlB", "frontoparietalControlC",
                   "dorsalAttentionA", "dorsalAttentionB",
                   "salienceVentralAttentionA", "salienceVentralAttentionB",
                   "limbicOrbitofrontal", "limbicTemporopolar", "temporoparietal",
                   "somatomotorA", "somatomotorBAuditory", "visualCentral", "visualPeripheral")
pretty_network <- c(
  defaultA = "Default Mode A", defaultB = "Default Mode B", defaultC = "Default Mode C",
  frontoparietalControlA = "Control A", frontoparietalControlB = "Control B", frontoparietalControlC = "Control C",
  dorsalAttentionA = "Dorsal Attn A", dorsalAttentionB = "Dorsal Attention B",
  salienceVentralAttentionA = "Salience/Ventral Attn A", salienceVentralAttentionB = "Salience/Ventral Attn B",
  limbicOrbitofrontal = "Limbic Orbitofrontal", limbicTemporopolar = "Limbic Temporal Pole",
  somatomotorA = "Somatomotor A", somatomotorBAuditory = "Somatomotor B/Auditory",
  visualCentral = "Visual Central", visualPeripheral = "Visual Peripheral", temporoparietal = "Temporal Parietal"
)

pretty_pair <- function(region) {
  nets <- strsplit(region, "_")[[1]]
  nets <- ifelse(nets %in% names(pretty_network), unname(pretty_network[nets]), nets)
  if (nets[1] == nets[2]) paste(nets[1], "within") else paste(nets[1], "with", nets[2])
}

# "17Networks_LH_DefaultA_pCunPCC_1" -> "LH default mode A\npCunPCC 1"
pretty_region <- function(region) {
  parts <- strsplit(sub("^17Networks_", "", region), "_")[[1]]
  hemi <- parts[1]
  net_label <- unname(pretty_network[region_to_network(region)])
  sub_region <- paste(parts[-(1:2)], collapse = " ")
  paste0(hemi, " ", net_label, "\n", sub_region)
}

plot_traj <- function(reg, df, title, y.limits) {
  ggplot(df %>% filter(region == reg), aes(x = age, y = fitted, color = sex)) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = c(Female = colorF, Male = colorM)) +
    coord_cartesian(ylim = y.limits) +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 7) +
    theme(legend.position = "none",
          plot.title = element_text(size = 7, hjust = 0.5, lineheight = 0.9),
          axis.text = element_text(size = 7, color = "black"),
          axis.line = element_line(linewidth = 0.3),
          axis.ticks = element_line(linewidth = 0.3),
          plot.margin = margin(0.1, 0.15, 0.1, 0.15, "cm"))
}

assemble_traj <- function(plots, ncol, y_label, legend_height = 0.05) {
  grid <- plot_grid(plotlist = plots, ncol = ncol, align = "hv")
  grid <- annotate_figure(grid,
                          left = text_grob(y_label, size = 10, rot = 90),
                          bottom = text_grob("Age (years)", size = 10))
  plot_grid(grid, traj_legend, nrow = 2, rel_heights = c(1, legend_height))
}

add_17 <- function(region) if_else(startsWith(region, "17"), region, paste0("17", region))
# networks ordered from association to sensorimotor (for panel order)
network_order_17 <- c("DefaultA", "DefaultB", "DefaultC", "ContA", "ContB", "ContC", "DorsAttnA", "DorsAttnB",
                      "SalVentAttnA", "SalVentAttnB", "LimbicA", "LimbicB", "TempPar",
                      "SomMotA", "SomMotB", "VisCent", "VisPeri")


# ---- one hemisphere: lateral + medial ----
plot_cortex_pc1 <- function(df, hemi, sex_name, horizontal = FALSE) {
  df <- df %>% filter(sex == sex_name)
  views <- if (hemi == "left") c("left lateral", "left medial") else c("right lateral", "right medial")
  
  one_view <- function(view) {
    ggplot() +
      geom_brain(data = df, atlas = schaefer17_200, mapping = aes(fill = pc1_mean),
                 hemi = hemi, position = position_brain(view), color = "grey40", size = 0.1,
                 show.legend = FALSE) +
      pc1_fill +
      theme_void() +
      theme(plot.margin = unit(c(0, -0.3, 0, -0.3), "cm"))
  }
  
  plot_grid(one_view(views[1]), one_view(views[2]),
            ncol = if (horizontal) 2 else 1, align = "v", axis = "lr")
}


################
# Supp Fig. 5
################
# ---- female / male block for one dataset ----
sex_block <- function(df, sex_name, horizontal = FALSE) {
  brains <- plot_grid(plot_cortex_pc1(df, "left", sex_name, horizontal),
                      plot_cortex_pc1(df, "right", sex_name, horizontal), ncol = 2)
  plot_grid(ggdraw() + draw_label(sex_name, size = 12), brains,
            ncol = 1, rel_heights = c(0.12, 1))
}

plot_pc1_dataset <- function(dataset, title, horizontal = FALSE) {
  df <- pc1_brain %>% filter(dataset == !!dataset)
  brains <- plot_grid(sex_block(df, "Female", horizontal), NULL, sex_block(df, "Male", horizontal),
                      ncol = 3, rel_widths = c(1, 0.05, 1))
  plot_grid(ggdraw() + draw_label(title, size = 12, fontface = "bold"), brains,
            ncol = 1, rel_heights = c(if (horizontal) 0.18 else 0.1, 1))
}

################
# Supp Fig. 6
################
plot_global_dispersion <- function(dataset, title) {
  
  outputs_root <- sprintf("%s/%s/gradient_dispersion", output_root, dataset)
  resi_outputs_dir <- paste0(outputs_root, "/resi")
  
  global_result <- read.csv(sprintf("%s/%s_resi_sex_dispersion_global.csv", resi_outputs_dir, dataset))
  gam_df <- read.csv(sprintf("%s/dispersion_global_%s.csv", outputs_root, dataset))
  
  gam_df <- gam_df %>% drop_na(sex)
  gam_df$sex <- factor(gam_df$sex, levels = c("Female", "Male"))
  
  p_val <- global_result$p
  p_lab <- if (p_val < 0.0001) {
    " < 0.0001"
  } else {
    paste0(" = ", format(signif(p_val, 2), scientific = FALSE))
  }
  message(sprintf("%s: RESI = %.3f, p = %.3g", dataset, global_result$RESI_Xinyu, p_val))  # exact values for the text
  resi_lab <- sprintf("%.2f", global_result$RESI_Xinyu)
  subtitle_expr <- bquote("RESI = " * .(resi_lab) * ", " * italic(p) * .(p_lab))
  
  # residualize dispersion for age and in-scanner motion
  gam_df$resid_adj <- residuals(mgcv::gam(dispersion ~ s(age, k = 3, fx = TRUE) + meanFD_avgSes,
                                          data = gam_df, method = "REML"))
  
  ggplot(gam_df, aes(x = sex, y = resid_adj, fill = sex, color = sex)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6, linewidth = 0.4) +
    geom_jitter(width = 0.1, alpha = 0.15, size = 0.3) +
    scale_fill_manual(values = c(Female = "#C53121", Male = "#0B8CACFF")) +
    scale_color_manual(values = c(Female = "#C53121", Male = "#0B8CACFF")) +
    labs(title = title, subtitle = subtitle_expr, x = NULL, y = NULL) +
    coord_cartesian(ylim = c(-600, 1500)) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 10, color = "grey30", hjust = 0.5),
          axis.text = element_text(size = 10, color = "black"),
          axis.line = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4),
          plot.margin = margin(0.1, 0.2, 0.1, 0.2, "cm"))
}
