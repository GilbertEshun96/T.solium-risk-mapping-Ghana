# ============================================================================== #
# 1. SETUP: LIBRARIES & PATH DEFINITIONS
# ============================================================================== #
library(sf); library(sp); library(terra); library(raster); library(tidyterra)
library(gstat); library(classInt)
library(dplyr); library(tidyr); library(e1071); library(irr); library(foreign); library(moments)
library(ggplot2); library(ggnewscale); library(ggpubr); library(gridExtra); library(grid); library(patchwork)
library(viridis); library(colorspace)
library(knitr); library(officer); library(flextable)

library(future); library(future.apply)
plan(multisession, workers = 4)

setwd("D:/My studies documents/Taenia Solium/R codes")

# ============================================================================== #
# 2. LOAD SHAPEFILES, LAKES & INITIAL GRID
# ============================================================================== #
shp2 <- st_read("D:/My studies documents/Taenia Solium/R codes/gadm41_GHA_shp/gadm41_GHA_2.shp")
shp_lakes <- st_read("D:/My studies documents/Taenia Solium/R codes/Ghana lakes/Lakes.shp")
shp_lakes <- st_transform(shp_lakes, st_crs(shp2))
sf_use_s2(FALSE)

shp2 <- st_make_valid(shp2)
shp_lakes <- st_make_valid(shp_lakes)
shp2 <- shp2[!st_is_empty(shp2), ]
shp_lakes <- shp_lakes[!st_is_empty(shp_lakes), ]

shp_subnational <- st_difference(shp2, st_union(shp_lakes))
shp_subnational <- st_transform(shp_subnational, 4326)

master_temp_rast <- rast(vect(shp_subnational), res = 0.0083)
crs(master_temp_rast) <- "EPSG:4326"
ghana_vect <- vect(shp_subnational)

# -----------------------------
# FIX: clip lakes/lagoons to Ghana
# -----------------------------
ghana_sf <- sf::st_transform(sf::st_as_sf(shp2), 4326)
lakes_sf <- sf::st_transform(sf::st_as_sf(shp_lakes), 4326)

ghana_sf <- sf::st_make_valid(ghana_sf)
lakes_sf <- sf::st_make_valid(lakes_sf)

if (sf::st_crs(lakes_sf) != sf::st_crs(ghana_sf)) {
  lakes_sf <- sf::st_transform(lakes_sf, sf::st_crs(ghana_sf))
}

ghana_union <- sf::st_union(ghana_sf)

lakes_sf_correct <- suppressWarnings(
  sf::st_intersection(lakes_sf, ghana_union)
)

# remove empty/problem geometries
lakes_sf_correct <- lakes_sf_correct[!sf::st_is_empty(lakes_sf_correct), ]

# keep only polygon geometries
lakes_sf_correct <- suppressWarnings(
  sf::st_collection_extract(lakes_sf_correct, "POLYGON")
)

# ============================================================================== #
# 3. GRID CREATION FUNCTION
# ============================================================================== #
create_grid_and_geodata_sf <- function(shp, cellsize = 0.0083, geodata_sf){
  grid <- st_make_grid(shp, cellsize = cellsize, what = "centers")
  grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)
  st_crs(grid_sf) <- st_crs(shp)
  grid_sf <- st_intersection(grid_sf, shp)
  if(!inherits(geodata_sf,"sf")) stop("geodata_sf must be an sf object")
  list(grid_sf = grid_sf, geodata_obj_sf = geodata_sf)
}

# ============================================================================== #
# 4. YEAR CONFIGURATION
# ============================================================================== #
analysis_years <- list(
  "2022" = list(
    data_path="D:/My studies documents/Taenia Solium/R codes/GHPR8CSV/GHPR8CFL.SAV",
    geo_path="D:/My studies documents/Taenia Solium/R codes/GHGE8AFL/GHGE8AFL.shp",
    pig_col="HV246I",
    glw_path="C:/Users/gilbe/Downloads/GLW4-2020.D-DA.PGS.tif",
    glw_name="Pig_Density_GLW4",
    wc1_codes=c(14,23,31,42,43,96)
  ),
  "2014" = list(
    data_path="D:/My studies documents/Taenia Solium/R codes/Datasets/2014/GHPR72SV/GHPR72FL.SAV",
    geo_path="D:/My studies documents/Taenia Solium/R codes/Datasets/2014/GHGE71FL/GHGE71FL.shp",
    pig_col="HV207G",
    glw_path="D:/My studies documents/Taenia Solium/R codes/Datasets/5_Pg_2015_Da.tif",
    glw_name="Pig_Density_GLW3",
    wc1_codes=c(14,23,31,42,43,96)
  ),
  "2008" = list(
    data_path="D:/My studies documents/Taenia Solium/R codes/Datasets/2008/GHPR5ASV/GHPR5AFL.SAV",
    geo_path="D:/My studies documents/Taenia Solium/R codes/Datasets/2008/GHGE5AFL/GHGE5AFL.shp",
    pig_col="HV207G",
    glw_path="D:/My studies documents/Taenia Solium/R codes/Datasets/5_Pg_2010_Da.tif",
    glw_name="Pig_Density_GLW2",
    wc1_codes=c(14,23,31,42,96)
  ),
  "2003" = list(
    data_path="D:/My studies documents/Taenia Solium/R codes/Datasets/2003/GHHR4BSV/GHHR4BFL.SAV",
    geo_path="D:/My studies documents/Taenia Solium/R codes/Datasets/2003/GHGE4BFL/GHGE4BFL.shp",
    pig_col=NULL,
    glw_path="C:/Users/gilbe/OneDrive/Desktop/My studies documents/Taenia Solium/R codes/Glb_Pigs_CC2006_AD.tif",
    glw_name="Pig_Density_GLW1",
    wc1_codes=c(21,31,32,96)
  )
)

# ============================================================================== #
# 5. MASTER ANALYSIS FUNCTION - CORRECTED FULL VERSION
# ============================================================================== #

run_full_master_analysis_parallel <- function(
    formula, var_column_name, var_label, year_val,
    current_clu, current_geodata_sf, current_grid_sf, current_doc){
  
  year_dir <- paste0("outputs_", year_val, "/")
  dir.create(paste0(year_dir, "raw_maps"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(year_dir, "rasterized_results"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(year_dir, "variograms"), recursive = TRUE, showWarnings = FALSE)
  
  safe_label <- gsub("[^A-Za-z0-9_]", "_", var_label)
  outfile <- paste0(year_dir, "rasterized_results/", safe_label, "_data.tif")
  
  p_raw <- ggplot() +
    geom_sf(data = shp_subnational, fill = "grey90", colour = "grey80") +
    geom_point(data = current_clu, aes(x = x, y = y, color = .data[[var_column_name]])) +
    scale_color_viridis_c(option = "viridis") +
    coord_sf() +
    theme_bw() +
    labs(title = paste(year_val, "Raw:", var_label),
         x = "Longitude", y = "Latitude", colour = var_label)
  
  ggsave(paste0(year_dir, "raw_maps/", safe_label, "_raw_points.png"),
         p_raw, width = 7, height = 6)
  
  geodata_sp <- as(current_geodata_sf, "Spatial")
  
  vgm_emp <- variogram(formula, geodata_sp)
  if (nrow(vgm_emp) == 0) vgm_emp <- data.frame(dist = 1, gamma = 0)
  
  nugget_init <- ifelse(length(vgm_emp$gamma) > 0, vgm_emp$gamma[1], 0)
  total_sill  <- ifelse(length(vgm_emp$gamma) > 0, max(vgm_emp$gamma, na.rm = TRUE), 1)
  psill_init  <- total_sill - nugget_init
  
  idx <- which(vgm_emp$gamma >= 0.95 * total_sill)
  range_init <- if (length(idx) > 0) vgm_emp$dist[min(idx)] else max(vgm_emp$dist)
  range_init <- max(range_init, 0.01)
  
  initial_models <- list(
    Spherical   = vgm(psill = psill_init, model = "Sph", nugget = nugget_init, range = range_init),
    Exponential = vgm(psill = psill_init, model = "Exp", nugget = nugget_init, range = range_init),
    Gaussian    = vgm(psill = psill_init, model = "Gau", nugget = nugget_init, range = range_init),
    Wave        = vgm(psill = psill_init, model = "Wav", nugget = nugget_init, range = range_init)
  )
  
  fit_models <- lapply(
    initial_models,
    function(m) tryCatch(fit.variogram(vgm_emp, model = m), error = function(e) NULL)
  )
  
  plot_variogram_gg <- function(emp_vgm, model_fit, title){
    if (is.null(model_fit)) return(ggplot() + labs(title = paste(title, "(Failed)")))
    
    df_emp <- as.data.frame(emp_vgm)
    model_line <- variogramLine(model_fit, maxdist = max(df_emp$dist))
    
    ggplot() +
      geom_point(data = df_emp, aes(x = dist, y = gamma), color = "blue") +
      geom_line(data = model_line, aes(x = dist, y = gamma), color = "red") +
      labs(x = "Distance", y = "Semivariance", title = title) +
      theme_minimal()
  }
  
  grid_img <- arrangeGrob(
    plot_variogram_gg(vgm_emp, fit_models$Spherical, "Spherical Fit"),
    plot_variogram_gg(vgm_emp, fit_models$Exponential, "Exponential Fit"),
    plot_variogram_gg(vgm_emp, fit_models$Gaussian, "Gaussian Fit"),
    plot_variogram_gg(vgm_emp, fit_models$Wave, "Wave Fit"),
    ncol = 2,
    nrow = 2
  )
  
  img_path <- paste0(year_dir, "variograms/", safe_label, "_variogram_grid.png")
  ggsave(filename = img_path, grid_img, width = 7, height = 6)
  
  sse_calc <- function(emp_vgm, model_fit){
    if (is.null(model_fit)) return(Inf)
    v_line <- variogramLine(model_fit, maxdist = max(emp_vgm$dist))
    preds <- approx(x = v_line$dist, y = v_line$gamma, xout = emp_vgm$dist)$y
    sum((emp_vgm$gamma - preds)^2, na.rm = TRUE)
  }
  
  sse_results <- data.frame(
    Model = names(fit_models),
    Nugget = sapply(fit_models, function(m) if (is.null(m)) NA else m$psill[1]),
    Partial_Sill = sapply(fit_models, function(m) if (is.null(m)) NA else m$psill[2]),
    Range = sapply(fit_models, function(m) if (is.null(m)) NA else m$range[2]),
    SSE = sapply(fit_models, function(m) sse_calc(vgm_emp, m))
  )
  
  best_idx <- which.min(sse_results$SSE)
  best_key <- names(fit_models)[best_idx]
  chosen_fit <- fit_models[[best_key]]
  
  ft <- flextable(sse_results) %>%
    colformat_double(j = c("Nugget", "Partial_Sill", "Range", "SSE"), digits = 6) %>%
    autofit() %>%
    bold(i = best_idx, part = "body")
  
  current_doc <- current_doc %>%
    body_add_par(paste0(year_val, " - Analysis for: ", var_label), style = "heading 1") %>%
    body_add_flextable(ft) %>%
    body_add_par(paste("Chosen variogram model:", best_key), style = "Normal") %>%
    body_add_par(paste("Variogram Model Comparison for", var_label), style = "heading 2") %>%
    body_add_img(src = img_path, width = 6, height = 5) %>%
    body_add_break()
  
  if (file.exists(outfile)) {
    message(paste("Raster already exists. Loading:", outfile))
    r_clipped <- rast(outfile)
  } else {
    grid_sp <- as(current_grid_sf, "Spatial")
    
    if (nrow(geodata_sp) < 2 || all(is.na(geodata_sp[[var_column_name]]))) {
      message("Not enough points for kriging; using IDW fallback")
      kriged <- gstat::idw(formula, geodata_sp, grid_sp, idp = 2.0)
    } else {
      kriged <- tryCatch({
        krige(formula, geodata_sp, grid_sp, model = chosen_fit, nmax = 50, maxdist = Inf)
      }, error = function(e){
        message("Kriging failed; using IDW instead")
        gstat::idw(formula, geodata_sp, grid_sp, idp = 2.0)
      })
    }
    
    v <- vect(st_as_sf(kriged))
    r_kriged <- rasterize(v, master_temp_rast, field = "var1.pred")
    crs(r_kriged) <- crs(master_temp_rast)
    r_clipped <- mask(crop(r_kriged, ghana_vect), ghana_vect)
    
    writeRaster(r_clipped, outfile, overwrite = TRUE)
  }
  
  raster_png <- paste0(year_dir, "rasterized_results/", safe_label, "_raster_map.png")
  
  png(filename = raster_png, width = 2100, height = 1800, res = 300)
  plot(r_clipped, main = paste("Interpolated:", var_label), col = viridis(100))
  plot(ghana_vect, add = TRUE, border = "black")
  dev.off()
  
  current_doc <- current_doc %>%
    body_add_par(paste("Rasterized Kriging Map (Land Only):", var_label), style = "heading 2") %>%
    body_add_img(src = raster_png, width = 6, height = 5) %>%
    body_add_break()
  
  list(raster = r_clipped, doc = current_doc)
}



# ============================================================================== #
# 6. FULL WORKFLOW EXECUTION 
# ============================================================================== #

doc <- read_docx()

all_rasters <- list()
percentiles <- list(`66th` = 0.666, `75th` = 0.75)

for (yr in names(analysis_years)) {
  
  config <- analysis_years[[yr]]
  all_rasters[[yr]] <- list()
  year_dir <- paste0("outputs_", yr, "/")
  dir.create(year_dir, showWarnings = FALSE, recursive = TRUE)
  
  data_yr <- read.spss(config$data_path, to.data.frame = TRUE, use.value.labels = FALSE)
  geo_yr <- st_read(config$geo_path)
  geo_yr <- st_transform(geo_yr, 4326)
  
  data_yr$wc1 <- ifelse(data_yr$HV205 %in% config$wc1_codes, 1, 0)
  data_yr$w1w2 <- ifelse(data_yr$HV270 %in% c(1, 2), 1, 0)
  
  if (yr == "2022") {
    data_yr$wc2 <- ifelse(data_yr$HV205 == 31, 1, 0)
    data_yr$w1  <- ifelse(data_yr$HV270 == 1, 1, 0)
  }
  
  if (!is.null(config$pig_col)) {
    data_yr$Pig_Mean_DHS <- data_yr[[config$pig_col]]
  }
  
  myClu <- geo_yr[geo_yr$LONGNUM != 0 & !is.na(geo_yr$LONGNUM) & !is.na(geo_yr$LATNUM), ]
  myClu <- myClu[, c("DHSCC", "DHSCLUST", "DHSREGCO", "DHSREGNA", "URBAN_RURA",
                     "ADM1NAME", "LONGNUM", "LATNUM", "geometry")]
  names(myClu) <- c("Cc", "HV001", "REGCODE", "HV024", "IsRural", "AdName", "x", "y", "geometry")
  
  vars_to_agg <- c("wc1", "w1w2")
  
  if (yr == "2022") {
    vars_to_agg <- c(vars_to_agg, "wc2", "w1")
  }
  
  if ("Pig_Mean_DHS" %in% names(data_yr)) {
    vars_to_agg <- c(vars_to_agg, "Pig_Mean_DHS")
  }
  
  for (v in vars_to_agg) {
    temp <- aggregate(x = data_yr[[v]], by = list(HV001 = data_yr$HV001), FUN = mean, na.rm = TRUE)
    names(temp)[2] <- v
    myClu <- left_join(myClu, temp, by = "HV001")
  }
  
  Grid_out <- create_grid_and_geodata_sf(shp_subnational, geodata_sf = myClu)
  
  glw_rast <- rast(config$glw_path)
  glw_rast <- project(glw_rast, master_temp_rast)
  glw_rast <- crop(glw_rast, ghana_vect)
  glw_rast <- mask(glw_rast, ghana_vect)
  
  all_rasters[[yr]][[config$glw_name]] <- glw_rast
  all_rasters[[yr]][["Pig_Density_Final"]] <- glw_rast
  
  glw_dir <- paste0("outputs_", yr, "/rasterized_results/")
  dir.create(glw_dir, recursive = TRUE, showWarnings = FALSE)
  
  tif_file <- paste0(glw_dir, config$glw_name, "_Ghana.tif")
  writeRaster(glw_rast, filename = tif_file, overwrite = TRUE)
  
  png_file <- paste0(glw_dir, config$glw_name, "_Ghana.png")
  png(filename = png_file, width = 2100, height = 1800, res = 300)
  plot(glw_rast, main = paste("GLW Pig Density", yr), col = viridis::viridis(100))
  plot(ghana_vect, add = TRUE, border = "black")
  dev.off()
  
  message(paste("Saved GLW raster for", yr, "->", tif_file, "and", png_file))
  
  vars_list <- list(
    list(formula = wc1 ~ 1, col = "wc1", label = "Unimproved San."),
    list(formula = w1w2 ~ 1, col = "w1w2", label = "Poverty W1+W2")
  )
  
  if (yr == "2022") {
    vars_list <- c(
      vars_list,
      list(list(formula = wc2 ~ 1, col = "wc2", label = "Open Defecation")),
      list(list(formula = w1 ~ 1, col = "w1", label = "Wealth Q1"))
    )
  }
  
  if ("Pig_Mean_DHS" %in% names(myClu)) {
    vars_list <- c(
      vars_list,
      list(list(formula = Pig_Mean_DHS ~ 1, col = "Pig_Mean_DHS", label = "Mean Pig (DHS)"))
    )
  }
  
  for (vv in vars_list) {
    res <- run_full_master_analysis_parallel(
      formula = vv$formula,
      var_column_name = vv$col,
      var_label = vv$label,
      year_val = yr,
      current_clu = myClu,
      current_geodata_sf = Grid_out$geodata_obj_sf,
      current_grid_sf = Grid_out$grid_sf,
      current_doc = doc
    )
    
    all_rasters[[yr]][[vv$label]] <- res$raster
    doc <- res$doc
  }
}

dir.create("outputs", showWarnings = FALSE, recursive = TRUE)

print(doc, target = "outputs/Ghana_Variogram_Analysis_FULL.docx")
message("✅ Word document saved: outputs/Ghana_Variogram_Analysis_FULL.docx")



# ============================================================================== #
# 7. HISTOGRAMS, PERCENTILE ANALYSES & HIGH-RISK RASTERS
# ============================================================================== #

variable_map <- c(
  "wc1" = "Unimproved San.",
  "wc2" = "Open Defecation",
  "w1"  = "Wealth Q1",
  "w1w2" = "Poverty W1+W2",
  "Mean Pig (DHS)" = "Mean Pig (DHS)",
  "Pig_Density_GLW1" = "Pig_Density_GLW1",
  "Pig_Density_GLW2" = "Pig_Density_GLW2",
  "Pig_Density_GLW3" = "Pig_Density_GLW3",
  "Pig_Density_GLW4" = "Pig_Density_GLW4"
)

variable_colors <- c(
  "wc1" = "#FFD700",
  "wc2" = "#FFD700",
  "w1"  = "#FF0000",
  "w1w2" = "#FF0000",
  "Mean Pig (DHS)" = "#9400D3",
  "Pig_Density_GLW1" = "#9400D3",
  "Pig_Density_GLW2" = "#9400D3",
  "Pig_Density_GLW3" = "#9400D3",
  "Pig_Density_GLW4" = "#9400D3"
)

threshold_list_66 <- list()
threshold_list_75 <- list()
dir.create("outputs", showWarnings = FALSE)

plot_raster_with_lakes <- function(rast_bin, land_vect, lakes_sf, high_col, main_title = NULL, legend = FALSE) {
  col_vals <- c("#D9D9D9", high_col)
  
  plot(rast_bin, col = col_vals, axes = FALSE, box = FALSE, main = main_title, legend = legend)
  plot(lakes_sf, add = TRUE, col = "lightblue", border = NA)
  plot(land_vect, add = TRUE, border = "black", lwd = 0.5)
}

for (pname in names(percentiles)) {
  
  q_val <- percentiles[[pname]]
  
  png(paste0("outputs/All_Variables_Histograms_", pname, "_FULL.png"),
      width = 4000, height = 4000, res = 300)
  par(mfrow = c(4, 6), mar = c(4, 4, 3, 1))
  
  for (yr in names(all_rasters)) {
    
    vars_this_year <- switch(
      yr,
      "2022" = c("wc1", "wc2", "w1", "w1w2", "Mean Pig (DHS)", "Pig_Density_GLW4"),
      "2014" = c("wc1", "w1w2", "Mean Pig (DHS)", "Pig_Density_GLW3"),
      "2008" = c("wc1", "w1w2", "Mean Pig (DHS)", "Pig_Density_GLW2"),
      "2003" = c("wc1", "w1w2", "Pig_Density_GLW1")
    )
    
    for (label in vars_this_year) {
      
      raster_name <- variable_map[[label]]
      
      if (!(raster_name %in% names(all_rasters[[yr]]))) next
      
      r <- all_rasters[[yr]][[raster_name]]
      if (is.null(r)) next
      
      vals <- values(r, mat = FALSE)
      vals <- vals[!is.na(vals) & vals >= 0]
      if (length(vals) == 0) next
      
      q_real <- quantile(vals, probs = q_val, na.rm = TRUE)
      key <- paste0(yr, "_", label)
      
      if (pname == "66th") {
        threshold_list_66[[key]] <- q_real
      } else {
        threshold_list_75[[key]] <- q_real
      }
      
      is_count <- grepl("Pig", label)
      plot_vals <- if (is_count) log1p(vals) else vals
      hist_xlab <- if (is_count) "Log(Count + 1)" else "Probability / Proportion"
      
      hist(
        plot_vals,
        main = paste(yr, label),
        xlab = hist_xlab,
        col = "gray95",
        border = "white",
        breaks = 40,
        cex.main = 0.7
      )
      
      abline(
        v = if (is_count) log1p(q_real) else q_real,
        col = "red",
        lwd = 2,
        lty = 2
      )
      
      mtext(
        paste(pname, ":", round(q_real, 2)),
        side = 3,
        line = -1,
        adj = 0.95,
        cex = 0.5,
        col = "red"
      )
      
      bin_rast <- r >= q_real
      high_col <- variable_colors[[label]]
      
      file_name_noleg <- paste0(
        "outputs/",
        yr, "_", gsub(" ", "_", label),
        "_HighRisk_", pname, "_FULL_noLegend.png"
      )
      
      png(file_name_noleg, width = 2100, height = 2400, res = 300)
      plot_raster_with_lakes(
        bin_rast,
        ghana_vect,
        lakes_sf_correct,
        high_col,
        main_title = NULL,
        legend = FALSE
      )
      dev.off()
      
      file_name_leg <- paste0(
        "outputs/",
        yr, "_", gsub(" ", "_", label),
        "_HighRisk_", pname, "_FULL_legend.png"
      )
      
      png(file_name_leg, width = 2100, height = 2400, res = 300)
      plot_raster_with_lakes(
        bin_rast,
        ghana_vect,
        lakes_sf_correct,
        high_col,
        main_title = NULL,
        legend = TRUE
      )
      dev.off()
    }
  }
  
  dev.off()
}

threshold_df_66 <- data.frame(
  Variable = names(threshold_list_66),
  Value = as.numeric(unlist(threshold_list_66))
)

threshold_df_75 <- data.frame(
  Variable = names(threshold_list_75),
  Value = as.numeric(unlist(threshold_list_75))
)

doc_thresh <- read_docx() %>%
  body_add_par("Summary of Risk Thresholds", style = "heading 1") %>%
  body_add_par("66.6th Percentile", style = "heading 2") %>%
  body_add_flextable(
    flextable(threshold_df_66) %>%
      colformat_double(j = "Value", digits = 6) %>%
      autofit()
  ) %>%
  body_add_par("75th Percentile", style = "heading 2") %>%
  body_add_flextable(
    flextable(threshold_df_75) %>%
      colformat_double(j = "Value", digits = 6) %>%
      autofit()
  )

print(doc_thresh, target = "outputs/Ghana_Risk_Thresholds_Summary_FULL.docx")

message("✅ Threshold summary saved: outputs/Ghana_Risk_Thresholds_Summary_FULL.docx")




# -----------------------------
# FIX: clip lakes/lagoons to Ghana
# -----------------------------
ghana_sf <- sf::st_transform(sf::st_as_sf(shp2), 4326)
lakes_sf <- sf::st_transform(sf::st_as_sf(shp_lakes), 4326)

ghana_sf <- sf::st_make_valid(ghana_sf)
lakes_sf <- sf::st_make_valid(lakes_sf)

if (sf::st_crs(lakes_sf) != sf::st_crs(ghana_sf)) {
  lakes_sf <- sf::st_transform(lakes_sf, sf::st_crs(ghana_sf))
}

ghana_union <- sf::st_union(ghana_sf)

lakes_sf_correct <- suppressWarnings(
  sf::st_intersection(lakes_sf, ghana_union)
)

# remove empty/problem geometries
lakes_sf_correct <- lakes_sf_correct[!sf::st_is_empty(lakes_sf_correct), ]

# keep only polygon geometries
lakes_sf_correct <- suppressWarnings(
  sf::st_collection_extract(lakes_sf_correct, "POLYGON")
)

# ==============================================================================
# FULL CORRECTED SCRIPT: SPATIOTEMPORAL COMBINATIONS & GRIDS
# Lakes/lagoons now appear in the legend
# ==============================================================================

library(terra)
library(ggplot2)
library(sf)
library(patchwork)
library(tidyterra)
library(grid)

risk_palette <- c(
  "0" = "#D9D9D9",
  "1" = "#FFD700",
  "2" = "#9400D3",
  "3" = "#00FA9A",
  "4" = "#FF0000",
  "5" = "#FF8C00",
  "6" = "#FFB6C1",
  "7" = "#5D4037",
  "8" = "lightblue"
)

labels_custom <- c(
  "0" = "all low",
  "1" = "A = poor sanitation",
  "2" = "B = high pig density",
  "3" = "AB",
  "4" = "C = high poverty",
  "5" = "AC",
  "6" = "BC",
  "7" = "ABC",
  "8" = "lakes & lagoons"
)

dirs <- c(
  "outputs/HighRes/Combinations_No_Axes_LakesBlue_v5",
  "outputs/HighRes/Combinations_With_Axes_LakesBlue_v5",
  "outputs/HighRes/Grids_NoAxes_LakesBlue_v5",
  "outputs/HighRes/Grids_WithAxes_LakesBlue_v5",
  "outputs/Combinations_No_Axes_LakesBlue_v5",
  "outputs/Combinations_With_Axes_LakesBlue_v5"
)

lapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE))

variable_map <- c(
  "wc1" = "Unimproved San.",
  "wc2" = "Open Defecation",
  "w1"  = "Wealth Q1",
  "w1w2" = "Poverty W1+W2",
  "Mean Pig (DHS)" = "Mean Pig (DHS)",
  "Pig_Density_GLW1" = "Pig_Density_GLW1",
  "Pig_Density_GLW2" = "Pig_Density_GLW2",
  "Pig_Density_GLW3" = "Pig_Density_GLW3",
  "Pig_Density_GLW4" = "Pig_Density_GLW4"
)

# -----------------------------
# Helper function for triple overlap plots
# -----------------------------
plot_triple_overlap_clean <- function(rA, rB, rC, shp, q_val) {
  
  valsA <- values(rA, mat = FALSE)
  valsB <- values(rB, mat = FALSE)
  valsC <- values(rC, mat = FALSE)
  
  valsA <- valsA[!is.na(valsA)]
  valsB <- valsB[!is.na(valsB)]
  valsC <- valsC[!is.na(valsC)]
  
  thA <- quantile(valsA, q_val, na.rm = TRUE)
  thB <- quantile(valsB, q_val, na.rm = TRUE)
  thC <- quantile(valsC, q_val, na.rm = TRUE)
  
  overlap <- mask(
    crop(
      (rA >= thA) * 1 +
        (rB >= thB) * 2 +
        (rC >= thC) * 4,
      shp
    ),
    shp
  )
  
  overlap_f <- as.factor(overlap)
  
  levels(overlap_f) <- data.frame(
    ID = 0:7,
    class = as.character(0:7)
  )
  
  ggplot() +
    geom_spatraster(data = overlap_f) +
    
    geom_sf(
      data = lakes_sf_correct,
      aes(fill = "8"),
      color = NA,
      inherit.aes = FALSE,
      show.legend = TRUE
    ) +
    
    geom_sf(
      data = st_as_sf(shp),
      fill = NA,
      color = "black",
      size = 0.3,
      inherit.aes = FALSE
    ) +
    
    scale_fill_manual(
      values = risk_palette,
      labels = labels_custom,
      name = "Risk\nFactor",
      breaks = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
      limits = c("0", "1", "2", "3", "4", "5", "6", "7", "8"),
      na.value = "white",
      drop = FALSE
    ) +
    
    guides(
      fill = guide_legend(
        ncol = 3,
        byrow = TRUE,
        title.position = "top",
        title.hjust = 0.5,
        override.aes = list(
          color = NA,
          alpha = 1
        )
      )
    )
}

# -----------------------------
# EVERY POSSIBLE COMBINATION
# -----------------------------
for (yr in names(all_rasters)) {
  
  A_opts <- intersect(c("wc1", "wc2"), names(variable_map))
  B_opts <- intersect(
    c("Mean Pig (DHS)", "Pig_Density_GLW1", "Pig_Density_GLW2",
      "Pig_Density_GLW3", "Pig_Density_GLW4"),
    names(variable_map)
  )
  C_opts <- intersect(c("w1", "w1w2"), names(variable_map))
  
  for (q_name in names(percentiles)) {
    
    q_val <- percentiles[[q_name]]
    
    for (a in A_opts) for (b in B_opts) for (c in C_opts) {
      
      rastA <- all_rasters[[yr]][[variable_map[[a]]]]
      rastB <- all_rasters[[yr]][[variable_map[[b]]]]
      rastC <- all_rasters[[yr]][[variable_map[[c]]]]
      
      if (is.null(rastA) || is.null(rastB) || is.null(rastC)) next
      
      p_base <- plot_triple_overlap_clean(
        rastA,
        rastB,
        rastC,
        ghana_vect,
        q_val
      )
      
      p_no_axes <- p_base +
        theme_void() +
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 15),
          legend.key.width = unit(1.1, "cm"),
          legend.key.height = unit(0.7, "cm"),
          plot.margin = margin(12, 12, 12, 12),
          plot.background = element_rect(fill = "white", color = NA)
        )
      
      ggsave(
        paste0(
          "outputs/Combinations_No_Axes_LakesBlue_v5/",
          yr, "_", q_name, "_",
          gsub(" ", "_", a), "_",
          gsub(" ", "_", b), "_",
          gsub(" ", "_", c),
          "_lakesblue_v5.png"
        ),
        p_no_axes,
        width = 8,
        height = 10,
        bg = "white"
      )
      
      ggsave(
        paste0(
          "outputs/HighRes/Combinations_No_Axes_LakesBlue_v5/",
          yr, "_", q_name, "_",
          gsub(" ", "_", a), "_",
          gsub(" ", "_", b), "_",
          gsub(" ", "_", c),
          "_HR_lakesblue_v5.png"
        ),
        p_no_axes,
        width = 12,
        height = 15,
        dpi = 300,
        bg = "white"
      )
      
      p_with_axes <- p_base +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 18, face = "bold"),
          legend.text = element_text(size = 15),
          legend.key.width = unit(1.1, "cm"),
          legend.key.height = unit(0.7, "cm"),
          plot.margin = margin(12, 12, 12, 12),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)
        )
      
      ggsave(
        paste0(
          "outputs/Combinations_With_Axes_LakesBlue_v5/",
          yr, "_", q_name, "_",
          gsub(" ", "_", a), "_",
          gsub(" ", "_", b), "_",
          gsub(" ", "_", c),
          "_Axes_lakesblue_v5.png"
        ),
        p_with_axes,
        width = 8,
        height = 10,
        bg = "white"
      )
      
      ggsave(
        paste0(
          "outputs/HighRes/Combinations_With_Axes_LakesBlue_v5/",
          yr, "_", q_name, "_",
          gsub(" ", "_", a), "_",
          gsub(" ", "_", b), "_",
          gsub(" ", "_", c),
          "_Axes_HR_lakesblue_v5.png"
        ),
        p_with_axes,
        width = 12,
        height = 15,
        dpi = 300,
        bg = "white"
      )
    }
  }
}

# -----------------------------
# Curated 2x2 spatiotemporal grids
# -----------------------------
glw_assignment <- list(
  "2003" = "Pig_Density_GLW1",
  "2008" = "Pig_Density_GLW2",
  "2014" = "Pig_Density_GLW3",
  "2022" = "Pig_Density_GLW4"
)

for (p_lab in names(percentiles)) {
  
  q_val <- percentiles[[p_lab]]
  plot_list_no <- list()
  plot_list_with <- list()
  
  for (yr in names(glw_assignment)) {
    
    pig_var <- glw_assignment[[yr]]
    
    rA <- all_rasters[[yr]][[variable_map[["wc1"]]]]
    rB <- all_rasters[[yr]][[variable_map[[pig_var]]]]
    rC <- all_rasters[[yr]][[variable_map[["w1w2"]]]]
    
    if (is.null(rA) || is.null(rB) || is.null(rC)) next
    
    p_grid <- plot_triple_overlap_clean(
      rA,
      rB,
      rC,
      ghana_vect,
      q_val
    ) +
      labs(title = yr)
    
    plot_list_no[[yr]] <- p_grid +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
        legend.position = "none",
        plot.margin = margin(14, 14, 14, 14),
        plot.background = element_rect(fill = "white", color = NA)
      )
    
    plot_list_with[[yr]] <- p_grid +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
        legend.position = "none",
        plot.margin = margin(14, 14, 14, 14),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      )
  }
  
  final_grid_no <- (
    plot_list_no[["2003"]] + plot_list_no[["2008"]]
  ) / (
    plot_list_no[["2014"]] + plot_list_no[["2022"]]
  ) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 17),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.75, "cm"),
      plot.margin = margin(16, 16, 16, 16),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  final_grid_with <- (
    plot_list_with[["2003"]] + plot_list_with[["2008"]]
  ) / (
    plot_list_with[["2014"]] + plot_list_with[["2022"]]
  ) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 20, face = "bold"),
      legend.text = element_text(size = 17),
      legend.key.width = unit(1.2, "cm"),
      legend.key.height = unit(0.75, "cm"),
      plot.margin = margin(16, 16, 16, 16),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(
    paste0("outputs/Ghana_Grid_NoAxes_", p_lab, "_lakesblue_v5.png"),
    final_grid_no,
    width = 12,
    height = 14,
    bg = "white"
  )
  
  ggsave(
    paste0("outputs/Ghana_Grid_WithAxes_", p_lab, "_lakesblue_v5.png"),
    final_grid_with,
    width = 12,
    height = 14,
    bg = "white"
  )
  
  ggsave(
    paste0(
      "outputs/HighRes/Grids_NoAxes_LakesBlue_v5/Ghana_Grid_NoAxes_",
      p_lab,
      "_HR_lakesblue_v5.png"
    ),
    final_grid_no,
    width = 18,
    height = 21,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    paste0(
      "outputs/HighRes/Grids_WithAxes_LakesBlue_v5/Ghana_Grid_WithAxes_",
      p_lab,
      "_HR_lakesblue_v5.png"
    ),
    final_grid_with,
    width = 18,
    height = 21,
    dpi = 300,
    bg = "white"
  )
}

cat("✅ Done: lakes & lagoons are now included in the legend.\n")





# -----------------------------
# FIX: clip lakes/lagoons to Ghana FULL boundary
# IMPORTANT: use shp2, NOT ghana_vect
# because ghana_vect is Ghana land with lakes removed
# -----------------------------
ghana_full_sf <- sf::st_transform(sf::st_as_sf(shp2), 4326)
lakes_sf <- sf::st_transform(sf::st_as_sf(shp_lakes), 4326)

ghana_full_sf <- sf::st_make_valid(ghana_full_sf)
lakes_sf <- sf::st_make_valid(lakes_sf)

if (sf::st_crs(lakes_sf) != sf::st_crs(ghana_full_sf)) {
  lakes_sf <- sf::st_transform(lakes_sf, sf::st_crs(ghana_full_sf))
}

ghana_union <- sf::st_union(ghana_full_sf)

lakes_sf_correct <- suppressWarnings(
  sf::st_intersection(lakes_sf, ghana_union)
)

lakes_sf_correct <- lakes_sf_correct[!sf::st_is_empty(lakes_sf_correct), ]

lakes_sf_correct <- suppressWarnings(
  sf::st_collection_extract(lakes_sf_correct, "POLYGON")
)



# ============================================================================== #
# FINAL 2 x 2 GRID PER YEAR
# top-left  = sanitation (wc1)
# top-right = pig (GLW)
# bottom-left = poverty (w1w2)
# bottom-right = combined
# - Uses the SAME binary logic as Section 7
# - Works for both 66th and 75th percentiles
# - No titles on maps
# - One legend at the bottom
# - Lakes/lagoons clipped to Ghana full boundary
# - New output names so old files are not replaced
# - Slightly increased map spacing
# - Slightly increased legend text
# ============================================================================== #

# -----------------------------
# Output directory
# -----------------------------
dir.create("outputs/HighRes/Grids_2x2_ByYear_Final_LakesBlue_v7_BaseR", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# GLW assignment by year
# -----------------------------
glw_assignment <- list(
  "2003" = "Pig_Density_GLW1",
  "2008" = "Pig_Density_GLW2",
  "2014" = "Pig_Density_GLW3",
  "2022" = "Pig_Density_GLW4"
)

# -----------------------------
# Colours
# -----------------------------
col_low   <- "#D9D9D9"
col_wc1   <- "#FFD700"
col_glw   <- "#9400D3"
col_w1w2  <- "#FF0000"
col_lakes <- "lightblue"

col_combined <- c(
  "0" = "#D9D9D9", # all low
  "1" = "#FFD700", # A
  "2" = "#9400D3", # B
  "3" = "#00FA9A", # AB
  "4" = "#FF0000", # C
  "5" = "#FF8C00", # AC
  "6" = "#FFB6C1", # BC
  "7" = "#5D4037"  # ABC
)

# -----------------------------
# Helper: EXACT same binary logic as Section 7
# -----------------------------
make_binary_section7 <- function(r, q_val) {
  vals <- values(r, mat = FALSE)
  vals <- vals[!is.na(vals) & vals >= 0]
  if(length(vals) == 0) return(NULL)
  
  q_real <- quantile(vals, probs = q_val, na.rm = TRUE)
  bin_rast <- r >= q_real
  return(bin_rast)
}

# -----------------------------
# Helper: binary panel
# -----------------------------
plot_binary_panel <- function(bin_rast, high_col) {
  par(mar = c(0.75, 0.75, 0.75, 0.75))
  plot(
    bin_rast,
    col = c(col_low, high_col),
    axes = FALSE,
    box = FALSE,
    legend = FALSE
  )
  plot(lakes_sf_correct, add = TRUE, col = col_lakes, border = NA)
  plot(ghana_vect, add = TRUE, border = "black", lwd = 0.5)
}

# -----------------------------
# Helper: combined panel
# -----------------------------
plot_combined_panel <- function(bin_wc1, bin_glw, bin_w1w2) {
  par(mar = c(0.75, 0.75, 0.75, 0.75))
  
  combined <- (bin_wc1 * 1) + (bin_glw * 2) + (bin_w1w2 * 4)
  
  plot(
    combined,
    col = unname(col_combined[as.character(0:7)]),
    type = "classes",
    axes = FALSE,
    box = FALSE,
    legend = FALSE
  )
  plot(lakes_sf_correct, add = TRUE, col = col_lakes, border = NA)
  plot(ghana_vect, add = TRUE, border = "black", lwd = 0.5)
}

# -----------------------------
# Helper: one bottom legend
# -----------------------------
draw_single_bottom_legend <- function() {
  par(mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  
  legend(
    "center",
    title = "Risk factor",
    legend = c(
      "all low",
      "A = poor sanitation",
      "B = high pig density",
      "AB",
      "C = high poverty",
      "AC",
      "BC",
      "ABC",
      "lakes & lagoons"
    ),
    fill = c(
      "#D9D9D9",
      "#FFD700",
      "#9400D3",
      "#00FA9A",
      "#FF0000",
      "#FF8C00",
      "#FFB6C1",
      "#5D4037",
      "lightblue"
    ),
    border = NA,
    bty = "n",
    horiz = FALSE,
    ncol = 3,
    cex = 1.35,
    pt.cex = 2.4,
    title.cex = 1.45
  )
}

# -----------------------------
# Build 2 x 2 grids for both thresholds
# -----------------------------
for (pname in names(percentiles)) {
  
  q_val <- percentiles[[pname]]
  
  for (yr in names(glw_assignment)) {
    
    pig_label <- glw_assignment[[yr]]
    
    r_wc1  <- all_rasters[[yr]][[ variable_map[["wc1"]] ]]
    r_glw  <- all_rasters[[yr]][[ variable_map[[pig_label]] ]]
    r_w1w2 <- all_rasters[[yr]][[ variable_map[["w1w2"]] ]]
    
    if (is.null(r_wc1) || is.null(r_glw) || is.null(r_w1w2)) next
    
    wc1_bin  <- make_binary_section7(r_wc1, q_val)
    glw_bin  <- make_binary_section7(r_glw, q_val)
    w1w2_bin <- make_binary_section7(r_w1w2, q_val)
    
    if (is.null(wc1_bin) || is.null(glw_bin) || is.null(w1w2_bin)) next
    
    out_file <- paste0(
      "outputs/HighRes/Grids_2x2_ByYear_Final_LakesBlue_v7_BaseR/",
      yr, "_", pname, "_2x2_lakesblue_v7_baser.png"
    )
    
    png(out_file, width = 4200, height = 4200, res = 300)
    
    layout(
      matrix(c(1, 2,
               3, 4,
               5, 5), nrow = 3, byrow = TRUE),
      heights = c(1, 1, 0.28)
    )
    
    par(oma = c(0.6, 0.6, 0.6, 0.6))
    
    # top-left
    plot_binary_panel(wc1_bin, col_wc1)
    
    # top-right
    plot_binary_panel(glw_bin, col_glw)
    
    # bottom-left
    plot_binary_panel(w1w2_bin, col_w1w2)
    
    # bottom-right
    plot_combined_panel(wc1_bin, glw_bin, w1w2_bin)
    
    # legend row
    draw_single_bottom_legend()
    
    dev.off()
  }
}

cat("✅ Done: 2 x 2 base-R grids created for 66th and 75th percentiles with corrected lakes, larger legend, and slightly more spacing.\n")



# ============================================================================== #
# FINAL 2 x 2 GRID PER YEAR (ggplot2 + ggsave)
# top-left     = sanitation (wc1)
# top-right    = pig (GLW)
# bottom-left  = poverty (w1w2)
# bottom-right = combined
#
# FEATURES
# - Uses the SAME binary logic as Section 7
# - Works for both 66th and 75th percentiles
# - No titles on maps
# - One large readable legend at the bottom
# - High-quality export with ggsave()
# - Outside Ghana is PURE WHITE
# - Legend title = "Risk Factor"
# - Standard spacing between panels
# - Lakes/lagoons clipped to Ghana full boundary and added as separate sf layer
# - New output names so old files are not replaced
# - Slightly increased map spacing
# - Slightly increased legend text
# ============================================================================== #

# -----------------------------
# Packages
# -----------------------------
library(terra)
library(sf)
library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

# -----------------------------
# Output directory
# -----------------------------
dir.create("outputs/HighRes/Grids_2x2_ByYear_Final_LakesBlue_v7_GG", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# GLW assignment by year
# -----------------------------
glw_assignment <- list(
  "2003" = "Pig_Density_GLW1",
  "2008" = "Pig_Density_GLW2",
  "2014" = "Pig_Density_GLW3",
  "2022" = "Pig_Density_GLW4"
)

# -----------------------------
# Colours
# -----------------------------
col_low   <- "#D9D9D9"
col_wc1   <- "#FFD700"
col_glw   <- "#9400D3"
col_w1w2  <- "#FF0000"
col_lakes <- "lightblue"

col_combined <- c(
  "0" = "#D9D9D9", # all low
  "1" = "#FFD700", # A
  "2" = "#9400D3", # B
  "3" = "#00FA9A", # AB
  "4" = "#FF0000", # C
  "5" = "#FF8C00", # AC
  "6" = "#FFB6C1", # BC
  "7" = "#5D4037"  # ABC
)

# -----------------------------
# Convert vectors to sf once
# -----------------------------
ghana_sf <- st_as_sf(ghana_vect)
lakes_sf <- lakes_sf_correct
ghana_vect_terra <- terra::vect(ghana_sf)

# -----------------------------
# Helper: EXACT same binary logic as Section 7
# -----------------------------
make_binary_section7 <- function(r, q_val) {
  vals <- values(r, mat = FALSE)
  vals <- vals[!is.na(vals) & vals >= 0]
  
  if (length(vals) == 0) return(NULL)
  
  q_real <- quantile(vals, probs = q_val, na.rm = TRUE)
  ifel(r >= q_real, 1, 0)
}

# -----------------------------
# Helper: crop + mask to Ghana
# This keeps outside Ghana pure white
# -----------------------------
crop_mask_to_ghana <- function(r) {
  if (is.null(r)) return(NULL)
  r2 <- terra::crop(r, ghana_vect_terra, snap = "out")
  r2 <- terra::mask(r2, ghana_vect_terra)
  r2
}

# -----------------------------
# Helper: raster to data.frame
# -----------------------------
rast_to_df <- function(r, value_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- value_name
  df
}

# -----------------------------
# Helper: theme for maps
# Standard spacing via plot.margin
# -----------------------------
theme_map_clean <- function() {
  theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(16, 24, 16, 24),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# -----------------------------
# Helper: binary panel
# -----------------------------
plot_binary_panel_gg <- function(bin_rast, high_col) {
  df <- rast_to_df(bin_rast, "class")
  
  df$class <- factor(
    df$class,
    levels = c(0, 1),
    labels = c("Low", "High")
  )
  
  ggplot() +
    geom_raster(
      data = df,
      aes(x = x, y = y, fill = class),
      na.rm = FALSE
    ) +
    geom_sf(
      data = lakes_sf,
      fill = col_lakes,
      colour = NA,
      inherit.aes = FALSE
    ) +
    geom_sf(
      data = ghana_sf,
      fill = NA,
      colour = "black",
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c("Low" = col_low, "High" = high_col),
      drop = FALSE,
      na.value = "white"
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean()
}

# -----------------------------
# Helper: combined panel
# -----------------------------
plot_combined_panel_gg <- function(bin_wc1, bin_glw, bin_w1w2) {
  combined <- (bin_wc1 * 1) + (bin_glw * 2) + (bin_w1w2 * 4)
  combined <- crop_mask_to_ghana(combined)
  
  df <- rast_to_df(combined, "group")
  df$group <- factor(df$group, levels = 0:7)
  
  ggplot() +
    geom_raster(
      data = df,
      aes(x = x, y = y, fill = group),
      na.rm = FALSE
    ) +
    geom_sf(
      data = lakes_sf,
      fill = col_lakes,
      colour = NA,
      inherit.aes = FALSE
    ) +
    geom_sf(
      data = ghana_sf,
      fill = NA,
      colour = "black",
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = col_combined,
      drop = FALSE,
      na.value = "white"
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean()
}

# -----------------------------
# Helper: one big bottom legend
# -----------------------------
make_bottom_legend <- function() {
  legend_labels <- c(
    "all low",
    "A = poor sanitation",
    "B = high pig density",
    "AB",
    "C = high poverty",
    "AC",
    "BC",
    "ABC",
    "lakes & lagoons"
  )
  
  legend_cols <- c(
    "all low" = "#D9D9D9",
    "A = poor sanitation" = "#FFD700",
    "B = high pig density" = "#9400D3",
    "AB" = "#00FA9A",
    "C = high poverty" = "#FF0000",
    "AC" = "#FF8C00",
    "BC" = "#FFB6C1",
    "ABC" = "#5D4037",
    "lakes & lagoons" = "lightblue"
  )
  
  legend_df <- data.frame(
    x = seq_along(legend_labels),
    y = 1,
    lab = factor(legend_labels, levels = legend_labels)
  )
  
  p_leg <- ggplot(legend_df, aes(x = x, y = y, fill = lab)) +
    geom_tile() +
    scale_fill_manual(
      values = legend_cols,
      breaks = legend_labels,
      drop = FALSE,
      name = "Risk Factor"
    ) +
    guides(
      fill = guide_legend(
        ncol = 3,
        byrow = TRUE,
        override.aes = list(colour = NA)
      )
    ) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 19),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.82, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      legend.spacing.y = unit(0.28, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
  
  cowplot::get_legend(p_leg)
}

# -----------------------------
# Main loop: build 2 x 2 grids for both thresholds
# -----------------------------
legend_grob <- make_bottom_legend()

for (pname in names(percentiles)) {
  
  q_val <- percentiles[[pname]]
  
  for (yr in names(glw_assignment)) {
    
    pig_label <- glw_assignment[[yr]]
    
    r_wc1  <- all_rasters[[yr]][[ variable_map[["wc1"]] ]]
    r_glw  <- all_rasters[[yr]][[ variable_map[[pig_label]] ]]
    r_w1w2 <- all_rasters[[yr]][[ variable_map[["w1w2"]] ]]
    
    if (is.null(r_wc1) || is.null(r_glw) || is.null(r_w1w2)) next
    
    wc1_bin  <- make_binary_section7(r_wc1, q_val)
    glw_bin  <- make_binary_section7(r_glw, q_val)
    w1w2_bin <- make_binary_section7(r_w1w2, q_val)
    
    if (is.null(wc1_bin) || is.null(glw_bin) || is.null(w1w2_bin)) next
    
    # keep outside Ghana white
    wc1_bin  <- crop_mask_to_ghana(wc1_bin)
    glw_bin  <- crop_mask_to_ghana(glw_bin)
    w1w2_bin <- crop_mask_to_ghana(w1w2_bin)
    
    # panels
    p1 <- plot_binary_panel_gg(wc1_bin, col_wc1)
    p2 <- plot_binary_panel_gg(glw_bin, col_glw)
    p3 <- plot_binary_panel_gg(w1w2_bin, col_w1w2)
    p4 <- plot_combined_panel_gg(wc1_bin, glw_bin, w1w2_bin)
    
    # standard patchwork spacing
    final_plot <-
      (
        (p1 + p2) /
          (p3 + p4) /
          wrap_elements(legend_grob)
      ) +
      plot_layout(heights = c(1, 1, 0.26)) &
      theme(
        plot.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(18, 18, 18, 18)
      )
    
    out_file <- paste0(
      "outputs/HighRes/Grids_2x2_ByYear_Final_LakesBlue_v7_GG/",
      yr, "_", pname, "_2x2_lakesblue_v7_gg.png"
    )
    
    ggsave(
      filename = out_file,
      plot = final_plot,
      width = 16,
      height = 16,
      units = "in",
      dpi = 400,
      bg = "white",
      limitsize = FALSE
    )
  }
}

cat("✅ Done: 2 x 2 ggplot grids created for 66th and 75th percentiles using corrected lakes, larger legends, and slightly increased map spacing.\n")




# -----------------------------
# FIX: clip lakes/lagoons to Ghana FULL boundary
# IMPORTANT: use shp2, NOT ghana_vect
# because ghana_vect is Ghana land with lakes removed
# -----------------------------
ghana_full_sf <- sf::st_transform(sf::st_as_sf(shp2), 4326)
lakes_sf_raw <- sf::st_transform(sf::st_as_sf(shp_lakes), 4326)

ghana_full_sf <- sf::st_make_valid(ghana_full_sf)
lakes_sf_raw <- sf::st_make_valid(lakes_sf_raw)

if (sf::st_crs(lakes_sf_raw) != sf::st_crs(ghana_full_sf)) {
  lakes_sf_raw <- sf::st_transform(lakes_sf_raw, sf::st_crs(ghana_full_sf))
}

ghana_union <- sf::st_union(ghana_full_sf)

lakes_sf_correct <- suppressWarnings(
  sf::st_intersection(lakes_sf_raw, ghana_union)
)

lakes_sf_correct <- lakes_sf_correct[!sf::st_is_empty(lakes_sf_correct), ]

lakes_sf_correct <- suppressWarnings(
  sf::st_collection_extract(lakes_sf_correct, "POLYGON")
)



# ============================================================================== #
# SENSITIVITY ANALYSIS MAPS - 2022 ONLY
# Builds maps for:
# 1. Full maps of all scenario tests
# 2. ABC only
# 3a. AB, BC, ABC separate (others grey)
# 3b. AB + BC + ABC combined as "Risk of Interest" (others grey)
# Uses only the 2022 rasters already stored in all_rasters
# Thresholds used: 66th and 75th
# Lakes/lagoons clipped to Ghana full boundary and added as separate sf layer
# Legend sizes and map spacing kept consistent
# ============================================================================== #

# ============================================================================== #
# 1. SETUP
# ============================================================================== #
library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(cowplot)
library(grid)

# --- DIRECT MAPPING: SYNCHRONIZING WITH ALL_RASTERS ---
wc1          <- all_rasters[["2022"]][["Unimproved San."]]
wc2          <- all_rasters[["2022"]][["Open Defecation"]]
w1w2         <- all_rasters[["2022"]][["Poverty W1+W2"]]
w1           <- all_rasters[["2022"]][["Wealth Q1"]]
pig_glw4     <- all_rasters[["2022"]][["Pig_Density_GLW4"]]
pig_dhs_mean <- all_rasters[["2022"]][["Mean Pig (DHS)"]]

ghana_vect <- vect(shp_subnational)
ghana_sf   <- st_as_sf(shp_subnational)
lakes_sf   <- lakes_sf_correct

area_rast  <- cellSize(wc1, unit = "km")
area_ghana <- mask(area_rast, ghana_vect)
total_gh_area <- global(area_ghana, sum, na.rm = TRUE)[1,1]

# Threshold labels and values
percentiles <- c("66th" = 0.666, "75th" = 0.75)

# ============================================================================== #
# 2. OUTPUT FOLDERS
# ============================================================================== #
dir.create("outputs/Sensitivity_2022_LakesBlue_v2", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/01_Full_Maps", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/02_ABC_Only", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/03_AB_BC_ABC_Separate", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/04_RiskOfInterest_Combined", recursive = TRUE, showWarnings = FALSE)

# ============================================================================== #
# 3. COLOURS
# Bitmask convention:
# A = sanitation = 1
# B = pig        = 2
# C = poverty    = 4
# ============================================================================== #
col_low   <- "#D9D9D9"
col_A     <- "#FFD700"
col_B     <- "#9400D3"
col_AB    <- "#00FA9A"
col_C     <- "#FF0000"
col_AC    <- "#FF8C00"
col_BC    <- "#FFB6C1"
col_ABC   <- "#5D4037"
col_lakes <- "lightblue"

# combined risk-of-interest colour for AB + BC + ABC
col_roi   <- "#74A586"

col_full <- c(
  "0" = col_low,
  "1" = col_A,
  "2" = col_B,
  "3" = col_AB,
  "4" = col_C,
  "5" = col_AC,
  "6" = col_BC,
  "7" = col_ABC
)

# ============================================================================== #
# 4. HELPERS
# ============================================================================== #

# Same threshold logic as your earlier code: mask to Ghana, then threshold
create_bin <- function(r, q) {
  if (!inherits(r, "SpatRaster")) r <- rast(r)
  r_masked <- mask(r, ghana_vect)
  vals <- values(r_masked, mat = FALSE)
  vals <- vals[!is.na(vals) & vals >= 0]
  if (length(vals) == 0) return(NULL)
  thresh <- quantile(vals, q, na.rm = TRUE)
  r_masked >= thresh
}

crop_mask_to_ghana <- function(r) {
  if (is.null(r)) return(NULL)
  r2 <- crop(r, ghana_vect, snap = "out")
  r2 <- mask(r2, ghana_vect)
  r2
}

rast_to_df <- function(r, value_name = "value") {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- value_name
  df
}

theme_map_clean <- function() {
  theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(12, 12, 12, 12),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA)
    )
}

# ------------------------------------------------------------------------------ #
# Plot: full 0-7 combination map
# ------------------------------------------------------------------------------ #
plot_full_combo <- function(r) {
  df <- rast_to_df(r, "group")
  df$group <- factor(df$group, levels = 0:7)
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = group), na.rm = FALSE) +
    geom_sf(data = lakes_sf, fill = col_lakes, colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(values = col_full, drop = FALSE, na.value = "white") +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean()
}

# ------------------------------------------------------------------------------ #
# Plot: two-class map (other vs selected risk)
# ------------------------------------------------------------------------------ #
plot_two_class <- function(r, high_col) {
  df <- rast_to_df(r, "group")
  df$group <- factor(df$group, levels = c(0, 1), labels = c("Other", "Risk"))
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = group), na.rm = FALSE) +
    geom_sf(data = lakes_sf, fill = col_lakes, colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(
      values = c("Other" = col_low, "Risk" = high_col),
      drop = FALSE,
      na.value = "white"
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean()
}

# ------------------------------------------------------------------------------ #
# Plot: four-class map (other, AB, BC, ABC)
# ------------------------------------------------------------------------------ #
plot_ab_bc_abc <- function(r) {
  df <- rast_to_df(r, "group")
  df$group <- factor(df$group, levels = c(0, 1, 2, 3),
                     labels = c("Other", "AB", "BC", "ABC"))
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = group), na.rm = FALSE) +
    geom_sf(data = lakes_sf, fill = col_lakes, colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(
      values = c(
        "Other" = col_low,
        "AB"    = col_AB,
        "BC"    = col_BC,
        "ABC"   = col_ABC
      ),
      drop = FALSE,
      na.value = "white"
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean()
}

# ------------------------------------------------------------------------------ #
# Legends
# ------------------------------------------------------------------------------ #
make_legend_full <- function() {
  labs <- c("all low", "A = poor sanitation", "B = high pig density",
            "AB", "C = high poverty", "AC", "BC", "ABC", "lakes & lagoons")
  cols <- c(
    "all low" = col_low,
    "A = poor sanitation" = col_A,
    "B = high pig density" = col_B,
    "AB" = col_AB,
    "C = high poverty" = col_C,
    "AC" = col_AC,
    "BC" = col_BC,
    "ABC" = col_ABC,
    "lakes & lagoons" = col_lakes
  )
  
  dd <- data.frame(x = seq_along(labs), y = 1, lab = factor(labs, levels = labs))
  
  p <- ggplot(dd, aes(x = x, y = y, fill = lab)) +
    geom_tile() +
    scale_fill_manual(values = cols, breaks = labs, drop = FALSE, name = "Risk Factor") +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(colour = NA))) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 19),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.82, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      legend.spacing.y = unit(0.28, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  cowplot::get_legend(p)
}

make_legend_two_class <- function(high_lab, high_col) {
  labs <- c("all low / other risks", high_lab, "lakes & lagoons")
  
  cols <- c(col_low, high_col, col_lakes)
  names(cols) <- labs
  
  dd <- data.frame(x = seq_along(labs), y = 1, lab = factor(labs, levels = labs))
  
  p <- ggplot(dd, aes(x = x, y = y, fill = lab)) +
    geom_tile() +
    scale_fill_manual(values = cols, breaks = labs, drop = FALSE, name = "Risk Factor") +
    guides(fill = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(colour = NA))) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 18),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.82, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      legend.spacing.y = unit(0.28, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  cowplot::get_legend(p)
}

make_legend_ab_bc_abc <- function() {
  labs <- c("all low / other risks", "AB", "BC", "ABC", "lakes & lagoons")
  cols <- c(
    "all low / other risks" = col_low,
    "AB" = col_AB,
    "BC" = col_BC,
    "ABC" = col_ABC,
    "lakes & lagoons" = col_lakes
  )
  
  dd <- data.frame(x = seq_along(labs), y = 1, lab = factor(labs, levels = labs))
  
  p <- ggplot(dd, aes(x = x, y = y, fill = lab)) +
    geom_tile() +
    scale_fill_manual(values = cols, breaks = labs, drop = FALSE, name = "Risk Factor") +
    guides(fill = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(colour = NA))) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 22, face = "bold"),
      legend.text = element_text(size = 19),
      legend.key.width = unit(1.3, "cm"),
      legend.key.height = unit(0.82, "cm"),
      legend.spacing.x = unit(0.5, "cm"),
      legend.spacing.y = unit(0.28, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(0, 0, 0, 0)
    )
  
  cowplot::get_legend(p)
}

legend_full  <- make_legend_full()

legend_abc   <- make_legend_two_class(
  high_lab = "ABC = triple risk",
  high_col = col_ABC
)

legend_roi   <- make_legend_two_class(
  high_lab = "ABC + AB + BC (High pig density with poor sanitation and/or high poverty)",
  high_col = col_roi
)

legend_three <- make_legend_ab_bc_abc()

# ============================================================================== #
# 5. SCENARIOS / TESTS
# "Model" here = scenario / risk map
# ============================================================================== #

build_combo <- function(S_rast, B_rast, C_rast, qS, qB, qC) {
  S_bin <- create_bin(S_rast, qS)
  B_bin <- create_bin(B_rast, qB)
  C_bin <- create_bin(C_rast, qC)
  
  if (is.null(S_bin) || is.null(B_bin) || is.null(C_bin)) return(NULL)
  
  combo <- S_bin * 1 + B_bin * 2 + C_bin * 4
  crop_mask_to_ghana(combo)
}

# We build the requested 5 tests:
# base   = wc1 + w1w2 + pig_glw4 at 66th
# test_A = wc1 + w1w2 + pig_glw4 at 75th
# test_B = wc1 + w1w2 + pig_dhs_mean at 66th
# test_C = wc2 + w1w2 + pig_glw4 at 66th
# test_D = wc1 + w1 + pig_glw4 at 66th

scenario_defs <- list(
  base   = list(S = wc1, B = pig_glw4,     C = w1w2, qS = 0.666, qB = 0.666, qC = 0.666),
  test_A = list(S = wc1, B = pig_glw4,     C = w1w2, qS = 0.75,  qB = 0.75,  qC = 0.75),
  test_B = list(S = wc1, B = pig_dhs_mean, C = w1w2, qS = 0.666, qB = 0.666, qC = 0.666),
  test_C = list(S = wc2, B = pig_glw4,     C = w1w2, qS = 0.666, qB = 0.666, qC = 0.666),
  test_D = list(S = wc1, B = pig_glw4,     C = w1,   qS = 0.666, qB = 0.666, qC = 0.666)
)

scenario_combos <- lapply(scenario_defs, function(z) {
  build_combo(
    S_rast = z$S,
    B_rast = z$B,
    C_rast = z$C,
    qS = z$qS,
    qB = z$qB,
    qC = z$qC
  )
})

# ============================================================================== #
# 6. EXPORT ALL REQUESTED MAP TYPES
# ============================================================================== #

for (sc_name in names(scenario_combos)) {
  
  combo <- scenario_combos[[sc_name]]
  if (is.null(combo)) next
  
  # label for folder/file
  # use 75th only for test_A, otherwise 66th
  p_lab <- if (sc_name == "test_A") "75th" else "66th"
  
  # --------------------------------------------------------------------------- #
  # 1. FULL MAP OF THE SCENARIO
  # --------------------------------------------------------------------------- #
  p_full <- plot_full_combo(combo)
  final_full <- (p_full / wrap_elements(legend_full)) +
    plot_layout(heights = c(1, 0.20)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(
    filename = paste0("outputs/Sensitivity_2022_LakesBlue_v2/01_Full_Maps/2022_", sc_name, "_", p_lab, "_Full_lakesblue_v2.png"),
    plot = final_full,
    width = 12,
    height = 13,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
  
  # --------------------------------------------------------------------------- #
  # 2. ABC ONLY
  # --------------------------------------------------------------------------- #
  abc_only <- ifel(combo == 7, 1, 0)
  abc_only <- crop_mask_to_ghana(abc_only)
  
  p_abc <- plot_two_class(abc_only, col_ABC)
  final_abc <- (p_abc / wrap_elements(legend_abc)) +
    plot_layout(heights = c(1, 0.18)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(
    filename = paste0("outputs/Sensitivity_2022_LakesBlue_v2/02_ABC_Only/2022_", sc_name, "_", p_lab, "_ABC_Only_lakesblue_v2.png"),
    plot = final_abc,
    width = 12,
    height = 13,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
  
  # --------------------------------------------------------------------------- #
  # 3a. AB, BC, ABC SEPARATE
  # recode:
  # 0 = other
  # 1 = AB
  # 2 = BC
  # 3 = ABC
  # --------------------------------------------------------------------------- #
  ab_bc_abc <- ifel(combo == 3, 1,
                    ifel(combo == 6, 2,
                         ifel(combo == 7, 3, 0)))
  ab_bc_abc <- crop_mask_to_ghana(ab_bc_abc)
  
  p_three <- plot_ab_bc_abc(ab_bc_abc)
  final_three <- (p_three / wrap_elements(legend_three)) +
    plot_layout(heights = c(1, 0.20)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(
    filename = paste0("outputs/Sensitivity_2022_LakesBlue_v2/03_AB_BC_ABC_Separate/2022_", sc_name, "_", p_lab, "_AB_BC_ABC_lakesblue_v2.png"),
    plot = final_three,
    width = 12,
    height = 13,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
  
  # --------------------------------------------------------------------------- #
  # 3b. COMBINED RISK OF INTEREST = AB + BC + ABC
  # --------------------------------------------------------------------------- #
  roi <- ifel((combo == 3) | (combo == 6) | (combo == 7), 1, 0)
  roi <- crop_mask_to_ghana(roi)
  
  p_roi <- plot_two_class(roi, col_roi)
  final_roi <- (p_roi / wrap_elements(legend_roi)) +
    plot_layout(heights = c(1, 0.18)) &
    theme(plot.background = element_rect(fill = "white", colour = NA))
  
  ggsave(
    filename = paste0("outputs/Sensitivity_2022_LakesBlue_v2/04_RiskOfInterest_Combined/2022_", sc_name, "_", p_lab, "_RiskOfInterest_lakesblue_v2.png"),
    plot = final_roi,
    width = 12,
    height = 13,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
}

cat("✅ Done: 2022 sensitivity analysis maps created for base, test_A, test_B, test_C, and test_D.\n")
cat("✅ Outputs saved under outputs/Sensitivity_2022_LakesBlue_v2/\n")




# ============================================================================== #
# 7. EXTRA COMPARISON PANELS
# 1x2: base vs test_A
# 1x2: base vs test_B
# 1x3: base vs test_C vs test_D
# For: Full maps, ABC only, AB/BC/ABC separate, Risk of Interest
# No titles, no subtitles, one legend below
# Lakes/lagoons clipped to Ghana full boundary and added as separate sf layer
# ============================================================================== #

# -----------------------------
# Output folders
# -----------------------------
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/Full_Maps", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/ABC_Only", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/AB_BC_ABC_Separate", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/RiskOfInterest", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper: blank-theme version of existing plotters
# -----------------------------
plot_full_combo_notitle <- function(r) {
  df <- rast_to_df(r, "group")
  df$group <- factor(df$group, levels = 0:7)
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = group), na.rm = FALSE) +
    geom_sf(data = lakes_sf, fill = col_lakes, colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(values = col_full, drop = FALSE, na.value = "white") +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean() +
    theme(plot.margin = margin(20, 20, 20, 20))
}

plot_two_class_notitle <- function(r, high_col) {
  df <- rast_to_df(r, "group")
  df$group <- factor(df$group, levels = c(0, 1), labels = c("Other", "Risk"))
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = group), na.rm = FALSE) +
    geom_sf(data = lakes_sf, fill = col_lakes, colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(
      values = c("Other" = col_low, "Risk" = high_col),
      drop = FALSE,
      na.value = "white"
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean() +
    theme(plot.margin = margin(20, 20, 20, 20))
}

plot_ab_bc_abc_notitle <- function(r) {
  df <- rast_to_df(r, "group")
  df$group <- factor(df$group, levels = c(0, 1, 2, 3),
                     labels = c("Other", "AB", "BC", "ABC"))
  
  ggplot() +
    geom_raster(data = df, aes(x = x, y = y, fill = group), na.rm = FALSE) +
    geom_sf(data = lakes_sf, fill = col_lakes, colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(
      values = c(
        "Other" = col_low,
        "AB"    = col_AB,
        "BC"    = col_BC,
        "ABC"   = col_ABC
      ),
      drop = FALSE,
      na.value = "white"
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_map_clean() +
    theme(plot.margin = margin(20, 20, 20, 20))
}

# -----------------------------
# Helper: derive layers from combo
# -----------------------------
make_abc_only <- function(combo) {
  crop_mask_to_ghana(ifel(combo == 7, 1, 0))
}

make_roi <- function(combo) {
  crop_mask_to_ghana(ifel((combo == 3) | (combo == 6) | (combo == 7), 1, 0))
}

make_ab_bc_abc <- function(combo) {
  out <- ifel(combo == 3, 1,
              ifel(combo == 6, 2,
                   ifel(combo == 7, 3, 0)))
  crop_mask_to_ghana(out)
}

# -----------------------------
# Helper: save comparison row with one shared legend
# Extra spacing added with plot_spacer() between maps
# -----------------------------
save_comparison_panel <- function(plot_list, legend_grob, out_file, n_maps) {
  
  if (n_maps == 2) {
    map_row <- plot_list[[1]] + plot_spacer() + plot_list[[2]] +
      plot_layout(ncol = 3, widths = c(1, 0.05, 1))
    out_width <- 15
  } else if (n_maps == 3) {
    map_row <- plot_list[[1]] + plot_spacer() + plot_list[[2]] + plot_spacer() + plot_list[[3]] +
      plot_layout(ncol = 5, widths = c(1, 0.05, 1, 0.05, 1))
    out_width <- 21
  } else {
    map_row <- wrap_plots(plotlist = plot_list, ncol = n_maps)
    out_width <- 14
  }
  
  combined <- map_row / wrap_elements(legend_grob) +
    plot_layout(heights = c(1, 0.18)) &
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(12, 12, 12, 12)
    )
  
  ggsave(
    filename = out_file,
    plot = combined,
    width = out_width,
    height = 8.8,
    units = "in",
    dpi = 400,
    bg = "white",
    limitsize = FALSE
  )
}

# -----------------------------
# Build comparison panels
# -----------------------------
scenario_order_12_A  <- c("base", "test_A")
scenario_order_12_B  <- c("base", "test_B")
scenario_order_13_CD <- c("base", "test_C", "test_D")

comparison_sets <- list(
  "base_vs_testA"       = scenario_order_12_A,
  "base_vs_testB"       = scenario_order_12_B,
  "base_vs_testC_testD" = scenario_order_13_CD
)

for (comp_name in names(comparison_sets)) {
  
  scs <- comparison_sets[[comp_name]]
  n_maps_use <- length(scs)
  
  # --------------------------------------------------
  # FULL MAPS
  # --------------------------------------------------
  p_full_list <- lapply(scs, function(sc) {
    combo <- scenario_combos[[sc]]
    plot_full_combo_notitle(combo)
  })
  
  save_comparison_panel(
    plot_list = p_full_list,
    legend_grob = legend_full,
    out_file = paste0("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/Full_Maps/", comp_name, "_Full_lakesblue_v2.png"),
    n_maps = n_maps_use
  )
  
  # --------------------------------------------------
  # ABC ONLY
  # --------------------------------------------------
  p_abc_list <- lapply(scs, function(sc) {
    combo <- scenario_combos[[sc]]
    abc_only <- make_abc_only(combo)
    plot_two_class_notitle(abc_only, col_ABC)
  })
  
  save_comparison_panel(
    plot_list = p_abc_list,
    legend_grob = legend_abc,
    out_file = paste0("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/ABC_Only/", comp_name, "_ABC_Only_lakesblue_v2.png"),
    n_maps = n_maps_use
  )
  
  # --------------------------------------------------
  # AB / BC / ABC SEPARATE
  # --------------------------------------------------
  p_three_list <- lapply(scs, function(sc) {
    combo <- scenario_combos[[sc]]
    ab_bc_abc <- make_ab_bc_abc(combo)
    plot_ab_bc_abc_notitle(ab_bc_abc)
  })
  
  save_comparison_panel(
    plot_list = p_three_list,
    legend_grob = legend_three,
    out_file = paste0("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/AB_BC_ABC_Separate/", comp_name, "_AB_BC_ABC_lakesblue_v2.png"),
    n_maps = n_maps_use
  )
  
  # --------------------------------------------------
  # RISK OF INTEREST
  # --------------------------------------------------
  p_roi_list <- lapply(scs, function(sc) {
    combo <- scenario_combos[[sc]]
    roi <- make_roi(combo)
    plot_two_class_notitle(roi, col_roi)
  })
  
  save_comparison_panel(
    plot_list = p_roi_list,
    legend_grob = legend_roi,
    out_file = paste0("outputs/Sensitivity_2022_LakesBlue_v2/05_Comparisons/RiskOfInterest/", comp_name, "_RiskOfInterest_lakesblue_v2.png"),
    n_maps = n_maps_use
  )
}

cat("✅ Extra comparison panels created for Full maps, ABC only, AB/BC/ABC separate, and Risk of Interest.\n")




# ============================================================================== #
# 8. CONSENSUS MAPS
# - ABC consensus across the 5 scenario tests
# - Risk of Interest consensus across the 5 scenario tests
# - With axes and without axes
# - 0 = grey
# - lakes & lagoons included in legend
# - legend on the side
# - Lakes/lagoons clipped to Ghana full boundary and added as separate sf layer
# ============================================================================== #

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(colorspace)

# -----------------------------
# Helper layers from combo
# -----------------------------
make_abc_only <- function(combo) {
  crop_mask_to_ghana(ifel(combo == 7, 1, 0))
}

make_roi <- function(combo) {
  crop_mask_to_ghana(ifel((combo == 3) | (combo == 6) | (combo == 7), 1, 0))
}

# -----------------------------
# Build binary layers for all 5 scenarios
# -----------------------------
abc_layers <- lapply(scenario_combos, make_abc_only)
roi_layers <- lapply(scenario_combos, make_roi)

abc_layers <- abc_layers[!sapply(abc_layers, is.null)]
roi_layers <- roi_layers[!sapply(roi_layers, is.null)]

# -----------------------------
# Consensus rasters
# Values 0 to 5
# -----------------------------
consensus_abc <- Reduce(`+`, abc_layers)
consensus_roi <- Reduce(`+`, roi_layers)

consensus_abc <- crop_mask_to_ghana(consensus_abc)
consensus_roi <- crop_mask_to_ghana(consensus_roi)

# -----------------------------
# Convert to factors
# Values 0 to 5 remain raster classes
# Lakes are added separately as sf layer
# -----------------------------
consensus_abc_f <- as.factor(consensus_abc)
levels(consensus_abc_f) <- data.frame(
  value = 0:5,
  label = c("0", "1", "2", "3", "4", "5")
)

consensus_roi_f <- as.factor(consensus_roi)
levels(consensus_roi_f) <- data.frame(
  value = 0:5,
  label = c("0", "1", "2", "3", "4", "5")
)

# -----------------------------
# Colours
# 0 = grey
# 1-5 = increasing intensity
# lakes & lagoons = lightblue as separate sf layer
# -----------------------------
abc_cols <- c(
  "0" = "#D9D9D9",
  "1" = lighten(col_ABC, 0.55),
  "2" = lighten(col_ABC, 0.40),
  "3" = lighten(col_ABC, 0.25),
  "4" = lighten(col_ABC, 0.12),
  "5" = col_ABC,
  "lakes & lagoons" = "lightblue"
)

roi_cols <- c(
  "0" = "#D9D9D9",
  "1" = lighten(col_roi, 0.55),
  "2" = lighten(col_roi, 0.40),
  "3" = lighten(col_roi, 0.25),
  "4" = lighten(col_roi, 0.12),
  "5" = col_roi,
  "lakes & lagoons" = "lightblue"
)

# -----------------------------
# Output folder
# -----------------------------
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/With_Axes", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/No_Axes", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Base theme builders
# -----------------------------
theme_consensus_with_axes <- function() {
  theme_minimal() +
    theme(
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 18, face = "plain"),
      legend.text = element_text(size = 16),
      legend.key.width = unit(1.1, "cm"),
      legend.key.height = unit(0.75, "cm"),
      plot.margin = margin(12, 12, 12, 12),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

theme_consensus_no_axes <- function() {
  theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 18, face = "plain"),
      legend.text = element_text(size = 16),
      legend.key.width = unit(1.1, "cm"),
      legend.key.height = unit(0.75, "cm"),
      plot.margin = margin(12, 12, 12, 12),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# -----------------------------
# ABC consensus - WITH AXES
# -----------------------------
p_consensus_abc_axes <- ggplot() +
  geom_spatraster(data = consensus_abc_f, aes(fill = label)) +
  geom_sf(data = lakes_sf, aes(fill = "lakes & lagoons"), color = NA, inherit.aes = FALSE) +
  geom_sf(data = ghana_sf, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(
    values = abc_cols,
    breaks = c("5", "4", "3", "2", "1", "0", "lakes & lagoons"),
    drop = FALSE,
    na.value = "white",
    name = "Consensus score"
  ) +
  coord_sf(
    xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
    ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_consensus_with_axes()

ggsave(
  filename = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/With_Axes/Consensus_ABC_WithAxes_lakesblue_v2.png",
  plot = p_consensus_abc_axes,
  width = 11,
  height = 9,
  units = "in",
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

# -----------------------------
# ABC consensus - NO AXES
# -----------------------------
p_consensus_abc_noaxes <- ggplot() +
  geom_spatraster(data = consensus_abc_f, aes(fill = label)) +
  geom_sf(data = lakes_sf, aes(fill = "lakes & lagoons"), color = NA, inherit.aes = FALSE) +
  geom_sf(data = ghana_sf, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(
    values = abc_cols,
    breaks = c("5", "4", "3", "2", "1", "0", "lakes & lagoons"),
    drop = FALSE,
    na.value = "white",
    name = "Consensus score"
  ) +
  coord_sf(
    xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
    ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_consensus_no_axes()

ggsave(
  filename = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/No_Axes/Consensus_ABC_NoAxes_lakesblue_v2.png",
  plot = p_consensus_abc_noaxes,
  width = 11,
  height = 9,
  units = "in",
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

# -----------------------------
# ROI consensus - WITH AXES
# -----------------------------
p_consensus_roi_axes <- ggplot() +
  geom_spatraster(data = consensus_roi_f, aes(fill = label)) +
  geom_sf(data = lakes_sf, aes(fill = "lakes & lagoons"), color = NA, inherit.aes = FALSE) +
  geom_sf(data = ghana_sf, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(
    values = roi_cols,
    breaks = c("5", "4", "3", "2", "1", "0", "lakes & lagoons"),
    drop = FALSE,
    na.value = "white",
    name = "Consensus score"
  ) +
  coord_sf(
    xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
    ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_consensus_with_axes()

ggsave(
  filename = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/With_Axes/Consensus_RiskOfInterest_WithAxes_lakesblue_v2.png",
  plot = p_consensus_roi_axes,
  width = 11,
  height = 9,
  units = "in",
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

# -----------------------------
# ROI consensus - NO AXES
# -----------------------------
p_consensus_roi_noaxes <- ggplot() +
  geom_spatraster(data = consensus_roi_f, aes(fill = label)) +
  geom_sf(data = lakes_sf, aes(fill = "lakes & lagoons"), color = NA, inherit.aes = FALSE) +
  geom_sf(data = ghana_sf, fill = NA, color = "black", linewidth = 0.4, inherit.aes = FALSE) +
  scale_fill_manual(
    values = roi_cols,
    breaks = c("5", "4", "3", "2", "1", "0", "lakes & lagoons"),
    drop = FALSE,
    na.value = "white",
    name = "Consensus score"
  ) +
  coord_sf(
    xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
    ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  theme_consensus_no_axes()

ggsave(
  filename = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/No_Axes/Consensus_RiskOfInterest_NoAxes_lakesblue_v2.png",
  plot = p_consensus_roi_noaxes,
  width = 11,
  height = 9,
  units = "in",
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

cat("✅ Consensus maps created with axes and without axes for ABC and risks of interest.\n")
cat("✅ Lakes/lagoons are clipped to Ghana full boundary and added as separate lightblue sf layer.\n")


# -----------------------------
# 1 x 2 CONSENSUS PANEL
# KEEP INDIVIDUAL LEGENDS
# -----------------------------

library(patchwork)

dir.create(
  "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/Comparison_1x2",
  recursive = TRUE,
  showWarnings = FALSE
)

# Keep maps EXACTLY as they are
# including their own legends

p_consensus_abc_panel <- p_consensus_abc_noaxes +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  )

p_consensus_roi_panel <- p_consensus_roi_noaxes +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  )

# Build 1 x 2 panel
combined_consensus_panel <-
  
  p_consensus_abc_panel +
  patchwork::plot_spacer() +
  p_consensus_roi_panel +
  
  patchwork::plot_layout(
    ncol = 3,
    widths = c(1, 0.05, 1)
  ) &
  
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )

# Save
ggsave(
  filename =
    "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/Comparison_1x2/Consensus_ABC_ROI_1x2_SeparateLegends_lakesblue_v2.png",
  
  plot = combined_consensus_panel,
  
  width = 22,
  height = 9,
  
  units = "in",
  
  dpi = 400,
  
  bg = "white",
  
  limitsize = FALSE
)

cat("✅ 1 x 2 consensus panel with separate legends exported successfully.\n")













# ============================================================================== #
# 9. CONSENSUS REPORTS FOR THE TWO TYPES
# - ABC consensus
# - Risks of interest consensus
#
# CORRECTED TO USE:
#   - scenario_defs
#   - scenario_combos
# already created earlier
#
# Lakes/lagoons:
# - Uses lakes_sf_correct from the solved lake clipping step
# - No shp_lakes used for plotting
# - Report images point to LakesBlue_v2 outputs
# ============================================================================== #

library(terra)
library(sf)
library(dplyr)
library(flextable)
library(officer)
library(colorspace)

# ------------------------------------------------------------------------------ #
# Helper: build scenario binary layers from scenario_defs
# ------------------------------------------------------------------------------ #
scenario_bins <- lapply(scenario_defs, function(z) {
  S_bin <- create_bin(z$S, z$qS)
  B_bin <- create_bin(z$B, z$qB)
  C_bin <- create_bin(z$C, z$qC)
  
  list(S = S_bin, B = B_bin, C = C_bin)
})

# ------------------------------------------------------------------------------ #
# 4. NATIONAL VARIABLE CONTRIBUTION ANALYSIS (LOGIC ALIGNMENT)
# ------------------------------------------------------------------------------ #
analyze_vars_consensus <- function(S, P, G, model_name) {
  S_num <- as.numeric(S); P_num <- as.numeric(P); G_num <- as.numeric(G)
  S_num[is.na(S_num)] <- 0
  P_num[is.na(P_num)] <- 0
  G_num[is.na(G_num)] <- 0
  
  # Bitmask: S=1, G=2, P=4
  combo <- (S_num * 1) + (G_num * 2) + (P_num * 4)
  
  get_area <- function(code, label) {
    class_mask <- combo == code
    class_area <- area_ghana * class_mask
    km2 <- global(class_area, sum, na.rm = TRUE)[1,1]
    
    data.frame(
      Model = model_name,
      Variable = label,
      Area_km2 = round(km2, 2),
      Percentage = round((km2 / total_gh_area) * 100, 2)
    )
  }
  
  rbind(
    get_area(0, "all low"),
    get_area(1, "A = poor sanitation"),
    get_area(2, "B = high pig density"),
    get_area(3, "AB= poor sanitation + high pig density"),
    get_area(4, "C = high poverty"),
    get_area(5, "AC= poor sanitation + high poverty"),
    get_area(6, "BC= high pig density + high poverty"),
    get_area(7, "ABC= triple risk")
  )
}

table_1_full <- rbind(
  analyze_vars_consensus(
    scenario_bins$base$S, scenario_bins$base$C, scenario_bins$base$B,
    "Baseline (0.66)"
  ),
  analyze_vars_consensus(
    scenario_bins$test_A$S, scenario_bins$test_A$C, scenario_bins$test_A$B,
    "Test A (0.75)"
  ),
  analyze_vars_consensus(
    scenario_bins$test_B$S, scenario_bins$test_B$C, scenario_bins$test_B$B,
    "Test B (DHS Pigs)"
  ),
  analyze_vars_consensus(
    scenario_bins$test_C$S, scenario_bins$test_C$C, scenario_bins$test_C$B,
    "Test C (OD)"
  ),
  analyze_vars_consensus(
    scenario_bins$test_D$S, scenario_bins$test_D$C, scenario_bins$test_D$B,
    "Test D (W1)"
  )
)

# ------------------------------------------------------------------------------ #
# 7. DISTRICT ZONING RISK ANALYSIS (SYNCHRONIZED BITMASK)
# Done separately for ALL five scenario maps
# ------------------------------------------------------------------------------ #
analyze_district_zoning <- function(S, P, G, model_name) {
  S_num <- as.numeric(S); S_num[is.na(S_num)] <- 0
  P_num <- as.numeric(P); P_num[is.na(P_num)] <- 0
  G_num <- as.numeric(G); G_num[is.na(G_num)] <- 0
  
  combo_rast <- (S_num * 1) + (G_num * 2) + (P_num * 4)
  combo_rast <- mask(combo_rast, ghana_vect)
  
  dist_list <- list()
  
  for(i in 1:nrow(ghana_vect)) {
    poly <- ghana_vect[i,]
    d_name <- poly$NAME_2
    r_name <- poly$NAME_1
    
    d_combo <- crop(combo_rast, poly) %>% mask(poly)
    d_area  <- crop(area_rast, poly) %>% mask(poly)
    total_area <- global(d_area, sum, na.rm = TRUE)[1,1]
    
    for(j in 1:7) {
      target_id <- j
      mask_cat  <- d_combo == target_id
      cat_area  <- global(d_area * mask_cat, sum, na.rm = TRUE)[1,1]
      
      if(!is.na(cat_area) && cat_area > 0) {
        dist_list[[length(dist_list) + 1]] <- data.frame(
          Model = model_name,
          Region = r_name,
          District = d_name,
          Variable_Option = labels_custom[as.character(target_id)],
          Area_km2 = round(cat_area, 2),
          Area_Percentage = round((cat_area / total_area) * 100, 2)
        )
      }
    }
  }
  
  do.call(rbind, dist_list) %>% arrange(Region, District, desc(Area_Percentage))
}

final_zoning_baseline <- analyze_district_zoning(
  scenario_bins$base$S, scenario_bins$base$C, scenario_bins$base$B,
  "Baseline (0.66)"
)

final_zoning_testA <- analyze_district_zoning(
  scenario_bins$test_A$S, scenario_bins$test_A$C, scenario_bins$test_A$B,
  "Test A (0.75)"
)

final_zoning_testB <- analyze_district_zoning(
  scenario_bins$test_B$S, scenario_bins$test_B$C, scenario_bins$test_B$B,
  "Test B (DHS Pigs)"
)

final_zoning_testC <- analyze_district_zoning(
  scenario_bins$test_C$S, scenario_bins$test_C$C, scenario_bins$test_C$B,
  "Test C (OD)"
)

final_zoning_testD <- analyze_district_zoning(
  scenario_bins$test_D$S, scenario_bins$test_D$C, scenario_bins$test_D$B,
  "Test D (W1)"
)

# ------------------------------------------------------------------------------ #
# 8. DISTRICT CONSENSUS ANALYSIS
# Separate table for:
#   - ABC consensus
#   - Risks of interest consensus
# ------------------------------------------------------------------------------ #
analyze_district_consensus <- function(consensus_rast, ghana_vect, area_rast) {
  cons_list <- list()
  
  for(i in 1:nrow(ghana_vect)) {
    poly <- ghana_vect[i,]
    d_cons <- crop(consensus_rast, poly) %>% mask(poly)
    d_area <- crop(area_rast, poly) %>% mask(poly)
    tot_a  <- global(d_area, sum, na.rm = TRUE)[1,1]
    
    for(v in 1:5) {
      mask_v <- d_cons == v
      area_v <- global(d_area * mask_v, sum, na.rm = TRUE)[1,1]
      
      if(!is.na(area_v) && area_v > 0) {
        cons_list[[length(cons_list) + 1]] <- data.frame(
          Region = poly$NAME_1,
          District = poly$NAME_2,
          Consensus_Score = v,
          Area_km2 = round(area_v, 2),
          Pct_of_District = round((area_v / tot_a) * 100, 2)
        )
      }
    }
  }
  
  do.call(rbind, cons_list) %>% arrange(Region, District, desc(Consensus_Score))
}

final_cons_abc_table <- analyze_district_consensus(consensus_abc, ghana_vect, area_rast)
final_cons_roi_table <- analyze_district_consensus(consensus_roi, ghana_vect, area_rast)

# ------------------------------------------------------------------------------ #
# Reusable function for table coloring (consistent with risk_palette)
# ------------------------------------------------------------------------------ #
color_risk_table <- function(ft, col_name) {
  ft %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'A = poor sanitation'")), j = col_name, bg = risk_palette["1"]) %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'B = high pig density'")), j = col_name, bg = risk_palette["2"]) %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'AB'")), j = col_name, bg = risk_palette["3"]) %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'C = high poverty'")), j = col_name, bg = risk_palette["4"]) %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'AC'")), j = col_name, bg = risk_palette["5"]) %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'BC'")), j = col_name, bg = risk_palette["6"]) %>%
    bg(i = as.formula(paste0("~ ", col_name, " == 'ABC'")), j = col_name, bg = risk_palette["7"]) %>%
    color(
      i = as.formula(paste0("~ ", col_name, " %in% c('B = high pig density', 'C = high poverty', 'ABC')")),
      j = col_name,
      color = "white"
    )
}

# ------------------------------------------------------------------------------ #
# Flextables: national and district zoning
# ------------------------------------------------------------------------------ #
ft_national <- color_risk_table(flextable(table_1_full), "Variable") %>%
  theme_vanilla() %>%
  merge_v(j = ~ Model) %>%
  autofit()

ft_zoning_baseline <- color_risk_table(flextable(final_zoning_baseline), "Variable_Option") %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  autofit()

ft_zoning_testA <- color_risk_table(flextable(final_zoning_testA), "Variable_Option") %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  autofit()

ft_zoning_testB <- color_risk_table(flextable(final_zoning_testB), "Variable_Option") %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  autofit()

ft_zoning_testC <- color_risk_table(flextable(final_zoning_testC), "Variable_Option") %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  autofit()

ft_zoning_testD <- color_risk_table(flextable(final_zoning_testD), "Variable_Option") %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  autofit()

# ------------------------------------------------------------------------------ #
# Flextables: consensus tables using the SAME colour scales as the consensus maps
# ------------------------------------------------------------------------------ #
ft_consensus_abc <- flextable(final_cons_abc_table) %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  bg(i = ~ Consensus_Score == 5, j = "Consensus_Score", bg = col_ABC) %>%
  bg(i = ~ Consensus_Score == 4, j = "Consensus_Score", bg = lighten(col_ABC, 0.12)) %>%
  bg(i = ~ Consensus_Score == 3, j = "Consensus_Score", bg = lighten(col_ABC, 0.25)) %>%
  bg(i = ~ Consensus_Score == 2, j = "Consensus_Score", bg = lighten(col_ABC, 0.40)) %>%
  bg(i = ~ Consensus_Score == 1, j = "Consensus_Score", bg = lighten(col_ABC, 0.55)) %>%
  color(i = ~ Consensus_Score >= 4, j = "Consensus_Score", color = "white") %>%
  autofit()

ft_consensus_roi <- flextable(final_cons_roi_table) %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region + District) %>%
  bg(i = ~ Consensus_Score == 5, j = "Consensus_Score", bg = col_roi) %>%
  bg(i = ~ Consensus_Score == 4, j = "Consensus_Score", bg = lighten(col_roi, 0.12)) %>%
  bg(i = ~ Consensus_Score == 3, j = "Consensus_Score", bg = lighten(col_roi, 0.25)) %>%
  bg(i = ~ Consensus_Score == 2, j = "Consensus_Score", bg = lighten(col_roi, 0.40)) %>%
  bg(i = ~ Consensus_Score == 1, j = "Consensus_Score", bg = lighten(col_roi, 0.55)) %>%
  color(i = ~ Consensus_Score >= 4, j = "Consensus_Score", color = "white") %>%
  autofit()

# ------------------------------------------------------------------------------ #
# Output folders
# ------------------------------------------------------------------------------ #
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/07_Consensus_Reports/ABC", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_2022_LakesBlue_v2/07_Consensus_Reports/RiskOfInterest", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------ #
# ABC consensus report
# ------------------------------------------------------------------------------ #
doc_consensus_abc <- read_docx() %>%
  body_add_par("Ghana 2022: ABC Consensus Report", style = "heading 1") %>%
  body_add_par("1. National Variable Contribution", style = "heading 2") %>%
  body_add_flextable(ft_national) %>%
  body_add_par("2. ABC consensus map with axes", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/With_Axes/Consensus_ABC_WithAxes_lakesblue_v2.png",
    width = 5.5,
    height = 6.5
  ) %>%
  body_add_par("3. ABC consensus map without axes", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/No_Axes/Consensus_ABC_NoAxes_lakesblue_v2.png",
    width = 5.5,
    height = 6.5
  ) %>%
  body_add_par("4. District Zoning: Baseline (0.66)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_baseline) %>%
  body_add_par("5. District Zoning: Test A (0.75)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testA) %>%
  body_add_par("6. District Zoning: Test B (DHS Pigs)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testB) %>%
  body_add_par("7. District Zoning: Test C (OD)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testC) %>%
  body_add_par("8. District Zoning: Test D (W1)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testD) %>%
  body_add_par("9. District Area by Consensus Score (ABC only)", style = "heading 2") %>%
  body_add_flextable(ft_consensus_abc)

print(
  doc_consensus_abc,
  target = "outputs/Sensitivity_2022_LakesBlue_v2/07_Consensus_Reports/ABC/Ghana_2022_ABC_Consensus_Report_lakesblue_v2.docx"
)

# ------------------------------------------------------------------------------ #
# Risks of interest consensus report
# ------------------------------------------------------------------------------ #
doc_consensus_roi <- read_docx() %>%
  body_add_par("Ghana 2022: Risks of Interest Consensus Report", style = "heading 1") %>%
  body_add_par("1. National Variable Contribution", style = "heading 2") %>%
  body_add_flextable(ft_national) %>%
  body_add_par("2. Risks of interest consensus map with axes", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/With_Axes/Consensus_RiskOfInterest_WithAxes_lakesblue_v2.png",
    width = 5.5,
    height = 6.5
  ) %>%
  body_add_par("3. Risks of interest consensus map without axes", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Sensitivity_2022_LakesBlue_v2/06_Consensus/No_Axes/Consensus_RiskOfInterest_NoAxes_lakesblue_v2.png",
    width = 5.5,
    height = 6.5
  ) %>%
  body_add_par("4. District Zoning: Baseline (0.66)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_baseline) %>%
  body_add_par("5. District Zoning: Test A (0.75)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testA) %>%
  body_add_par("6. District Zoning: Test B (DHS Pigs)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testB) %>%
  body_add_par("7. District Zoning: Test C (OD)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testC) %>%
  body_add_par("8. District Zoning: Test D (W1)", style = "heading 2") %>%
  body_add_flextable(ft_zoning_testD) %>%
  body_add_par("9. District Area by Consensus Score (Risks of interest)", style = "heading 2") %>%
  body_add_flextable(ft_consensus_roi)

print(
  doc_consensus_roi,
  target = "outputs/Sensitivity_2022_LakesBlue_v2/07_Consensus_Reports/RiskOfInterest/Ghana_2022_RiskOfInterest_Consensus_Report_lakesblue_v2.docx"
)

cat("\n--- CONSENSUS REPORTS COMPLETE: ABC AND RISKS OF INTEREST EXPORTED SUCCESSFULLY ---\n")








# ============================================================================== #
# 5. FINAL ROBUST EXPORT: PNG HEATMAPS + WORD TABLES (STATISTICAL CALCULATION)
# For:
#   1. Full maps
#   2. ABC-only maps
#   3. Risks-of-interest maps
# Uses existing scenario_combos already created earlier
# New output names so old files are not replaced
# ============================================================================== #

library(tibble)
library(scales)
library(flextable)
library(officer)
library(dplyr)
library(irr)
library(tidyr)
library(ggplot2)
library(terra)

# ------------------------------------------------------------------------------ #
# Helper layers from combo
# ------------------------------------------------------------------------------ #
make_full_map <- function(combo) {
  crop_mask_to_ghana(combo)
}

make_abc_only <- function(combo) {
  crop_mask_to_ghana(ifel(combo == 7, 1, 0))
}

make_roi <- function(combo) {
  crop_mask_to_ghana(ifel((combo == 3) | (combo == 6) | (combo == 7), 1, 0))
}

# ------------------------------------------------------------------------------ #
# 5.1. Setup Comparison Lists
# ------------------------------------------------------------------------------ #
all_models_full <- list(
  "Base"   = make_full_map(scenario_combos$base),
  "Test A" = make_full_map(scenario_combos$test_A),
  "Test B" = make_full_map(scenario_combos$test_B),
  "Test C" = make_full_map(scenario_combos$test_C),
  "Test D" = make_full_map(scenario_combos$test_D)
)

all_models_abc <- list(
  "Base"   = make_abc_only(scenario_combos$base),
  "Test A" = make_abc_only(scenario_combos$test_A),
  "Test B" = make_abc_only(scenario_combos$test_B),
  "Test C" = make_abc_only(scenario_combos$test_C),
  "Test D" = make_abc_only(scenario_combos$test_D)
)

all_models_roi <- list(
  "Base"   = make_roi(scenario_combos$base),
  "Test A" = make_roi(scenario_combos$test_A),
  "Test B" = make_roi(scenario_combos$test_B),
  "Test C" = make_roi(scenario_combos$test_C),
  "Test D" = make_roi(scenario_combos$test_D)
)

m_names <- names(all_models_full)
n <- length(m_names)

# ------------------------------------------------------------------------------ #
# Helper: calculate matrices
# ------------------------------------------------------------------------------ #
calc_similarity_mats <- function(all_models, m_names, ghana_vect, binary_mode = TRUE) {
  n <- length(m_names)
  
  kap_mat <- matrix(NA, n, n, dimnames = list(m_names, m_names))
  jac_mat <- matrix(NA, n, n, dimnames = list(m_names, m_names))
  spr_mat <- matrix(NA, n, n, dimnames = list(m_names, m_names))
  
  for(i in 1:n) {
    for(j in 1:n) {
      v1 <- values(mask(all_models[[i]], ghana_vect), mat = FALSE)
      v2 <- values(mask(all_models[[j]], ghana_vect), mat = FALSE)
      
      valid <- complete.cases(v1, v2)
      v1_c <- v1[valid]
      v2_c <- v2[valid]
      
      kap_mat[i, j] <- kappa2(data.frame(v1_c, v2_c))$value
      
      if (binary_mode) {
        union_sum <- sum(v1_c == 1 | v2_c == 1)
        jac_mat[i, j] <- ifelse(
          union_sum == 0,
          NA,
          sum(v1_c == 1 & v2_c == 1) / union_sum
        )
      } else {
        # Multi-class exact match proportion for full maps
        jac_mat[i, j] <- sum(v1_c == v2_c) / length(v1_c)
      }
      
      spr_mat[i, j] <- suppressWarnings(cor(v1_c, v2_c, method = "spearman"))
    }
  }
  
  list(kappa = kap_mat, jaccard = jac_mat, spearman = spr_mat)
}

# ------------------------------------------------------------------------------ #
# 5.2. Calculate matrices
# ------------------------------------------------------------------------------ #
stats_full <- calc_similarity_mats(all_models_full, m_names, ghana_vect, binary_mode = FALSE)
stats_abc  <- calc_similarity_mats(all_models_abc,  m_names, ghana_vect, binary_mode = TRUE)
stats_roi  <- calc_similarity_mats(all_models_roi,  m_names, ghana_vect, binary_mode = TRUE)

# ------------------------------------------------------------------------------ #
# 5.3. GGPLOT2 HEATMAP FUNCTION (For PNGs)
# ------------------------------------------------------------------------------ #
save_heatmap_ggplot <- function(mat, title, filename) {
  df_long <- as.data.frame(mat) %>%
    rownames_to_column(var = "Model_A") %>%
    pivot_longer(-Model_A, names_to = "Model_B", values_to = "Value")
  
  df_long$Model_A <- factor(df_long$Model_A, levels = m_names)
  df_long$Model_B <- factor(df_long$Model_B, levels = m_names)
  
  p <- ggplot(df_long, aes(x = Model_A, y = Model_B, fill = Value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(Value), "", sprintf("%.3f", Value)), color = Value > 0.7), size = 4) +
    scale_fill_gradientn(
      colors = c("#ffffff", "#fee5d9", "#fb6a4a", "#cb181d"),
      limits = c(0, 1),
      na.value = "grey90"
    ) +
    scale_color_manual(values = c("black", "white"), guide = "none", na.translate = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename, p, width = 8, height = 7, dpi = 300)
}

# ------------------------------------------------------------------------------ #
# 5.4. WORD TABLE FUNCTION (Values Only)
# ------------------------------------------------------------------------------ #
make_word_table <- function(mat, title) {
  df <- as.data.frame(mat) %>%
    rownames_to_column(var = "Model")
  
  flextable(df) %>%
    theme_vanilla() %>%
    set_caption(title) %>%
    colformat_double(digits = 3) %>%
    autofit()
}

# ============================================================================== #
# EXECUTE EXPORTS
# ============================================================================== #
dir.create("outputs/Sensitivity_Analysis_LakesBlue_v2/Full", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_Analysis_LakesBlue_v2/ABC", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Sensitivity_Analysis_LakesBlue_v2/RiskOfInterest", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------ #
# 1. SAVE PNGs - FULL MAPS
# ------------------------------------------------------------------------------ #
save_heatmap_ggplot(
  stats_full$kappa,
  "Kappa Matrix - Full Maps",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/Full/Heatmap_Kappa_Full_lakesblue_v2.png"
)

save_heatmap_ggplot(
  stats_full$jaccard,
  "Exact Match Matrix - Full Maps",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/Full/Heatmap_ExactMatch_Full_lakesblue_v2.png"
)

save_heatmap_ggplot(
  stats_full$spearman,
  "Spearman Matrix - Full Maps",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/Full/Heatmap_Spearman_Full_lakesblue_v2.png"
)

# ------------------------------------------------------------------------------ #
# 2. SAVE PNGs - ABC
# ------------------------------------------------------------------------------ #
save_heatmap_ggplot(
  stats_abc$kappa,
  "Kappa Matrix - ABC",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/ABC/Heatmap_Kappa_ABC_lakesblue_v2.png"
)

save_heatmap_ggplot(
  stats_abc$jaccard,
  "Jaccard Matrix - ABC",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/ABC/Heatmap_Jaccard_ABC_lakesblue_v2.png"
)

save_heatmap_ggplot(
  stats_abc$spearman,
  "Spearman Matrix - ABC",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/ABC/Heatmap_Spearman_ABC_lakesblue_v2.png"
)

# ------------------------------------------------------------------------------ #
# 3. SAVE PNGs - Risk of Interest
# ------------------------------------------------------------------------------ #
save_heatmap_ggplot(
  stats_roi$kappa,
  "Kappa Matrix - Risks of Interest",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/RiskOfInterest/Heatmap_Kappa_RiskOfInterest_lakesblue_v2.png"
)

save_heatmap_ggplot(
  stats_roi$jaccard,
  "Jaccard Matrix - Risks of Interest",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/RiskOfInterest/Heatmap_Jaccard_RiskOfInterest_lakesblue_v2.png"
)

save_heatmap_ggplot(
  stats_roi$spearman,
  "Spearman Matrix - Risks of Interest",
  "outputs/Sensitivity_Analysis_LakesBlue_v2/RiskOfInterest/Heatmap_Spearman_RiskOfInterest_lakesblue_v2.png"
)

# ------------------------------------------------------------------------------ #
# 4. SAVE WORD DOCUMENT - FULL MAPS
# ------------------------------------------------------------------------------ #
doc_stats_full <- read_docx() %>%
  body_add_par("T. solium Sensitivity Analysis: Statistical Tables for Full Maps", style = "heading 1") %>%
  
  body_add_par("1. Kappa Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_full$kappa, "Kappa Matrix - Full Maps")) %>%
  
  body_add_par("2. Exact Match Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_full$jaccard, "Exact Match Matrix - Full Maps")) %>%
  
  body_add_par("3. Spearman Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_full$spearman, "Spearman Matrix - Full Maps")) %>%
  
  body_end_section_landscape()

print(
  doc_stats_full,
  target = "outputs/Sensitivity_Analysis_LakesBlue_v2/Full/Ghana_Sensitivity_Tables_Full_lakesblue_v2.docx"
)

# ------------------------------------------------------------------------------ #
# 5. SAVE WORD DOCUMENT - ABC
# ------------------------------------------------------------------------------ #
doc_stats_abc <- read_docx() %>%
  body_add_par("T. solium Sensitivity Analysis: Statistical Tables for ABC", style = "heading 1") %>%
  
  body_add_par("1. Kappa Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_abc$kappa, "Kappa Matrix - ABC")) %>%
  
  body_add_par("2. Jaccard Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_abc$jaccard, "Jaccard Matrix - ABC")) %>%
  
  body_add_par("3. Spearman Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_abc$spearman, "Spearman Matrix - ABC")) %>%
  
  body_end_section_landscape()

print(
  doc_stats_abc,
  target = "outputs/Sensitivity_Analysis_LakesBlue_v2/ABC/Ghana_Sensitivity_Tables_ABC_lakesblue_v2.docx"
)

# ------------------------------------------------------------------------------ #
# 6. SAVE WORD DOCUMENT - Risks of Interest
# ------------------------------------------------------------------------------ #
doc_stats_roi <- read_docx() %>%
  body_add_par("T. solium Sensitivity Analysis: Statistical Tables for Risks of Interest", style = "heading 1") %>%
  
  body_add_par("1. Kappa Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_roi$kappa, "Kappa Matrix - Risks of Interest")) %>%
  
  body_add_par("2. Jaccard Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_roi$jaccard, "Jaccard Matrix - Risks of Interest")) %>%
  
  body_add_par("3. Spearman Data", style = "heading 2") %>%
  body_add_flextable(make_word_table(stats_roi$spearman, "Spearman Matrix - Risks of Interest")) %>%
  
  body_end_section_landscape()

print(
  doc_stats_roi,
  target = "outputs/Sensitivity_Analysis_LakesBlue_v2/RiskOfInterest/Ghana_Sensitivity_Tables_RiskOfInterest_lakesblue_v2.docx"
)

cat("\n--- DONE: PNGs and Landscape Word Tables Saved for Full maps, ABC, and Risks of Interest ---\n")








# ============================================================================== #
# PIG DISAGREEMENT MAPS + PERCENTAGE TABLES (WORD)
# ============================================================================== #

library(terra)
library(sf)
library(ggplot2)
library(tidyterra)
library(dplyr)
library(flextable)
library(officer)

# ------------------------------------------------------------------------------ #
# OUTPUT FOLDER
# ------------------------------------------------------------------------------ #
dir.create("outputs/Disagreement_Analysis_LakesBlue_v2", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------ #
# EXACT COLOURS
# ------------------------------------------------------------------------------ #
col_low_dis     <- "#D9D9D9"   # all low & other risks
col_dhs_only    <- "#000080"   # navy blue
col_glw_only    <- "#E30B5C"   # raspberry

# For base/test B comparisons, keep the same one-only colours consistently:
col_base_only   <- "#000080"   # navy blue
col_testB_only  <- "#E30B5C"   # raspberry

# Agreement colours use the exact colours already in your workflow:
# binary pigs agreement = exact pig colour
# ABC agreement = exact ABC brown
# ROI agreement = exact ROI colour
# expects these objects already exist:
# col_B, col_ABC, col_roi

# ------------------------------------------------------------------------------ #
# HELPER: disagreement raster
# 0 = all low & other risks
# 1 = map 1 only
# 2 = map 2 only
# 3 = both / agreement
# ------------------------------------------------------------------------------ #
make_disagreement <- function(r1, r2) {
  r1n <- as.numeric(r1)
  r2n <- as.numeric(r2)
  
  r1n[is.na(r1n)] <- 0
  r2n[is.na(r2n)] <- 0
  
  out <- (r1n == 1 & r2n == 0) * 1 +
    (r1n == 0 & r2n == 1) * 2 +
    (r1n == 1 & r2n == 1) * 3
  
  crop_mask_to_ghana(out)
}

# ------------------------------------------------------------------------------ #
# HELPER: percentage table
# IMPORTANT:
# percentages are calculated only from:
#   1 = map 1 only
#   2 = map 2 only
#   3 = both
# class 0 is excluded from the denominator
# ------------------------------------------------------------------------------ #
calc_disagreement_pct <- function(rast, label_map1, label_map2, label_both) {
  v <- values(rast, mat = FALSE)
  v <- v[!is.na(v)]
  
  active <- v[v %in% c(1, 2, 3)]
  total_active <- length(active)
  
  if (total_active == 0) {
    return(data.frame(
      Category = c(label_map1, label_map2, label_both),
      Percentage = c(0, 0, 0)
    ))
  }
  
  data.frame(
    Category = c(label_map1, label_map2, label_both),
    Percentage = c(
      round(sum(active == 1) / total_active * 100, 2),
      round(sum(active == 2) / total_active * 100, 2),
      round(sum(active == 3) / total_active * 100, 2)
    )
  )
}

# ------------------------------------------------------------------------------ #
# HELPER: plot disagreement map
# no axes, no title, legend included
# Lakes added as separate sf layer, not rasterized
# ------------------------------------------------------------------------------ #
plot_disagreement_map <- function(rast, cols, labels, breaks_order, legend_title = "Conflict") {
  r_fact <- as.factor(rast)
  levels(r_fact) <- data.frame(
    value = 0:3,
    label = labels[1:4]
  )
  
  ggplot() +
    geom_spatraster(data = r_fact, aes(fill = label)) +
    geom_sf(data = lakes_sf, aes(fill = "lakes & lagoons"), colour = NA, inherit.aes = FALSE) +
    geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.4, inherit.aes = FALSE) +
    scale_fill_manual(
      values = cols,
      breaks = breaks_order,
      drop = FALSE,
      na.value = "white",
      name = legend_title
    ) +
    coord_sf(
      xlim = st_bbox(ghana_sf)[c("xmin", "xmax")],
      ylim = st_bbox(ghana_sf)[c("ymin", "ymax")],
      expand = FALSE
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 18, face = "plain"),
      legend.text = element_text(size = 16),
      legend.key.width = unit(1.1, "cm"),
      legend.key.height = unit(0.75, "cm"),
      plot.margin = margin(12, 12, 12, 12),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# ------------------------------------------------------------------------------ #
# 1. BINARY MAP OF GLW4 AND PIG DHS 2022
# agreement colour = exact pig colour (col_B)
# map1 = GLW4
# map2 = DHS
# ------------------------------------------------------------------------------ #
glw_bin <- create_bin(pig_glw4, 0.666)
dhs_bin <- create_bin(pig_dhs_mean, 0.666)

dis_glw_dhs <- make_disagreement(glw_bin, dhs_bin)

labels_glw_dhs <- c(
  "all low & other risks",
  "GLW4 only",
  "DHS only",
  "Both",
  "lakes & lagoons"
)

breaks_glw_dhs <- c(
  "Both",
  "GLW4 only",
  "DHS only",
  "all low & other risks",
  "lakes & lagoons"
)

cols_glw_dhs <- c(
  "all low & other risks" = col_low_dis,
  "GLW4 only" = col_glw_only,
  "DHS only" = col_dhs_only,
  "Both" = col_B,
  "lakes & lagoons" = "lightblue"
)

p_glw_dhs <- plot_disagreement_map(
  dis_glw_dhs,
  cols = cols_glw_dhs,
  labels = labels_glw_dhs,
  breaks_order = breaks_glw_dhs,
  legend_title = "Pig data source"
)

ggsave(
  filename = "outputs/Disagreement_Analysis_LakesBlue_v2/Disagreement_GLW4_vs_DHS_lakesblue_v2.png",
  plot = p_glw_dhs,
  width = 8,
  height = 10,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

pct_glw_dhs <- calc_disagreement_pct(
  dis_glw_dhs,
  label_map1 = "GLW4 only",
  label_map2 = "DHS only",
  label_both = "Both"
)

# ------------------------------------------------------------------------------ #
# 2. ABC OF BASE MAP AND ABC OF TEST B
# agreement colour = exact ABC brown (col_ABC)
# map1 = Base
# map2 = Test B
# ------------------------------------------------------------------------------ #
abc_base  <- make_abc_only(scenario_combos$base)
abc_testB <- make_abc_only(scenario_combos$test_B)

dis_abc <- make_disagreement(abc_base, abc_testB)

labels_abc <- c(
  "all low & other risks",
  "GLW4 only",
  "DHS only",
  "Both",
  "lakes & lagoons"
)

breaks_abc <- c(
  "Both",
  "GLW4 only",
  "DHS only",
  "all low & other risks",
  "lakes & lagoons"
)

cols_abc <- c(
  "all low & other risks" = col_low_dis,
  "GLW4 only" = col_glw_only,
  "DHS only" = col_dhs_only,
  "Both" = col_ABC,
  "lakes & lagoons" = "lightblue"
)

p_abc <- plot_disagreement_map(
  dis_abc,
  cols = cols_abc,
  labels = labels_abc,
  breaks_order = breaks_abc,
  legend_title = "Pig data source"
)

ggsave(
  filename = "outputs/Disagreement_Analysis_LakesBlue_v2/Disagreement_ABC_Base_vs_TestB_lakesblue_v2.png",
  plot = p_abc,
  width = 8,
  height = 10,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

pct_abc <- calc_disagreement_pct(
  dis_abc,
  label_map1 = "GLW4 only",
  label_map2 = "DHS only",
  label_both = "Both"
)

# ------------------------------------------------------------------------------ #
# 3. RISK OF INTEREST OF BASE MAP AND RISK OF INTEREST OF TEST B
# agreement colour = exact ROI colour (col_roi)
# map1 = Base
# map2 = Test B
# ------------------------------------------------------------------------------ #
roi_base  <- make_roi(scenario_combos$base)
roi_testB <- make_roi(scenario_combos$test_B)

dis_roi <- make_disagreement(roi_base, roi_testB)

labels_roi <- c(
  "all low & other risks",
  "GLW4 only",
  "DHS only",
  "Both",
  "lakes & lagoons"
)

breaks_roi <- c(
  "Both",
  "GLW4 only",
  "DHS only",
  "all low & other risks",
  "lakes & lagoons"
)

cols_roi <- c(
  "all low & other risks" = col_low_dis,
  "GLW4 only" = col_glw_only,
  "DHS only" = col_dhs_only,
  "Both" = col_roi,
  "lakes & lagoons" = "lightblue"
)

p_roi <- plot_disagreement_map(
  dis_roi,
  cols = cols_roi,
  labels = labels_roi,
  breaks_order = breaks_roi,
  legend_title = "Pig data source"
)

ggsave(
  filename = "outputs/Disagreement_Analysis_LakesBlue_v2/Disagreement_RiskOfInterest_Base_vs_TestB_lakesblue_v2.png",
  plot = p_roi,
  width = 8,
  height = 10,
  dpi = 400,
  bg = "white",
  limitsize = FALSE
)

pct_roi <- calc_disagreement_pct(
  dis_roi,
  label_map1 = "GLW4 only",
  label_map2 = "DHS only",
  label_both = "Both"
)

# ------------------------------------------------------------------------------ #
# WORD TABLES
# ------------------------------------------------------------------------------ #
ft_glw_dhs <- flextable(pct_glw_dhs) %>%
  theme_vanilla() %>%
  colformat_double(j = "Percentage", digits = 2) %>%
  autofit()

ft_abc <- flextable(pct_abc) %>%
  theme_vanilla() %>%
  colformat_double(j = "Percentage", digits = 2) %>%
  autofit()

ft_roi <- flextable(pct_roi) %>%
  theme_vanilla() %>%
  colformat_double(j = "Percentage", digits = 2) %>%
  autofit()

# ------------------------------------------------------------------------------ #
# WORD DOCUMENT
# ------------------------------------------------------------------------------ #
doc_disagreement <- read_docx() %>%
  body_add_par("Ghana 2022: Disagreement Analysis", style = "heading 1") %>%
  
  body_add_par("1. Binary pig map: GLW4 vs DHS 2022", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Disagreement_Analysis_LakesBlue_v2/Disagreement_GLW4_vs_DHS_lakesblue_v2.png",
    width = 5.8,
    height = 7.2
  ) %>%
  body_add_par("Percentages", style = "heading 3") %>%
  body_add_flextable(ft_glw_dhs) %>%
  
  body_add_par("2. ABC of base map vs ABC of test B", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Disagreement_Analysis_LakesBlue_v2/Disagreement_ABC_Base_vs_TestB_lakesblue_v2.png",
    width = 5.8,
    height = 7.2
  ) %>%
  body_add_par("Percentages", style = "heading 3") %>%
  body_add_flextable(ft_abc) %>%
  
  body_add_par("3. Risk of interest of base map vs risk of interest of test B", style = "heading 2") %>%
  body_add_img(
    src = "outputs/Disagreement_Analysis_LakesBlue_v2/Disagreement_RiskOfInterest_Base_vs_TestB_lakesblue_v2.png",
    width = 5.8,
    height = 7.2
  ) %>%
  body_add_par("Percentages", style = "heading 3") %>%
  body_add_flextable(ft_roi)

print(
  doc_disagreement,
  target = "outputs/Disagreement_Analysis_LakesBlue_v2/Ghana_2022_Disagreement_Analysis_lakesblue_v2.docx"
)

cat("\n--- DISAGREEMENT MAPS AND WORD DOCUMENT EXPORTED SUCCESSFULLY ---\n")

# ------------------------------------------------------------------------------ #
# 1 x 3 DISAGREEMENT PANEL
# KEEP INDIVIDUAL LEGENDS
# ------------------------------------------------------------------------------ #

library(patchwork)

dir.create(
  "outputs/Disagreement_Analysis_LakesBlue_v2/Comparison_1x3",
  recursive = TRUE,
  showWarnings = FALSE
)

# Keep maps EXACTLY as they are
# including their own legends

p_glw_dhs_panel <- p_glw_dhs +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  )

p_abc_panel <- p_abc +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  )

p_roi_panel <- p_roi +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  )

# Build 1 x 3 panel
combined_disagreement_panel <-
  
  p_glw_dhs_panel +
  patchwork::plot_spacer() +
  p_abc_panel +
  patchwork::plot_spacer() +
  p_roi_panel +
  
  patchwork::plot_layout(
    ncol = 5,
    widths = c(1, 0.05, 1, 0.05, 1)
  ) &
  
  theme(
    plot.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    plot.margin = margin(
      12,
      12,
      12,
      12
    )
  )

# Save
ggsave(
  filename =
    "outputs/Disagreement_Analysis_LakesBlue_v2/Comparison_1x3/Disagreement_1x3_SeparateLegends_lakesblue_v2.png",
  
  plot = combined_disagreement_panel,
  
  width = 30,
  height = 10,
  
  units = "in",
  
  dpi = 400,
  
  bg = "white",
  
  limitsize = FALSE
)

cat("✅ 1 x 3 disagreement panel with separate legends exported successfully.\n")







# ============================================================================== #
# FINAL TABLE: AREA (km²) + PERCENTAGE + CHANGE (ABC & ROI)
# LAND ONLY (EXCLUDES LAKES)
# ============================================================================== #

library(terra)
library(dplyr)
library(flextable)
library(officer)

dir.create("outputs/Sensitivity_2022_LakesBlue_v2", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------ #
# 1. LAND AREA (km²)
# ------------------------------------------------------------------------------ #

area_rast <- terra::cellSize(scenario_combos$base, unit = "km")
area_land <- terra::mask(area_rast, ghana_vect)

total_land_area <- terra::global(area_land, sum, na.rm = TRUE)[1, 1]

# ------------------------------------------------------------------------------ #
# 2. FUNCTION: AREA + %
# ------------------------------------------------------------------------------ #

calc_area_pct <- function(combo_rast, target_classes) {
  
  r <- terra::crop(combo_rast, ghana_vect, snap = "out")
  r <- terra::mask(r, ghana_vect)
  
  target_mask <- r %in% target_classes
  
  area_km2 <- terra::global(area_land * target_mask, sum, na.rm = TRUE)[1, 1]
  pct <- (area_km2 / total_land_area) * 100
  
  return(list(area = area_km2, pct = pct))
}

# ------------------------------------------------------------------------------ #
# 3. SCENARIOS
# ------------------------------------------------------------------------------ #

scenario_labels <- c(
  "base"   = "Base",
  "test_A" = "Test A",
  "test_B" = "Test B",
  "test_C" = "Test C",
  "test_D" = "Test D"
)

# ------------------------------------------------------------------------------ #
# 4. BUILD TABLE
# ------------------------------------------------------------------------------ #

results <- lapply(names(scenario_labels), function(sc) {
  
  abc <- calc_area_pct(scenario_combos[[sc]], 7)
  roi <- calc_area_pct(scenario_combos[[sc]], c(3,6,7))
  
  data.frame(
    TEST = scenario_labels[[sc]],
    
    ABC_AREA_km2 = abc$area,
    ABC_PERCENT  = abc$pct,
    
    ROI_AREA_km2 = roi$area,
    ROI_PERCENT  = roi$pct
  )
})

coverage_table <- bind_rows(results)

# ------------------------------------------------------------------------------ #
# 5. BASE VALUES
# ------------------------------------------------------------------------------ #

base_abc_area <- coverage_table$ABC_AREA_km2[coverage_table$TEST == "Base"]
base_roi_area <- coverage_table$ROI_AREA_km2[coverage_table$TEST == "Base"]

base_abc_pct  <- coverage_table$ABC_PERCENT[coverage_table$TEST == "Base"]
base_roi_pct  <- coverage_table$ROI_PERCENT[coverage_table$TEST == "Base"]

# ------------------------------------------------------------------------------ #
# 6. ADD CHANGES
# ------------------------------------------------------------------------------ #

coverage_table <- coverage_table %>%
  mutate(
    # AREA CHANGE (km²)
    ABC_AREA_CHANGE_km2 = ifelse(TEST == "Base", NA, ABC_AREA_km2 - base_abc_area),
    ROI_AREA_CHANGE_km2 = ifelse(TEST == "Base", NA, ROI_AREA_km2 - base_roi_area),
    
    # % CHANGE (relative)
    ABC_PERCENT_CHANGE = ifelse(TEST == "Base", NA,
                                ((ABC_PERCENT - base_abc_pct) / base_abc_pct) * 100),
    
    ROI_PERCENT_CHANGE = ifelse(TEST == "Base", NA,
                                ((ROI_PERCENT - base_roi_pct) / base_roi_pct) * 100)
  ) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 2))
  )

# ------------------------------------------------------------------------------ #
# 7. FORMAT FOR OUTPUT
# ------------------------------------------------------------------------------ #

coverage_table_print <- coverage_table
coverage_table_print[is.na(coverage_table_print)] <- "-"

print(coverage_table_print)

# ------------------------------------------------------------------------------ #
# 8. EXPORT CSV
# ------------------------------------------------------------------------------ #

write.csv(
  coverage_table_print,
  "outputs/Sensitivity_2022_LakesBlue_v2/ABC_ROI_Area_Percentage_Table.csv",
  row.names = FALSE
)

# ------------------------------------------------------------------------------ #
# 9. EXPORT WORD
# ------------------------------------------------------------------------------ #

ft <- flextable(coverage_table_print) %>%
  theme_vanilla() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("ABC and Risk of Interest: Area and Percentage (Land Only)", style = "heading 1") %>%
  body_add_par("Includes absolute area (km²), percentage of land, and change relative to Base.", style = "Normal") %>%
  body_add_flextable(ft)

print(
  doc,
  target = "outputs/Sensitivity_2022_LakesBlue_v2/ABC_ROI_Area_Percentage_Table.docx"
)

cat("✅ Done: Area + percentage + change table created (land only, lakes excluded).\n")


# ==============================================================================
# STANDALONE REPORT: GHANA DISTRICT MODAL RISK (BASE MAP)
# Only legend corrected; results and image quality unchanged
# ==============================================================================

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(officer)
library(flextable)
library(grid)

out_dir <- "outputs/Sensitivity_2022_LakesBlue_v2/Modal_Base"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

combo_rast <- scenario_combos$base
combo_rast <- crop_mask_to_ghana(combo_rast)
combo_rast <- terra::mask(combo_rast, ghana_vect)

risk_labels <- c(
  "0" = "all low",
  "1" = "A = poor sanitation",
  "2" = "B = high pig density",
  "3" = "AB = poor sanitation + high pig density",
  "4" = "C = high poverty",
  "5" = "AC = poor sanitation + high poverty",
  "6" = "BC = high pig density + high poverty",
  "7" = "ABC = poor sanitation + high pig density + high poverty"
)

custom_colors <- c(
  "all low" = "#D9D9D9",
  "A = poor sanitation" = "#FFD700",
  "B = high pig density" = "#9400D3",
  "AB = poor sanitation + high pig density" = "#00FA9A",
  "C = high poverty" = "#FF0000",
  "AC = poor sanitation + high poverty" = "#FF8C00",
  "BC = high pig density + high poverty" = "#FFB6C1",
  "ABC = poor sanitation + high pig density + high poverty" = "#5D4037",
  "lakes & lagoons" = "lightblue"
)

modal_fun <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

district_modes <- terra::extract(
  combo_rast,
  ghana_vect,
  fun = modal_fun,
  na.rm = TRUE
)

modal_report_df <- data.frame(
  Region   = ghana_vect$NAME_1,
  District = ghana_vect$NAME_2,
  Mode_ID  = district_modes[, 2]
) %>%
  mutate(
    Mode_ID = as.integer(Mode_ID),
    Average_Modal_Risk = risk_labels[as.character(Mode_ID)]
  )

modal_report_df$Average_Modal_Risk[is.na(modal_report_df$Average_Modal_Risk)] <- "all low"

ghana_modal_sf <- sf::st_as_sf(ghana_vect) %>%
  left_join(
    modal_report_df,
    by = c("NAME_1" = "Region", "NAME_2" = "District")
  )

ghana_modal_sf$Average_Modal_Risk <- factor(
  ghana_modal_sf$Average_Modal_Risk,
  levels = names(custom_colors)
)

lakes_sf_modal <- lakes_sf_correct
lakes_sf_modal$Average_Modal_Risk <- factor(
  "lakes & lagoons",
  levels = names(custom_colors)
)

district_modal_vect <- terra::vect(ghana_modal_sf)

modal_raster <- terra::rasterize(
  district_modal_vect,
  combo_rast,
  field = "Mode_ID",
  touches = TRUE
)

modal_raster <- terra::crop(modal_raster, ghana_vect, snap = "out")
modal_raster <- terra::mask(modal_raster, ghana_vect)

modal_raster_path <- file.path(
  out_dir,
  "Ghana_District_Modal_Risk_Base_2022_lakesblue_v2.tif"
)

terra::writeRaster(
  modal_raster,
  modal_raster_path,
  overwrite = TRUE
)

p_base <- ggplot() +
  geom_sf(
    data = ghana_modal_sf,
    aes(fill = Average_Modal_Risk),
    color = "black",
    linewidth = 0.15
  ) +
  geom_sf(
    data = lakes_sf_modal,
    aes(fill = Average_Modal_Risk),
    color = NA,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = custom_colors,
    breaks = names(custom_colors),
    limits = names(custom_colors),
    name = "District Modal\nRisk Factor",
    drop = FALSE
  ) +
  coord_sf(
    xlim = st_bbox(ghana_modal_sf)[c("xmin", "xmax")],
    ylim = st_bbox(ghana_modal_sf)[c("ymin", "ymax")],
    expand = FALSE
  ) +
  guides(
    fill = guide_legend(
      ncol = 2,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      override.aes = list(color = NA, alpha = 1)
    )
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 8.2),
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.x = unit(0.15, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.margin = margin(t = 4, r = 4, b = 4, l = 4),
    axis.title = element_blank(),
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_with_axes <- p_base

p_no_axes <- p_base +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 8.2),
    legend.key.width = unit(0.75, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.x = unit(0.15, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.margin = margin(t = 4, r = 4, b = 4, l = 4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

path_with_axes <- file.path(
  out_dir,
  "Ghana_Modal_Risk_Base_With_Axes_lakesblue_v2.png"
)

path_no_axes <- file.path(
  out_dir,
  "Ghana_Modal_Risk_Base_No_Axes_lakesblue_v2.png"
)

ggsave(path_with_axes, p_with_axes, width = 8, height = 10, dpi = 400, bg = "white")
ggsave(path_no_axes, p_no_axes, width = 8, height = 10, dpi = 400, bg = "white")

ft_modal <- flextable(
  modal_report_df %>%
    select(Region, District, Average_Modal_Risk)
) %>%
  theme_vanilla() %>%
  merge_v(j = ~ Region) %>%
  set_header_labels(
    Average_Modal_Risk = "Average (Modal) Risk"
  ) %>%
  bg(i = ~ Average_Modal_Risk == "all low",
     j = "Average_Modal_Risk", bg = "#D9D9D9") %>%
  bg(i = ~ Average_Modal_Risk == "A = poor sanitation",
     j = "Average_Modal_Risk", bg = "#FFD700") %>%
  bg(i = ~ Average_Modal_Risk == "B = high pig density",
     j = "Average_Modal_Risk", bg = "#9400D3") %>%
  bg(i = ~ Average_Modal_Risk == "AB = poor sanitation + high pig density",
     j = "Average_Modal_Risk", bg = "#00FA9A") %>%
  bg(i = ~ Average_Modal_Risk == "C = high poverty",
     j = "Average_Modal_Risk", bg = "#FF0000") %>%
  bg(i = ~ Average_Modal_Risk == "AC = poor sanitation + high poverty",
     j = "Average_Modal_Risk", bg = "#FF8C00") %>%
  bg(i = ~ Average_Modal_Risk == "BC = high pig density + high poverty",
     j = "Average_Modal_Risk", bg = "#FFB6C1") %>%
  bg(i = ~ Average_Modal_Risk == "ABC = poor sanitation + high pig density + high poverty",
     j = "Average_Modal_Risk", bg = "#5D4037") %>%
  color(
    i = ~ Average_Modal_Risk %in% c(
      "B = high pig density",
      "C = high poverty",
      "ABC = poor sanitation + high pig density + high poverty"
    ),
    j = "Average_Modal_Risk",
    color = "white"
  ) %>%
  autofit()

doc_modal_report <- read_docx() %>%
  body_add_par(
    "Ghana T. solium Mapping: District Modal Risk Report (Base Map)",
    style = "heading 1"
  ) %>%
  body_add_par(
    "Modal risk was calculated from the base combination map. Lakes and lagoons are excluded from land calculations and shown as a separate light-blue map layer.",
    style = "Normal"
  ) %>%
  body_add_par(
    paste("Report generated:", Sys.Date()),
    style = "Normal"
  ) %>%
  body_add_par("1. Modal Risk Visualization With Axes", style = "heading 2") %>%
  body_add_img(src = path_with_axes, width = 5.8, height = 7.2) %>%
  body_add_break() %>%
  body_add_par("2. Modal Risk Visualization Without Axes", style = "heading 2") %>%
  body_add_img(src = path_no_axes, width = 5.8, height = 7.2) %>%
  body_add_break() %>%
  body_add_par("3. District Modal Risk Classification Table", style = "heading 2") %>%
  body_add_flextable(ft_modal)

out_file <- file.path(
  out_dir,
  "Ghana_District_Modal_Risk_Base_Report_lakesblue_v2.docx"
)

print(doc_modal_report, target = out_file)

cat("\n--- MODAL BASE REPORT COMPLETE ---\n")
cat("Folder:", out_dir, "\n")
cat("Raster:", modal_raster_path, "\n")
cat("PNG with axes:", path_with_axes, "\n")
cat("PNG no axes:", path_no_axes, "\n")
cat("Word report:", out_file, "\n")







         # ==============================================================================
         # FULL CODE: GHANA BASE RISK + MODAL RISK MAPS — V6 FINAL PRODUCTION READY
         # ==============================================================================
         
         library(terra)
         library(sf)
         library(ggplot2)
         library(dplyr)
         library(ggnewscale)
         library(ggspatial)
         library(patchwork)
         library(grid)
         library(tibble)
         library(tidyterra)
         library(ragg)
         library(ggrepel)
         library(shadowtext)
         
         # Suppress spherical geometry warnings/errors by forcing planar calculations
         sf::sf_use_s2(FALSE)
         
         out_dir <- "outputs/Prevalence_On_Risk_Maps_V6"
         dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
         out_dir_donor <- file.path(out_dir, "Donor_Presentation")
         dir.create(out_dir_donor, showWarnings = FALSE, recursive = TRUE)
         
         # ==============================================================================
         # FONT SIZE CONSTANTS (SCALED UP FOR MAXIMUM READABILITY)
         # ==============================================================================
         
         SZ_MAP_TITLE   <- 18
         SZ_DIST_LABEL  <- 6.2  
         SZ_INSET_LABEL <- 5.5  
         
         # Legend Typography Scale-Up
         SZ_LEG_HEADING <- 8.0  
         SZ_LEG_SUBHEAD <- 6.8  
         SZ_LEG_ITEM    <- 6.2  
         SZ_LEG_MAG     <- 5.6  
         SZ_LEG_SOURCE  <- 6.2  
         
         # ==============================================================================
         # 1. PREVALENCE DATA (SYNCHRONIZED AND CORRECTED)
         # ==============================================================================
         
         prev_districts <- tribble(
           ~Disease, ~District, ~Village_Town, ~Latitude, ~Longitude, ~prevalence, ~prevalence_source, ~estimate_type,
           
           "PCC", "Kumasi",               "Kumasi (Abbatoir)",               6.660167,  -1.604278,  4.2,  "Informed", "Point estimate",
           "PCC", "Ablekuma West",        "(Shiabu, Chemunaa, Glefe)",       5.527320,  -0.254650,  1.8,  "Informed", "Point estimate",
           "PCC", "Weija Gbawe",          "(Gbawe, Pambros)",                5.556882,  -0.306150,  4.5,  "Informed", "Point estimate",
           "PCC", "Korle Klottey",        "Bola Beach",                      5.550300,  -0.186220, 11.6,  "Informed", "Point estimate",
           "PCC", "Ledzokuku",            "Teshie",                          5.582234,  -0.113510,  5.9,  "Informed", "Point estimate",
           "PCC", "Bawku West",           "Bawku / Zebilla",                10.939946,  -0.482999, 56.25, "Informed", "Median",
           "PCC", "Kasena Nankana East",  "Navrongo",                       10.892812,  -1.088373, 59.7,  "Informed", "Median",
           "PCC", "Bolgatanga",           "Bolgatanga",                     10.785396,  -0.848732, 39.7,  "Informed", "Median",
           "PCC", "Bongo",                "(Bongo, Bougosonge, Namoo)",     10.910450,  -0.825440, 21.4,  "Informed", "Point estimate",
           "PCC", "East Mamprusi",        "Nalerigu",                       10.525509,  -0.370179, 38.6,  "Informed", "Point estimate",
           "PCC", "Garu",                 "Garu",                           10.853225,  -0.175937, 33.4,  "Informed", "Point estimate",
           "PCC", "Kasena Nankana West",  "Chiana",                         10.855414,  -1.272788, 41.7,  "Informed", "Point estimate",
           "PCC", "Pusiga",               "Widana",                         11.085833,  -0.020833, 49.4,  "Informed", "Point estimate",
           "PCC", "Wa Municipal",         "Wa",                             10.062326,  -2.505427, 17.1,  "Informed", "Point estimate",
           "PCC", "West Gonja",           "Damango",                         9.085065,  -1.820076, 29.8,  "Informed", "Point estimate",
           "PCC", "Yendi",                "Yendi",                           9.446648,  -0.003507, 30.0,  "Informed", "Point estimate",
           "PCC", "Bunkpurugu & Yunyoo",  "Bunkpurugu & Yunyoo",             10.547874,  -0.029556, 57.8,  "Informed", "Point estimate",
           
           "HCC", "Bolgatanga",           "Bolgatanga",                     10.785396,  -0.848732,  6.5,  "Informed", "Point estimate",
           "HCC", "Bawku West",           "Zebilla",                        10.939946,  -0.482999,  3.1,  "Informed", "Point estimate",
           "HCC", "Kasena Nankana East",  "Navrongo",                       10.892812,  -1.088373,  4.3,  "Informed", "Point estimate",
           "HCC", "Ablekuma West",        "(Shiabu, Chemunaa, Glefe)",       5.527320,  -0.254650,  8.4,  "Informed", "Point estimate",
           "HCC", "Accra",                "(James Town, Korle Gonno)",       5.532410,  -0.220660,  3.5,  "Informed", "Point estimate",
           "HCC", "Ayawaso Central",      "Pig Farm",                        5.600200,  -0.201620,  2.4,  "Informed", "Point estimate",
           "HCC", "Korle Klottey",        "(Osu-Kajaanor, Bola Beach)",      5.550300,  -0.186220,  3.2,  "Informed", "Point estimate",
           "HCC", "Ledzokuku",            "Teshie",                          5.582234,  -0.113510,  8.1,  "Informed", "Point estimate",
           "HCC", "Weija Gbawe",          "(Gbawe, Pambros, Mallam)",        5.556882,  -0.306150,  3.9,  "Informed", "Point estimate",
           
           "HTT", "Kintampo North",  "(Gulumpe, Jato Akura, Chiranda)",      8.192963,  -1.677778,  4.6,  "Informed", "Median",
           "HTT", "Accra",           "Accra",                                5.550000,  -0.200000,  0.7,  "Observed", "Point estimate",
           "HTT", "Kumasi",          "KATH",                                 6.697479,  -1.631690,  0.4,  "Observed", "Median",
           "HTT", "Bunkpurugu & Yunyoo",  "Bunkpurugu & Yunyoo",             10.547874,  -0.029556, 44.6,  "Informed", "Point estimate",
           "HTT", "Kintampo South",  "Kintampo North & South",               7.981857,  -1.759460,  1.5,  "Observed", "Point estimate",
           
           "HCC-HTT", "Kintampo North", "Kintampo",                          8.053714,  -1.728919,  4.6,  "Observed", "Point estimate",
           "NCC",     "Wenchi",          "Wenchi",                           7.738011,  -2.103519,   NA,  "Observed", "Cases only"
         )
         
         # ==============================================================================
         # 2. MAP OBJECTS
         # ==============================================================================
         
         ghana_sf <- st_as_sf(ghana_vect)
         
         # ==============================================================================
         # 3. DISTRICT POLYGON BOUNDARY JOIN
         # ==============================================================================
         
         district_polygons <- ghana_sf %>%
           mutate(
             District_match = case_when(
               NAME_2 == "Bunkpurugu Nakpanduri" ~ "Bunkpurugu & Yunyoo",
               TRUE                              ~ NAME_2
             )
           ) %>%
           group_by(District_match) %>%
           summarise(geometry = st_union(geometry), .groups = "drop")
         
         district_inside_points <- district_polygons %>%
           st_point_on_surface() %>%
           mutate(
             inside_lon = st_coordinates(.)[, 1],
             inside_lat = st_coordinates(.)[, 2]
           ) %>%
           st_drop_geometry() %>%
           select(District_match, inside_lon, inside_lat)
         
         # ==============================================================================
         # 4. CLEAN DATA
         # ==============================================================================
         
         prev_clean0 <- prev_districts %>%
           mutate(
             District_match    = District,
             Disease           = factor(Disease, levels = c("PCC","HCC","HTT","NCC","HCC-HTT")),
             prevalence_source = factor(prevalence_source, levels = c("Observed","Informed")),
             estimate_type     = factor(estimate_type, levels = c("Point estimate","Median","Cases only")),
             prevalence_group  = case_when(
               Disease == "PCC" & prevalence <  10                       ~ "PCC: 0-9.99",
               Disease == "PCC" & prevalence >= 10 & prevalence < 25    ~ "PCC: 10-24.99",
               Disease == "PCC" & prevalence >= 25 & prevalence < 40    ~ "PCC: 25-39.99",
               Disease == "PCC" & prevalence >= 40                       ~ "PCC: >=40",
               Disease == "HCC" & prevalence <  4                       ~ "HCC: 0-3.99",
               Disease == "HCC" & prevalence >= 4  & prevalence < 8     ~ "HCC: 4-7.99",
               Disease == "HCC" & prevalence >= 8                       ~ "HCC: >=8",
               Disease == "HTT" & prevalence <= 0.7                     ~ "HTT: 0-0.7",
               Disease == "HTT" & prevalence >= 1.5 & prevalence <= 5.5 ~ "HTT: 1.5-5.5",
               Disease == "HTT" & prevalence >= 44.6                    ~ "HTT: >=44.6",
               Disease == "NCC"                                         ~ "NCC: cases only",
               Disease == "HCC-HTT"                                     ~ "HCC-HTT",
               TRUE ~ NA_character_
             ),
             point_size = case_when(
               prevalence_group == "PCC: 0-9.99"   ~ 3.0,
               prevalence_group == "PCC: 10-24.99" ~ 4.4,
               prevalence_group == "PCC: 25-39.99" ~ 5.8,
               prevalence_group == "PCC: >=40"     ~ 7.2,
               prevalence_group == "HCC: 0-3.99"   ~ 3.0,
               prevalence_group == "HCC: 4-7.99"   ~ 4.4,
               prevalence_group == "HCC: >=8"      ~ 5.8, 
               prevalence_group == "HTT: 0-0.7"    ~ 3.0,
               prevalence_group == "HTT: 1.5-5.5"  ~ 4.4,
               prevalence_group == "HTT: >=44.6"   ~ 7.2,
               Disease %in% c("NCC","HCC-HTT")     ~ 5.0, 
               TRUE ~ 4.4
             )
           ) %>%
           left_join(district_inside_points, by = "District_match") %>%
           mutate(
             nudge_strength = 0.42,
             Longitude_base = ifelse(!is.na(inside_lon), Longitude + nudge_strength*(inside_lon - Longitude), Longitude),
             Latitude_base  = ifelse(!is.na(inside_lat), Latitude  + nudge_strength*(inside_lat - Latitude), Latitude)
           )
         
         # ==============================================================================
         # 5. GEOMETRICALLY CONSTRAINED ENGINE (FORCE WEIJA AND KORLE KLOTTEY SNAP)
         # ==============================================================================
         
         safe_offset_one_district <- function(dat, polygon_row, name) {
           dat <- dat %>% arrange(Disease, prevalence_group)
           n   <- nrow(dat)
           if (n == 1) {
             dat$Longitude_plot <- dat$Longitude_base
             dat$Latitude_plot  <- dat$Latitude_base
             return(dat)
           }
           
           angles      <- 2 * pi * (seq_len(n) - 1) / n
           base_radius <- 0.042 
           
           # Flags to bypass correction for specific districts
           skip_correction <- FALSE 
           
           # 1. Apply standard cluster adjustments
           if (name %in% c("Korle Klottey", "Ledzokuku", "Ablekuma West", "Accra")) {
             angles <- pi * (seq_len(n) - 1) / (n - 1) + 0.15 
           }
           
           # 2. HARD-CODED DE-CONFLICTION FOR WEIJA GBAWE
           if (name == "Weija Gbawe") {
             pcc_idx <- which(dat$Disease == "PCC")
             hcc_idx <- which(dat$Disease == "HCC")
             if(length(pcc_idx) > 0) angles[pcc_idx] <- 0.5 
             if(length(hcc_idx) > 0) angles[hcc_idx] <- 2.5 
             base_radius <- 0.08 
           }
           
           candidate   <- dat
           candidate$Longitude_plot <- dat$Longitude_base + base_radius * cos(angles)
           candidate$Latitude_plot  <- dat$Latitude_base  + base_radius * sin(angles)
           
           # 3. OVERRIDE: FORCE KORLE KLOTTEY PCC INSIDE POLYGON
           if (name == "Korle Klottey") {
             pcc_idx <- which(candidate$Disease == "PCC")
             if(length(pcc_idx) > 0) {
               candidate$Longitude_plot[pcc_idx] <- -0.19 
               candidate$Latitude_plot[pcc_idx]  <- 5.56
               skip_correction <- TRUE # Force-bypass the nearest-neighbor loop
             }
           }
           
           # Geometric Constraint Logic
           candidate_sf <- st_as_sf(candidate, coords = c("Longitude_plot","Latitude_plot"), crs = st_crs(ghana_sf), remove = FALSE)
           inside       <- as.logical(st_within(candidate_sf, polygon_row, sparse = FALSE)[,1])
           
           # Only run correction if skip_correction is FALSE
           if (!skip_correction) {
             shrink_steps <- c(0.034, 0.026, 0.018, 0.012, 0.006, 0.000)
             for (i in which(!inside)) {
               for (r in shrink_steps) {
                 test_lon <- candidate$Longitude_plot[i] # Uses current calculated pos
                 test_lat <- candidate$Latitude_plot[i]
                 test_sf  <- st_as_sf(data.frame(x = test_lon, y = test_lat), coords = c("x","y"), crs = st_crs(ghana_sf))
                 
                 if (as.logical(st_within(test_sf, polygon_row, sparse = FALSE)[1,1])) {
                   candidate$Longitude_plot[i] <- test_lon
                   candidate$Latitude_plot[i]  <- test_lat
                   break
                 }
                 if (r == 0.000) {
                   nearest_pt <- st_nearest_points(test_sf, polygon_row)
                   coords_m   <- st_coordinates(st_cast(nearest_pt, "POINT"))
                   candidate$Longitude_plot[i] <- coords_m[2, 1]
                   candidate$Latitude_plot[i]  <- coords_m[2, 2]
                 }
               }
             }
           }
           return(candidate)
         }
         # ==============================================================================
         # 6. SYMBOL AND COLOUR SETTINGS
         # ==============================================================================
         
         disease_cols <- c("PCC" = "#FFFF00", "HCC" = "#0057FF", "HTT" = "#FFFF00", "NCC" = "#FFFF00", "HCC-HTT" = "#FFFF00")
         source_cols  <- c("Observed" = "red", "Informed" = "black")
         risk_labels  <- c("0" = "all low", "1" = "A = poor sanitation", "2" = "B = high pig density", "3" = "AB", "4" = "C = high poverty", "5" = "AC", "6" = "BC", "7" = "ABC")
         
         custom_colors <- c(
           "all low"              = "#D9D9D9",
           "A = poor sanitation"  = "#FFD700",
           "B = high pig density" = "#9400D3",
           "AB"                   = "#00FA9A",
           "C = high poverty"     = "#FF0000",
           "AC"                   = "#FF8C00",
           "BC"                   = "#FFB6C1",
           "ABC"                  = "#5D4037",
           "lakes & lagoons"      = "lightblue"
         )
         
         # ==============================================================================
         # 7. ACCRA INSET DISTRICTS
         # ==============================================================================
         
         accra_inset_districts <- c("Ablekuma West", "Weija Gbawe", "Korle Klottey", "Ledzokuku", "Accra", "Ayawaso Central")
         
         # ==============================================================================
         # 8. MAIN MAP CANVAS BOUNDARIES & LABEL ANCHORS (MODIFIED FOR SPACING & ARROWS)
         # ==============================================================================
         
         coord_xmin <- -5.50
         coord_xmax <-  2.80
         coord_ymin <-  4.30
         coord_ymax <- 11.85
         
         # MODIFIED: text_y coordinates updated to create specific spacing on the right
         district_label_anchors <- tribble(
           ~District,              ~text_x,  ~text_y,  ~hjust, ~vjust,
           # Top Row (Keep as is)
           "Kasena Nankana East",  -2.25,    11.42,    0.5,    0.0,  
           "Bongo",                -1.15,    11.42,    0.5,    0.0,  
           "Bolgatanga",           -0.05,    11.42,    0.5,    0.0,  
           "Pusiga",                0.85,    11.42,    0.5,    0.0,
           
           # Left Column (Keep as is)
           "Kasena Nankana West",  -3.30,    11.15,    1.0,    0.5,
           "Wa Municipal",         -3.30,     9.80,    1.0,    0.5,
           "West Gonja",           -3.30,     8.90,    1.0,    0.5,
           "Kintampo North",       -3.30,     8.10,    1.0,    0.5,
           "Kintampo South",       -3.30,     7.50,    1.0,    0.5,
           "Wenchi",               -3.30,     6.90,    1.0,    0.5,
           "Kumasi",               -3.30,     6.20,    1.0,    0.5,
           
           # MODIFIED: Adjusted text_y to increase vertical gaps on the right
           "Bawku West",            1.15,    10.80,    0.0,    0.5,
           "Garu",                  1.15,    10.10,    0.0,    0.5, # Pushed down
           "Bunkpurugu & Yunyoo",   1.15,     9.40,    0.0,    0.5, # Pushed down
           "East Mamprusi",         1.15,     8.70,    0.0,    0.5, # Pushed down
           "Yendi",                 1.15,     8.00,    0.0,    0.5  # Pushed down
         )
         
         # Label-aware engine: Automatically calculates arrow tips based on new anchor positions
         main_arrows_df <- prev_clean %>%
           filter(!District %in% accra_inset_districts) %>%
           left_join(district_label_anchors, by = "District") %>%
           filter(!is.na(text_x)) %>%
           mutate(
             # Dynamic padding: Adjusts based on anchor alignment (hjust/vjust)
             x_pad = (hjust - 0.5) * 0.4, 
             y_pad = (vjust - 0.5) * 0.2,
             
             # Target boundary
             target_x = text_x + x_pad,
             target_y = text_y + y_pad,
             
             # Distance to target
             dx = target_x - Longitude_plot,
             dy = target_y - Latitude_plot,
             dist = sqrt(dx^2 + dy^2),
             
             # Final arrow tip: Stops 0.15 units before reaching the text
             offset_dist = pmax(0.15, dist - 0.15),
             
             arrow_xend = Longitude_plot + (dx / dist) * offset_dist,
             arrow_yend = Latitude_plot + (dy / dist) * offset_dist
           )
         # ==============================================================================
         # 9. ACCRA REGIONAL INSET SETUP (CORRECTED ROBUST FILTER)
         # ==============================================================================
         
         accra_xlim <- c(-0.35, -0.05) 
         accra_ylim <- c(5.45,  5.65)  
         
         accra_rect_sf      <- st_as_sf(st_as_sfc(st_bbox(c(xmin = accra_xlim[1], xmax = accra_xlim[2], ymin = accra_ylim[1], ymax = accra_ylim[2]), crs = 4326)))
         accra_districts_sf <- district_polygons %>% st_filter(accra_rect_sf, .predicate = st_intersects)
         
         # FIXED: Filter by District list to prevent nudged points from disappearing
         prev_sf_accra      <- prev_sf %>% filter(District %in% accra_inset_districts)
         
         accra_corner_dots  <- data.frame(lon = c(accra_xlim[2], accra_xlim[2]), lat = c(accra_ylim[1], accra_ylim[2]))
         
         accra_label_anchors <- tribble(
           ~District,         ~text_x,  ~text_y,  ~hjust, ~vjust,
           "Weija Gbawe",     -0.30,    5.50,     0.5,    0.5,
           "Ablekuma West",   -0.24,    5.48,     0.5,    0.5,
           "Accra",           -0.18,    5.52,     0.5,    0.5,
           "Ayawaso Central", -0.12,    5.62,     0.5,    0.5,
           "Korle Klottey",   -0.12,    5.50,     0.5,    0.5,
           "Ledzokuku",       -0.08,    5.58,     0.5,    0.5
         )
         
         accra_arrows_df <- prev_sf_accra %>%
           st_drop_geometry() %>%
           left_join(accra_label_anchors, by = "District") %>%
           filter(!is.na(text_x)) %>%
           mutate(
             dx = text_x - Longitude_plot,
             dy = text_y - Latitude_plot,
             dist = sqrt(dx^2 + dy^2),
             offset_dist = 0.005,
             arrow_xend = text_x - (dx / dist) * offset_dist,
             arrow_yend = text_y - (dy / dist) * offset_dist
           )
         
         
        # ==============================================================================
         # 10. MAIN MAP FUNCTION
         # ==============================================================================
         
         make_main_map <- function(p_risk, map_title = "Ghana") {
           p_risk +
             ggnewscale::new_scale_fill() +
             geom_sf(data = prev_sf %>% filter(Disease %in% c("PCC", "HCC")), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 21, stroke = 0.4, alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(Disease == "HTT"), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 22, stroke = 0.4, alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(Disease == "HCC-HTT"), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 23, stroke = 0.4, alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(Disease == "NCC"), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 8, stroke = 0.4, alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(estimate_type == "Point estimate"), shape = 16, size = 1.0, colour = "black", show.legend = FALSE) +
             geom_segment(data = main_arrows_df, aes(x = Longitude_plot, y = Latitude_plot, xend = arrow_xend, yend = arrow_yend), arrow = arrow(length = unit(0.18, "cm"), type = "closed"), linewidth = 0.75, colour = "black", alpha = 0.95) +
             geom_text(data = district_label_anchors, aes(x = text_x, y = text_y, label = District), size = SZ_DIST_LABEL, fontface = "bold", colour = "black") +
             scale_fill_manual(values = disease_cols) +
             scale_colour_manual(values = source_cols) +
             scale_size_identity() +
             coord_sf(xlim = c(coord_xmin, coord_xmax), ylim = c(coord_ymin, coord_ymax), expand = FALSE) +
             labs(title = map_title) +
             theme_void() +
             theme(plot.title = element_text(size = SZ_MAP_TITLE, face = "bold", hjust = 0.5), plot.margin = margin(6, 6, 6, 6))
         }
         
         # ==============================================================================
         # 11. REGIONAL INSET MAP FUNCTION (CLEAN TEXT STYLE)
         # ==============================================================================
         
         make_inset_map <- function(p_risk) {
           p_risk +
             ggnewscale::new_scale_fill() +
             geom_sf(data = accra_districts_sf, fill = NA, colour = "grey30", linewidth = 0.8, inherit.aes = FALSE) +
             # Consistent symbol rendering
             geom_sf(data = prev_sf_accra %>% filter(Disease %in% c("PCC", "HCC")), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 21, stroke = 0.5, alpha = 0.9, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(Disease == "HTT"), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 22, stroke = 0.5, alpha = 0.9, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(Disease == "HCC-HTT"), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 23, stroke = 0.5, alpha = 0.9, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(Disease == "NCC"), aes(fill = Disease, colour = prevalence_source, size = point_size), shape = 8, stroke = 0.5, alpha = 0.9, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(estimate_type == "Point estimate"), shape = 16, size = 1.1, colour = "black", show.legend = FALSE) +
             # Retaining your arrow infrastructure
             geom_segment(data = accra_arrows_df, aes(x = Longitude_plot, y = Latitude_plot, xend = arrow_xend, yend = arrow_yend), 
                          arrow = arrow(length = unit(0.08, "cm"), type = "closed"), linewidth = 0.4, colour = "black", inherit.aes = FALSE) +
             # UPDATED LABELING: Standard bold, centered, black font (No Shadow)
             geom_text(data = accra_label_anchors, aes(x = text_x, y = text_y, label = District), 
                       size = 4.2, fontface = "bold", colour = "black", hjust = 0.5, vjust = 0.5, inherit.aes = FALSE) +
             scale_fill_manual(values = disease_cols) +
             scale_colour_manual(values = source_cols) +
             scale_size_identity() +
             coord_sf(xlim = c(-0.37, -0.05), ylim = c(5.45, 5.65), expand = FALSE) +
             labs(title = "Greater Accra") + 
             theme_void() +
             theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b=5)),
                   panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
                   plot.background = element_rect(fill = "white", colour = NA))
         }
         # ==============================================================================
         # 12. LEGEND PANEL
         # ==============================================================================
         
         risk_legend_df <- tribble(
           ~x, ~y, ~label, ~fill,
           1, 28.0, "all low",              "#D9D9D9",
           1, 27.0, "A = poor sanitation",  "#FFD700",
           1, 26.0, "B = high pig density", "#9400D3",
           1, 25.0, "C = high poverty",     "#FF0000",
           1, 24.0, "AB",                   "#00FA9A",
           1, 23.0, "AC",                   "#FF8C00",
           1, 22.0, "BC",                   "#FFB6C1",
           1, 21.0, "ABC",                  "#5D4037",
           1, 20.0, "lakes & lagoons",      "lightblue"
         )
         disease_legend_df <- tribble(
           ~x, ~y, ~Disease, ~shape, ~fill,
           1, 17.6, "PCC",     21, "#FFFF00",
           1, 16.6, "HCC",     21, "#0057FF",
           1, 15.6, "HTT",     22, "#FFFF00",
           1, 14.6, "NCC",      8, "#FFFF00",
           1, 13.6, "HCC-HTT", 23, "#FFFF00"
         )
         magnitude_legend_df <- tribble(
           ~x, ~y,   ~heading, ~label,       ~shape, ~fill,      ~size,
           1, 10.8,  "PCC",    "0-9.99",      21,    "#FFFF00",   3.0,
           1,  9.9,  "PCC",    "10-24.99",    21,    "#FFFF00",   4.4,
           1,  9.0,  "PCC",    "25-39.99",    21,    "#FFFF00",   5.8,
           1,  8.1,  "PCC",    ">=40",        21,    "#FFFF00",   7.2,
           1,  6.7,  "HCC",    "0-3.99",      21,    "#0057FF",   3.0,
           1,  5.8,  "HCC",    "4-7.99",      21,    "#0057FF",   4.4,
           1,  4.9,  "HCC",    ">=8",         21,    "#0057FF",   5.8, 
           1,  3.5,  "HTT",    "0-0.7",       22,    "#FFFF00",   3.0,
           1,  2.6,  "HTT",    "1.5-5.5",     22,    "#FFFF00",   4.4,
           1,  1.7,  "HTT",    ">=44.6",      22,    "#FFFF00",   7.2
         )
         source_legend_df <- tribble(
           ~x, ~y,   ~label,     ~outline,
           1, -0.2,  "Informed", "black",
           1, -1.1,  "Observed", "red"
         )
         type_legend_df <- tribble(
           ~x, ~y,   ~label,            ~dot,
           1, -3.2,  "Point estimate",  TRUE,
           1, -4.1,  "Median of range", FALSE
         )
         
         p_right_legend <- ggplot() +
           annotate("text", x = 0.4, y = 29.0, label = "Risk factors", hjust = 0, fontface = "bold", size = SZ_LEG_HEADING) +
           geom_tile(data = risk_legend_df, aes(x = x, y = y, fill = fill), width = 0.38, height = 0.65, colour = "grey35") +
           geom_text(data = risk_legend_df, aes(x = 1.35, y = y, label = label), hjust = 0, size = SZ_LEG_ITEM) +
           
           annotate("text", x = 0.4, y = 18.7, label = "Disease", hjust = 0, fontface = "bold", size = SZ_LEG_HEADING) +
           geom_point(data = disease_legend_df, aes(x = x, y = y, shape = factor(shape), fill = fill), size = 6.0, colour = "black", stroke = 0.5) +
           geom_text(data = disease_legend_df, aes(x = 1.35, y = y, label = Disease), hjust = 0, size = SZ_LEG_ITEM) +
           
           annotate("text", x = 0.4,  y = 11.9,  label = "Magnitude", hjust = 0, fontface = "bold", size = SZ_LEG_HEADING) +
           annotate("text", x = 0.48, y = 11.25, label = "PCC", hjust = 0, fontface = "bold", size = SZ_LEG_SUBHEAD) +
           annotate("text", x = 0.48, y =  7.15, label = "HCC", hjust = 0, fontface = "bold", size = SZ_LEG_SUBHEAD) +
           annotate("text", x = 0.48, y =  3.95, label = "HTT", hjust = 0, fontface = "bold", size = SZ_LEG_SUBHEAD) +
           geom_point(data = magnitude_legend_df, aes(x = x, y = y, shape = factor(shape), fill = fill, size = size * 1.1), colour = "black", stroke = 0.45) +
           geom_text(data = magnitude_legend_df, aes(x = 1.45, y = y, label = label), hjust = 0, size = SZ_LEG_MAG) +
           
           annotate("text", x = 0.4, y = 0.75, label = "Prevalence source", hjust = 0, fontface = "bold", size = SZ_LEG_HEADING) +
           geom_point(data = source_legend_df, aes(x = x, y = y, colour = outline), shape = 21, size = 6.0, fill = "#FFFF00", stroke = 0.5) +
           geom_text(data = source_legend_df, aes(x = 1.35, y = y, label = label), hjust = 0, size = SZ_LEG_SOURCE) +
           
           annotate("text", x = 0.4, y = -2.25, label = "Prevalence type", hjust = 0, fontface = "bold", size = SZ_LEG_HEADING) +
           geom_point(data = type_legend_df, aes(x = x, y = y), shape = 21, size = 6.0, fill = "#FFFF00", colour = "black", stroke = 0.5) +
           geom_point(data = type_legend_df %>% filter(dot), aes(x = x, y = y), shape = 16, size = 1.5, colour = "black") +
           geom_text(data = type_legend_df, aes(x = 1.35, y = y, label = label), hjust = 0, size = SZ_LEG_SOURCE) +
           
           scale_fill_identity() +
           scale_colour_identity() +
           scale_shape_manual(values = c("8"=8,"21"=21,"22"=22,"23"=23)) +
           scale_size_identity() +
           coord_cartesian(xlim = c(0.3, 4.8), ylim = c(-4.8, 29.6), expand = FALSE) +
           theme_void() +
           theme(
             plot.background = element_rect(fill = "white", colour = NA),
             plot.margin      = margin(8, 2, 8, 4),
             legend.position  = "none"
           )
         
         # ==============================================================================
         # 13. LAYOUT ASSEMBLY (MAIN MAP EXPANDED)
         # ==============================================================================
         
         build_layout <- function(p_main, p_inset) {
           p_inset_box <- p_inset + theme(plot.margin = margin(5, 5, 5, 5))
           (p_main | (plot_spacer() / p_inset_box) | p_right_legend) + 
             plot_layout(widths = c(6.5, 2, 1.5))
         }
         
         # ==============================================================================
         # 14. CONNECTOR (CENTROID TARGET TRACKING ASSEMBLY)
         # ==============================================================================
         
         pw_main_w  <- 6.5
         pw_mid_w   <- 1.8
         pw_leg_w   <- 1.7
         pw_total_w <- pw_main_w + pw_mid_w + pw_leg_w
         
         main_frac  <- pw_main_w / pw_total_w
         inset_x_lo <- main_frac - 0.015  
         
         map_x_range <- coord_xmax - coord_xmin
         map_y_range <- coord_ymax - coord_ymin
         title_frac  <- 0.050
         bot_frac    <- 0.020
         plot_h_frac <- 1.0 - title_frac - bot_frac
         
         rect_cx_frac <- (mean(accra_xlim) - coord_xmin) / map_x_range
         rect_cy_frac <- (mean(accra_ylim) - coord_ymin) / map_y_range
         
         dot1_x <- main_frac * (0.012 + rect_cx_frac * (1 - 2*0.012))
         dot1_y <- bot_frac  + rect_cy_frac * plot_h_frac
         
         accra_frac_from_bot <- (mean(accra_ylim) - coord_ymin) / map_y_range
         spacer_h_u  <- round((1 - accra_frac_from_bot) * 12)
         inset_h_u   <- max(round(accra_frac_from_bot * 18), 1)
         total_h_u   <- spacer_h_u + inset_h_u
         inset_frac  <- inset_h_u / total_h_u
         inset_y_ctr <- bot_frac + (inset_frac * 0.40) * plot_h_frac 
         
         draw_connector_npc <- function() {
           grid::grid.lines(
             x  = unit(c(dot1_x, inset_x_lo), "npc"),
             y  = unit(c(dot1_y, inset_y_ctr), "npc"),
             gp = grid::gpar(col = "red", lty = "dashed", lwd = 1.8)
           )
           grid::grid.points(
             x    = unit(dot1_x, "npc"), y    = unit(dot1_y, "npc"),
             pch  = 19, size = unit(0.008, "npc"), gp = grid::gpar(col = "red", fill = "red")
           )
           grid::grid.points(
             x    = unit(inset_x_lo,  "npc"), y    = unit(inset_y_ctr, "npc"),
             pch  = 19, size = unit(0.008, "npc"), gp = grid::gpar(col = "red", fill = "red")
           )
         }
         
         # ==============================================================================
         # 15. SAVE ENGINE CHANNELS
         # ==============================================================================
         
         save_png <- function(obj, path, width, height, dpi, bg = "white") {
           ragg::agg_png(path, width = width, height = height, units = "in", res = dpi, background = bg)
           print(obj); draw_connector_npc(); dev.off()
           invisible(path)
         }
         
         save_tiff <- function(obj, path, width, height, dpi, bg = "white") {
           ragg::agg_tiff(path, width = width, height = height, units = "in", res = dpi, background = bg, compression = "lzw")
           print(obj); draw_connector_npc(); dev.off()
           invisible(path)
         }
         
         save_pdf <- function(obj, path, width, height, bg = "white") {
           cairo_pdf(path, width = width, height = height, bg = bg)
           print(obj); draw_connector_npc(); dev.off()
           invisible(path)
         }
         
         # ==============================================================================
         # 16. BASE RISK MAP GENERATION
         # ==============================================================================
         
         combo_rast <- scenario_combos$base
         combo_rast <- crop_mask_to_ghana(combo_rast)
         combo_rast <- terra::mask(combo_rast, ghana_vect)
         combo_f    <- as.factor(combo_rast)
         levels(combo_f) <- data.frame(value = 0:7, risk_label = unname(risk_labels[as.character(0:7)]))
         
         lakes_sf_base <- lakes_sf_correct %>% mutate(risk_label = factor("lakes & lagoons", levels = names(custom_colors)))
         
         p_base_risk <- ggplot() +
           geom_spatraster(data = combo_f, aes(fill = risk_label), show.legend = FALSE) +
           geom_sf(data = lakes_sf_base, aes(fill = risk_label), color = NA, inherit.aes = FALSE, show.legend = FALSE) +
           geom_sf(data = ghana_sf, fill = NA, colour = "black", linewidth = 0.40) +
           scale_fill_manual(values = custom_colors, breaks = names(custom_colors), limits = names(custom_colors), drop = FALSE, na.value = "white")
         
         final_base <- build_layout(make_main_map(p_base_risk, "Ghana"), make_inset_map(p_base_risk))
         gc()
         
         # ==============================================================================
         # 17. MODAL RISK MAP GENERATION
         # ==============================================================================
         
         modal_fun <- function(x, na.rm = TRUE) {
           if (na.rm) x <- x[!is.na(x)]
           if (length(x) == 0) return(NA)
           ux <- unique(x)
           ux[which.max(tabulate(match(x, ux)))]
         }
         
         district_modes  <- terra::extract(combo_rast, ghana_vect, fun = modal_fun, na.rm = TRUE)
         modal_report_df <- data.frame(Region = ghana_vect$NAME_1, District = ghana_vect$NAME_2, Mode_ID  = district_modes[, 2]) %>%
           mutate(Mode_ID = as.integer(Mode_ID), Average_Modal_Risk = risk_labels[as.character(Mode_ID)])
         modal_report_df$Average_Modal_Risk[is.na(modal_report_df$Average_Modal_Risk)] <- "all low"
         
         ghana_modal_sf <- st_as_sf(ghana_vect) %>% left_join(modal_report_df, by = c("NAME_1" = "Region", "NAME_2" = "District"))
         ghana_modal_sf$Average_Modal_Risk <- factor(ghana_modal_sf$Average_Modal_Risk, levels = names(custom_colors))
         
         lakes_sf_modal <- lakes_sf_correct %>% mutate(Average_Modal_Risk = factor("lakes & lagoons", levels = names(custom_colors)))
         
         p_modal_risk <- ggplot() +
           geom_sf(data = ghana_modal_sf, aes(fill = Average_Modal_Risk), color = "black", linewidth = 0.12, show.legend = FALSE) +
           geom_sf(data = lakes_sf_modal, aes(fill = Average_Modal_Risk), color = NA, inherit.aes = FALSE, show.legend = FALSE) +
           scale_fill_manual(values = custom_colors, breaks = names(custom_colors), limits = names(custom_colors), drop = FALSE)
         
         final_modal <- build_layout(make_main_map(p_modal_risk, "Ghana"), make_inset_map(p_modal_risk))
         gc()
         
         # ==============================================================================
         # 18. EXPORT CHANNELS EXECUTION
         # ==============================================================================
         
         map_exports <- list(list(name = "Base", obj = final_base), list(name = "Modal", obj = final_modal))
         
         for (m in map_exports) {
           nm  <- m$name; obj <- m$obj
           
           save_png(obj, file.path(out_dir, paste0("Ghana_", nm, "_RiskMap_V12345667.png")), width = 22, height = 14, dpi = 900)
           save_pdf(obj, file.path(out_dir, paste0("Ghana_", nm, "_RiskMap_V12345667.pdf")), width = 22, height = 14)
           save_png(obj, file.path(out_dir_donor, paste0("Ghana_", nm, "_Donor_HighRes_V12345667.png")), width = 30, height = 18, dpi = 600)
           save_pdf(obj, file.path(out_dir_donor, paste0("Ghana_", nm, "_Donor_Print_V12345667.pdf")), width = 30, height = 18)
           save_tiff(obj, file.path(out_dir_donor, paste0("Ghana_", nm, "_Donor_Archival_V12345667.tiff")), width = 30, height = 18, dpi = 400)
           gc()
         }
         
         cat("Done. V6 outputs successfully compiled with scale modifications active.\n")
         
         
        
         
          # ==============================================================================
         # GHANA T. SOLIUM STUDY DISTRICT + PREVALENCE MAP
         # ==============================================================================
         
         library(sf)
         library(ggplot2)
         library(dplyr)
         library(tibble)
         library(ggnewscale)
         library(patchwork)
         library(grid)
         library(ragg)
         
         sf::sf_use_s2(FALSE)
         
         out_dir <- "outputs/Study_District_Prevalence_Maps"
         dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
         
         SZ_MAP_TITLE    <- 18
         SZ_MAIN_LABEL   <- 4.2
         SZ_INSET_LABEL  <- 3.1
         SZ_LEG_HEADING  <- 5.0
         SZ_LEG_ITEM     <- 4.5
         SZ_LEG_MAG      <- 4.2
         
         # ==============================================================================
         # 1. PREVALENCE DATA
         # ==============================================================================
         
         prev_districts <- tribble(
           ~Disease, ~District, ~Village_Town, ~Latitude, ~Longitude, ~prevalence, ~data_type,
           
           "PCC", "Bolgatanga",            "Bolgatanga", 10.785396, -0.848732, 51.8,  "Range (13.6-90)",
           "PCC", "Kasena Nankana East",   "Navrongo",   10.892812, -1.088373, 43.05, "Range (19.4-66.7)",
           "PCC", "Bawku West",            "Zebilla",    10.939946, -0.482999, 77.8,  "Point",
           "PCC", "Kumasi Metropolitan",   "Kumasi",      6.660167, -1.604278,  2.3,  "Point",
           "PCC", "Wa",                    "Wa",         10.062326, -2.505427,  8.3,  "Point",
           "PCC", "West Gonja",            "Damango",     9.085065, -1.820076, 10.8,  "Point",
           "PCC", "Yendi",                 "Yendi",       9.446648, -0.003507, 11.0,  "Point",
           "PCC", "East Mamprusi",         "Nalerigu",   10.525509, -0.370179, 12.1,  "Point",
           "PCC", "Bongo",                 "Bongo",      10.909080, -0.807506, 10.7,  "Point",
           "PCC", "Bongo",                 "Bougosonge", 10.850000, -0.816670, 15.4,  "Point",
           "PCC", "Kasena Nankana West",   "Chiana",     10.855414, -1.272788, 14.0,  "Point",
           "PCC", "Garu",                  "Garu",       10.853225, -0.175937, 10.5,  "Point",
           "PCC", "Bongo",                 "Namoo",      10.972280, -0.852142,  2.3,  "Point",
           "PCC", "Pusiga",                "Widana",     11.085833, -0.020833, 16.1,  "Point",
           "PCC", "Bunkpurugu Nakpanduri", "Bunkpurugu", 10.547874, -0.029556, 19.0,  "Point",
           
           "PCC", "Ledzokuku",             "Teshie",                                      NA, NA, 0.0, "Point",
           "PCC", "Weija Gbawe",           "Gbawe; Pambros",                              NA, NA, 0.0, "Point",
           "PCC", "Ablekuma West",         "Glefe; Shaibu; Chemunaa; Bola Beach/Labadi",   NA, NA, 0.0, "Point",
           "PCC", "Korle Klottey",         "Korle Korley",                                NA, NA, 0.0, "Point",
           "PCC", "Accra Metropolitan",    "Accra",                                       NA, NA, 0.0, "Point",
           
           "NCC", "Wenchi",                "Wenchi",      7.738011, -2.103519, 2.0, "Point",
           
           "HTT", "Bunkpurugu Nakpanduri", "Bunkpurugu", 10.547874, -0.029556, 13.2, "Point",
           "HTT", "Accra Metropolitan",    "Accra",       5.550000, -0.200000, 0.7, "Point",
           "HTT", "Kumasi",                "Kumasi",      6.697479, -1.631690, 0.4, "Range (0.1-0.7)",
           "HTT", "Kintampo North",        "Kintampo",    8.053714, -1.728919, 1.1, "Point",
           
           "HCC", "Bolgatanga",            "Bolgatanga", 10.785396, -0.848732, 6.7, "Point",
           "HCC", "Bawku West",            "Zebilla",    10.939946, -0.482999, 0.0, "Point",
           "HCC", "Kasena Nankana East",   "Navrongo",   10.892812, -1.088373, 0.0, "Point",
           "HCC", "Korle Klottey",         "Osu",         5.560593, -0.182439, 0.0, "Point",
           "HCC", "Weija Gbawe",           "Gbawe",       5.584882, -0.312821, 3.8, "Point",
           "HCC", "Ablekuma West",         "Bola Beach",  5.540000, -0.190000, 12.5, "Point",
           "HCC", "Ablekuma West",         "Chemunaa",    5.526600, -0.247200, 2.4, "Point",
           "HCC", "Accra",                 "James Town",  5.533000, -0.212000, 0.0, "Point",
           "HCC", "Ablekuma West",         "Kajaano",     5.560593, -0.182439, 33.3, "Point",
           "HCC", "Accra",                 "Korle Gonno", 5.531827, -0.229322, 0.0, "Point",
           "HCC", "Ayawaso Central",       "Pig Farm",    5.600200, -0.201620, 0.0, "Point",
           "HCC", "Ledzokuku",             "Teshie",      5.582234, -0.113510, 0.0, "Point",
           
           "HCC-HTT", "Kintampo North",     "Kintampo",    8.053714, -1.728919, 4.6, "Point"
         )
         
         # ==============================================================================
         # 2. DISTRICT NAME CLEANING
         # ==============================================================================
         
         clean_district <- function(x) {
           case_when(
             x %in% c("Kumasi Metropolitan", "KMA") ~ "Kumasi",
             x %in% c("Wa") ~ "Wa Municipal",
             x %in% c("Accra Metropolitan") ~ "Accra",
             x %in% c("Kasena-Nankana East") ~ "Kasena Nankana East",
             x %in% c("Kasena-Nankana West") ~ "Kasena Nankana West",
             x %in% c("Weija-Gbawe") ~ "Weija Gbawe",
             x %in% c("Korle-Klottey") ~ "Korle Klottey",
             TRUE ~ x
           )
         }
         
         # ==============================================================================
         # 3. MAP OBJECTS
         # ==============================================================================
         
         ghana_sf <- st_as_sf(ghana_vect) %>%
           mutate(District_match = clean_district(NAME_2))
         
         district_polygons <- ghana_sf %>%
           group_by(District_match) %>%
           summarise(geometry = st_union(geometry), .groups = "drop")
         
         district_inside_points <- district_polygons %>%
           st_point_on_surface() %>%
           mutate(
             inside_lon = st_coordinates(.)[, 1],
             inside_lat = st_coordinates(.)[, 2]
           ) %>%
           st_drop_geometry() %>%
           select(District_match, inside_lon, inside_lat)
         
         # ==============================================================================
         # 4. CLEAN PREVALENCE DATA
         # ==============================================================================
         
         prev_clean <- prev_districts %>%
           mutate(
             District_match = clean_district(District),
             Disease = factor(Disease, levels = c("PCC", "HCC", "HTT", "NCC", "HCC-HTT")),
             
             prevalence_group = case_when(
               Disease == "PCC" & prevalence > 0 & prevalence < 10       ~ "PCC: 0.1-9.99",
               Disease == "PCC" & prevalence >= 10 & prevalence < 25     ~ "PCC: 10-24.99",
               Disease == "PCC" & prevalence >= 25 & prevalence < 40     ~ "PCC: 25-39.99",
               Disease == "PCC" & prevalence >= 40                       ~ "PCC: >=40",
               
               Disease == "HCC" & prevalence > 0 & prevalence < 4        ~ "HCC: 0.1-3.99",
               Disease == "HCC" & prevalence >= 4 & prevalence < 8       ~ "HCC: 4-7.99",
               Disease == "HCC" & prevalence >= 8                        ~ "HCC: >=8",
               
               Disease == "HTT" & prevalence > 0 & prevalence <= 0.7     ~ "HTT: 0.1-0.7",
               Disease == "HTT" & prevalence > 0.7 & prevalence <= 5.5   ~ "HTT: 1.5-5.5",
               Disease == "HTT" & prevalence > 5.5                       ~ "HTT: >=44.6",
               
               Disease == "NCC"                                         ~ "NCC: cases only",
               Disease == "HCC-HTT"                                     ~ "HCC-HTT",
               TRUE ~ NA_character_
             ),
             
             point_size = case_when(
               prevalence_group == "PCC: 0.1-9.99" ~ 3.0,
               prevalence_group == "PCC: 10-24.99" ~ 4.4,
               prevalence_group == "PCC: 25-39.99" ~ 5.8,
               prevalence_group == "PCC: >=40"     ~ 7.2,
               prevalence_group == "HCC: 0.1-3.99" ~ 3.0,
               prevalence_group == "HCC: 4-7.99"   ~ 4.4,
               prevalence_group == "HCC: >=8"      ~ 5.8,
               prevalence_group == "HTT: 0.1-0.7"  ~ 3.0,
               prevalence_group == "HTT: 1.5-5.5"  ~ 4.4,
               prevalence_group == "HTT: >=44.6"   ~ 7.2,
               Disease %in% c("NCC", "HCC-HTT")    ~ 5.0,
               TRUE ~ 4.4
             )
           ) %>%
           left_join(district_inside_points, by = "District_match") %>%
           mutate(
             nudge_strength = 0.42,
             
             Longitude_plot = ifelse(
               !is.na(Longitude) & !is.na(inside_lon),
               Longitude + nudge_strength * (inside_lon - Longitude),
               Longitude
             ),
             Latitude_plot = ifelse(
               !is.na(Latitude) & !is.na(inside_lat),
               Latitude + nudge_strength * (inside_lat - Latitude),
               Latitude
             ),
             
             Longitude_plot = case_when(
               District_match == "Bunkpurugu Nakpanduri" & Disease == "PCC" ~ Longitude_plot - 0.035,
               District_match == "Bunkpurugu Nakpanduri" & Disease == "HTT" ~ Longitude_plot + 0.035,
               District_match == "Kintampo North" & Disease == "HTT"       ~ Longitude_plot - 0.035,
               District_match == "Kintampo North" & Disease == "HCC-HTT"   ~ Longitude_plot + 0.035,
               TRUE ~ Longitude_plot
             ),
             Latitude_plot = case_when(
               District_match == "Bunkpurugu Nakpanduri" & Disease == "PCC" ~ Latitude_plot + 0.020,
               District_match == "Bunkpurugu Nakpanduri" & Disease == "HTT" ~ Latitude_plot - 0.020,
               District_match == "Kintampo North" & Disease == "HTT"       ~ Latitude_plot + 0.020,
               District_match == "Kintampo North" & Disease == "HCC-HTT"   ~ Latitude_plot - 0.020,
               TRUE ~ Latitude_plot
             )
           )
         
         study_districts <- unique(prev_clean$District_match)
         
         ghana_study_sf <- ghana_sf %>%
           mutate(
             study_status = ifelse(
               District_match %in% study_districts,
               "Study district",
               "No study district"
             )
           )
         
         lakes_sf_study <- lakes_sf_correct %>%
           mutate(study_status = "lakes & lagoons")
         
         # ==============================================================================
         # 5. ZERO-PREVALENCE ASTERISKS
         # ==============================================================================
         
         district_centroids <- district_polygons %>%
           st_point_on_surface() %>%
           mutate(
             Longitude_star = st_coordinates(.)[, 1],
             Latitude_star  = st_coordinates(.)[, 2]
           ) %>%
           st_drop_geometry()
         
         pcc_zero_districts <- prev_clean %>%
           filter(Disease == "PCC", prevalence == 0) %>%
           distinct(District_match) %>%
           pull(District_match)
         
         hcc_zero_districts <- prev_clean %>%
           filter(Disease == "HCC", prevalence == 0) %>%
           distinct(District_match) %>%
           pull(District_match)
         
         zero_stars <- bind_rows(
           district_centroids %>%
             filter(District_match %in% pcc_zero_districts) %>%
             mutate(zero_type = "0% PCC prevalence in study district"),
           
           district_centroids %>%
             filter(District_match %in% hcc_zero_districts) %>%
             mutate(zero_type = "0% HCC prevalence in study district")
         )
         
         # ==============================================================================
         # 6. NON-ZERO PREVALENCE POINTS
         # ==============================================================================
         
         prev_points <- prev_clean %>%
           filter(prevalence > 0, !is.na(Latitude_plot), !is.na(Longitude_plot))
         
         prev_sf <- st_as_sf(
           prev_points,
           coords = c("Longitude_plot", "Latitude_plot"),
           crs = st_crs(ghana_sf),
           remove = FALSE
         )
         
         # ==============================================================================
         # 7. COLOURS AND SYMBOLS
         # ==============================================================================
         
         district_cols <- c(
           "Study district"    = "#8B5A2B",
           "No study district" = "#D9D9D9",
           "lakes & lagoons"   = "lightblue"
         )
         
         disease_cols <- c(
           "PCC"     = "#FFFF00",
           "HCC"     = "#0057FF",
           "HTT"     = "#FFFF00",
           "NCC"     = "#FFFF00",
           "HCC-HTT" = "#FFFF00"
         )
         
         zero_star_cols <- c(
           "0% PCC prevalence in study district" = "darkgreen",
           "0% HCC prevalence in study district" = "blue"
         )
         
         # ==============================================================================
         # 8. MAP EXTENT
         # ==============================================================================
         
         coord_xmin <- -5.50
         coord_xmax <-  2.80
         coord_ymin <-  4.30
         coord_ymax <- 11.85
         
         # ==============================================================================
         # 9. LABEL ANCHORS
         # ==============================================================================
         
         accra_inset_districts <- c(
           "Ablekuma West", "Weija Gbawe", "Korle Klottey",
           "Ledzokuku", "Accra", "Ayawaso Central"
         )
         
         main_label_anchors <- tribble(
           ~Village_Town, ~text_x, ~text_y, ~hjust, ~vjust,
           "Navrongo",    -2.65, 11.52, 0.5, 0.5,
           "Chiana",      -3.45, 11.18, 1.0, 0.5,
           "Namoo",       -1.95, 11.72, 0.5, 0.5,
           "Bongo",       -1.25, 11.52, 0.5, 0.5,
           "Bougosonge",  -0.65, 11.67, 0.5, 0.5,
           "Bolgatanga",   0.05, 11.50, 0.5, 0.5,
           "Widana",       0.95, 11.70, 0.5, 0.5,
           "Zebilla",      1.15, 10.65, 0.0, 0.5,
           "Garu",         1.15, 10.20, 0.0, 0.5,
           "Bunkpurugu",   1.15,  9.55, 0.0, 0.5,
           "Nalerigu",     1.15,  8.95, 0.0, 0.5,
           "Yendi",        1.15,  8.25, 0.0, 0.5,
           "Wa",          -3.35,  9.90, 1.0, 0.5,
           "Damango",     -3.35,  8.90, 1.0, 0.5,
           "Kintampo",    -3.35,  7.80, 1.0, 0.5,
           "Wenchi",      -3.35,  6.95, 1.0, 0.5,
           "Kumasi",      -3.35,  6.20, 1.0, 0.5
         )
         
         main_label_df <- prev_points %>%
           filter(!District_match %in% accra_inset_districts) %>%
           left_join(main_label_anchors, by = "Village_Town") %>%
           filter(!is.na(text_x)) %>%
           mutate(
             dx = text_x - Longitude_plot,
             dy = text_y - Latitude_plot,
             dist = sqrt(dx^2 + dy^2),
             arrow_xend = text_x - (dx / dist) * 0.12,
             arrow_yend = text_y - (dy / dist) * 0.12
           )
         
         # ==============================================================================
         # 10. ACCRA INSET SETUP
         # ==============================================================================
         
         accra_xlim <- c(-0.35, -0.05)
         accra_ylim <- c(5.45, 5.65)
         
         accra_rect_sf <- st_as_sf(
           st_as_sfc(
             st_bbox(
               c(xmin = accra_xlim[1], xmax = accra_xlim[2],
                 ymin = accra_ylim[1], ymax = accra_ylim[2]),
               crs = st_crs(ghana_sf)
             )
           )
         )
         
         accra_districts_sf <- ghana_study_sf %>%
           st_filter(accra_rect_sf, .predicate = st_intersects)
         
         lakes_sf_accra <- lakes_sf_study %>%
           st_filter(accra_rect_sf, .predicate = st_intersects)
         
         prev_sf_accra <- prev_sf %>%
           filter(District_match %in% accra_inset_districts)
         
         zero_stars_accra <- zero_stars %>%
           filter(District_match %in% accra_inset_districts)
         
         accra_label_anchors <- tribble(
           ~Village_Town, ~text_x, ~text_y, ~hjust, ~vjust,
           "Gbawe",       -0.315, 5.535, 0.5, 0.5,
           "Chemunaa",    -0.245, 5.468, 0.5, 0.5,
           "Bola Beach",  -0.165, 5.470, 0.5, 0.5,
           "Kajaano",     -0.125, 5.540, 0.5, 0.5,
           "Accra",       -0.165, 5.505, 0.5, 0.5
         )
         
         accra_label_df <- prev_points %>%
           filter(District_match %in% accra_inset_districts) %>%
           left_join(accra_label_anchors, by = "Village_Town") %>%
           filter(!is.na(text_x)) %>%
           mutate(
             dx = text_x - Longitude_plot,
             dy = text_y - Latitude_plot,
             dist = sqrt(dx^2 + dy^2),
             arrow_xend = text_x - (dx / dist) * 0.006,
             arrow_yend = text_y - (dy / dist) * 0.006
           )
         
         # ==============================================================================
         # 11. LEGEND DATA
         # ==============================================================================
         
         district_legend_df <- tribble(
           ~x, ~y, ~label, ~fill,
           1, 21.0, "Study district",    "#8B5A2B",
           1, 20.0, "No study district", "#D9D9D9",
           1, 19.0, "lakes & lagoons",   "lightblue"
         )
         
         disease_legend_df <- tribble(
           ~x, ~y, ~Disease, ~shape, ~fill,
           1, 16.6, "PCC",     21, "#FFFF00",
           1, 15.6, "HCC",     21, "#0057FF",
           1, 14.6, "HTT",     22, "#FFFF00",
           1, 13.6, "NCC",      8, "#FFFF00",
           1, 12.6, "HCC-HTT", 23, "#FFFF00"
         )
         
         magnitude_legend_df <- tribble(
           ~x, ~y,   ~heading, ~label,       ~shape, ~fill,      ~size,
           1, 10.0,  "PCC",    "0.1-9.99",    21,    "#FFFF00",   3.0,
           1,  9.1,  "PCC",    "10-24.99",    21,    "#FFFF00",   4.4,
           1,  8.2,  "PCC",    "25-39.99",    21,    "#FFFF00",   5.8,
           1,  7.3,  "PCC",    ">=40",        21,    "#FFFF00",   7.2,
           1,  5.9,  "HCC",    "0.1-3.99",    21,    "#0057FF",   3.0,
           1,  5.0,  "HCC",    "4-7.99",      21,    "#0057FF",   4.4,
           1,  4.1,  "HCC",    ">=8",         21,    "#0057FF",   5.8,
           1,  2.7,  "HTT",    "0.1-0.7",     22,    "#FFFF00",   3.0,
           1,  1.8,  "HTT",    "1.5-5.5",     22,    "#FFFF00",   4.4,
           1,  0.9,  "HTT",    ">=44.6",      22,    "#FFFF00",   7.2
         )
         
         zero_legend_df <- tribble(
           ~x, ~y, ~label, ~colour,
           1, -0.8, "PCC", "darkgreen",
           1, -1.8, "HCC", "blue"
         )
         
         type_legend_df <- tribble(
           ~x, ~y, ~label, ~dot,
           1, -3.6, "Point estimate", TRUE,
           1, -4.5, "Midpoint of range", FALSE
         )
         
         # ==============================================================================
         # 12. RIGHT LEGEND
         # ==============================================================================
         
         p_right_legend <- ggplot() +
           annotate("text", x = 0.4, y = 22.0, label = "District status",
                    hjust = 0, fontface = "plain", size = SZ_LEG_HEADING) +
           geom_tile(data = district_legend_df, aes(x = x, y = y, fill = fill),
                     width = 0.38, height = 0.65, colour = "grey35") +
           geom_text(data = district_legend_df, aes(x = 1.35, y = y, label = label),
                     hjust = 0, size = SZ_LEG_ITEM, fontface = "plain") +
           
           annotate("text", x = 0.4, y = 17.7, label = "Disease",
                    hjust = 0, fontface = "plain", size = SZ_LEG_HEADING) +
           geom_point(data = disease_legend_df,
                      aes(x = x, y = y, shape = factor(shape), fill = fill),
                      size = 6.0, colour = "black", stroke = 0.5) +
           geom_text(data = disease_legend_df,
                     aes(x = 1.35, y = y, label = Disease),
                     hjust = 0, size = SZ_LEG_ITEM, fontface = "plain") +
           
           annotate("text", x = 0.4, y = 11.1, label = "Prevalence magnitude (%)",
                    hjust = 0, fontface = "plain", size = SZ_LEG_HEADING) +
           annotate("text", x = 0.48, y = 10.45, label = "PCC",
                    hjust = 0, fontface = "plain", size = SZ_LEG_ITEM) +
           annotate("text", x = 0.48, y = 6.35, label = "HCC",
                    hjust = 0, fontface = "plain", size = SZ_LEG_ITEM) +
           annotate("text", x = 0.48, y = 3.15, label = "HTT",
                    hjust = 0, fontface = "plain", size = SZ_LEG_ITEM) +
           geom_point(data = magnitude_legend_df,
                      aes(x = x, y = y, shape = factor(shape), fill = fill, size = size * 1.1),
                      colour = "black", stroke = 0.45) +
           geom_text(data = magnitude_legend_df,
                     aes(x = 1.45, y = y, label = label),
                     hjust = 0, size = SZ_LEG_MAG, fontface = "plain") +
           
           annotate("text", x = 0.4, y = 0.0, label = "Site(s) with 0% prevalence",
                    hjust = 0, fontface = "plain", size = SZ_LEG_HEADING) +
           geom_text(data = zero_legend_df,
                     aes(x = x, y = y, label = "*", colour = colour),
                     size = 9, fontface = "bold") +
           geom_text(data = zero_legend_df,
                     aes(x = 1.35, y = y, label = label),
                     hjust = 0, size = SZ_LEG_ITEM, fontface = "plain") +
           
           annotate("text", x = 0.4, y = -2.7, label = "Prevalence type",
                    hjust = 0, fontface = "plain", size = SZ_LEG_HEADING) +
           geom_point(data = type_legend_df, aes(x = x, y = y),
                      shape = 21, size = 6.0, fill = "#FFFF00",
                      colour = "black", stroke = 0.5) +
           geom_point(data = type_legend_df %>% filter(dot),
                      aes(x = x, y = y),
                      shape = 16, size = 1.5, colour = "black") +
           geom_text(data = type_legend_df,
                     aes(x = 1.35, y = y, label = label),
                     hjust = 0, size = SZ_LEG_ITEM, fontface = "plain") +
           
           scale_fill_identity() +
           scale_colour_identity() +
           scale_shape_manual(values = c("8" = 8, "21" = 21, "22" = 22, "23" = 23)) +
           scale_size_identity() +
           coord_cartesian(xlim = c(0.3, 6.6), ylim = c(-5.0, 22.8), expand = FALSE) +
           theme_void() +
           theme(
             plot.background = element_rect(fill = "white", colour = NA),
             plot.margin = margin(8, 2, 8, 4),
             legend.position = "none"
           )
         
         # ==============================================================================
         # 13. MAIN MAP
         # ==============================================================================
         
         make_main_map <- function(map_title = "Ghana") {
           ggplot() +
             geom_sf(data = ghana_study_sf, aes(fill = study_status),
                     colour = "black", linewidth = 0.20, show.legend = FALSE) +
             geom_sf(data = lakes_sf_study, aes(fill = study_status),
                     colour = NA, inherit.aes = FALSE, show.legend = FALSE) +
             scale_fill_manual(values = district_cols, breaks = names(district_cols),
                               limits = names(district_cols), drop = FALSE) +
             ggnewscale::new_scale_fill() +
             geom_sf(data = prev_sf %>% filter(Disease %in% c("PCC", "HCC")),
                     aes(fill = Disease, size = point_size),
                     shape = 21, colour = "black", stroke = 0.45,
                     alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(Disease == "HTT"),
                     aes(fill = Disease, size = point_size),
                     shape = 22, colour = "black", stroke = 0.45,
                     alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(Disease == "HCC-HTT"),
                     aes(fill = Disease, size = point_size),
                     shape = 23, colour = "black", stroke = 0.45,
                     alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(Disease == "NCC"),
                     aes(fill = Disease, size = point_size),
                     shape = 8, colour = "black", stroke = 0.45,
                     alpha = 0.98, show.legend = FALSE) +
             geom_sf(data = prev_sf %>% filter(data_type == "Point"),
                     shape = 16, size = 1.0, colour = "black", show.legend = FALSE) +
             geom_segment(data = main_label_df,
                          aes(x = Longitude_plot, y = Latitude_plot,
                              xend = arrow_xend, yend = arrow_yend),
                          arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
                          linewidth = 0.6, colour = "black", inherit.aes = FALSE) +
             geom_text(data = main_label_df,
                       aes(x = text_x, y = text_y, label = Village_Town,
                           hjust = hjust, vjust = vjust),
                       size = SZ_MAIN_LABEL, fontface = "plain",
                       colour = "black", inherit.aes = FALSE) +
             geom_text(data = zero_stars,
                       aes(x = Longitude_star, y = Latitude_star,
                           label = "*", colour = zero_type),
                       size = 9, fontface = "bold", show.legend = FALSE) +
             scale_fill_manual(values = disease_cols) +
             scale_colour_manual(values = zero_star_cols) +
             scale_size_identity() +
             coord_sf(xlim = c(coord_xmin, coord_xmax),
                      ylim = c(coord_ymin, coord_ymax),
                      expand = FALSE) +
             labs(title = map_title) +
             theme_void() +
             theme(
               plot.title = element_text(size = SZ_MAP_TITLE, face = "bold", hjust = 0.5),
               plot.margin = margin(6, 6, 6, 6),
               legend.position = "none"
             )
         }
         
         # ==============================================================================
         # 14. ACCRA INSET MAP
         # ==============================================================================
         
         make_inset_map <- function() {
           ggplot() +
             geom_sf(data = accra_districts_sf, aes(fill = study_status),
                     colour = "black", linewidth = 0.35, show.legend = FALSE) +
             geom_sf(data = lakes_sf_accra, aes(fill = study_status),
                     colour = NA, inherit.aes = FALSE, show.legend = FALSE) +
             scale_fill_manual(values = district_cols, breaks = names(district_cols),
                               limits = names(district_cols), drop = FALSE) +
             ggnewscale::new_scale_fill() +
             geom_sf(data = prev_sf_accra %>% filter(Disease %in% c("PCC", "HCC")),
                     aes(fill = Disease, size = point_size),
                     shape = 21, colour = "black", stroke = 0.5,
                     alpha = 0.95, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(Disease == "HTT"),
                     aes(fill = Disease, size = point_size),
                     shape = 22, colour = "black", stroke = 0.5,
                     alpha = 0.95, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(Disease == "HCC-HTT"),
                     aes(fill = Disease, size = point_size),
                     shape = 23, colour = "black", stroke = 0.5,
                     alpha = 0.95, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(Disease == "NCC"),
                     aes(fill = Disease, size = point_size),
                     shape = 8, colour = "black", stroke = 0.5,
                     alpha = 0.95, show.legend = FALSE) +
             geom_sf(data = prev_sf_accra %>% filter(data_type == "Point"),
                     shape = 16, size = 1.0, colour = "black", show.legend = FALSE) +
             geom_segment(data = accra_label_df,
                          aes(x = Longitude_plot, y = Latitude_plot,
                              xend = arrow_xend, yend = arrow_yend),
                          arrow = arrow(length = unit(0.08, "cm"), type = "closed"),
                          linewidth = 0.4, colour = "black", inherit.aes = FALSE) +
             geom_text(data = accra_label_df,
                       aes(x = text_x, y = text_y, label = Village_Town,
                           hjust = hjust, vjust = vjust),
                       size = SZ_INSET_LABEL, fontface = "plain",
                       colour = "black", inherit.aes = FALSE) +
             geom_text(data = zero_stars_accra,
                       aes(x = Longitude_star, y = Latitude_star,
                           label = "*", colour = zero_type),
                       size = 8, fontface = "bold", show.legend = FALSE) +
             scale_fill_manual(values = disease_cols) +
             scale_colour_manual(values = zero_star_cols) +
             scale_size_identity() +
             coord_sf(xlim = c(-0.37, -0.05),
                      ylim = c(5.45, 5.65),
                      expand = FALSE) +
             labs(title = "Greater Accra") +
             theme_void() +
             theme(
               plot.title = element_text(size = 14, face = "plain", hjust = 0.5, margin = margin(b = 5)),
               panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
               plot.background = element_rect(fill = "white", colour = NA),
               legend.position = "none"
             )
         }
         
         # ==============================================================================
         # 15. LAYOUT
         # ==============================================================================
         
         build_layout <- function(p_main, p_inset) {
           p_inset_box <- p_inset + theme(plot.margin = margin(5, 5, 5, 5))
           
           (p_main | (plot_spacer() / p_inset_box) | p_right_legend) +
             plot_layout(widths = c(6.5, 2, 1.8))
         }
         
         final_map <- build_layout(
           make_main_map("Ghana"),
           make_inset_map()
         )
         
         # ==============================================================================
         # 16. ACCRA INSET CONNECTOR
         # ==============================================================================
         
         pw_main_w  <- 6.5
         pw_mid_w   <- 2
         pw_leg_w   <- 1.8
         pw_total_w <- pw_main_w + pw_mid_w + pw_leg_w
         
         main_frac  <- pw_main_w / pw_total_w
         inset_x_lo <- main_frac - 0.015
         
         map_x_range <- coord_xmax - coord_xmin
         map_y_range <- coord_ymax - coord_ymin
         
         title_frac  <- 0.050
         bot_frac    <- 0.020
         plot_h_frac <- 1.0 - title_frac - bot_frac
         
         rect_cx_frac <- (mean(accra_xlim) - coord_xmin) / map_x_range
         rect_cy_frac <- (mean(accra_ylim) - coord_ymin) / map_y_range
         
         dot1_x <- main_frac * (0.012 + rect_cx_frac * (1 - 2 * 0.012))
         dot1_y <- bot_frac + rect_cy_frac * plot_h_frac
         
         inset_y_ctr <- 0.18
         
         draw_connector_npc <- function() {
           grid::grid.lines(
             x = unit(c(dot1_x, inset_x_lo), "npc"),
             y = unit(c(dot1_y, inset_y_ctr), "npc"),
             gp = grid::gpar(col = "red", lty = "dashed", lwd = 1.8)
           )
           
           grid::grid.points(
             x = unit(dot1_x, "npc"),
             y = unit(dot1_y, "npc"),
             pch = 19,
             size = unit(0.008, "npc"),
             gp = grid::gpar(col = "red", fill = "red")
           )
           
           grid::grid.points(
             x = unit(inset_x_lo, "npc"),
             y = unit(inset_y_ctr, "npc"),
             pch = 19,
             size = unit(0.008, "npc"),
             gp = grid::gpar(col = "red", fill = "red")
           )
         }
         
         # ==============================================================================
         # 17. EXPORT PNG ONLY
         # ==============================================================================
         
         save_png <- function(obj, path, width, height, dpi, bg = "white") {
           ragg::agg_png(
             filename = path,
             width = width,
             height = height,
             units = "in",
             res = dpi,
             background = bg
           )
           
           print(obj)
           draw_connector_npc()
           
           dev.off()
           
           invisible(path)
         }
         
         save_png(
           final_map,
           file.path(out_dir, "Ghana_StudyDistrict_Prevalence_Map_V5.png"),
           width = 14,
           height = 9,
           dpi = 300
         )
         
         cat("Done. PNG exported successfully.\n")
         
         
         
         
         
         