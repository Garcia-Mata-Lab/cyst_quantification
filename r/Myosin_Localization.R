# Script name: Myosin_Localization.R
# Compiles radial localization data from "Rotating_Cyst_Profile.ijm" in Fiji/ImageJ
# Outputs CSV files summarizing the normalized data
# Outputs line graphs with shaded SEM showing a normalized profile of protein intensity from the center of the cyst outwards
# Author: Dr. Madeline Lovejoy and Dr. GPT
# Date: 2026-08-24
# R version 4.1.2 (2021-11-01)

# tidyverse version 2.0.0

#### Load Libraries ####
# if you have never used this package in your R console, you will need to install it first with the code: install.packages("tidyverse")
library(tidyverse)

#### Sets Working Directory to the folder that has this R code ####
folder <- getwd()

#### Load CSV Files ending in "Myo.csv" ####
Files <- list.files(
  path = folder,
  pattern = "Myo.csv",
  recursive = TRUE,
  full.names = TRUE
)

#### Names the wanted parameters for normalization as variables ####
target_length <- 100      # Number of interpolated points
distance_max <- 10        # Normalize distance from 0 to 10
drop_fraction <- 0.15     # Threshold to detect cyst boundary (15% of peak)
buffer_rows <- 10         # Rows to retain beyond detected boundary of cyst

#### Outlier Threshold ####
outlier_z_threshold <- 2.5

#### Initialize List #### WILL BE 0 AFTER THIS STEP
cyst_list <- list()

# ===============================
# PROCESS EACH CYST
# ===============================
for (f in Files) {
 
  df <- read_csv(f, show_col_types = FALSE)
 
  # Skip files that do not have measurements for all 360 degrees (columns)
  if (ncol(df) < 360) {
    warning(paste("Skipping file:", basename(f)))
    next
  }
 
  # Convert first 360 columns to a matrix
  data_matrix <- as.matrix(df[, 1:360])
 
  # Compute mean across angles (radial profile)
  radial_profile <- rowMeans(data_matrix, na.rm = TRUE)
 
  # Identify center of the cyst
  pixels <- length(radial_profile)
  center_pixel <- floor(pixels / 2)
# Extract half profile (center to the edge of the cyst)
  half_profile_raw <- radial_profile[(center_pixel + 1):pixels]
 
  # Smooth profile to reduce noise
  half_profile_smooth <- stats::filter(
    half_profile_raw,
    rep(1/5, 5),
    sides = 2
  )
  half_profile_smooth <- as.numeric(half_profile_smooth)
 
  # Replace NA values introduced by smoothing
  half_profile_smooth[is.na(half_profile_smooth)] <-
    half_profile_raw[is.na(half_profile_smooth)]
 
  # Detect cyst boundary based on intensity drop
  peak_intensity <- max(half_profile_smooth, na.rm = TRUE)
  threshold <- drop_fraction * peak_intensity
 
  drop_index <- which(half_profile_smooth < threshold)[1]
 
  # Determine cutoff index to discard measurements after detection of the cyst edge
  if (!is.na(drop_index)) {
    cutoff_index <- min(drop_index + buffer_rows,
                        length(half_profile_raw))
  } else {
    cutoff_index <- length(half_profile_raw)
  }
 
  # Trim data to remove measurements outside the cyst
  half_profile_trimmed <- half_profile_raw[1:cutoff_index]
 
  # Compute raw peak intensity from trimmed data
  raw_peak <- max(half_profile_trimmed, na.rm = TRUE)
 
  # Skip if peak is invalid
  if (is.na(raw_peak) || raw_peak == 0) {
    warning(paste("Invalid peak in file:", basename(f)))
    next
  }
 
  # Normalize each cyst to its own peak
  half_profile_norm <- half_profile_trimmed / raw_peak
 
  # Interpolate to standardized distance of 100 points
  interp_profile <- approx(
    x = seq_along(half_profile_norm),
    y = half_profile_norm,
    xout = seq(1, length(half_profile_norm),
               length.out = target_length)
  )$y
 
  # Standardized distance from center
  distance_norm <- seq(0, distance_max,
                       length.out = target_length)
 
  # Extract condition from filename
  condition <- case_when(
    str_detect(f, "CTRL") ~ "CTRL",
    str_detect(f, "SGEFKD") ~ "SGEF KD",
    str_detect(f, "WTResc") ~ "WT Rescue",
    TRUE ~ "Other"
  )
 
  # Store data
  cyst_list[[f]] <- tibble(
    Distance = distance_norm,
    IntensityNorm = interp_profile,
    RawPeak = raw_peak,
    Cyst = basename(f),
    Condition = condition
  )
}

#### Combine All Cysts for each condition ####
profile_df <- bind_rows(cyst_list)

# Ensure consistent condition order
profile_df$Condition <- factor(
  profile_df$Condition,
  levels = c("CTRL", "SGEF KD", "WT Rescue", "Other")
)

# ===============================
# CALCULATE RAW PEAK STATISTICS
# ===============================
peak_df <- profile_df %>%
  distinct(Cyst, Condition, RawPeak)

condition_peaks <- peak_df %>%
  group_by(Condition) %>%
  summarise(
    MeanRawPeak = mean(RawPeak, na.rm = TRUE),
    .groups = "drop"
  )

# Extract CTRL reference peak (will normalize the intensity for the other conditions to CTRL)
ctrl_mean_peak <- condition_peaks %>%
  filter(Condition == "CTRL") %>%
  pull(MeanRawPeak)

if (length(ctrl_mean_peak) == 0) {
  stop("No CTRL samples found. Check file naming.")
}

# Compute scaling factors relative to CTRL
condition_peaks <- condition_peaks %>%
  mutate(
    ScalingFactor = MeanRawPeak / ctrl_mean_peak
  )

# Will show in the console the peak for each condition
print(condition_peaks)

# ===============================
# APPLY SCALING FACTORS
# ===============================
profile_df <- profile_df %>%
  left_join(
    condition_peaks %>% dplyr::select(Condition, ScalingFactor),
    by = "Condition"
  ) %>%
  mutate(
    IntensityScaled = IntensityNorm * ScalingFactor
  )

# ===============================
# REMOVE EXTREME PROFILE OUTLIERS
# ===============================

# Compute provisional mean profile
temp_summary <- profile_df %>%
  group_by(Condition, Distance) %>%
  summarise(
    MeanProfile = mean(IntensityScaled, na.rm = TRUE),
    .groups = "drop"
  )

# Compare each cyst to its condition mean to normalize values
deviation_df <- profile_df %>%
  left_join(
    temp_summary,
    by = c("Condition", "Distance")
  ) %>%
  group_by(Cyst, Condition) %>%
  summarise(
    ProfileDeviation = mean(
      abs(IntensityScaled - MeanProfile),
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Compute Z-scores within each condition
outlier_df <- deviation_df %>%
  group_by(Condition) %>%
  mutate(
    DeviationZ = (
      ProfileDeviation -
        mean(ProfileDeviation, na.rm = TRUE)
    ) / sd(ProfileDeviation, na.rm = TRUE)
  ) %>%
  ungroup()

# Identify outliers
outliers_to_remove <- outlier_df %>%
  filter(DeviationZ > outlier_z_threshold)

# Print removed cysts in the Console
cat("\n===============================\n")
cat("REMOVED OUTLIER CYSTS\n")
cat("===============================\n")

print(outliers_to_remove)

# Remove outlier cysts
profile_df <- profile_df %>%
  filter(!Cyst %in% outliers_to_remove$Cyst)

# ===============================
# RECALCULATE MEAN AND SEM after removing outliers
# ===============================
summary_df <- profile_df %>%
  group_by(Condition, Distance) %>%
  summarise(
    n = sum(!is.na(IntensityScaled)),
    Mean = mean(IntensityScaled, na.rm = TRUE),
    SEM = sd(IntensityScaled, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

# ===============================
# GENERATE MEAN ± SEM PLOT
# ===============================
condition_colors <- c(
  "CTRL" = "blue",
  "SGEF KD" = rgb(255, 128, 0, maxColorValue = 255), # orange
  "WT Rescue" = rgb(0, 215, 149, maxColorValue = 255), # turquoise
  "Other" = "gray50"
)

# Graphing
p <- ggplot(
  summary_df,
  aes(
    x = Distance,
    y = Mean,
    color = Condition,
    fill = Condition
  )
) +
 
 # Line thickness  
geom_line(linewidth = 1.5) +
 
# SEM ribbon  
geom_ribbon(
    aes(
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
  # Transparency
    alpha = 0.2,
    color = NA
  ) +

# Uses the colors assigned for each condition
  scale_color_manual(values = condition_colors) +
 
  scale_fill_manual(values = condition_colors) +

# Axes labels
  labs(
    x = "Normalized distance from cyst center",
    y = "Relative myosin intensity"
  ) +
 
 # Font size 
theme_classic(base_size = 16) +
 
  # Transparent graph background
theme(
    legend.title = element_blank(),
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    )
  )

# Display plot
p

# ===============================
# GLOBAL Y-AXIS LIMITS
# ===============================

global_y_max <- max(
  profile_df$IntensityScaled,
  na.rm = TRUE
)

global_y_min <- min(
  profile_df$IntensityScaled,
  na.rm = TRUE
)
# ===============================
# INDIVIDUAL CYST PLOTS - Average line and a transparent line for each cyst - one graph per condition
# ===============================
plots_by_condition <- lapply(
  levels(profile_df$Condition),
  function(cond) {
   
    df_cond <- profile_df %>%
      filter(Condition == cond)
   
    if (nrow(df_cond) == 0) return(NULL)
   
    mean_df <- df_cond %>%
      group_by(Distance) %>%
      summarise(
        Mean = mean(IntensityScaled, na.rm = TRUE),
        .groups = "drop"
      )
   
    ggplot() +
     
      geom_line(
        data = df_cond,
        aes(
          x = Distance,
          y = IntensityScaled,
          group = Cyst
        ),
        color = condition_colors[cond],
        linewidth = 0.6,
        alpha = 0.35
      ) +
     
      geom_line(
        data = mean_df,
        aes(
          x = Distance,
          y = Mean
        ),
        color = condition_colors[cond],
        linewidth = 1.6
      ) +
     
      labs(
        title = cond,
        x = "Normalized distance from cyst center",
        y = "Relative myosin intensity"
      ) +
     
      coord_cartesian(
        ylim = c(0, 1.25)
      ) +
      theme_classic(base_size = 16) +
     
      theme(
        legend.position = "none",
        plot.title = element_text(
          hjust = 0.5,
          face = "bold"
        )
      )
  }
)

names(plots_by_condition) <- levels(profile_df$Condition)

# Display plots
plots_by_condition

# ===============================
# SAVE OUTPUTS
# ===============================
ggsave(
  filename = file.path(
    folder,
# ***CHANGE THIS NAME DEPENDING ON HOW YOU WOULD LIKE YOUR GRAPH SAVED - THIS GRAPH IS ALL CONDITIONS COMBINED***
    "Exp3_Combined_Myo_Intensity_Plot.pdf"
  ),
  plot = p,
  width = 6,
  height = 4,
  units = "in",
  dpi = 300
)

for (cond in names(plots_by_condition)) {
 
  plot_obj <- plots_by_condition[[cond]]
 
  if (is.null(plot_obj)) next
 
  safe_name <- gsub(" ", "_", cond)
 
  ggsave(
    filename = file.path(
      folder,
      paste0(
        "Exp3_",
        safe_name,
# ***CHANGE THIS NAME DEPENDING ON HOW YOU WOULD LIKE YOUR GRAPHS SAVED - THIS WILL SAVE ONE GRAPH PER CONDITION WITH THE INDIVIDUAL RADIAL PROFILES***
        "_Individual_Cyst_Profiles.pdf"
      )
    ),
    plot = plot_obj,
    width = 6,
    height = 4,
    units = "in",
    dpi = 300
  )
}

write_csv(
  profile_df,
  file.path(
    folder,
# ***CHANGE THE EXPERIMENT NUMBER, BUT KEEP THE REST OF THE NAME THE SAME***
    "Exp3_Normalized_Cyst_Profiles.csv"
  )
)

write_csv(
  summary_df,
  file.path(
    folder,
# ***CHANGE THIS NAME***
    "Exp3_Mean_SEM_Cyst_Profiles.csv"
  )
)

write_csv(
  condition_peaks,
  file.path(
    folder,
# ***CHANGE THIS NAME***
    "Exp3_Condition_Scaling_Factors.csv"
  )
)

write_csv(
  outlier_df,
  file.path(
    folder,
    "Exp3_Profile_Outlier_Analysis.csv"
  )
)


# THE FOLLOWING CODE SHOULD BE USED AFTER ALL YOUR EXPERIMENTS ARE QUANTIFIED - ***OUTPUT CSVS FROM THE PREVIOUS CODE SHOULD BE MOVED INTO ONE FOLDER***
# ***MANUALLY CHANGE WORKING DIRECTORY TO FOLDER WITH ALL EXPERIMENTS*** (Session -> Change working directory)
#----------------------------------
library(tidyverse)

folder_combined <- getwd()

files <- list.files(
  path = folder_combined,
  pattern = "Normalized_Cyst_Profiles\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if(length(files)==0){
  stop("No normalized profile CSV files found.")
}

# ===============================
# LOAD ALL EXPERIMENTS
# ===============================
combined_df <- lapply(files, function(f){
 
  df <- read_csv(f, show_col_types = FALSE)
 
  exp_name <- str_extract(
    basename(f),
    "Exp\\d+"
  )
 
  if(is.na(exp_name)){
    exp_name <- basename(f)
  }
 
  df %>%
    mutate(
      Experiment = exp_name
    )
 
}) %>%
  bind_rows()

combined_df$Condition <- factor(
  combined_df$Condition,
  levels=c(
    "CTRL",
    "SGEF KD",
    "WT Rescue",
    "Other"
  )
)

condition_colors <- c(
  "CTRL"="blue",
  "SGEF KD"=rgb(255, 128, 0, maxColorValue = 255),
  "WT Rescue"=rgb(0, 215, 149, maxColorValue = 255),
  "Other"="gray50"
)

# ===============================
# OVERALL MEAN ± SEM COMPARING ALL EXPERIMENTS
# ===============================
overall_summary <- combined_df %>%
  group_by(
    Condition,
    Distance
  ) %>%
  summarise(
    n = sum(!is.na(IntensityScaled)),
    Mean = mean(
      IntensityScaled,
      na.rm=TRUE
    ),
    SEM = sd(
      IntensityScaled,
      na.rm=TRUE
    )/sqrt(n),
    .groups="drop"
  )

# ===============================
# MEAN OF EACH EXPERIMENT
# ===============================
experiment_summary <- combined_df %>%
  group_by(
    Condition,
    Experiment,
    Distance
  ) %>%
  summarise(
    Mean = mean(
      IntensityScaled,
      na.rm=TRUE
    ),
    .groups="drop"
  )

# Line types for experiments for showing experiment averages
linetype_values <- c(
  "dashed",
  "dotdash",
  "twodash",
  "longdash",
  "dotted"
)

experiment_levels <- sort(
  unique(experiment_summary$Experiment)
)

experiment_linetypes <- setNames(
  linetype_values[
    seq_along(experiment_levels)
  ],
  experiment_levels
)

# ===============================
# COMBINED PLOT
# ===============================
combined_plot <- ggplot() +
 
  geom_ribbon(
    data = overall_summary,
    aes(
      x = Distance,
      ymin = Mean-SEM,
      ymax = Mean+SEM,
      fill = Condition
    ),
    alpha = 0.15,
    color = NA
  ) +
 
  geom_line(
    data = experiment_summary,
    aes(
      x = Distance,
      y = Mean,
      color = Condition,
      linetype = Experiment,
      group = interaction(
        Condition,
        Experiment
      )
    ),
    linewidth = 0.7,
    alpha = 0.45
  ) +
 
  geom_line(
    data = overall_summary,
    aes(
      x = Distance,
      y = Mean,
      color = Condition
    ),
    linewidth = 2
  ) +
 
  scale_color_manual(
    values = condition_colors
  ) +
 
  scale_fill_manual(
    values = condition_colors
  ) +
 
  scale_linetype_manual(
    values = experiment_linetypes
  ) +
 
  guides(
    linetype = "none",
    fill = "none",
    color = guide_legend(
      title = NULL,
      override.aes = list(
        linewidth = 2.5,
        linetype = "solid"
      )
    )
  ) +
 
  labs(
    x = "Normalized distance from cyst center",
# ***CHANGE Y-AXIS NAME DEPENDING ON WHICH PROTEIN YOU'RE QUANTIFYING***
    y = "Relative pMLC intensity"
  ) +
 
  theme_classic(base_size = 16) +
 
  theme(
    panel.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    plot.background = element_rect(
      fill = "transparent",
      color = NA
    ),
    legend.position = "right",
    # Axis titles
    axis.title.x = element_text(
      face = "bold",
      color = "black"
    ),
   
    axis.title.y = element_text(
      face = "bold",
      color = "black"
    ),
   
    # Axis tick labels (numbers)
    axis.text.x = element_text(
      face = "bold",
      color = "black"
    ),
   
    axis.text.y = element_text(
      face = "bold",
      color = "black"
    )
  )

combined_plot

# ===============================
# SAVE
# ===============================
output_dir <- file.path(
  folder_combined,
  "Combined_Results"
)

dir.create(
  output_dir,
  showWarnings = FALSE
)

ggsave(
  filename = file.path(
    output_dir,
    "Combined_Experiments_All_Conditions.pdf"
  ),
  plot = combined_plot,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "transparent"
)


# ***THIS CODE IS JUST ME TRYING TO FIGURE OUT IF I CAN TEST FOR STATISTICAL SIGNIFICANCE AT DIFFERENT POINTS IN THE GRAPH - YOU DO NOT NEED TO RUN ANY OF THE FOLLOWING CODE***
# ===============================
# CALCULATE SLOPE FOR EACH CYST
# (Distance = 3 to 6)
# ===============================

slope_df <- combined_df %>%
  filter(Distance >= 4,
         Distance <= 9) %>%
  group_by(
    Condition,
    Experiment,
    Cyst
  ) %>%
  summarise(
    Slope = coef(
      lm(IntensityScaled ~ Distance)
    )[2],
    .groups = "drop"
  )

# View results
print(slope_df)
# ===============================
# SAVE SLOPES
# ===============================

write_csv(
  slope_df,
  file.path(
    output_dir,
    "Radial_Profile_Slopes_4_to_9_apical.csv"
  )
)

# ===============================
# CALCULATE SLOPE FOR EACH CYST
# (Distance = 7.5 to 8.5)
# ===============================

slope_df <- combined_df %>%
  filter(Distance >= 7.5,
         Distance <= 8.5) %>%
  group_by(
    Condition,
    Experiment,
    Cyst
  ) %>%
  summarise(
    Slope = coef(
      lm(IntensityScaled ~ Distance)
    )[2],
    .groups = "drop"
  )

# View results
print(slope_df)

# ===============================
# SAVE SLOPES
# ===============================

write_csv(
  slope_df,
  file.path(
    output_dir,
    "Radial_Profile_Slopes_7_to_8.5_basal.csv"
  )
)
