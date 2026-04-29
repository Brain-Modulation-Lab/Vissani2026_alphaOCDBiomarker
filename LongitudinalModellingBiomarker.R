# =============================================================================
# 1. LIBRARIES & DIRECTORIES
# =============================================================================
library(tidyverse)
library(cowplot)
library(scales)
library(patchwork)
library(Hmisc)
library(readxl)
library(here)
library(boot)
library(viridis) # For the Turbo palette
library(lme4) 
library(lmerTest) 
library(sjPlot)
library(emmeans)
library(broom.mixed)

ROOT_NEXUS  <- "/Volumes/Nexus4"
PATH_OCD     <- file.path(ROOT_NEXUS, "OCD")
PATH_PROJECT <- file.path(ROOT_NEXUS, "your_folderpath")
PATH_SCRIPTS <- file.path(PATH_PROJECT, "your_script_path")

PATH_RESULTS <- file.path(PATH_PROJECT,"Results")
PATH_CLINICAL         <- file.path(PATH_RESULTS, "output/clinical")
PATH_FIGURES          <- file.path(PATH_RESULTS, "figures")
PATH_BRAINSENSESURVEY <- file.path(PATH_RESULTS, "output/brainsensesurvey")
PATH_THERAPYHISTORY   <- file.path(PATH_OCD, "groupanalyses/therapy_history")
PATH_OUTPUT <- file.path(PATH_FIGURES, "brainsensesurvey","all","Rstats")

# set working directory

# add PROJECT-SPECIFIC Utils.R
source(here("utils","utils.R"))

# Constants
FONT_FAMILY    <- "Helvetica"
nBoots         <- 10000
COLORS_COMPARE <- c("First Session (S1)" = "gray40", "Best Y-BOCS" = "#d7191c")
aligned_freq_vector <- seq(-32, 32, by = 0.5) # Based on your MATLAB setup
DAYS_IN_YEAR <- 365.25

# =============================================================================
# 2. SCRIPT-SPECIFIC HELPER FUNCTIONS
# =============================================================================



plot_trajectory <- function(input_data, title_label) {
  clim_val <- quantile(abs(input_data$ephys_delta), 0.98, na.rm = TRUE)
  
  # 1. Recalculate gapless integer rank
  plot_df <- input_data %>%
    arrange(dy_imp, patient_id, session_id, channel_id) %>%
    mutate(x_idx = as.numeric(factor(obs_id, levels = unique(obs_id))))
  
  # 2. Find threshold index
  v_line <- plot_df %>% filter(dy_imp >= 35) %>% arrange(dy_imp) %>% slice(1) %>% pull(x_idx)
  if(length(v_line) == 0) v_line <- NA
  
  # 3. Define EXACT X-limits for both plots
  # This ensures the 0.5 offset of geom_rect is accounted for in the ramp plot
  x_lims <- c(0.5, max(plot_df$x_idx) + 0.5)
  
  # Heatmap
  p_heat <- ggplot(plot_df) +
    geom_rect(aes(xmin=x_idx-0.5, xmax=x_idx+0.5, ymin=ymin, ymax=ymax, fill=ephys_delta)) +
    geom_vline(xintercept = v_line, linetype = "dashed", color = "black", linewidth = 0.8) +
    scale_y_log10(breaks = c(2, 4, 8, 12, 21, 35), expand = c(0,0)) +
    # Force X-axis limits
    scale_x_continuous(limits = x_lims, expand = c(0,0)) +
    coord_cartesian(ylim = c(2, 35)) +
    scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", 
                         limits = c(-clim_val, clim_val), oob = scales::squish) +
    theme_cowplot(font_size = 12) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank(), axis.line.x = element_blank(),
          plot.margin = margin(5, 5, 0, 5)) +
    labs(y = "Frequency [Hz]", title = title_label, fill = "Δ Power")
  
  # Ramp
  p_ramp <- ggplot(plot_df %>% distinct(x_idx, dy_imp), aes(x = x_idx, y = dy_imp)) +
    geom_area(fill = "gray85", color = "black", linewidth = 0.5) +
    geom_point(shape = 21, fill = "white", color = "black", size = 2.2, stroke = 0.4) +
    geom_hline(yintercept = 35, linetype = "dashed", color = "darkred", linewidth = 0.8) +
    geom_vline(xintercept = v_line, linetype = "dashed", color = "black", linewidth = 0.8) +
    # Use the same exact X-limits as the heatmap
    scale_x_continuous(limits = x_lims, expand = c(0,0)) +
    theme_cowplot(font_size = 12) +
    theme(plot.margin = margin(0, 5, 5, 5)) +
    labs(y = "Y-BOCS imp [%S1]", x = "# recording")
  
  # 4. Use patchwork to collect and align
  # 'align = "v"' is the key to matching the plot areas regardless of Y-axis label length
  return(wrap_plots(p_heat, p_ramp, ncol = 1, heights = c(4, 1)) & 
           theme(plot.margin = margin(5, 5, 5, 5)))
}

# Data loading
DB <- load_file(file.path(PATH_BRAINSENSESURVEY, "OCD_survey-all_bipolar-alpha_cons.csv"))
STIM_PERIODS <- load_file(file.path(PATH_THERAPYHISTORY, "platform-refined-20260223","sub-all_ses-postop_therapyhistory.csv"))

target_vector_vars <- c("Frequency", "FrequencyFlat", "MeanPower", "MeanPowerNorm", 
                        "MeanPowerFlat", "MeanPowerFlatAlphaAligned", "MeanPowerFlatBetaAligned")
target_pattern <- paste0("^(", paste(target_vector_vars, collapse = "|"), ")_")

DF_EPHYS_LONG <- DB %>%
  select(-matches("preVisit|postVisit")) %>%
  rename_with(~make.unique(.), everything()) %>% 
  pivot_longer(
    cols = matches(target_pattern),
    names_to = c(".value", "bin_index"),
    names_pattern = "(.*)_(\\d+)$"
  ) %>%
  mutate(bin_index = as.numeric(bin_index),
         patient_id = as.character(patient_id),
         alignedFrequency = aligned_freq_vector[bin_index]) %>%
  filter(!is.na(Frequency) | !is.na(MeanPowerFlatAlphaAligned))

# Pivot to Long, Standardize Hemisphere, and Fill Asymmetric NAs
STIM_PREPARED <- STIM_PERIODS %>%
  mutate(patient_id = as.character(patient_id)) %>%
  filter(days_end_after_stimon >= days_start_after_stimon) %>%
  pivot_longer(
    cols = matches("^[LR]_"),
    names_to = c("Hemisphere", ".value"),
    names_pattern = "(L|R)_(.*)"
  ) %>%
  mutate(Hemisphere = ifelse(Hemisphere == "L", "Left", "Right")) %>%
  rename(Amp = Amp, Freq = Freq, PW = PW) %>%
  # Fill settings forward for each side (handles the NA issue in OP07)
  group_by(patient_id, Hemisphere) %>%
  arrange(days_start_after_stimon) %>%
  fill(Amp, Freq, PW, Contacts, Group, .direction = "down") %>%
  # Extend last period to Infinity for late sessions
  mutate(days_end_after_stimon = if_else(
    days_end_after_stimon == max(days_end_after_stimon, na.rm = TRUE), 
    Inf, 
    days_end_after_stimon
  )) %>%
  ungroup()

DF_EPHYS_STIM <- DF_EPHYS_LONG %>%
  # Use left_join so Session 1 is preserved even if it doesn't match a stim period
  left_join(
    STIM_PREPARED, 
    by = join_by(
      patient_id, Hemisphere,
      days_after_stimon >= days_start_after_stimon,
      days_after_stimon <= days_end_after_stimon
    )
  ) #%>%
# # Handle overlaps for sessions that DID match multiple stim periods
# group_by(patient_id, session_id, Hemisphere, bin_index) %>%
# # We use slice_max but allow for NAs (Session 1)
# arrange(desc(duration_days)) %>%
# slice(1) %>%
# ungroup()

# Filter for patients with at least 2 sessions
valid_patients <- DF_EPHYS_STIM %>%
  group_by(patient_id) %>%
  summarise(n_sess = n_distinct(session_id)) %>%
  filter(n_sess >= 2) %>%
  pull(patient_id)

SELECTED_SESSIONS <- DF_EPHYS_STIM %>%
  group_by(patient_id, Hemisphere) %>%
  summarise(
    session_first = 1, 
    session_best  = get_best_session(session_id, days_after_stimon, ybocs),
    .groups = "drop"
  ) %>%
  filter(!is.na(session_best))

DF_EPHYS_STIM <- DF_EPHYS_STIM %>% filter(patient_id %in% valid_patients) %>%
  mutate(channel_id = paste0("N", Sensing_ChannelbipN, "_P", Sensing_ChannelbipP),
         years_after_stimon = days_after_stimon / DAYS_IN_YEAR) %>% inner_join(SELECTED_SESSIONS, by = c("patient_id", "Hemisphere"))

# Prepare data for plot
AUDIT_EPHYS <- DF_EPHYS_STIM %>%
  distinct(patient_id, session_id, Hemisphere, days_after_stimon, Amp, Group) %>%
  mutate(years_after_stimon = days_after_stimon / DAYS_IN_YEAR)

# Cap the Infinite horizon for plotting purposes

AUDIT_STIM_LOG <- STIM_PREPARED %>%
  filter(patient_id %in% unique(AUDIT_EPHYS$patient_id)) %>%
  mutate(
    years_start = days_start_after_stimon / DAYS_IN_YEAR,
    years_end   = if_else(is.infinite(days_end_after_stimon), 
                          max(AUDIT_EPHYS$years_after_stimon, na.rm = TRUE) + 0.1, 
                          days_end_after_stimon / DAYS_IN_YEAR)
  )

# Create the Audit Plot
plot_audit <- ggplot() +
  geom_rect(data = AUDIT_STIM_LOG,
            aes(xmin = years_start, xmax = years_end,
                ymin = Amp - 0.2, ymax = Amp + 0.2, fill = Group), 
            alpha = 0.3) +
  geom_point(data = AUDIT_EPHYS,
             aes(x = years_after_stimon, y = Amp), 
             color = "black", size = 2, shape = 4, stroke = 0.8) +
  facet_grid(patient_id ~ Hemisphere, scales = "free_x") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = pretty_breaks(n = 5)) +
  theme_cowplot(font_size = 10, font_family = FONT_FAMILY) +
  theme(strip.text.y = element_text(angle = 0), legend.position = "bottom") +
  labs(title = "Stimulation Assignment Audit",
       subtitle = "X = Joined Ephys Sessions | Shaded = Massaged Stim Log",
       x = "Time after DBS ON [years]", y = "Amplitude (mA)")

# Output Results
print(plot_audit)


# Plot overall heatmap of psd with ybocs change and ephysio change w.r.t. S1 ----
DF_YBOCS_S1 <- DF_EPHYS_STIM %>%
  filter(session_id == 1) %>%
  distinct(patient_id, ybocs) %>%
  rename(ybocs_s1_global = ybocs)

DF_SPEC_CHANGES <- DF_EPHYS_STIM %>%
  left_join(DF_YBOCS_S1, by = "patient_id") %>%
  group_by(patient_id, Hemisphere, channel_id, bin_index) %>%
  arrange(session_id) %>%
  mutate(
    power_s1 = first(MeanPowerFlat),
    ephys_delta = MeanPowerFlat - power_s1,
    dy_imp = -100 * (ybocs - ybocs_s1_global) / ybocs_s1_global,
    obs_id = paste(patient_id, session_id, Hemisphere, channel_id, sep = "_"),
    ymin = FrequencyFlat - 0.5, ymax = FrequencyFlat + 0.5
  ) %>%
  # Crucial: Filter NAs first to prevent white strips in the heatmap
  filter(session_id > 1, !is.na(ephys_delta)) %>%
  ungroup()

# To see the Pooled version:
fig_pooled <- plot_trajectory(DF_SPEC_CHANGES, "Pooled Clinical Trajectory (L+R)")
print(fig_pooled)

# To see the Left side version:
fig_left <- plot_trajectory(DF_SPEC_CHANGES %>% filter(Hemisphere == "Left"), "Left Hemisphere Only")
print(fig_left)

# To see the Right side version:
fig_right <- plot_trajectory(DF_SPEC_CHANGES %>% filter(Hemisphere == "Right"), "Right Hemisphere Only")
print(fig_right)

# 2] Plotting ybocs trajectory with stim settings and spectra/psd

# Power fields to summarise (keep for completeness)
ephysio_fields <- c(
  "MeanPowerAlphaToT", "MeanPowerAlpha", "MeanPowerAlphaFlat",
  "MaxPowerAlphaToT",  "MaxPowerAlpha",  "MaxPowerAlphaFlat",
  "MeanPowerBetaToT",  "MeanPowerBeta",  "MeanPowerBetaFlat",
  "MaxPowerBetaToT",   "MaxPowerBeta",   "MaxPowerBetaFlat"
)
# Define the stim fields you want to keep
stim_fields <- c("period_id", "Group", "days_start_after_stimon", 
                 "days_end_after_stimon", "duration_days", 
                 "Amp", "Freq", "PW", "Contacts")

patients_to_highlight <- tibble(
  patient_id =  c("OP12", "OP03", "OP04", "OP06"),
  ybocs = c(31, 36, 33, 30),
  years_after_stimon = -150 / DAYS_IN_YEAR
)

# Process the Longitudinal Data
DF_CLIN_STIM <- DF_EPHYS_STIM %>% 
  left_join(DF_YBOCS_S1, by = "patient_id") %>%
  mutate(
    ybocs_s1 = ybocs_s1_global,
    delta_improvement = ybocs_s1 - ybocs,
    pct_improvement = 100 * (delta_improvement) / ybocs_s1
  ) %>%
  group_by(patient_id, session_id, Hemisphere) %>%
  summarise(
    ybocs = first(ybocs),
    days_after_stimon = first(days_after_stimon),
    # Ensure timing is in years for alignment
    years_after_stimon = first(days_after_stimon) / DAYS_IN_YEAR, 
    delta_improvement = first(delta_improvement),
    
    pct_improvement = first(pct_improvement),
    session_first = first(session_first),
    session_best = first(session_best),
    across(any_of(stim_fields), first),
    
    across(
      any_of(ephysio_fields),
      list(
        hemiMax  = ~ safe_max(.x),
        hemiMean = ~ safe_mean(.x)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )



# =============================================================================
# 1. CALCULATE SCALED DELTA FOR BOTH HEMISPHERES
# =============================================================================
# This avoids hardcoding and just looks for what you already computed
computed_fields <- colnames(DF_CLIN_STIM)[grep("_hemiMax$|_hemiMean$", colnames(DF_CLIN_STIM))]

# 3. Compute Baselines and Scaled Deltas
DF_OVERLAY_DATA <- DF_CLIN_STIM %>%
  group_by(patient_id, Hemisphere) %>%
  mutate(
    # Step A: Get S1 Baseline for every computed field
    across(all_of(computed_fields), 
           ~ .[session_id == session_first][1], 
           .names = "{.col}_s1"),
    
    # Step B: Calculate Raw Deltas
    across(all_of(computed_fields), 
           ~ . - get(paste0(cur_column(), "_s1")), 
           .names = "{.col}_delta") 
    ) %>%
  ungroup()

# get scaled value for exemplary plots
DF_EXAMPLES <- DF_OVERLAY_DATA %>% 
  filter(Hemisphere == "Left",patient_id %in% patients_to_highlight$patient_id) %>%
  mutate (
  # Step C: Scale Deltas (Standardized to +/- 12 units centered on 25)
  across(ends_with("_delta"), 
         ~ {
           m <- max(abs(.), na.rm = TRUE)
           if(is.na(m) || m == 0) rep(0, n()) else (. / m) * 12
         }, 
         .names = "{gsub('_delta', '_scaled', .col)}")
  )

# 1. Update Stim Data coordinates (Rug sits between -10 and -2)
DF_STIM_FINAL <- STIM_PERIODS %>%
  mutate(patient_id = as.character(patient_id),
         x_start = days_start_after_stimon / DAYS_IN_YEAR,
         x_end = days_end_after_stimon / DAYS_IN_YEAR) %>%
  filter(x_end > x_start) %>%
  pivot_longer(cols = c(L_Amp, R_Amp), names_to = "side", values_to = "amp") %>%
  mutate(
    y_min = if_else(side == "L_Amp", -5, -10),
    y_max = if_else(side == "L_Amp", -2, -7)     
  )


final_faceted_plot <- ggplot() +
  
  # A. Bilateral Stim Rug
  geom_rect(data = DF_STIM_FINAL %>%
              filter(patient_id %in% patients_to_highlight$patient_id), 
            aes(xmin = x_start, xmax = x_end, ymin = y_min, ymax = y_max, fill = amp), 
            color = "white", linewidth = 0.1) + 
  
  # B. Margin Labels
  geom_text(data = DF_STIM_FINAL %>%
              filter(patient_id %in% patients_to_highlight$patient_id) %>% 
              group_by(patient_id, side) %>% slice(1),
            aes(x = -0.6, y = (y_min + y_max)/2, label = if_else(side == "L_Amp", "L", "R")),
            size = 3, fontface = "bold", hjust = 1) +
  
  # C. NEURAL OVERLAYS (Anchored to Initial Y-BOCS)
  geom_line(data = DF_EXAMPLES, 
            aes(x = years_after_stimon, 
                # Calculation: Start at first Y-BOCS + (Scaled Power 0-10)
                y = ave(ybocs, patient_id, FUN = \(x) x[1]) + 
                  ( (MeanPowerAlphaFlat_hemiMax - ave(MeanPowerAlphaFlat_hemiMax, patient_id, FUN = min)) / 
                      (ave(MeanPowerAlphaFlat_hemiMax, patient_id, FUN = max) - ave(MeanPowerAlphaFlat_hemiMax, patient_id, FUN = min)) * 10 )), 
            color = "#1CA9C9", linewidth = 1, alpha = 0.6) +
  
  geom_point(data = DF_EXAMPLES, 
             aes(x = years_after_stimon, 
                 y = ave(ybocs, patient_id, FUN = \(x) x[1]) + 
                   ( (MeanPowerAlphaFlat_hemiMax - ave(MeanPowerAlphaFlat_hemiMax, patient_id, FUN = min)) / 
                       (ave(MeanPowerAlphaFlat_hemiMax, patient_id, FUN = max) - ave(MeanPowerAlphaFlat_hemiMax, patient_id, FUN = min)) * 10 )), 
             color = "#1CA9C9", size = 4, alpha = 0.6) +
  
  
  
  geom_line(data = DF_EXAMPLES, 
            aes(x = years_after_stimon, 
                # Calculation: Start at first Y-BOCS + (Scaled Power 0-10)
                y = ave(ybocs, patient_id, FUN = \(x) x[1]) + 
                  ( (MaxPowerAlphaFlat_hemiMax - ave(MaxPowerAlphaFlat_hemiMax, patient_id, FUN = min)) / 
                      (ave(MaxPowerAlphaFlat_hemiMax, patient_id, FUN = max) - ave(MaxPowerAlphaFlat_hemiMax, patient_id, FUN = min)) * 10 )), 
            color = "#23a", linewidth = 1, alpha = 0.6) +
  
  geom_point(data = DF_EXAMPLES, 
             aes(x = years_after_stimon, 
                 y = ave(ybocs, patient_id, FUN = \(x) x[1]) + 
                   ( (MaxPowerAlphaFlat_hemiMax - ave(MaxPowerAlphaFlat_hemiMax, patient_id, FUN = min)) / 
                       (ave(MaxPowerAlphaFlat_hemiMax, patient_id, FUN = max) - ave(MaxPowerAlphaFlat_hemiMax, patient_id, FUN = min)) * 10 )), 
             color = "#23a", size = 4, alpha = 0.6) +
  
  
  # D. Main Y-BOCS Trajectory
  geom_line(data = DF_CLIN_STIM %>% filter(patient_id %in% patients_to_highlight$patient_id) %>% distinct(patient_id, session_id, .keep_all = TRUE), 
            aes(x = years_after_stimon, y = ybocs), 
            color = "black", linewidth = 1.2, alpha = 0.8) +
  
  # E. Highlights & Points
  geom_point(data = DF_CLIN_STIM %>% filter(patient_id %in% patients_to_highlight$patient_id) %>% distinct(patient_id, session_id, .keep_all = TRUE), 
             aes(x = years_after_stimon, y = ybocs), 
             shape = 21, fill = "white", color = "black", size = 4, stroke = 0.6) +
  geom_point(data = patients_to_highlight, aes(x = years_after_stimon, y = ybocs), 
             color = "black", size = 5.5) +
  
  # F. Reference Lines
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  
  # G. Final Styling
  facet_wrap(~patient_id, scales = "free_x", ncol = 4) +
  
  scale_y_continuous(
    name = "Y-BOCS Score (Black) & Scaled Alpha (Blue)",
    breaks = seq(0, 40, 10),
    limits = c(-12, 50), # Expanded slightly to accommodate the +10 swing
    sec.axis = sec_axis(~ ., name = "Neural Shift Anchored to Baseline")
  ) +
  
  scale_fill_viridis_c(option = "magma", name = "Amp [mA]", direction = -1) +
  theme_cowplot(font_size = 11) +
  theme(axis.title.y.right = element_text(size = 9, color = "gray40"),
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank(),
        legend.position = "bottom",
        coord_cartesian(clip = "off"))

print(final_faceted_plot)

save_fig(final_faceted_plot, "sub-all_ses-postop_device-all_ephysio_clinic_trajectories-examples", file.path(PATH_FIGURES,"brainsensesurvey","all"), w = 12, h = 6)

# Exemplary plots clinic + power (OP12) ---------------
OFFSET <- 10
SCALE  <- 15 # Adjust this if Power units are much smaller/larger than Y-BOCS units

exemplary_patient <- ggplot() + 
  # A. Power Data - Left Hemisphere (Solid)
  # We plot it as (Power * SCALE) + OFFSET so it sits in the Y-BOCS 25+ zone
  geom_line(data = DF_EPHYS_STIM %>% filter(patient_id == "OP12", Hemisphere == "Left", channel_id == "N3_P1"), 
            aes(x = years_after_stimon, y = (MaxPowerAlphaFlat * SCALE) + OFFSET), 
            color = "#1CA9C9", linewidth = 1, alpha = 0.6) +
  geom_point(data = DF_EPHYS_STIM %>% filter(patient_id == "OP12", Hemisphere == "Left", channel_id == "N3_P1"), 
             aes(x = years_after_stimon, y = (MaxPowerAlphaFlat * SCALE) + OFFSET), 
             color = "#1CA9C9", size = 4, alpha = 0.6) +
  
  # B. Main Y-BOCS Trajectory (Primary Axis)
  geom_line(data = DF_CLIN_STIM %>% filter(patient_id == "OP12") %>% distinct(session_id, .keep_all = TRUE), 
            aes(x = years_after_stimon, y = ybocs), 
            color = "black", linewidth = 1.2, alpha = 0.8) +
  geom_point(data = DF_CLIN_STIM %>% filter(patient_id == "OP12") %>% distinct(session_id, .keep_all = TRUE), 
             aes(x = years_after_stimon, y = ybocs), 
             shape = 21, fill = "white", color = "black", size = 4, stroke = 0.6) +
  
  # C. Reference Lines
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = OFFSET, color = "gray60", linetype = "dotted", alpha = 0.5) +
  
  # D. Dual Axis Scaling
  # The formula in sec_axis must be the INVERSE of the formula in aes(y=...)
  scale_y_continuous(
    name = "Y-BOCS Total Score",
    breaks = seq(0, 40, 10),
    limits = c(0, 42), 
    sec.axis = sec_axis(~ (. - OFFSET) / SCALE, 
                        name = "Alpha Power (Relative Shift)")
  ) +
  
  labs(x = "Time after DBS ON [years]") +
  theme_cowplot(font_size = 11) +
  theme(
    axis.title.y.right = element_text(size = 10, color = "#1CA9C9", face = "bold"),
    axis.text.y.right = element_text(color = "#1CA9C9"),
    legend.position = "bottom"
  )

print(exemplary_patient)





# (Spearman) Correlation between Y-BOCS and Ephysio changes ---------------

suffixes <- c("_hemiMax", "_hemiMean")
clinical_vars <- c("pct_improvement", "delta_improvement")

# Generate all neural delta combinations
neural_vars <- expand.grid(base = ephysio_fields, suffix = suffixes) %>%
  mutate(full = paste0(base, suffix, "_delta")) %>%
  pull(full)

# Define Hemisphere scopes
target_hemis <- list(
  "Both"  = c("Left", "Right"),
  "Left" = c("Left"),
  "Right" = c("Right")
)



# PRE-CALCULATE GLOBAL AXIS LIMITS
x_limits_map <- DF_OVERLAY_DATA %>%
  summarise(across(all_of(neural_vars), ~ max(abs(.), na.rm = TRUE) * 1.1))

y_limits_map <- DF_OVERLAY_DATA %>%
  summarise(across(all_of(clinical_vars), ~ max(abs(.), na.rm = TRUE) * 1.1))


stats_summary <- data.frame()
baseline_options <- list("No_S1" = TRUE) 


for (bl_label in names(baseline_options)) {
  exclude_s1 <- baseline_options[[bl_label]]
  
  for (h_label in names(target_hemis)) {
    for (clin_v in clinical_vars) {
      for (neur_v in neural_vars) {
        
        # Prepare Data
        df_all <- DF_OVERLAY_DATA %>%
          {if(exclude_s1) filter(., session_id != session_first, days_after_stimon > 0) else .} %>%
          filter(Hemisphere %in% target_hemis[[h_label]]) %>%
          select(patient_id, session_id, Hemisphere, x = !!sym(neur_v), y = !!sym(clin_v)) %>%
          drop_na()
        
        if (nrow(df_all) < 5) next 
        
        # Stats Logic (Spearman + Bootstrap)
        boot_res <- boot(data = df_all, statistic = spearman_boot_fn, R = nBoots)
        boot_ci  <- tryCatch(boot.ci(boot_res, type = "perc"), error = function(e) NULL)
        rho      <- boot_res$t0
        ci_low   <- if(!is.null(boot_ci)) boot_ci$percent[4] else NA
        ci_high  <- if(!is.null(boot_ci)) boot_ci$percent[5] else NA
        p_val    <- cor.test(df_all$x, df_all$y, method = "spearman")$p.value
        
        # Format p-value (Scientific if < 0.001)
        p_label <- if(p_val < 0.001) formatC(p_val, format = "e", digits = 2) else round(p_val, 3)
        stat_text <- sprintf("rho = %.2f\n95%% CI [%.2f, %.2f]\np = %s", 
                             rho, ci_low, ci_high, p_label)
        
        stats_summary <- rbind(stats_summary, data.frame(
          baseline_status = bl_label,
          hemi_scope = h_label,
          neural_metric = neur_v,
          clinical_metric = clin_v,
          rho = rho,
          ci_lower = ci_low,
          ci_upper = ci_high,
          p_value = p_val,
          n_obs = nrow(df_all)
        ))
        
        
        # Plotting Limits
        lim_x <- x_limits_map[[neur_v]]
        lim_y <- y_limits_map[[clin_v]]
        
        p <- ggplot(df_all, aes(x = x, y = y)) +
          # 1. Background Grid & Origin Lines
          geom_vline(xintercept = 0, color = "gray90", linewidth = 0.8) +
          geom_hline(yintercept = 0, color = "gray90", linewidth = 0.8) +
          
          # 2. GLOBAL LINEAR FIT
          # group = 1 ensures it pools all patients for the fit
          geom_smooth(aes(group = 1), method = "lm", color = "gray20", 
                      fill = "gray80", alpha = 0.15, linewidth = 1) +
          
          # 3. DATA POINTS
          geom_point(aes(color = patient_id, shape = Hemisphere), 
                     size = 3.5, alpha = 0.8) +
          
          # 4. STATISTICAL ANNOTATION
          annotate("text", x = -lim_x * 0.95, y = lim_y * 0.9, 
                   label = stat_text, hjust = 0, vjust = 1, 
                   size = 3.5, fontface = "italic", color = "gray30") +
          
          # 5. AESTHETICS & TURBO PALETTE (Fixes the Legend)
          scale_shape_manual(values = c("Left" = 16, "Right" = 17)) +
          scale_color_viridis_d(option = "turbo", name = "Patient ID") +
          
          coord_cartesian(xlim = c(-lim_x, lim_x), ylim = c(-lim_y, lim_y)) +
          
          labs(
            title = paste("Analysis:", gsub("_", " ", bl_label)),
            subtitle = paste("Hemi:", h_label, "| Neural Metric:", neur_v),
            x = paste("Delta", neur_v),
            y = gsub("_", " ", clin_v)
          ) +
          theme_cowplot(font_size = 12) +
          theme(
            legend.position = "right",
            plot.title = element_text(face = "bold")
          )
        
        fig_name <- sprintf("sub-all_status-%s_hemi-%s_%s_vs_%s", 
                            bl_label, h_label, neur_v, clin_v)
        save_fig(p, fig_name, file.path(PATH_FIGURES, "brainsensesurvey","all","Rstats"), w = 9, h = 9)
        
      }
    }
  }
}

# EXPORT FINAL MASTER STATS
write.csv(stats_summary, file.path(PATH_FIGURES, "brainsensesurvey","all","Rstats","sub-all_datatype-stats_comprehensive_summary.csv"), row.names = FALSE)


# (Mixed-Model interaction model ) Correlation between Y-BOCS and Ephysio taking into account confounding factors ---------------


# 1. Setup Variables
absolute_neural_vars <- expand.grid(base = ephysio_fields, suffix = suffixes) %>%
  mutate(full = paste0(base, suffix)) %>%
  pull(full)

absolute_clinical_vars <- c("ybocs")

# Pre-calculate axis limits
x_abs_limits <- DF_OVERLAY_DATA %>%
  summarise(across(all_of(absolute_neural_vars), ~ max(., na.rm = TRUE) * 1.05))

y_abs_limits <- DF_OVERLAY_DATA %>%
  summarise(across(all_of(absolute_clinical_vars), ~ max(., na.rm = TRUE) * 1.05))

# --- PART A: INTERACTION MODELS (Full Dataset) ---
# 1. Setup Variables
absolute_neural_vars <- expand.grid(base = ephysio_fields, suffix = suffixes) %>%
  mutate(full = paste0(base, suffix)) %>%
  pull(full)

absolute_clinical_vars <- c("ybocs")

# Pre-calculate axis limits for absolute values
x_abs_limits <- DF_OVERLAY_DATA %>%
  summarise(across(all_of(absolute_neural_vars), ~ max(., na.rm = TRUE) * 1.05))

y_abs_limits <- DF_OVERLAY_DATA %>%
  summarise(across(all_of(absolute_clinical_vars), ~ max(., na.rm = TRUE) * 1.05))

# -------------------------------------------------------------------------
# 2. Analysis Loop: Interaction Models (Full Dataset)
# -------------------------------------------------------------------------
abs_lmm_interaction_coefs <- data.frame()
abs_lmm_interaction_anova <- data.frame()

for (clin_v in absolute_clinical_vars) {
  for (neur_v in absolute_neural_vars) {
    
    # A. Clean Data
    df_abs <- DF_OVERLAY_DATA %>%
      select(patient_id, Hemisphere, ybocs = !!sym(clin_v), 
             NeuralVal = !!sym(neur_v), days_after_stimon, Amp) %>%
      mutate(across(c(patient_id, Hemisphere), as.factor),
             across(c(ybocs, NeuralVal, days_after_stimon, Amp), as.numeric)) %>%
      filter(!is.na(NeuralVal), !is.na(ybocs))
    
    if (nrow(df_abs) < 10) next
    
    # B. Fit Mixed Model
    # ybocs ~ Hemisphere * Neural + Time + Amp + Random Intercept
    model_abs <- lmer(ybocs ~ Hemisphere * NeuralVal + days_after_stimon + Amp + (1|patient_id), 
                      data = df_abs, REML = TRUE)
    
    # C. Extract Coefficients (t-table)
    coefs <- as.data.frame(summary(model_abs)$coefficients)
    coefs$term <- rownames(coefs)
    coefs <- coefs %>% filter(term != "(Intercept)") # <--- REMOVE INTERCEPT
    coefs$neural_metric <- neur_v
    coefs$clinical_metric <- clin_v
    abs_lmm_interaction_coefs <- rbind(abs_lmm_interaction_coefs, coefs)
    
    # D. Extract ANOVA (F-table)
    anova_res <- as.data.frame(anova(model_abs, type = "3", ddf = "Satterthwaite"))
    anova_res$term <- rownames(anova_res)
    anova_res$neural_metric <- neur_v
    anova_res$clinical_metric <- clin_v
    abs_lmm_interaction_anova <- rbind(abs_lmm_interaction_anova, anova_res)
    
    # E. Formatting for Plot
    f_val  <- anova_res["NeuralVal", "F value"]
    p_val  <- anova_res["NeuralVal", "Pr(>F)"]
    df_num <- anova_res["NeuralVal", "NumDF"]
    df_den <- anova_res["NeuralVal", "DenDF"]
    
    p_label <- if(p_val < 0.001) formatC(p_val, format = "e", digits = 2) else round(p_val, 3)
    stat_text <- sprintf("LMM Main Effect (Neural):\nF(%0.1f, %0.1f) = %0.2f\np = %s", 
                         df_num, df_den, f_val, p_label)
    
    lim_x <- x_abs_limits[[neur_v]]
    lim_y <- y_abs_limits[[clin_v]]
    
    # F. Create Plot
    p <- ggplot(df_abs, aes(x = NeuralVal, y = ybocs)) +
      geom_smooth(aes(group = Hemisphere, linetype = Hemisphere), 
                  method = "lm", color = "gray20", fill = "gray80", alpha = 0.15) +
      geom_point(aes(color = patient_id, shape = Hemisphere), size = 3.5, alpha = 0.8) +
      annotate("text", x = lim_x * 0.95, y = lim_y * 0.95, 
               label = stat_text, hjust = 1, vjust = 1, size = 3.5, fontface = "italic", color = "gray30") +
      scale_shape_manual(values = c("Left" = 16, "Right" = 17)) +
      scale_color_viridis_d(option = "turbo", name = "Patient ID") +
      coord_cartesian(xlim = c(0, lim_x), ylim = c(0, lim_y)) +
      labs(title = paste("Absolute LMM:", neur_v), x = neur_v, y = clin_v,
           subtitle = "Fixed: Hemi * Neural + Days + Amp | Random: (1|Patient)") +
      theme_cowplot(font_size = 12)
    
    save_fig(p, sprintf("sub-all_model-LMM_absolute_%s_vs_%s", neur_v, clin_v), 
             file.path(PATH_FIGURES, "brainsensesurvey","all","Rstats"), w = 9, h = 9)
  }
}

# -------------------------------------------------------------------------
# 3. Analysis Loop: Simple Models (By Hemisphere Subset)
# -------------------------------------------------------------------------
abs_lmm_simple_coefs <- data.frame()
abs_lmm_simple_anova <- data.frame()

for (h_label in names(target_hemis)) {
  for (clin_v in absolute_clinical_vars) {
    for (neur_v in absolute_neural_vars) {
      
      df_abs <- DF_OVERLAY_DATA %>%
        filter(Hemisphere %in% target_hemis[[h_label]]) %>%
        select(patient_id, ybocs = !!sym(clin_v), NeuralVal = !!sym(neur_v), days_after_stimon, Amp) %>%
        mutate(patient_id = as.factor(patient_id),
               across(c(ybocs, NeuralVal, days_after_stimon, Amp), as.numeric)) %>%
        filter(!is.na(NeuralVal), !is.na(ybocs))
      
      if (nrow(df_abs) < 10) next
      
      model_abs <- lmer(ybocs ~ NeuralVal + days_after_stimon + Amp + (1|patient_id), 
                        data = df_abs, REML = TRUE)
      
      # Extract Stats (t-table)
      coefs <- as.data.frame(summary(model_abs)$coefficients)
      coefs$term <- rownames(coefs)
      coefs <- coefs %>% filter(term != "(Intercept)") # <--- REMOVE INTERCEPT
      coefs$hemi_scope <- h_label
      coefs$neural_metric <- neur_v
      coefs$clinical_metric <- clin_v
      abs_lmm_simple_coefs <- rbind(abs_lmm_simple_coefs, coefs)
      
      # Extract ANOVA (F-table)
      anova_res <- as.data.frame(anova(model_abs, type = "3", ddf = "Satterthwaite"))
      anova_res$term <- rownames(anova_res)
      anova_res$hemi_scope <- h_label
      anova_res$neural_metric <- neur_v
      anova_res$clinical_metric <- clin_v
      abs_lmm_simple_anova <- rbind(abs_lmm_simple_anova, anova_res)
      
      # Plot Labeling
      p_label <- if(anova_res["NeuralVal", "Pr(>F)"] < 0.001) formatC(anova_res["NeuralVal", "Pr(>F)"], format = "e", digits = 2) else round(anova_res["NeuralVal", "Pr(>F)"], 3)
      stat_text <- sprintf("LMM Main Effect (Neural):\nF(%0.1f, %0.1f) = %0.2f\np = %s", 
                           anova_res["NeuralVal", "NumDF"], anova_res["NeuralVal", "DenDF"], anova_res["NeuralVal", "F value"], p_label)
      
      # Plotting
      lim_x <- x_abs_limits[[neur_v]]
      lim_y <- y_abs_limits[[clin_v]]
      
      p <- ggplot(df_abs, aes(x = NeuralVal, y = ybocs)) +
        geom_smooth(method = "lm", color = "gray20", fill = "gray80", alpha = 0.15) +
        geom_point(aes(color = patient_id), size = 3.5, alpha = 0.8) +
        annotate("text", x = lim_x * 0.95, y = lim_y * 0.95, label = stat_text, 
                 hjust = 1, vjust = 1, size = 3.5, fontface = "italic", color = "gray30") +
        scale_color_viridis_d(option = "turbo", name = "Patient ID") +
        coord_cartesian(xlim = c(0, lim_x), ylim = c(0, lim_y)) +
        labs(title = paste("Absolute LMM:", neur_v, "| Hemi:", h_label), x = neur_v, y = clin_v) +
        theme_cowplot(font_size = 12)
      
      save_fig(p, sprintf("sub-all_model-LMM-simple_absolute_hemi-%s_%s_vs_%s", h_label, neur_v, clin_v), 
               file.path(PATH_FIGURES, "brainsensesurvey","all","Rstats"), w = 9, h = 9)
    }
  }
}

# -------------------------------------------------------------------------
# 4. Final Exports
# -------------------------------------------------------------------------

# Interaction Exports
write.csv(abs_lmm_interaction_coefs, file.path(PATH_OUTPUT, "sub-all_LMM_interaction_COEFFICIENTS.csv"), row.names = FALSE)
write.csv(abs_lmm_interaction_anova, file.path(PATH_OUTPUT, "sub-all_LMM_interaction_ANOVA.csv"), row.names = FALSE)

# Simple Subset Exports
write.csv(abs_lmm_simple_coefs, file.path(PATH_OUTPUT, "sub-all_LMM_simple_COEFFICIENTS.csv"), row.names = FALSE)
write.csv(abs_lmm_simple_anova, file.path(PATH_OUTPUT, "sub-all_LMM_simple_ANOVA.csv"), row.names = FALSE)



# Plot exemplary abseline psds ------

#patients_to_plot <- c("OP05","OP06")
# 1. Summarise and calculate YBOCS per patient
DF_EPHYS_EXEMPLARY <- DF_EPHYS_LONG %>%
  filter(session_id == 1, patient_id %in% c("OP05","OP06","OP09","OP12")) %>%
  mutate(FrequencyFlat = as.numeric(FrequencyFlat)) %>%
  group_by(patient_id, Hemisphere, FrequencyFlat) %>%
  summarise(
    n_obs     = n(),
    avg_power = mean(.data$MeanPowerFlat, na.rm = TRUE),
    sd_power  = sd(.data$MeanPowerFlat, na.rm = TRUE),
    se_power  = sd_power / sqrt(n_obs),
    ci_lower  = avg_power - (1.96 * se_power),
    ci_upper  = avg_power + (1.96 * se_power),
    # Keep YBOCS here for the sorting/labeling
    ybocs     = first(ybocs), 
    .groups   = "drop"
  ) %>%
  # Apply sorting: Lowest YBOCS at top (desc = FALSE) or Highest at top (desc = TRUE)
  mutate(patient_id = fct_reorder(patient_id, ybocs, .desc = FALSE))

# 2. Extract labels for the inner-plot text
DF_LABELS <- DF_EPHYS_EXEMPLARY %>%
  group_by(patient_id, Hemisphere) %>%
  summarise(
    ybocs = first(ybocs),
    x_pos = 28, # Positioned near the right (30Hz mark)
    y_pos = max(ci_upper, na.rm = TRUE) * 0.9,
    .groups = "drop"
  )

# 3. Build the Plot
final_ephys_plot_clean <- ggplot(DF_EPHYS_EXEMPLARY, aes(x = FrequencyFlat, y = avg_power)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Hemisphere), alpha = 0.2) +
  geom_line(aes(color = Hemisphere), linewidth = 0.7) +
  
  # YBOCS text INSIDE the plot area
  geom_text(data = DF_LABELS, 
            aes(x = x_pos, y = y_pos, label = paste0("YBOCS: ", round(ybocs, 1))),
            size = 3, family = FONT_FAMILY, fontface = "italic", hjust = 1) +
  
  facet_grid(patient_id ~ Hemisphere, scales = "fixed") + 
  
  scale_x_continuous(
    name = "Frequency (Hz)",
    trans = "log10",
    limits = c(2, 35), 
    breaks = c(4, 8, 12, 30),
    labels = c("4", "8", "12", "30"),
    expand = expansion(mult = c(0, 0.02)) 
  ) +
  
  scale_y_continuous(
    name = "periodic PSD [dB]", 
    expand = expansion(mult = c(0.1, 0.2)) 
  ) +
  
  scale_color_manual(values = c("Left" = "#3182bd", "Right" = "#de2d26")) +
  scale_fill_manual(values = c("Left" = "#3182bd", "Right" = "#de2d26")) +
  theme_cowplot(font_size = 10, font_family = FONT_FAMILY) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )

print(final_ephys_plot_clean)

save_fig(final_ephys_plot_clean, "sub=all_ses-postop-first_bipolar-flat-spectral_hemi-mean_individual", 
         file.path(PATH_FIGURES, "brainsensesurvey","first"), w = 4, h = 8)



DF_EPHYS_EXEMPLARY <- DF_EPHYS_LONG %>%
  filter(session_id == 1, patient_id %in% c("OP05","OP06","OP09","OP12")) %>%
  mutate(FrequencyFlat = as.numeric(Frequency)) %>%
  group_by(patient_id, Hemisphere, Frequency) %>%
  summarise(
    n_obs     = n(),
    avg_power = mean(.data$MeanPower, na.rm = TRUE),
    sd_power  = sd(.data$MeanPower, na.rm = TRUE),
    se_power  = sd_power / sqrt(n_obs),
    ci_lower  = avg_power - (1.96 * se_power),
    ci_upper  = avg_power + (1.96 * se_power),
    # Keep YBOCS here for the sorting/labeling
    ybocs     = first(ybocs), 
    .groups   = "drop"
  ) %>%
  # Apply sorting: Lowest YBOCS at top (desc = FALSE) or Highest at top (desc = TRUE)
  mutate(patient_id = fct_reorder(patient_id, ybocs, .desc = FALSE))



# 3. Build the Plot
final_ephys_plot_clean <- ggplot(DF_EPHYS_EXEMPLARY, aes(x = Frequency, y = avg_power)) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = Hemisphere), alpha = 0.2) +
  geom_line(aes(color = Hemisphere), linewidth = 0.7) +
  
  # YBOCS text INSIDE the plot area
  geom_text(data = DF_LABELS, 
            aes(x = x_pos, y = y_pos, label = paste0("YBOCS: ", round(ybocs, 1))),
            size = 3, family = FONT_FAMILY, fontface = "italic", hjust = 1) +
  
  facet_grid(patient_id ~ Hemisphere, scales = "free_y") + 
  
  scale_x_continuous(
    name = "Frequency (Hz)",
    trans = "log10",
    limits = c(2, 35), 
    breaks = c(4, 8, 12, 30),
    labels = c("4", "8", "12", "30"),
    expand = expansion(mult = c(0, 0.02)) 
  ) +
  
  scale_y_log10(
    name = "periodic PSD [dB]", 
    expand = expansion(mult = c(0.1, 0.2)) 
  ) +
  
  scale_color_manual(values = c("Left" = "#3182bd", "Right" = "#de2d26")) +
  scale_fill_manual(values = c("Left" = "#3182bd", "Right" = "#de2d26")) +
  theme_cowplot(font_size = 10, font_family = FONT_FAMILY) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )

print(final_ephys_plot_clean)

save_fig(final_ephys_plot_clean, "sub=all_ses-postop-first_bipolar-absolute-spectral_hemi-mean_individual", 
         file.path(PATH_FIGURES, "brainsensesurvey","first"), w = 4, h = 8)


# Assuming 'neural' is your MeanPowerAlphaFlat_hemiMax variable
model <- lmer(ybocs ~ MaxPowerAlphaFlat_hemiMax + days_after_stimon + (1|patient_id), 
              data = DF_OVERLAY_DATA)


# ------- plot mixed model ----


# Assuming 'model' is: lmer(ybocs ~ Hemisphere * MeanPowerAlphaFlat_hemiMax + stim_amp + days_after_stimon + (1|patient_id))

# 1. Fit the model
# Assuming 'neural' is your MeanPowerAlphaFlat_hemiMax variable
model <- lmer(ybocs ~ Hemisphere * MaxPowerAlphaFlat_hemiMax + Amp + days_after_stimon + (1|patient_id), 
              data = DF_OVERLAY_DATA)

# 2. Extract Coefficients for the LEFT Hemisphere (The reference level)
intercept_lmm <- fixef(model)["(Intercept)"]
slope_left    <- fixef(model)["MaxPowerAlphaFlat_hemiMax"]
# Extract Fixed Effects
cf <- fixef(model)

# Specifically for the 'Left' (Reference) effect:
# intercept_lmm = Intercept + (Mean Stim Amp * Beta_Stim) + (Mean Days * Beta_Days)
# slope_lmm     = MeanPowerAlphaFlat_hemiMax (The base slope)

avg_amp  <- mean(DF_OVERLAY_DATA$Amp, na.rm = TRUE)
avg_days <- mean(DF_OVERLAY_DATA$days_after_stimon, na.rm = TRUE)

lmm_intercept <- cf["(Intercept)"] + (cf["Amp"] * avg_amp) + (cf["days_after_stimon"] * avg_days)
lmm_slope     <- cf["MaxPowerAlphaFlat_hemiMax"]


# 1. Create the Forest Plot
# 'prob.inner' and 'prob.outer' create the shaded CI regions
plot_model(model, 
           type = "est", 
           show.values = TRUE, 
           value.offset = .3,
           vline.color = "gray90",
           title = "Fixed Effects on Y-BOCS Total Score",
           axis.labels = c(
             "HemisphereRight:MaxPowerAlphaFlat_hemiMax" = "Alpha * Hemisphere (Interaction)",
             "days_after_stimon" = "Time (Days Post-Op)",
             "stim_amp" = "Stimulation Amplitude (mA)",
             "MaxPowerAlphaFlat_hemiMax" = "Alpha Power (Left/Ref)",
             "HemisphereRight" = "Hemisphere (Right)"
           )) +
  COW_THEME +
  scale_color_manual(values = c("pos" = "#1CA9C9", "neg" = "#de2d26")) # Blue for increase, Red for decrease







# 1. Extract ANOVA (Type III) for F-stats
# This handles the F-significance logic
anova_df <- as.data.frame(anova(model, type = 3)) %>%
  tibble::rownames_to_column("anova_term") %>%
  rename(F_stat = `F value`, p_f = `Pr(>F)`)

# 2. Extract Fixed Effects (Global Model)
df_fixef <- tidy(model, conf.int = TRUE) %>%
  filter(effect == "fixed", term != "(Intercept)") %>%
  mutate(
    # Mapping model terms to ANOVA rows for the F-stat join
    anova_term = case_when(
      term == "MaxPowerAlphaFlat_hemiMax" ~ "MaxPowerAlphaFlat_hemiMax",
      term == "Amp"                        ~ "Amp",
      term == "days_after_stimon"           ~ "days_after_stimon",
      term == "HemisphereRight"            ~ "Hemisphere",
      grepl(":", term)                      ~ "Hemisphere:MaxPowerAlphaFlat_hemiMax",
      TRUE ~ term
    )
  ) %>%
  left_join(anova_df, by = "anova_term") %>%
  mutate(
    term_label = case_when(
      term == "MaxPowerAlphaFlat_hemiMax" ~ "Alpha Power (Left/Ref)",
      term == "Amp"                        ~ "Stimulation Amplitude",
      term == "days_after_stimon"           ~ "Time (Days Post-Op)",
      term == "HemisphereRight"            ~ "Hemisphere (Right Offset)",
      grepl(":", term)                      ~ "Interaction: Alpha × Hemi",
      TRUE ~ term
    ),
    Source = "Model Terms (Fixed Effects)",
    # Significance based on F-test p-value (p_f)
    is_sig = p_f < 0.05,
    stars = case_when(p_f < 0.001 ~ "***", p_f < 0.01 ~ "**", p_f < 0.05 ~ "*", TRUE ~ "")
  )

# 3. Extract Simple Slopes (No F-stat for derived slopes, so we use t-stat p-values)
slopes_res <- emtrends(model, ~ Hemisphere, var = "MaxPowerAlphaFlat_hemiMax")
df_slopes <- as.data.frame(summary(slopes_res, infer = TRUE)) %>%
  mutate(
    term_label = paste("Slope:", Hemisphere),
    Source = "Hemisphere-Specific Results",
    is_sig = p.value < 0.05,
    stars = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""),
    F_stat = NA  # ANOVA F doesn't apply to individual slope extractions
  ) %>%
  rename(estimate = 2, conf.low = lower.CL, conf.high = upper.CL)

# 4. Final Merge
df_final <- bind_rows(df_fixef, df_slopes) %>%
  mutate(term_label = factor(term_label, levels = rev(unique(term_label))))

# 5. The High-Impact Plot
plot_LMEstats<- ggplot(df_final, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  
  # Error bars
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = is_sig), 
                 height = 0, linewidth = 1.2) +
  
  # Points
  geom_point(aes(fill = estimate < 0, shape = Source), size = 5, color = "white", stroke = 1) +
  
  # LABEL: Adding F-statistic and Stars
  geom_text(aes(label = ifelse(!is.na(F_stat), 
                               paste0("F=", round(F_stat, 1), stars), 
                               stars)), 
            vjust = -1.2, size = 3.5, fontface = "bold", family = "sans") +
  
  # Aesthetic tweaks
  scale_shape_manual(values = c("Model Terms (Fixed Effects)" = 21, "Hemisphere-Specific Results" = 23)) +
  scale_fill_manual(values = c("TRUE" = "#e74c3c", "FALSE" = "#3498db"), 
                    labels = c("Symptom Increase", "Clinical Improvement"), name = "Effect") +
  scale_color_manual(values = c("TRUE" = "#2c3e50", "FALSE" = "gray80"), guide = "none") +
  
  labs(
    title = "Statistical Drivers of Y-BOCS Recovery",
    subtitle = "Significance and Stars based on ANOVA F-tests. Red = Improvement.",
    x = "Beta Coefficient (Estimate)", y = NULL
  ) +
  COW_THEME +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(face = "bold", color = "#2c3e50"),
    legend.position = "bottom"
  )

save_fig(plot_LMEstats, paste0("sub-all_LMM-","MaxPowerAlphaFlat_hemiMax","_interaction_forest-plot"), 
         file.path(PATH_FIGURES, "brainsensesurvey","all"), w = 4, h = 3)

