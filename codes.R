# Load packages ----------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)
library(lmerTest)
library(car)

options(contrasts = c("contr.sum", "contr.poly" ))


# Plot settings ----------------------------------------------------------

custom_colors <- c(
  "A"    = "red4",
  "B"    = "midnightblue",
  "C"    = "darkorange2",
  "D"    = "goldenrod",
  "E"    = "darkgreen",
  "WT_1" = "grey70",
  "WT_2" = "grey34",
  "SW"   = "black"
)

custom_shapes <- c(
  "A"    = 1,
  "B"    = 2,
  "C"    = 3,
  "D"    = 4,
  "E"    = 5,
  "WT_1" = 6,
  "WT_2" = 7,
  "SW"   = 8
)


# Data import ------------------------------------------------------------

# Select larvae.csv
larval_pupal_data <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
) %>%
  mutate(
    across(
      c(Genetic_line, Familiar_diet, Treatment, Replicate, Container_ID),
      factor
    ),
    Familiar_diet = factor(Familiar_diet, levels = c("CF", "FW"))
  )


# Select survival.csv
adult_survival_data <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(Date.of.Death), !is.na(Sex), Date.of.Death != "") %>%
  mutate(
    Date.of.Death    = as.Date(Date.of.Death, "%d/%m/%Y"),
    Date.of.Eclosion = as.Date(Date.of.Eclosion, "%d/%m/%Y"),
    lifespan         = as.numeric(Date.of.Death - Date.of.Eclosion)
  ) %>%
  filter(lifespan != 0) %>%
  mutate(
    across(
      c(ID, Container_ID, Sex, Genetic_line, Treatment, Familiar_diet, Replicate, Group),
      factor
    ),
    Familiar_diet = factor(Familiar_diet, levels = c("CF", "FW"))
  )


# Select repro_morphs.csv
adult_repro_data <- read.csv(
  file.choose(),
  stringsAsFactors = FALSE
) %>%
  mutate(
    across(
      c(ID, Container_ID, Genetic_line, Familiar_diet, Treatment, Sex, Replicate, Operator),
      factor
    ),
    Familiar_diet = factor(Familiar_diet, levels = c("CF", "FW")),
    
    ow = as.numeric(as.character(ow)),
    tw = as.numeric(as.character(tw)),
    ag = as.numeric(as.character(ag)),
    
    GSI = case_when(
      as.character(GSI) %in% c("", "NA", "N/A", "n/a", "-", "#DIV/0!", "Inf", "NaN") ~ NA_character_,
      TRUE ~ as.character(GSI)
    ),
    ag_gonad_ratio = case_when(
      as.character(ag_gonad_ratio) %in% c("", "NA", "N/A", "n/a", "-", "#DIV/0!", "Inf", "NaN") ~ NA_character_,
      TRUE ~ as.character(ag_gonad_ratio)
    ),
    
    GSI = as.numeric(GSI),
    ag_gonad_ratio = as.numeric(ag_gonad_ratio)
  )


# Models -----------------------------------------------------------------

data_female <- adult_repro_data %>% filter(Sex == "F")
data_male   <- adult_repro_data %>% filter(Sex == "M")

model_specs <- list(
  
  "Specific growth rate" = list(
    data = larval_pupal_data,
    model_1 = SGR_per_day ~ Genetic_line * Treatment,
    model_2 = SGR_per_day ~ Familiar_diet * Treatment + (1 | Genetic_line)
  ),
  
  "Prepupal weight" = list(
    data = larval_pupal_data,
    model_1 = avg_pupal_weight ~ Genetic_line * Treatment,
    model_2 = avg_pupal_weight ~ Familiar_diet * Treatment + (1 | Genetic_line)
  ),
  
  "Larval duration" = list(
    data = larval_pupal_data,
    model_1 = larval_duration ~ Genetic_line * Treatment,
    model_2 = larval_duration ~ Familiar_diet * Treatment + (1 | Genetic_line)
  ),
  
  "Pupation rate" = list(
    data = larval_pupal_data,
    model_1 = pupation_rate ~ Genetic_line * Treatment,
    model_2 = pupation_rate ~ Familiar_diet * Treatment + (1 | Genetic_line)
  ),
  
  "Pupal duration" = list(
    data = larval_pupal_data,
    model_1 = pupal_duration ~ Genetic_line * Treatment,
    model_2 = pupal_duration ~ Familiar_diet * Treatment + (1 | Genetic_line)
  ),
  
  "Eclosion rate" = list(
    data = larval_pupal_data,
    model_1 = eclosion_rate ~ Genetic_line * Treatment,
    model_2 = eclosion_rate ~ Familiar_diet * Treatment + (1 | Genetic_line)
  ),
  
  "Adult body weight" = list(
    data = adult_repro_data,
    model_1 = Body.weight ~ Genetic_line * Treatment * Sex +
      (1 | Container_ID) +
      (1 | Operator),
    model_2 = Body.weight ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID) +
      (1 | Operator)
  ),
  
  "Adult lifespan" = list(
    data = adult_survival_data,
    model_1 = lifespan ~ Genetic_line * Treatment * Sex +
      (1 | Container_ID),
    model_2 = lifespan ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Ovaries" = list(
    data = data_female,
    model_1 = ow ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = ow ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Testes" = list(
    data = data_male,
    model_1 = tw ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = tw ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Female accessory gland" = list(
    data = data_female,
    model_1 = ag ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = ag ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Male accessory gland" = list(
    data = data_male,
    model_1 = ag ~ Genetic_line * Treatment +
      (1 | Container_ID) +
      (1 | Operator),
    model_2 = ag ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID) +
      (1 | Operator)
  ),
  
  "Female GSI" = list(
    data = data_female,
    model_1 = GSI ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = GSI ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Male GSI" = list(
    data = data_male,
    model_1 = GSI ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = GSI ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Female AG to gonad ratio" = list(
    data = data_female,
    model_1 = ag_gonad_ratio ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = ag_gonad_ratio ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  ),
  
  "Male AG to gonad ratio" = list(
    data = data_male,
    model_1 = ag_gonad_ratio ~ Genetic_line * Treatment +
      (1 | Container_ID),
    model_2 = ag_gonad_ratio ~ Familiar_diet * Treatment +
      (1 | Genetic_line) +
      (1 | Container_ID)
  )
)


# Fit models -------------------------------------------------------------

fit_model <- function(formula, data) {
  if (length(lme4::findbars(formula)) == 0) {
    lm(formula, data = data)
  } else {
    lmer(formula, data = data)
  }
}

is_mixed_model <- function(model) {
  inherits(model, "merMod")
}

anova_type3 <- function(model) {
  if (is_mixed_model(model)) {
    anova(model, type = 3)
  } else {
    car::Anova(model, type = 3)
  }
}

models_1 <- purrr::map(
  model_specs,
  ~ fit_model(.x$model_1, .x$data)
)

models_2 <- purrr::map(
  model_specs,
  ~ fit_model(.x$model_2, .x$data)
)

anova_model_1 <- purrr::map(
  models_1,
  anova_type3
)

anova_model_2 <- purrr::map(
  models_2,
  anova_type3
)


# Model diagnostics ------------------------------------------------------

plot_model_diagnostics <- function(model, model_name) {
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow = c(1, 2))
  
  plot(
    fitted(model),
    resid(model),
    main = paste(model_name, "- residuals vs fitted"),
    xlab = "Fitted values",
    ylab = "Residuals"
  )
  abline(h = 0, lty = 2)
  
  qqnorm(
    resid(model),
    main = paste(model_name, "- QQ plot")
  )
  qqline(resid(model))
}

purrr::iwalk(
  models_1,
  ~ plot_model_diagnostics(.x, paste("Model 1:", .y))
)

purrr::iwalk(
  models_2,
  ~ plot_model_diagnostics(.x, paste("Model 2:", .y))
)


# Check singular fits and random-effect variance -------------------------

check_singular <- function(model) {
  if (is_mixed_model(model)) {
    lme4::isSingular(model)
  } else {
    NA
  }
}

get_random_effects <- function(model) {
  if (is_mixed_model(model)) {
    VarCorr(model)
  } else {
    NA
  }
}

singular_model_1 <- purrr::map_lgl(
  models_1,
  check_singular
)

singular_model_2 <- purrr::map_lgl(
  models_2,
  check_singular
)

random_effects_model_1 <- purrr::map(
  models_1,
  get_random_effects
)

random_effects_model_2 <- purrr::map(
  models_2,
  get_random_effects
)

singular_model_1
singular_model_2

random_effects_model_1
random_effects_model_2


# P-value helper functions -----------------------------------------------

p_to_stars <- function(p) {
  if (is.na(p)) {
    return(" n.s.")
  } else if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(" n.s.")
  }
}

format_p <- function(p) {
  if (is.na(p)) {
    return(NA_character_)
  } else if (p < 0.001) {
    return("<0.001***")
  } else {
    return(paste0(format(round(p, 3), nsmall = 3), p_to_stars(p)))
  }
}

extract_p <- function(anova_result, term_name) {
  
  anova_df <- as.data.frame(anova_result)
  anova_df$Term <- rownames(anova_df)
  
  p_col <- grep("Pr", names(anova_df), value = TRUE)
  
  if (!term_name %in% anova_df$Term) {
    return(NA_real_)
  }
  
  anova_df %>%
    filter(Term == term_name) %>%
    pull(all_of(p_col)) %>%
    as.numeric()
}


# Model 2 P-value table --------------------------------------------------

model2_table <- tibble::tibble(
  Response_variable = names(anova_model_2),
  `Dietary background` = purrr::map_chr(
    anova_model_2,
    ~ format_p(extract_p(.x, "Familiar_diet"))
  ),
  `Experimental diet` = purrr::map_chr(
    anova_model_2,
    ~ format_p(extract_p(.x, "Treatment"))
  ),
  `Dietary background × Experimental diet` = purrr::map_chr(
    anova_model_2,
    ~ format_p(extract_p(.x, "Familiar_diet:Treatment"))
  )
)


View(model2_table)



# Plotting functions -----------------------------------------------------

summarize_reaction <- function(df, traits, include_sex = FALSE, sex_filter = NULL) {
  
  if (!is.null(sex_filter)) {
    df <- df %>% filter(Sex == sex_filter)
  }
  
  grouping_vars <- c("Trait", "Treatment", "Genetic_line", "Familiar_diet")
  
  if (include_sex) {
    grouping_vars <- c(grouping_vars, "Sex")
  }
  
  df %>%
    select(
      Genetic_line,
      Familiar_diet,
      Treatment,
      any_of("Sex"),
      all_of(traits)
    ) %>%
    pivot_longer(
      cols = all_of(traits),
      names_to = "Trait",
      values_to = "Value"
    ) %>%
    group_by(across(all_of(grouping_vars))) %>%
    summarise(
      Mean = mean(Value, na.rm = TRUE),
      SD   = sd(Value, na.rm = TRUE),
      N    = sum(!is.na(Value)),
      SE   = SD / sqrt(N),
      CI_lower = Mean - qt(0.975, df = N - 1) * SE,
      CI_upper = Mean + qt(0.975, df = N - 1) * SE,
      .groups = "drop"
    )
}


plot_reaction <- function(df, trait_name, y_lab, facet_sex = FALSE,
                          title = NULL, remove_legend = FALSE,
                          remove_x_label = FALSE) {
  
  dodge_width <- 0.2
  plot_title <- if (is.null(title) || is.na(title)) NULL else title
  
  p <- ggplot(
    filter(df, Trait == trait_name),
    aes(
      x = Treatment,
      y = Mean,
      colour = Genetic_line,
      shape = Genetic_line,
      linetype = Familiar_diet,
      group = Genetic_line
    )
  ) +
    geom_line(
      position = position_dodge(width = dodge_width),
      linewidth = 0.9
    ) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper),
      position = position_dodge(width = dodge_width),
      width = 0.25,
      linewidth = 0.6
    ) +
    geom_point(
      position = position_dodge(width = dodge_width),
      size = 2,
      stroke = 1
    ) +
    scale_color_manual(name = "Genetic line", values = custom_colors) +
    scale_shape_manual(name = "Genetic line", values = custom_shapes) +
    scale_linetype_manual(
      name = "Dietary background",
      values = c("CF" = "solid", "FW" = "dashed")
    ) +
    scale_x_discrete(
      expand = expansion(add = c(0.20, 0.20))
    ) +
    labs(
      x = if (remove_x_label) NULL else "Experimental diet",
      y = y_lab,
      title = plot_title
    ) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, color = "black", size = 10),
      axis.text.y     = element_text(color = "black", size = 10),
      axis.title      = element_text(size = 12, color = "black"),
      panel.border    = element_rect(color = "black", fill = NA, linewidth = 1),
      panel.grid      = element_blank(),
      strip.text      = element_text(size = 12, face = "bold"),
      plot.title      = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.title    = element_text(size = 12),
      legend.text     = element_text(size = 10),
      legend.position = if (remove_legend) "none" else "right"
    )
  
  if (facet_sex) {
    p <- p + facet_wrap(~ Sex, ncol = 2)
  }
  
  p
}


# Model 1 labels for plots -----------------------------------------------

make_model1_label_basic <- function(trait_name, show_ns = TRUE) {
  
  anova_result <- anova_model_1[[trait_name]]
  
  label_df <- tibble::tibble(
    short_term = c("G", "E", "G×E"),
    model_term = c(
      "Genetic_line",
      "Treatment",
      "Genetic_line:Treatment"
    )
  ) %>%
    mutate(
      p_value = purrr::map_dbl(model_term, ~ extract_p(anova_result, .x)),
      stars   = purrr::map_chr(p_value, p_to_stars),
      label   = paste0(short_term, stars)
    )
  
  if (!show_ns) {
    label_df <- label_df %>%
      filter(stars != " n.s.")
  }
  
  if (nrow(label_df) == 0) {
    return("n.s.")
  }
  
  paste(label_df$label, collapse = "\n")
}


make_model1_label_sex <- function(trait_name, show_ns = TRUE) {
  
  anova_result <- anova_model_1[[trait_name]]
  
  label_df <- tibble::tibble(
    short_term = c("G", "E", "S", "G×E", "G×S", "E×S", "G×E×S"),
    model_term = c(
      "Genetic_line",
      "Treatment",
      "Sex",
      "Genetic_line:Treatment",
      "Genetic_line:Sex",
      "Treatment:Sex",
      "Genetic_line:Treatment:Sex"
    )
  ) %>%
    mutate(
      p_value = purrr::map_dbl(model_term, ~ extract_p(anova_result, .x)),
      stars   = purrr::map_chr(p_value, p_to_stars),
      label   = paste0(short_term, stars)
    )
  
  if (!show_ns) {
    label_df <- label_df %>%
      filter(stars != " n.s.")
  }
  
  if (nrow(label_df) == 0) {
    return("n.s.")
  }
  
  paste(label_df$label, collapse = "\n")
}


add_model_label_basic <- function(p, label_text) {
  
  p +
    annotate(
      "text",
      x = -Inf,
      y = Inf,
      label = label_text,
      hjust = -0.08,
      vjust = 1.08,
      size = 3.5,
      fontface = "bold",
      family = "mono",
      lineheight = 0.9
    )
}


add_model_label_female_facet <- function(p, label_text) {
  
  label_df <- data.frame(
    Sex = factor("F"),
    label = label_text
  )
  
  p +
    geom_text(
      data = label_df,
      aes(
        x = -Inf,
        y = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = -0.08,
      vjust = 1.08,
      size = 3.5,
      fontface = "bold",
      family = "mono",
      lineheight = 0.9
    )
}


# Plot summaries ---------------------------------------------------------

summary_larvae <- summarize_reaction(
  larval_pupal_data,
  c("SGR_per_day", "avg_pupal_weight", "larval_duration")
)

summary_pupal <- summarize_reaction(
  larval_pupal_data,
  c("pupation_rate", "pupal_duration", "eclosion_rate")
)

summary_bw <- summarize_reaction(
  adult_repro_data,
  "Body.weight",
  include_sex = TRUE
)

summary_lifespan <- summarize_reaction(
  adult_survival_data,
  "lifespan",
  include_sex = TRUE
)

summary_ow_f <- summarize_reaction(
  adult_repro_data,
  "ow",
  sex_filter = "F"
)

summary_tw_m <- summarize_reaction(
  adult_repro_data,
  "tw",
  sex_filter = "M"
)

summary_ag_f <- summarize_reaction(
  adult_repro_data,
  "ag",
  sex_filter = "F"
)

summary_ag_m <- summarize_reaction(
  adult_repro_data,
  "ag",
  sex_filter = "M"
)

summary_GSI_f <- summarize_reaction(
  adult_repro_data,
  "GSI",
  sex_filter = "F"
)

summary_GSI_m <- summarize_reaction(
  adult_repro_data,
  "GSI",
  sex_filter = "M"
)

summary_ag_ratio_f <- summarize_reaction(
  adult_repro_data,
  "ag_gonad_ratio",
  sex_filter = "F"
)

summary_ag_ratio_m <- summarize_reaction(
  adult_repro_data,
  "ag_gonad_ratio",
  sex_filter = "M"
)


# Generate annotated plots -----------------------------------------------

plot_specs <- tibble::tibble(
  plot_name = c(
    "p_SGR",
    "p_pupal_weight",
    "p_larval_duration",
    "p_pupation_rate",
    "p_pupal_duration",
    "p_eclosion_rate",
    "p_bw",
    "p_lifespan",
    "p_ow_f",
    "p_tw_m",
    "p_ag_f",
    "p_ag_m",
    "p_GSI_f",
    "p_GSI_m",
    "p_ag_ratio_f",
    "p_ag_ratio_m"
  ),
  df = list(
    summary_larvae,
    summary_larvae,
    summary_larvae,
    summary_pupal,
    summary_pupal,
    summary_pupal,
    summary_bw,
    summary_lifespan,
    summary_ow_f,
    summary_tw_m,
    summary_ag_f,
    summary_ag_m,
    summary_GSI_f,
    summary_GSI_m,
    summary_ag_ratio_f,
    summary_ag_ratio_m
  ),
  trait_name = c(
    "SGR_per_day",
    "avg_pupal_weight",
    "larval_duration",
    "pupation_rate",
    "pupal_duration",
    "eclosion_rate",
    "Body.weight",
    "lifespan",
    "ow",
    "tw",
    "ag",
    "ag",
    "GSI",
    "GSI",
    "ag_gonad_ratio",
    "ag_gonad_ratio"
  ),
  y_lab = c(
    "Specific growth rate (per day)",
    "Mean prepupal weight (g)",
    "Larval duration (days)",
    "Pupation success (%)",
    "Pupal duration (days)",
    "Eclosion success (%)",
    "Body weight (g)",
    "Lifespan (days)",
    "Ovaries (mg)",
    "Testes (mg)",
    "Female accessory gland (mg)",
    "Male accessory gland (mg)",
    "Female GSI",
    "Male GSI",
    "Female AG to gonad ratio",
    "Male AG to gonad ratio"
  ),
  facet_sex = c(
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE,
    TRUE, TRUE,
    FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE
  ),
  title = c(
    NA, NA, NA,
    NA, NA, NA,
    NA, NA,
    "Female", "Male", NA, NA,
    "Female", "Male", "Female", "Male"
  ),
  remove_legend = c(
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE,
    FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE, FALSE
  ),
  remove_x_label = c(
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE,
    FALSE, FALSE,
    TRUE, TRUE, FALSE, FALSE,
    TRUE, TRUE, FALSE, FALSE
  ),
  model_name = c(
    "Specific growth rate",
    "Prepupal weight",
    "Larval duration",
    "Pupation rate",
    "Pupal duration",
    "Eclosion rate",
    "Adult body weight",
    "Adult lifespan",
    "Ovaries",
    "Testes",
    "Female accessory gland",
    "Male accessory gland",
    "Female GSI",
    "Male GSI",
    "Female AG to gonad ratio",
    "Male AG to gonad ratio"
  ),
  label_type = c(
    "basic", "basic", "basic",
    "basic", "basic", "basic",
    "sex", "sex",
    "basic", "basic", "basic", "basic",
    "basic", "basic", "basic", "basic"
  )
)

plots <- purrr::pmap(
  list(
    plot_specs$df,
    plot_specs$trait_name,
    plot_specs$y_lab,
    plot_specs$facet_sex,
    plot_specs$title,
    plot_specs$remove_legend,
    plot_specs$remove_x_label,
    plot_specs$model_name,
    plot_specs$label_type
  ),
  function(df, trait_name, y_lab, facet_sex, title,
           remove_legend, remove_x_label, model_name, label_type) {
    
    p <- plot_reaction(
      df = df,
      trait_name = trait_name,
      y_lab = y_lab,
      facet_sex = facet_sex,
      title = title,
      remove_legend = remove_legend,
      remove_x_label = remove_x_label
    )
    
    label_text <- if (label_type == "sex") {
      make_model1_label_sex(model_name)
    } else {
      make_model1_label_basic(model_name)
    }
    
    if (label_type == "sex") {
      add_model_label_female_facet(p, label_text)
    } else {
      add_model_label_basic(p, label_text)
    }
  }
)

names(plots) <- plot_specs$plot_name
list2env(plots, envir = .GlobalEnv)


# View annotated plots ---------------------------------------------------

p_SGR
p_pupal_weight
p_larval_duration

p_pupation_rate
p_pupal_duration
p_eclosion_rate

p_bw
p_lifespan

p_ow_f
p_tw_m
p_ag_f
p_ag_m

p_GSI_f
p_GSI_m
p_ag_ratio_f
p_ag_ratio_m

### Figure 1 - Growth rate ###

# Read the data(select: Growth rate data.csv)
data_growth <- read.csv(
  file.choose()
)

# Check required columns
if (!all(c("line", "diet", "day", "weight_g") %in% colnames(data_growth))) {
  stop("The dataset must contain: 'line', 'diet', 'day', 'weight_g'.")
}

# Exclude missing
if (any(is.na(data_growth))) {
  warning("Missing values detected; they will be removed before summary.")
}

# ----- Set ORDER of facet panels -----
# Put SW and WT_2 in the last row
data_growth$line <- factor(
  data_growth$line,
  levels = c("A","B","C","D","E","WT_1","WT_2","SW")
)

# Summaries
summary_data_growth <- data_growth %>% 
  group_by(line, diet, day) %>% 
  summarise(
    mean_biomass = mean(weight_g, na.rm = TRUE), 
    se_biomass   = sd(weight_g, na.rm = TRUE) / sqrt(n()),
    .groups = "drop_last"
  )

# Plot
ggplot(summary_data_growth, aes(
  x = day, 
  y = mean_biomass, 
  linetype = diet, 
  group = interaction(diet, line)
)) +
  geom_line(linewidth = 0.6, color = "black") +
  geom_point(size = 1, color = "black") +
  geom_errorbar(aes(ymin = mean_biomass - se_biomass,
                    ymax = mean_biomass + se_biomass),
                width = 0.5, size = 0.5, color = "black") +
  facet_wrap(~ line, scales = "free_y") +
  labs(
    x = "Time (Days)",
    y = "Biomass (g)"
  ) +
  scale_linetype_manual(values = c("CF" = "dashed", "FW" = "solid")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid   = element_blank(),
    strip.text   = element_text(size = 12, face = "bold"),
    axis.text    = element_text(size = 10, color = "black"),
    axis.title   = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    legend.position     = c(0.95, 0.1),
    legend.justification = c(1.1, 0.6)
  )


#### Table 1: GSI and accessory gland to gonad ratio ---------------------

# Clean and format data --------------------------------------------------

newdata <- newdata %>%
  mutate(
    Genetic.line = factor(Genetic.line),
    Treatment    = factor(Treatment),
    Sex          = factor(Sex),
    
    GSI = case_when(
      as.character(GSI) %in% c("", "NA", "N/A", "n/a", "-", "#DIV/0!", "Inf", "NaN") ~ NA_character_,
      TRUE ~ as.character(GSI)
    ),
    
    ag_gonad_ratio = case_when(
      as.character(ag_gonad_ratio) %in% c("", "NA", "N/A", "n/a", "-", "#DIV/0!", "Inf", "NaN") ~ NA_character_,
      TRUE ~ as.character(ag_gonad_ratio)
    ),
    
    GSI = as.numeric(GSI),
    ag_gonad_ratio = as.numeric(ag_gonad_ratio)
  )


# Summary and post-hoc function -----------------------------------------

get_summary <- function(df, measure) {
  
  df <- df %>%
    filter(!is.na(.data[[measure]]))
  
  if (nrow(df) == 0) {
    warning(paste("Column", measure, "is empty or has only NAs."))
    return(data.frame(
      group = character(),
      mean = numeric(),
      se = numeric(),
      letter = character()
    ))
  }
  
  summ <- df %>%
    group_by(group) %>%
    summarise(
      mean = mean(.data[[measure]], na.rm = TRUE),
      N    = sum(!is.na(.data[[measure]])),
      se   = sd(.data[[measure]], na.rm = TRUE) / sqrt(N),
      .groups = "drop"
    )
  
  if (length(unique(df$group)) < 2) {
    summ$letter <- "a"
    return(summ)
  }
  
  group_n <- df %>%
    group_by(group) %>%
    summarise(
      N = sum(!is.na(.data[[measure]])),
      .groups = "drop"
    )
  
  if (any(group_n$N < 2)) {
    warning(paste("Some groups have fewer than 2 observations for", measure))
  }
  
  bart <- bartlett.test(
    as.formula(paste(measure, "~ group")),
    data = df
  )
  
  if (bart$p.value > 0.05) {
    
    aov_model <- aov(
      as.formula(paste(measure, "~ group")),
      data = df
    )
    
    test_res <- agricolae::LSD.test(
      aov_model,
      "group",
      group = TRUE
    )
    
    letters <- test_res$groups
    
  } else {
    
    test_res <- agricolae::kruskal(
      df[[measure]],
      df$group,
      group = TRUE,
      p.adj = "bonferroni"
    )
    
    letters <- test_res$groups
  }
  
  letters <- data.frame(
    group = rownames(letters),
    letter = letters$groups,
    stringsAsFactors = FALSE
  )
  
  left_join(summ, letters, by = "group")
}


