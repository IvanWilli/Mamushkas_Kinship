# Sweden 1970-2017 parity kin with accumulated deaths.
# This script reads:
#   fltper_1x1.txt: female HMD period life table
#   mltper_1x1.txt: male HMD period life table
#   SWEpft.txt: period fertility table by birth order
# It prepares the full 1970-2017 list inputs
library(tidyverse)
source("R/kin_multi_stage_time_variant_2sex_with_death.R")

years <- 1970:2017
model_years <- 1970:2017
ns <- 3 # parity states: 0, 1, 2+
selected_kin <- c("d", "gd", "ggd", "ys", "os", "m", "gm", "ggm")
flt <- read.table(
  "data/fltper_1x1.txt",
  skip = 2,
  header = TRUE,
  na.strings = ".",
  stringsAsFactors = FALSE
)
mlt <- read.table(
  "data/mltper_1x1.txt",
  skip = 2,
  header = TRUE,
  na.strings = ".",
  stringsAsFactors = FALSE
)
fert <- read.table(
  "data/SWEpft.txt",
  skip = 2,
  header = TRUE,
  na.strings = ".",
  stringsAsFactors = FALSE
)
flt$age_num <- as.integer(gsub("[^0-9]", "", flt$Age))
mlt$age_num <- as.integer(gsub("[^0-9]", "", mlt$Age))
fert$age_num <- as.integer(gsub("[^0-9]", "", fert$x))
max_age <- max(flt$age_num[flt$Year %in% years], mlt$age_num[mlt$Year %in% years])
ages <- 0:max_age
na <- length(ages)
U_list_females <- vector("list", length(years))
U_list_males <- vector("list", length(years))
F_list_females <- vector("list", length(years))
T_list_females <- vector("list", length(years))
for (i in seq_along(years)) {
  y <- years[i]
  f_life <- flt[flt$Year == y, c("age_num", "qx")]
  m_life <- mlt[mlt$Year == y, c("age_num", "qx")]
  f_fert <- fert[
    fert$Year == y,
    c("age_num", "m1x", "m2x", "m3x", "m4x", "m5px", "l1x", "l2x", "l3x", "l4x")
  ]
  f_surv <- rep(0, na)
  m_surv <- rep(0, na)
  f_surv[match(f_life$age_num, ages)] <- 1 - f_life$qx
  m_surv[match(m_life$age_num, ages)] <- 1 - m_life$qx
  f_surv[is.na(f_surv)] <- 0
  m_surv[is.na(m_surv)] <- 0
  U_list_females[[i]] <- matrix(f_surv, nrow = na, ncol = ns)
  U_list_males[[i]] <- matrix(m_surv, nrow = na, ncol = ns)
  fert_mat <- matrix(0, nrow = na, ncol = ns)
  higher_order_rates <- as.matrix(f_fert[, c("m2x", "m3x", "m4x", "m5px")])
  higher_order_survivors <- as.matrix(f_fert[, c("l1x", "l2x", "l3x", "l4x")])
  higher_order_rates[is.na(higher_order_rates)] <- 0
  higher_order_survivors[is.na(higher_order_survivors)] <- 0
  higher_order_survivors_sum <- rowSums(higher_order_survivors)
  m2x <- ifelse(
    higher_order_survivors_sum > 0,
    rowSums(higher_order_rates * higher_order_survivors) / higher_order_survivors_sum,
    0
  )
  fert_rates <- cbind(m1x = f_fert$m1x, m2x = m2x)
  fert_rates[is.na(fert_rates)] <- 0
  fert_mat[match(f_fert$age_num, ages), ] <- cbind(fert_rates, fert_rates[, "m2x"])
  F_list_females[[i]] <- fert_mat
  T_list_females[[i]] <- vector("list", na)
  for (a in seq_len(na)) {
    T_age <- diag(ns)
    rates <- pmin(pmax(fert_mat[a, 1:(ns - 1)], 0), 1)
    for (stage in seq_len(ns - 1)) {
      T_age[stage, stage] <- 1 - rates[stage]
      T_age[stage + 1, stage] <- rates[stage]
    }
    T_age[ns, ns] <- 1
    T_list_females[[i]][[a]] <- T_age
  }
}
saveRDS(
  list(
    years = years,
    ages = ages,
    ns = ns,
    U_list_females = U_list_females,
    U_list_males = U_list_males,
    F_list_females = F_list_females,
    F_list_males = F_list_females,
    T_list_females = T_list_females,
    T_list_males = T_list_females,
    birth_female = 0.49,
    parity = TRUE
  ),
  "data/sweden_1970_2017_kin_multistage_inputs.rds"
)

run_idx <- match(1970:2017, years)
out <- kin_multi_stage_time_variant_2sex_with_death(
  U_list_females = U_list_females[run_idx],
  U_list_males = U_list_males[run_idx],
  F_list_females = F_list_females[run_idx],
  F_list_males = F_list_females[run_idx],
  T_list_females = T_list_females[run_idx],
  T_list_males = T_list_females[run_idx],
  birth_female = 0.49,
  parity = TRUE,
  summary_kin = FALSE,
  sex_Focal = "Female",
  initial_stage_Focal = 1, output_years = 2017
)
save(out, file = "data/out.rda")

# compute mamushka
load("data/out.rda")
library(tidyverse)

kin_country_year <- out %>%
    mutate(parity_kin = ifelse(stage_kin == 1, 0, 1)) %>%
    summarise(living = sum(living), dead = sum(dead), .by = c(kin, parity_kin, age_focal))

# M_a
M <- kin_country_year %>%
  filter(kin == "Focal") %>%
  summarise(M = 1-sum(living[parity_kin==0]), .by = c(age_focal)) %>%
  as.data.frame() %>% pull(M)
# M %>% plot()

# m_d
m_d <- kin_country_year %>%
  filter(kin == "d") %>%
  summarise(living = sum(living), dead = sum(dead), .by = c(age_focal, parity_kin)) %>%
  pivot_wider(names_from = parity_kin, values_from = living:dead) %>%
  mutate(m_d = (living_1 + dead_1)/(living_1 + dead_1 + living_0 + dead_0)) %>%
  replace(is.na(.), 0) %>%
  select(age_focal, m_d)
m_d <- m_d$m_d
# m_d %>% plot()

# gm_m
gm_m <- 1-(1-m_d)^1.78
# gm_m2 <- 1 - px(x = 0, n = cumsum(asfr_matrix[,1]), p = m_d, k = 1)
# plot(gm_m); lines(gm_m2)

  # GM_a
GM <- gm_m * M
# plot(GM)

# Kolk comparison, using saved Mamushka results only.
kolk_2017 <- readxl::read_xlsx("data/kolk_sweden_2017.xlsx")
mamushka_age <- seq_along(M) - 1
kolk_age_range <- range(
  kolk_2017$Edad[!is.na(kolk_2017$M) | !is.na(kolk_2017$GM)],
  na.rm = TRUE
)

comp_plot <- bind_rows(
  data.frame(age = mamushka_age, living = M, kin = "Mother", Type = "Mamushka"),
  data.frame(age = mamushka_age, living = GM, kin = "Grandmother", Type = "Mamushka"),
  kolk_2017 |>
    select(age = Edad, living = M) |>
    mutate(kin = "Mother", Type = "Register"),
  kolk_2017 |>
    select(age = Edad, living = GM) |>
    mutate(kin = "Grandmother", Type = "Register")
) |>
  mutate(
    kin = factor(kin, levels = c("Mother", "Grandmother")),
    Type = factor(Type, levels = c("Register", "Mamushka"))
  ) |>
  filter(age >= kolk_age_range[1], age <= kolk_age_range[2])

comparison_plot <- comp_plot |> 
  ggplot() +
  geom_line(
    aes(age, living, color = Type, linetype = Type),
    linewidth = 1.5) +
  facet_grid(
    ~ kin,
    scales = "free_y",
    switch = "y"
  ) + scale_color_brewer(palette="Set2") +
  labs(y = "", x = "Age") +
  scale_x_continuous(breaks = seq(0, 110, 10), expand = expansion(mult = c(0, 0.01))) +
  scale_y_continuous(breaks = seq(0, 1, .1)) +
  theme_bw(base_size = 13)

ggsave(
  "plots/sweden_2017_mamushka_kolk_comparison.pdf",
  comparison_plot,
  width = 9,
  height = 6
)



