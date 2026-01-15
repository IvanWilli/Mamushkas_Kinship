# compute mamushka distribution for ARG in 2001 and 2022

library(tidyverse)
library(devtools)
library(foreign)
library(DemoKin)
library(ggrepel)
source("R/funs.R")

# load VR and census data
load("data/ARG_inputs.rda")

# comparison againss hfd countries
HFD_countries <- read.table("data/mi.txt", skip = 2) %>%
  select(Country = 1, Year = 2, Age = 3, f1 = 4) %>%
  filter(Age %in% 15:49, Year >= 1980) %>%
  mutate(Age = as.integer(Age), f1 = as.numeric(f1))

# plot f1 for both years
ARG_f1_years <-
  bind_rows(
    asfr_matrix_2001 %>% as.data.frame() %>% select(Value = m1x) %>%
      mutate(Age = 0:100, Indicator ="f1(x)", Year = "2001"),
    asfr_matrix_2022 %>% as.data.frame() %>% select(Value = m1x) %>%
      mutate(Age = 0:100, Indicator ="f1(x)", Year = "2022")) %>%
  filter(Age %in% 10:50) %>%
  ggplot(aes(Age, Value, col = Year)) +
  geom_line(data = HFD_countries, aes(Age, f1, group = interaction(Country,Year)),
            col = "lightgrey", linewidth = .2) +
  geom_step(size = 1) +
  labs(y = "", col = "") +
  scale_x_continuous(limits = c(14, 50), labels = seq(10,50,5), breaks = seq(10,50,5)) +
  theme_bw()
ARG_f1_years
ggsave(plot = ARG_f1_years, filename = "plots/ARG_f1_years.pdf")

# get mamushka results
results_2001 <- get_momushka(asfr_matrix_2001, ltf_matrix_2001)
results_2022 <- get_momushka(asfr_matrix_2022, ltf_matrix_2022)

# save results
save(results_2001, results_2022, file = "data/ARG_momushka.rda")
load("data/ARG_momushka.rda")

# plot results
ARG_results_counts <- bind_rows(results_2001$results %>% mutate(Year = 2001),
          results_2022$results %>% mutate(Year = 2022)) %>% 
  select(Age, Year, M, D, GM, GGM, `M-GM`, `GM-GGM`) %>% 
  pivot_longer(M:`GM-GGM`) %>% 
  ggplot(aes(Age, value, col = name)) +
  geom_line(size = 2) +
  labs(col = "", y = "") +
  scale_x_continuous(breaks = seq(0, 100, 10), labels =  seq(0, 100, 10), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 1, .1), labels =  seq(0, 1, .1)) +
  theme_bw() +
  scale_color_brewer(palette="Set2") +
  facet_wrap(~Year)
ARG_results_counts
ggsave(plot = ARG_results_counts, filename = "plots/ARG_results_counts.pdf")

ARG_results_distr <- bind_rows(results_2001$results %>% mutate(Year = 2001),
          results_2022$results %>% mutate(Year = 2022)) %>% 
  select(Age, Year, 
         `Not mothers` = D, 
         `Only mothers` = `M-GM`, 
         `Only grandmothers` = `GM-GGM`,
         `At least great-grandmothers` = GGM) %>%
  pivot_longer(`Not mothers`:`At least great-grandmothers`) %>%
  mutate(name = factor(name, 
                       levels = c("Not mothers",
                                  "Only mothers",
                                  "Only grandmothers",
                                  "At least great-grandmothers"
                                  ))) %>% 
  ggplot(aes(Age, value, fill = name)) +
  geom_area() +
  labs(col = "", y = "", fill = "") +
  scale_x_continuous(expand=c(0,0), breaks = seq(0, 100, 10), labels =  seq(0, 100, 10)) +
  scale_y_continuous(expand=c(0,0), breaks = seq(0, 1, .1), labels =  seq(0, 1, .1)) +
  theme_bw() +
  scale_fill_brewer(palette="Set2") +
  facet_wrap(~Year)
ARG_results_distr
ggsave(plot = ARG_results_distr, filename = "plots/ARG_results_distr.pdf")

# table results
ARG_table_results <- bind_rows(results_2001$results %>% mutate(Year = 2001, .before = 1),
          results_2022$results %>% mutate(Year = 2022, .before = 1)) %>% 
  select(Year, Age, D, M, GM, m_d, m_gd, gm_m, GGM, ggm_gm) %>% 
  filter(Age %in% seq(15, 100, 5)) %>% 
  mutate(across(3:10, round, 2))

# piramyd under stable population
ARG_results_age_distr <- bind_rows(
          results_2001$results %>% mutate(Year = 2001, .before = 1),
          results_2022$results %>% mutate(Year = 2022, .before = 1)) %>%
  filter(!is.na(Age)) %>% 
  mutate(D = c*D,
        `M-GM` = c*`M-GM`,
        `GM-GGM` = c*`GM-GGM`,
        GGM = c*GGM) %>% 
  select(Year, Age, D, `M-GM`, `GM-GGM`, GGM) %>% 
  pivot_longer(D:GGM) %>% 
  ggplot(aes(Age, value, fill = name)) +
  labs(fill = "", y = "", x = "") +
  geom_area() + coord_flip() + scale_fill_brewer(palette="Set2") +
  theme_bw() +
  facet_grid(rows = vars(Year), switch = "y")
ARG_results_age_distr
ggsave(plot = ARG_results_age_distr, filename = "plots/ARG_results_age_distr.pdf")

# totals with stable population
bind_rows(results_2001$results %>% mutate(Year = 2001, .before = 1),
          results_2022$results %>% mutate(Year = 2022, .before = 1)) %>% 
  summarise(D = sum(c*D, na.rm = T),
            GM = sum(c*GM, na.rm = T),
            GGM = sum(c*GGM, na.rm = T),
            `M-GM` = sum(c*`M-GM`, na.rm = T),
            `GM-GGM` = sum(c*`GM-GGM`, na.rm = T),
            GGM = sum(c*GGM, na.rm = T),
            .by = Year)

# specific ages with stable population, 
bind_rows(results_2001$results %>% mutate(Year = 2001, .before = 1),
          results_2022$results %>% mutate(Year = 2022, .before = 1)) %>% 
  summarise(`GM-GGM_ma` = sum(`GM-GGM`*c*Age, na.rm = T)/sum(c*`GM-GGM`, na.rm = T),
            GGM_ma = sum(GGM*c*Age, na.rm = T)/sum(c*GGM, na.rm = T),
            `GM-GGM65plus` = sum(c[Age>=65]*`GM-GGM`[Age>=65], na.rm = T)/sum(c[Age>=65], na.rm = T),
            GGM65plus = sum(c[Age>=65]*GGM[Age>=65], na.rm = T)/sum(c[Age>=65], na.rm = T),
            D505plus = sum(c[Age>=50]*D[Age>=50], na.rm = T)/sum(c[Age>=50], na.rm = T),
            .by = Year)

# totals with c2022
load(file = "data/ARG_pop_2022.rda")
bind_rows(results_2001$results %>% mutate(Year = 2001, .before = 1),
          results_2022$results %>% mutate(Year = 2022, .before = 1)) %>% 
  select(-c) %>% 
  left_join(pop_2022 %>% filter(Sexo == "M") %>% select(Age = Edad, c = pop)) %>% 
  summarise(D = sum(c*D, na.rm = T)/sum(c, na.rm = T),
            GM = sum(c*GM, na.rm = T)/sum(c, na.rm = T),
            GGM = sum(c*GGM, na.rm = T)/sum(c, na.rm = T),
            `M-GM` = sum(c*`M-GM`, na.rm = T)/sum(c, na.rm = T),
            `GM-GGM` = sum(c*`GM-GGM`, na.rm = T)/sum(c, na.rm = T),
            GGM = sum(c*GGM, na.rm = T)/sum(c, na.rm = T),
            .by = Year)

# export results
write.csv(results, file = "ARG_momushka.csv", row.names = F)

# childleness plot --------------------------------------------------------

plot_childless <- readxl::read_xlsx("data/Childleness.xlsx") %>% 
  mutate(Country = ifelse(Country %in% c("Japan", "Sweden"), paste0(Country," ",Cohort), NA)) %>% 
  bind_rows(data.frame(CFR = c(2.5, 1.5), 
                       Childless = c(10, 30), 
                       Country = c("Arg 2001", "Arg 2022"))) %>% 
  filter(!is.na(CFR)) %>% 
  ggplot(aes(CFR, Childless, label=Country)) +
  geom_point() + geom_text_repel(
                                  segment.color = 'grey50',
                                  min.segment.length = 0,     # always draw segment
                                  point.padding = 0.6,        # keep label further from point (key)
                                  box.padding   = 0.8,        # keep labels further apart
                                  force         = 5,          # stronger repulsion
                                  force_pull    = 0,          # don't pull labels back to the point
                                  max.overlaps  = Inf) + 
  labs(y = "% Childless") +
  geom_smooth(span = .99, linetype = 2, se = F, col = 1) +
  geom_point(data = . %>% filter(Country %in% c("Arg 2001", "Arg 2020")), col = 2, size = 6) +
  theme_bw()
plot_childless
ggsave(plot = plot_childless, filename = "plots/plot_childless.pdf")


# siblings ----------------------------------------------------------------

# arg case
load("data/ARG_inputs.rda")

# Example fertility and survival rates

tfr_2001 <- cumsum(asfr_matrix_2001[,1])
Sib_x_2001 <- sapply(0:99, function(a){
  # a = 30
  tfr_a <- c(tfr_2001[-c(1:a)], rep(max(tfr_2001), a))
  sum(dist_mat_2001 * (1-exp(-tfr_a)*(1+tfr_a)))
})
tfr_2022 <- cumsum(asfr_matrix_2022[,1])
Sib_x_2022 <- sapply(0:99, function(a){
  # a = 30
  tfr_a <- c(tfr_2022[-c(1:a)], rep(max(tfr_2022), a))
  sum(dist_mat_2022 * (1-exp(-tfr_a)*(1+tfr_a)))
})

plot_siblings <- tibble(Age = 0:99, 
           "2001" = Sib_x_2001,
           "2022" = Sib_x_2022) %>% 
  pivot_longer(`2001`:`2022`, names_to = "Year", values_to = "S(x)") %>% 
  ggplot(aes(Age, `S(x)`, col = Year)) +
  geom_line(linewidth = 1) +
  labs(y  ="") +
  scale_color_brewer(palette="Set2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1), labels = seq(0,1,.1), breaks =  seq(0,1,.1)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 30)) +
  theme_bw()
plot_siblings
ggsave(plot = plot_siblings, filename = "plots/plot_sibling.pdf")

# comparison against poisson ----------------------------------------------

# c, gd and ggd ever had by age

total_kin_2001 <- results_2001$kin_country_year %>% 
  filter(kin %in% c("d", "gd", "ggd")) %>% 
  summarise(living = sum(living), 
            dead = sum(dead), 
            total = living + dead, 
            .by = c(age_focal, kin)) 
poisson_2001 <- total_kin_2001 %>% 
  summarise(M = 1-exp(-total[kin == "d"]*2), # same sex
            GM = 1-exp(-total[kin == "gd"]*2),
            GGM = 1-exp(-total[kin == "ggd"]*2), 
            GM = M * GM,
            GGM = GM *GGM, .by = age_focal)
total_kin_2022 <- results_2022$kin_country_year %>% 
  filter(kin %in% c("d", "gd", "ggd")) %>% 
  summarise(living = sum(living), 
            dead = sum(dead), 
            total = living + dead, 
            .by = c(age_focal, kin)) 
poisson_2022 <- total_kin_2022 %>% 
  summarise(M = 1-exp(-total[kin == "d"]*2), # same sex
            GM = 1-exp(-total[kin == "gd"]*2),
            GGM = 1-exp(-total[kin == "ggd"]*2), 
            GM = M * GM,
            GGM = GM *GGM, .by = age_focal)

# comparison plot
plot_poissson_app <- bind_rows(
  poisson_2001 %>% 
    rename(Age = age_focal) %>% 
    mutate(Model = "Poisson", Year = 2001),
  poisson_2022 %>% 
    rename(Age = age_focal) %>% 
    mutate(Model = "Poisson", Year = 2022),
  results_2001$results %>%
    select(Age, M, GM, GGM) %>% 
    mutate(Model = "Binomial", Year = 2001),
  results_2022$results %>%
    select(Age, M, GM, GGM) %>% 
    mutate(Model = "Binomial", Year = 2022)) %>% 
  pivot_longer(M:GGM, names_to = "kin", values_to = "total") %>% 
    ggplot(aes(Age, total, col = kin, linetype = Model)) + 
  geom_line() +
  scale_y_continuous(breaks = seq(0,1,.1),labels = seq(0,1,.1))+
  scale_x_continuous(limits = c(0, 99)) +
  theme_bw() + 
  facet_wrap(~Year)
plot_poissson_app
ggsave(plot = plot_poissson_app, filename = "plots/plot_poissson_app.pdf")

