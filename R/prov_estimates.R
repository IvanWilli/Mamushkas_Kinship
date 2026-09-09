library(tidyverse)
library(devtools)
library(foreign)
library(DemoTools)
source("funs.R")
load_all("C:/Proyectos/DemoKin")

# nacs deis
nacs_eevv <- read.csv("https://www.argentina.gob.ar/sites/default/files/2021/03/nacweb22.csv", sep = ";", header = F) %>% 
  select(Prov = 1, Edad = 4, B = 8) %>% 
  summarise(B = sum(B), .by = c(Prov, Edad)) %>% 
  mutate(Edad = case_when(substr(Edad, 1, 1) == 1 ~ 10,
                          substr(Edad, 1, 1) == 2 ~ 15,
                          substr(Edad, 1, 1) == 3 ~ 20,
                          substr(Edad, 1, 1) == 4 ~ 25,
                          substr(Edad, 1, 1) == 5 ~ 30,
                          substr(Edad, 1, 1) == 6 ~ 35,
                          substr(Edad, 1, 1) == 7 ~ 40,
                          substr(Edad, 1, 1) == 8 ~ 45))
sum(nacs_eevv$B)

# nacs deis
defs_eevv <- read.csv("https://www.argentina.gob.ar/sites/default/files/2021/03/defweb22.csv", sep = ";", header = F) %>% 
  select(Prov = 1, Sexo = 2, Edad = 5, D = 6) %>% 
  summarise(D = sum(D), .by = c(Prov, Sexo, Edad)) %>% 
  mutate(Edad = case_when(Edad == "01_Menor de 1 año" ~ "0",
                           Edad == "02_1 a 9" ~ "1",
                           Edad == "99_Sin especificar" ~ NA,
                          T ~ substr(Edad, 4, 6)),
         Edad = as.integer(Edad),
         Sexo = case_when(Sexo == 1 ~ "Varon",
                          Sexo == 2 ~ "Mujer",
                          T ~ NA))
defs_eevv <- bind_rows(
  defs_eevv %>% filter(Edad == 1) %>% mutate(D = D/2, Edad = 1),
  defs_eevv %>% filter(Edad == 1) %>% mutate(D = D/2, Edad = 5),
  defs_eevv %>% filter(Edad != 1),
)
sum(defs_eevv$D)

# pobl
provincias <- data.frame(
  cod_Prov = c(6, 2, 10, 22, 26, 14, 18, 30, 34, 38, 42, 46, 
               50, 54, 58, 62, 66, 70, 74, 78, 82, 86, 94, 90),
  Prov = c("Buenos Aires", "Caba", "Catamarca", "Chaco", "Chubut",
                "Córdoba", "Corrientes", "Entre Ríos", "Formosa", "Jujuy",
                "La Pampa", "La Rioja", "Mendoza", "Misiones", "Neuquén",
                "Río Negro", "Salta", "San Juan", "San Luis", "Santa Cruz",
                "Santa Fe", "Santiago del Estero", "Tierra del Fuego", "Tucumán"),
  stringsAsFactors = FALSE
)
pob <- readxl::read_xlsx("data/C2022_prov.xlsx", 
                         sheet = "Pob", range = "J14:N2806") %>% 
  mutate(Edad = ifelse(Edad == "edad", 0, as.integer(Edad))) %>%
  filter(Edad %in% 0:100) %>% 
  mutate(Edad = case_when(Edad == 0 ~ 0,
                          Edad %in% 1:4 ~ 1,
                          T ~ trunc(Edad/5)*5),  
         Mujer = as.integer(Mujer),
         Varon = as.integer(Varon)) %>% 
  summarise(Mujer = sum(as.integer(Mujer)),
            Varon = sum(as.integer(Varon)), .by = c(Prov, Edad)) %>% 
  replace(is.na(.), 0) %>% 
  left_join(provincias)

# data
provs_pob_def_nac <- pob %>% 
  left_join(defs_eevv %>% 
              pivot_wider(names_from = Sexo, values_from = D) %>% 
              rename(Mujer_D = Mujer, Varon_D = Varon, cod_Prov = Prov) %>% 
              select(-`NA`)) %>% 
  left_join(nacs_eevv %>% 
              rename(cod_Prov = Prov)) %>% 
  replace(is.na(.), 0) %>% 
  filter(Prov != "Total")

# TGF
TGF <- provs_pob_def_nac %>% 
  summarise(TGF = sum(B/Mujer, na.rm = T)*5, .by = Prov)

# Tablas
tm <- provs_pob_def_nac %>% 
  filter(Edad <= 80) %>% 
  arrange(Prov, Edad) %>% 
  split(.$Prov) %>% 
  map_df(function(r){
    print(r$Prov[1])
    bind_rows(
      lt_abridged2single(nMx = r$Mujer_D/r$Mujer, Age = r$Edad, Sex = "f", OAnew = 100) %>% mutate(Sexo = "Mujer"),
      lt_abridged2single(nMx = r$Varon_D/r$Varon, Age = r$Edad, Sex = "m", OAnew = 100) %>% mutate(Sexo = "Varon")
    ) %>% mutate(Prov = r$Prov[1], cod_Prov = r$cod_Prov[1])
  })

# plot
plot_e0_tgf <- inner_join(
  tm %>% 
    filter(Age == 0, Sexo == "Mujer") %>% select(Prov, `e(0)` = ex),
  provs_pob_def_nac %>% 
    summarise(TGF = sum(B/Mujer, na.rm = T)*5, Mujer = sum(Mujer)/sum(provs_pob_def_nac$Mujer)*100, 
              .by = Prov)) %>% 
  ggplot(aes(TGF, `e(0)`)) +
  geom_text(data = . %>% filter(TGF < 1.1 | TGF > 1.7), aes(label = Prov), 
            position = "jitter", vjust="inward",hjust="inward") +
  geom_point(alpha = .5, aes(size = Mujer)) +
  scale_x_continuous(labels = seq(.9, 2, .1), breaks = seq(.9, 2, .1)) +
  labs(size = "% Pobl") +
  theme_bw()
ggsave(plot = plot_e0_tgf, filename = "plots/plot_e0_tgf.pdf")

# hijos tenidos
HNVT <- readxl::read_xlsx("data/C2022_prov.xlsx", sheet = "HNVT", range = "w14:AP239") %>% 
  select(-Total) %>% 
  mutate_at(vars(`0`:`16`), as.integer) %>%
  replace(is.na(.), 0) %>% 
  pivot_longer(`0`:`16`, names_to = "HNVT", values_to = "W") %>% 
  mutate(HNVT = as.integer(HNVT)) %>% 
  filter(Edad != "Total") 

# check parides
HNVT %>%
  mutate(B = HNVT * W) %>% 
  summarise(Fx = sum(B)/sum(W), .by = c(Edad, Prov)) %>% 
  ggplot() +
  geom_boxplot(aes(Edad, Fx)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# último año
HNVUA <- readxl::read_xlsx("data/C2022_prov.xlsx", sheet = "HNVUA", range = "m14:q214") %>% 
  filter(Edad != "Total") 
  
# hijos tenidos y último año
HNVUA_HNVT <- readxl::read_xlsx("data/C2022_prov.xlsx", sheet = "HNVT_HNVUA", range = "w15:ap914") %>% 
  filter(HNUA == "Si") %>% 
  select(Prov, Edad, `1`) %>% 
  filter(Edad != "Total") 

HNV_prov <- HNVUA_HNVT %>% 
  left_join(HNVUA) %>% 
  left_join(HNVT %>% summarise(B = sum(HNVT*W), 
                               M = sum(W[HNVT>0])/sum(W),
                               W = sum(W), .by = c(Prov, Edad))) %>% 
  filter(Edad != "Total") %>% 
  mutate(fx = `Si`/W,
         pB1 = as.integer(`1`)/Si)

# check pB1
HNV_prov %>% 
  ggplot() +
  geom_boxplot(aes(Edad, pB1)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# add fx eevv
HNV_EEVV_prov <- HNV_prov %>%
  mutate(Edad = as.integer(substr(Edad, 1, 2))) %>% 
  left_join(provs_pob_def_nac %>% 
              select(Prov, Edad, Mujer, cod_Prov, B_eevv = B)) %>% 
  mutate(fx_eevv = B_eevv/Mujer)

# check TFR
HNV_EEVV_prov %>% 
  summarise(TFR_censo = sum(fx)*5, 
            TFR_eevv = sum(fx_eevv)*5, .by = Prov) %>% 
  ggplot(aes(TFR_censo, TFR_eevv)) +
  geom_point() +
  geom_text(aes(label = Prov)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")

# calcular f1
HNV_EEVV_prov <- HNV_EEVV_prov %>% 
  mutate(f1 = pB1/(1-M)*fx,
         f2 = (fx-(1-M)*f1)/M)

# f1 by prov
plot_f1_provs <- HNV_EEVV_prov %>% 
  filter(Prov != "Total", Edad > 10) %>% 
  ggplot(aes(Edad, f1)) +
  geom_smooth(aes(Edad + 2.5, f1), se = F, alpha = .2, col = "grey", span = .95) +
  geom_step() +
  facet_wrap(~ Prov) +
  theme_minimal()
ggsave(plot = plot_f1_provs, filename = "plots/plot_f1_provs.pdf")

# f1 y f2 caba y formosa
plot_f1_caba_misiones <- HNV_EEVV_prov %>% 
  filter(Prov %in% c("Caba", "Misiones"), Edad > 10) %>%
  select(Edad, Prov, "f1(x)" = f1, "f2+(x)" = f2) %>% 
  pivot_longer(`f1(x)`:`f2+(x)`) %>% 
  ggplot(aes(Edad, value, col = Prov)) +
  geom_step(linewidth = 1) +
  theme_minimal() +
  labs(y = "", x = "Age") +
  facet_wrap(~name)
ggsave(plot = plot_f1_caba_misiones, filename = "plots/plot_f1_caba_misiones.pdf")

# momushka
provs <- unique(tm$Prov)
prov_results_count <- map_df(provs, function(prov){
  print(prov)
  asfr_matrix_2022_prov <- HNV_EEVV_prov %>%
    filter(Prov == prov) %>% 
    select(Edad, fx, f1, f2) %>% 
    right_join(data.frame(Edad = 0:100) %>% 
                 mutate(Edad = trunc(Edad/5)*5), by = "Edad") %>% 
    arrange(Edad) %>% 
    replace(is.na(.), 0) %>% 
    select(-Edad) %>% as.matrix()
  ltf_matrix_2022_prov <- tm %>% 
    filter(Prov == prov, Sexo == "Mujer") %>% 
    mutate(p1x = 1-nqx, p2x = p1x) %>% 
    select(p1x, p2x)
  results_2022_prov <- get_momushka(asfr_matrix_2022_prov, 
                                    ltf_matrix_2022_prov)
  cod_prov = HNV_EEVV_prov %>% filter(Prov == prov) %>% slice(1) %>% pull(cod_Prov)
  return(results_2022_prov$results %>% mutate(Prov = prov, cod_Prov = cod_prov))
})
write.csv(prov_results_count, file = "data/prov_results_count.csv")

# plot
plot_prov_results_count <- prov_results_count %>% 
  select(Prov, Age, M, GM, GGM) %>% 
  pivot_longer(M:GGM) %>% 
  filter((name == "M" & Age %in% 15:50) |
         (name == "GM" & Age %in% 30:99) |
         (name == "GGM" & Age %in% 30:99)) %>% 
  mutate(name = factor(name, levels = c("M", "GM", "GGM"))) %>% 
  ggplot(aes(Age, value)) +
  geom_line(col = "grey", aes(group = Prov), alpha = .8) +
  geom_line(data = . %>% filter(Prov %in% c("Caba", "Misiones")), aes(col = Prov), size = 2) +
  geom_text(aes(label = Prov), data = . %>% 
              filter(Prov %in% c("Caba", "Misiones"), 
                     (Age == 30 & name == "M") | 
                     (Age == 70 & name == "GM") |
                     (Age == 90 & name == "GGM") )) +
  labs(col = "", y = "") +
  scale_x_continuous(breaks = seq(0, 100, 10), labels =  seq(0, 100, 10), expand = c(0,0)) +
  scale_y_continuous(breaks = seq(0, 1, .1), labels =  seq(0, 1, .1)) +
  theme_bw() +
  theme(legend.position = "") +
  facet_wrap(~name, scales = "free_x")
ggsave(plot = plot_prov_results_count, filename = "plots/prov_results_count.pdf")

# table 
prov_results_count %>% 
  filter(Age %in% c(30, 60, 90)) %>% 
  select(Prov, Age, M, D, GM, GGM) %>% 
  dput()

# poblacionalmente
table_provs_pobl <- prov_results_count %>% 
  select(Prov, Age, M, D, GM, GGM, `M-GM`, `GM-GGM`) %>% 
  left_join(readxl::read_xlsx("data/C2022_prov.xlsx", 
                              sheet = "Pob", range = "J14:N2806") %>% 
              mutate(Edad = ifelse(Edad == "edad", 0, as.integer(Edad)),
                     Mujer = as.integer(Mujer)) %>%
              filter(Edad %in% 0:100, Prov == "Total") %>% 
              select(Age = Edad, Mujer)) %>% 
  filter(!is.na(Age)) %>% 
  summarise(M = round(sum(M*Mujer)/sum(Mujer)*100, 1), # xM = sum(M*Mujer*Edad)/sum(M*Mujer),
            GM = round(sum(GM*Mujer)/sum(Mujer)*100, 1), # xGM = sum(GM*Mujer*Edad)/sum(GM*Mujer),
            GGM = round(sum(GGM*Mujer)/sum(Mujer)*100, 1), # xGGM = sum(GGM*Mujer*Edad)/sum(GGM*Mujer)
            .by = Prov) %>% 
  arrange(-M) %>% 
  as.data.frame()
