# microsimulation

#### 0. PARÁMETROS BÁSICOS ----
max_age  <- 100        # edad máxima
T_years  <- 300        # años de simulación
N0       <- 5000       # población inicial de mujeres
ages <- 0:max_age

#### 1. PERFILES DEMOGRÁFICOS (EJEMPLO) ----
load("data/ARG_inputs.rda")
fert <- asfr_matrix_2001[,1]; sum(fert)
fert1 <- asfr_matrix_2001[,2]
fert2p <- asfr_matrix_2001[,3]
qx <- 1-ltm_matrix_2001[,1]

## Función auxiliar: obtener probabilidad de muerte por edad
get_qx <- function(age_vec) {
  age_vec[age_vec > max_age] <- max_age
  qx[age_vec + 1]  # +1 porque ages parte en 0 y se indexa desde 1
}

#### 2. POBLACIÓN INICIAL ----
pop <- data.frame(
  id        = 1:N0,
  edad      = sample(0:max_age, N0, replace = TRUE),
  sexo      = "F",                  # solo mujeres al inicio
  viva      = TRUE,
  paridad   = 0L,                   # hijos tenidos (total)
  id_madre  = NA_integer_,          # NA para la 1ª generación
  id_padre  = NA_integer_           # NUEVO: columna para padre
)
N0m <- N0  # same size
pop_m <- data.frame(
  id        = (max(pop$id) + 1):(max(pop$id) + N0m),
  edad      = sample(0:max_age, N0m, replace = TRUE),
  sexo      = "M",
  viva      = TRUE,
  paridad   = 0L,
  id_madre  = NA_integer_,
  id_padre  = NA_integer_
)
pop <- rbind(pop, pop_m)

# --- TRACK SETUP (add once, before the loop) ---
track <- data.frame(
  t          = integer(),
  sex        = character(),
  age        = integer(),
  population = integer(),  # live population at BEGINNING of year t (by age, sex)
  deaths     = integer(),  # deaths during year t (by age at death, sex)
  newborns   = integer(),  # newborns during year t (by mother's age, newborn sex)
  stringsAsFactors = FALSE
)

#### 3. BUCLE DE MICROSIMULACIÓN ----
set.seed(1234)
  for (t in 1:T_years) {
    print(t)
    
    # --- TRACK: population at beginning of year (before ageing/deaths/births) ---
    pop_beg <- pop[pop$viva, c("sexo", "edad")]
    pop_beg_tab <- as.data.frame(table(
      sex = pop_beg$sexo,
      age = pop_beg$edad
    ))
    names(pop_beg_tab)[3] <- "population"
    
    track <- rbind(
      track,
      data.frame(
        t          = t,
        sex        = as.character(pop_beg_tab$sex),
        age        = as.integer(as.character(pop_beg_tab$age)),
        population = as.integer(pop_beg_tab$population),
        deaths     = 0L,
        newborns   = 0L
      )
    )
    
    ## 3.2 nuevas muertes
    vivos_idx <- which(pop$viva)
    if (length(vivos_idx) > 0) {
      death_prob <- get_qx(pop$edad[vivos_idx])
      mueren     <- runif(length(vivos_idx)) < death_prob
      
      # --- TRACK: deaths during year (by sex and age at death) ---
      if (any(mueren)) {
        dead_df <- pop[vivos_idx[mueren], c("sexo", "edad")]
        dead_tab <- as.data.frame(table(
          sex = dead_df$sexo,
          age = dead_df$edad
        ))
        names(dead_tab)[3] <- "deaths"
        
        track <- rbind(
          track,
          data.frame(
            t          = t,
            sex        = as.character(dead_tab$sex),
            age        = as.integer(as.character(dead_tab$age)),
            population = 0L,
            deaths     = as.integer(dead_tab$deaths),
            newborns   = 0L
          )
        )
      }
      
      pop$viva[vivos_idx[mueren]] <- FALSE
    }
    
    ## 3.3 Nacimientos (madres = mujeres vivas 15–49)
    potenciales_madres <- which(pop$viva & pop$sexo == "F" &
                                  pop$edad >= 15 & pop$edad <= 49)
    
    if (length(potenciales_madres) > 0) {
      
      edad_m <- pop$edad[potenciales_madres]
      par_m  <- pop$paridad[potenciales_madres]
      
      # Probabilidad de parto por orden
      prob_nac <- ifelse(par_m == 0,
                         fert1[edad_m + 1],   # orden 1
                         fert2p[edad_m + 1])  # orden 2+
      
      # A lo sumo un nacimiento por año por mujer (Bernoulli)
      hay_nac <- runif(length(prob_nac)) < prob_nac
      madres_que_paren <- potenciales_madres[hay_nac]
      
      n_nac <- length(madres_que_paren)
      if (n_nac > 0) {
        
        #### SECCIÓN CLAVE: ELECCIÓN DE PADRES ----
        # Regla:
        #  - varones vivos
        #  - edad en [edad_madre + 1, edad_madre + 4]
        #  - un mismo varón no puede ser padre de >1 hijo en ESTE año
        #  - pero sí puede ser padre nuevamente en años futuros
        
        padres_asignados <- rep(NA_integer_, n_nac)
        
        for (k in seq_len(n_nac)) {
          idx_madre <- madres_que_paren[k]
          edad_m_k  <- pop$edad[idx_madre]
          
          # candidatos: varones vivos, edad 1–4 años mayor, no usados ya este año
          candidatos <- which(
            pop$viva &
              pop$sexo == "M" &
              pop$edad >= (edad_m_k + 0) &
              pop$edad <= (edad_m_k + 0) &
              !(1:nrow(pop) %in% padres_asignados[1:(k - 1)])
          )
          
          if (length(candidatos) > 0) {
            padres_asignados[k] <- sample(candidatos, 1)
          } else {
            padres_asignados[k] <- NA_integer_  # sin padre admisible este año
          }
        }
        
        #### Crear nuevos individuos (recién nacidos) ----
        new_ids   <- (nrow(pop) + 1):(nrow(pop) + n_nac)
        sexos_nac <- ifelse(runif(n_nac) < 0.5, "F", "M")
        
        # --- TRACK: newborns by mother's age and newborn sex ---
        mother_age <- pop$edad[madres_que_paren]
        birth_tab <- as.data.frame(table(
          sex = sexos_nac,
          age = mother_age
        ))
        names(birth_tab)[3] <- "newborns"
        
        track <- rbind(
          track,
          data.frame(
            t          = t,
            sex        = as.character(birth_tab$sex),
            age        = as.integer(as.character(birth_tab$age)), # mother's age
            population = 0L,
            deaths     = 0L,
            newborns   = as.integer(birth_tab$newborns)
          )
        )
        
        nuevos <- data.frame(
          id        = new_ids,
          edad      = 0L,
          sexo      = sexos_nac,
          viva      = TRUE,
          paridad   = 0L,
          id_madre  = pop$id[madres_que_paren],
          id_padre  = ifelse(is.na(padres_asignados),
                             NA_integer_,
                             pop$id[padres_asignados])
        )
        
        #### Actualizar paridad de madres y padres ----
        pop$paridad[madres_que_paren] <- pop$paridad[madres_que_paren] + 1L
        
        padres_validos <- !is.na(padres_asignados)
        if (any(padres_validos)) {
          pop$paridad[padres_asignados[padres_validos]] <-
            pop$paridad[padres_asignados[padres_validos]] + 1L
        }
      }
    }
    
    ## 3.1 Envejecimiento (AL FINAL, siempre; solo sobre población existente)
    pop$edad[pop$viva] <- pop$edad[pop$viva] + 1L
    pop$viva[pop$edad > max_age] <- FALSE
    
    ## Agregar nacimientos al final del año (para que queden edad 0)
    if (!is.null(nuevos)) {
      pop <- rbind(pop, nuevos)
    }
  }

save(pop, track, file = "data/microsimulation.rda")

# nacs sin padres por falta de casos en edades [-1, +4]
recent <- which(pop$edad <= 100)  # everyone alive or dead in the last 100 cohorts
mean(is.na(pop$id_padre[recent & !is.na(pop$id_madre)]))

# indicadores -------------------------------------------------------------

load("data/microsimulation.rda")

#### CLASIFICACIÓN FINAL: MADRES, ABUELOS/AS, BISABUELOS/AS ----

pop2 <- pop %>%
  mutate(
    is_mother = (sexo == "F" & paridad > 0L),
    is_parent = (paridad > 0L)
  )

# --- helper: given a child-level indicator, flag parents who have ≥1 child with that indicator
flag_has_child_with <- function(df, child_flag_col) {
  child_flag <- rlang::ensym(child_flag_col)
  
  # maternal side
  m_side <- df %>%
    filter(!is.na(id_madre)) %>%
    group_by(id_madre) %>%
    summarise(has = any(!!child_flag, na.rm = TRUE), .groups = "drop") %>%
    rename(id = id_madre)
  
  # paternal side
  p_side <- df %>%
    filter(!is.na(id_padre)) %>%
    group_by(id_padre) %>%
    summarise(has = any(!!child_flag, na.rm = TRUE), .groups = "drop") %>%
    rename(id = id_padre)
  
  bind_rows(m_side, p_side) %>%
    group_by(id) %>%
    summarise(has = any(has), .groups = "drop")
}

# 4.2 Grandparents (total): ≥1 child who is a parent
has_child_parent <- flag_has_child_with(pop2, is_parent)

pop2 <- pop2 %>%
  left_join(has_child_parent, by = c("id")) %>%
  rename(is_grandparent_total = has) %>%
  mutate(
    is_grandparent_total = coalesce(is_grandparent_total, FALSE),
    is_grandmother_total = (sexo == "F" & is_grandparent_total)
  )

# 4.3 Great-grandparents (total): ≥1 child who is a grandparent
has_child_grandparent <- flag_has_child_with(pop2, is_grandparent_total)

pop2 <- pop2 %>%
  left_join(has_child_grandparent, by = c("id")) %>%
  rename(is_greatgrandparent_total = has) %>%
  mutate(
    is_greatgrandparent_total = coalesce(is_greatgrandparent_total, FALSE),
    is_greatgrandmother_total = (sexo == "F" & is_greatgrandparent_total)
  )

#### PROPORCIONES POR EDAD (mujeres vivas) ----

res <- pop2 %>%
  filter(viva, sexo == "F") %>%
  group_by(edad) %>%
  summarise(
    n_mujeres_vivas     = n(),
    prop_madres         = mean(is_mother),
    prop_abuelas_tot    = mean(is_grandmother_total),
    prop_bisabuelas_tot = mean(is_greatgrandmother_total),
    .groups = "drop"
  ) %>%
  right_join(tibble(edad = 0:max_age), by = "edad") %>%
  arrange(edad) %>%
  mutate(
    n_mujeres_vivas     = coalesce(n_mujeres_vivas, 0L),
    prop_madres         = coalesce(prop_madres, NA_real_),
    prop_abuelas_tot    = coalesce(prop_abuelas_tot, NA_real_),
    prop_bisabuelas_tot = coalesce(prop_bisabuelas_tot, NA_real_)
  )

# comparación con mamushka ------------------------------------------------

load("data/ARG_momushka.rda")
microsim_mamushka <- results_2001 %>% 
  select(Age, M, GM, GGM) %>% 
  mutate(Type = "Mamushka") %>% 
  bind_rows(res %>% 
              select(Age = edad, 
                     M = prop_madres, 
                     GM = prop_abuelas_tot, 
                     GGM = prop_bisabuelas_tot) %>% 
              mutate(Type = "Microsimulation"))

microsim_mamushka %>% 
  filter(!is.na(Age)) %>%
  pivot_longer(M:GGM, names_to = "kin", values_to = "Proportion") %>% 
  summarise(dif = Proportion[Type == "Mamushka"] - Proportion[Type != "Mamushka"], 
            .by = c(Age, kin)) %>% 
  ggplot(aes(Age, dif, col = kin)) +
  geom_point()

plot_microsim_momushka <- microsim_mamushka %>% 
  pivot_longer(M:GGM, names_to = "kin", values_to = "Proportion") %>% 
  ggplot(aes(Age, Proportion, col = kin)) +
  geom_line(data = . %>% filter(Type == "Mamushka"), size = 1.5) +
  geom_point(data = . %>% filter(Type != "Mamushka"), size = 2) +
  # geom_smooth(data = . %>% filter(Type != "Mamushka"), span = .2) +
  labs(col = "", y ="") +
  scale_x_continuous(expand=c(0,0), breaks = seq(0, 100, 10), labels =  seq(0, 100, 10)) +
  scale_y_continuous(expand=c(0,.05), breaks = seq(0, 1.1, .1), labels =  seq(0, 1.1, .1)) +
  scale_color_brewer(palette="Set2") +
  theme_bw()
plot_microsim_momushka
ggsave(plot = plot_microsim_momushka, filename = "plots/plot_microsim_momushka.pdf")

# tabla
table_microsim_momushka <- results_2001 %>% 
  select(Age, M, GM, GGM) %>% 
  left_join(res %>% 
              select(Age = edad, 
                     n_mujeres_vivas, 
                     M_micr = prop_madres, 
                     GM_micr = prop_abuelas_tot, 
                     GGM_micr = prop_bisabuelas_tot)) %>% 
  mutate(Age10 = trunc(Age/10)*10) %>% 
  filter(Age10 %in% 10:90, n_mujeres_vivas > 0) %>% 
  mutate(Age10 = paste0(Age10, "-", Age10+9)) %>% 
  replace(is.na(.), 0)
table_microsim_momushka <- bind_rows(
  table_microsim_momushka %>% 
    summarise(M_Mam = sum(M*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GM_Mam = sum(GM*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GGM_Mam = sum(GGM*n_mujeres_vivas)/sum(n_mujeres_vivas),
              M_Micr = sum(M_micr*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GM_Micr = sum(GM_micr*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GGM_Micr = sum(GGM_micr*n_mujeres_vivas)/sum(n_mujeres_vivas), .by = Age10),
  table_microsim_momushka %>% 
    summarise(M_Mam = sum(M*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GM_Mam = sum(GM*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GGM_Mam = sum(GGM*n_mujeres_vivas)/sum(n_mujeres_vivas),
              M_Micr = sum(M_micr*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GM_Micr = sum(GM_micr*n_mujeres_vivas)/sum(n_mujeres_vivas),
              GGM_Micr = sum(GGM_micr*n_mujeres_vivas)/sum(n_mujeres_vivas)) %>% 
    mutate(Age10 = "Total")) %>% 
  mutate(M_dif = M_Mam - M_Micr,
         GM_dif = GM_Mam - GM_Micr,
         GGM_dif = GGM_Mam - GGM_Micr) %>% 
  mutate(across(-Age10, ~ round(.x, 2))) %>% 
  as.data.frame()
table_microsim_momushka %>% dput()

# comparación overdispersion ----------------------------------------------

# finding over-dispersion in the microsimulation
library(dplyr)
library(aod)

df_counts <- pop %>% 
  filter(viva, sexo == "F") %>% 
  mutate(is_mother = paridad > 0) %>%
  group_by(edad) %>%
  summarise(y = sum(is_mother), n = n(), .groups = "drop")
fit <- betabin(cbind(y, n - y) ~ 1, ~ 1, data = df_counts)
kappa <- unname(coef(fit)["rho"])
kappa

# beta-binomial
px <- function(x, n, p, k = Inf) {
  # Binomial case: k = Inf (or very large)
  if (is.infinite(k)) {
    return(dbinom(x, size = n, prob = p))
  }
  # Beta-binomial case
  alpha <- k * p
  beta  <- k * (1 - p)
  # Stable computation:
  # P(X=x) = choose(n,x) * B(x+alpha, n-x+beta) / B(alpha, beta)
  log_p <- lchoose(n, x) + lbeta(x + alpha, n - x + beta) - lbeta(alpha, beta)
  exp(log_p)
}

# do beta binomial of probabilities
load("data/ARG_momushka.rda")
gm_m <- 1-(1-results_2001$m_d[-102])^cumsum(asfr_matrix_2001[,1])
gm_m2 <- 1 - px(x = 0, n = cumsum(asfr_matrix_2001[,1]), p = results_2001$m_d[-102], k = 5)
gm_m2 <- ifelse(is.na(gm_m2), 0, gm_m2)
plot(gm_m); lines(gm_m2)
GM2 <- gm_m2 * results_2001$M[-102]
ggm_gm <- 1-(1-results_2001$m_gd[-102])^results_2001$gd[-102]
ggm_gm2 <- 1 - px(x = 0, n = results_2001$gd[-102], p = results_2001$m_gd[-102], 
                  k = 5)
ggm_gm2 <- ifelse(is.na(ggm_gm2), 0, ggm_gm2)
plot(ggm_gm); lines(ggm_gm2)
GGM2 <- ggm_gm2 * GM2

# add to plot
plot_microsim_momushka_overdisp <- plot_microsim_momushka +
  geom_line(linetype = 2, linewidth = 1,
            data = data.frame(Age = 0:100, M = M, GM = GM2, GGM = GGM2) %>% 
              mutate(Type = "Beta-Binomial") %>% 
              pivot_longer(M:GGM, names_to = "kin", values_to = "Proportion")
              )
plot_microsim_momushka_overdisp
ggsave(plot = plot_microsim_momushka_overdisp, 
       filename = "plots/plot_microsim_momushka_overdisp.pdf")

# comparación de inputs ---------------------------------------------------

load("data/microsimulation.rda")

# get last 50 years
rates_check <- track %>%
  filter(t %in% 251:300) %>% 
  group_by(age) %>%
  summarise(
    Exposed_feamles = sum(population[sex == "F"]),
    Exposed = sum(population),
    Dx      = sum(deaths),
    qx_hat  = Dx / Exposed,
    Bx      = sum(newborns),
    fx_hat  = Bx / Exposed_feamles,
    .groups = "drop"
  )

# plot check
library(tidyverse)
plot_check_micro_inputs <- bind_rows(rates_check %>% 
                                       select(age, qx = qx_hat, fx = fx_hat) %>% 
                                       mutate(Type = "Microsimulation"),
                                     data.frame(age = 0:100, qx = qx, fx = fert, Type = "Input")) %>% 
  pivot_longer(qx:fx, names_to = "Indicator", values_to = "Value") %>% 
  ggplot(aes(age, Value)) +
  geom_point(data = . %>% filter(Type == "Microsimulation")) +
  geom_line(data = . %>% filter(Type != "Microsimulation")) +
  labs(col = "", y ="") +
  scale_x_continuous(expand=c(0,0), breaks = seq(0, 100, 10), labels =  seq(0, 100, 10)) +
  scale_y_log10(expand=c(0,0)) +
  scale_color_brewer(palette="Set2") +
  theme_bw() + labs(x = "Age", y = "") +
  facet_wrap(~Indicator)
plot_check_micro_inputs
ggsave(plot = plot_check_micro_inputs, 
       filename = "plots/plot_check_micro_inputs.pdf")

# # versión estando vivas ---------------------------------------------------
# 
# #### 4’. ESTATUS BASADOS EN DESCENDENCIA VIVA ----
# 
# n <- nrow(pop)
# 
# ## 4’.1 Hijos vivos (para todos, hombres y mujeres)
# 
# child_alive <- pop$viva  # cada individuo está vivo o no
# 
# # Rama materna: ¿algún hijo vivo?
# has_alive_child_m <- tapply(child_alive, pop$id_madre, any, na.rm = TRUE)
# 
# # Rama paterna: ¿algún hijo vivo?
# has_alive_child_p <- tapply(child_alive, pop$id_padre, any, na.rm = TRUE)
# 
# # Vector alineado con pop$id
# has_alive_child <- rep(FALSE, n)
# 
# if (!is.null(has_alive_child_m)) {
#   idx_m <- match(names(has_alive_child_m), pop$id)
#   has_alive_child[idx_m] <- has_alive_child[idx_m] | has_alive_child_m
# }
# 
# if (!is.null(has_alive_child_p)) {
#   idx_p <- match(names(has_alive_child_p), pop$id)
#   has_alive_child[idx_p] <- has_alive_child[idx_p] | has_alive_child_p
# }
# 
# ## Madres con al menos un hijo vivo
# is_mother_alive_child <- (pop$sexo == "F" & has_alive_child)
# 
# 
# ## 4’.2 Nietos vivos: hijos que tienen al menos un hijo vivo
# 
# child_has_alive_child <- has_alive_child   # propiedad de cada individuo: ¿tiene hijos vivos?
# 
# # Por la rama materna: ¿algún hijo que tenga hijos vivos?
# has_alive_grandchild_m <- tapply(child_has_alive_child, pop$id_madre, any, na.rm = TRUE)
# 
# # Por la rama paterna
# has_alive_grandchild_p <- tapply(child_has_alive_child, pop$id_padre, any, na.rm = TRUE)
# 
# has_alive_grandchild <- rep(FALSE, n)
# 
# if (!is.null(has_alive_grandchild_m)) {
#   idx_m2 <- match(names(has_alive_grandchild_m), pop$id)
#   has_alive_grandchild[idx_m2] <-
#     has_alive_grandchild[idx_m2] | has_alive_grandchild_m
# }
# 
# if (!is.null(has_alive_grandchild_p)) {
#   idx_p2 <- match(names(has_alive_grandchild_p), pop$id)
#   has_alive_grandchild[idx_p2] <-
#     has_alive_grandchild[idx_p2] | has_alive_grandchild_p
# }
# 
# ## Abuelas con al menos un nieto vivo
# is_grandmother_alive_grandchild <- (pop$sexo == "F" & has_alive_grandchild)
# 
# 
# ## 4’.3 Bisnietos vivos: hijos que tienen al menos un nieto vivo
# 
# child_has_alive_grandchild <- has_alive_grandchild
# 
# has_alive_greatgrandchild_m <- tapply(child_has_alive_grandchild, pop$id_madre, any, na.rm = TRUE)
# has_alive_greatgrandchild_p <- tapply(child_has_alive_grandchild, pop$id_padre, any, na.rm = TRUE)
# 
# has_alive_greatgrandchild <- rep(FALSE, n)
# 
# if (!is.null(has_alive_greatgrandchild_m)) {
#   idx_m3 <- match(names(has_alive_greatgrandchild_m), pop$id)
#   has_alive_greatgrandchild[idx_m3] <-
#     has_alive_greatgrandchild[idx_m3] | has_alive_greatgrandchild_m
# }
# 
# if (!is.null(has_alive_greatgrandchild_p)) {
#   idx_p3 <- match(names(has_alive_greatgrandchild_p), pop$id)
#   has_alive_greatgrandchild[idx_p3] <-
#     has_alive_greatgrandchild[idx_p3] | has_alive_greatgrandchild_p
# }
# 
# ## Bisabuelas con al menos un bisnieto vivo
# is_greatgrandmother_alive_ggchild <- (pop$sexo == "F" & has_alive_greatgrandchild)
# 
# #### 5’. PROPORCIONES POR EDAD (solo con descendencia viva) ----
# 
# mujeres_vivas <- (pop$viva & pop$sexo == "F")
# ages_final    <- 0:max_age
# 
# res_alive <- data.frame(
#   edad                        = ages_final,
#   n_mujeres_vivas             = sapply(ages_final, function(a)
#     sum(mujeres_vivas & pop$edad == a)),
#   prop_madres_con_hijo_vivo   = NA_real_,
#   prop_abuelas_con_nieto_vivo = NA_real_,
#   prop_bisabuelas_con_bn_vivo = NA_real_
# )
# 
# for (i in seq_along(ages_final)) {
#   a   <- ages_final[i]
#   idx <- which(mujeres_vivas & pop$edad == a)
#   if (length(idx) > 0) {
#     res_alive$prop_madres_con_hijo_vivo[i]   <- mean(is_mother_alive_child[idx])
#     res_alive$prop_abuelas_con_nieto_vivo[i] <- mean(is_grandmother_alive_grandchild[idx])
#     res_alive$prop_bisabuelas_con_bn_vivo[i] <- mean(is_greatgrandmother_alive_ggchild[idx])
#   }
# }
# 
# head(res_alive, 20)
# 
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# 
# ## 1) Unimos res (ever) y res_alive (living) en una sola tabla ancha
# 
# both <- res %>%
#   select(
#     edad,
#     mothers_ever      = prop_madres,
#     grandmothers_ever = prop_abuelas_tot,
#     ggms_ever         = prop_bisabuelas_tot
#   ) %>%
#   left_join(
#     res_alive %>%
#       select(
#         edad,
#         mothers_alive      = prop_madres_con_hijo_vivo,
#         grandmothers_alive = prop_abuelas_con_nieto_vivo,
#         ggms_alive         = prop_bisabuelas_con_bn_vivo
#       ),
#     by = "edad"
#   )
# 
# ## 2) Construimos explícitamente el long: 3 tipos de kin × 2 status
# 
# combined <- bind_rows(
#   # Mothers
#   both %>%
#     transmute(
#       edad,
#       kin     = "Mothers",
#       status  = "Ever",
#       proportion = mothers_ever
#     ),
#   both %>%
#     transmute(
#       edad,
#       kin     = "Mothers",
#       status  = "With living descendant",
#       proportion = mothers_alive
#     ),
#   # Grandmothers
#   both %>%
#     transmute(
#       edad,
#       kin     = "Grandmothers",
#       status  = "Ever",
#       proportion = grandmothers_ever
#     ),
#   both %>%
#     transmute(
#       edad,
#       kin     = "Grandmothers",
#       status  = "With living descendant",
#       proportion = grandmothers_alive
#     ),
#   # Great-grandmothers
#   both %>%
#     transmute(
#       edad,
#       kin     = "Great-grandmothers",
#       status  = "Ever",
#       proportion = ggms_ever
#     ),
#   both %>%
#     transmute(
#       edad,
#       kin     = "Great-grandmothers",
#       status  = "With living descendant",
#       proportion = ggms_alive
#     )
# )
# 
# ## (Opcional) chequeo rápido de duplicados
# combined %>%
#   count(edad, kin, status) %>%
#   filter(n > 1)
# # debería devolver 0 filas
# 
# ## 3) Gráfico superpuesto
# 
# ggplot(combined,
#        aes(x = edad, y = proportion,
#            color = kin, linetype = status)) +
#   geom_line(size = 1) +
#   scale_y_continuous(limits = c(0, 1),
#                      name   = "Proportion among women alive") +
#   scale_x_continuous(name = "Age") +
#   scale_color_discrete(name = "") +
#   scale_linetype_discrete(name = "") +
#   labs(
#     title    = "Age-specific proportions of mothers, grandmothers and great-grandmothers",
#     subtitle = "Ever vs with ≥1 living child, grandchild or great-grandchild"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     legend.position = "top",
#     legend.box      = "vertical"
#   )
# 
# ######### dif
# ombined %>% filter(kin == "Mothers", edad > 60) %>% 
#   pivot_wider(names_from = status, values_from = proportion) %>% 
#   mutate(dif = Ever - `With living descendant`) %>% as.data.frame()
# 
# 
# 
# head(res, 20)
# res$prop_madres %>% plot()
# mean(res$prop_madres[60:88])
# res$prop_abuelas %>% plot()
# mean(res$prop_abuelas[80:88])
# res$prop_bisabuelas %>% plot()
# mean(res$prop_bisabuelas[90:99])
# 
# res_plot <- res %>%
#   select(edad, prop_madres, prop_abuelas_tot, prop_bisabuelas_tot) %>%
#   pivot_longer(cols = -edad,
#                names_to = "tipo",
#                values_to = "proporcion")
# 
# # Etiquetas más amigables
# res_plot$tipo <- factor(res_plot$tipo,
#                         levels = c("prop_madres",
#                                    "prop_abuelas_tot",
#                                    "prop_bisabuelas_tot"),
#                         labels = c("Madres",
#                                    "Abuelas (totales)",
#                                    "Bisabuelas (totales)"))
# 
# # Gráfico
# (ggplot(res_plot, aes(x = edad, y = proporcion, color = tipo)) +
#     geom_line(size = 1) +
#     scale_y_continuous(limits = c(0,1), name = "Proporción") +
#     scale_x_continuous(name = "Edad") +
#     labs(title = "Proporción de madres, abuelas y bisabuelas por edad",
#          color = "") +
#     theme_minimal(base_size = 14) +
#     theme(legend.position = "top")) %>% plotly::ggplotly()