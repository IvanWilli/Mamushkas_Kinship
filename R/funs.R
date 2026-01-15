# main function
get_momushka <- function(asfr_matrix, ltf_matrix, pop = NULL){
  
  
  # asfr_matrix = asfr_matrix_2001; ltf_matrix = ltf_matrix_2001
  
  # prepare for DemoKin
  list_matrices <- make_mulstistate_parity_matrices(asfr_matrix[,-1], 
                                                    ltf_matrix, 
                                                    birth_female = 1)
  # run DemoKin
  kin_country_year <- 
    DemoKin::kin_multi_stage(
      U = list_matrices$U,
      f = list_matrices$F.,
      D = list_matrices$D,
      H = list_matrices$H, 
      birth_female=1,
      parity = TRUE)
  kin_country_year <- kin_country_year %>% 
    mutate(parity_kin = ifelse(stage_kin == 1, 0, 1)) %>% 
    summarise(living = sum(living), dead = sum(dead), .by = c(kin, parity_kin, age_focal))
  
  # check
    # kin_country_year %>% 
    #   filter(age_focal == 60) %>% 
    #   summarise(sum(living), .by = kin)
  
  # get results
  
  # M_a
  M <- kin_country_year %>% 
    filter(kin == "focal") %>% 
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
  gm_m <- 1-(1-m_d)^cumsum(asfr_matrix[,1])
  # gm_m2 <- 1 - px(x = 0, n = cumsum(asfr_matrix[,1]), p = m_d, k = 1)
  # plot(gm_m); lines(gm_m2)
  
  # GM_a 
  GM <- gm_m * M
  # plot(GM)
  
  # M
  M_only <- M - GM
  # plot(M_only)
  
  # m_gd
  m_gd <- kin_country_year %>% 
    filter(kin == "gd") %>% 
    summarise(living = sum(living), dead = sum(dead), .by = c(age_focal, parity_kin)) %>% 
    pivot_wider(names_from = parity_kin, values_from = living:dead) %>% 
    mutate(gd = living_1 + dead_1 + living_0 + dead_0,
           m_gd = (living_1 + dead_1)/gd) %>% 
    replace(is.na(.), 0)
  gd <- m_gd$gd
  m_gd <- m_gd$m_gd
  # m_gd %>% plot()
  # gd %>% plot()
  
  # ggm_gm
  ggm_gm <- 1-(1-m_gd)^gd
  # plot(ggm_gm)
  
  # GGM
  GGM <- ggm_gm * GM
  # plot(GGM)
  
  # GM only
  GM_only <- GM - GGM
  # plot(GM_only)
  
  # D only
  D <- 1 - M
  # plot(D)
  
  # stable distribution
  if(is.null(pop)){
    ages <- 101
    Ut = matrix(0, nrow=ages, ncol=ages)
    Ut[row(Ut)-1 == col(Ut)] <- ltf_matrix$p1x[-ages]
    Ut[ages, ages] = ltf_matrix$p1x[ages]
    ft = matrix(0, nrow=ages, ncol=ages)
    ft[1,] = asfr_matrix[,1]
    A = Ut + ft
    A_decomp = eigen(A)
    lambda = as.double(A_decomp$values[1])
    c <- as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
  }else{
    c <- pop/sum(pop) 
  }
  
  # results
  results <- tibble(Age = 0:100, M, m_d, m_gd, gm_m, ggm_gm, gd, 
                    D , GM, GGM, `M-GM` = M-GM, `GM-GGM` = GM-GGM, 
                    c)
  
  # # plot results
  # results %>% 
  #   select(Age, M, D, GM, GGM, `M-GM`, `GM-GGM`) %>% 
  #   pivot_longer(M:`GM-GGM`) %>% 
  #   ggplot(aes(Age, value, col = name)) +
  #   geom_line(size = 2) +
  #   labs(col = "", y = "") +
  #   scale_x_continuous(breaks = seq(0, 100, 10), labels =  seq(0, 100, 10), expand = c(0,0)) +
  #   scale_y_continuous(breaks = seq(0, 1, .1), labels =  seq(0, 1, .1)) +
  #   theme_bw()
  
  # add population results
  results <- results %>% 
    bind_rows(
      results %>%  
        summarise_at(vars(M:`GM-GGM`), ~ sum(.x * results$c)))
  
  return(list(results = results, kin_country_year = kin_country_year))
  
}

# function to create lists for the parity case given a set of coniditonal rates and survival probabilities with stages in columns and ages in rows
make_mulstistate_parity_matrices <- function(f_parity, p_parity, birth_female=.5){
  # f_parity = asfr_matrix[,-1]; p_parity = ltf_matrix; birth_female = 1
  
  ages <- nrow(f_parity)
  stages <- ncol(f_parity) + 1
  F_list <- U_list <- D_list <- H_list <- list()
  for(x in 1:ages){
    cond_probs <- as.numeric(f_parity[x,]/(1+f_parity[x,]/2))
    U_age <- matrix(0,stages,stages)
    diag(U_age) <-  c(1 - cond_probs, 1)
    U_age[row(U_age)-1==col(U_age)] <- cond_probs
    U_list[[x]] <- U_age
    F_age <- matrix(0,stages,stages)
    F_age[1,] <- c(cond_probs,cond_probs[stages-1])*birth_female
    F_list[[x]] <- F_age
  }
  p_parity <- as.matrix(cbind(p_parity, p_parity[,"p2x"]))
  for(s in 1:stages){
    H_age <- D_age <- matrix(0,ages,ages)
    H_age[1,] <- 1
    D_age[row(D_age)-1==col(D_age)] <- p_parity[-ages,s]
    D_age[stages, stages] <- p_parity[ages,s]
    H_list[[s]] <- H_age
    D_list[[s]] <- D_age
  }
  return(list(U = U_list, F. = F_list, H = H_list, D = D_list))
}
