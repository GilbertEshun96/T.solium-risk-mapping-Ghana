# ==========================================
# BAYESIAN PREVALENCE ANALYSIS: GHANA
# ==========================================
set.seed(123)

library(prevalence)
library(officer)
library(flextable)

# --- 1. BOLGATANGA (Ab-ELISA) ---
res_bolga_ab <- truePrev(x = 52, n = 137, SE = 0.95, SP = 0.95)

# --- 2. NAVRONGO (Ab-ELISA) ---
res_nav_ab <- truePrev(x = 2, n = 3, SE = 0.95, SP = 0.95)

# --- 3. BAWKU WEST / ZEBILLA (Ab-ELISA) ---
res_bawku_w_ab <- truePrev(x = 7, n = 9, SE = 0.95, SP = 0.95)

# --- 4. LEKMA (Ab-ELISA) ---
res_lekma <- truePrev(x = 0, n = 22, SE = 0.95, SP = 0.95)

# --- 5. AMA (Ab-ELISA) ---
res_ama <- truePrev(x = 0, n = 58, SE = 0.95, SP = 0.95)

# --- 6. WEIJA-GBAWE (Ab-ELISA) ---
res_weija <- truePrev(x = 0, n = 16, SE = 0.95, SP = 0.95)

# --- 7. KUMASI (Ante-mortem) ---
res_kumasi <- truePrev(x = 95, n = 4121,
                       SE = ~dunif(0.221, 0.387),
                       SP = ~dunif(0.90, 1.00))

# --- 8. WA (Post-mortem) ---
res_wa <- truePrev(x = 25, n = 303,
                   SE = ~dunif(0.221, 0.387),
                   SP = ~dunif(0.90, 1.00))

# --- 9. DAMANGO (Post-mortem) ---
res_damango <- truePrev(x = 9, n = 83,
                        SE = ~dunif(0.221, 0.387),
                        SP = ~dunif(0.90, 1.00))

# --- 10. YENDI (Post-mortem) ---
res_yendi <- truePrev(x = 8, n = 73,
                      SE = ~dunif(0.221, 0.387),
                      SP = ~dunif(0.90, 1.00))

# --- 11. NALERIGU (Post-mortem) ---
res_nalerigu <- truePrev(x = 4, n = 33,
                         SE = ~dunif(0.221, 0.387),
                         SP = ~dunif(0.90, 1.00))

# --- 12. BAWKU (Post-mortem) ---
res_bawku_pm <- truePrev(x = 7, n = 56,
                         SE = ~dunif(0.221, 0.387),
                         SP = ~dunif(0.90, 1.00))

# --- 13. BOLGATANGA (Post-mortem) ---
res_bolga_pm <- truePrev(x = 11, n = 81,
                         SE = ~dunif(0.221, 0.387),
                         SP = ~dunif(0.90, 1.00))

# --- 14. BONGO (Post-mortem) ---
res_bongo_pm <- truePrev(x = 13, n = 122,
                         SE = ~dunif(0.221, 0.387),
                         SP = ~dunif(0.90, 1.00))

# --- 15. BOUGOSONG (Post-mortem) ---
res_boug_pm <- truePrev(x = 4, n = 26,
                        SE = ~dunif(0.221, 0.387),
                        SP = ~dunif(0.90, 1.00))

# --- 16. CHIANA (Post-mortem) ---
res_chiana_pm <- truePrev(x = 8, n = 57,
                          SE = ~dunif(0.221, 0.387),
                          SP = ~dunif(0.90, 1.00))

# --- 17. GARU (Post-mortem) ---
res_garu_pm <- truePrev(x = 4, n = 38,
                        SE = ~dunif(0.221, 0.387),
                        SP = ~dunif(0.90, 1.00))

# --- 18. NAMOO (Post-mortem) ---
res_namoo_pm <- truePrev(x = 1, n = 43,
                         SE = ~dunif(0.221, 0.387),
                         SP = ~dunif(0.90, 1.00))

# --- 19. NAVRONGO (Post-mortem) ---
res_nav_pm <- truePrev(x = 12, n = 62,
                       SE = ~dunif(0.221, 0.387),
                       SP = ~dunif(0.90, 1.00))

# --- 20. WIDANA (Post-mortem) ---
res_widana_pm <- truePrev(x = 5, n = 31,
                          SE = ~dunif(0.221, 0.387),
                          SP = ~dunif(0.90, 1.00))

# --- 21. UPPER EAST (Post-mortem Multi) ---
res_ue_multi <- truePrev(x = 50, n = 495,
                         SE = ~dunif(0.221, 0.387),
                         SP = ~dunif(0.90, 1.00))

# --- 22. UER (Full carcass) ---
res_uer_full <- truePrev(x = 7, n = 60, SE = 1.00, SP = 1.00)

# --- 23. BURUNKURUGU (Post-mortem) ---
res_burunkurugu_pm <- truePrev(x = 22, n = 117,
                               SE = ~dunif(0.221, 0.387),
                               SP = ~dunif(0.90, 1.00))

# --- 24. TESHIE (Ab-ELISA) ---
res_teshie <- truePrev(x = 0, n = 22, SE = 0.95, SP = 0.95)

# --- 25. GBAWE (Ab-ELISA) ---
res_gbawe <- truePrev(x = 0, n = 6, SE = 0.95, SP = 0.95)

# --- 26. PAMBROS (Ab-ELISA) ---
res_pambros <- truePrev(x = 0, n = 10, SE = 0.95, SP = 0.95)

# --- 27. GLEFE (Ab-ELISA) ---
res_glefe <- truePrev(x = 0, n = 6, SE = 0.95, SP = 0.95)

# --- 28. SHAIBU (Ab-ELISA) ---
res_shaibu <- truePrev(x = 0, n = 12, SE = 0.95, SP = 0.95)

# --- 29. CHEMUNAA (Ab-ELISA) ---
res_chemunaa <- truePrev(x = 0, n = 39, SE = 0.95, SP = 0.95)

# --- 30. BOLA BEACH (Ab-ELISA) ---
res_bola <- truePrev(x = 0, n = 7, SE = 0.95, SP = 0.95)

# --- 31. BOLGATANGA (Ab-ELISA 2) ---
res_bolga_ab2 <- truePrev(x = 9, n = 10, SE = 0.95, SP = 0.95)

# --- 32. NAVRONGO (Ab-ELISA 2) ---
res_nav_ab2 <- truePrev(x = 2, n = 3, SE = 0.95, SP = 0.95)

# --- 33. ZEBILLA (Ab-ELISA) ---
res_zebilla <- truePrev(x = 7, n = 9, SE = 0.95, SP = 0.95)

# --- 34. OTHER UPPER (Ab-ELISA) ---
res_other <- truePrev(x = 68, n = 121, SE = 0.95, SP = 0.95)


# ==========================================
# EXTRACTION AND TABULATION
# ==========================================

all_res <- list(
  res_bolga_ab, res_nav_ab, res_bawku_w_ab, res_lekma, res_ama, res_weija,
  res_kumasi, res_wa, res_damango, res_yendi, res_nalerigu, res_bawku_pm,
  res_bolga_pm, res_bongo_pm, res_boug_pm, res_chiana_pm, res_garu_pm,
  res_namoo_pm, res_nav_pm, res_widana_pm, res_ue_multi, res_uer_full,
  res_burunkurugu_pm, res_teshie, res_gbawe, res_pambros, res_glefe,
  res_shaibu, res_chemunaa, res_bola, res_bolga_ab2, res_nav_ab2,
  res_zebilla, res_other
)

names_list <- c(
  "Bolgatanga", "Navrongo", "Bawku West/Zebilla", "LEKMA", "AMA", "Weija-Gbawe",
  "Kumasi", "Wa", "Damango", "Yendi", "Nalerigu", "Bawku",
  "Bolgatanga", "Bongo", "Bougosong", "Chiana", "Garu", "Namoo",
  "Navrongo", "Widana", "Upper East", "UER Full carcass",
  "Burunkurugu", "Teshie", "Gbawe", "Pambros", "Glefe",
  "Shaibu", "Chemunaa", "Bola Beach", "Bolgatanga", "Navrongo",
  "Zebilla", "Other Upper"
)

x_vals <- c(
  52, 2, 7, 0, 0, 0,
  95, 25, 9, 8, 4, 7,
  11, 13, 4, 8, 4, 1,
  12, 5, 50, 7,
  22, 0, 0, 0, 0,
  0, 0, 0, 9, 2,
  7, 68
)

n_vals <- c(
  137, 3, 9, 22, 58, 16,
  4121, 303, 83, 73, 33, 56,
  81, 122, 26, 57, 38, 43,
  62, 31, 495, 60,
  117, 22, 6, 10, 6,
  12, 39, 7, 10, 3,
  9, 121
)

extract_all_data <- function(mod_obj, loc_name, x_raw, n_raw) {
  tp_chain <- as.numeric(mod_obj@mcmc[[1]][[1]])
  
  inf_mean <- round(mean(tp_chain, na.rm = TRUE) * 100, 1)
  inf_lb   <- round(quantile(tp_chain, 0.025, na.rm = TRUE) * 100, 1)
  inf_ub   <- round(quantile(tp_chain, 0.975, na.rm = TRUE) * 100, 1)
  
  b_test <- binom.test(x_raw, n_raw)
  act_p  <- round((x_raw / n_raw) * 100, 1)
  act_lb <- round(b_test$conf.int[1] * 100, 1)
  act_ub <- round(b_test$conf.int[2] * 100, 1)
  
  data.frame(
    Location = as.character(loc_name),
    Actual = act_p,
    Actual_LB = act_lb,
    Actual_UB = act_ub,
    Informed = inf_mean,
    Informed_LB = inf_lb,
    Informed_UB = inf_ub,
    stringsAsFactors = FALSE
  )
}

final_df <- do.call(
  rbind,
  mapply(
    extract_all_data,
    all_res,
    names_list,
    x_vals,
    n_vals,
    SIMPLIFY = FALSE
  )
)

ft <- flextable(final_df) %>%
  set_header_labels(
    Actual_LB = "Actual LB",
    Actual_UB = "Actual UB",
    Informed_LB = "Informed LB",
    Informed_UB = "Informed UB"
  ) %>%
  colformat_double(digits = 1) %>%
  theme_vanilla() %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Ghana Porcine Cysticercosis: Final Prevalence Report", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "Ghana_Prevalence_Final_Ordered_1dp.docx")