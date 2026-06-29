# ==========================================
# BAYESIAN PREVALENCE ANALYSIS: HUM
# ==========================================
set.seed(123) 

library(prevalence)
library(officer)
library(flextable)

# --- 1. ACCRA (Ab-ELISA) ---
res_accra_ab <- truePrev(x = 8, n = 238, SE = 0.94, SP = 0.95)
res_accra_ab@model
res_accra_ab@diagnostics
plot(res_accra_ab@mcmc[[1]], main = "Accra: True Prev")

# --- 2. UPPER EAST (Ab-ELISA) ---
res_ue_ab <- truePrev(x = 2, n = 84, SE = 0.94, SP = 0.95)
res_ue_ab@model
res_ue_ab@diagnostics
plot(res_ue_ab@mcmc[[1]], main = "Upper East: True Prev")

# --- 3. WENCHI (Contrast-enhanced) ---
res_wenchi <- truePrev(x = 2, n = 2, SE = 0.936, SP = 0.821)
res_wenchi@model
res_wenchi@diagnostics
plot(res_wenchi@mcmc[[1]], main = "Wenchi: True Prev")

# --- 4. GOROGO (Stool microscopy) ---
res_gorogo <- truePrev(x = 27, n = 172, SE = 0.496, SP = 1)
res_gorogo@model
res_gorogo@diagnostics
plot(res_gorogo@mcmc[[1]], main = "Gorogo: True Prev")

# --- 5. BUNKPURUGU (Kato-Katz) ---
res_bunkpurugu <- truePrev(x = 65, n = 494, SE = ~dunif(0.07, 0.94), SP = ~dunif(0.96, 1.00))
res_bunkpurugu@model
res_bunkpurugu@diagnostics
plot(res_bunkpurugu@mcmc[[1]], main = "Bunkpurugu: True Prev")

# --- 6. JATO AKURA (Kato-Katz) ---
res_jato_akura <- truePrev(x = 3, n = 258, SE = ~dunif(0.07, 0.94), SP = ~dunif(0.96, 1.00))
res_jato_akura@model
res_jato_akura@diagnostics
plot(res_jato_akura@mcmc[[1]], main = "Jato Akura: True Prev")

# --- 7. ASHANTI (Multiplex PCR) ---
res_ashanti <- truePrev(x = 0, n = 2046, SE = ~dunif(0.88, 0.98), SP = ~dunif(0.94, 1.00))
res_ashanti@model
res_ashanti@diagnostics
plot(res_ashanti@mcmc[[1]], main = "Ashanti: True Prev")

# --- 8. ACCRA METRO (Stool microscopy) ---
res_accra_metro <- truePrev(x = 2, n = 300, SE = 0.496, SP = 1)
res_accra_metro@model
res_accra_metro@diagnostics
plot(res_accra_metro@mcmc[[1]], main = "Accra Metro: True Prev")

# --- 9. KUMASI 1 (Real-time PCR) ---
res_kumasi1 <- truePrev(x = 1, n = 905, SE = ~dunif(0.88, 0.98), SP = ~dunif(0.94, 1.00))
res_kumasi1@model
res_kumasi1@diagnostics
plot(res_kumasi1@mcmc[[1]], main = "Kumasi 1: True Prev")

# --- 10. KUMASI 33 (Western blot) ---
res_kumasi33 <- truePrev(x = 1, n = 905, SE = 0.97, SP = 0.99)
res_kumasi33@model
res_kumasi33@diagnostics
plot(res_kumasi33@mcmc[[1]], main = "Kumasi 33: True Prev")

# --- 11. KINTAMPO (RES33 antigen) ---
res_kintampo <- truePrev(x = 26, n = 571, SE = 0.99, SP = 0.93)
res_kintampo@model
res_kintampo@diagnostics
plot(res_kintampo@mcmc[[1]], main = "Kintampo: True Prev")

# --- 12. KUMASI 2 (Stool microscopy) ---
res_kumasi2 <- truePrev(x = 13, n = 2000, SE = 0.496, SP = 1)
res_kumasi2@model
res_kumasi2@diagnostics
plot(res_kumasi2@mcmc[[1]], main = "Kumasi 2: True Prev")

# --- 13. KINTAMPO N (Kato-Katz) ---
res_kintampo_n <- truePrev(x = 8, n = 443, SE = ~dunif(0.07, 0.94), SP = ~dunif(0.96, 1.00))
res_kintampo_n@model
res_kintampo_n@diagnostics
plot(res_kintampo_n@mcmc[[1]], main = "Kintampo N: True Prev")

# --- 14. BOLGATANGA (Ab-ELISA) ---
res_bolga_ab <- truePrev(x = 2, n = 30, SE = 0.94, SP = 0.95)
res_bolga_ab@model
res_bolga_ab@diagnostics
plot(res_bolga_ab@mcmc[[1]], main = "Bolgatanga: True Prev")

# --- 15. ZEBILLA (Ab-ELISA) ---
res_zebilla <- truePrev(x = 0, n = 31, SE = 0.94, SP = 0.95)
res_zebilla@model
res_zebilla@diagnostics
plot(res_zebilla@mcmc[[1]], main = "Zebilla: True Prev")

# --- 16. NAVRONGO (Ab-ELISA) ---
res_nav_ab <- truePrev(x = 0, n = 23, SE = 0.94, SP = 0.95)
res_nav_ab@model
res_nav_ab@diagnostics
plot(res_nav_ab@mcmc[[1]], main = "Navrongo: True Prev")

# --- 17. MALLAM (Ab-ELISA) ---
res_mallam <- truePrev(x = 0, n = 9, SE = 0.94, SP = 0.95)
res_mallam@model
res_mallam@diagnostics
plot(res_mallam@mcmc[[1]], main = "Mallam: True Prev")

# --- 18. OSU (Ab-ELISA) ---
res_osu <- truePrev(x = 0, n = 15, SE = 0.94, SP = 0.95)
res_osu@model
res_osu@diagnostics
plot(res_osu@mcmc[[1]], main = "Osu: True Prev")

# --- 19. GBAWE (Ab-ELISA) ---
res_gbawe <- truePrev(x = 1, n = 26, SE = 0.94, SP = 0.95)
res_gbawe@model
res_gbawe@diagnostics
plot(res_gbawe@mcmc[[1]], main = "Gbawe: True Prev")

# --- 20. BOLA BEACH (Ab-ELISA) ---
res_bola <- truePrev(x = 5, n = 40, SE = 0.94, SP = 0.95)
res_bola@model
res_bola@diagnostics
plot(res_bola@mcmc[[1]], main = "Bola Beach: True Prev")

# --- 21. CHEMUNAA (Ab-ELISA) ---
res_chemunaa <- truePrev(x = 1, n = 42, SE = 0.94, SP = 0.95)
res_chemunaa@model
res_chemunaa@diagnostics
plot(res_chemunaa@mcmc[[1]], main = "Chemunaa: True Prev")

# --- 22. GLEFE (Ab-ELISA) ---
res_glefe <- truePrev(x = 0, n = 2, SE = 0.94, SP = 0.95)
res_glefe@model
res_glefe@diagnostics
plot(res_glefe@mcmc[[1]], main = "Glefe: True Prev")

# --- 23. JAMES TOWN (Ab-ELISA) ---
res_jamestown <- truePrev(x = 0, n = 9, SE = 0.94, SP = 0.95)
res_jamestown@model
res_jamestown@diagnostics
plot(res_jamestown@mcmc[[1]], main = "James Town: True Prev")

# --- 24. KAJAANO (Ab-ELISA) ---
res_kajaano <- truePrev(x = 3, n = 9, SE = 0.94, SP = 0.95)
res_kajaano@model
res_kajaano@diagnostics
plot(res_kajaano@mcmc[[1]], main = "Kajaano: True Prev")

# --- 25. KORLE GONNO (Ab-ELISA) ---
res_korle <- truePrev(x = 0, n = 19, SE = 0.94, SP = 0.95)
res_korle@model
res_korle@diagnostics
plot(res_korle@mcmc[[1]], main = "Korle Gonno: True Prev")

# --- 26. PAMBROS (Ab-ELISA) ---
res_pambros <- truePrev(x = 0, n = 1, SE = 0.94, SP = 0.95)
res_pambros@model
res_pambros@diagnostics
plot(res_pambros@mcmc[[1]], main = "Pambros: True Prev")

# --- 27. PIG FARM (Ab-ELISA) ---
res_pigfarm <- truePrev(x = 0, n = 44, SE = 0.94, SP = 0.95)
res_pigfarm@model
res_pigfarm@diagnostics
plot(res_pigfarm@mcmc[[1]], main = "Pig farm: True Prev")

# --- 28. SHAIBU (Ab-ELISA) ---
res_shaibu <- truePrev(x = 0, n = 11, SE = 0.94, SP = 0.95)
res_shaibu@model
res_shaibu@diagnostics
plot(res_shaibu@mcmc[[1]], main = "Shaibu: True Prev")

# --- 29. TESHIE (Ab-ELISA) ---
res_teshie <- truePrev(x = 0, n = 11, SE = 0.94, SP = 0.95)
res_teshie@model
res_teshie@diagnostics
plot(res_teshie@mcmc[[1]], main = "Teshie: True Prev")


# ==========================================
# EXTRACTION AND TABULATION
# ==========================================

all_res <- list(res_accra_ab, res_ue_ab, res_wenchi, res_gorogo, res_bunkpurugu, res_jato_akura,
                res_ashanti, res_accra_metro, res_kumasi1, res_kumasi33, res_kintampo, res_kumasi2,
                res_kintampo_n, res_bolga_ab, res_zebilla, res_nav_ab, res_mallam, res_osu,
                res_gbawe, res_bola, res_chemunaa, res_glefe, res_jamestown, res_kajaano,
                res_korle, res_pambros, res_pigfarm, res_shaibu, res_teshie)

names_list <- c("Accra", "Upper East", "Wenchi", "Gorogo", "Bunkpurugu", "Jato Akura",
                "Ashanti", "Accra Metro", "Kumasi 1", "Kumasi 33", "Kintampo", "Kumasi 2",
                "Kintampo N", "Bolgatanga", "Zebilla", "Navrongo", "Mallam", "Osu",
                "Gbawe", "Bola Beach", "Chemunaa", "Glefe", "James Town", "Kajaano",
                "Korle Gonno", "Pambros", "Pig farm", "Shaibu", "Teshie")

x_vals <- c(8, 2, 2, 27, 65, 3, 0, 2, 1, 1, 26, 13, 8, 2, 0, 0, 0, 0, 1, 5, 1, 0, 0, 3, 0, 0, 0, 0, 0)
n_vals <- c(238, 84, 2, 172, 494, 258, 2046, 300, 905, 905, 571, 2000, 443, 30, 31, 23, 9, 15, 26, 40, 42, 2, 9, 9, 19, 1, 44, 11, 11)

extract_all_data <- function(mod_obj, loc_name, x_raw, n_raw) {
  tp_chain <- as.numeric(mod_obj@mcmc[[1]][[1]])
  inf_mean <- round(mean(tp_chain, na.rm = TRUE) * 100, 1)
  inf_lb   <- round(quantile(tp_chain, 0.025, na.rm = TRUE) * 100, 1)
  inf_ub   <- round(quantile(tp_chain, 0.975, na.rm = TRUE) * 100, 1)
  
  b_test <- binom.test(x_raw, n_raw)
  act_p  <- round((x_raw / n_raw) * 100, 1)
  act_lb <- round(b_test$conf.int[1] * 100, 1)
  act_ub <- round(b_test$conf.int[2] * 100, 1)
  
  data.frame(Location = loc_name, Actual = act_p, Actual_LB = act_lb, Actual_UB = act_ub,
             Informed = inf_mean, Informed_LB = inf_lb, Informed_UB = inf_ub)
}

final_df <- do.call(rbind, mapply(extract_all_data, all_res, names_list, x_vals, n_vals, SIMPLIFY = FALSE))

ft <- flextable(final_df) %>%
  set_header_labels(Actual_LB = "Actual LB", Actual_UB = "Actual UB",
                    Informed_LB = "Informed LB", Informed_UB = "Informed UB") %>%
  colformat_double(digits = 1) %>% theme_vanilla() %>% autofit()

doc <- read_docx() %>% 
  body_add_par("HUM: Final Bayesian Prevalence Analysis", style = "heading 1") %>%
  body_add_flextable(ft)

print(doc, target = "hum_prevalence_final.docx")