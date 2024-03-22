library(data.table)

t0 <- c(0.2, 0.4, 0.7, 0.8)

## {M_n} distribution
bounded_uniform <- function(N, lambda, p = 0.2){
  sample(
    x = seq(floor(lambda * (1 - p)), floor(lambda * (1 + p)), by = 1),
    size = N,
    replace = TRUE
  )
}


Hlogistic_d3_bis <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

get_real_data_mean_bis <- function(t = seq(0.1, 0.9, len = 10)){
  # \eta(t)
  cost_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * cos(2 * pi * k * ti))
  sint_mat <- outer(X = t, Y = 1:50, function(ti, k) sqrt(2) * sin(2 * pi * k * ti))
  eta <- cbind(1, cost_mat, sint_mat)

  # Basis coeffient
  basis_coef <- c(
    240.851203112216, 0.511305802240387, 0.0686061390530684, 0.404870711185904,
    0.163861146404104, 0.114790691976136, 0.104308592101677, 0.142469296029503,
    0.130499379756486, -0.0361260524404399, -0.0678188776098126, 0.00380999048511739,
    -0.0503325208862268, -0.0224294395747781, -0.0577702125832923, -0.0392160482687889,
    0.0287914974650662, -0.0088367407619057, 0.0342607237900621, 0.0273805020027171,
    0.0316142323237759, -0.0041035647844982, -0.00674473387407733, -0.013464230018548,
    -0.0534193098888712, 0, -0.0314081783398105, -0.00635579389377677, -0.0109138999402073,
    -0.00258921655373403, 0.0210582431037979, 0.00270039632612643, 0.0220059327400361,
    0.00498016602040829, 0.00612809662966225, -0.000468228505337897, -0.015169843922895,
    -0.0128459609384212, -0.0106467952350991, -0.0127090790046628, -0.00165413277187792,
    0.00462608335044012, 0.00161552211498331, 0.00387471684545869, 0.0168292422759318,
    0.00724332519473612, -0.00489353583983273, -0.0127769082071871, -0.0178778867020446,
    -0.0194110273298402, -0.0188371035035502, 0.422453270228039, 1.00539572630755,
    -0.24482290499743, -0.261068571761477, 0.00307387023311774, -0.116983058709441,
    -0.232838350406328, -0.0329785119290798, -0.0816897661769816, -0.0647624678180593,
    -0.0220324199101271, 0.0885506433535071, 0.0146671099049825, -0.037714370452663,
    0.0436172347334241, 0.0274110323188496, 0.00692447566082562, -0.0230790920539398, 0,
    -0.0055879729086357, -0.024539268711286, -0.00612357817533665, 0.0330311138598858,
    -0.0349369137221206, 0.00811058206224394, 0.0086712004621395, 0.00485477680861141,
    0.000146907988786585, 0.017293445884018, 0.00678226738589246, -0.0237728236554015,
    -0.00555310365060365, 0.00790951023830742, -0.0134462843070566, -0.0194592780694928,
    -0.0110838300364163, 0.00437116487075195, 0, 0.0165232688198793, 0.000983634631283174,
    -0.00550366281650836, 0.00681777170946617, 0, -0.00151544500544163, -0.0129957746034038,
    0.0164014779142208, 0, -0.0210844330910173, 0.00363142721126929, -0.00200062166486549
  )

  # mean function estimation
  muhat <- eta %*% basis_coef

  # Remove objects
  rm(cost_mat, sint_mat, eta, basis_coef)
  gc()

  return(muhat[, 1])
}

## mfBm
ker_d3_bis <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.2)

simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3_bis, Hvec = Hvec, design = "d3_new_3011")

simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_bis,
              process_mean = get_real_data_mean_bis, white_noise = "mfBm",
              hurst = Hlogistic_d3_bis, Hvec = Hvec, design = "d3_new_3011")

estim_locreg_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_3011")
estim_locreg_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_3011")

ggpubr::ggarrange(
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_3011", param = "Ht", Hfun = Hlogistic_d3_bis),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_3011", param = "Ht", Hfun = Hlogistic_d3_bis),
  ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_3011", param = "Lt", Hfun = Hlogistic_d3_bis),
  ggplot_locreg(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_3011", param = "Lt", Hfun = Hlogistic_d3_bis),
  nrow = 2, ncol = 2)

# ggplot_locreg(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d1", param = "Ht", Hfun = Hlogistic_d3_bis)


# simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")
#
# simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")
#
# simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")
#
# simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")
#
# simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")
#
# simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")
#
# simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211")
#
# simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
#               process = "FAR", process_ker = ker_d3_bis,
#               process_mean = get_real_data_mean_bis, white_noise = "mfBm",
#               hurst = Hlogistic_d3, Hvec = Hvec, design = "d3_new_2211_2")


## Merge data simulated data
merge_data_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211"){
  # design = "d3_new_2211"
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/to_delete/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(file_title)
  id_mc_max <- max(dt[, unique(id_mc)])

  # design = "d3_new_2211_2"
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/to_delete/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_2.RDS")
  dt_2 <- readRDS(file_title)

  # Merge
  dt_2[, id_mc := id_mc + id_mc_max]
  dt_merge <- rbind(dt, dt_2)

  # Save
  saveRDS(object = dt_merge,
          file = paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                                           process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_d3.RDS"))
}

merge_data_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
merge_data_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
merge_data_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")
merge_data_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", design = "d3_new_2211")


# If we use All real data ----
Hlogistic_d3_all <- function(t){
  hurst_logistic(t, h_left = 0.4, h_right = 0.6,
                 change_point_position = 0.5, slope = 50)
}

ker_d3_all <- function(s,t) get_real_data_far_kenel(s = s, t = t, operator_norm = 0.7)

#
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")

#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")

#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 75, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 75, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
#
gc() ; gc()
simulate_data(Nmc = 25, Ni = 150, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")
simulate_data(Nmc = 25, Ni = 1000, lambdai = 40, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")
simulate_data(Nmc = 25, Ni = 200, lambdai = 150, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")
simulate_data(Nmc = 25, Ni = 400, lambdai = 300, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")

#
gc() ; gc()
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_1")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_2")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_3")
simulate_data(Nmc = 75, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_4")
simulate_data(Nmc = 25, Ni = 1000, lambdai = 1000, t0, sig = 0.25,
              process = "FAR", process_ker = ker_d3_all,
              process_mean = get_real_data_mean, white_noise = "mfBm", Lt = 4,
              hurst = Hlogistic_d3_all, Hvec = Hvec, design = "d3_all_5")

merge_data_fun <- function(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all"){
  # design = "d3_new_2211"
  file_title <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                       process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,".RDS")
  dt <- readRDS(file_title)
  id_mc_max <- dt[, max(unique(id_mc))]

  for(i in 1:n_sup_basis){
    file_name <- paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_", design,"_", i, ".RDS")
    dt_2 <- readRDS(file_name)

    # Merge
    dt_2[, id_mc := id_mc + id_mc_max]
    dt <- rbind(dt, dt_2)
    id_mc_max <- dt[, max(unique(id_mc))]
  }

  # Save
  saveRDS(object = dt,
          file = paste0("./inst/12_mc_simulate_data/", process, "/data/dt_mc_",
                        process,"_", white_noise, "_", "N=", N, "_lambda=", lambda, "_d3.RDS"))
}

merge_data_fun(N = 150, lambda = 40, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 1000, lambda = 40, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 200, lambda = 150, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 400, lambda = 300, process = "FAR", white_noise = "mfBm", n_sup_basis = 5, design = "d3_all")
merge_data_fun(N = 1000, lambda = 1000, process = "FAR", white_noise = "mfBm", n_sup_basis = 2, design = "d3_all")
