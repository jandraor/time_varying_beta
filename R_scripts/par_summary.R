
summarise_pars <- function(par_names, cutoff_list) {
  map_df(par_names, function(par_name) {
    cutoff_CI <- cutoff_list[[par_name]]
    
    profile_par <- prof_est_df %>% 
      filter(profile_over == par_name) %>% 
      filter(loglik > cutoff_CI)
    
    lims <- range(profile_par[[par_name]] %>% as.numeric)
    
    data.frame(Parameter = par_name,
               lower_limit = lims[[1]], upper_limit = lims[[2]])
  })
}

## Taken from 'add reference'
## Weighted quantile function
wquant <- function (x, weights, probs = c(0.025,0.25, 0.5, 0.75, 0.975)) {
  idx <- order(x)
  x <- x[idx]
  weights <- weights[idx]
  w <- cumsum(weights)/sum(weights)
  rval <- approx(w,x,probs,rule=1)
  rval$y
}

var_quantiles_by_order <- function(var_df, rnd) {
  
  var_df |> 
    group_by(order) |> 
    summarise(mean   = mean(value),
              q_val  = quantile(value, c(0.025, 0.25, 0.5, 0.75, 0.975)),
              q_type = c("q2.5", "q25", "q50", "q75", "q97.5")) |> 
    ungroup() |> 
    mutate(q_val = round(q_val, rnd),
           mean  = round(mean, rnd)) |> 
    pivot_wider(names_from = q_type, values_from = q_val)
}

summarise_predicted_incidence <- function(inc_df) {
  
  inc_df <- inc_df |>  select(-order, -variable)
  
  wkl_df <- inc_df |>  filter(time <= 77) |> 
    mutate(week = ((time - 1) %/% 7) + 1) |>  
    group_by(iter, week) |> 
    summarise(value = sum(value))
  
  qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  wkl_df |>  group_by(week) |> 
    summarise(median = median(value), vals = quantile(value, qs),
              lims = c("q2.5", "q25", "q50", "q75", "q97.5")) |>  
    pivot_wider(names_from = lims, values_from = vals)
}

summarise_predicted_Z <- function(Z_df) {
  
  Z_df <- Z_df %>% select(-stock) %>% 
    mutate(week = time / 7)
  
  qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  Z_df |>  group_by(week) |> 
    summarise(median = median(value), vals = quantile(value, qs),
              lims = c("q2.5", "q25", "q50", "q75", "q97.5")) |>  
    pivot_wider(names_from = lims, values_from = vals)
}

summarise_predicted_Re <- function(pst_df, Z_df, dly_o) {
  
  Z_df <- Z_df |>  select(-stock) |> 
    mutate(week = time / 7)
  
  zeta_df <- pst_df |> select(zeta) |> 
    mutate(iter = row_number())
  
  Z_df |>  rename(Z = value) |>  left_join(zeta_df, by = "iter") |>  
    mutate(beta = zeta * Z,
           R    = estimate_r(beta)) -> R_df
  
  mdl_path  <- str_glue("./models/SEI3R_order_{dly_o}.stmx")
  mdl       <- read_xmile(mdl_path)
  stocks    <- sd_stocks(mdl)
  S_df      <- extract_timeseries_stock("S", pst_df, stocks, "x")
  
  S_df |>  mutate(s = value / 4937796) -> s_df
  
  bind_cols(R_df[, c("iter", "week", "R")], select(s_df, s)) |>  
    mutate(Re = R * s) -> Re_df
  
  qs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  Re_df |>  group_by(week) |> 
    summarise(median = median(R), vals = quantile(R, qs),
              lims = c("q2.5", "q25", "q50", "q75", "q97.5")) |> 
    pivot_wider(names_from = lims, values_from = vals)
}

summarise_hindcast <- function (var_name, filter_lik = FALSE) {
  
  summary_file <- str_glue("{folder}/summary_{var_name}.rds")
  
  if(!file.exists(summary_file)) {
    
    gc()
    file_name <- str_glue("{folder}/{var_name}_thesis.rds")
    df <- readRDS(file_name)
    
    if(filter_lik) {
      df <- df |> filter(loglik > max(loglik) - 20)
    }
      
    df <- df |> mutate(weight = exp(loglik - mean(loglik)))
    
    summary_hindcast <- df |>   group_by(time) |> 
      summarise(value = wquant(value, weight),
                type = c("q2.5", "q25", "q50", "q75", "q97.5")) |> ungroup() |>
      pivot_wider(names_from = type, values_from = value)
    
    saveRDS(summary_hindcast, summary_file)
    
    summary_hindcast
  } else {
    readRDS(summary_file)
  }
}

# Sample from the filtering distribution
sample_filt_dist <- function(var_input, folder, seed_val, filter_lik = FALSE) {
  
  set.seed(seed_val)
  
  fn_var <- file.path(folder, str_glue("{var_input}_thesis.rds"))
  
  var_df <- readRDS(fn_var)
  
  if(filter_lik) {
    var_df <- var_df |> filter(loglik > max(loglik) - 20)
  }
    
  var_df <- var_df |> mutate(weight = exp(loglik - mean(loglik)))
  
  map_df(1:11, function(i) {
    
    time_df <- filter(var_df, time == i) |>  
      mutate(rel_weight = weight / sum(weight))
    
    samples <- sample(time_df$value, 1e4, replace = TRUE,
                      prob = time_df$rel_weight)
    gc()
    
    data.frame(time = i, value = samples)
  }) |>  mutate(week = time)
}
  