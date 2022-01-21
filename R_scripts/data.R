get_data <- function() {
  
  irish_data   <- get_irish_data()
  
  drv_data_obj <- get_driving_data()
  
  daily_df <- data.frame(time = 1:79, 
                         y1 = irish_data$y,
                         y2 = drv_data_obj$df$y2)
  
  wkl_inc <- get_wkl_incidence(irish_data)
  
  imp <- drv_data_obj$df$y2
  
  wkl_df <- wkl_inc |>  rename(y1 = y) |>   
    mutate(y2 = imp[1:length(imp) %% 7 == 0])
  
  data_list <- list(Daily  = daily_df,
                    Weekly = wkl_df)
  
  data_list
}

get_wkl_incidence <- function(irish_data) {
  
  irish_data |> slice(1:77) |> 
    mutate(week = ((time - 1) %/% 7) + 1) |> 
    group_by(week) |> summarise(y = sum(y)) |>  
    mutate(time = week * 7)
}

get_irish_data <- function() {
  raw_data <- read_excel("./Data/irish_data.xlsx")
  
  total_cases <- slice(raw_data, nrow(raw_data)) |>  
    pull(CumulativeCases)
  
  irish_data <- raw_data |>  
    rename(time = Time,y = ReportedCases) |>  
    select(-CumulativeCases)
  
  first_date      <- as_date("2020-02-29")
  last_date       <- first_date + nrow(irish_data) - 1
  irish_data$date <- seq(first_date, last_date, 1)
  
  irish_data
}

get_apple_data <- function(start_date = "2020-02-28", 
                           end_date = "2020-05-17") {
  filepath        <- "./Data/Apple/applemobilitytrends-2020-06-03.csv"
  raw_data        <- read_csv(filepath)
  
  all_transp_data <- raw_data |>  filter(region == "Ireland") |>  
    dplyr::select(-geo_type, -alternative_name, -`sub-region`, -country) |>  
    pivot_longer(c(-region, -transportation_type), 
                 names_to = "date", values_to = "index") |>  
    mutate(date = as_date(date)) |>  
    filter(date >= start_date, date <= end_date) 
}

get_driving_data <- function() {
  
  apple_data <- get_apple_data(start_date = "2020-02-28",
                               end_date = "2020-05-17")
  
  driving_data <- filter(apple_data, transportation_type == "driving")
  first_val    <- driving_data[1, "index"] %>% as.numeric()
  
  driving_data <- driving_data |>  mutate(index = index / first_val,
                                          time = row_number() - 1) |>  
    slice(-1)
  raw_data     <- driving_data$index
  imp          <- na_interpolation(raw_data)
  
  df <- data.frame(date = seq(ymd('2020-02-29'), ymd('2020-05-17'), "day"),
                   time = 1:length(imp),
                   y2   = imp)
  
  list(raw_data = raw_data,
       imputed_data = imp,
       df = df)
}