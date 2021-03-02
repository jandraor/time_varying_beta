get_apple_data <- function(start_date = "2020-02-28", 
                           end_date = "2020-05-17") {
  filepath        <- "./Data/Apple/applemobilitytrends-2020-06-03.csv"
  raw_data        <- read_csv(filepath)
  
  all_transp_data <- raw_data %>% filter(region == "Ireland") %>% 
    dplyr::select(-geo_type, -alternative_name, -`sub-region`, -country) %>% 
    pivot_longer(c(-region, -transportation_type), 
                 names_to = "date", values_to = "index") %>% 
    mutate(date = as_date(date)) %>% 
    filter(date >= start_date, date <= end_date) 
}

get_driving_data <- function() {
  
  apple_data <- get_apple_data(start_date = "2020-02-28",
                               end_date = "2020-05-17")
  
  driving_data <- filter(apple_data, transportation_type == "driving")
  first_val    <- driving_data[1, "index"] %>% as.numeric()
  
  driving_data <- driving_data %>% mutate(index = index / first_val,
                                          time = row_number() - 1) %>% 
    slice(-1)
  raw_data     <- driving_data$index
  imp          <- na_interpolation(raw_data)
  
  list(raw_data = raw_data,
       imputed_data = imp)
}