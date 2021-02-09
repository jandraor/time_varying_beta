get_irish_data <- function() {
  raw_data <- read_excel("./Data/irish_data.xlsx")
  
  total_cases <- slice(raw_data, nrow(raw_data)) %>% 
    pull(CumulativeCases)
  
 irish_data <- raw_data %>% 
   rename(time = Time,y = ReportedCases) %>% 
   select(-CumulativeCases)
 
 first_date      <- as_date("2020-02-29")
 last_date       <- first_date + nrow(irish_data) - 1
 irish_data$date <- seq(first_date, last_date, 1)
 
 irish_data
}