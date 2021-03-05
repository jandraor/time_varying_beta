get_weekly_df <- function() {
  source("./R_scripts/irish_data.R")
  
  irish_data  <- get_irish_data()
  total_cases <- sum(irish_data$y)
  
  wkl_inc <- irish_data %>% slice(1:77)%>%
    mutate(week = ((time - 1) %/% 7) + 1) %>% 
    group_by(week) %>% summarise(y = sum(y)) %>% 
    mutate(time = week * 7)
  
  source("./R_scripts/apple.R")
  
  drv_data_obj <- get_driving_data()
  
  imp <- drv_data_obj$imputed_data
  
  wkl_df <- wkl_inc %>% rename(y1 = y) %>% 
    mutate(y2 = imp[1:length(imp) %% 7 == 0])
}
