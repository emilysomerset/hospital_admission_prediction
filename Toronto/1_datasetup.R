library(readxl)
library(magrittr)
library(tidyverse)
library(janitor)
library(lubridate)
library(forcats)
# Check whether total = sum (cases by age)
load("./Toronto/Data/work_d_toronto_2024_08_08.RData")
case_age <- read_excel("./Toronto/Data/Laboratory confirmed weekly case counts or rates of COVID-19 by age group in Toronto Public Health.xlsx",
                       col_types = c("date", "numeric", "date", "date", "text"))
case_age <- case_age %>% 
  clean_names() %>% 
  mutate(age_group_report = factor(age_group_report)) %>% 
  rename(case = number_of_cases_or_episodes,
         age=age_group_report) %>% 
  mutate_at(vars(grep("_date", colnames(.), value = TRUE)), ymd)


case_age <- case_age %>% 
  filter(earliest_week_end_date <= "2024-06-01" & earliest_week_end_date>="2020-03-07")

which(!((case_age %>% 
  group_by(earliest_week_end_date) %>% 
  summarise(case = sum(case, na.rm = TRUE)) %>% 
  filter(earliest_week_end_date >= '2020-03-07') %$% case) == work_d_toronto$number_of_cases))

## good. 

### Look into the unknown age groups
case_age %>% filter(age=="Unknown") %>% 
  filter(!is.na(case)) #21 cases with unknown age groups.

case_age <- case_age %>% 
  filter(!(is.na(case) & age == "Unknown"))

case_age<- case_age %>% 
  mutate(case = ifelse(is.na(case),0, case))

##### Hospital beds
hosp_age <- read_excel("./Toronto/Data/COVID-19 bed occupancy by age group in Toronto Public Health.xlsx",
                       col_types = c("date", "numeric","text", "date", "date"))

hosp_age <- hosp_age %>% 
  clean_names()
colnames(hosp_age) <- c("week_start_date","average_bed_occupancy","age","earliest_week_start_date","earliest_week_end_date")

hosp_age <- hosp_age %>% 
  filter(earliest_week_end_date <= "2024-06-01" & earliest_week_end_date>="2020-03-07")%>% 
  mutate_at(vars(grep("_date", colnames(.), value = TRUE)), ymd) 

hosp_age$age %>% unique()
case_age$age %>% unique()

case_age <- case_age %>% 
  mutate(age = fct_collapse(age, `0 to 4` = c("<1","1 to 4"))) %>% 
  group_by(age, earliest_week_end_date,earliest_week_start_date, week_start_date) %>% 
  summarise(case = sum(case))

hosp_age$age <- factor(hosp_age$age, levels = c("0 to 4","5 to 11",  "12 to 17", "18 to 39", "40 to 59", "60 to 79", "80+" ))

case_age$age <- factor(case_age$age, levels = c("0 to 4","5 to 11",  "12 to 19", "20 to 39", "40 to 59", "60 to 79", "80+","Unknown" ))

hosp_age$age2 = hosp_age$age
hosp_age$age2 = factor(hosp_age$age2, label = c("0 to 4","5 to 11",  "12 to 19", "20 to 39", "40 to 59", "60 to 79", "80+"))
  

case_age <- case_age %>% 
  left_join(hosp_age %>% rename("age_orig_hosp"="age") %>% 
              dplyr::select("earliest_week_end_date","age2","age_orig_hosp","average_bed_occupancy"), by = c("earliest_week_end_date","age"="age2")) 


case_age <- case_age %>% 
  arrange(-desc(earliest_week_end_date)) %>% 
  group_by(earliest_week_end_date) %>% 
  mutate(totcase = sum(case)) 

plausible_total <- function(avg, n = 7, digits = 1) {
  # Compute lower and upper bounds for the true average before rounding
  lower_bound <- avg - 0.5 * 10^(-digits)
  upper_bound <- avg + 0.5 * 10^(-digits)
  
  # Pick the midpoint of the interval as a plausible true average
  plausible_avg <- (lower_bound + upper_bound) / 2
  
  # Compute total and round to integer
  total <- round(plausible_avg * n)
  
  return(total)
}

case_age <- case_age %>% 
  mutate(integer_number = plausible_total(average_bed_occupancy)) 

case_age <- case_age %>% 
  ungroup() %>% 
  slice(592:nrow(.))


## total hospital admissions
hosp_adm <- read_excel("Toronto/Data/CityofToronto_COVID-19_Data.xlsx", 
                       sheet = "10- Severity by Epi Week", col_types = c("date", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "numeric"))

hosp_adm<- hosp_adm %>% 
  clean_names() %>% 
  mutate(episode_week = ymd(episode_week)) %>% 
  rename(earliest_week_start_date = "episode_week")

case_age<- case_age %>% 
  left_join(hosp_adm, by = "earliest_week_start_date")

save(file = './Toronto/work_d.RData',  list= "case_age")
