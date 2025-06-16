rm(list=ls())
library(dplyr)
library(splines)
# library(OSplines)
library(TMB)
library(aghq)
library(bayesplot)
library(lemon)
library(openxlsx)
library(lubridate)
library(magrittr)
library(reshape2)
library(janitor)
library(readr)
library(ggplot2)

# raw_case <- read_csv("~/Wastewater/reproduction_number/data_leighton/raw/covid-case-counts-moh.csv")
cases <- read_csv("./New Zealand/data/cases.csv")

# raw_case <- raw_case %>% clean_names()
# 
# raw_case %>% 
#   filter(report_date == '2020-02-26')
# 
# raw_case %>% 
#   filter(report_date == '2020-03-04')
# 
# raw_case %>% 
#   filter(report_date == '2020-03-16')
# 
# emily_process = raw_case %>% 
#   filter(case_status=="Confirmed") %>% 
#   group_by(report_date) %>% 
#   summarise(ncases = length(case_status),
#             nborders = length(which(district =="At the border")))

# cases %>% 
#   left_join(emily_process, by = c('date'='report_date'))
# 
# 
# cases$date[which(!(cases$date %in% emily_process$report_date))]
# emily_process$report_date[which(!(emily_process$report_date %in% cases$date ))]


ggplot(cases %>% filter(date >= '2022-03-01'&date <= '2022-04-01'), aes(date, nCases))+ 
  geom_bar(stat="identity")+ 
  scale_x_date(labels = function(x)wday(x,label = TRUE), breaks = scales::pretty_breaks(20))


cases_weekly <- cases%>% 
  mutate(number_of_cases = c(rep(NA,6), zoo::rollapply(nCases,7,sum))) %>% 
  filter(date %in% seq(ymd("2020-03-13"),max(cases$date),7)) 

cases_weekly2 <- 
  data.frame(date = seq(min(cases$date), max(cases$date),1)) %>% 
  left_join(cases, by = "date") %>% 
  mutate(dow = wday(date,label = TRUE),
         dow_num = wday(date)) %>% 
  mutate(starting_mon = date[which(dow=='Mon')[1]]) %>% 
  filter(date >= starting_mon) %>% 
  mutate(number_of_cases = c(zoo::rollapply(nCases,7,sum,na.rm=TRUE),rep(NA,6))) %>% 
  filter(date %in% seq(ymd("2020-03-02"),max(cases$date),7)) 
  
cases_epiweekly <- 
    data.frame(date = seq(min(cases$date), max(cases$date),1)) %>% 
    left_join(cases, by = "date") %>% 
    mutate(dow = wday(date,label = TRUE),
           dow_num = wday(date)) %>% 
    mutate(starting_sun = date[which(dow=='Sun')[1]]) %>% 
    filter(date >= starting_sun) %>% 
  mutate(number_of_cases = c(zoo::rollapply(nCases,7,sum,na.rm=TRUE),rep(NA,6))) %>% 
  filter(date %in% seq(ymd("2020-03-01"),max(cases$date),7)) 

# ggplot(cases, aes(date, nCases*10))+
#   geom_line()+ 
#   geom_point()+ 
#   geom_line(data = cases_weekly2, aes(date, number_of_cases), col = "red")+
#   geom_point(data = cases_weekly2, aes(date, number_of_cases), col = "red")+
#   geom_line(data = cases_weekly3, aes(date, number_of_cases), col = "blue")+
#   geom_point(data = cases_weekly3, aes(date, number_of_cases), col = "blue")+
#   scale_x_date(limits = c(ymd("2022-01-01"),ymd("2023-03-31")))

ggplot(cases, aes(date, nCases*10))+
  geom_line()+ 
  geom_line(data = cases_weekly, aes(date, number_of_cases), col = "red")+
  scale_x_date(limits = c(ymd("2022-01-01"),ymd("2023-03-31")))

cases_weekly <- cases_weekly %>% 
  dplyr::select(-'nCases',-'nBorderCases',-'nLocalCases',-'nCasesMovingAverage') %>% 
  rename('earliest_week_end_date' = date)

cases_epiweekly <- cases_epiweekly %>% 
  dplyr::select(-'nCases',-'nBorderCases',-'nLocalCases',-'nCasesMovingAverage',-'starting_sun',-'dow_num') %>% 
  rename('earliest_week_start_date' = date) %>% 
  mutate(earliest_week_end_date = earliest_week_start_date + 6) %>% 
  mutate(dow2 = wday(earliest_week_end_date, label = TRUE))


ggplot(cases, aes(date, nCases))+
  # geom_line()+
  # geom_point()+
  geom_line(aes(date, nCasesMovingAverage))+
  geom_point(aes(date, nCasesMovingAverage))+
  # geom_line(data = cases_epiweekly , aes(earliest_week_end_date, number_of_cases), col = "red")+
  # geom_point(data = cases_epiweekly , aes(earliest_week_end_date, number_of_cases), col = "red")+
  # geom_line(data = cases_epiweekly , aes(earliest_week_start_date, number_of_cases), col = "blue")+
  # geom_point(data = cases_epiweekly , aes(earliest_week_start_date, number_of_cases), col = "blue")+
  scale_x_date(limits = c(ymd("2022-01-01"),ymd("2023-03-31")))

zoo::rollmean(cases$nCases, 7, fill=NA, align="center")[1:10]
# ggplot()+ 
#   geom_line(data = cases_weekly, aes(date, number_of_cases), col = "red")+
#   geom_point(data = cases_weekly, aes(date, number_of_cases), col = "red")

save(file = "~/Wastewater/reproduction_number/cases_weekly_new_zealand.RData", list = "cases_epiweekly")

