# synthetic bacterial data

tbl <- "~/Local/gitdata/asms_2023/bac_data/Sheet 2-Table 1.csv" %>% read_csv()

fuzz <- function(x, mean = 1, sd = 0.1){
  x * rnorm(length(x), mean, sd)
}

tbl_new <- tbl %>%
  pivot_longer(matches('time'),
               names_to = 'time',
               values_to = 'abundance') %>%
  separate(bug,
           into = c('bug',"dose"),
           sep = "dose\\:") %>%
  group_by(Metabolite) %>%
  mutate(abundance = fuzz(abundance)) %>%
  ungroup() %>%
  mutate(time = str_extract(time, "[0-9]+")) %>%
  mutate(abundance = ifelse(bug == "e coli ", fuzz(abundance, 10, 5), abundance)) %>%
  mutate(abundance = ifelse(bug == "p aeruginosa ", fuzz(abundance, 10, 10), abundance)) %>%
  mutate(abundance = round(abundance)) %>%
  mutate(abundance = abs(abundance))

tbl_new %>%
  mutate(time = ifelse(time == 2, 70, time)) %>%
  ggplot(aes(time, abundance)) +
  geom_point(aes(color=dose)) +
  facet_grid(bug ~ Metabolite) +
  scale_y_log10()

tbl_new %>%
  mutate(time = glue::glue("Time_{time}min")) %>%
  mutate(time = ifelse(time == 'Time_2min', 'Time_2hr', time)) %>%
  pivot_wider(names_from = c('Metabolite','time'),
              values_from = 'abundance') %>%
  unite(Culture, bug, dose, sep="dose:") %>%
  writexl::write_xlsx(path = "data/bacterial-metabolites_dose-simicillin_wide.xlsx")

tbl_new %>%
  mutate(time = glue::glue("Time_{time}min")) %>%
  mutate(time = ifelse(time == 'Time_2min', 'Time_2hr', time)) %>%
  pivot_wider(names_from = c('time'),
              values_from = 'abundance') %>%
  unite(Culture, bug, dose, sep="dose:") %>%
  writexl::write_xlsx(path = "data/bacterial-metabolites_dose-simicillin_long.xlsx")

tbl_new %>%
  mutate(time = ifelse(time == '2', '120', time)) %>%
  mutate(time = as.numeric(time)) %>%
  mutate(dose = str_extract(dose, "[0-9]+")) %>%
  mutate(bug = trimws(bug)) %>%
  rename(time_min = time) %>%
  rename(dose_mg = dose) %>%
  rename(culture = bug) %>%
  rename(metabolite = Metabolite) %>%
  write_csv("data/bacterial-metabolites_dose-simicillin_tidy.csv")
