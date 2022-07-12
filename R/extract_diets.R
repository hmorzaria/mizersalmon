

diets <- dietCheck %>%
  mutate(ts = Time/365)  %>% # given in days, but other outputs are in years
  left_join(age_mature, by = c("Predator" = "Code")) %>%
  mutate(Cohort = Cohort + 1) %>%  # starts at zero in dietCheck
  mutate("Stage" = case_when(
    Cohort < Cohort_mature ~ "Juvenile",
    TRUE ~ "Adult"
  )) %>%
  filter(Predator %in% groups_df$Code) %>%
  select(!contains(c("Stock","Updated", "Cohort", "Time", "age_mat", "Cohort_mature"))) %>%
  pivot_longer(cols = !contains(c("Time","Predator","Stage","ts"))) %>%
  rename("Prey" = name) %>%
  group_by(Predator, Stage, Prey, ts) %>%
  summarise(Consumed = sum(value)) %>%
  ungroup() %>%
  group_by(Predator, Stage, Prey) %>%
  summarise(Diet = mean(Consumed)) %>%
  filter(Diet>0)

diets_simplified <- diets %>%
  mutate(Prey = case_when(Prey %in% groups_df$Code ~ Prey,
                          TRUE ~ "Other")) %>%
  group_by(Predator, Stage, Prey) %>%
  summarise(Diet = sum(Diet, na.rm=TRUE))

# as proportions
diets_proportions <- diets_simplified %>%
  left_join(
    diets_simplified %>%
      group_by(Predator, Stage) %>%
      summarise(total_diet = sum(Diet, na.rm=TRUE)),
    by = c("Predator", "Stage")
  ) %>%
  mutate(Diets_proportion = Diet/total_diet)
