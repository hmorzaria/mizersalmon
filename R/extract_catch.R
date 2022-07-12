#' @title extract species group spatial overlap by cohort
#' @param ncfile, name of output file from AMPS
#' @param func.groups, a data frame of salmon groups in the model
#'
#' @return spatial overlap
#' @export
#'
#' @descriptions for each pair of species groups, how much of their biomass overlaps spatially and temporally
#' @descriptions approximating with correlation over space and time
#' @author Hem Nalini Morzaria-Luna, hmorzarialuna_gmail.com July 2022


catches <- read.csv(file.path(modelOutputsPath, "outputCatch.txt"), sep=" ") %>%  as_tibble()


groups_catch <- catches %>%
  mutate(ts = Time/365+1) %>%
  select(c(ts, groups_df$Code)) %>%
  pivot_longer(cols = groups_df$Code) %>%
  rename(Catch = value, Code = name)
