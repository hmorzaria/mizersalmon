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


extract_spatialoverlap <- function(ncfile){


  output <- nc_open(file.path(modelOutputsPath, "output.nc"))

  biomass_list <- NULL;
  for(g in 1:ngroups){
    this_numCohorts <- groups_df$NumCohorts[g]
    for(c in 1:this_numCohorts){
      thisVar <- groups_df$Name[g] %>%  str_trim(side="both") %>%  paste0(c)

      xx <- ncvar_get(nc = output, paste0(thisVar, "_Nums")) * (
        ncvar_get(nc = output, paste0(thisVar, "_ResN")) + ncvar_get(nc = output, paste0(thisVar, "_StructN"))) * mg_2_tonne * X_CN

      dimnames(xx) <- dimList
      this_biomass <- xx %>%  as_tibble() %>%
        mutate(layer = paste0("l",1:nlayers)) %>%
        pivot_longer(cols=!contains("layer")) %>%
        separate(name, into=c("box","ts")) %>%
        mutate("Name" = groups_df$Name[g], "Code" = groups_df$Code[g], "Cohort" = c)

      biomass_list[[(c-1)*ngroups + g]] <- this_biomass
    }
  }


  ## without stage, just species
  groups_biomassNOstage <- bind_rows(biomass_list) %>%
    mutate(ts = as.integer(ts)) %>%
    # rename("Biomass" = value) %>%
    left_join(age_mature, by = "Code") %>%
    group_by(layer, box, ts, Name, Code) %>%
    summarise("Biomass" = sum(value, na.rm=TRUE))

  biomass4spatialOverlap <- groups_biomassNOstage %>%
    filter(ts %in% timestep_4speicesOverlap)


  speciesCorr <- biomass4spatialOverlap %>%
    ungroup() %>%
    select(layer, box, ts, Code,  Biomass) %>%
    pivot_wider(values_from = Biomass, names_from = c(Code))  %>%
    select(groups_df$Code) %>%
    correlate(use = "everything")
  speciesOverlap <- speciesCorr %>%
    pivot_longer(cols = !contains("term")) %>%
    rename(corr = value) %>%
    mutate(Overlap = case_when(is.na(corr) ~ 1,
                               corr <= -0.9 ~ 0,
                               TRUE ~ (corr + 0.9)/(2 + 0.9))) %>% # 90% negative correlation will lead to zero overlap
    separate(term, into = c("Code_A")) %>%
    separate(name, into = c("Code_B"))
  ### if using stage in overlap (adul/juvenile) ##
  # speciesCorr <- biomass4spatialOverlap %>%
  #   ungroup() %>%
  #   select(layer, box, ts, Code, Stage, Biomass) %>%
  #   pivot_wider(values_from = Biomass, names_from = c(Code, Stage))  %>%
  #   select(contains("_")) %>%
  #   correlate(use = "everything")
  # speciesOverlap <- speciesCorr %>%
  #   pivot_longer(cols = !contains("term")) %>%
  #   rename(corr = value) %>%
  #   mutate(Overlap = case_when(is.na(corr) ~ 1,
  #                              corr <= -0.9 ~ 0,
  #                              TRUE ~ (corr + 0.9)/(2 + 0.9))) %>% # 80% negative correlation will lead to zero overlap
  #   separate(term, into = c("Code_A","Stage_A")) %>%
  #   separate(name, into = c("Code_B", "Stage_B")) %>%
  #   mutate(Stage_A = as.character(Stage_A),
  #          Stage_B = as.character(Stage_B))


}
