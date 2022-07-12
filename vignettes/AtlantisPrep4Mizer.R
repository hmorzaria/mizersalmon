# read in Atlantis model and outputs and create files for setting up Mizer model
#
# files read in: 
# model inputs: biology prm, groups.csv
# model outputs: output.nc, outputCATCH.txt, outputDietCheck.txt
#
# files created: 
# species_params (Species, Name,  w_inf,  w_mat, interaction_ASQ, interaction_BAL,...)
# diet_proportions (Predator, Stage, Prey, Diet, total_diet, Diets_proportion)
# groups_biomass_catch (Timestep, Code,  Catch, Biomass, Effort)
#
source(file.path(DIR$AssessmentDevelopment, "get_first_number.R"))

modelInputsPath <- file.path(DIR$Base, "FromCRAM2017Repository","ArchivedModels", "Base")
modelOutputsPath <- file.path(DIR$Base, "FromCRAM2017Repository","ArchivedModels", "Base", "outputFish")
mizerModelPath <- file.path(DIR$Base, "Models","Mizer")

timestep_4initial <- 1 # use this/these timesteps for B0 (takes mean if more than one timestep)
timestep_4speicesOverlap <- 1:5 # use this/these for species overlap (taken from correlation wrt space and time)
timesteps_4diets <- 1:50 # summarise diets over these timesteps

biol_prm <- readLines(file.path(modelInputsPath, "CRAM_BH_hybrid_biol.prm"))
output <- nc_open(file.path(modelOutputsPath, "output.nc"))
catches <- read.csv(file.path(modelOutputsPath, "outputCatch.txt"), sep=" ") %>%  as_tibble()
dietCheck <- read.csv(file.path(modelOutputsPath, "outputDietCheck.txt"), sep=" ") %>%  as_tibble()

groups_df <- read.csv(file.path(modelInputsPath, "CRAM_groups.csv")) %>%  as_tibble() %>% 
  filter(NumCohorts >1)
ngroups <- dim(groups_df)[1]

thisVol <- ncvar_get(output, "volume") 
nlayers <- dim(thisVol)[1]; nboxes <- dim(thisVol)[2]; ntimesteps <- dim(thisVol)[3]
dimList <- NULL; dimList[[1]] <- paste0("l",1:nlayers); dimList[[2]] <- paste0("b",1:nboxes); dimList[[3]] <- 1:ntimesteps

####################################
## get age of maturity from biol.prm file (by age, I mean 'cohort')
#####################################
age_mature_prep <- biol_prm[grep("age_mat", biol_prm)] %>% 
  as_tibble() %>% 
  separate(value, into = c("Code", "x1","x2"), sep = "_", remove=FALSE)
age_mature_prep$age_mat <- unlist(lapply(age_mature_prep$value, get_first_number))
age_mature <- age_mature_prep %>%  
  select(c("Code","age_mat")) %>% 
  mutate(Cohort_mature = age_mat + 1) # allow for starting at zero in age_mat

########################################
## get biomass from model outputs ##
######################################
## keep spatial structure
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

groups_biomass <- bind_rows(biomass_list) %>% 
  mutate(Cohort = as.character(Cohort),
         ts = as.integer(ts)) %>% 
  # rename("Biomass" = value) %>% 
  left_join(age_mature, by = "Code") %>% 
  mutate("Stage" = case_when(
    Cohort < Cohort_mature ~ "Juvenile",
    TRUE ~ "Adult"
  )) %>% 
  group_by(layer, box, ts, Name, Code, Stage) %>% 
  summarise("Biomass" = sum(value, na.rm=TRUE))

## without stage, just species
groups_biomassNOstage <- bind_rows(biomass_list) %>% 
  mutate(ts = as.integer(ts)) %>% 
  # rename("Biomass" = value) %>% 
  left_join(age_mature, by = "Code") %>% 
  group_by(layer, box, ts, Name, Code) %>% 
  summarise("Biomass" = sum(value, na.rm=TRUE))

# for catches, need all timesteps but not the spatial detail
biomassByTimestep <- groups_biomass %>% 
  group_by(Name, Code, Stage, ts) %>% 
  summarise(Biomass = sum(Biomass, na.rm=TRUE))

# summ over cohorts, layers and boxes, to get total biomass by species and timestep
biomass_bySpeciesTimestep <-  groups_biomass %>% 
  group_by(Code, ts) %>% 
  summarise(Biomass = sum(Biomass, na.rm=TRUE))
# B0 approx using timestep_4initial defined near the start
initial_byGroup <- biomass_bySpeciesTimestep %>% 
  filter(ts %in% timestep_4initial) %>% 
  group_by(Code) %>% 
  summarise(B0 = mean(Biomass, na.rm=TRUE))

biomass4spatialOverlap <- groups_biomassNOstage %>% 
  filter(ts %in% timestep_4speicesOverlap) 

########################################
## get species group spatial overlap by cohort ##
######################################
# for each pair of species groups, how much of their biomass overlaps spatially and temporally
# approximating with correlation over space and time...
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

########################################
## get catches, then combine with biomass and approximate effort as catch/biomass ##
######################################
groups_catch <- catches %>% 
  mutate(ts = Time/365+1) %>% 
  select(c(ts, groups_df$Code)) %>% 
  pivot_longer(cols = groups_df$Code) %>% 
  rename(Catch = value, Code = name)

# join with biomass, summarise over boxes, so its species, cohort, timestep
groups_biomass_catch <- groups_catch %>% 
  left_join(
    biomass_bySpeciesTimestep, by = c("Code","ts")
  ) %>% 
  mutate("Effort" = Catch/Biomass) %>% 
  rename(Timestep = ts)

##########################################
## diets ##
########################################
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


##########################################
## # average weight of individuals by Stage, weight at maturity, also w_inf - take as 95th percentile ##
########################################
weights_list <- NULL
for(g in 1:ngroups){
  this_numCohorts <- groups_df$NumCohorts[g]
  for(c in 1:this_numCohorts){
    thisVar <- groups_df$Name[g] %>%  str_trim(side="both") %>%  paste0(c)
    
    weight <-   (ncvar_get(nc = output, paste0(thisVar, "_ResN")) + ncvar_get(nc = output, paste0(thisVar, "_StructN"))) * mg_2_tonne * X_CN  
    
    nums <- ncvar_get(nc = output, paste0(thisVar, "_Nums")) 
    
    dimnames(weight) <- dimList; dimnames(nums) <- dimList
    avg_weights <- weight %>%  as_tibble() %>% 
      mutate(layer = paste0("l",1:nlayers)) %>% 
      pivot_longer(cols=!contains("layer")) %>% 
      separate(name, into=c("box","ts")) %>% 
      mutate("Name" = groups_df$Name[g], "Code" = groups_df$Code[g], "Cohort" = c) %>% 
      filter(value > 0) %>% 
      rename(Weight = value)
    
    total_nums <- nums %>%  as_tibble() %>% 
      mutate(layer = paste0("l",1:nlayers)) %>% 
      pivot_longer(cols=!contains("layer")) %>% 
      separate(name, into=c("box","ts")) %>% 
      mutate("Name" = groups_df$Name[g], "Code" = groups_df$Code[g], "Cohort" = c) %>% 
      filter(value > 0) %>% 
      rename(Number = value) %>% 
      left_join(avg_weights, by = c("layer", "box", "ts", "Name", "Code", "Cohort")) %>% 
      filter(Number > 1e-16)
    
    weights_list[[(c-1)*ngroups + g]] <- total_nums 
  }
}

w_infs <- bind_rows(weights_list) %>% 
  group_by(Code, Name) %>% 
  summarise(w_inf = quantile(Weight, 0.95))

avg_weights <-  bind_rows(weights_list) %>% 
  mutate(Cohort = as.character(Cohort),
         ts = as.integer(ts)) %>% 
  left_join(age_mature, by = "Code") %>% 
  mutate("Stage" = case_when(
    Cohort < Cohort_mature ~ "Juvenile",
    TRUE ~ "Adult"
  )) 
w_mature <- avg_weights %>% 
  filter(Cohort == Cohort_mature & Number > 1e-12) %>% 
  group_by(Code) %>% 
  summarise(w_mat = mean(Weight, na.rm=TRUE))
w_mature %>% 
  ggplot(aes(x =  Code, y = w_mat))+
  geom_bar(stat="identity")

#### bring together species_params to output ###
interaction_matrix <- speciesOverlap %>% 
  mutate(resource = paste0("interaction_", Code_B)) %>% 
  select(Code_A, resource, Overlap) %>% 
  pivot_wider(names_from = resource, values_from = Overlap)

species_params <- groups_df %>% 
  select(Code) %>% 
  left_join(
    w_infs, by = "Code"
  ) %>% 
  rename(Species = Code) %>% 
  left_join(w_mature, by = c("Species" = "Code")) %>% 
  left_join(interaction_matrix, by = c("Species" = "Code_A")) 

#######################################
## save out ##
save(list=c("species_params", "diets_proportions", "groups_biomass_catch"), file= file.path(mizerModelPath, "AtlantisPars4Mizer"))


#######################################################################
#######################################
## plots for testing
# diets_proportions %>% 
#   ggplot(aes(x = Predator, y = Diets_proportion, fill = Prey))+
#   geom_bar(stat="identity")+
#   facet_wrap(~Stage)+
#   coord_flip()
# 
# speciesOverlap %>% 
#   filter(Code_A=="ASQ") %>% 
#   ggplot(aes(x = Code_B, y =  Overlap, fill = Stage_B)) +
#   geom_bar(stat="identity")
# 
# 
# speciesOverlap %>% 
#   ggplot(aes(x = Code_A, y = Overlap, fill = Stage_A))+
#   geom_boxplot(outlier.shape = NA)+
#   theme(axis.text.x = element_text(angle=90))+
#   ylim(0.2,0.8)
# 
# avg_weights %>% 
#   filter(Code=="HOK") %>% 
#   ggplot(aes(x = Weight, color = Cohort))+
#   geom_freqpoly()
# avg_weights %>% 
#   filter(Code=="HOK") %>% 
#   mutate(ts = as.character(ts)) %>% 
#   ggplot(aes(x = ts, y = Weight))+
#   geom_boxplot()+
#   facet_wrap(~Cohort)
# 
# avg_weights %>% 
#   filter(Code=="HOK" & Cohort == "10") %>% 
#   ggplot(aes(x = Weight)) +
#   geom_histogram()
# 
# avg_weights %>% 
#   filter(Code=="HOK" & Cohort == "10" & ts>120) %>% 
#   ggplot(aes(x = Weight)) +
#   geom_histogram()
# 
# groups_biomass_catch %>% 
#   filter(Code=="HOK") %>% 
#   group_by(Stage, ts) %>% 
#   filter(ts>120) %>% 
#   summarise(Catch = mean(Catch),
#             Biomass = mean(Biomass))
# groups_biomass_catch %>% 
#   filter(Code=="HOK") %>%  
#   ggplot(aes(x = ts, y = Effort))+
#   geom_line()
#
# biomass_bySpeciesTimestep %>%
#   ggplot(aes(x = ts, y = Biomass, color = Code)) +
#   geom_line(stat="identity")
# groups_biomass %>%
# filter(grepl("^A|^B", Code)) %>%
# ggplot(aes(x = ts, y = Biomass, fill = Cohort)) +
# geom_bar(stat="identity") +
# facet_wrap(~Code)
# 
# # check spatial overlap makes sense
# spatial2plot <- biomass4spatialOverlap %>% 
#   filter(ts %in% timestep_4speicesOverlap) %>% 
#   select(layer, box, ts, Code, Cohort, Biomass) %>% 
#   group_by(layer, box, Code, Cohort) %>% 
#   summarise(Biomass = mean(Biomass, na.rm=TRUE)) %>% 
#   mutate(Cohort = as.integer(Cohort),
#          box_number = as.integer(gsub( "b","", box)))
# spatial2plot %>% 
#   filter(Code %in% c("BEE","LIN") & Biomass >0) %>% 
#   ggplot(aes(x = box_number, y = layer, fill = log(Biomass)))+
#   geom_tile() +
#   facet_grid( Cohort ~ Code, scales="free")
# 
# speciesOverlap %>% 
#   filter(Overlap<0.45 & Code_A !="ASQ")
#   
#   filter(Code_A=="SPD" & Code_B == "HOK")
# 
############
