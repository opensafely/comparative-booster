
library('tidyverse')
library('arrow')
library('here')
library('glue')

#source(here("analysis", "lib", "utility_functions.R"))

remotes::install_github("https://github.com/wjchulme/dd4d")
library('dd4d')


population_size <- 20000

# get nth largest value from list
nthmax <- function(x, n=1){
  dplyr::nth(sort(x, decreasing=TRUE), n)
}

nthmin <- function(x, n=1){
  dplyr::nth(sort(x, decreasing=FALSE), n)
}


source(here("lib", "design", "design.R"))


studystart_date <- as.Date(study_dates$studystart_date)
studyend_date <- as.Date(study_dates$studyend_date)
followupend_date <- as.Date(study_dates$followup_date)
index_date <- studystart_date

firstpfizer_date <- as.Date(study_dates$firstpfizer_date)
firstaz_date <- as.Date(study_dates$firstaz_date)
firstmoderna_date <- as.Date(study_dates$firstmoderna_date)

index_day <- 0L
studystart_day <- as.integer(studystart_date - index_date)
studyend_day <- as.integer(studyend_date - index_date)
firstpfizer_day <- as.integer(firstpfizer_date - index_date)
firstaz_day <- as.integer(firstaz_date - index_date)
firstmoderna_day <- as.integer(firstmoderna_date - index_date)


known_variables <- c(
  "index_date", "studystart_date", "studyend_date", "firstpfizer_date", "firstaz_date", "firstmoderna_date",
  "index_day",  "studystart_day", "studyend_day", "firstpfizer_day", "firstaz_day", "firstmoderna_day"
)

sim_list = lst(

  dereg_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+120)),
    missing_rate = ~0.99
  ),

  has_follow_up_previous_6weeks = bn_node(
    ~rbernoulli(n=..n, p=0.999)
  ),

  hscworker = bn_node(
    ~rbernoulli(n=..n, p=0.01)
  ),

  age = bn_node(
    ~as.integer(rnorm(n=..n, mean=60, sd=14))
  ),

  age_august2021 = bn_node(~age),

  sex = bn_node(
    ~rfactor(n=..n, levels = c("F", "M"), p = c(0.51, 0.49)),
    missing_rate = ~0.001 # this is shorthand for ~(rbernoulli(n=..n, p = 0.2))
  ),

  bmi = bn_node(
    ~rfactor(n=..n, levels = c("Not obese", "Obese I (30-34.9)", "Obese II (35-39.9)", "Obese III (40+)"), p = c(0.5, 0.2, 0.2, 0.1)),
  ),

  ethnicity = bn_node(
    ~rfactor(n=..n, levels = c(1,2,3,4,5), p = c(0.8, 0.05, 0.05, 0.05, 0.05)),
    missing_rate = ~ 0.25
  ),

  ethnicity_6_sus = bn_node(
    ~rfactor(n=..n, levels = c(0,1,2,3,4,5), p = c(0.1, 0.7, 0.05, 0.05, 0.05, 0.05)),
    missing_rate = ~ 0
  ),

  practice_id = bn_node(
    ~as.integer(runif(n=..n, 1, 200))
  ),

  msoa = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 100)), levels=1:100),
    missing_rate = ~ 0.005
  ),

  stp = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 36)), levels=1:36)
  ),

  region = bn_node(
    variable_formula = ~rfactor(n=..n, levels=c(
      "North East",
      "North West",
      "Yorkshire and The Humber",
      "East Midlands",
      "West Midlands",
      "East",
      "London",
      "South East",
      "South West"
    ), p = c(0.2, 0.2, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05))
  ),

  imd = bn_node(
    ~factor(plyr::round_any(runif(n=..n, 1, 32000), 100), levels=seq(0,32000,100)),
    missing_rate = ~0.02
  ),

  imd_integer = bn_node(
    ~as.integer(as.character(imd)),
    keep=FALSE
  ),

  imd_Q5 = bn_node(
    ~factor(
      case_when(
        (imd_integer >= 0) & (imd_integer < 32844*1/5) ~ "1 (most deprived)",
        (imd_integer >= 32844*1/5) & (imd_integer < 32844*2/5) ~ "2",
        (imd_integer >= 32844*2/5) & (imd_integer < 32844*3/5) ~ "3",
        (imd_integer >= 32844*3/5) & (imd_integer < 32844*4/5) ~ "4",
        (imd_integer >= 32844*4/5) & (imd_integer <= 32844*5/5) ~ "5 (least deprived)",
        TRUE ~ "Unknown"
      ),
      levels= c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown")
    ),
    missing_rate = ~0
  ),

  rural_urban = bn_node(
    ~rfactor(n=..n, levels = 1:9, p = rep(1/9, 9)),
    missing_rate = ~ 0.1
  ),

  care_home_tpp = bn_node(
    ~rbernoulli(n=..n, p = 0.01)
  ),

  care_home_code = bn_node(
    ~rbernoulli(n=..n, p = 0.01)
  ),

  ## vaccination variables

  first_vax_type = bn_node(~rcat(n=..n, c("pfizer","az"), c(0.5,0.5)), keep=FALSE),
  third_vax_type = bn_node(~rcat(n=..n, c("pfizer","moderna"), c(0.5,0.5)), keep=FALSE),
  covid_vax_pfizer_1_day = bn_node(
    ~as.integer(runif(n=..n, firstpfizer_day, firstpfizer_day+60)),
    missing_rate = ~1-(first_vax_type=="pfizer")
  ),
  covid_vax_pfizer_2_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_pfizer_1_day+30, covid_vax_pfizer_1_day+60)),
    needs = c("covid_vax_pfizer_1_day"),
  ),
  covid_vax_pfizer_3_day = bn_node(
    ~as.integer(runif(n=..n, pmax(pmin(covid_vax_pfizer_2_day+15, covid_vax_az_2_day+15), studystart_day), pmax(pmin(covid_vax_pfizer_2_day+15, covid_vax_az_2_day+15), studystart_day)+100)),
    needs = c("covid_vax_pfizer_1_day", "covid_vax_pfizer_2_day"),
    missing_rate = ~1-(third_vax_type=="pfizer")
  ),
  covid_vax_pfizer_4_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_pfizer_3_day+120, covid_vax_pfizer_3_day+200)),
    needs = c("covid_vax_pfizer_1_day", "covid_vax_pfizer_2_day", "covid_vax_pfizer_3_day"),
    missing_rate = ~0.99
  ),
  covid_vax_az_1_day = bn_node(
    ~as.integer(runif(n=..n, firstaz_day, firstaz_day+60)),
    missing_rate = ~1-(first_vax_type=="az")
  ),
  covid_vax_az_2_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_az_1_day+30, covid_vax_az_1_day+60)),
    needs = c("covid_vax_az_1_day"),
  ),
  covid_vax_az_3_day = bn_node(
    ~as.integer(runif(n=..n, pmax(pmin(covid_vax_pfizer_2_day+15, covid_vax_az_2_day+15), studystart_day), pmax(pmin(covid_vax_pfizer_2_day+15, covid_vax_az_2_day+15), studystart_day)+100)),
    needs = c("covid_vax_az_1_day", "covid_vax_az_2_day")
  ),
  covid_vax_az_4_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_az_3_day+120, covid_vax_az_3_day+200)),
    needs = c("covid_vax_az_1_day", "covid_vax_az_2_day", "covid_vax_az_3_day"),
    missing_rate = ~0.99
  ),
  covid_vax_moderna_1_day = bn_node(
    ~as.integer(runif(n=..n, studystart_day, studystart_day+60)),
    missing_rate = ~1-(third_vax_type=="moderna")
  ),
  covid_vax_moderna_2_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_moderna_1_day+30, covid_vax_moderna_1_day+60)),
    needs = c("covid_vax_moderna_1_day")
  ),
  covid_vax_moderna_3_day = bn_node(
    ~as.integer(runif(n=..n, pmax(covid_vax_moderna_2_day+15, studystart_day), pmax(covid_vax_moderna_2_day, studystart_day)+100)),
    needs = c("covid_vax_moderna_1_day", "covid_vax_moderna_2_day")
  ),
  covid_vax_moderna_4_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_moderna_3_day+120, covid_vax_moderna_3_day+200)),
    needs = c("covid_vax_moderna_1_day", "covid_vax_moderna_2_day", "covid_vax_moderna_3_day"),
    missing_rate = ~0.99

  ),

  # covid vax any
  anycovidvax_1_day = bn_node(
    ~pmin(covid_vax_pfizer_1_day, covid_vax_az_1_day, na.rm=TRUE),
  ),
  anycovidvax_2_day = bn_node(
    ~pmin(covid_vax_pfizer_2_day, covid_vax_az_2_day, na.rm=TRUE),
  ),
  anycovidvax_3_day = bn_node(
    ~pmin(covid_vax_pfizer_3_day, covid_vax_az_3_day, covid_vax_moderna_1_day, covid_vax_moderna_3_day, na.rm=TRUE),
  ),
  anycovidvax_4_day = bn_node(
    ~pmin(covid_vax_pfizer_4_day, covid_vax_az_4_day, covid_vax_moderna_4_day, na.rm=TRUE),
  ),



  ## baseline clinical variables

  asthma = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_neuro_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_resp_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  sev_obesity = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  diabetes = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  sev_mental = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_heart_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_kidney_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_liver_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),

  immunosuppressed = bn_node( ~immdx | immrx),
  immdx = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  immrx = bn_node( ~rbernoulli(n=..n, p = 0.02)),

  asplenia = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  cancer_haem = bn_node( ~rbernoulli(n=..n, p = 0.01)),
  cancer_nonhaem = bn_node( ~rbernoulli(n=..n, p = 0.01)),
  solid_organ_transplant = bn_node( ~rbernoulli(n=..n, p = 0.01)),
  hiv_aids = bn_node( ~rbernoulli(n=..n, p = 0.01)),

  learndis = bn_node( ~rbernoulli(n=..n, p = 0.02)),

  cev_ever = bn_node( ~rbernoulli(n=..n, p = 0.05)),
  cev = bn_node( ~rbernoulli(n=..n, p = 0.05)),

  endoflife = bn_node( ~rbernoulli(n=..n, p = 0.001)),
  housebound = bn_node( ~rbernoulli(n=..n, p = 0.001)),

  prior_covid_test_frequency = bn_node(
    ~as.integer(rpois(n=..n, lambda=3)),
    missing_rate = ~0
  ),

  inhospital = bn_node( ~rbernoulli(n=..n, p = 0.01)),

  ## pre-baseline events where event date is relevant

  primary_care_covid_case_0_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day-100, anycovidvax_3_day-1)),
    missing_rate = ~0.7
  ),

  covid_test_0_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day-100, anycovidvax_3_day-1)),
    missing_rate = ~0.7
  ),

  postest_0_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day-100, anycovidvax_3_day-1)),
    missing_rate = ~0.9
  ),

  covidemergency_0_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day-100, anycovidvax_3_day-1)),
    missing_rate = ~0.99
  ),


  covidadmitted_0_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day-100, anycovidvax_3_day-1)),
    missing_rate = ~0.99
  ),

  ## post-baseline events (outcomes)

  primary_care_covid_case_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.7
  ),

  covid_test_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.7
  ),

  postest_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.7
  ),

  emergency_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.8
  ),
  emergencyhosp_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.85
  ),


  covidemergency_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.8
  ),

  covidemergencyhosp_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.85
  ),

  # respemergency_day = bn_node(
  #   ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
  #   missing_rate = ~0.95
  # ),
  #
  # respemergencyhosp_day = bn_node(
  #   ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
  #   missing_rate = ~0.95
  # ),

  covidadmitted_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.7
  ),

  covidcritcare_day = bn_node(
    ~covidadmitted_day,
    needs = "covidadmitted_day",
    missing_rate = ~0.7
  ),


  admitted_unplanned_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.7
  ),
  admitted_planned_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.7
  ),

  coviddeath_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),

  death_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.90
  ),

  # fractures
  fractureemergency_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.95
  ),

  fractureadmitted_day = bn_node(
    ~as.integer(runif(n=..n, anycovidvax_3_day, anycovidvax_3_day+100)),
    missing_rate = ~0.97
  ),

  fracturedeath_day = bn_node(
    ~death_day,
    missing_rate = ~0.95,
    needs = "death_day"
  ),

  test_count = bn_node(
    ~ rpois(n = ..n, 1)
  ),

  postest_count = bn_node(
    ~ rpois(n = ..n, 0.1)
  )


)
bn <- bn_create(sim_list, known_variables = known_variables)

bn_plot(bn)
bn_plot(bn, connected_only=TRUE)

set.seed(10)

dummydata <-bn_simulate(bn, pop_size = population_size, keep_all = FALSE, .id="patient_id")


dummydata_processed <- dummydata %>%
  # change index date from fixed to third / booster dose date

  mutate(
    across(c(

    ))
  ) %>%

  mutate(

  ) %>%
  #convert logical to integer as study defs output 0/1 not TRUE/FALSE
  #mutate(across(where(is.logical), ~ as.integer(.))) %>%
  #convert integer days to dates since index date and rename vars
  mutate(across(ends_with("_day"), ~ as.Date(as.character(index_date + .)))) %>%
  rename_with(~str_replace(., "_day", "_date"), ends_with("_day"))


fs::dir_create(here("lib", "dummydata"))
write_feather(dummydata_processed, sink = here("lib", "dummydata", "dummyinput.feather"))


