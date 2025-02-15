
from cohortextractor import (
  StudyDefinition,
  patients,
  codelist_from_csv,
  codelist,
  filter_codes_by_category,
  combine_codelists,
)

# Import codelists from codelists.py
import codelists

# import json module
import json

# import study dates defined in "./lib/design/study-dates.R" script
with open("./lib/design/study-dates.json") as f:
  study_dates = json.load(f)

# change these in design.R if necessary
studystart_date = study_dates["studystart_date"] 
studyend_date = study_dates["studyend_date"]
followupend_date = study_dates["followupend_date"]
postestfollowupend_date = study_dates["postestfollowupend_date"]
firstpfizer_date = study_dates["firstpfizer_date"]
firstaz_date = study_dates["firstaz_date"]
firstmoderna_date = study_dates["firstmoderna_date"]


## Functions for extracting a series of time dependent variables
# These define study defintion variable signatures such that
# variable_1_date is the the first event date on or after the index date
# variable_2_date is the first event date strictly after variable_2_date
# ...
# variable_n_date is the first event date strictly after variable_n-1_date


def vaccination_date_X(name, index_date, n, product_name_matches=None, target_disease_matches=None):
  # vaccination date, given product_name
  def var_signature(
    name,
    on_or_after,
    product_name_matches,
    target_disease_matches
  ):
    return {
      name: patients.with_tpp_vaccination_record(
        product_name_matches=product_name_matches,
        target_disease_matches=target_disease_matches,
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
      ),
    }
    
  variables = var_signature(f"{name}_1_date", index_date, product_name_matches, target_disease_matches)
  for i in range(2, n+1):
    variables.update(var_signature(
      f"{name}_{i}_date", 
      f"{name}_{i-1}_date + 1 days",
      # pick up subsequent vaccines occurring one day or later -- people with unrealistic dosing intervals are later excluded
      product_name_matches,
      target_disease_matches
    ))
  return variables



# Specify study defeinition
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "2020-01-01", "latest": studyend_date},
    "rate": "uniform",
    "incidence": 0.2,
    "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
    "float": {"distribution": "normal", "mean": 25, "stddev": 5},
  },
  
  #index_date = index_date,
  
  # This line defines the study population
  population=patients.satisfying(
    f"""
      registered
      AND
      age_august2021 >= 18
      AND
      NOT has_died
      AND 
      anycovidvax_3_date >= startdate
      AND
      anycovidvax_3_date <= enddate
    """,
    # we define baseline variables on the day _before_ the study date (start date = day of first possible booster vaccination)
    registered=patients.registered_as_of(
      "anycovidvax_3_date - 1 day",
    ),
    has_died=patients.died_from_any_cause(
      on_or_before="anycovidvax_3_date - 1 day",
      returning="binary_flag",
    ), 
   
    startdate = patients.fixed_value(studystart_date),
    enddate = patients.fixed_value(studyend_date),
  
    
  ),
  
  
  #################################################################
  ## Covid vaccine dates
  #################################################################
  
  # pfizer
  **vaccination_date_X(
    name = "covid_vax_pfizer",
    # use 1900 to capture all possible recorded covid vaccinations, including date errors
    # any vaccines occurring before national rollout are later excluded
    index_date = "1900-01-01", 
    n = 4,
    product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)"
  ),
  
  # az
  **vaccination_date_X(
    name = "covid_vax_az",
    index_date = "1900-01-01",
    n = 4,
    product_name_matches="COVID-19 Vaccine Vaxzevria 0.5ml inj multidose vials (AstraZeneca)"
  ),
  
  # moderna
  **vaccination_date_X(
    name = "covid_vax_moderna",
    index_date = "1900-01-01",
    n = 4,
    product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)"
  ),
  
  # any covid vaccine
  **vaccination_date_X(
    name = "anycovidvax",
    index_date = "1900-01-01",
    n = 4,
    target_disease_matches="SARS-2 CORONAVIRUS"
  ),
  
  
  
  
  ###############################################################################
  ## Admin and demographics
  ###############################################################################
  
  has_follow_up_previous_6weeks=patients.registered_with_one_practice_between(
    start_date="anycovidvax_3_date - 42 days",
    end_date="anycovidvax_3_date",
  ),
  
  dereg_date=patients.date_deregistered_from_all_supported_practices(
    on_or_after="anycovidvax_3_date",
    date_format="YYYY-MM-DD",
  ),

  
  age=patients.age_as_of( 
    "anycovidvax_3_date - 1 day",
  ),
  
  # for jcvi group definitions
  age_august2021=patients.age_as_of( 
    "2020-08-31",
  ),
  
  sex=patients.sex(
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
      "incidence": 1,
    }
  ),
  
  # https://github.com/opensafely/risk-factors-research/issues/51
  bmi=patients.categorised_as(
    {
      "Not obese": "DEFAULT",
      "Obese I (30-34.9)": """ bmi_value >= 30 AND bmi_value < 35""",
      "Obese II (35-39.9)": """ bmi_value >= 35 AND bmi_value < 40""",
      "Obese III (40+)": """ bmi_value >= 40 AND bmi_value < 100""",
      # set maximum to avoid any impossibly extreme values being classified as obese
    },
    bmi_value=patients.most_recent_bmi(
      on_or_after="anycovidvax_3_date - 5 years",
      minimum_age_at_measurement=16
    ),
    return_expectations={
      "rate": "universal",
      "category": {
        "ratios": {
          "Not obese": 0.7,
          "Obese I (30-34.9)": 0.1,
          "Obese II (35-39.9)": 0.1,
          "Obese III (40+)": 0.1,
        }
      },
    },
  ),
  

  # Ethnicity in 6 categories
  ethnicity = patients.with_these_clinical_events(
    codelists.ethnicity,
    returning="category",
    find_last_match_in_period=True,
    include_date_of_match=False,
    return_expectations={
      "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
      "incidence": 0.75,
    },
  ),
  
  # ethnicity variable that takes data from SUS
  ethnicity_6_sus = patients.with_ethnicity_from_sus(
    returning="group_6",  
    use_most_frequent_code=True,
    return_expectations={
      "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
      "incidence": 0.8,
    },
  ),
  
  ################################################################################################
  ## Practice and patient ID variables
  ################################################################################################
  # practice pseudo id
  practice_id=patients.registered_practice_as_of(
    "anycovidvax_3_date - 1 day",
    returning="pseudo_id",
    return_expectations={
      "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
      "incidence": 1,
    },
  ),
  
  # msoa
  msoa=patients.address_as_of(
    "anycovidvax_3_date - 1 day",
    returning="msoa",
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {"E02000001": 0.0625, "E02000002": 0.0625, "E02000003": 0.0625, "E02000004": 0.0625,
        "E02000005": 0.0625, "E02000007": 0.0625, "E02000008": 0.0625, "E02000009": 0.0625, 
        "E02000010": 0.0625, "E02000011": 0.0625, "E02000012": 0.0625, "E02000013": 0.0625, 
        "E02000014": 0.0625, "E02000015": 0.0625, "E02000016": 0.0625, "E02000017": 0.0625}},
    },
  ),    
  
  # stp is an NHS administration region based on geography
  stp=patients.registered_practice_as_of(
    "anycovidvax_3_date - 1 day",
    returning="stp_code",
    return_expectations={
      "rate": "universal",
      "category": {
        "ratios": {
          "STP1": 0.1,
          "STP2": 0.1,
          "STP3": 0.1,
          "STP4": 0.1,
          "STP5": 0.1,
          "STP6": 0.1,
          "STP7": 0.1,
          "STP8": 0.1,
          "STP9": 0.1,
          "STP10": 0.1,
        }
      },
    },
  ),
  # NHS administrative region
  region=patients.registered_practice_as_of(
    "anycovidvax_3_date - 1 day",
    returning="nuts1_region_name",
    return_expectations={
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.2,
          "East Midlands": 0.1,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.1,
          "South East": 0.1,
          "South West": 0.1
          #"" : 0.01
        },
      },
    },
  ),
  
  ## IMD - quintile
  
  imd=patients.address_as_of(
    "anycovidvax_3_date - 1 day",
    returning="index_of_multiple_deprivation",
    round_to_nearest=100,
    return_expectations={
      "category": {"ratios": {c: 1/320 for c in range(100, 32100, 100)}}
    }
  ),
  
  imd_Q5=patients.categorised_as(
    {
      "Unknown": "DEFAULT",
      "1 (most deprived)": "imd >= 0 AND imd < 32844*1/5",
      "2": "imd >= 32844*1/5 AND imd < 32844*2/5",
      "3": "imd >= 32844*2/5 AND imd < 32844*3/5",
      "4": "imd >= 32844*3/5 AND imd < 32844*4/5",
      "5 (least deprived)": "imd >= 32844*4/5 AND imd <= 32844",
    },
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {"Unknown": 0.02, "1 (most deprived)": 0.18, "2": 0.2, "3": 0.2, "4": 0.2, "5 (least deprived)": 0.2}},
    },
  ),
  
  #rurality
  rural_urban=patients.address_as_of(
    "anycovidvax_3_date - 1 day",
    returning="rural_urban_classification",
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {1: 0.125, 2: 0.125, 3: 0.125, 4: 0.125, 5: 0.125, 6: 0.125, 7: 0.125, 8: 0.125}},
    },
  ),



  ################################################################################################
  ## occupation / residency
  ################################################################################################


  # health or social care worker  
  hscworker = patients.with_healthcare_worker_flag_on_covid_vaccine_record(returning="binary_flag"),
  
  care_home_tpp=patients.satisfying(
    "care_home='1'",
    
    care_home = patients.care_home_status_as_of(
      "anycovidvax_3_date - 1 day",
      categorised_as={
          "1": "IsPotentialCareHome",
          "": "DEFAULT",  # use empty string
      },
    ),
  ),
  
  
  # Patients in long-stay nursing and residential care
  care_home_code=patients.with_these_clinical_events(
      codelists.carehome,
      on_or_before="anycovidvax_3_date - 1 day",
      returning="binary_flag",
      return_expectations={"incidence": 0.01},
  ),
  
  

  ################################################################################################
  ## Pre-baseline events where event date is of interest
  ################################################################################################


  # Positive case identification prior to study start date
  primary_care_covid_case_0_date=patients.with_these_clinical_events(
    combine_codelists(
      codelists.covid_primary_care_code,
      codelists.covid_primary_care_positive_test,
      codelists.covid_primary_care_sequelae,
    ),
    returning="date",
    date_format="YYYY-MM-DD",
    on_or_before="anycovidvax_3_date - 1 day",
    find_last_match_in_period=True,
  ),
  
  # covid PCR test dates from SGSS
  covid_test_0_date=patients.with_test_result_in_sgss(
    pathogen="SARS-CoV-2",
    test_result="any",
    on_or_before="anycovidvax_3_date - 1 day",
    returning="date",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
    restrict_to_earliest_specimen_date=False,
  ),

  
  # positive covid test
  postest_0_date=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="positive",
      returning="date",
      date_format="YYYY-MM-DD",
      on_or_before="anycovidvax_3_date - 1 day",
      find_last_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
  ),
  
  # emergency attendance for covid
  covidemergency_0_date=patients.attended_emergency_care(
    returning="date_arrived",
    on_or_before="anycovidvax_3_date - 1 day",
    with_these_diagnoses = codelists.covid_emergency,
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),
  
    # Positive covid admission prior to study start date
  covidadmitted_0_date=patients.admitted_to_hospital(
    returning="date_admitted",
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_these_diagnoses=codelists.covid_icd10,
    on_or_before="anycovidvax_3_date - 1 day",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),
  

  


  ############################################################
  ## Clinical information as at (day before) 3rd / booster dose date
  ############################################################
  
  
  # From PRIMIS

  asthma = patients.satisfying(
    """
      astadm OR
      (ast AND astrxm1 AND astrxm2 AND astrxm3)
      """,
    # Asthma Admission codes
    astadm=patients.with_these_clinical_events(
      codelists.astadm,
      returning="binary_flag",
      on_or_before="anycovidvax_3_date - 1 day",
    ),
    # Asthma Diagnosis code
    ast = patients.with_these_clinical_events(
      codelists.ast,
      returning="binary_flag",
      on_or_before="anycovidvax_3_date - 1 day",
    ),
    # Asthma systemic steroid prescription code in month 1
    astrxm1=patients.with_these_medications(
      codelists.astrx,
      returning="binary_flag",
      between=["anycovidvax_3_date - 30 days", "anycovidvax_3_date - 1 day"],
    ),
    # Asthma systemic steroid prescription code in month 2
    astrxm2=patients.with_these_medications(
      codelists.astrx,
      returning="binary_flag",
      between=["anycovidvax_3_date - 60 days", "anycovidvax_3_date - 31 days"],
    ),
    # Asthma systemic steroid prescription code in month 3
    astrxm3=patients.with_these_medications(
      codelists.astrx,
      returning="binary_flag",
      between= ["anycovidvax_3_date - 90 days", "anycovidvax_3_date - 61 days"],
    ),

  ),

  # Chronic Neurological Disease including Significant Learning Disorder
  chronic_neuro_disease=patients.with_these_clinical_events(
    codelists.cns_cov,
    returning="binary_flag",
    on_or_before="anycovidvax_3_date - 1 day",
  ),

  # Chronic Respiratory Disease
  chronic_resp_disease = patients.satisfying(
    "asthma OR resp_cov",
    resp_cov=patients.with_these_clinical_events(
      codelists.resp_cov,
      returning="binary_flag",
      on_or_before="anycovidvax_3_date - 1 day",
    ),
  ),

  sev_obesity = patients.satisfying(
    """
      sev_obesity_date > bmi_date OR
      bmi_value1 >= 40
      """,

    bmi_stage_date=patients.with_these_clinical_events(
      codelists.bmi_stage,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    sev_obesity_date=patients.with_these_clinical_events(
      codelists.sev_obesity,
      returning="date",
      find_last_match_in_period=True,
      ignore_missing_values=True,
      between= ["bmi_stage_date", "anycovidvax_3_date - 1 day"],
      date_format="YYYY-MM-DD",
    ),

    bmi_date=patients.with_these_clinical_events(
      codelists.bmi,
      returning="date",
      ignore_missing_values=True,
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    bmi_value1=patients.with_these_clinical_events(
      codelists.bmi,
      returning="numeric_value",
      ignore_missing_values=True,
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
    ),

  ),

  diabetes = patients.satisfying(
    "(dmres_date < diab_date) OR (diab_date AND (NOT dmres_date))",
    
    diab_date=patients.with_these_clinical_events(
      codelists.diab,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    dmres_date=patients.with_these_clinical_events(
      codelists.dmres,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
  ),

  sev_mental=patients.satisfying(
    "(smhres_date < sev_mental_date) OR (sev_mental_date AND (NOT smhres_date))",

    # Severe Mental Illness codes
    sev_mental_date=patients.with_these_clinical_events(
      codelists.sev_mental,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
    # Remission codes relating to Severe Mental Illness
    smhres_date=patients.with_these_clinical_events(
      codelists.smhres,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),
  ),


  # Chronic heart disease codes
  chronic_heart_disease=patients.with_these_clinical_events(
    codelists.chd_cov,
    returning="binary_flag",
    on_or_before="anycovidvax_3_date - 1 day",
  ),

  chronic_kidney_disease=patients.satisfying(
    """
      ckd OR
      (ckd15_date AND ckd35_date >= ckd15_date)
      """,

    # Chronic kidney disease codes - all stages
    ckd15_date=patients.with_these_clinical_events(
      codelists.ckd15,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    # Chronic kidney disease codes-stages 3 - 5
    ckd35_date=patients.with_these_clinical_events(
      codelists.ckd35,
      returning="date",
      find_last_match_in_period=True,
      on_or_before="anycovidvax_3_date - 1 day",
      date_format="YYYY-MM-DD",
    ),

    # Chronic kidney disease diagnostic codes
    ckd=patients.with_these_clinical_events(
      codelists.ckd_cov,
      returning="binary_flag",
      on_or_before="anycovidvax_3_date - 1 day",
    ),
  ),


  # Chronic Liver disease codes
  chronic_liver_disease=patients.with_these_clinical_events(
    codelists.cld,
    returning="binary_flag",
    on_or_before="anycovidvax_3_date - 1 day",
  ),

  ## immunosupressive conditions / medications ####
  #################################################

  # Immunosuppression diagnosis codes
  immdx=patients.with_these_clinical_events(
    codelists.immdx_cov,
    returning="binary_flag",
    on_or_before="anycovidvax_3_date - 1 day",
  ),
    
  # Immunosuppression medication codes
  immrx=patients.with_these_medications(
    codelists.immrx,
    returning="binary_flag",
    between=["anycovidvax_3_date - 182 days", "anycovidvax_3_date - 1 day"]
  ),
  
  immunosuppressed=patients.satisfying(
    "immrx OR immdx",
  ),
  
  # Asplenia or Dysfunction of the Spleen codes
  asplenia=patients.with_these_clinical_events(
    codelists.spln_cov,
    returning="binary_flag",
    on_or_before="anycovidvax_3_date - 1 day",
  ),
  
  # haematological cancer
  cancer_haem = patients.with_these_clinical_events(
    codelists.cancer_haem_snomed,
    returning="binary_flag",
    between=["anycovidvax_3_date - 3 years", "anycovidvax_3_date - 1 day"],
  ),
  
  # solid cancer
  cancer_nonhaem = patients.with_these_clinical_events( 
    codelists.cancer_nonhaem_snomed,
    between=["anycovidvax_3_date - 3 years", "anycovidvax_3_date - 1 day"],
    returning="binary_flag",
  ), 
  
  # solid organ transplant
  solid_organ_transplant = patients.with_these_clinical_events( 
    codelists.solid_organ_transplant,
    on_or_before="anycovidvax_3_date - 1 day",
    returning="binary_flag",
  ), 
  
  # HIV / AIDS
  hiv_aids = patients.with_these_clinical_events( 
    codelists.hiv_aids,
    on_or_before="anycovidvax_3_date - 1 day",
    returning="binary_flag",
  ), 
  
  
  
  # NHSD Haematological diseases
  
  # haematological_disease_nhsd=patients.satisfying(
  #   """
  #   haematopoietic_stem_cell_transplant_nhsd_snomed OR
  #   haematological_malignancies_nhsd_snomed OR
  #   sickle_cell_disease_nhsd_snomed OR
  #   """,
  #   return_expectations={
  #     "incidence": 0.05,
  #   },
  #   haematopoietic_stem_cell_transplant_nhsd_snomed=patients.with_these_clinical_events(
  #     haematopoietic_stem_cell_transplant_nhsd_snomed_codes,
  #     on_or_before="anycovidvax_3_date - 1 day",
  #     returning="binary_flag",
  #   ),
  #   haematological_malignancies_nhsd_snomed=patients.with_these_clinical_events(
  #     haematological_malignancies_nhsd_snomed_codes,
  #     on_or_before="anycovidvax_3_date - 1 day",
  #     returning="binary_flag",
  #   ),
  #   sickle_cell_disease_nhsd_snomed=patients.with_these_clinical_events(
  #     sickle_cell_disease_nhsd_snomed_codes,
  #     on_or_before="anycovidvax_3_date - 1 day",
  #     returning="binary_flag",
  #   ),
  # ),
  


  

  # Wider Learning Disability
  learndis=patients.with_these_clinical_events(
    codelists.learndis,
    returning="binary_flag",
    on_or_before="anycovidvax_3_date - 1 day",
  ),


  # to represent household contact of shielding individual
  # hhld_imdef_dat=patients.with_these_clinical_events(
  #   codelists.hhld_imdef,
  #   returning="date",
  #   find_last_match_in_period=True,
  #   on_or_before="anycovidvax_3_date - 1 day",
  #   date_format="YYYY-MM-DD",
  # ),
  #
  # #####################################
  # # primis employment codelists
  # #####################################
  #
  # # Carer codes
  # carer_date=patients.with_these_clinical_events(
  #   codelists.carer,
  #   returning="date",
  #   find_last_match_in_period=True,
  #   on_or_before="anycovidvax_3_date - 1 day",
  #   date_format="YYYY-MM-DD",
  # ),
  # # No longer a carer codes
  # notcarer_date=patients.with_these_clinical_events(
  #   codelists.notcarer,
  #   returning="date",
  #   find_last_match_in_period=True,
  #   on_or_before="anycovidvax_3_date - 1 day",
  #   date_format="YYYY-MM-DD",
  # ),
  # # Employed by Care Home codes
  # carehome_date=patients.with_these_clinical_events(
  #   codelists.carehomeemployee,
  #   returning="date",
  #   find_last_match_in_period=True,
  #   on_or_before="anycovidvax_3_date - 1 day",
  #   date_format="YYYY-MM-DD",
  # ),
  # # Employed by nursing home codes
  # nursehome_date=patients.with_these_clinical_events(
  #   codelists.nursehomeemployee,
  #   returning="date",
  #   find_last_match_in_period=True,
  #   on_or_before="anycovidvax_3_date - 1 day",
  #   date_format="YYYY-MM-DD",
  # ),
  # # Employed by domiciliary care provider codes
  # domcare_date=patients.with_these_clinical_events(
  #   codelists.domcareemployee,
  #   returning="date",
  #   find_last_match_in_period=True,
  #   on_or_before="anycovidvax_3_date - 1 day",
  #   date_format="YYYY-MM-DD",
  # ),
  
  
  
  #####################################
  # JCVI groups
  #####################################
  
  cev_ever = patients.with_these_clinical_events(
    codelists.shield,
    returning="binary_flag",
    on_or_before = "anycovidvax_3_date - 1 day",
    find_last_match_in_period = True,
  ),

  cev = patients.satisfying(
    """severely_clinically_vulnerable AND NOT less_vulnerable""",
    ##### The shielded patient list was retired in March/April 2021 when shielding ended
    ##### so it might be worth using that as the end date instead of index_date, as we're not sure
    ##### what has happened to these codes since then, e.g. have doctors still been adding new
    ##### shielding flags or low-risk flags? Depends what you're looking for really. Could investigate separately.

    ### SHIELDED GROUP - first flag all patients with "high risk" codes
    severely_clinically_vulnerable = patients.with_these_clinical_events(
      codelists.shield,
      returning="binary_flag",
      on_or_before = "anycovidvax_3_date - 1 day",
      find_last_match_in_period = True,
    ),

    # find date at which the high risk code was added
    date_severely_clinically_vulnerable = patients.date_of(
      "severely_clinically_vulnerable",
      date_format="YYYY-MM-DD",
    ),

    ### NOT SHIELDED GROUP (medium and low risk) - only flag if later than 'shielded'
    less_vulnerable = patients.with_these_clinical_events(
      codelists.nonshield,
      between=["date_severely_clinically_vulnerable + 1 day", "anycovidvax_3_date - 1 day",],
    ),

  ),
  
  endoflife = patients.satisfying(
    """
    midazolam OR
    endoflife_coding
    """,
  
    midazolam = patients.with_these_medications(
      codelists.midazolam,
      returning="binary_flag",
      on_or_before = "anycovidvax_3_date - 1 day",
    ),
    
    endoflife_coding = patients.with_these_clinical_events(
      codelists.eol,
      returning="binary_flag",
      on_or_before = "anycovidvax_3_date - 1 day",
      find_last_match_in_period = True,
    ),
        
  ),
    
  housebound = patients.satisfying(
    """housebound_date
    AND NOT no_longer_housebound
    AND NOT moved_into_care_home
    """,
        
    housebound_date=patients.with_these_clinical_events( 
      codelists.housebound, 
      on_or_before="anycovidvax_3_date - 1 day",
      find_last_match_in_period = True,
      returning="date",
      date_format="YYYY-MM-DD",
    ),   
    no_longer_housebound=patients.with_these_clinical_events( 
      codelists.no_longer_housebound, 
      between=["housebound_date", "anycovidvax_3_date - 1 day"],
    ),
    moved_into_care_home=patients.with_these_clinical_events(
      codelists.carehome,
      between=["housebound_date", "anycovidvax_3_date - 1 day"],
    ),
  ),
  
  prior_covid_test_frequency=patients.with_test_result_in_sgss(
    pathogen="SARS-CoV-2",
    test_result="any",
    between=["anycovidvax_3_date - 182 days", "anycovidvax_3_date - 1 day"], # 182 days = 26 weeks
    returning="number_of_matches_in_period", 
    date_format="YYYY-MM-DD",
    restrict_to_earliest_specimen_date=False,
  ),
  
  # overnight hospital admission at time of 3rd / booster dose
  inhospital = patients.satisfying(
  
    "discharged_0_date >= anycovidvax_3_date",
    
    discharged_0_date=patients.admitted_to_hospital(
      returning="date_discharged",
      on_or_before="anycovidvax_3_date", # this is the admission date
      # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
      # see https://docs.opensafely.org/study-def-variables/#sus for more info
      with_admission_method = ['11', '12', '13', '21', '2A', '22', '23', '24', '25', '2D', '28', '2B', '81'],
      with_patient_classification = ["1"], # ordinary admissions only
      date_format="YYYY-MM-DD",
      find_last_match_in_period=True,
    ), 
  ),
  

  ############################################################
  ## Post-baseline variables (outcomes)
  ############################################################


  # Positive case identification after study start date
  primary_care_covid_case_date=patients.with_these_clinical_events(
    combine_codelists(
      codelists.covid_primary_care_code,
      codelists.covid_primary_care_positive_test,
      codelists.covid_primary_care_sequelae,
    ),
    returning="date",
    date_format="YYYY-MM-DD",
    on_or_after="anycovidvax_3_date",
    find_first_match_in_period=True,
  ),
  
  
  # covid PCR test dates from SGSS
  covid_test_date=patients.with_test_result_in_sgss(
    pathogen="SARS-CoV-2",
    test_result="any",
    on_or_after="anycovidvax_3_date",
    find_first_match_in_period=True,
    restrict_to_earliest_specimen_date=False,
    returning="date",
    date_format="YYYY-MM-DD",
  ),
  
  # positive covid test
  postest_date=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="positive",
      returning="date",
      date_format="YYYY-MM-DD",
      on_or_after="anycovidvax_3_date",
      find_first_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
  ),
  
  # emergency attendance for covid, as per discharge diagnosis
  covidemergency_date=patients.attended_emergency_care(
    returning="date_arrived",
    date_format="YYYY-MM-DD",
    on_or_after="anycovidvax_3_date",
    with_these_diagnoses = codelists.covid_emergency,
    find_first_match_in_period=True,
  ),
  
  # emergency attendance for covid, as per discharge diagnosis, resulting in discharge to hospital
  covidemergencyhosp_date=patients.attended_emergency_care(
    returning="date_arrived",
    date_format="YYYY-MM-DD",
    on_or_after="anycovidvax_3_date",
    find_first_match_in_period=True,
    with_these_diagnoses = codelists.covid_emergency,
    discharged_to = codelists.discharged_to_hospital,
  ),
  
  # emergency attendance for respiratory illness
  # FIXME -- need to define codelist
  # respemergency_date=patients.attended_emergency_care(
  #   returning="date_arrived",
  #   date_format="YYYY-MM-DD",
  #   on_or_after="anycovidvax_3_date",
  #   with_these_diagnoses = codelists.resp_emergency,
  #   find_first_match_in_period=True,
  # ),
  
  # emergency attendance for respiratory illness, resulting in discharge to hospital
  # FIXME -- need to define codelist
  # respemergencyhosp_date=patients.attended_emergency_care(
  #   returning="date_arrived",
  #   date_format="YYYY-MM-DD",
  #   on_or_after="anycovidvax_3_date",
  #   find_first_match_in_period=True,
  #   with_these_diagnoses = codelists.resp_emergency,
  #   discharged_to = codelists.discharged_to_hospital,
  # ),
  
  # any emergency attendance
  emergency_date=patients.attended_emergency_care(
    returning="date_arrived",
    on_or_after="anycovidvax_3_date",
    date_format="YYYY-MM-DD",
    find_first_match_in_period=True,
  ),
  
  # emergency attendance resulting in discharge to hospital
  emergencyhosp_date=patients.attended_emergency_care(
    returning="date_arrived",
    on_or_after="anycovidvax_3_date",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
    discharged_to = codelists.discharged_to_hospital,
  ),
  
  
  
  
  # unplanned hospital admission
  admitted_unplanned_date=patients.admitted_to_hospital(
    returning="date_admitted",
    on_or_after="anycovidvax_3_date",
    # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_patient_classification = ["1"], # ordinary admissions only
    date_format="YYYY-MM-DD",
    find_first_match_in_period=True,
  ),
  
  # planned hospital admission
  admitted_planned_date=patients.admitted_to_hospital(
    returning="date_admitted",
    on_or_after="anycovidvax_3_date",
    # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["11", "12", "13", "81"],
    with_patient_classification = ["1"], # ordinary admissions only 
    date_format="YYYY-MM-DD",
    find_first_match_in_period=True,
  ),
  
  # Positive covid admission prior to study start date
  covidadmitted_date=patients.admitted_to_hospital(
    returning="date_admitted",
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_these_diagnoses=codelists.covid_icd10,
    on_or_after="anycovidvax_3_date",
    date_format="YYYY-MM-DD",
    find_first_match_in_period=True,
  ),
  
  covidcritcare_date=patients.admitted_to_hospital(
      returning="date_admitted",
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_these_diagnoses=codelists.covid_icd10,
      with_at_least_one_day_in_critical_care=True,
      on_or_after="anycovidvax_3_date",
      date_format="YYYY-MM-DD",
      find_first_match_in_period=True,
    ),
  
  # Covid-related death
  coviddeath_date=patients.with_these_codes_on_death_certificate(
    codelists.covid_icd10,
    returning="date_of_death",
    date_format="YYYY-MM-DD",
  ),
  
  # All-cause death
  death_date=patients.died_from_any_cause(
    returning="date_of_death",
    date_format="YYYY-MM-DD",
  ),
  
  # a+e attendance due to fractures
  fractureemergency_date=patients.attended_emergency_care(
    returning="date_arrived",
    date_format="YYYY-MM-DD",
    on_or_after="anycovidvax_3_date",
    with_these_diagnoses = codelists.fractures_snomedECDS,
    find_first_match_in_period=True,
  ),
  
  # admission due to fractures
  fractureadmitted_date=patients.admitted_to_hospital(
    returning="date_admitted",
    on_or_after="anycovidvax_3_date",
    with_these_diagnoses=codelists.fractures_icd10,
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    date_format="YYYY-MM-DD",
    find_first_match_in_period=True,
  ),
  
  # death due to fractures
  fracturedeath_date=patients.with_these_codes_on_death_certificate(
    codelists.fractures_icd10,
    returning="date_of_death",
    date_format="YYYY-MM-DD",
  ),
  
  test_count = patients.with_test_result_in_sgss(
      pathogen = "SARS-CoV-2",
      test_result = "any",
      returning = "number_of_matches_in_period",
      between = ["anycovidvax_3_date", f"{postestfollowupend_date}"],
      restrict_to_earliest_specimen_date=False
    ),
  
  postest_count = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "number_of_matches_in_period",
    between = ["anycovidvax_3_date", f"{postestfollowupend_date}"],
    restrict_to_earliest_specimen_date=False
  ),


)
