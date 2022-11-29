## ----setup, message=FALSE-------------------------------------------------------------------------------------------
# CRAN package, please using install.packages() to install
library(dplyr)
library(admiral)
library(metacore)
library(metatools)
library(stringr)


## -------------------------------------------------------------------------------------------------------------------
qs <- haven::read_xpt(file.path(path$sdtm, "qs.xpt"))
adsl <- haven::read_xpt(file.path(path$adam, "adsl.xpt"))
adas_qc <- haven::read_xpt(file.path(path$adam, "adadas.xpt"))


## -------------------------------------------------------------------------------------------------------------------
qs <- convert_blanks_to_na(qs)


## -------------------------------------------------------------------------------------------------------------------
## placeholder for origin=predecessor, use metatool::build_from_derived()

## derive ADT/ADY
# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P)
adas1 <- qs %>%
  # subset to interested QSTESTCD
  filter(QSTESTCD %in%
    c(str_c("ACITM", str_pad(1:14, 2, pad = "0")), "ACTOT")) %>%
  # Join ADSL with QS (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = QSDTC
  ) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

## derive AVISIT/AVISITN
adas2 <- adas1 %>%
  mutate(
    AVISIT = case_when(
      ADY <= 1 ~ "Baseline",
      ADY >= 2 & ADY <= 84 ~ "Week 8",
      ADY >= 85 & ADY <= 140 ~ "Week 16",
      ADY > 140 ~ "Week 24",
      TRUE ~ NA_character_
    ),
    PARAMCD = QSTESTCD
  )


## -------------------------------------------------------------------------------------------------------------------
## placeholder: replace derive_vars_merged_lookup by metatool::create_var_from_codelist()/create_cat_var()
## Add PARAMCD PARAM and PARAMN - from LOOK-UP table ----
# Replace derive_vars_merged_lookup() with PARAMCD lookup function


## -------------------------------------------------------------------------------------------------------------------
## derive AWRANGE/AWTARGET/AWTDIFF/AWLO/AWHI/AWU
aw_lookup <- tribble(
  ~AVISIT, ~AWRANGE, ~AWTARGET, ~AWLO, ~AWHI,
  "Baseline", "<=1", 1, NA_integer_, 1,
  "Week 8", "2-84", 56, 2, 84,
  "Week 16", "85-140", 112, 85, 140,
  "Week 24", ">140", 168, 141, NA_integer_
)

adas3 <- derive_vars_merged(
  adas2,
  dataset_add = aw_lookup,
  by_vars = vars(AVISIT)
) %>%
  mutate(
    AWTDIFF = abs(AWTARGET - ADY),
    AWU = "DAYS"
  )


## baseline information ---- ABLFL/AVAL/BASE/CHG/PCHG
adas4 <- adas3 %>%
  mutate(
    ABLFL = QSBLFL,
    AVAL = QSSTRESN
  ) %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()



## -------------------------------------------------------------------------------------------------------------------
adas5 <- adas4 %>%
  mutate(diff = AWTARGET - ADY) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, PARAMCD, AVISIT),
      order = vars(AWTDIFF, diff),
      new_var = ANL01FL,
      mode = "first"
    ),
    filter = !is.na(AVISIT)
  )


## -------------------------------------------------------------------------------------------------------------------
## placeholder for using metacore/metatools
