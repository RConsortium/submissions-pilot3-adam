###########################################################################
#' developers : Kangjie Zhang/
#' date: 29NOV2022
#' modification History:
#' program ADADAS
###########################################################################

## setup
library(dplyr)
library(admiral)
library(metacore)
library(metatools)
library(stringr)
library(xportr)

qs <- haven::read_xpt(file.path("sdtm", "qs.xpt"))
adsl <- haven::read_xpt(file.path("adam", "adsl.xpt"))

qs <- convert_blanks_to_na(qs)


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

## derive AVISIT/AVISITN/ABLFL/AVAL
avisitn_lookup <- tibble::tribble(
  ~AVISIT, ~AVISITN,
  "Baseline", 0,
  "Week 8", 8,
  "Week 16", 16,
  "Week 24", 24
)

adas2 <- adas1 %>%
  mutate(
    AVISIT = case_when(
      ADY <= 1 ~ "Baseline",
      ADY >= 2 & ADY <= 84 ~ "Week 8",
      ADY >= 85 & ADY <= 140 ~ "Week 16",
      ADY > 140 ~ "Week 24",
      TRUE ~ NA_character_
    ),
    PARAMCD = QSTESTCD,
    ABLFL = QSBLFL,
    AVAL = QSSTRESN
  ) %>%
  # Add AVISITN
  derive_vars_merged(
    dataset_add = avisitn_lookup,
    by_vars = vars(AVISIT)
  )

## placeholder: replace look up table by
## metatool::create_var_from_codelist()/create_cat_var()

# derive PARAMCD=ACTOT, DTYPE=LOCF
# A dataset with combinations of PARAMCD, AVISIT which are expected.
actot_expected_obsv <- tibble::tribble(
  ~PARAMCD, ~AVISITN, ~AVISIT,
  "ACTOT", 0, "Baseline",
  "ACTOT", 8, "Week 8",
  "ACTOT", 16, "Week 16",
  "ACTOT", 24, "Week 24"
)

adas_locf <-derive_locf_records(
  data = adas2,
  dataset_expected_obs = actot_expected_obsv,
  by_vars = vars(STUDYID, USUBJID, PARAMCD),
  order = vars(AVISITN, AVISIT)
)

## derive AWRANGE/AWTARGET/AWTDIFF/AWLO/AWHI/AWU
aw_lookup <- tribble(
  ~AVISIT, ~AWRANGE, ~AWTARGET, ~AWLO, ~AWHI,
  "Baseline", "<=1", 1, NA_integer_, 1,
  "Week 8", "2-84", 56, 2, 84,
  "Week 16", "85-140", 112, 85, 140,
  "Week 24", ">140", 168, 141, NA_integer_
)

adas3 <- derive_vars_merged(
  adas_locf,
  dataset_add = aw_lookup,
  by_vars = vars(AVISIT)
) %>%
  mutate(
    AWTDIFF = abs(AWTARGET - ADY),
    AWU = "DAYS"
  )


## baseline information BASE/CHG/PCHG
adas4 <- adas3 %>%
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


## ANL01FL
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


## remaining work:
## 1. derive PARAMCD=ACTOT, DTYPE=LOCF
## 2. pending define to use metacore/metatools + output XPT
## 3. QC in qc_adadas.R program


## placeholder for derive PARAMCD=ACTOT, DTYPE=LOCF
## placeholder for using metacore/metatools
## out to an XPT
adas5 %>%
  #  pending define: xportr_type(adsl_spec, "ADADAS") %>%
  #  pending define: xportr_length(adsl_spec, "ADADAS") %>%
  xportr_write("submission/datasets/adadas.xpt",
    label = "ADAS-COG Analysis Dataset"
  )
