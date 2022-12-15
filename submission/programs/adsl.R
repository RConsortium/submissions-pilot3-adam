###########################################################################
#' developers : Steven Haesendonckx/Declan Hodges
#' date: 09DEC2022
#' modification History:
#'
###########################################################################

# Set up ------------------------------------------------------------------

library(haven)
library(admiral)
library(dplyr)
library(tidyr)
library(pilot3)

# read source -------------------------------------------------------------
# When SAS datasets are imported into R using read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values


dm <- convert_blanks_to_na(read_xpt(file.path("sdtm", "dm.xpt")))
ds <- convert_blanks_to_na(read_xpt(file.path("sdtm", "ds.xpt")))
ex <- convert_blanks_to_na(read_xpt(file.path("sdtm", "ex.xpt")))
qs <- convert_blanks_to_na(read_xpt(file.path("sdtm", "qs.xpt")))
sv <- convert_blanks_to_na(read_xpt(file.path("sdtm", "sv.xpt")))
vs <- convert_blanks_to_na(read_xpt(file.path("sdtm", "vs.xpt")))
sc <- convert_blanks_to_na(read_xpt(file.path("sdtm", "sc.xpt")))
mh <- convert_blanks_to_na(read_xpt(file.path("sdtm", "mh.xpt")))

# adsl_prod <- convert_blanks_to_na(read_xpt(file.path(adam[2], "adsl.xpt")))

# toprogram <- setdiff(colnames(adsl_prod), colnames(dm))

# Formats -----------------------------------------------------------------

# site groups - see CSR
# dm %>%
#   group_by(SITEID, ACTARMCD) %>%
#   summarise(n=n_distinct(USUBJID)) %>%
#   filter(n<3, ACTARMCD != "Scrnfail") %>%
#   summarise(n_distinct(SITEID))

# Disposition information -------------------------------------------------
# unique(ds[order(ds[["DSCAT"]]) , c("DSCAT", "DSDECOD")])

ds00 <- ds %>%
  filter(DSCAT == "DISPOSITION EVENT", DSDECOD != "SCREEN FAILURE") %>%
  derive_vars_dt(
    dtc = DSSTDTC,
    new_vars_prefix = "EOS",
    highest_imputation = "n",
  ) %>%
  mutate(
    DISCONFL = ifelse(!is.na(EOSDT) & DSDECOD != "COMPLETED", "Y", NA),
    DSRAEFL = ifelse(DSTERM == "ADVERSE EVENT", "Y", NA),
    DCDECOD = DSDECOD
  ) %>%
  select(STUDYID, USUBJID, EOSDT, DISCONFL, DSRAEFL, DSDECOD, DSTERM, DCDECOD)

# Treatment information ---------------------------------------------------

ex_dt <- ex %>%
  derive_vars_dt(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    highest_imputation = "n",
  ) %>%
  # treatment end is imputed by discontinuation if subject discontinued after visit 3 = randomization as per protocol
  derive_vars_merged(
    dataset_add = ds00,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(EOSDT = EOSDT),
    filter_add = DCDECOD != "COMPLETED"
  ) %>%
  derive_vars_dt(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    highest_imputation = "Y",
    min_dates = vars(EXSTDT),
    max_dates = vars(EOSDT),
    date_imputation = "last",
    flag_imputation = "none"
  ) %>%
  mutate(DOSE = EXDOSE * (EXENDT - EXSTDT + 1))

ex_dose <- ex_dt %>%
  group_by(STUDYID, USUBJID, EXTRT) %>%
  summarise(cnt = n_distinct(EXTRT), CUMDOSE = sum(DOSE))

ex_dose[which(ex_dose[["cnt"]] > 1), "USUBJID"] # are there subjects with mixed treatments?

adsl00 <- dm %>%
  select(-DOMAIN) %>%
  filter(ACTARMCD != "Scrnfail") %>%
  # planned treatment
  mutate(
    TRT01P = ARM,
    TRT01PN = case_when(
      ARM == "Placebo" ~ 0,
      ARM == "Xanomeline High Dose" ~ 81,
      ARM == "Xanomeline Low Dose" ~ 54
    )
  ) %>%
  # actual treatment - It is assumed TRT01A=TRT01P which is not really true.
  mutate(
    TRT01A = TRT01P,
    TRT01AN = TRT01PN
  ) %>%
  # treatment start
  derive_vars_merged(
    dataset_add = ex_dt,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        grepl("PLACEBO", EXTRT, fixed = TRUE))) &
      !is.na(EXSTDT),
    new_vars = vars(TRTSDT = EXSTDT),
    order = vars(EXSTDT, EXSEQ),
    mode = "first",
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  # treatment end
  derive_vars_merged(
    dataset_add = ex_dt,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        grepl("PLACEBO", EXTRT, fixed = TRUE))) &
      !is.na(EXENDT),
    new_vars = vars(TRTEDT = EXENDT),
    order = vars(EXENDT, EXSEQ),
    mode = "last",
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  # treatment duration
  derive_var_trtdurd() %>%
  # dosing
  left_join(ex_dose, by = c("STUDYID", "USUBJID")) %>%
  select(-cnt) %>%
  mutate(AVGDD = round(CUMDOSE / TRTDURD, digits = 1))

# Demographic grouping ----------------------------------------------------
# distinct(adsl_prod[which(adsl_prod$SITEGR1 == "900"), c("SITEID", "SITEGR1")])

adsl01 <- adsl00 %>%
  mutate(
    SITEGR1 = format_sitegr1(SITEID),
    AGEGR1 = format_agegr1(AGE),
    AGEGR1N = format_agegr1n(AGE),
    RACEN = format_racen(RACE)
  )

# Population flag ---------------------------------------------------------
# SAFFL - Y if ITTFL='Y' and TRTSDT ne missing. N otherwise
# ITTFL - Y if ARMCD ne ' '. N otherwise
# EFFFL - Y if SAFFL='Y AND at least one record in QS for ADAS-Cog and for CIBIC+ with VISITNUM>3, N otherwise
# these variables are also in suppdm, but define said derived

qstest <- distinct(qs[, c("QSTESTCD", "QSTEST")])

eff <- qs %>%
  filter(VISITNUM > 3, QSTESTCD %in% c("CIBIC", "ACTOT")) %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(effcnt = n_distinct(QSTESTCD))

adsl02 <- adsl01 %>%
  left_join(eff, by = c("STUDYID", "USUBJID")) %>%
  mutate(
    SAFFL = case_when(
      ARMCD != "Scrnfail" & ARMCD != "" & !is.na(TRTSDT) ~ "Y",
      ARMCD == "Scrnfail" ~ NA_character_,
      TRUE ~ "N"
    ),
    ITTFL = case_when(
      ARMCD != "Scrnfail" & ARMCD != "" ~ "Y",
      ARMCD == "Scrnfail" ~ NA_character_,
      TRUE ~ "N"
    ),
    EFFFL = case_when(
      ARMCD != "Scrnfail" & ARMCD != "" & !is.na(TRTSDT) & effcnt == 2 ~ "Y",
      ARMCD == "Scrnfail" ~ NA_character_,
      TRUE ~ "N"
    )
  )

# Study Visit compliance --------------------------------------------------
# these variables are also in suppdm, but define said derived

sv00 <- sv %>%
  select(STUDYID, USUBJID, VISIT, VISITDY, SVSTDTC) %>%
  mutate(
    FLG = "Y",
    VISITCMP = case_when(
      VISIT == "WEEK 8" ~ "COMP8FL",
      VISIT == "WEEK 16" ~ "COMP16FL",
      VISIT == "WEEK 24" ~ "COMP24FL",
      TRUE ~ "ZZZ" # ensures every subject with one visit will get a row with minimally 'N'
    )
  ) %>%
  arrange(STUDYID, USUBJID, VISITDY) %>%
  distinct(STUDYID, USUBJID, VISITCMP, FLG) %>%
  pivot_wider(names_from = VISITCMP, values_from = FLG, values_fill = "N") %>%
  select(-ZZZ)

adsl03 <- adsl02 %>%
  left_join(sv00, by = c("STUDYID", "USUBJID"))

# Disposition -------------------------------------------------------------

adsl04 <- adsl03 %>%
  left_join(ds00, by = c("STUDYID", "USUBJID")) %>%
  select(-DSDECOD) %>%
  derive_var_disposition_status(
    dataset_ds = ds00,
    new_var = EOSSTT,
    status_var = DSDECOD, # this variable is removed after reformat
    filter_ds = !is.na(USUBJID)
  ) %>%
  derive_vars_disposition_reason(
    dataset_ds = ds00,
    new_var = DCSREAS,
    reason_var = DSDECOD,
    filter_ds = !is.na(USUBJID),
    format_new_vars = format_dcsreas # could not include dsterm in formatting logic
  ) %>%
  mutate(DCSREAS = ifelse(DSTERM == "PROTOCOL ENTRY CRITERIA NOT MET", "I/E Not Met", DCSREAS))

# Baseline variables ------------------------------------------------------
# selection definition from define

vs00 <- vs %>%
  filter((VSTESTCD == "HEIGHT" & VISITNUM == 1) | (VSTESTCD == "WEIGHT" & VISITNUM == 3)) %>%
  mutate(AVAL = round(VSSTRESN, digits = 1)) %>%
  select(STUDYID, USUBJID, VSTESTCD, AVAL) %>%
  pivot_wider(names_from = VSTESTCD, values_from = AVAL, names_glue = "{VSTESTCD}BL") %>%
  mutate(
    BMIBL = round(WEIGHTBL / (HEIGHTBL / 100)^2, digits = 1),
    BMIBLGR1 = format_bmiblgr1(BMIBL)
  )

sc00 <- sc %>%
  filter(SCTESTCD == "EDLEVEL") %>%
  select(STUDYID, USUBJID, SCTESTCD, SCSTRESN) %>%
  pivot_wider(names_from = SCTESTCD, values_from = SCSTRESN, names_glue = "EDUCLVL")

adsl05 <- adsl04 %>%
  left_join(vs00, by = c("STUDYID", "USUBJID")) %>%
  left_join(sc00, by = c("STUDYID", "USUBJID"))

# Disease information -----------------------------------------------------

visit1dt <- sv %>%
  filter(VISITNUM == 1) %>%
  derive_vars_dt(
    dtc = SVSTDTC,
    new_vars_prefix = "VISIT1",
  ) %>%
  select(STUDYID, USUBJID, VISIT1DT)

visnumen <- sv %>%
  filter(VISITNUM < 100) %>%
  arrange(STUDYID, USUBJID, SVSTDTC) %>%
  group_by(STUDYID, USUBJID) %>%
  slice(n()) %>%
  ungroup() %>%
  mutate(VISNUMEN = ifelse(round(VISITNUM, digits = 0) == 13, 12, round(VISITNUM, digits = 0))) %>%
  select(STUDYID, USUBJID, VISNUMEN)

disonsdt <- mh %>%
  filter(MHCAT == "PRIMARY DIAGNOSIS") %>%
  derive_vars_dt(
    dtc = MHSTDTC,
    new_vars_prefix = "DISONS",
  ) %>%
  select(STUDYID, USUBJID, DISONSDT)

adsl06 <- adsl05 %>%
  left_join(visit1dt, by = c("STUDYID", "USUBJID")) %>%
  left_join(visnumen, by = c("STUDYID", "USUBJID")) %>%
  left_join(disonsdt, by = c("STUDYID", "USUBJID")) %>%
  derive_vars_duration(
    new_var = DURDIS,
    start_date = DISONSDT,
    end_date = VISIT1DT,
    out_unit = "months",
    add_one = TRUE
  ) %>%
  mutate(
    DURDIS = round(DURDIS, digits = 1),
    DURDSGR1 = format_durdsgr1(DURDIS)
  ) %>%
  derive_vars_dt(
    dtc = RFENDTC,
    new_vars_prefix = "RFEN",
  )

mmsetot <- qs %>%
  filter(QSCAT == "MINI-MENTAL STATE") %>%
  group_by(STUDYID, USUBJID) %>%
  summarise(MMSETOT = sum(as.numeric(QSORRES), na.rm = TRUE)) %>%
  select(STUDYID, USUBJID, MMSETOT)

adsl07 <- adsl06 %>%
  left_join(mmsetot, by = c("STUDYID", "USUBJID"))

# Add Labels --------------------------------------------------------------

adsl <- adsl07[, c(
  "STUDYID", "USUBJID", "SUBJID", "SITEID", "SITEGR1", "ARM",
  "TRT01P", "TRT01PN", "TRT01A", "TRT01AN", "TRTSDT", "TRTEDT",
  "TRTDURD", "AVGDD", "CUMDOSE", "AGE", "AGEGR1", "AGEGR1N", "AGEU",
  "RACE", "RACEN", "SEX", "ETHNIC", "SAFFL", "ITTFL", "EFFFL",
  "COMP8FL", "COMP16FL", "COMP24FL", "DISCONFL", "DSRAEFL", "DTHFL",
  "BMIBL", "BMIBLGR1", "HEIGHTBL", "WEIGHTBL", "EDUCLVL", "DISONSDT",
  "DURDIS", "DURDSGR1", "VISIT1DT", "RFSTDTC", "RFENDTC", "VISNUMEN",
  "RFENDT", "DCDECOD", "EOSSTT", "DCSREAS", "MMSETOT"
)]

# labs_prod <- sapply(colnames(adsl_prod), FUN = function(x) attr(adsl_prod[[x]], "label"))
# labs <- sapply(colnames(adsl), FUN = function(x) attr(adsl[[x]], "label"))

# setdiff(labs_prod, labs)

# labs[unlist(lapply(labs, is.null))]


adsl[["AVGDD"]] <- as.numeric(adsl[["AVGDD"]])
adsl[["CUMDOSE"]] <- as.numeric(adsl[["CUMDOSE"]])

attr(adsl[["TRT01P"]], "label") <- "Planned Treatment for Period 01"
attr(adsl[["TRT01PN"]], "label") <- "Planned Treatment for Period 01 (N)"
attr(adsl[["TRT01A"]], "label") <- "Actual Treatment for Period 01"
attr(adsl[["TRT01AN"]], "label") <- "Actual Treatment for Period 01 (N)"
attr(adsl[["DTHFL"]], "label") <- "Subject Died?"
attr(adsl[["SITEGR1"]], "label") <- "Pooled Site Group 1"
attr(adsl[["TRTSDT"]], "label") <- "Date of First Exposure to Treatment"
attr(adsl[["TRTEDT"]], "label") <- "Date of Last Exposure to Treatment"
attr(adsl[["TRTDURD"]], "label") <- "Total Treatment Duration (Days)"
attr(adsl[["AVGDD"]], "label") <- "Avg Daily Dose (as planned)"
attr(adsl[["CUMDOSE"]], "label") <- "Cumulative Dose (as planned)"
attr(adsl[["AGEGR1"]], "label") <- "Pooled Age Group 1"
attr(adsl[["AGEGR1N"]], "label") <- "Pooled Age Group 1 (N)"
attr(adsl[["RACEN"]], "label") <- "Race (N)"
attr(adsl[["SAFFL"]], "label") <- "Safety Population Flag"
attr(adsl[["ITTFL"]], "label") <- "Intent-To-Treat Population Flag"
attr(adsl[["EFFFL"]], "label") <- "Efficacy Population Flag"
attr(adsl[["COMP8FL"]], "label") <- "Completers of Week 8 Population Flag"
attr(adsl[["COMP16FL"]], "label") <- "Completers of Week 16 Population Flag"
attr(adsl[["COMP24FL"]], "label") <- "Completers of Week 24 Population Flag"
attr(adsl[["DISCONFL"]], "label") <- "Did the Subject Discontinue the Study?"
attr(adsl[["DSRAEFL"]], "label") <- "Discontinued due to AE?"
attr(adsl[["BMIBL"]], "label") <- "Baseline BMI (kg/m^2)"
attr(adsl[["BMIBLGR1"]], "label") <- "Pooled Baseline BMI Group 1"
attr(adsl[["HEIGHTBL"]], "label") <- "Baseline Height (cm)"
attr(adsl[["WEIGHTBL"]], "label") <- "Baseline Weight (kg)"
attr(adsl[["EDUCLVL"]], "label") <- "Years of Education"
attr(adsl[["DISONSDT"]], "label") <- "Date of Onset of Disease"
attr(adsl[["DURDIS"]], "label") <- "Duration of Disease (Months)"
attr(adsl[["DURDSGR1"]], "label") <- "Pooled Disease Duration Group 1"
attr(adsl[["VISIT1DT"]], "label") <- "Date of Visit 1"
attr(adsl[["RFENDT"]], "label") <- "Date of Discontinuation/Completion"
attr(adsl[["VISNUMEN"]], "label") <- "End of Trt Visit (Vis 12 or Early Term.)"
attr(adsl[["EOSSTT"]], "label") <- "End of Study Status"
attr(adsl[["DCSREAS"]], "label") <- "Reason for Discontinuation from Study"
attr(adsl[["MMSETOT"]], "label") <- "MMSE Total"

labsupdated <- sapply(colnames(adsl), FUN = function(x) attr(adsl[[x]], "label"))
labsupdated[unlist(lapply(labsupdated, is.null))]

# Output ------------------------------------------------------------------

write_xpt(adsl, file.path("submission/datasets/adsl.xpt"))

# END of Code -------------------------------------------------------------
