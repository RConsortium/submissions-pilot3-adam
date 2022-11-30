###########################################################################
#' developers : Steven Haesendonckx/
#' date: 28NOV2022
#' modification History:
#' 
###########################################################################

# Set up ------------------------------------------------------------------

fcts <- c("eff_models.R", "fmt.R", "helpers.R", "Tplyr_helpers.R")
invisible(sapply(fcts, FUN = function(x) source(file.path("R/", x), )))

library(haven)
library(admiral)
library(dplyr)

# read source -------------------------------------------------------------
# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values


dm <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "dm.xpt")))
ds <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "ds.xpt")))
ex <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "ex.xpt")))
qs <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "qs.xpt")))
sv <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "sv.xpt")))
vs <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "vs.xpt")))
sc <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "sc.xpt")))
mh <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "mh.xpt")))

adsl_prod <- admiral::convert_blanks_to_na(haven::read_xpt(file.path(adam[2], "adsl.xpt")))

toprogram <- setdiff(colnames(adsl_prod), colnames(dm))

# Formats -----------------------------------------------------------------

# site groups - if not pooled then SITEGR1=SITEID. If pooled, SITEGR1 will be 900 - no SAP available
format_siteid <- function(x) {
  dplyr::case_when(
    x %in% c("702", "706", "707", "711", "714", "715", "717") ~ "900",
    TRUE ~ x
  )
}

# AGEGR1 - AGEGR1 = 1 if AGE <65. AGEGR1 = 2 if AGE 65-80. AGEGR1 = 3 if AGE >80
format_agegr1 <- function(x) {
  dplyr::case_when(
    x < 65 ~ "<65",
    dplyr::between(x, 65, 80) ~ "65-80",
    x > 80 ~ ">80",
  )
}

format_agegr1n <- function(x) {
  dplyr::case_when(
    x < 65 ~ 1,
    dplyr::between(x, 65, 80) ~ 2,
    x > 80 ~ 3,
  )
}

# Race group numbering
format_racen <- function(x) {
  dplyr::case_when(
    x == "WHITE" ~ 1,
    x == "BLACK OR AFRICAN AMERICAN" ~ 2,
    x == "AMERICAN INDIAN OR ALASKA NATIVE" ~ 6
  )
}
# BMI group
format_bmi <- function(x) {
  dplyr::case_when(
    !is.na(x) & x < 25 ~ "<25",
    25 <= x & x < 30 ~ "25-<30",
    30 <= x ~ ">=30"
  )
}

# Disease duration group
format_dis <- function(x) {
  dplyr::case_when(
    !is.na(x) & x < 12 ~ "<12",
    12 <= x ~ ">=12"
  )
}



# Disposition information -------------------------------------------------
unique(ds[order(ds[["DSCAT"]]) , c("DSCAT", "DSDECOD")])

ds00 <- ds %>%
  dplyr::filter(DSCAT == "DISPOSITION EVENT", DSDECOD != "SCREEN FAILURE") %>%
  admiral::derive_vars_dt(
    dtc = DSSTDTC,
    new_vars_prefix = "EOS",
    highest_imputation = "n",
  ) %>%
  dplyr::mutate(DISCONFL = ifelse(!is.na(EOSDT) & DSDECOD != "COMPLETED", "Y", NA),
                DSRAEFL = ifelse(DSTERM == "ADVERSE EVENT", "Y", NA),
                DCDECOD = DSDECOD
  ) %>%
  dplyr::select(STUDYID, USUBJID, EOSDT, DISCONFL, DSRAEFL, DSDECOD, DSTERM, DCDECOD) 

# Treatment information ---------------------------------------------------

ex_dt <- ex %>%
  admiral::derive_vars_dt(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST",
    highest_imputation = "n",
  ) %>%
  # treatment end is imputed by discontinuation if subject discontinued after visit 3 = randomization as per protocol
  admiral::derive_vars_merged(
    dataset_add = ds00,
    by_vars = vars(STUDYID, USUBJID),
    new_vars = vars(EOSDT = EOSDT),
    filter_add = DCDECOD != "COMPLETED"
  ) %>%
  admiral::derive_vars_dt(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    highest_imputation = "Y",
    min_dates = vars(EXSTDT), 
    max_dates = vars(EOSDT),
    date_imputation = "last"
  ) %>%
  dplyr::mutate(DOSE = EXDOSE* (EXENDT-EXSTDT + 1))

ex_dose <- ex_dt %>%
  dplyr::group_by(STUDYID, USUBJID, EXTRT) %>%
  dplyr::summarise(cnt = dplyr::n_distinct(EXTRT), CUMDOSE = sum(DOSE)) 

ex_dose[which(ex_dose[["cnt"]] > 1), "USUBJID"] # are there subjects with mixed treatments?

adsl00 <- dm %>%
  dplyr::select(-DOMAIN) %>%
  dplyr::filter(ACTARMCD != "Scrnfail") %>%
  
  # planned treatment
  dplyr::mutate(TRT01P = ARM, 
                TRT01PN = dplyr::case_when(ARM == "Placebo" ~ 0,
                                           ARM == "Xanomeline High Dose" ~ 81,
                                           ARM =="Xanomeline Low Dose" ~ 54)) %>%
  
  # actual treatment - It is assumed TRT01A=TRT01P which is not really true.
  dplyr::mutate(TRT01A = TRT01P, 
                TRT01AN = TRT01PN) %>%
  
  # treatment start
  admiral::derive_vars_merged(
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
  admiral::derive_vars_merged(
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
  admiral::derive_var_trtdurd() %>%
  
  # dosing
  dplyr::left_join(ex_dose, by = c("STUDYID", "USUBJID")) %>%
  dplyr::select(-cnt) %>%
  dplyr::mutate(AVGDD = round(CUMDOSE/TRTDURD, digits = 1))

# Demographic grouping ----------------------------------------------------
#distinct(adsl_prod[which(adsl_prod$SITEGR1 == "900"), c("SITEID", "SITEGR1")])

adsl01 <- adsl00 %>%
  dplyr::mutate(
    SITEGR1 = format_siteid(SITEID),
    AGEGR1 = format_agegr1(AGE),
    AGEGR1N = format_agegr1n(AGE),
    RACEN = format_racen(RACE)
  )

# Population flag ---------------------------------------------------------
# SAFFL - Y if ITTFL='Y' and TRTSDT ne missing. N otherwise
# ITTFL - Y if ARMCD ne ' '. N otherwise
# EFFFL - Y if SAFFL='Y AND at least one record in QS for ADAS-Cog and for CIBIC+ with VISITNUM>3, N otherwise
# these variables are also in suppdm, but define said derived

qstest <- distinct(qs[,c("QSTESTCD", "QSTEST")])

eff <- qs %>%
  dplyr::filter(VISITNUM>3, QSTESTCD %in% c("CIBIC", "ACTOT")) %>%
  dplyr::group_by(STUDYID, USUBJID) %>%
  dplyr::summarise(effcnt = dplyr::n_distinct(QSTESTCD))

adsl02 <- adsl01 %>%
  dplyr::left_join(eff, by = c("STUDYID", "USUBJID")) %>%
  dplyr::mutate(SAFFL = dplyr::case_when(
                          ARMCD != "Scrnfail" & ARMCD != "" & !is.na(TRTSDT) ~ "Y",
                          ARMCD == "Scrnfail" ~ NA_character_,
                          TRUE ~ "N"
                        ),
                
                ITTFL = dplyr::case_when(
                  ARMCD != "Scrnfail" & ARMCD != "" ~ "Y",
                  ARMCD == "Scrnfail" ~ NA_character_,
                  TRUE ~ "N"
                        ),
                
                EFFFL = dplyr::case_when(
                  ARMCD != "Scrnfail" & ARMCD != "" & !is.na(TRTSDT) & effcnt == 2 ~ "Y",
                  ARMCD == "Scrnfail" ~ NA_character_,
                  TRUE ~ "N"
                )
  )

# Study Visit compliance --------------------------------------------------
# these variables are also in suppdm, but define said derived
sv00 <- sv %>%
  dplyr::select(STUDYID, USUBJID, VISIT, VISITDY, SVSTDTC) %>%
  dplyr::mutate(FLG = "Y",
                VISITCMP = dplyr::case_when(
                          VISIT == "WEEK 8" ~ "COMP8FL",
                          VISIT == "WEEK 16" ~ "COMP16FL",
                          VISIT == "WEEK 24" ~ "COMP24FL",
                          TRUE ~ "ZZZ" # ensures every subject with one visit will get a row with minimally 'N'
                )) %>%
  dplyr::arrange(STUDYID, USUBJID, VISITDY) %>%
  dplyr::distinct(STUDYID, USUBJID, VISITCMP, FLG) %>%
  tidyr::pivot_wider(names_from = VISITCMP, values_from = FLG, values_fill = "N") %>%
  dplyr::select(-ZZZ)

adsl03 <- adsl02 %>%
  dplyr::left_join(sv00, by = c("STUDYID", "USUBJID")) 

# Disposition -------------------------------------------------------------

format_dcsreas <- function(dsdecod) {
  dplyr::case_when(
    dsdecod == "ADVERSE EVENT" ~ "Adverse Event",
    dsdecod == "STUDY TERMINATED BY SPONSOR" ~ "Sponsor Decision",
    dsdecod == "DEATH" ~ "Death",
    dsdecod == "WITHDRAWAL BY SUBJECT" ~ "Withdrew Consent",
    dsdecod == "PHYSICIAN DECISION" ~ "Physician Decision",
    dsdecod == "PROTOCOL VIOLATION" ~ "Protocol Violation",
    dsdecod == "LOST TO FOLLOW-UP" ~ "Lost to Follow-up",
    dsdecod == "LACK OF EFFICACY" ~ "Lack of Efficacy"
  )
}

adsl04 <- adsl03 %>%
  dplyr::left_join(ds00, by = c("STUDYID", "USUBJID")) %>%
  dplyr::select(-DSDECOD) %>%
  admiral::derive_var_disposition_status(
    dataset_ds = ds00,
    new_var = EOSSTT,
    status_var = DSDECOD, #this variable is removed after reformat
    filter_ds = !is.na(USUBJID)
  ) %>%
  admiral::derive_vars_disposition_reason(
    dataset_ds = ds00,
    new_var = DCSREAS,
    reason_var = DSDECOD,
    filter_ds = !is.na(USUBJID),
    format_new_vars = format_dcsreas #could not include dsterm in formatting logic
  ) %>%
  dplyr::mutate(DCSREAS = ifelse(DSTERM == "PROTOCOL ENTRY CRITERIA NOT MET", "I/E Not Met", DCSREAS))

# Baseline variables ------------------------------------------------------
# selection definition from define

vs00 <- vs %>%
  dplyr::filter((VSTESTCD == "HEIGHT" & VISITNUM == 1) | (VSTESTCD == "WEIGHT" & VISITNUM == 3)) %>%
  dplyr::mutate(AVAL = round(VSSTRESN, digits = 1)) %>%
  dplyr::select(STUDYID, USUBJID, VSTESTCD, AVAL) %>%
  tidyr::pivot_wider(names_from = VSTESTCD, values_from = AVAL, names_glue = "{VSTESTCD}BL") %>%
  dplyr::mutate(BMIBL = round(WEIGHTBL/(HEIGHTBL/100)^2, digits = 1),
                BMIBLGR1 = format_bmi(BMIBL)
               ) 

sc00 <- sc %>%
  dplyr::filter(SCTESTCD == "EDLEVEL") %>%
  dplyr::select(STUDYID, USUBJID, SCTESTCD, SCSTRESN) %>%
  tidyr::pivot_wider(names_from = SCTESTCD, values_from = SCSTRESN, names_glue = "EDUCLVL")

adsl05 <- adsl04 %>%
  dplyr::left_join(vs00, by = c("STUDYID", "USUBJID")) %>%
  dplyr::left_join(sc00, by = c("STUDYID", "USUBJID"))

# Disease information -----------------------------------------------------

visit1dt <- sv %>%
  dplyr::filter(VISITNUM == 1) %>%
  admiral::derive_vars_dt(
    dtc = SVSTDTC,
    new_vars_prefix = "VISIT1",
  ) %>%
  dplyr::select(STUDYID, USUBJID, VISIT1DT)

visnumen <- sv %>%
  dplyr::filter(VISITNUM < 100) %>%
  dplyr::arrange(STUDYID, USUBJID, SVSTDTC) %>%
  dplyr::group_by(STUDYID, USUBJID) %>%
  dplyr::slice(n()) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(VISNUMEN = ifelse(round(VISITNUM, digits = 0) == 13, 12, round(VISITNUM, digits = 0))) %>%
  dplyr::select(STUDYID, USUBJID, VISNUMEN)

disonsdt <- mh %>%
  dplyr::filter(MHCAT == 'PRIMARY DIAGNOSIS') %>%
  admiral::derive_vars_dt(
    dtc = MHSTDTC,
    new_vars_prefix = "DISONS",
  ) %>%
  dplyr::select(STUDYID, USUBJID, DISONSDT)

adsl06 <- adsl05 %>%
  dplyr::left_join(visit1dt, by = c("STUDYID", "USUBJID")) %>%
  dplyr::left_join(visnumen, by = c("STUDYID", "USUBJID")) %>%
  dplyr::left_join(disonsdt, by = c("STUDYID", "USUBJID")) %>%
  admiral::derive_vars_duration(new_var=DURDIS,
                              start_date = DISONSDT,
                              end_date = VISIT1DT, 
                              out_unit = "months",
                              add_one = TRUE) %>%
  dplyr::mutate(DURDIS = round(DURDIS, digits = 1),
                DURDSGR1 = format_dis(DURDIS)) %>%
  admiral::derive_vars_dt(
    dtc = RFENDTC,
    new_vars_prefix = "RFEN",
  ) 

mmsetot <- qs %>%
  dplyr::filter(QSCAT == "MINI-MENTAL STATE") %>%
  dplyr::group_by(STUDYID, USUBJID) %>%
  dplyr::summarise(MMSETOT = sum(as.numeric(QSORRES), na.rm = TRUE)) %>%
  dplyr::select(STUDYID, USUBJID, MMSETOT)

adsl07 <- adsl06 %>%
  dplyr::left_join(mmsetot, by = c("STUDYID", "USUBJID"))
  
# Add Labels --------------------------------------------------------------

#dput(colnames(adsl_prod))

adsl <- adsl07[ , c("STUDYID", "USUBJID", "SUBJID", "SITEID", "SITEGR1", "ARM", 
"TRT01P", "TRT01PN", "TRT01A", "TRT01AN", "TRTSDT", "TRTEDT", 
"TRTDURD", "AVGDD", "CUMDOSE", "AGE", "AGEGR1", "AGEGR1N", "AGEU", 
"RACE", "RACEN", "SEX", "ETHNIC", "SAFFL", "ITTFL", "EFFFL", 
"COMP8FL", "COMP16FL", "COMP24FL", "DISCONFL", "DSRAEFL", "DTHFL", 
"BMIBL", "BMIBLGR1", "HEIGHTBL", "WEIGHTBL", "EDUCLVL", "DISONSDT", 
"DURDIS", "DURDSGR1", "VISIT1DT", "RFSTDTC", "RFENDTC", "VISNUMEN", 
"RFENDT", "DCDECOD", "EOSSTT", "DCSREAS", "MMSETOT")]

labs_prod <- sapply(colnames(adsl_prod), FUN = function(x) attr(adsl_prod[[x]], "label"))
labs <- sapply(colnames(adsl), FUN = function(x) attr(adsl[[x]], "label"))

setdiff(labs_prod, labs)

labs[unlist(lapply(labs,is.null))]


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
labsupdated[unlist(lapply(labsupdated,is.null))]

# QC dev vs prod ----------------------------------------------------------

## Metadata compare (labels)

# difflabels <- dplyr::setdiff(labs_prod, labsupdated)
# discr_labels <- unlist(labs_prod)[which(unlist(labs_prod) %in% difflabels)]

## Content check using in-house package

 # dfcompare(
 #      file = "tlcompare"
 #     ,left = prod
 #     ,right = adtte
 #     ,keys = c("STUDYID", "USUBJID")
 #     ,showdiffs = 1000
 #     ,debug = FALSE
 # )

# Output ------------------------------------------------------------------

haven::write_xpt(adsl, file.path("submission/datasets/adsl.xpt"))

# END of Code -------------------------------------------------------------

# QC ----------------------------------------------------------------------


dfcompare(
  file = "adsl_compare"
  ,left = adsl_prod
  ,right = adsl
  ,keys = c("STUDYID", "USUBJID")
  ,showdiffs = 10000
  ,debug = F
)
