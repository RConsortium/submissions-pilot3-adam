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

# ae <- haven::read_xpt(file.path("sdtm", "ae.xpt"))
# lb <- haven::read_xpt(file.path("sdtm", "lb.xpt"))
# vs <- haven::read_xpt(file.path("sdtm", "vs.xpt"))

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
    x== "BLACK OR AFRICAN AMERICAN" ~ 2,
    x== "AMERICAN INDIAN OR ALASKA NATIVE" ~ 6
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

# "BMIBL" VS
# "BMIBLGR1" VS
# "HEIGHTBL" VS
# "WEIGHTBL" VS
# "EDUCLVL" SC

# Disease information -----------------------------------------------------

# "DISONSDT"
# "DURDIS"
# "DURDSGR1"
# "MMSETOT" 

# More visit information --------------------------------------------------

# "VISIT1DT"
# "VISNUMEN"
# "RFENDT"  

# QC ----------------------------------------------------------------------

adsl <- adsl04

## Content check using in-house package
adsl[["AVGDD"]] <- as.numeric(adsl[["AVGDD"]])
adsl[["CUMDOSE"]] <- as.numeric(adsl[["CUMDOSE"]])

dfcompare(
  file = "adsl_compare"
  ,left = adsl_prod
  ,right = adsl
  ,keys = c("STUDYID", "USUBJID")
  ,showdiffs = 10000
  ,debug = F
)

  
  
## in define.xml following derivation should use TRTDURD: AVGDD = CUMDOSE/TRTDUR
## in define.xml we should mention that screen failures are removed from DM

















# Add Labels --------------------------------------------------------------



# QC dev vs prod ----------------------------------------------------------

## Metadata compare (labels)

# difflabels <- dplyr::setdiff(labsprod, labsupdated)
# discr_labels <- unlist(labsprod)[which(unlist(labsprod) %in% difflabels)]

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
