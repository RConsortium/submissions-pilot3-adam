#' Derive SITEGR1
#' 
#' `format_sitegr1` derives the SITEGR1 variable, which pools sites with fewer than 3 subjects in any one treatment group.
format_sitegr1 <- function(x) {
  case_when(
    x %in% c("702", "706", "707", "711", "714", "715", "717") ~ "900",
    TRUE ~ x
  )
}

#' Derive AGEGR1
#' 
#' `format_agegr1` derives the AGEGR1 variable, which pools subjects into age groups.
format_agegr1 <- function(x) {
  case_when(
    x < 65 ~ "<65",
    between(x, 65, 80) ~ "65-80",
    x > 80 ~ ">80",
  )
}

#' Derive AGEGR1N
#' 
#' `format_agegr1n` derives the AGEGR1N variable, which is the numeric code for AGEGR1.
format_agegr1n <- function(x) {
  case_when(
    x < 65 ~ 1,
    between(x, 65, 80) ~ 2,
    x > 80 ~ 3,
  )
}

#' Derive RACEN
#' 
#' `format_racen` derives the RACEN variable, which is the numeric code for RACE.
format_racen <- function(x) {
  case_when(
    x == "WHITE" ~ 1,
    x == "BLACK OR AFRICAN AMERICAN" ~ 2,
    x == "AMERICAN INDIAN OR ALASKA NATIVE" ~ 6
  )
}

#' Derive BMIBLGR1
#' 
#' `format_bmiblgr1` derives the BMIBLGR1 variable, which pools subjects into baseline BMI groups.
format_bmiblgr1 <- function(x) {
  case_when(
    !is.na(x) & x < 25 ~ "<25",
    25 <= x & x < 30 ~ "25-<30",
    30 <= x ~ ">=30"
  )
}

#' Derive DURDSGR1
#' 
#' `format_durdsgr1` derives the DURDSGR1 variable, which pools subjects into disease duration groups.
format_durdsgr1 <- function(x) {
  case_when(
    !is.na(x) & x < 12 ~ "<12",
    12 <= x ~ ">=12"
  )
}

#' Derive DCSREAS
#' 
#' `format_dcsreas` derives the DCSREAS variable, which contains subjects' disposition reason(s).
format_dcsreas <- function(dsdecod) {
  case_when(
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