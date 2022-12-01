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


lb <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "lb.xpt")))
supplb <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "supplb.xpt")))
sv <- admiral::convert_blanks_to_na(haven::read_xpt(file.path("sdtm", "sv.xpt")))

adsl <- admiral::convert_blanks_to_na(haven::read_xpt(envsetup::read_path(adam, "adsl.xpt")))

prodh <- admiral::convert_blanks_to_na(haven::read_xpt(file.path(adam[2], "adlbh.xpt")))
prodc <- admiral::convert_blanks_to_na(haven::read_xpt(file.path(adam[2], "adlbc.xpt")))
prodhy <- admiral::convert_blanks_to_na(haven::read_xpt(file.path(adam[2], "adlbhy.xpt")))
prodhpv <- admiral::convert_blanks_to_na(haven::read_xpt(file.path(adam[2], "adlbhpv.xpt")))
prodcpv <- admiral::convert_blanks_to_na(haven::read_xpt(file.path(adam[2], "adlbcpv.xpt")))



toprogram <- setdiff(colnames(prodh), c(colnames(lb), unique(supplb[["QNAM"]])))

# Formats -----------------------------------------------------------------

format_siteid <- function(x) {
  dplyr::case_when(
    x %in% c("702", "706", "707", "711", "714", "715", "717") ~ "900",
    TRUE ~ x
  )
}

# Add supplemental information --------------------------------------------

sup <- supplb %>%
  dplyr::select(STUDYID, USUBJID, IDVAR, IDVARVAL, QNAM, QLABEL, QVAL) %>%
  tidyr::pivot_wider(id_cols  = c(STUDYID, USUBJID, IDVARVAL),
                     names_from = QNAM,
                     values_from = QVAL) %>%
  dplyr::mutate(LBSEQ = as.numeric(IDVARVAL)) %>%
  dplyr::select(-IDVARVAL)

adlb00 <- lb %>%
  dplyr::left_join(sup, by = c("STUDYID", "USUBJID", "LBSEQ")) 

# ADSL information --------------------------------------------------------

adsl <- adsl %>%
  dplyr::select(STUDYID, USUBJID, TRT01PN, TRT01P, TRT01AN, TRT01A, TRTSDT, TRTEDT, AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX,
                COMP24FL, DSRAEFL, SAFFL) %>%
  dplyr::mutate(SUBJID = sub(".*-", "", USUBJID))


adlb01 <- adlb00 %>%
  dplyr::left_join(adsl, by = c("STUDYID", "USUBJID"))

# Dates -------------------------------------------------------------------
# x <- sapply(lb$LBDTC, FUN = nchar)
# x[x!=16]

adlb02 <- adlb01 %>%
  admiral::derive_vars_dtm(
    new_vars_prefix = "A",
    dtc = LBDTC,
    highest_imputation = "s", # admiral assumes seconds are present before populating ADTM
    ignore_seconds_flag = T
  ) %>%
  admiral::derive_vars_dt(
    new_vars_prefix = "A",
    dtc = LBDTC,
    highest_imputation = "n"
  ) %>%
  admiral::derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))


# adlb02[which(nchar(adlb02$LBDTC) == 16), c("LBDTC", "ADTM")]
# adlb02[which(nchar(adlb02$LBDTC) == 16), c("LBDTC", "ADTM", "ATMF")]
# adlb02[which(nchar(adlb02$LBDTC) == 10), c("LBDTC", "ADTM", "ADT")]
# adlb02[which(nchar(adlb02$LBDTC) != 10 & nchar(adlb02$LBDTC) != 16), c("LBDTC", "ADTM", "ADT")]

# AVAL(C) -----------------------------------------------------------------
# No imputations are done for values below LL or above UL

adlb03 <- adlb02 %>%
  dplyr::mutate(AVAL = LBSTRESN,
                AVALC = ifelse(is.na(AVAL), LBSTRESC, NA))

# Parameter ---------------------------------------------------------------

h <- unique(prodh[, c("PARAMCD", "PARAM", "PARAMN")])
c <- unique(prodc[, c("PARAMCD", "PARAM", "PARAMN")])
hpv <- unique(prodhpv[, c("PARAMCD", "PARAM", "PARAMN")])
hy <- unique(prodhy[, c("PARAMCD", "PARAM", "PARAMN")])
cpv <- unique(prodcpv[, c("PARAMCD", "PARAM", "PARAMN")])

params <- h %>%
  dplyr::bind_rows(c, hpv, hy, cpv) %>%
  arrange(PARAMN)

dput(params)

PARAM For PARAMN<100: same as LB.LBTEST concatenated with LB.LBSTRESU. 
      For PARAMN>100: LB.LBTEST concatenated with LB.LBSTRESU concatenated with change from previous visit, relative to normal range

PARAMCD For PARAMN<100: same as LBTESTCD
        For PARAMN>100: _ + LBTESTCD
        
PARAMN Numeric code for Parameter

# Baseline ----------------------------------------------------------------

bsl <- adlb02 %>%
  dplyr::filter(ADT <= TRTSDT) %>%
  dplyr::arrange(STUDYID, USUBJID, LBTESTCD, ADT, ADTM, LBSEQ) %>%
  dplyr::group_by(STUDYID, USUBJID, LBTESTCD) %>%
  dplyr::slice(n()) %>%
  dplyr::select(STUDYID, USUBJID, LBTESTCD, LBSEQ) %>%
  dplyr::mutate(BASE = AVAL)


"AVISIT" Baseline, Week x, End of Treatment (duplicate where endpoint = "Y")
"AVISITN" 


"PARCAT1"  

"AVAL"     "BASE"     "CHG"      "A1LO"     "A1HI"     "R2A1LO"   "R2A1HI"   "BR2A1LO"  "BR2A1HI" 
"ANL01FL"  "ALBTRVAL" "ANRIND"   "BNRIND"   "ABLFL"    "AENTMTFL"

PARCAT1 Parameter
Category 1
text 5 ADLBCAT harcoded to 'CHEM'
AVAL Analysis Value float 8
BASE Baseline Value float 8
CHG Change from
Baseline
float 8
A1LO Analysis Range 1
Lower Limit
float 8 LB.LBSTNRLO
Variable Label Type Length Display
Format
Code List /
Controlled Terms
Source/Derivation/Comments
A1HI Analysis Range 1
Upper Limit
float 8 LB.LBSTNRHI
R2A1LO Ratio to Analysis
Range 1 Lower
Limit
float 8 AVAL / A1HI
R2A1HI Ratio to Analysis
Range 1 Upper
Limit
float 8 AVAL / A1LO
BR2A1LO Base Ratio to
Analysis Range 1
Lower Lim
float 8 AVAL / A1HI at baseline
BR2A1HI Base Ratio to
Analysis Range 1
Upper Lim
float 8 AVAL / A1LO at baseline
ANL01FL Analysis Record
Flag 1
text 1 Y_BLANK
ALBTRVAL Amount
Threshold Range
float 8 Maximum of [LBSTRESN-(1.5*ULN)] and [(.5*LLN) -
LBSTRESN]
ANRIND Analysis
Reference Range
Indicator
text 1 LNH
Variable Label Type Length Display
Format
Code List /
Controlled Terms
Source/Derivation/Comments
BNRIND Baseline
Reference Range
Indicator
text 1 LNH
ABLFL Baseline Record
Flag
text 1 Y_BLANK
AENTMTFL Last value in
treatment visit
text 1 Y_BLANK Last observed value for this lab parameter during
treatment phase: 'Y' if VISITNUM=12, if subject
discontinues prior to VISIT 12, then this variable is set to
'Y' if this is the last assessment of this analyte for the
subject
LBSEQ Sequence
Number
integer 8 LB.LBSEQ
LBNRIND Reference Range
Indicator
text 8 LBNRIND LB.LBNRIND
LBSTRESN Numeric
Result/Finding in
Standard Units
float 8 LB.LBSTRESN

