###########################################################################
#' developers : Thomas Neitmann/
#' date: 16DEC2022
#' modification History:
#' QC ADAE
###########################################################################
library(haven)
library(diffdf)

adsl <- read_xpt(file.path("submission", "datasets", "adsl.xpt"))
qc_adsl <- read_xpt(file.path("adam", "adsl.xpt"))

diffdf(adsl, qc_adsl, keys = c("STUDYID", "USUBJID"))


# -------------------------------------------#
# QC/Check against original TDF ADAE dataset #
# -------------------------------------------#
# adsl  <- haven::read_xpt("https://github.com/RConsortium/submissions-pilot3-adam/blob/main/adam/adsl.xpt?raw=true")
# adae_orig <- haven::read_xpt("https://github.com/RConsortium/submissions-pilot3-adam/blob/main/adam/adae.xpt?raw=true") %>% 
#   convert_blanks_to_na() #blanks to NA


adae <- read_xpt(file.path("submission", "datasets", "adae.xpt"))
qc_adsl <- read_xpt(file.path("adam", "adae.xpt")) %>% 
  convert_blanks_to_na() #blanks to NA
#---------#
# Compare #
#---------#
diffdf(adae, adae_orig, keys = c("STUDYID","USUBJID","AESEQ")   ) 
