if (Sys.getenv("GITHUB_ACTIONS") == "") {
  source("renv/activate.R")
  envsetup::rprofile(config::get(file = file.path(getwd(), "_envsetup.yml"), config = "prod"))
  source("inst/startup.R")
} else {
  options(repos = "https://cran.microsoft.com/snapshot/2022-11-01")
}
