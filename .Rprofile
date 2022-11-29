if (Sys.getenv("GITHUB_ACTIONS") == "") {
  source("renv/activate.R")
  envsetup::rprofile(config::get(file = file.path(getwd(), "_envsetup.yml")))
  source("inst/startup.R")
}
