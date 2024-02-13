#List out packages used for R analysis in the adrg.pdf

library(tidyr)
library(data.table)

pkgloaded <- sessionInfo()$loadedOnly #get intial list of packages that were loaded
pkgother <- sessionInfo()$otherPkgs #get initla list of packages that were also listed in the Session but may not have been used

loaded <- data.frame(rbindlist(pkgloaded, idcol = TRUE, fill = TRUE)) %>%
            select(Package, Title, Version, Description) %>%
            mutate(loaded = "Y")
other <- data.frame(rbindlist(pkgother, idcol = TRUE, fill = TRUE)) %>%
            select(Package, Title, Version, Description) %>%
            mutate(loaded = "N")

pkgdesc <- bind_rows(loaded, other) # stacks all package data frames.
# NOTE column 'loaded', from this data frame can be used to subset out packages not used and may not be needed for the adrg.pdf

# update path to save CSV in your local area.
write.csv(pkgdesc, "/cloud/project/submission/pkgs4adrg.csv", row.names = FALSE)
