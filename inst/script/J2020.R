library("JASPAR2020")

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "vertebrates" # Change Tax Group Based on the organism of study 

J2020 <- getMatrixSet(JASPAR2020, opts) 

# For JASPAR2022

#library("JASPAR2022")

#opts <- list()
#opts[["collection"]] <- "CORE"
#opts[["tax_group"]] <- "vertebrates"

#J2022 <- getMatrixSet(JASPAR2022, opts)

# For JASPAR2024

#library("JASPAR2024")
#library("RSQLite")

#opts <- list()
#opts[["collection"]] <- "CORE"
#opts[["tax_group"]] <- "vertebrates"

#library("httr")
#httr::set_config(config(ssl_verifypeer = 0L))
#JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))
#J2024 <- getMatrixSet(JASPARConnect, opts)