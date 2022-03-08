## Avoid R CMD check NOTEs on unknown words in package description
##
## Usage:
## saveRDS(readLines(".aspell/WORDLIST"), file = ".aspell/WORDLIST.rds", version = 2L)
##
## References:
## 1. help("aspell-utils", package = "utils")
## 2. https://dirk.eddelbuettel.com/blog/2017/08/10/
Rd_files <- vignettes <- R_files <- description <-
    list(encoding = "UTF-8",
         language = "en",
         dictionaries = c("en_stats", "WORDLIST"))
