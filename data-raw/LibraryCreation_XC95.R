# Helper Script when working on creating a Library
# moves files to proper location
# creates documentation and installs the new library
# Erik.Leppo@tetratech.com
# 20170223
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NEWS
# Render then Copy NEWS so picked up in help
rmarkdown::render("NEWS.rmd", "all")
file.copy("NEWS.md", "NEWS", overwrite = TRUE)
#file.remove("NEWS.html")
#file.remove("NEWS.md")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Library Name
myLibrary <- "XC95"
# Load Library
library(devtools)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create Package
# create(myLibrary)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Document, Install, and Reload Library
## Generate Documentation
setwd(paste0("./", myLibrary))
devtools::document()
## Install New Package (locally)
setwd("..") # return to root directory first
devtools::install(myLibrary)
## Reload library
library(myLibrary, character.only = TRUE)
# change wd back to package
setwd(paste0("./", myLibrary))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## Restart R within RStudio:  Ctrl + Shift + F10
library("XC95")
help(package="XC95")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Upload to Github via GitHub Desktop utility
# 0. download from web via "clone or download" via "Open in Desktop" (GitHub Desktop) if not already in GitHub Desktop
# 1. Make changes in download/clone folder. (done above)
# 3. Open GH Desktop commit changes then sync.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# install from GitHub (via devtools)
devtools::install_github(paste0("leppott/",myLibrary))
#


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remove installed packages (if needed for troubleshooting)
search() # find
#detach(3) # remove by number
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# to build package
#https://thepoliticalmethodologist.com/2014/08/14/building-and-maintaining-r-packages-with-devtools-and-roxygen2/
# To build the package as a compressed file in your working directory, run build(current.code, path=getwd()).

# to save internal data for examples
# example
#http://r-pkgs.had.co.nz/data.html#data-sysdata
# have to be at root directory (above package)
#devtools::use_data(NV.predictors,NV.bugs,pkg="MMIcalcNV",internal=TRUE,overwrite=TRUE)
## verify with data()

# To save RMD files
# http://stackoverflow.com/questions/30377213/how-to-include-rmarkdown-file-in-r-package
# /pkg/inst/rmd/
# system.file("rmd/file.Rmd", package="packagename")
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
# Create Package
# create(myLibrary)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Travis_CI ####
# https://juliasilge.com/blog/beginners-guide-to-travis/

# Add yaml file
devtools::use_travis()
