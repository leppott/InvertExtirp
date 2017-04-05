NEWS-XC95
================

<!-- NEWS.md is generated from NEWS.Rmd. Please edit that file -->
    #> Last Update: 2017-04-05 15:43:41

Version history.

Planned Updates
===============

-   Update several of the functions and add documentation for the package.

-   Add (working) examples along with data.

-   Functions left to create working examples: fish.gam, fish.wt.cdf, taxon.response, taxon.response.sort, tolerance

v0.0.0.9002
===========

2017-04-05

-   Rename package from InvertExtirp to XC95. Update DESCRIPTION, NEWS, and README.

-   Add data to use with examples.

-   Added raw data for weightcdf(). Included script for creating the RDA files.

-   Added documentation for the above data files (ss and bio.sample).

v0.0.0.9001
===========

2017-04-04

-   Update the R files for readability, documentation, and dependant packages for the functions. "tolerance\_cov.R" requires the "reshape" package (update one call to reshape::cast and add to DESCRIPTION).

-   Rename "tolerance\_cov.R" to "tolerance.R" to match the function.

v0.0.0.9000
===========

2017-03-29

-   Created GitHub repository.

-   Added ReadME.RMD and NEWS.RMD

-   7 functions; curve.shape, fish.gam, fish.wt.cdf, taxon.response, taxon.response.sort, tolerance\_cov, weightedcdf2.
