

<div id="badges"><!-- pkgdown markup -->
<a href="https://CRAN.R-project.org/web/checks/check_results_calmate.html"><img border="0" src="https://www.r-pkg.org/badges/version/calmate" alt="CRAN check status"/></a> <a href="https://github.com/HenrikBengtsson/calmate/actions?query=workflow%3AR-CMD-check"><img border="0" src="https://github.com/HenrikBengtsson/calmate/actions/workflows/R-CMD-check.yaml/badge.svg?branch=develop" alt="R CMD check status"/></a>     <a href="https://app.codecov.io/gh/HenrikBengtsson/calmate"><img border="0" src="https://codecov.io/gh/HenrikBengtsson/calmate/branch/develop/graph/badge.svg" alt="Coverage Status"/></a> 
</div>

# calmate: Improved Allele-Specific Copy Number of SNP Microarrays for Downstream Segmentation 

The CalMaTe method calibrates preprocessed allele-specific copy number estimates (ASCNs) from DNA microarrays by controlling for single-nucleotide polymorphism-specific allelic crosstalk. The resulting ASCNs are on average more accurate, which increases the power of segmentation methods for detecting changes between copy number states in tumor studies including copy neutral loss of heterozygosity. CalMaTe applies to any ASCNs regardless of preprocessing method and microarray technology, e.g. Affymetrix and Illumina.


## Citing calmate

Whenever using the **calmate** package or the CalMaTe method, please cite:

* M. Ortiz-Estevez, A. Aramburu, H. Bengtsson, P. Neuvial, and A. Rubio, CalMaTe: A method and software to improve allele-specific copy number of SNP arrays for downstream segmentation, Bioinformatics, Volume 28, Issue 13, 2012, Pages 1793-1794, [doi:10.1093/bioinformatics/bts248](https://doi.org/10.1093/bioinformatics/bts248)

## Installation
R package calmate is available on [CRAN](https://cran.r-project.org/package=calmate) and can be installed in R as:
```r
install.packages("calmate")
```


### Pre-release version

To install the pre-release version that is available in Git branch `develop` on GitHub, use:
```r
remotes::install_github("HenrikBengtsson/calmate", ref="develop")
```
This will install the package from source.  

