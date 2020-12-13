


# calmate: Improved Allele-Specific Copy Number of SNP Microarrays for Downstream Segmentation


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

## Contributions

This Git repository uses the [Git Flow](https://nvie.com/posts/a-successful-git-branching-model/) branching model (the [`git flow`](https://github.com/petervanderdoes/gitflow-avh) extension is useful for this).  The [`develop`](https://github.com/HenrikBengtsson/calmate/tree/develop) branch contains the latest contributions and other code that will appear in the next release, and the [`master`](https://github.com/HenrikBengtsson/calmate) branch contains the code of the latest release, which is exactly what is currently on [CRAN](https://cran.r-project.org/package=calmate).

Contributing to this package is easy.  Just send a [pull request](https://help.github.com/articles/using-pull-requests/).  When you send your PR, make sure `develop` is the destination branch on the [calmate repository](https://github.com/HenrikBengtsson/calmate).  Your PR should pass `R CMD check --as-cran`, which will also be checked by <a href="https://travis-ci.org/HenrikBengtsson/calmate">Travis CI</a> and <a href="https://ci.appveyor.com/project/HenrikBengtsson/calmate">AppVeyor CI</a> when the PR is submitted.

We abide to the [Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/) of Contributor Covenant.


## Software status

| Resource      | CRAN        | GitHub Actions      | Travis CI       | AppVeyor CI      |
| ------------- | ------------------- | ------------------- | --------------- | ---------------- |
| _Platforms:_  | _Multiple_          | _Multiple_          | _Linux & macOS_ | _Windows_        |
| R CMD check   | <a href="https://cran.r-project.org/web/checks/check_results_calmate.html"><img border="0" src="http://www.r-pkg.org/badges/version/calmate" alt="CRAN version"></a> |        | <a href="https://travis-ci.org/HenrikBengtsson/calmate"><img src="https://travis-ci.org/HenrikBengtsson/calmate.svg" alt="Build status"></a>   | <a href="https://ci.appveyor.com/project/HenrikBengtsson/calmate"><img src="https://ci.appveyor.com/api/projects/status/github/HenrikBengtsson/calmate?svg=true" alt="Build status"></a> |
| Test coverage |                     |                     | <a href="https://codecov.io/gh/HenrikBengtsson/calmate"><img src="https://codecov.io/gh/HenrikBengtsson/calmate/branch/develop/graph/badge.svg" alt="Coverage Status"/></a>     |                  |
