# CRAN submission calmate 0.13.0

on 2022-03-08

Thanks in advance


## Notes not sent to CRAN

### R CMD check validation

The package has been verified using `R CMD check --as-cran` on:

| R version     | GitHub | R-hub    | win-builder |
| ------------- | ------ | -------- | ----------- |
| 3.4.x         | L      |          |             |
| 3.6.x         | L      |          |             |
| 4.0.x         | L      | L        |             |
| 4.1.x         | L M W  | L M M1 W | W           |
| devel         | L M W  | L        | W           |

*Legend: OS: L = Linux, M = macOS, M1 = macOS M1, W = Windows*


R-hub checks:

```r
res <- rhub::check(platform = c(
  "debian-clang-devel", "debian-gcc-patched", "linux-x86_64-centos-epel",
  "macos-highsierra-release-cran", "macos-m1-bigsur-release",
  "windows-x86_64-release"))
print(res)
```

passed with all OK.
