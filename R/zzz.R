# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE


.First.lib <- function(libname, pkgname) {
  pd <- packageDescription(pkgname);

  packageStartupMessage(pkgname, " v", pd$Version, " (", 
    pd$Date, ") successfully loaded. See ?", pkgname, " for help.");
}
