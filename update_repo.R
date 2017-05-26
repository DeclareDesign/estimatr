# derived from the drat package
# License: GPL (>= 2)
# https://github.com/eddelbuettel/drat

getPathForPackage <- function(file) {
  pkgtype <- identifyPackageType(file)
  fields <- getPackageInfo(file)
  rversion <- unname(fields["Rmajor"])

  if (pkgtype == "source") {
    ret <- file.path("src", "contrib")
  } else if (pkgtype == "win.binary") {
    ret <- file.path("bin", "windows", "contrib", rversion)
  } else if (pkgtype == "mac.binary") {
    if (fields["OSflavour"] == "") {
      # non-binary package, treated as Mavericks
      message("Note: Non-binary OS X package will be installed in Mavericks path.")
      fields["Mavericks"] <- "yes"
    }
    if (unname(fields["Mavericks"]) == "yes") {
      ret <- file.path("bin", "macosx", "mavericks", "contrib", rversion)
    } else {
      ret <- file.path("bin", "macosx", "contrib", rversion)
    }
  }
  return(ret)
}

identifyPackageType <- function(file) {
  ##from src/library/tools/R/packages.R
  ret <- if (grepl("_.*\\.tar\\..*$", file)) {
    "source"
  } else if (grepl("_.*\\.tgz$", file)) {
    "mac.binary"
  } else if (grepl("_.*\\.zip$", file)) {
    "win.binary"
  } else {
    stop("Unknown package type", call. = FALSE)
  }
  return(ret)
}

path <- ifelse(
  .Platform$OS.type == 'windows',
  file.path('..', '${APPVEYOR_PROJECT_NAME:-$PKG_REPO}'),
  file.path('..')
)

files <-
  dir(path, pattern = ifelse(.Platform$OS.type == 'windows', '.zip', '.t*z'))

print(paste('processing', pkg))

for (f in files) {
  write_PACKAGES(pkgdir, type = identifyPackageType(f), ...)
}
