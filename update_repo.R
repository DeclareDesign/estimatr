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
      fields["El-Capitan"] <- "yes"
    }
    if (unname(fields["El-Capitan"]) == "yes") {
      ret <- file.path("bin", "macosx", "el-capitan", "contrib", rversion)
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
  pkgdir <- file.path(repodir, reldir)

  if (!file.exists(pkgdir)) {
    ## TODO: this could be in a git branch, need checking
    if (!dir.create(pkgdir, recursive = TRUE)) {
      stop("Directory ", pkgdir, " couldn't be created\n", call. = FALSE)
    }
  }

  ## copy file into repo
  if (!file.copy(f, pkgdir, overwrite = TRUE)) {
    stop("File ", f, " can not be copied to ", pkgdir, call. = FALSE)
  }

  write_PACKAGES(pkgdir, type = identifyPackageType(f), ...)

}
