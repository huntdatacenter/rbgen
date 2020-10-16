## Default behavior according to R extension manual
files <- Sys.glob(paste0("*", SHLIB_EXT))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
files <- Sys.glob(paste0("*", ".a"))
dest <- file.path(R_PACKAGE_DIR, paste0('libs', R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(files, dest, overwrite = TRUE)
if(file.exists("symbols.rds"))
    file.copy("symbols.rds", dest, overwrite = TRUE)

## Add headers
headers <- list.files("include", pattern = "\\.(hpp|h|a)$", recursive = TRUE, full.names = TRUE, include.dirs = TRUE, no.. = TRUE)
if (any(file.exists(headers))){
    for (header in headers) {
      dest <- dirname(file.path(R_PACKAGE_DIR, header))
      dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      file.copy(from = header, to = dest, overwrite = TRUE, recursive = TRUE)
    }
}
