# Project-local R library; analysis uses ChemmineR, without a Python/Java backend.
prepare_R_chemistry_runtime <- function(install_missing=FALSE) {
 lib <- file.path('data','raw','chemistry24','R_library')
 dir.create(lib,recursive=TRUE,showWarnings=FALSE)
 .libPaths(c(normalizePath(lib),.libPaths()))
 if(!requireNamespace('ChemmineR',quietly=TRUE)) {
  if(!install_missing) stop('ChemmineR missing. Run 24b with its authorised install step.')
  repos <- c(BioC='https://bioconductor.org/packages/3.23/bioc',CRAN='https://cloud.r-project.org')
  meta <- available.packages(repos=repos,type='source')
  deps <- tools::package_dependencies('ChemmineR',db=meta,
   which=c('Depends','Imports','LinkingTo'),recursive=TRUE)[[1]]
  missing <- deps[!vapply(deps,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))]
  # Official binaries include compiled graphics libraries needed by R imports.
  if(Sys.info()[['sysname']]=='Darwin' && length(missing)) {
   binary <- available.packages(repos=repos[['CRAN']],type='binary')
   available <- intersect(missing,rownames(binary))
   if(length(available)) install.packages(available,lib=lib,type='binary',
    dependencies=NA,repos=repos[['CRAN']])
  }
  install.packages('ChemmineR',lib=lib,type='source',dependencies=NA,repos=repos)
 }
 if(!requireNamespace('ChemmineR',quietly=TRUE)) stop('ChemmineR installation failed; inspect install messages.')
 invisible(lib)
}
