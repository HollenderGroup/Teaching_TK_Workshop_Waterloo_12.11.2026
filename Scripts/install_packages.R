## TK Workshop - install the required R packages (run once)
## Open this file in RStudio and click "Source".

pkgs <- c("readxl", "writexl", "deSolve", "rstudioapi")

if (getRversion() < "4.3.0") {
  warning("R ", getRversion(), " is older than recommended (4.3 or newer). ",
          "Please update R from https://cran.r-project.org if you run into problems.")
}

missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
}

ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (all(ok)) {
  message("✔ All set! R ", getRversion(), " with ", paste(pkgs, collapse = ", "), ".")
} else {
  stop("Could not install: ", paste(pkgs[!ok], collapse = ", "),
       ". Please contact the instructor before the workshop.")
}
