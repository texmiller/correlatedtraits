#' @export
pdf_doc = function() {
  devtools::document()

  pkg <- "correlatedtraits"

  path <- find.package(pkg)

  if (paste0(pkg, ".pdf") %in% list.files(path)) {
   file.remove(paste0(path, paste0("\\", pkg, ".pdf")))
  }

  system(paste(shQuote(file.path(R.home("bin"), "R")),
               "CMD", "Rd2pdf", shQuote(path)))
}
