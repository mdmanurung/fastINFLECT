## Remove malformed empty-path records emitted by some pkgdown/search-index
## combinations. Run after pkgdown::build_site().

search_path <- file.path("docs", "search.json")
if (file.exists(search_path)) {
  records <- jsonlite::fromJSON(search_path, simplifyVector = FALSE)
  keep <- vapply(records, function(record) {
    is.list(record) &&
      is.character(record$path) &&
      length(record$path) == 1L &&
      nzchar(record$path)
  }, logical(1))

  if (!all(keep)) {
    jsonlite::write_json(
      records[keep],
      search_path,
      auto_unbox = TRUE,
      null = "null",
      pretty = FALSE
    )
  }
}
