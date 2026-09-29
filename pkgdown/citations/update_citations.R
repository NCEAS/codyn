# Compile the list of publications that cite codyn and write it as BibTeX.
#
# Reads the DOIs and filters from pkgdown/citations/citation-sources.yml, asks
# OpenAlex (https://openalex.org) for every work citing any of those DOIs, adds
# the papers describing codyn themselves (but not the software releases), drops
# false matches and preprints of papers that were later published, and writes the
# result to the BibTeX file named in the config. The file is only rewritten
# when the set of publications changes, so scheduled runs produce no diff when
# nothing new has been published.
#
# Run from the package root:
#   Rscript pkgdown/citations/update_citations.R
#
# Optional environment variables:
#   OPENALEX_API_KEY   OpenAlex API key, sent with each request if set
#   OPENALEX_MAILTO    contact email for the OpenAlex polite pool
#   CITATION_SUMMARY   path of a file to write a plain-text summary of changes
#                      to, used as the body of the workflow's commit message

library(yaml)
library(jsonlite)

config_path <- "pkgdown/citations/citation-sources.yml"
api <- "https://api.openalex.org"

# Types are ranked when choosing one record among works with the same title,
# so that a published article wins over its preprint.
type_rank <- c("article", "review", "book-chapter", "book", "conference-paper",
               "dissertation", "report", "dataset", "other", "preprint")

openalex_get <- function(path, query = list(), tries = 4) {
  key <- Sys.getenv("OPENALEX_API_KEY")
  mailto <- Sys.getenv("OPENALEX_MAILTO")
  if (nzchar(key)) query$api_key <- key
  if (nzchar(mailto)) query$mailto <- mailto
  url <- paste0(api, path)
  if (length(query) > 0) {
    url <- paste0(url, "?", paste(names(query),
                                  vapply(query, utils::URLencode, "", reserved = TRUE),
                                  sep = "=", collapse = "&"))
  }
  for (i in seq_len(tries)) {
    res <- tryCatch(fromJSON(url, simplifyVector = FALSE), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    if (i < tries) Sys.sleep(2^i)
  }
  stop("OpenAlex request failed: ", url, "\n", conditionMessage(res))
}

work_fields <- paste(c("id", "doi", "title", "publication_year", "type",
                       "authorships", "primary_location", "biblio"), collapse = ",")

fetch_doi <- function(doi) {
  openalex_get(paste0("/works/doi:", doi), list(select = work_fields))
}

fetch_citing <- function(ids) {
  works <- list()
  cursor <- "*"
  while (!is.null(cursor)) {
    page <- openalex_get("/works", list(
      filter = paste0("cites:", paste(ids, collapse = "|")),
      select = work_fields,
      `per-page` = "200",
      cursor = cursor
    ))
    works <- c(works, page$results)
    cursor <- page$meta$next_cursor
  }
  works
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || identical(x, "")) y else x

# Convert an OpenAlex title (which may contain HTML markup and entities) to
# BibTeX-safe text, keeping italics.
clean_text <- function(x) {
  if (is.null(x)) return("")
  x <- gsub("<(i|em)>(.*?)</\\1>", "\u0001\\2\u0002", x, perl = TRUE, ignore.case = TRUE)
  x <- gsub("<[^>]+>", "", x)
  entities <- c("&amp;" = "&", "&lt;" = "<", "&gt;" = ">", "&quot;" = "\"",
                "&#39;" = "'", "&apos;" = "'", "&nbsp;" = " ")
  for (e in names(entities)) x <- gsub(e, entities[[e]], x, fixed = TRUE)
  x <- gsub("[{}\\\\]", "", x)
  x <- gsub("([&%$#_])", "\\\\\\1", x)
  x <- gsub("\u0001", "\\\\emph{", x)
  x <- gsub("\u0002", "}", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

title_key <- function(x) gsub("[^a-z0-9]", "", tolower(x %||% ""))

ascii_name <- function(x) {
  x <- iconv(x, "UTF-8", "ASCII//TRANSLIT", sub = "")
  gsub("[^A-Za-z]", "", x)
}

work_to_bibtex <- function(w, corrections = list()) {
  id <- sub("^https://openalex.org/", "", w$id)
  authors <- vapply(w$authorships, function(a) {
    clean_text(a$raw_author_name %||% a$author$display_name %||% "")
  }, "")
  authors <- authors[nzchar(authors)]
  first_last <- if (length(authors)) ascii_name(tail(strsplit(authors[1], " ")[[1]], 1)) else ""
  key <- paste0(first_last %||% "anon", w$publication_year, "_", id)

  source <- clean_text(w$primary_location$source$display_name)
  b <- w$biblio
  pages <- if (!is.null(b$first_page)) {
    if (!is.null(b$last_page) && b$last_page != b$first_page) {
      paste0(b$first_page, "--", b$last_page)
    } else {
      b$first_page
    }
  }
  entry_type <- switch(w$type,
    article = , review = if (nzchar(source)) "article" else "misc",
    preprint = if (nzchar(source)) "article" else "misc",
    "book-chapter" = "incollection",
    book = "book",
    "conference-paper" = "inproceedings",
    dissertation = "phdthesis",
    report = "techreport",
    "misc"
  )
  container <- switch(entry_type,
    article = "journal", incollection = , inproceedings = "booktitle",
    phdthesis = "school", techreport = "institution", "howpublished")

  fields <- list(
    author = paste(authors, collapse = " and "),
    title = paste0("{", clean_text(w$title), "}"),
    year = w$publication_year,
    volume = b$volume,
    number = b$issue,
    pages = pages,
    doi = sub("^https://doi.org/", "", w$doi %||% ""),
    url = if (is.null(w$doi)) w$primary_location$landing_page_url
  )
  fields[[container]] <- source
  fix <- corrections[[id]]
  if (!is.null(fix$title)) fix$title <- paste0("{", clean_text(fix$title), "}")
  fields[names(fix)] <- fix
  fields <- Filter(function(v) !is.null(v) && nzchar(v), lapply(fields, as.character))
  body <- paste0("  ", names(fields), " = {", unlist(fields), "}", collapse = ",\n")
  list(key = key, year = w$publication_year, title = clean_text(fix$title %||% w$title),
       text = paste0("@", entry_type, "{", key, ",\n", body, "\n}"))
}

main <- function() {
  config <- read_yaml(config_path)
  dois <- vapply(config$dois, function(d) d$doi, "")
  kinds <- vapply(config$dois, function(d) d$kind %||% "", "")
  if (!all(kinds %in% c("paper", "software"))) {
    stop("Each DOI in ", config_path, " needs `kind: paper` or `kind: software`: ",
         paste(dois[!kinds %in% c("paper", "software")], collapse = ", "))
  }

  message("Resolving ", length(dois), " DOIs in OpenAlex")
  sources <- lapply(dois, fetch_doi)
  source_ids <- sub("^https://openalex.org/", "", vapply(sources, function(w) w$id, ""))
  for (i in seq_along(dois)) message("  ", dois[i], " (", kinds[i], ") -> ", source_ids[i])

  # The papers describing codyn are listed along with the works citing them;
  # software releases are only used to find citing works.
  citing <- fetch_citing(source_ids)
  message("Retrieved ", length(citing), " citing works")
  citing_ids <- sub("^https://openalex.org/", "", vapply(citing, function(w) w$id, ""))
  works <- c(sources[kinds == "paper"], citing[!citing_ids %in% source_ids])

  ids <- sub("^https://openalex.org/", "", vapply(works, function(w) w$id, ""))
  type <- vapply(works, function(w) w$type %||% "other", "")
  year <- vapply(works, function(w) as.integer(w$publication_year %||% NA), 1L)
  keep <- !(ids %in% unlist(config$exclude_works)) &
    !(type %in% unlist(config$exclude_types)) &
    !is.na(year) & year >= config$min_year
  message("Dropped ", sum(!keep), " works that are excluded or too old")
  works <- works[keep]
  type <- type[keep]

  # Keep one record per title, preferring published versions over preprints
  # and records with a DOI.
  rank <- match(type, type_rank, nomatch = length(type_rank) - 1)
  has_doi <- vapply(works, function(w) !is.null(w$doi), TRUE)
  tkey <- vapply(works, function(w) title_key(w$title), "")
  ord <- order(tkey, rank, !has_doi)
  works <- works[ord]
  dup <- duplicated(tkey[ord]) & nzchar(tkey[ord])
  message("Dropped ", sum(dup), " duplicate records (mostly preprints of published papers)")
  works <- works[!dup]

  entries <- lapply(works, work_to_bibtex, corrections = config$corrections %||% list())
  entries <- entries[order(-vapply(entries, `[[`, 1, "year"),
                           vapply(entries, `[[`, "", "key"))]
  keys <- vapply(entries, `[[`, "", "key")
  body <- paste(vapply(entries, `[[`, "", "text"), collapse = "\n\n")

  out <- config$output
  old_keys <- character(0)
  if (file.exists(out)) {
    old <- readLines(out, encoding = "UTF-8")
    old_keys <- sub("^@[a-z]+\\{(.*),$", "\\1", grep("^@", old, value = TRUE))
    old_body <- paste(old[!grepl("^%", old)], collapse = "\n")
    if (identical(trimws(old_body), trimws(body))) {
      message("No changes; ", out, " is up to date with ", length(entries), " publications")
      return(invisible())
    }
  }

  header <- c(
    "% Publications citing codyn, compiled from OpenAlex (https://openalex.org).",
    "% Generated by pkgdown/citations/update_citations.R from pkgdown/citations/citation-sources.yml;",
    "% do not edit by hand.",
    paste0("% Sources: ", paste(dois, collapse = ", ")),
    paste0("% Retrieved: ", format(Sys.Date())),
    paste0("% Entries: ", length(entries)),
    ""
  )
  dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
  con <- file(out, open = "w", encoding = "UTF-8")
  writeLines(c(header, body), con)
  close(con)

  added <- entries[!keys %in% old_keys]
  removed <- setdiff(old_keys, keys)
  message("Wrote ", length(entries), " publications to ", out, " (",
          length(added), " added, ", length(removed), " removed)")

  summary_path <- Sys.getenv("CITATION_SUMMARY")
  if (nzchar(summary_path)) {
    lines <- c(
      paste0("Total publications: ", length(entries), " (", length(added),
             " added, ", length(removed), " removed)"),
      ""
    )
    if (length(added)) {
      lines <- c(lines, "Added:", vapply(added, function(e) {
        paste0("- ", e$year, ": ", e$title, " (", sub(".*_", "", e$key), ")")
      }, ""), "")
    }
    if (length(removed)) {
      lines <- c(lines, "Removed:", paste0("- ", removed), "")
    }
    lines <- c(lines, "To remove a false match, add its OpenAlex ID to exclude_works",
               "in pkgdown/citations/citation-sources.yml.")
    writeLines(lines, summary_path)
  }
}

main()
