args <- commandArgs(trailingOnly = TRUE)

value_after <- function(flag, default = NULL) {
  hit <- which(args == flag)
  if (!length(hit)) {
    return(default)
  }
  idx <- hit[[1]] + 1
  if (idx > length(args)) {
    stop("Missing value after ", flag)
  }
  args[[idx]]
}

has_flag <- function(flag) {
  any(args == flag)
}

input <- value_after("--input")
input_format <- value_after("--input-format", "auto")
object_name <- value_after("--object-name", "")
list_policy <- value_after("--list-policy", "first_table")
output <- value_after("--output")
summary_path <- value_after("--summary")
include_rownames <- has_flag("--include-rownames")

if (is.null(input) || is.null(output) || is.null(summary_path)) {
  stop("Required arguments: --input, --output, --summary")
}

detect_format <- function(path, requested) {
  if (requested != "auto") {
    return(requested)
  }
  lower <- tolower(basename(path))
  if (grepl("\\.rds$", lower)) {
    return("rds")
  }
  if (grepl("\\.(rdata|rda)$", lower)) {
    return("rdata")
  }
  "auto"
}

read_rdata <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  names_loaded <- load(path, envir = env)
  if (!length(names_loaded)) {
    stop("No objects were loaded from RData input")
  }
  if (object_name != "") {
    if (!exists(object_name, envir = env, inherits = FALSE)) {
      stop("Object '", object_name, "' was not found in the RData input")
    }
    return(list(value = get(object_name, envir = env, inherits = FALSE), selected_path = object_name))
  }
  selected_name <- names_loaded[[1]]
  list(value = get(selected_name, envir = env, inherits = FALSE), selected_path = selected_name)
}

read_input <- function(path, fmt, object_name) {
  if (fmt == "rds") {
    return(list(value = readRDS(path), format = "rds", selected_path = NULL))
  }
  if (fmt == "rdata") {
    loaded <- read_rdata(path, object_name)
    return(list(value = loaded$value, format = "rdata", selected_path = loaded$selected_path))
  }

  rds_attempt <- tryCatch(
    list(ok = TRUE, value = readRDS(path)),
    error = function(e) list(ok = FALSE, error = e)
  )
  if (rds_attempt$ok) {
    return(list(value = rds_attempt$value, format = "rds", selected_path = NULL))
  }

  rdata_attempt <- tryCatch(
    list(ok = TRUE, value = read_rdata(path, object_name)),
    error = function(e) list(ok = FALSE, error = e)
  )
  if (rdata_attempt$ok) {
    loaded <- rdata_attempt$value
    return(list(value = loaded$value, format = "rdata", selected_path = loaded$selected_path))
  }

  stop(
    "Input could not be read as RDS or RData. RDS error: ",
    conditionMessage(rds_attempt$error),
    "; RData error: ",
    conditionMessage(rdata_attempt$error)
  )
}

is_rectangular <- function(x) {
  is.data.frame(x) || is.matrix(x) || methods::is(x, "DataFrame")
}

coerce_rectangular <- function(x) {
  if (is.data.frame(x)) {
    return(as.data.frame(x, check.names = FALSE))
  }
  if (is.matrix(x)) {
    return(as.data.frame(x, check.names = FALSE))
  }
  if (methods::is(x, "DataFrame")) {
    out <- tryCatch(as.data.frame(x, check.names = FALSE), error = function(e) NULL)
    if (!is.null(out) && is.data.frame(out)) {
      return(out)
    }
  }
  NULL
}

atomic_list_to_df <- function(x) {
  if (!is.list(x) || !length(x)) {
    return(NULL)
  }
  keep <- vapply(x, function(y) is.atomic(y) && length(dim(y)) == 0, logical(1))
  if (!any(keep)) {
    return(NULL)
  }
  cols <- x[keep]
  lens <- vapply(cols, length, integer(1))
  if (length(unique(lens)) != 1) {
    return(NULL)
  }
  as.data.frame(cols, check.names = FALSE, stringsAsFactors = FALSE)
}

find_first_table <- function(x, path = "root") {
  tab <- coerce_rectangular(x)
  if (!is.null(tab)) {
    attr(tab, "selected_path") <- path
    return(tab)
  }
  if (!is.list(x)) {
    return(NULL)
  }
  x_names <- names(x)
  for (i in seq_along(x)) {
    nm <- if (is.null(x_names)) "" else x_names[[i]]
    child_path <- if (!is.na(nm) && nzchar(nm)) {
      paste0(path, "$", nm)
    } else {
      paste0(path, "[[", i, "]]")
    }
    candidate <- find_first_table(x[[i]], child_path)
    if (!is.null(candidate)) {
      return(candidate)
    }
  }
  NULL
}

select_from_object <- function(x, name) {
  if (name == "") {
    return(x)
  }
  if (is.list(x) && name %in% names(x)) {
    return(x[[name]])
  }
  if (!is.null(names(x)) && name %in% names(x)) {
    return(x[[name]])
  }
  stop("Object/list element '", name, "' was not found")
}

fmt <- detect_format(input, input_format)
input_data <- read_input(input, fmt, object_name)
fmt <- input_data$format

if (fmt == "rdata") {
  selected <- input_data$value
  selected_path <- input_data$selected_path
} else {
  selected <- select_from_object(input_data$value, object_name)
  selected_path <- if (object_name == "") "root" else object_name
}

table <- coerce_rectangular(selected)

if (is.null(table) && is.list(selected) && list_policy == "atomic_columns") {
  table <- atomic_list_to_df(selected)
  selected_path <- if (object_name == "") "root atomic-list columns" else paste(object_name, "atomic-list columns")
}

if (is.null(table) && is.list(selected) && list_policy == "first_table") {
  table <- find_first_table(selected, selected_path)
  selected_path <- attr(table, "selected_path", exact = TRUE)
}

if (is.null(table)) {
  stop("Selected object cannot be converted to a rectangular table")
}

if (include_rownames) {
  rn <- rownames(table)
  if (!is.null(rn) && !identical(rn, as.character(seq_len(nrow(table))))) {
    table <- cbind(rowname = rn, table)
  }
}

write.table(table, file = output, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "")

summary_lines <- c(
  paste("input_format", fmt, sep = "\t"),
  paste("selected_path", selected_path, sep = "\t"),
  paste("object_class", paste(class(selected), collapse = ","), sep = "\t"),
  paste("rows", nrow(table), sep = "\t"),
  paste("columns", ncol(table), sep = "\t"),
  paste("column_names", paste(names(table), collapse = ","), sep = "\t")
)
writeLines(summary_lines, con = summary_path)
