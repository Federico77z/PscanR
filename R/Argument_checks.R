# Assisted-by: OpenAI Codex and Claude Code (code refactoring, review and
# documentation). All changes were reviewed and tested by the authors.

# Shared argument checks for the exported functions. Each one stops with the
# same kind of message, "'<argument>' must be ...", and without the internal
# call, so the user sees the argument they passed rather than the helper that
# checked it. Checks specific to one input (a background table, a background
# file, a full background) stay next to the code that needs them.

.ps_assert_class <- function(x, classes, arg) {
    if (!any(vapply(classes, function(cl) is(x, cl), logical(1)))) {
        msg <- sprintf("'%s' must be a %s", arg,
            paste(classes, collapse = " or "))
        stop(msg, call. = FALSE)
    }
    invisible(x)
}

.ps_assert_character <- function(x, arg) {
    if (!is.character(x)) {
        stop(sprintf("'%s' must be a character vector", arg), call. = FALSE)
    }
    invisible(x)
}

.ps_assert_string <- function(x, arg) {
    if (!is.character(x) || length(x) != 1L || is.na(x)) {
        stop(sprintf("'%s' must be a single string", arg), call. = FALSE)
    }
    invisible(x)
}

.ps_assert_flag <- function(x, arg) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(sprintf("'%s' must be TRUE or FALSE", arg), call. = FALSE)
    }
    invisible(x)
}

.ps_assert_number <- function(x, arg, min = -Inf, max = Inf) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x < min ||
            x > max) {
        bounds <- if (is.finite(min) && is.finite(max)) {
            sprintf(" between %s and %s", min, max)
        } else if (is.finite(min)) {
            sprintf(" of at least %s", min)
        } else if (is.finite(max)) {
            sprintf(" of at most %s", max)
        } else {
            ""
        }
        msg <- sprintf("'%s' must be a single number%s", arg, bounds)
        stop(msg, call. = FALSE)
    }
    invisible(x)
}

# The output file can be created, or overwritten if it exists.
.ps_assert_writable_file <- function(file, arg) {
    .ps_assert_string(file, arg)
    if (file.exists(file)) {
        if (file.access(file, mode = 2) != 0) {
            stop(sprintf("Cannot write to the existing file %s", file),
                call. = FALSE)
        }
    } else if (file.access(dirname(file), mode = 2) != 0) {
        stop(sprintf("Cannot write to the directory %s", dirname(file)),
            call. = FALSE)
    }
    invisible(file)
}
