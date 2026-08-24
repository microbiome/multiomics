################################################################################
# integration with other packages

.require_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("'", pkg, "' package not found. Please install the '", pkg, "' package ",
            "to use this function.",
            call. = FALSE
        )
    }
}

################################### TESTING ###################################
# Methods for testing

# This function unifies input testing. The message will always be in same format
# also it makes the code simpler in main function since testing is done here.
# Borrowed from HoloFoodR.
.check_input <- function(variable, supported_class, supported_values = NULL, limits = NULL,
                         length = NULL, variable_name = .get_name_in_parent(variable)) {
    # Convert supported classes to character
    classes_char <- lapply(supported_class, function(class) {
        if (is.null(class)) {
            class <- "NULL"
        }
        return(class)
    })
    classes_char <- unlist(classes_char)
    # Based on number of acceptable classes, the msg is different
    class_txt <- .create_msg_from_list(classes_char)
    # Create a message
    msg <- paste0("'", variable_name, "' must be ", class_txt, ".")

    # If supported values were provided
    if (!is.null(supported_values)) {
        # Convert supported values to character
        values_char <- lapply(supported_values, function(value) {
            if (is.null(value)) {
                value <- "NULL"
            }
            value <- as.character(value)
            return(value)
        })
        values_char <- unlist(values_char)
        # Collapse into text
        values_txt <- paste0("'", paste(values_char, collapse = "', '"), "'")
        msg <- paste0(
            msg, " It must be one of the following options: ", values_txt
        )
    }

    # If limits were provided
    if (!is.null(limits)) {
        msg <- paste0(msg, " (Numeric constrains: ")
        # Add thresholds to message
        if (!is.null(limits$upper)) {
            msg <- paste0(msg, limits$upper, ">x")
        } else if (!is.null(limits$include_upper)) {
            msg <- paste0(msg, limits$upper, ">=x")
        }
        if (!is.null(limits$lower)) {
            msg <- paste0(msg, "x>", limits$lower)
        } else if (!is.null(limits$include_lower)) {
            msg <- paste0(msg, "x>=", limits$include_lower)
        }
        msg <- paste0(msg, ")")
    }

    # If length was provided
    if (!is.null(length)) {
        msg <- paste0(
            msg, " The length must be ",
            paste0(length, collapse = " or "), "."
        )
    }

    # List all the input types. Run the check if the variable must be that type.
    # If correct type was found, change the result to TRUE.
    input_correct <- FALSE
    if ("NULL" %in% classes_char && is.null(variable)) {
        input_correct <- TRUE
    }
    if ("logical scalar" %in% classes_char && .is_a_bool(variable)) {
        input_correct <- TRUE
    }
    if ("logical vector" %in% classes_char && is.logical(variable)) {
        input_correct <- TRUE
    }
    if ("character scalar" %in% classes_char && .is_non_empty_string(
        variable
    )) {
        input_correct <- TRUE
    }
    if ("character vector" %in% classes_char && .is_non_empty_character(
        variable
    )) {
        input_correct <- TRUE
    }
    if ("numeric scalar" %in% classes_char && .is_a_numeric(variable)) {
        input_correct <- TRUE
    }
    if ("numeric vector" %in% classes_char && is.numeric(variable)) {
        input_correct <- TRUE
    }
    if ("integer vector" %in% classes_char && .is_integer(variable)) {
        input_correct <- TRUE
    }
    if ("integer scalar" %in% classes_char && .is_an_integer(variable)) {
        input_correct <- TRUE
    }
    if ("list" %in% classes_char && is.list(variable) && !is.data.frame(
        variable
    )) {
        input_correct <- TRUE
    }
    if ("data.frame" %in% classes_char && is.data.frame(variable)) {
        input_correct <- TRUE
    }
    if ("matrix" %in% classes_char && is.matrix(variable)) {
        input_correct <- TRUE
    }
    # If supported values were provided. Check these only if the variable
    # is not numeric, NULL or list.
    if (!is.null(supported_values) && !is.null(variable) &&
        !is.numeric(variable) && !is.list(variable)) {
        # Test that if variable is in supported values
        values_correct <- lapply(supported_values, function(value) {
            res <- FALSE
            if (is.null(value) && is.null(variable) || value %in% variable) {
                res <- TRUE
            }
            return(res)
        })
        values_correct <- unlist(values_correct)
        # If not, then give FALSE even though class checks were correct
        if (!any(values_correct)) {
            input_correct <- FALSE
        }
    }
    # If limits were provided
    if (!is.null(limits) && !is.null(variable)) {
        if (!is.null(limits$upper) && variable >= limits$upper) {
            input_correct <- FALSE
        } else if (!is.null(
            limits$include_upper
        ) && variable > limits$include_upper) {
            input_correct <- FALSE
        }

        if (!is.null(limits$lower) && variable <= limits$lower) {
            input_correct <- FALSE
        } else if (!is.null(
            limits$include_lower
        ) && variable < limits$include_lower) {
            input_correct <- FALSE
        }
    }
    # Check length if provided
    if (!is.null(variable) && !is.null(length) &&
        !length(variable) %in% length) {
        input_correct <- FALSE
    }
    # Give error if variable was not correct type
    if (!input_correct) {
        stop(msg, call. = FALSE)
    }
    return(input_correct)
}

# This function creates a string from character values provided. The string
# can be used to messages. It creates a tidy list from list of values.
.create_msg_from_list <- function(classes_char, and_or = "or", ...) {
    if (length(classes_char) > 2) {
        class_txt <- paste0(
            paste(
                classes_char[seq_len(length(classes_char) - 1)],
                collapse = ", "
            ),
            " ", and_or, " ", classes_char[length(classes_char)]
        )
    } else if (length(classes_char) == 2) {
        class_txt <- paste0(
            classes_char[[1]], " ", and_or, " ", classes_char[[2]]
        )
    } else {
        class_txt <- classes_char
    }
    return(class_txt)
}

################################################################################
# testing

.is_a_bool <- function(x) {
    is.logical(x) && length(x) == 1L && !is.na(x)
}

.is_non_empty_character <- function(x) {
    is.character(x) && all(nzchar(x))
}

.is_non_empty_string <- function(x) {
    .is_non_empty_character(x) && length(x) == 1L
}

.is_a_string <- function(x) {
    is.character(x) && length(x) == 1L
}

.is_integer <- function(x) {
    is.numeric(x) && all(x %% 1 == 0)
}

.is_an_integer <- function(x) {
    .is_integer(x) && length(x) == 1L
}

.are_whole_numbers <- function(x) {
    tol <- 100 * .Machine$double.eps
    abs(x - round(x)) <= tol && !is.infinite(x)
}

.is_a_numeric <- function(x) {
    is.numeric(x) && length(x) == 1L
}

.is_numeric_string <- function(x) {
    x <- as.character(x)
    suppressWarnings({
        x <- as.numeric(x)
    })
    !is.na(x)
}

.is_function <- function(x) {
    typeof(x) == "closure" && is(x, "function")
}

.all_are_existing_files <- function(x) {
    all(file.exists(x))
}

.get_name_in_parent <- function(x) {
    .safe_deparse(do.call(substitute, list(substitute(x), parent.frame())))
}

.safe_deparse <- function(expr, ...) {
    paste0(deparse(expr, width.cutoff = 500L, ...), collapse = "")
}

################################################################################
# checks

#' @importFrom SummarizedExperiment assays
.check_assay_present <- function(assay.type, x, name = .get_name_in_parent(assay.type)) {
    if (!.is_non_empty_string(assay.type)) {
        stop("'", name, "' must be a single non-empty character value.",
            call. = FALSE
        )
    }
    if (!(assay.type %in% names(assays(x)))) {
        stop("'", name, "' must be a valid name of assays(x)", call. = FALSE)
    }
}

.check_rowTree_present <- function(tree.name, x, name = .get_name_in_parent(tree.name)) {
    if (!.is_non_empty_string(tree.name)) {
        stop("'", name, "' must be a single non-empty character value.",
            call. = FALSE
        )
    }
    if (!(tree.name %in% rowTreeNames(x))) {
        stop("'", name, "' must specify a tree from 'rowTreeNames(x)'.",
            call. = FALSE
        )
    }
}

.check_colTree_present <- function(tree.name, x, name = .get_name_in_parent(tree.name)) {
    if (!.is_non_empty_string(tree.name)) {
        stop("'", name, "' must be a single non-empty character value.",
            call. = FALSE
        )
    }
    if (!(tree.name %in% colTreeNames(x))) {
        stop("'", name, "' must specify a tree from 'colTreeNames(x)'.",
            call. = FALSE
        )
    }
}

# Check if alternative experiment can be found from altExp slot.
.check_altExp_present <- function(altexp, tse, altExpName = .get_name_in_parent(altexp),
                                  tse_name = .get_name_in_parent(tse), .disable.altexp = FALSE, ...) {
    # Disable altExp if specified
    if (!.is_a_bool(.disable.altexp)) {
        stop("'.disable.altexp' must be TRUE or FALSE.", call. = FALSE)
    }
    if (.disable.altexp) {
        altexp <- NULL
    }
    # Check that altexp.name must be an integer or name
    if (!(.is_a_string(altexp) || .is_an_integer(altexp) || is.null(altexp))) {
        stop(
            "'", altExpName, "' must be a string or an integer.",
            call. = FALSE
        )
    }
    # If is not NULL, but the object does not have altExp slot
    if (!is.null(altexp) && !is(tse, "SingleCellExperiment")) {
        stop(
            "'", altExpName, "', is specified but '", tse_name, "' does not ",
            "have altExp slot.",
            call. = FALSE
        )
    }
    # Then check that altExp can be found; name or index.
    if (!is.null(altexp) && !altexp %in% c(
        altExpNames(tse), seq_len(length(altExps(tse)))
    )) {
        stop(
            "'", altExpName, "', does not specify an experiment from altExp ",
            "slot of '", tse_name, "'.",
            call. = FALSE
        )
    }
}

# Check if metadata has the specified data.
.check_metadata_present <- function(data.type, x, name = .get_name_in_parent(data.type)) {
    if (!.is_non_empty_string(data.type)) {
        stop("'", name, "' must be a single non-empty character value.",
            call. = FALSE
        )
    }
    if (!(data.type %in% names(metadata(x)))) {
        stop("'", name, "' must be a valid name of metadata(x)", call. = FALSE)
    }
    return(data.type)
}

################################################################################

# This function retrieves specific tables from MAE
#' @importFrom MultiAssayExperiment intersectColumns
.get_shared_samples_from_mae <- function(x, experiments, assay.types, ...) {
    # Select experiments from MAE
    x <- .select_experiments(x, experiments)
    # Select samples that are shared between experiments
    x <- intersectColumns(x)
    # Select alternative experiments
    x <- .select_altexps(x, ...)
    # Select assays of each experiment
    x <- .select_assays_from_mae(x, assay.types)
    # Get tables as a list
    mat_list <- MultiAssayExperiment::assays(x) |> as.list()
    return(mat_list)
}

# This function retrieves specific tables from TreeSE
.get_shared_samples_from_tse <- function(x, experiments, assay.types) {
    # Select experiments from TreeSE
    tse_list <- .select_experiments_from_tse(x, experiments)
    # Select assays of each experiment
    tse_list <- .select_assays_from_tse_list(tse_list, assay.types)
    # Get tables as a list
    mat_list <- lapply(tse_list, assay) |> as.list()
    return(mat_list)
}

# These following functions are for extracting more than one table from SCE or
# MAE.
# Select experiments from MAE
#' @importFrom MultiAssayExperiment experiments
.select_experiments <- function(mae, experiments) {
    # Check that the value is correct
    is_name <- is.character(experiments) && length(experiments) > 0 &&
        length(experiments) <= length(experiments(mae)) &&
        all(experiments %in% names(mae))
    is_index <- is.numeric(experiments) && all(experiments %% 1 == 0) &&
        length(experiments) > 0 &&
        length(experiments) <= length(experiments(mae)) &&
        all(experiments > 0 & experiments <= length(experiments(mae)))
    if (!(is_name || is_index)) {
        stop("'experiments' must specify names or index of experiments of ",
            "'x'",
            call. = FALSE
        )
    }
    # Subset experiments
    mae <- mae[, , experiments]
    # Check that all objects are SE
    all_SE <- lapply(experiments(mae), function(x) {
        is(x, "SummarizedExperiment")
    }) |>
        unlist() |>
        all()
    if (!all_SE) {
        stop("All experiments must be SummarizedExperiment objects.",
            call. = FALSE
        )
    }
    return(mae)
}

# Select optionally alternative experiments
#' @importFrom SingleCellExperiment altExps
.select_altexps <- function(mae, altexps = NULL, ...) {
    # Check that value is correct
    is_name <- is.character(altexps) && length(altexps) > 0 &&
        length(altexps) <= length(experiments(mae))
    is_index <- is.numeric(altexps) && all(altexps %% 1 == 0) &&
        length(altexps) > 0 &&
        length(altexps) <= length(experiments(mae))
    if (!(is.null(altexps) || is_name || is_index)) {
        stop("'altexps' must be NULL or specify alternative experiments for ",
            "each experiment.",
            call. = FALSE
        )
    }
    # If specified, select altExps from experiments
    if (!is.null(altexps)) {
        if (!require("SingleCellExperiment")) {
            stop("To enable 'altexps' option, 'SingleCellExperiment' package ",
                "must be installed.",
                call. = FALSE
            )
        }
        names(altexps) <- names(mae)
        for (exp in names(mae)) {
            # Get altExp if it is not NA, which disables altExp for single
            # experiment
            if (!is.na(altexps[[exp]])) {
                if (!is(mae[[exp]], "SingleCellExperiment")) {
                    stop("Experiment '", exp, "' must be SingleCellExperiment ",
                        "object.",
                        call. = FALSE
                    )
                }
                # Check that alExp can be found
                is_name <- is.character(altexps[[exp]]) &&
                    altexps[[exp]] %in% altExpNames(mae[[exp]])
                is_index <- is.numeric(altexps[[exp]]) &&
                    all(altexps[[exp]] %% 1 == 0) &&
                    altexps[[exp]] > 0 &&
                    altexps[[exp]] <= length(altExps(mae[[exp]]))
                if (!(is_name || is_index)) {
                    stop("'", altexps[[exp]], "' does not specify altExp from ",
                        "experiment '", exp, "'.",
                        call. = FALSE
                    )
                }
                mae[[exp]] <- altExp(mae[[exp]], altexps[[exp]])
            }
        }
    }
    return(mae)
}

# Select assays to be included in MOFA
#' @importFrom MultiAssayExperiment experiments
#' @importFrom SummarizedExperiment assay assays assayNames
.select_assays_from_mae <- function(mae, assay.types) {
    # Check that value is correct
    is_name <- is.character(assay.types) && length(assay.types) > 0 &&
        length(assay.types) <= length(experiments(mae))
    if (!is_name) {
        stop("'assay.type' must specify name of assays. The lenght must equal ",
            "to 'experiments'.",
            call. = FALSE
        )
    }
    # Give corresponding experiment names to assay.types
    names(assay.types) <- names(mae)
    # For every experiment in MAE
    for (exp in names(mae)) {
        # Check that assay exists
        if (!assay.types[[exp]] %in% assayNames(mae[[exp]])) {
            stop("Cannot find assay '", assay.types[[exp]], "' from ",
                "experiment '", exp, "'.",
                call. = FALSE
            )
        }
        # Keep only selected assay.type from a given experiment
        assays(mae[[exp]]) <- assays(mae[[exp]])[assay.types[[exp]]]
    }
    return(mae)
}

#' @importFrom SingleCellExperiment altExps
.select_experiments_from_tse <- function(tse, experiments) {
    tse_list <- altExps(tse)
    main_name <- "main"

    all_names <- make.unique(c(names(tse_list), main_name))
    main_name <- all_names[[length(all_names)]]
    tse_list[[main_name]] <- tse

    if (!is.character(experiments)) {
        main_name <- length(tse_list)
    }
    experiments[is.na(experiments)] <- main_name

    # Check that the value is correct
    is_name <- is.character(experiments) && length(experiments) > 0 &&
        length(experiments) <= length(tse_list) &&
        all(experiments %in% names(tse_list))
    is_index <- is.numeric(experiments) && all(experiments %% 1 == 0) &&
        length(experiments) > 0 &&
        length(experiments) <= length(tse_list) &&
        all(experiments > 0 & experiments <= length(tse_list))
    if (!(is_name || is_index)) {
        stop("'experiments' must specify names or index of alternative ",
            "experiments of 'x'",
            call. = FALSE
        )
    }

    tse_list <- tse_list[experiments]
    return(tse_list)
}

.select_assays_from_tse_list <- function(tse_list, assay.types) {
    # Check that value is correct
    is_name <- is.character(assay.types) && length(assay.types) > 0 &&
        length(assay.types) <= length(tse_list)
    if (!is_name) {
        stop("'assay.type' must specify name of assays. The lenght must equal ",
            "to 'experiments'.",
            call. = FALSE
        )
    }
    # Give corresponding experiment names to assay.types
    names(assay.types) <- names(tse_list)
    # For every experiment in MAE
    for (exp in names(tse_list)) {
        # Check that assay exists
        if (!assay.types[[exp]] %in% assayNames(tse_list[[exp]])) {
            stop("Cannot find assay '", assay.types[[exp]], "' from ",
                "experiment '", exp, "'.",
                call. = FALSE
            )
        }
        # Keep only selected assay.type from a given experiment
        assays(tse_list[[exp]]) <- assays(tse_list[[exp]])[assay.types[[exp]]]
    }
    return(tse_list)
}

.retrieve_sample_metadata_from_mae <- function(x, experiments, ...) {
    # Select experiments from MAE
    x <- .select_experiments(x, experiments)
    # Select samples that are shared between experiments
    x <- intersectColumns(x)
    # Select alternative experiments
    x <- .select_altexps(x, ...)
    # Get all colDatas
    dfs <- lapply(x |> length() |> seq_len(), function(i) {
        getWithColData(x, i) |> colData()
    })
    df <- do.call(cbind, dfs)
    df <- df[, !duplicated(colnames(df)), drop = FALSE]
    return(df)
}

.retrieve_sample_metadata_from_tse <- function(x, experiments) {
    # Select experiments from TreeSE
    tse_list <- .select_experiments_from_tse(x, experiments)
    # Get all colDatas
    dfs <- lapply(tse_list, colData)
    df <- do.call(cbind, dfs)
    df <- df[, !duplicated(colnames(df)), drop = FALSE]
    return(df)
}
