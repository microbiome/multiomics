# All generic methods are listed here

#' @rdname getMantel
#' @export
setGeneric("getMantel", signature = "x", function(x, ...) {
    standardGeneric("getMantel")
})

#' @rdname getCrossAssociation
#' @export
setGeneric(
    "getCrossAssociation",
    signature = "x",
    function(x, method = "pearson", ...) {
        standardGeneric("getCrossAssociation")
    }
)

#' @rdname getRDA
#' @export
setGeneric("getRDA", signature = "x", function(x, ...) {
    standardGeneric("getRDA")
})

#' @rdname getMiRKAT
#' @export
setGeneric("getMiRKAT", signature = "x", function(x, ...) {
    standardGeneric("getMiRKAT")
})

#' @rdname getLASSO
#' @export
setGeneric("getLASSO", signature = "x", function(x, ...) {
    standardGeneric("getLASSO")
})
