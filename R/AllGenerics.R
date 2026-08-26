# All generic methods are listed here

#' @rdname getMantel
#' @export
setGeneric("getMantel", signature = "x", function(x, ...) {
    standardGeneric("getMantel")
})

#' @rdname getJointRDA
#' @export
setGeneric("getJointRDA", signature = "x", function(x, ...) {
    standardGeneric("getJointRDA")
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

#' @rdname getPairwiseAssociation
#' @export
setGeneric("getPairwiseAssociation", signature = "x", function(x, ...) {
    standardGeneric("getPairwiseAssociation")
})
