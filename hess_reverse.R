#' Inverser une reaction
#'
#' Echange reactifs et produits et change le signe de l'enthalpie.
#'
#' @param reaction Liste avec `reactants`, `products` (vecteurs nommes) et
#'   `delta_h` (nombre).
#' @return La reaction inversee (liste).
#' @examples
#' r <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -50)
#' hess_reverse(r)
#' @export
hess_reverse <- function(reaction) {

  if (!is.list(reaction) ||
      is.null(reaction$reactants) ||
      is.null(reaction$products) ||
      is.null(reaction$delta_h)) {
    stop("La reaction doit contenir reactants, products et delta_h.")
  }

  if (!is.numeric(reaction$delta_h) ||
      length(reaction$delta_h) != 1 ||
      !is.finite(reaction$delta_h)) {
    stop("delta_h doit etre un nombre fini.")
  }

  list(
    reactants = reaction$products,
    products = reaction$reactants,
    delta_h = -reaction$delta_h
  )
}
