#' Multiplier une reaction par un facteur
#'
#' Multiplie les coefficients stoechiometriques et l'enthalpie par un facteur
#' strictement positif (pour un facteur negatif, utiliser [hess_reverse()]).
#'
#' @param reaction Liste avec `reactants`, `products` et `delta_h`.
#' @param factor Nombre strictement positif.
#' @return La reaction multipliee (liste).
#' @examples
#' r <- list(reactants = c(A = 1), products = c(B = 1), delta_h = -50)
#' hess_multiply(r, 2)
#' @export
hess_multiply <- function(reaction, factor) {

  if (!is.list(reaction) ||
      is.null(reaction$reactants) ||
      is.null(reaction$products) ||
      is.null(reaction$delta_h)) {
    stop("Format de reaction invalide.")
  }

  if (!is.numeric(factor) ||
      length(factor) != 1 ||
      !is.finite(factor) ||
      factor <= 0) {
    stop("Le facteur doit etre un nombre positif.")
  }

  if (!is.numeric(reaction$delta_h) ||
      length(reaction$delta_h) != 1 ||
      !is.finite(reaction$delta_h)) {
    stop("delta_h doit etre un nombre fini.")
  }

  list(
    reactants = reaction$reactants * factor,
    products = reaction$products * factor,
    delta_h = reaction$delta_h * factor
  )
}
