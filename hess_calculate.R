#' Enthalpie totale a partir de coefficients et d'enthalpies
#'
#' Calcule la somme ponderee `sum(coefficients * enthalpies)`.
#'
#' @param coefficients Vecteur numerique des coefficients de combinaison.
#' @param enthalpies Vecteur numerique des enthalpies (kJ/mol), meme longueur.
#' @return Un nombre : l'enthalpie totale.
#' @examples
#' hess_calculate(c(1, 2), c(-100, -50))
#' @export
hess_calculate <- function(coefficients, enthalpies) {

  if (!is.numeric(coefficients) || !is.numeric(enthalpies)) {
    stop("Les coefficients et les enthalpies doivent etre numeriques.")
  }

  if (length(coefficients) != length(enthalpies)) {
    stop("Les deux vecteurs doivent avoir la meme longueur.")
  }

  if (length(coefficients) == 0 ||
      any(!is.finite(coefficients)) ||
      any(!is.finite(enthalpies))) {
    stop("Les vecteurs doivent contenir des valeurs finies.")
  }

  sum(coefficients * enthalpies)
}
