#' Verifier l'equilibrage des atomes d'une reaction
#'
#' Analyse les formules chimiques (noms des especes, ex. `"CH4(g)"`,
#' `"Ca(OH)2(s)"`) et compare le nombre d'atomes de chaque element a gauche
#' et a droite. Si une formule n'est pas reconnue, le resultat est `NA`.
#'
#' @param reaction Liste avec `reactants` et `products`.
#' @param tolerance Tolerance numerique.
#' @return Liste : `balanced` (TRUE/FALSE/NA), `difference` (atomes produits
#'   moins atomes reactifs, par element) et `message`.
#' @examples
#' r <- list(reactants = c("CH4(g)" = 1, "O2(g)" = 2),
#'           products = c("CO2(g)" = 1, "H2O(l)" = 2), delta_h = -890.3)
#' hess_check_balance(r)$balanced
#' @export
hess_check_balance <- function(reaction, tolerance = 1e-8) {

  if (!is.list(reaction) ||
      is.null(reaction$reactants) ||
      is.null(reaction$products)) {
    stop("La reaction doit contenir reactants et products.")
  }

  side_atoms <- function(x) {
    total <- numeric(0)
    for (sp in names(x)) {
      atoms <- .hess_parse_formula(sp)
      for (el in names(atoms)) {
        current <- if (el %in% names(total)) total[[el]] else 0
        total[el] <- current + x[[sp]] * atoms[[el]]
      }
    }
    total
  }

  get_or_zero <- function(vec, el) {
    if (el %in% names(vec)) vec[[el]] else 0
  }

  tryCatch({
    left <- side_atoms(reaction$reactants)
    right <- side_atoms(reaction$products)
    elements <- union(names(left), names(right))

    diff <- vapply(elements, function(el) {
      get_or_zero(right, el) - get_or_zero(left, el)
    }, numeric(1))

    balanced <- all(abs(diff) <= tolerance)

    if (balanced) {
      msg <- "Equation equilibree (atomes conserves)"
    } else {
      bad <- diff[abs(diff) > tolerance]
      msg <- paste0(
        "DESEQUILIBRE (produits - reactifs) : ",
        paste0(names(bad), " ", format(bad, digits = 4, trim = TRUE),
               collapse = ", ")
      )
    }

    list(balanced = balanced, difference = diff, message = msg)
  }, error = function(e) {
    list(
      balanced = NA,
      difference = NULL,
      message = paste0("Non verifiable (", conditionMessage(e), ")")
    )
  })
}
