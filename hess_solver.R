#' Trouver les coefficients de combinaison
#'
#' Resout `reactions %*% x = target` par pseudo-inverse (SVD). Fonctionne pour
#' les systemes carres, sur-determines et sous-determines (solution de norme
#' minimale). Le champ `exact` indique si la cible est reproduite.
#'
#' @param reactions Matrice numerique (une colonne par reaction).
#' @param target Vecteur cible (une valeur par ligne de la matrice).
#' @param tolerance Tolerance sur le residu (nombre positif).
#' @return Liste : `coefficients`, `residual`, `exact`, `rank`, `unique`.
#' @examples
#' hess_solver(matrix(c(1, 0, 0, 1), nrow = 2), c(2, 3))$coefficients
#' @export
hess_solver <- function(reactions, target, tolerance = 1e-8) {

  if (!is.matrix(reactions) || !is.numeric(reactions)) {
    stop("reactions doit etre une matrice numerique.")
  }

  if (!is.numeric(target) ||
      length(target) != nrow(reactions)) {
    stop("target doit avoir une valeur par ligne de la matrice.")
  }

  if (any(!is.finite(reactions)) ||
      any(!is.finite(target))) {
    stop("La matrice et la cible doivent contenir des valeurs finies.")
  }

  if (!is.numeric(tolerance) ||
      length(tolerance) != 1 ||
      !is.finite(tolerance) ||
      tolerance <= 0) {
    stop("tolerance doit etre un nombre positif.")
  }

  sv <- svd(reactions)
  rank_tol <- max(dim(reactions)) * .Machine$double.eps * max(sv$d, 0)
  keep <- sv$d > rank_tol
  rank <- sum(keep)

  if (rank == 0) {
    solution <- numeric(ncol(reactions))
  } else {
    u <- sv$u[, keep, drop = FALSE]
    v <- sv$v[, keep, drop = FALSE]
    solution <- as.vector(v %*% ((t(u) %*% as.numeric(target)) / sv$d[keep]))
  }

  solution[abs(solution) < 1e-10] <- 0
  solution <- round(solution, 10)

  residual <- as.vector(reactions %*% solution - as.numeric(target))
  exact <- sqrt(sum(residual^2)) <= tolerance

  if (!exact) {
    warning("La cible n'est pas reproduite exactement.")
  }

  names(solution) <- colnames(reactions)

  list(
    coefficients = solution,
    residual = residual,
    exact = exact,
    rank = rank,
    unique = (rank == ncol(reactions))
  )
}
