#' Analyse complete avec la loi de Hess (point d'entree du package)
#'
#' Une seule fonction pour tout faire. A partir des reactions (et
#' eventuellement d'une reaction cible), `hess_run()` active automatiquement
#' toutes les fonctions du package :
#'
#' 1. [hess_matrix()] et [hess_solver()] : recherche des coefficients si une
#'    cible est donnee et que `coefficients` est `NULL` ;
#' 2. [hess_reverse()] et [hess_multiply()] : mise en forme de chaque
#'    reaction (inversion si coefficient negatif, multiplication) ;
#' 3. [hess_combine()] : equation bilan ;
#' 4. [hess_calculate()] : enthalpie totale (verifiee contre la combinaison) ;
#' 5. [hess_check_balance()] : equilibrage des atomes ;
#' 6. [hess_classify()] : exothermique / endothermique / thermoneutre ;
#' 7. [hess_summary()] : resume ;
#' 8. [hess_diagram()] : diagramme energetique colore.
#'
#' @param reactions Liste de reactions (`reactants`, `products`, `delta_h`).
#'   Une seule reaction peut etre passee directement. Les noms de la liste
#'   servent d'etiquettes aux etapes.
#' @param coefficients Coefficients de combinaison (optionnel). Si `NULL` et
#'   `target` fourni : calcules automatiquement ; sinon 1 pour chacune.
#' @param target Reaction cible (`reactants` et `products`, `delta_h`
#'   optionnel comme valeur de reference).
#' @param reactants_energy Niveau d'energie des reactifs pour le diagramme.
#' @param show_diagram Afficher le diagramme (`TRUE` par defaut).
#' @param check_balance Verifier l'equilibrage des atomes (`TRUE` par defaut).
#' @param verbose Afficher le resume dans la console (`TRUE` par defaut).
#' @return De facon invisible, un objet `hess_run` (liste) avec `reactants`,
#'   `products`, `delta_h`, `classification`, `balance`, `coefficients`,
#'   `steps`, `equation`, `atom_check`, `target_reached` et `diagram`.
#' @examples
#' # Combustion du methane
#' ex <- hess_example("methane")
#' res <- hess_run(ex$reactions, show_diagram = FALSE)
#'
#' # Les coefficients sont trouves automatiquement a partir de la cible
#' ex <- hess_example("formation_methane")
#' res <- hess_run(ex$reactions, target = ex$target, show_diagram = FALSE)
#' res$delta_h
#'
#' \donttest{
#' hess_run(ex$reactions, target = ex$target)
#' }
#' @export
hess_run <- function(reactions, coefficients = NULL, target = NULL,
                     reactants_energy = 0, show_diagram = TRUE,
                     check_balance = TRUE, verbose = TRUE) {

  .hess_flag(show_diagram, "show_diagram")
  .hess_flag(check_balance, "check_balance")
  .hess_flag(verbose, "verbose")

  if (.hess_is_reaction(reactions)) {
    reactions <- list(reactions)
  }

  if (!is.list(reactions) || length(reactions) == 0) {
    stop("Fournissez une liste de reactions.", call. = FALSE)
  }

  n <- length(reactions)

  for (i in seq_len(n)) {
    .hess_validate_reaction(reactions[[i]], paste("Reaction", i))
  }

  step_names <- names(reactions)
  if (is.null(step_names)) step_names <- rep("", n)
  empty <- is.na(step_names) | step_names == ""
  step_names[empty] <- paste("Etape", which(empty))

  # 1. Coefficients : fournis ou calcules par le solveur -----------------
  m <- NULL
  if (!is.null(target)) {
    .hess_validate_reaction(target, "Cible", need_h = FALSE)
    m <- hess_matrix(reactions, target)
  }

  if (is.null(coefficients)) {
    if (is.null(target)) {
      coefficients <- rep(1, n)
    } else {
      sol <- suppressWarnings(hess_solver(m$matrix, m$target))
      if (!sol$exact) {
        stop("La reaction cible ne peut pas etre obtenue par combinaison ",
             "des reactions fournies.", call. = FALSE)
      }
      coefficients <- as.numeric(sol$coefficients)
    }
  } else {
    if (!is.numeric(coefficients) || length(coefficients) != n ||
        any(!is.finite(coefficients))) {
      stop("coefficients doit etre numerique, fini et de longueur ", n, ".",
           call. = FALSE)
    }
    coefficients <- as.numeric(coefficients)
  }

  target_reached <- NA
  if (!is.null(target)) {
    residual <- as.vector(m$matrix %*% coefficients - m$target)
    target_reached <- sqrt(sum(residual^2)) <= 1e-8
    if (!target_reached) {
      warning("Les coefficients fournis ne reproduisent pas exactement ",
              "la reaction cible.", call. = FALSE)
    }
  }

  names(coefficients) <- step_names

  # 2. Inversion / multiplication de chaque reaction ----------------------
  adjusted <- list()
  for (i in seq_len(n)) {
    ci <- coefficients[[i]]
    if (abs(ci) < 1e-12) next
    r <- reactions[[i]]
    if (ci < 0) r <- hess_reverse(r)
    adjusted[[length(adjusted) + 1L]] <- hess_multiply(r, abs(ci))
  }

  if (length(adjusted) == 0) {
    stop("Tous les coefficients sont nuls : aucune reaction a combiner.",
         call. = FALSE)
  }

  # 3. Equation bilan -----------------------------------------------------
  combined <- hess_combine(adjusted)

  # 4. Enthalpie totale (et verification croisee) -------------------------
  enthalpies <- unname(vapply(reactions, function(r) as.numeric(r$delta_h),
                              numeric(1)))
  delta_h <- hess_calculate(as.numeric(coefficients), enthalpies)

  if (abs(delta_h - combined$delta_h) > 1e-6 * max(1, abs(delta_h))) {
    stop("Incoherence interne entre hess_calculate() et hess_combine().",
         call. = FALSE)
  }

  steps <- coefficients * enthalpies
  steps <- steps[abs(coefficients) >= 1e-12]

  # 5. Atomes, 6. classification -----------------------------------------
  atom_check <- if (check_balance) hess_check_balance(combined) else NULL
  classification <- hess_classify(delta_h)

  output <- list(
    reactants = combined$reactants,
    products = combined$products,
    delta_h = delta_h,
    classification = classification,
    balance = combined$balance,
    coefficients = coefficients,
    steps = steps,
    equation = .hess_equation(combined$reactants, combined$products),
    atom_check = atom_check,
    target_reached = target_reached
  )

  if (!is.null(target) && !is.null(target$delta_h) &&
      is.numeric(target$delta_h) && is.finite(target$delta_h)) {
    output$reference_delta_h <- target$delta_h
  }

  class(output) <- "hess_run"

  # 7. Resume, 8. diagramme ----------------------------------------------
  if (verbose) {
    hess_summary(output)
  }

  if (show_diagram) {
    p <- hess_diagram(
      delta_h = delta_h,
      reactants_energy = reactants_energy,
      steps = if (length(steps) >= 2) steps else NULL
    )
    output$diagram <- p
    print(p)
  }

  invisible(output)
}

#' @export
print.hess_run <- function(x, ...) {
  hess_summary(x)
  invisible(x)
}
