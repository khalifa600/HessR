# HessR

Package R pour appliquer la **loi de Hess** aux réactions chimiques (combustion, formation, etc.).

**Une seule fonction suffit : `hess_run()`.** Dès que l'ingénieur fournit ses réactions, toutes les fonctions s'activent ensemble : résolution des coefficients, inversion, multiplication, combinaison, calcul de ΔH, vérification des atomes, classification, résumé et diagramme coloré.

## Installation (RStudio)

1. Dézipper le dossier `HessR`, puis double-cliquer sur `HessR.Rproj`.
2. Dans la console :

```r
install.packages(c("devtools", "testthat", "ggplot2", "rlang"))
devtools::document()   # génère l'aide (man/)
devtools::test()       # lance les tests
devtools::install()    # installe le package
library(HessR)
```

## Utilisation rapide

### 1. Une réaction (combustion du méthane)

```r
library(HessR)

methane <- list(
  reactants = c("CH4(g)" = 1, "O2(g)" = 2),
  products  = c("CO2(g)" = 1, "H2O(l)" = 2),
  delta_h   = -890.3            # kJ/mol
)

hess_run(methane)
```

### 2. Plusieurs réactions + réaction cible (coefficients trouvés automatiquement)

```r
ex <- hess_example("formation_methane")   # C + 2 H2 -> CH4

res <- hess_run(ex$reactions, target = ex$target)
res$delta_h        # -74.8 kJ/mol
res$coefficients   # 1, 2, -1
```

### 3. Coefficients choisis à la main

```r
hess_run(ex$reactions, coefficients = c(1, 2, -1))
```

## Ce que fait `hess_run()`

| Étape | Fonction |
|---|---|
| Matrice stœchiométrique | `hess_matrix()` |
| Coefficients (si cible donnée) | `hess_solver()` |
| Inversion (coefficient < 0) | `hess_reverse()` |
| Multiplication | `hess_multiply()` |
| Équation bilan | `hess_combine()` |
| Enthalpie totale | `hess_calculate()` |
| Vérification des atomes | `hess_check_balance()` |
| Classification | `hess_classify()` |
| Résumé console | `hess_summary()` |
| Diagramme coloré | `hess_diagram()` |

Options utiles : `show_diagram = FALSE`, `verbose = FALSE`, `check_balance = FALSE`, `reactants_energy = 0`.

Le résultat contient `delta_h`, `classification`, `equation`, `coefficients`, `steps`, `atom_check`, `target_reached` et `diagram` (objet `ggplot`, modifiable ou sauvegardable avec `ggplot2::ggsave()`).

## Format d'une réaction

```r
list(
  reactants = c("A" = 1, "B" = 2),   # vecteur nommé, coefficients > 0
  products  = c("C" = 1),
  delta_h   = -100                    # nombre fini, kJ/mol
)
```

Écrire les formules avec leur état physique (`"H2O(l)"`, `"CO2(g)"`, `"C(s)"`) permet la vérification automatique des atomes. Les noms de la liste de réactions servent d'étiquettes dans le résumé et le diagramme.

## Unités

Les enthalpies doivent être exprimées dans une unité cohérente (kJ/mol de réaction telle qu'écrite).

## Limites

HessR calcule des variations d'enthalpie. Il ne calcule pas le rendement réel d'un moteur ou d'une installation (pertes thermiques, conditions de fonctionnement). Les états physiques doivent être cohérents avec les données utilisées.
