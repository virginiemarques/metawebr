#' Impute missing traits with phylogenetic random forests
#'
#' Fills missing trait values with random-forest imputation (missForest) using
#' the position of each species in the phylogeny as extra predictors, and
#' averages over several trees and runs to account for phylogenetic and
#' imputation uncertainty.
#'
#' For each tree:
#' 1. the tree is pruned to the species of `data_traits`;
#' 2. the first `n_eigen` phylogenetic eigenvectors are computed: principal
#'    coordinates of the patristic distance matrix, as in `PVR::PVRdecomp()`;
#' 3. [missForest::missForest()] is run `n_runs` times on the traits plus the
#'    eigenvectors.
#'
#' Numeric traits are then averaged over all runs and trees; other traits take
#' their most frequent imputed value. Observed values are never changed.
#'
#' Only species in the phylogeny can be imputed. The others keep their missing
#' values, with a warning. By default the trees are the 100 complete fish trees
#' of Rabosky et al. (2018), downloaded with
#' `fishtree::fishtree_complete_phylogeny()`. They cover ray-finned fishes only,
#' so pass your own `trees` to include sharks and rays.
#'
#' Runtime grows with the number of trees x `n_runs`: with the defaults, about
#' 5 minutes on one core for 150 species. Call [set.seed()] first for
#' reproducible results, whatever `mc.cores`.
#'
#' @param data_traits Traits with species names (`Genus_species`) as row names.
#' @param traits Trait columns used as predictors and imputed. Numeric columns
#'   are imputed as numbers, the others as categories.
#' @param trees `NULL` (default) to download the fishtree complete trees, or a
#'   `phylo` or `multiPhylo` object whose tip labels are species names (with
#'   spaces or underscores).
#' @param n_eigen Number of phylogenetic eigenvectors used as predictors.
#' @param n_runs Number of missForest runs per tree.
#' @param maxiter,ntree Passed to [missForest::missForest()].
#' @param mc.cores Number of cores; trees are processed in parallel (forked
#'   processes, not available on Windows).
#'
#' @return `data_traits` with the missing values of `traits` filled for species
#'   in the phylogeny, and a column `imputed` listing the traits imputed for
#'   each taxon (`NA` when none).
#'
#' @references
#' Diniz-Filho, J.A.F., de Sant'Ana, C.E.R. & Bini, L.M. (1998). An eigenvector
#' method for estimating phylogenetic inertia. *Evolution*, 52, 1247-1262.
#'
#' Rabosky, D.L. et al. (2018). An inverse latitudinal gradient in speciation
#' rate for marine fishes. *Nature*, 559, 392-395.
#'
#' Stekhoven, D.J. & Buhlmann, P. (2012). MissForest: non-parametric missing
#' value imputation for mixed-type data. *Bioinformatics*, 28, 112-118.
#'
#' @seealso [infer_traits()] with `method = "100matrix"`, [clean_traits()] for a
#'   faster genus-based filling.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' traits <- impute_traits_phylo(df_traits_obis, mc.cores = 4)
#' traits[!is.na(traits$imputed), ]
#' }
#'
#' @export
impute_traits_phylo <- function(data_traits,
                                traits = c("CommonLengthEstim", "TrophicLevel", "DemersPelag"),
                                trees = NULL,
                                n_eigen = 20,
                                n_runs = 20,
                                maxiter = 20,
                                ntree = 300,
                                mc.cores = 1) {
  for (pkg in c("ape", "missForest")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is needed for this function: install.packages(\"%s\").", pkg, pkg),
           call. = FALSE)
    }
  }
  check_traits(data_traits, traits)
  for (arg in c("n_eigen", "n_runs", "maxiter", "ntree")) check_count(get(arg), arg)

  data_traits$imputed <- NA_character_
  missing <- is.na(data_traits[, traits, drop = FALSE])
  if (!any(missing)) {
    message("No missing trait values to impute.")
    return(data_traits)
  }

  species <- rownames(data_traits)
  trees <- prepare_trees(trees, species)
  tips <- Reduce(intersect, lapply(trees, function(tree) tree$tip.label))
  tree_species <- species[species %in% tips]

  not_in_tree <- species[rowSums(missing) > 0 & !species %in% tree_species]
  if (length(not_in_tree) > 0) {
    warning(sprintf("%d taxa with missing traits are not in the phylogeny and were not imputed: %s.",
                    length(not_in_tree), format_taxa(not_in_tree)), call. = FALSE)
  }
  if (length(tree_species) < 3) {
    stop("At least 3 species of `data_traits` must be in the phylogeny.", call. = FALSE)
  }
  if (!any(missing[tree_species, ])) {
    message("No missing trait values for species in the phylogeny.")
    return(data_traits)
  }

  x <- data_traits[tree_species, traits, drop = FALSE]
  categorical <- !vapply(x, is.numeric, logical(1))
  x[categorical] <- lapply(x[categorical], factor)

  # One seed per tree and run, drawn in the main process
  seeds <- matrix(sample.int(.Machine$integer.max, length(trees) * n_runs), nrow = length(trees))
  old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", old_seed, envir = globalenv()))

  runs <- run_parallel(seq_along(trees), function(k) {
    tree <- ape::keep.tip(trees[[k]], tree_species)
    predictors <- data.frame(x, phylo_eigenvectors(tree, n_eigen)[tree_species, , drop = FALSE])
    lapply(seq_len(n_runs), function(r) {
      set.seed(seeds[k, r])
      # randomForest warns about regressions on few unique values; harmless here
      imp <- suppressWarnings(missForest::missForest(predictors, maxiter = maxiter, ntree = ntree))
      imp$ximp[, traits, drop = FALSE]
    })
  }, mc.cores)
  runs <- unlist(runs, recursive = FALSE)

  # Average numeric traits, take the most frequent category for the others
  n_imputed <- 0
  for (trait in traits) {
    rows <- which(is.na(x[[trait]]))
    if (length(rows) == 0) next
    values <- vapply(runs, function(run) as.character(run[rows, trait]), character(length(rows)))
    values <- matrix(values, nrow = length(rows))
    filled <- if (categorical[[trait]]) {
      apply(values, 1, mode_or_na)
    } else {
      rowMeans(matrix(as.numeric(values), nrow = length(rows)))
    }

    taxa <- tree_species[rows]
    if (is.factor(data_traits[[trait]])) data_traits[[trait]] <- as.character(data_traits[[trait]])
    data_traits[taxa, trait] <- filled
    previous <- data_traits[taxa, "imputed"]
    data_traits[taxa, "imputed"] <- ifelse(is.na(previous), trait, paste(previous, trait, sep = ", "))
    n_imputed <- n_imputed + length(rows)
  }

  message(sprintf("Imputed %d missing values for %d taxa (%d trees x %d runs).",
                  n_imputed, sum(!is.na(data_traits$imputed)), length(trees), n_runs))
  data_traits
}

# List of phylo objects with underscores in tip labels
prepare_trees <- function(trees, species) {
  if (is.null(trees)) {
    if (!requireNamespace("fishtree", quietly = TRUE)) {
      stop("Package 'fishtree' is needed to download trees: install.packages(\"fishtree\"), or pass `trees`.",
           call. = FALSE)
    }
    binomials <- gsub("_", " ", species[grepl("^[A-Z][a-z]+_[a-z]", species)])
    trees <- withCallingHandlers(
      fishtree::fishtree_complete_phylogeny(species = binomials),
      # species absent from the trees are reported later
      missing_species = function(w) invokeRestart("muffleWarning")
    )
  }
  if (inherits(trees, "phylo")) trees <- list(trees)
  is_tree_list <- inherits(trees, "multiPhylo") ||
    (is.list(trees) && length(trees) > 0 && all(vapply(trees, inherits, logical(1), what = "phylo")))
  if (!is_tree_list) {
    stop("`trees` must be a phylo or multiPhylo object.", call. = FALSE)
  }
  lapply(seq_along(trees), function(i) {
    tree <- trees[[i]]
    tree$tip.label <- gsub(" ", "_", tree$tip.label)
    tree
  })
}

# First phylogenetic eigenvectors (principal coordinates of the patristic
# distances), computed as in PVR::PVRdecomp(); rows named by species
phylo_eigenvectors <- function(tree, n_eigen = 20) {
  D <- ape::cophenetic.phylo(tree)
  n <- nrow(D)
  J <- diag(n) - 1 / n
  G <- J %*% (-0.5 * D) %*% J
  e <- eigen(G, symmetric = TRUE)
  keep <- seq_len(min(n_eigen, sum(e$values > 0)))
  vectors <- e$vectors[, keep, drop = FALSE]
  dimnames(vectors) <- list(rownames(D), paste0("phylo_eigen_", keep))
  vectors
}
