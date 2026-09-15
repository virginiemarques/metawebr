# Trophic level and omnivory index of every node of a food web.
#
# `web` is an adjacency matrix with prey in rows and predators in columns.
# Gives the same result as `NetIndices::TrophInd(Tij = t(web))` (to ~1e-10) but
# solves the linear system with solve() instead of a generalised inverse, which
# is several times faster on large webs. As in NetIndices, self-feeding counts
# in the diet but does not raise the trophic level, and the generalised inverse
# is used as a fallback when the system is singular.
trophic_levels <- function(web) {
  Tij <- t(web)                       # rows = consumers, columns = resources
  p <- Tij / rowSums(Tij)             # diet fractions
  p[!is.finite(p)] <- 0
  A <- -p
  diag(A) <- 1
  n <- nrow(A)
  TL <- tryCatch(
    solve(A, rep(1, n)),
    error = function(e) drop(MASS::ginv(A) %*% rep(1, n))
  )
  OI <- rowSums(p * (outer(-TL, TL, "+") + 1)^2)
  data.frame(TL = unname(TL), OI = unname(OI), row.names = colnames(web))
}
