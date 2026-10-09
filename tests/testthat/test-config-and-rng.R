test_that("config merge replaces YAML sequences and preserves null", {
  base <- list(
    subsets = list(list(name = "template")),
    singler = list(labels = list(list(name = "main"))),
    groupings_main = "clusters"
  )
  override <- list(
    subsets = list(),
    singler = list(labels = list()),
    groupings_main = NULL
  )

  merged <- SPARKS:::.merge_config(base, override)

  expect_length(merged$subsets, 0L)
  expect_length(merged$singler$labels, 0L)
  expect_null(merged$groupings_main)
})

test_that("seed helper restores the caller RNG state", {
  set.seed(8137)
  expected_seed <- get(".Random.seed", envir = .GlobalEnv)
  restore_rng <- SPARKS:::.set_seed_preserving_state(2861L)
  seeded_draw <- runif(1)
  restore_rng()

  expect_identical(get(".Random.seed", envir = .GlobalEnv), expected_seed)
  expect_type(seeded_draw, "double")
})

test_that("seed helper produces the same stream for the same seed", {
  draw_once <- function() {
    restore_rng <- SPARKS:::.set_seed_preserving_state(73L)
    on.exit(restore_rng(), add = TRUE)
    runif(3)
  }

  expect_identical(draw_once(), draw_once())
})

test_that("seed helper restores an uninitialized RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)
  on.exit(if (had_seed) set.seed(8137), add = TRUE)

  restore_rng <- SPARKS:::.set_seed_preserving_state(2861L)
  runif(1)
  restore_rng()

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("seed helper rejects invalid values", {
  expect_error(SPARKS:::.validate_seed(NA_integer_), "seed must be")
  expect_error(SPARKS:::.validate_seed(-1L), "seed must be")
  expect_error(SPARKS:::.validate_seed(1.5), "seed must be")
  expect_error(SPARKS:::.validate_seed(c(1L, 2L)), "seed must be")
})
