create_pfm <- function(s, a, dt = 1, bins = 4) {
  Nt <- length(s)
  V <- data.frame(
    s = s[(dt + 1):Nt],
    as.data.frame(lapply(a, function(x) x[1:(Nt - dt)]))
  )
  Vb <- as.data.frame(lapply(V, function(x) {
    cut(x, breaks = bins, labels = FALSE, include.lowest = TRUE)
  }))
  f <- xtabs(~ ., data = Vb)
  p <- prop.table(f)
  array(p, dim = rep(bins, ncol(Vb)))
}

set.seed(10)
N <- 1e7   # ⚠️ R不如NumPy快，建议先试1e6
dt <- 1
bins <- 2

# Generate random binary variables
q1 <- round(runif(N))
q2 <- round(runif(N))

# Example 1: duplicated input (same agent twice)
target <- c(tail(q1, -dt), rep(NA, dt))  # shift like np.roll(q1, dt)
target <- ifelse(is.na(target), q1[1:dt], target)  # wrap around
agents <- list(q1, q1)

# Compute histogram / PMF
hist <- create_pfm(target, agents, dt, bins)