set.seed(10)

mediator <- function(N) {
  q1 <- numeric(N)
  q2 <- numeric(N)
  q3 <- numeric(N)
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)
  W3 <- rnorm(N, 0, 1)

  for (n in seq_len(N - 1)) {
    q1[n + 1] <- sin(q2[n]) + 0.001 * W1[n]
    q2[n + 1] <- cos(q3[n]) + 0.01 * W2[n]
    q3[n + 1] <- 0.5 * q3[n] + 0.1 * W3[n]
  }

  data.frame(q1 = q1, q2 = q2, q3 = q3)
}


confounder <- function(N) {
  q1 <- numeric(N)
  q2 <- numeric(N)
  q3 <- numeric(N)
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)
  W3 <- rnorm(N, 0, 1)

  for (n in seq_len(N - 1)) {
    q1[n + 1] <- sin(q1[n] + q3[n]) + 0.01 * W1[n]
    q2[n + 1] <- cos(q2[n] - q3[n]) + 0.01 * W2[n]
    q3[n + 1] <- 0.5 * q3[n] + 0.1 * W3[n]
  }

  data.frame(q1 = q1, q2 = q2, q3 = q3)
}


synergistic_collider <- function(N) {
  q1 <- numeric(N)
  q2 <- numeric(N)
  q3 <- numeric(N)
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)
  W3 <- rnorm(N, 0, 1)

  for (n in seq_len(N - 1)) {
    q1[n + 1] <- sin(q2[n] * q3[n]) + 0.001 * W1[n]
    q2[n + 1] <- 0.5 * q2[n] + 0.1 * W2[n]
    q3[n + 1] <- 0.5 * q3[n] + 0.1 * W3[n]
  }

  data.frame(q1 = q1, q2 = q2, q3 = q3)
}


redundant_collider <- function(N) {
  q1 <- numeric(N)
  q2 <- numeric(N)
  q3 <- numeric(N)
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)
  W3 <- rnorm(N, 0, 1)

  for (n in seq_len(N - 1)) {
    # original alternative: q1[n+1] = 0.3*q1[n] + 0.7*(sin(q2[n]*q3[n]) + 0.1*W1[n])
    q1[n + 1] <- 0.3 * q1[n] + (sin(q2[n] * q3[n]) + 0.001 * W1[n])
    q2[n + 1] <- 0.5 * q2[n] + 0.1 * W2[n]
    q3[n + 1] <- q2[n + 1]
  }

  data.frame(q1 = q1, q2 = q2, q3 = q3)
}
