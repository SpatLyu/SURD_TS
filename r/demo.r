setwd('./project/SURD_TS/')
source('./r/surd.r')

# ---------- Example ----------
set.seed(123)
df <- data.frame(
  s  = sin(1:200) + rnorm(200,0,0.1),
  a1 = cos(1:200) + rnorm(200,0,0.1),
  a2 = runif(200)
)

p <- create_pfm(df, target="s", agents=c("a1","a2"), bins=4, dt=1)

result <- surd(p)
str(result)