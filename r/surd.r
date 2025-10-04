
# ===============================
# R implementation of surd (Synergy, Unique, Redundancy Decomposition)
# ===============================

# ---------- Entropy helpers ----------
entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

# Conditional entropy H(T|A)
cond_entropy <- function(p, target_dim, cond_dims) {
  margin_target <- apply(p, target_dim, sum)
  margin_cond <- apply(p, cond_dims, sum)
  joint <- p
  H <- 0
  # flatten to iterate
  idx <- which(joint > 0, arr.ind = TRUE)
  for (i in seq_len(nrow(idx))) {
    prob <- joint[idx[i, , drop = FALSE]]
    t_val <- idx[i, target_dim]
    cond_val <- idx[i, cond_dims, drop = FALSE]
    p_ta <- prob
    p_a <- margin_cond[matrix(cond_val, nrow=1)]
    if (p_a > 0) {
      H <- H - p_ta * log(p_ta / p_a)
    }
  }
  return(H)
}

# ---------- surd core (partial information decomposition) ----------
surd_core <- function(p) {
  # normalize
  p <- p + 1e-14
  p <- p / sum(p)
  
  Ntot <- length(dim(p))
  Nvars <- Ntot - 1
  Nt <- dim(p)[1]
  inds <- 2:Ntot
  
  # Info leak
  H <- entropy(apply(p, 1, sum))
  Hc <- cond_entropy(p, 1, inds)
  info_leak <- Hc / H
  
  # Target marginal
  p_s <- apply(p, 1, sum)
  
  # Compute specific MI for each subset of variables
  combs <- list()
  Is <- list()
  
  all_inds <- inds
  for (k in seq_along(all_inds)) {
    cmb <- combn(all_inds, k, simplify = FALSE)
    for (j in cmb) {
      combs <- append(combs, list(j))
      noj <- setdiff(all_inds, j)
      
      p_a <- apply(p, c(1, j), sum)
      p_as <- apply(p, c(j), sum)
      
      p_a_s <- sweep(p_as, 1, p_s, "/")
      p_s_a <- sweep(p_as, j, apply(p, j, sum), "/")
      
      # compute specific MI
      tmp <- array(0, dim = c(Nt))
      for (t in 1:Nt) {
        rows <- p_a_s[t, , drop=FALSE]
        ps_val <- p_s[t]
        val <- rows * (log(p_s_a[t, , drop=FALSE] + 1e-14) - log(ps_val + 1e-14))
        tmp[t] <- sum(val)
      }
      Is[[paste(j, collapse = "_")]] <- tmp
    }
  }
  
  # MI
  MI <- lapply(Is, function(v) sum(v * p_s))
  
  I_R <- list()
  I_S <- list()
  for (nm in names(Is)) I_R[[nm]] <- 0
  
  # process redundancy/synergy
  for (t in 1:Nt) {
    I1 <- sapply(Is, function(ii) ii[t])
    order_idx <- order(I1)
    labs <- names(Is)[order_idx]
    lens <- sapply(strsplit(labs, "_"), length)
    I1_sorted <- I1[order_idx]
    
    # eliminate dominated info
    for (l in 1:max(lens)) {
      idx_l2 <- which(lens == l+1)
      if (length(idx_l2) > 0) {
        Il1max <- max(I1_sorted[lens==l])
        bad <- idx_l2[I1_sorted[idx_l2] < Il1max]
        I1_sorted[bad] <- 0
      }
    }
    
    Di <- diff(c(0, I1_sorted))
    red_vars <- inds
    
    for (i in seq_along(labs)) {
      info <- Di[i] * p_s[t]
      vars <- strsplit(labs[i], "_")[[1]]
      if (length(vars) == 1) {
        key <- paste(red_vars, collapse = "_")
        I_R[[key]] <- I_R[[key]] + info
        red_vars <- setdiff(red_vars, as.integer(vars))
      } else {
        I_S[[labs[i]]] <- (I_S[[labs[i]]] %||% 0) + info
      }
    }
  }
  
  return(list(I_R=I_R, I_S=I_S, MI=MI, info_leak=info_leak))
}

# helper for default null coalesce
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------- wrapper function ----------
surd <- function(data, target, agents, bins=5, dt=1) {
  Nt <- nrow(data)
  df <- data.frame(
    T = data[[target]][(dt+1):Nt]
  )
  for (v in agents) {
    df[[v]] <- data[[v]][1:(Nt-dt)]
  }
  
  df_disc <- as.data.frame(lapply(df, function(x) {
    cut(x, breaks=bins, labels=FALSE, include.lowest=TRUE)
  }))
  
  tbl <- as.data.frame(table(df_disc))
  prob <- tbl$Freq / sum(tbl$Freq)
  dims <- rep(bins, ncol(df_disc))
  
  p <- array(0, dim=dims)
  idx <- as.matrix(tbl[,1:ncol(df_disc)])
  for (i in seq_len(nrow(tbl))) {
    p[matrix(idx[i,], nrow=1)] <- prob[i]
  }
  
  return(surd_core(p))
}

# ---------- Example ----------
# set.seed(123)
# df <- data.frame(
#   s = sin(1:100) + rnorm(100,0,0.1),
#   a1 = cos(1:100) + rnorm(100,0,0.1),
#   a2 = runif(100)
# )
# result <- surd(df, target="s", agents=c("a1","a2"), bins=4, dt=1)
# str(result)
