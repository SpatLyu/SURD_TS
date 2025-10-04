# R implementation of SURD (Synergy, Unique, Redundancy Decomposition)

# ---------- Utilities ----------

mylog <- function(x) {
  y <- x
  y[y <= 0] <- 1e-14
  log(y)
}

entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

cond_entropy <- function(p, target_dim, cond_dims) {
  # H(T|A) = H(T,A) - H(A)
  joint <- margin.table(p, c(target_dim, cond_dims))
  pa <- margin.table(p, cond_dims)
  Hjoint <- entropy(c(joint))
  Ha <- entropy(c(pa))
  Hjoint - Ha
}

# ---------- Specific MI computation ----------

compute_Is <- function(p, inds) {
  Nt <- dim(p)[1]
  p_s <- apply(p, 1, sum) # marginal target distribution
  combs <- list()
  Is <- list()
  
  for (k in 1:length(inds)) {
    cmbs <- combn(inds, k, simplify=FALSE)
    for (j in cmbs) {
      combs[[length(combs)+1]] <- j
      
      # joint distribution P(S, Aj)
      p_as <- apply(p, c(1,j), sum)
      
      # --- 保证 p_as 一定是 (Nt × M) 的矩阵 ---
      if (is.null(dim(p_as))) {
        p_as <- array(p_as, dim=c(Nt,1))
      } else if (length(dim(p_as)) == 1) {
        p_as <- matrix(p_as, nrow=Nt)
      }
      
      # Specific MI for each target value
      val <- numeric(Nt)
      for (t in 1:Nt) {
        psa <- p_as[t, ]
        psa <- psa / sum(psa)
        ps  <- p_s[t]
        if (ps > 0) {
          val[t] <- sum(psa * (mylog(psa) - log(ps)))
        }
      }
      
      Is[[paste(j, collapse="-")]] <- val
    }
  }
  
  list(combs=combs, Is=Is, p_s=p_s)
}

# ---------- Compute the joint PMF of variables ----------

create_pfm <- function(df, target, agents, dt=1, bins=4) {
  Nt <- nrow(df)
  s <- df[[target]]
  a <- lapply(agents, function(x) df[[x]])
  
  V <- data.frame(
    s = s[(dt+1):Nt],
    as.data.frame(lapply(a, function(ai) ai[1:(Nt-dt)]))
  )
  
  # Discretize into bins
  df_disc <- as.data.frame(lapply(V, function(x) {
    as.integer(cut(x, breaks=bins, labels=FALSE, include.lowest=TRUE))
  }))
  
  tbl <- as.data.frame(table(df_disc))
  prob <- tbl$Freq / sum(tbl$Freq)
  dims <- rep(bins, ncol(df_disc))
  
  p <- array(0, dim=dims)
  idx <- as.matrix(tbl[,1:ncol(df_disc)])
  for (i in seq_len(nrow(tbl))) {
    p[matrix(as.numeric(idx[i,]), nrow=1)] <- prob[i]
  }
  
  p
}

# ---------- SURD main ----------

surd <- function(p) {
  # normalize
  p <- p + 1e-14
  p <- p / sum(p)
  
  Ntot <- length(dim(p))
  Nvars <- Ntot - 1
  Nt <- dim(p)[1]
  inds <- 2:Ntot
  
  # info leak
  H <- entropy(apply(p, 1, sum))
  Hc <- cond_entropy(p, 1, inds)
  info_leak <- Hc / H
  
  # specific MI
  tmp <- compute_Is(p, inds)
  combs <- tmp$combs
  Is <- tmp$Is
  p_s <- tmp$p_s
  
  # MI for each combination
  MI <- list()
  for (k in names(Is)) {
    MI[[k]] <- sum(Is[[k]] * p_s)
  }
  
  I_R <- list()
  for (cc in combs) {
    I_R[[paste(cc, collapse="-")]] <- 0
  }
  I_S <- list()
  for (cc in combs[(Nvars+1):length(combs)]) {
    I_S[[paste(cc, collapse="-")]] <- 0
  }
  
  for (t in 1:Nt) {
    I1 <- sapply(Is, function(v) v[t])
    ord <- order(I1)
    lab <- names(Is)[ord]
    lens <- sapply(strsplit(lab,"-"), length)
    I1s <- I1[ord]
    
    # update based on max lower-order
    for (l in 1:max(lens)) {
      inds_l2 <- which(lens == l+1)
      if (length(inds_l2)>0) {
        Il1max <- max(I1s[lens==l])
        bad <- inds_l2[I1s[inds_l2] < Il1max]
        I1s[bad] <- 0
      }
    }
    
    ord2 <- order(I1s)
    lab2 <- lab[ord2]
    Di <- diff(c(0,I1s[ord2]))
    red_vars <- inds
    
    for (i in seq_along(lab2)) {
      info <- Di[i] * p_s[t]
      vars <- strsplit(lab2[i],"-")[[1]]
      if (length(vars)==1) {
        nm <- paste(red_vars, collapse="-")
        I_R[[nm]] <- I_R[[nm]] + info
        red_vars <- setdiff(red_vars, as.integer(vars))
      } else {
        I_S[[lab2[i]]] <- I_S[[lab2[i]]] + info
      }
    }
  }
  
  list(I_R=I_R, I_S=I_S, MI=MI, info_leak=info_leak)
}
