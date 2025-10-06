# R implementation of SURD (Synergy, Unique, Redundancy Decomposition)

# ---------- Utilities ----------

#' Compute the logarithm in base 2 avoiding singularities
#' 
#' @param x Numeric vector or matrix input data
#' @return Numeric vector or matrix with logarithm in base 2 of the input
mylog <- function(x) {
  # 确保输入是数值类型
  if (is.list(x)) {
    x <- unlist(x)
  }
  
  # 转换为数值向量/矩阵
  x <- as.numeric(x)
  
  valid_indices <- (x != 0) & (!is.na(x)) & (!is.infinite(x))
  
  log_values <- numeric(length(x))
  dim(log_values) <- dim(x)
  log_values[valid_indices] <- log2(x[valid_indices])
  
  return(log_values)
}

#' Compute the entropy of a discrete probability distribution function
#' 
#' @param p Probability distribution of the signal
#' @return Entropy of the given distribution
entropy <- function(p) {
  # 确保 p 是数值向量
  if (is.list(p)) {
    p <- unlist(p)
  }
  p <- as.numeric(p)
  
  return(-sum(p * mylog(p), na.rm = TRUE))
}

#' Compute the joint entropy for specific dimensions of a probability distribution
#' 
#' @param p N-dimensional joint probability distribution (array)
#' @param indices Dimensions over which the entropy is to be computed
#' @return Joint entropy for specified dimensions
entropy_nvars <- function(p, indices) {
  all_dims <- seq_len(length(dim(p)))
  excluded_indices <- setdiff(all_dims, indices)
  
  # Marginalize over excluded dimensions
  if (length(excluded_indices) > 0) {
    # 使用更安全的方式计算边缘分布
    marginalized_distribution <- apply(p, indices, sum, simplify = FALSE)
    
    # 确保结果是一个数值数组
    if (is.list(marginalized_distribution)) {
      marginalized_distribution <- unlist(marginalized_distribution)
    }
    
    # 确保结果是一个数组
    if (!is.array(marginalized_distribution)) {
      # 如果结果不是数组，手动构建
      kept_dims <- dim(p)[indices]
      marginalized_distribution <- array(
        as.numeric(marginalized_distribution), 
        dim = kept_dims
      )
    }
  } else {
    marginalized_distribution <- p
  }
  
  return(entropy(marginalized_distribution))
}

#' Compute the conditional entropy between two sets of variables
#' 
#' @param p N-dimensional joint probability distribution (array)
#' @param target_indices Variables for which entropy is to be computed
#' @param conditioning_indices Conditioning variables
#' @return Conditional entropy
cond_entropy <- function(p, target_indices, conditioning_indices) {
  joint_indices <- union(target_indices, conditioning_indices)
  joint_entropy <- entropy_nvars(p, joint_indices)
  conditioning_entropy <- entropy_nvars(p, conditioning_indices)
  
  return(joint_entropy - conditioning_entropy)
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
  
  return(p)
}

# ---------- SURD main ----------

surd <- function(p) {
  #' 分解目标变量与一组代理变量之间的互信息为冗余、协同和独特信息
  #'
  #' @param p array - 多维直方图数组，第一维代表目标变量，后续维度代表代理变量
  #' @return list - 包含冗余信息、协同信息、互信息和信息泄漏的列表
  
  # 确保概率分布中没有零值
  p <- p + 1e-14
  p <- p / sum(p)
  
  # 维度信息
  dims <- dim(p)
  Ntot <- length(dims)  # 总维度数
  Nvars <- Ntot - 1     # 代理变量数
  Nt <- dims[1]         # 目标变量状态数
  inds <- 2:Ntot        # 代理变量索引
  
  # 计算信息泄漏
  H <- entropy_nvars(p, 1)
  Hc <- cond_entropy(p, 1, inds)
  info_leak <- Hc / H
  
  # 计算目标变量的边缘分布
  p_s <- apply(p, 1, sum)
  
  # 准备特定互信息计算
  combs <- list()
  Is <- list()
  
  # 遍历所有代理变量组合
  for(i in inds) {
    other_inds <- setdiff(inds, i)
    for(j_size in 1:length(other_inds)) {
      combs_j <- combn(other_inds, j_size, simplify = FALSE)
      for(j in combs_j) {
        comb_key <- c(i, j)
        comb_key_str <- paste(sort(comb_key), collapse = ",")
        
        # 计算当前组合的联合和条件分布
        sum_axes <- setdiff(inds, comb_key)
        if(length(sum_axes) == 0) {
          p_as <- p
        } else {
          p_as <- apply(p, c(1, comb_key), sum, simplify = FALSE)
          # 确保结果是数值数组
          if (is.list(p_as)) {
            p_as <- unlist(p_as)
          }
          if (!is.array(p_as)) {
            kept_dims <- dim(p)[c(1, comb_key)]
            p_as <- array(as.numeric(p_as), dim = kept_dims)
          }
        }
        
        p_a <- apply(p_as, comb_key, sum)
        p_s_a <- sweep(p_as, comb_key, p_a, "/")
        p_a_s <- sweep(p_as, 1, p_s, "/")
        
        # 计算特定互信息
        log_ratio <- log(p_s_a) - log(p_s)
        Is[[comb_key_str]] <- apply(p_a_s * log_ratio, 1, sum)
      }
    }
  }
  
  # 计算每个组合的互信息
  MI <- list()
  for(key in names(Is)) {
    MI[[key]] <- sum(Is[[key]] * p_s)
  }
  
  # 初始化冗余和协同项
  I_R <- list()
  I_S <- list()
  
  # 处理目标变量的每个值
  for(t in 1:Nt) {
    # 提取当前目标值的特定互信息
    I1 <- numeric()
    lab <- list()
    
    for(key in names(Is)) {
      I1 <- c(I1, Is[[key]][t])
      lab <- c(lab, list(as.numeric(strsplit(key, ",")[[1]])))
    }
    
    # 排序特定互信息
    i1 <- order(I1)
    I1_sorted <- I1[i1]
    lab_sorted <- lab[i1]
    lens <- sapply(lab_sorted, length)
    
    # 基于现有最大值更新特定互信息
    for(l in 1:(max(lens)-1)) {
      inds_l2 <- which(lens == l + 1)
      Il1max <- max(I1_sorted[lens == l])
      for(idx in inds_l2) {
        if(I1_sorted[idx] < Il1max) {
          I1_sorted[idx] <- 0
        }
      }
    }
    
    # 重新排序更新后的特定互信息值
    i1_new <- order(I1_sorted)
    I1_final <- I1_sorted[i1_new]
    lab_final <- lab_sorted[i1_new]
    
    # 计算排序后特定互信息值的差异
    Di <- diff(c(0, I1_final))
    red_vars <- inds
    
    # 将互信息分配到冗余和协同项
    for(i in seq_along(lab_final)) {
      info <- Di[i] * p_s[t]
      ll <- lab_final[[i]]
      
      if(length(ll) == 1) {
        key <- paste(sort(red_vars), collapse = ",")
        if(is.null(I_R[[key]])) I_R[[key]] <- 0
        I_R[[key]] <- I_R[[key]] + info
        red_vars <- setdiff(red_vars, ll)
      } else {
        key <- paste(sort(ll), collapse = ",")
        if(is.null(I_S[[key]])) I_S[[key]] <- 0
        I_S[[key]] <- I_S[[key]] + info
      }
    }
  }
  
  return(list(I_R = I_R, I_S = I_S, MI = MI, info_leak = info_leak))
}