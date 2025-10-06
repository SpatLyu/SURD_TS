import numpy as np
from itertools import combinations as icmb
from typing import Tuple, Dict, Optional


def mylog(x):
    """
    Compute the logarithm in base 2 avoiding singularities.
    
    Parameters:
    - x (np.array): Input data.

    Returns:
    - np.array: Logarithm in base 2 of the input.
    """
    valid_indices = (x != 0) & (~np.isnan(x)) & (~np.isinf(x))
    
    log_values = np.zeros_like(x)
    log_values[valid_indices] = np.log2(x[valid_indices])
    
    return log_values


def entropy(p):
    """
    Compute the entropy of a discrete probability distribution function.

    Parameters:
    - p (np.array): Probability distribution of the signal.

    Returns:
    - float: Entropy of the given distribution.
    """
    return -np.sum(p * mylog(p))


def entropy_nvars(p, indices):
    """
    Compute the joint entropy for specific dimensions of a probability distribution.

    Parameters:
    - p (np.array): N-dimensional joint probability distribution.
    - indices (tuple): Dimensions over which the entropy is to be computed.

    Returns:
    - float: Joint entropy for specified dimensions.

    Example: compute the joint entropy H(X0,X3,X7)
    >>> entropy_nvars(p, (0,3,7))
    """
    
    excluded_indices = tuple(set(range(p.ndim)) - set(indices))
    marginalized_distribution = p if not excluded_indices else p.sum(axis=excluded_indices)

    return entropy(marginalized_distribution)


def cond_entropy(p, target_indices, conditioning_indices):
    """
    Compute the conditional entropy between two sets of variables.

    Parameters:
    - p (np.array): N-dimensional joint probability distribution.
    - target_indices (tuple): Variables for which entropy is to be computed.
    - conditioning_indices (tuple): Conditioning variables.

    Returns:
    - float: Conditional entropy.

    Example: compute the conditional entropy H(X0,X2|X7)
    >>> cond_entropy(p, (0, 2), (7,))
    """
    joint_entropy = entropy_nvars(p, set(target_indices) | set(conditioning_indices))
    conditioning_entropy = entropy_nvars(p, conditioning_indices)

    return joint_entropy - conditioning_entropy


def mutual_info(p, set1_indices, set2_indices):
    """
    Compute the mutual information between two sets of variables.

    Parameters:
    - p (np.array): N-dimensional joint probability distribution.
    - set1_indices (tuple): Indices of the first set of variables.
    - set2_indices (tuple): Indices of the second set of variables.

    Returns:
    - float: Mutual information.

    Example: compute the mutual information I(X0,X5;X4,X2)
    >>> mutual_info(p, (0, 5), (4, 2))
    """
    entropy_set1 = entropy_nvars(p, set1_indices)
    conditional_entropy = cond_entropy(p, set1_indices, set2_indices)

    return entropy_set1 - conditional_entropy


def cond_mutual_info(p, ind1, ind2, ind3):
    """
    Compute the conditional mutual information between two sets of variables 
    conditioned to a third set.

    Parameters:
    - p (np.array): N-dimensional joint probability distribution.
    - ind1 (tuple): Indices of the first set of variables.
    - ind2 (tuple): Indices of the second set of variables.
    - ind3 (tuple): Indices of the conditioning variables.

    Returns:
    - float: Conditional mutual information.

    Example: compute the conditional mutual information I(X0,X5;X4,X2|X1)
    cond_mutual_info(p, (0, 5), (4, 2), (1,)))
    """
    # Merge indices of ind2 and ind3
    combined_indices = tuple(set(ind2) | set(ind3))
    
    # Compute conditional mutual information
    return cond_entropy(p, ind1, ind3) - cond_entropy(p, ind1, combined_indices)


def transfer_entropy(p, target_var):
    """
    Calculate the transfer entropy from each input variable to the target variable.

    Parameters:
    - p (np.array): Multi-dimensional array containing the pdfs of the variables.
      The first dimension corresponds to the index of the variable:
          p[0]  -> target variable (in future)
          p[1:] -> input variables (at present time)

    Returns:
    - np.array: Transfer entropy values for each input variable.
    """
    num_vars = len(p.shape) - 1  # Excluding the future variable
    TE = np.zeros(num_vars)
    
    for i in range(1, num_vars + 1):
        # The indices for the present variables
        present_indices = tuple(range(1, num_vars + 1))
        
        # The indices for the present variables excluding the i-th variable
        # conditioning_indices = tuple([target_var] + [j for j in range(1, num_vars + 1) if j != i])
        conditioning_indices = tuple([target_var] + [j for j in range(1, num_vars + 1) if j != i and j != target_var])
        
        # Conditional entropy of the future state of the target variable given its own past
        cond_ent_target_given_past = cond_entropy(p, (0,), conditioning_indices)
        
        # Conditional entropy of the future state of the target variable given its own past and the ith input variable
        cond_ent_target_given_past_and_input = cond_entropy(p, (0,), present_indices)
        
        # Transfer entropy calculation
        TE[i-1] = cond_ent_target_given_past - cond_ent_target_given_past_and_input
    
    return TE


def surd(p: np.ndarray, max_combs: Optional[int] = None) -> Tuple[Dict, Dict, Dict, float]:
    '''
    Unified SURD (Synergy, Unique, Redundancy Decomposition)
    for both standard and high-dimensional information decomposition.

    This function decomposes the mutual information between a target variable
    (the first dimension of p) and a set of agent variables (remaining dimensions)
    into redundancy, synergy, and unique information terms.

    Behavior:
    ----------
    - If max_combs is None:
        Perform the standard SURD decomposition.
    - If max_combs is specified:
        Perform the high-dimensional extension of SURD,
        computing redundancy and synergy up to the specified combination order.

    Parameters:
    -----------
    p : np.ndarray
        N-dimensional joint probability distribution. The first dimension
        corresponds to the target variable, and the remaining ones correspond
        to agent variables.
    max_combs : int or None, optional
        Maximum combination order for high-dimensional synergy computations.
        If None, standard SURD behavior is used.

    Returns:
    --------
    I_R : dict
        Redundancy and unique information for each variable combination.
    I_S : dict
        Synergy information for each variable combination.
    MI : dict
        Mutual information for each combination of agent variables.
    info_leak : float
        Information leak ratio Hc/H for the target variable.

    Example:
    --------
    >>> I_R, I_S, MI, info_leak = surd(p)               # Standard SURD
    >>> I_R, I_S, MI, info_leak = surd(p, max_combs=3)  # High-dimensional SURD
    '''

    # Avoid log singularities
    p = np.maximum(p, 1e-14)
    p = p / p.sum()

    Ntot = p.ndim
    Nvars = Ntot - 1
    Nt = p.shape[0]
    inds = range(1, Ntot)

    # --- Information leak computation (shared for both modes)
    H = entropy_nvars(p, (0,))
    Hc = cond_entropy(p, (0,), range(1, Ntot))
    info_leak = Hc / H

    # --- Standard SURD behavior
    if max_combs is None:
        p_s = p.sum(axis=(*inds,), keepdims=True)
        combs, Is = [], {}

        for i in inds:
            for j in list(icmb(inds, i)):
                combs.append(j)
                noj = tuple(set(inds) - set(j))

                p_a = p.sum(axis=(0, *noj), keepdims=True)
                p_as = p.sum(axis=noj, keepdims=True)

                p_a_s = p_as / p_s
                p_s_a = p_as / p_a

                Is[j] = (p_a_s * (mylog(p_s_a) - mylog(p_s))).sum(axis=j).ravel()

        MI = {k: (Is[k] * p_s.squeeze()).sum() for k in Is.keys()}
        I_R = {cc: 0 for cc in combs}
        I_S = {cc: 0 for cc in combs[Nvars:]}

        for t in range(Nt):
            I1 = np.array([ii[t] for ii in Is.values()])
            i1 = np.argsort(I1)
            lab = [combs[i_] for i_ in i1]
            lens = np.array([len(l) for l in lab])

            I1 = I1[i1]
            for l in range(1, lens.max()):
                inds_l2 = np.where(lens == l + 1)[0]
                Il1max = I1[lens == l].max()
                inds_ = inds_l2[I1[inds_l2] < Il1max]
                I1[inds_] = 0

            i1 = np.argsort(I1)
            lab = [lab[i_] for i_ in i1]
            Di = np.diff(I1[i1], prepend=0.0)
            red_vars = list(inds)

            for i_, ll in enumerate(lab):
                info = Di[i_] * p_s.squeeze()[t]
                if len(ll) == 1:
                    I_R[tuple(red_vars)] += info
                    red_vars.remove(ll[0])
                else:
                    I_S[ll] += info

        return I_R, I_S, MI, info_leak

    # --- High-dimensional SURD behavior
    else:
        max_inds = range(1, max_combs + 1)
        tot_inds = range(1, Ntot)

        p_s = p.sum(axis=(*tot_inds,), keepdims=True)
        combs, Is = [], {}
        red_combs = []

        for i in max_inds:
            for j in list(icmb(tot_inds, i)):
                combs.append(j)

                p_a = p.sum(axis=tuple(set(tot_inds) - set(j)), keepdims=True)
                p_as = p.sum(axis=tuple(set(tot_inds) - set(j)), keepdims=True)

                p_a_s = p_as / p_s
                p_s_a = p_as / p_a

                Is[j] = (p_a_s * (mylog(p_s_a) - mylog(p_s))).sum(axis=j).ravel()

        MI = {k: (Is[k] * p_s.squeeze()).sum() for k in Is.keys()}

        for i in tot_inds:
            for j in list(icmb(tot_inds, i)):
                red_combs.append(j)
        I_R = {cc: 0 for cc in red_combs}
        I_S = {cc: 0 for cc in combs[Nvars:]}

        for t in range(Nt):
            I1 = np.array([ii[t] for ii in Is.values()])
            i1 = np.argsort(I1)
            lab = [combs[i_] for i_ in i1]
            lens = np.array([len(l) for l in lab])

            I1 = I1[i1]
            for l in range(1, lens.max()):
                inds_l2 = np.where(lens == l + 1)[0]
                Il1max = I1[lens == l].max()
                inds_ = inds_l2[I1[inds_l2] < Il1max]
                I1[inds_] = 0

            i1 = np.argsort(I1)
            lab = [lab[i_] for i_ in i1]
            Di = np.diff(I1[i1], prepend=0.0)
            red_vars = list(tot_inds)

            for i_, ll in enumerate(lab):
                info = Di[i_] * p_s.squeeze()[t]
                if len(ll) == 1:
                    I_R[tuple(red_vars)] += info
                    red_vars.remove(ll[0])
                else:
                    I_S[ll] += info

        return I_R, I_S, MI, info_leak