# Case 1|4|567

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# This program numerically attempts to rule out solutions to
# case 1|4|567 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_14567(mu, nu, a4, a6)
    
    # We ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3),
    # cannot be less than the upper bound
    # (1+sqrt(2))/2, or cannot satisfy the 
    # inequality mu + sqrt(mu-mu^2) >= spread
    # Some helper variables are defined
    
    u = mu-nu
    v = mu+nu
    mn = mu*nu

    if !mu_nu_feasible(mu, nu, u)
        return false
    end
    
    # We ignore cases where the weight sum exceeds 1
    # Apply the relevant formulas for a_i, f_i, g_i
    
    asum = intersect(a4+a6, UNIT_INT)
    
    if isempty(asum)
        return false
    end

    a5 = a2_assume234(a6, mn, v)
    asum = intersect(asum + a5, UNIT_INT)
    if isempty(asum)
        return false
    end

    a7 = a4_assume1234(a6, mn, v)
    asum = intersect(asum + a7, UNIT_INT)
    if isempty(asum)
        return false
    end

    f6, g6 = fg3_assume23(mu, nu, a6, mn, v, false)
    f5, g5 = fg2_assume23(mu, nu, a6, f6, g6, false)
    f7, g7 = fg4_assume234(mu, nu, a5, f5, f6, g5, g6, false)
    f1, g1 = fg1_assume124(mu, nu, a7, f5, f7, g5, g7)
    
    if isempty(f1) || isempty(g1)
        return false
    end

    f4, g4 = fg3_assume23(mu, nu, a4, mn, v)
    f1top, g1top = fg2_assume23(mu, nu, a4, f4, g4)
    f1 = intersect(f1, f1top)
    g1 = intersect(g1, g1top)
    
    if isempty(f1) || isempty(g1)
        return false
    end

    a1 = intersect(1 - asum, UNIT_INT)
    

    # Edge density equals the sum of the squares of eigenvalues

    avec = [a1, 0, 0, a4, a5, a6, a7]
    if !density_feasible(mu, nu, avec)
        return false
    end
    
    # Check the inequalities and eigenvector equations
    
    fvec = [f1, nothing, nothing, f4, f5, f6, f7]
    gvec = [g1, nothing, nothing, g4, g5, g6, g7]

    if !fg_ineq_feasible(fvec,gvec)
        return false
    end
    if !fg_row_feasible1(mu, nu, fvec, gvec, avec)
        return false
    end
    
    # Check the norm and ellipse equations
    
    if !norm_feasible1(fvec, gvec, avec)
        return false
    end
    if !ellipse_feasible(mu, nu, fvec, gvec, u)
        return false
    end

    return true
end

