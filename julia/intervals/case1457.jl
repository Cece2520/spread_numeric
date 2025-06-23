# Case 1|4|57

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# This program numerically attempts to rule out solutions to
# case 1|4|57 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_1457(mu, nu, a4, a5, c)
    
    # We ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3),
    # cannot be less than the upper bound
    # (1+sqrt(2))/2, or cannot satisfy the 
    # inequality mu + sqrt(mu-mu^2) >= spread
    # Some helper variables are defined
    
    u = mu-nu
    v = mu+nu
    mn = mu*nu
    s = c * v - nu

    if !mu_nu_feasible(mu, nu, u, s, c)
        return false
    end
    
    # We ignore cases where the weight sum exceeds 1
    # Apply the relevant formulas for a_i, f_i, g_i
    
    asum = intersect(a4+a5, UNIT_INT)
    if isempty(asum)
        return false
    end

    f4, g4 = fg3_assume23(mu, nu, a4, mn, c, v, s)
    f1, g1 = fg2_assume23(mu, nu, a4, f4, g4)
    
    if isempty(f1) || isempty(g1)
        return false
    end
    
    f5, g5 = fg2_assume2N4(mu, nu, a5, mn, c, v, s, false)
    f7, g7 = fg4_assume2N4(mu, nu, a5, f5, g5, false)
    
    if isempty(f7) || isempty(g7)
        return false
    end
    
    a7 = a4_assume12N4(a5, mn, v)
    asum = intersect(asum+a7, UNIT_INT)
    a1 = intersect(1-asum, UNIT_INT)
    
    if isempty(a1)
        return false
    end
    
    # Edge density equals the sum of the squares of eigenvalues

    avec = [a1, 0, 0, a4, a5, 0, a7]
    if !density_feasible(mu, nu, avec)
        return false
    end
    
    # Check the inequalities and eigenvector equations
    
    fvec = [f1, nothing, nothing, f4, f5, nothing, f7]
    gvec = [g1, nothing, nothing, g4, g5, nothing, g7]

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
    if !ellipse_feasible(mu, nu, fvec, gvec, c, s)
        return false
    end

    return true
end

