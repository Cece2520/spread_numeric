# Case 1|4|7

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")

# This program numerically attempts to rule out solutions to
# case 1|4|7 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_147(mu, nu, a4, a7)
    
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
    
    asum = intersect(a4+a7, UNIT_INT)
    
    if isempty(asum)
        return false
    end
    
    f4, g4 = fg3_assume23(mu, nu, a4, mn, v)
    f1, g1 = fg2_assume23(mu, nu, a4, f4, g4)
    
    if isempty(f4) || isempty(g4)
        return false
    end
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, false)
    f1bot, g1bot = fg2_assume23(mu, nu, a7, f7, g7)
    
    f1 = intersect(f1, f1bot)
    g1 = intersect(g1, g1bot)
    
    if isempty(f1) || isempty(g1)
        return false
    end
    
    a1 = intersect(interval(1)-(a4+a7), UNIT_INT)
    
    # Edge density equals the sum of the squares of eigenvalues

    avec = [a1, 0, 0, a4, 0, 0, a7]
    if !density_feasible(mu, nu, avec)
        return false
    end
    
    # Check the inequalities and eigenvector equations
    
    fvec = [f1, nothing, nothing, f4, nothing, nothing, f7]
    gvec = [g1, nothing, nothing, g4, nothing, nothing, g7]

    if !fg_ineq_feasible(fvec,gvec)
        return false
    end
    
    if !fg_row_feasible(mu, nu, fvec, gvec, avec)
        return false
    end
    
    # Check the norm and ellipse equations

    if !norm_feasible(fvec, gvec, avec)
        return false
    end
    
    if !ellipse_feasible(mu, nu, fvec, gvec, u)
        return false
    end
    
    return true
end

