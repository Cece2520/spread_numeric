# Case 1|234|7

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# This program numerically attempts to rule out solutions to
# case 1|234|7 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_12347(mu, nu, a3, a7)
    
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
    
    asum = intersect(a3+a7, UNIT_INT)
    
    if isempty(asum)
        return false
    end

    a2 = a2_assume234(a3, mn, v)
    asum = intersect(asum + a2, UNIT_INT)
    if isempty(asum)
        return false
    end

    a4 = a4_assume1234(a3, mn, v)
    asum = intersect(asum + a4, UNIT_INT)
    if isempty(asum)
        return false
    end

    f3, g3 = fg3_assume23(mu, nu, a3, mn, v)
    f2, g2 = fg2_assume23(mu, nu, a3, f3, g3)
    f4, g4 = fg4_assume234(mu, nu, a2, f2, f3, g2, g3)
    f1, g1 = fg1_assume124(mu, nu, a4, f2, f4, g2, g4)
    
    if isempty(f1) || isempty(g1)
        return false
    end

    a1 = intersect(1-asum, UNIT_INT)
    avec = [a1, a2, a3, a4, 0, 0, a7]
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, false)
    f1bot, g1bot = fg2_assume23(mu, nu, a3, f3, g3)
    f1 = intersect(f1, f1bot)
    g1 = intersect(g1, g1bot)
    

    # Edge density equals the sum of the squares of eigenvalues

    if !density_feasible(mu, nu, avec)
        return false
    end
    
    # Check the inequalities and eigenvector equations
    
    fvec = [f1, f2, f3, f4, nothing, nothing, f7]
    gvec = [g1, g2, g3, g4, nothing, nothing, g7]

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

