# Case 1|234|57

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# This program numerically attempts to rule out solutions to
# case 1|234|57 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_123457(mu, nu, a3, a5)
    
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
    
    asum = intersect(a3+a5, UNIT_INT)
    
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

    f5, g5 = fg2_assume2N4(mu, nu, a5, mn, v, false)
    f7, g7 = fg4_assume2N4(mu, nu, a5, f5, g5, false)
    
    if isempty(f7) || isempty(g7)
        return false
    end

    f3, g3 = fg3_assume23(mu, nu, a3, mn, v)
    f2, g2 = fg2_assume23(mu, nu, a3, f3, g3)
    f4, g4 = fg4_assume234(mu, nu, a2, f2, f3, g2, g3)
    f1, g1 = fg1_assume124(mu, nu, a4, f2, f4, g2, g4)
    
    if isempty(f1) || isempty(g1)
        return false
    end

    a7 = a4_assume12N4(a5, mn, v)

    asum = intersect(asum + a7, UNIT_INT)
    
    a1 = intersect(1 - asum, UNIT_INT)
    if isempty(a1)
        return false
    end

    # Use alternate definition for a_4, a_7

    fvec = [f1, f2, f3, f4, f5, nothing, f7]
    gvec = [g1, g2, g3, g4, g5, nothing, g7]
    
    a4n, a7n = a47_assume2457(a1,a2,a3,a5,0,mu,nu,fvec,gvec)

    a4 = intersect(a4, a4n)
    if isempty(a4)
        return false
    end

    a7 = intersect(a7, a7n)
    if isempty(a7)
        return false
    end

    if isempty( intersect(a4 + a7, 1 - a1 - a2 - a3 - a5))
        return false             
    end

    # Edge density equals the sum of the squares of eigenvalues
    
    avec = [a1, a2, a3, a4, a5, 0, a7]
    if !density_feasible(mu, nu, avec)
        return false
    end

    # Check the inequalities and eigenvector equations

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

