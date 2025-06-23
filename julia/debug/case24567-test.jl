# Case 24|567

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# This program numerically attempts to rule out solutions to
# case 24|567 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_24567_test(mu, nu, a2, a6)
    
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
        return (false, "mu-nu")
    end

    # We ignore cases where the weight sum exceeds 1
    # Apply the relevant formulas for a_i, f_i, g_i
    
    asum = intersect(a2+a6, UNIT_INT)
    
    if isempty(asum)
        return (false, "asum")
    end

    a5 = a2_assume234(a6, mn, v)
    asum = intersect(asum + a5, UNIT_INT)
    if isempty(asum)
        return (false, "asum")
    end

    f2, g2 = fg2_assume2N4(mu, nu, a2, mn, v)
    f4, g4 = fg4_assume2N4(mu, nu, a2, f2, g2)
    
    if isempty(f4) || isempty(g4)
        return (false, "24")
    end

    f6, g6 = fg3_assume23(mu, nu, a6, mn, v, false)
    f5, g5 = fg2_assume23(mu, nu, a6, f6, g6, false)
    f7, g7 = fg4_assume234(mu, nu, a5, f5, f6, g5, g6, false)
    
    if isempty(f7) || isempty(g7)
        return (false, "567")
    end

    a7 = a4_assumeN2347(a6, mn, v)
    asum = intersect(asum + a7, UNIT_INT)
    
    a4 = intersect(1-asum, UNIT_INT)
    if isempty(a4)
        return (false, "47")
    end

    # Update each defined a_i based on sum a_i = 1

    a5 = intersect(a5, (1-a2-a4-a6-a7))
    if isempty(a5)
        return (false, "457")
    end

    a7 = intersect(a7, (1-a2-a4-a5-a6))
    if isempty(a5)
        return (false, "457")
    end
    
    a4 = intersect(a4, (1-a2-a5-a6-a7))
    if isempty(a5)
        return (false, "457")
    end

    # Use alternate definition for a_4, a_7

    fvec = [nothing, f2, nothing, f4, f5, f6, f7]
    gvec = [nothing, g2, nothing, g4, g5, g6, g7]
    
    a4n, a7n = a47_assume2457(0,a2,0,a5,a6,mu,nu,fvec,gvec)

    a4 = intersect(a4, a4n)
    if isempty(a4)
        return (false, "47n")
    end

    a7 = intersect(a7, a7n)
    if isempty(a7)
        return (false, "47n")
    end
    if isempty( intersect(a4 + a7, interval(1) - a2 - a5 - a6))
        return (false, "47n")             
    end

    # Edge density equals the sum of the squares of eigenvalues
    
    avec = [0, a2, 0, a4, a5, a6, a7]
    if !density_feasible(mu, nu, avec)
        return (false, "dens")
    end
    
    # Check the inequalities and eigenvector equations

    if !fg_ineq_feasible(fvec,gvec)
        return (false, "ineq")
    end
    if !fg_row_feasible47(mu, nu, fvec, gvec, avec)
        return (false, "row47")
    end
    
    # Check the norm and ellipse equations
    
    if !norm_feasible47(fvec, gvec, avec)
        return (false, "norm47")
    end
    if !ellipse_feasible(mu, nu, fvec, gvec, u)
        return (false, "ell")
    end

    return (true, "isfeas")
end

