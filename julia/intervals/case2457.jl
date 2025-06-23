# Case 24|57

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# This program numerically attempts to rule out solutions to
# case 24|57 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible_2457(mu, nu, a2, a5, c)
    
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
    
    asum = intersect(a2+a5, UNIT_INT)
    
    if isempty(asum)
        return false
    end

    f2, g2 = fg2_assume2N4(mu, nu, a2, mn, c, v, s)
    f4, g4 = fg4_assume2N4(mu, nu, a2, f2, g2)
    
    if isempty(f4) || isempty(g4)
        return false
    end
    
    f5, g5 = fg2_assume2N4(mu, nu, a5, mn, c, v, s, false)
    f7, g7 = fg4_assume2N4(mu, nu, a5, f5, g5, false)
    
    if isempty(f7) || isempty(g7)
        return false
    end
    
    # The equation mu*f4 = a5*f5 + a7*f7 implies a7 = (mu*f4 - a5*f5)/f7
    # The equation nu*g4 = a5*g5 + a7*g7 implies a7 = (nu*g4 - a5*g5)/g7
    # Similar equations hold for a4
    
    a7 = intersect( intersect((mu*f4-a5*f5)/f7, (nu*g4-a5*g5)/g7), UNIT_INT)
    a4 = intersect( intersect((mu*f7-a2*f2)/f4, (nu*g7-a2*g2)/g4), intersect(UNIT_INT, 1-(a2+a5+a7)))

    a7 = intersect(a7, 1-(a2+a4+a5))
    
    if isempty(a4) || isempty(a7)
        return false
    end

    avec = [0, a2, 0, a4, a5, 0, a7]


    # Edge density equals the sum of the squares of eigenvalues

    if !density_feasible(mu, nu, avec)
        return false
    end

    # Check the inequalities and eigenvector equations
    
    fvec = [nothing, f2, nothing, f4, f5, nothing, f7]
    gvec = [nothing, g2, nothing, g4, g5, nothing, g7]

    if !fg_ineq_feasible(fvec,gvec)
        return false
    end
    if !fg_row_feasible47(mu, nu, fvec, gvec, avec)
        return false
    end
    
    # Check the norm and ellipse equations
    
    if !norm_feasible47(fvec, gvec, avec)
        return false
    end
    if !ellipse_feasible(mu, nu, fvec, gvec, c, s)
        return false
    end

    return true
end

