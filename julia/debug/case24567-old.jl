# 24|567

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")


# numerically attempts to rule out solutions to 
# the constraints, falling into the given intervals

function is_feasible_24567_old(mu, nu, a2, a6)
    
    # first, we ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3)
    # some helper variables are also used
    
    u = mu-nu
    if isempty( intersect(u - SPR_MAX, POS))
        return false
    end
    v = mu+nu
    mn = mu*nu
    
    
    # ignore cases where weight sum exceeds 1
    
    asum = intersect(a2+a6, UNIT_INT)
    
    if isempty(asum)
        return false
    end

    # again, weight sum cannot exceed 1; 
    # apply formulas for other ai's, fi's, gi's
    
    a5 = a2_assume234(a6, mn, v)
    asum = intersect(asum + a5, UNIT_INT)
    if isempty(asum)
        return false
    end

    f2, g2 = fg2_assume2N4(mu, nu, a2, mn, v)
    f4, g4 = fg4_assume2N4(mu, nu, a2, f2, g2)
    
    if isempty(f4) || isempty(g4)
        return false
    end

    f6, g6 = fg3_assume23(mu, nu, a6, mn, v, false)
    f5, g5 = fg2_assume23(mu, nu, a6, f6, g6, false)
    f7, g7 = fg4_assume234(mu, nu, a5, f5, f6, g5, g6, false)
    
    if isempty(f7) || isempty(g7)
        return false
    end

    a7 = a4_assumeN2347(a6, mn, v)
    asum = intersect(asum + a7, UNIT_INT)
    
    a4 = intersect(1-asum, UNIT_INT)
    if isempty(a4)
        return false
    end

    # double-check the eigenvector eq'ns
    
    avec = [0, a2, 0, a4, a5, a6, a7]
    
    
    # check eigenvector equations
    
    fvec = [nothing, f2, nothing, f4, f5, f6, f7]
    gvec = [nothing, g2, nothing, g4, g5, g6, g7]
    
    if !fg_row_feasible(mu, nu, fvec, gvec, avec)
        return false
    end
    
    # might as well also check the norms and ellipse equations
    
    if !norm_feasible(fvec, gvec, avec)
        return false
    end
    if !ellipse_feasible(mu, nu, fvec, gvec, u)
        return false
    end

    return true
end

