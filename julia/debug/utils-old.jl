using IntervalArithmetic

SPR_MAX = sqrt(interval(4)/3)

NULL_INT = emptyinterval()
UNIT_INT = interval(0,1)
GEQ_ONE = interval(1, Inf)
POS = interval(0,Inf)

MAX_DEPTH = 50

# these are helper methods meant to be used when 
# handling cases that share formulas 

# every method begins with the variable to be returned 
# and ends with a string indicating assumptions made
# most often, the string is a list of indices w/ weight
# assumed to be positive
# when weights are assumed 0, the string corresponds to 
# an interval in increasing order; in the interval: 
    # any number in {1, ..., 7} indicates a positive weight
    # N indicates the weight is assumed 0
    # X (wildcard) indicates no assumption on the weight

# always, we assume that:
    # u = mu - nu, 
    # v = mu + nu, and 
    # mn = mu*nu

# all other variables are self-evident
# formulas are written to minimize FLOPs and accumulated 
# error when possible

function a2_assume234(a3, mn, v)
    
    a2num = 2 * a3 * (mn)^2
    a2denom = 2 * (mn + a3*v)^2 + a3^3*v
    
    return intersect(a2num / a2denom, UNIT_INT)
end  
    
function a4_assume1234(a3, mn, v)
    
    a4num = a3 * (2*(mn + a3*v)^2 + a3^3*v)^2
    a4denom = 4 * mn^2 * (mn + a3*v)^2
    a4denom -= 2*a3^3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom -= a3^5 * v * (2*mn + a3*v)
    
    return intersect(a4num / a4denom, UNIT_INT)
end

# can we avoid using this one? so complicated... 

function a4_assumeN2347(a3, mn, v)
    a4num = 4*((3*a3*v + mn)*(2*a3*v + mn) - a3*mn*v) * mn^2 * a3 
    a4num += 4*v* a3^4 * ((mn + a3*v)^2 + v^2 * (4*mn + a3*v))
    a4num += v^2 * a3^7

    a4denom = 4 * mn^2 * (mn + a3*v)^2
    a4denom -= 2 * a3^3 * (a3*v + mn) * (a3*v + 3*mn)*v
    a4denom -= a3^5 * v * (2*mn + a3*v)
    
    return intersect(a4num / a4denom, UNIT_INT)
end

function a4_assume12N4(a2, mn, v)
    
    a4num = 2*a2*mn^2
    a4denom = 2*(a2*v - mn)^2 - a2^3*v
    
    return intersect(a4num / a4denom, UNIT_INT)
end

function fg3_assume23(mu, nu, a3, mn, v, g_pos=true)
    
    f3num = intersect((a3+2*nu)*mu, -POS)
    g3num = intersect((a3+2*mu)*nu, -POS)
    denom = intersect(a3*v + 2*mn, -POS)
    
    f3 = intersect( sqrt( intersect(f3num / denom, UNIT_INT)), UNIT_INT)
    g3 = intersect( sqrt( intersect(g3num / denom, GEQ_ONE)), GEQ_ONE)
    
    if g_pos
        return f3, g3
    end
    return f3, -g3
end


function fg2_assume23(mu, nu, a3, f3, g3, g_pos=true)
    
    f2 = intersect((1+a3/mu)*f3, GEQ_ONE)
    g2 = ((1+a3/nu)*g3)
    
    if g_pos
        return f2, intersect(g2, UNIT_INT)
    end
    return f2, intersect(g2, -UNIT_INT)
end

    
function fg4_assume234(mu, nu, a2, f2, f3, g2, g3, g_pos = true)
    
    f4 = intersect(f3-a2*f2/mu, UNIT_INT)
    g4 = (g3-a2*g2/nu)
    
    if g_pos
        return f4, intersect(g4, GEQ_ONE)
    end
    return f4, intersect(g4, -GEQ_ONE)
end

function fg1_assume124(mu, nu, a4, f2, f4, g2, g4)
    
    f1 = intersect(f2+a4*f4/mu, GEQ_ONE)
    g1 = intersect(g2+a4*g4/nu, UNIT_INT)
    
    return f1, g1
end

function fg2_assume2N4(mu, nu, a2, mn, v, g_pos = true)
    
    f2num = intersect((a2-2*nu)*mu, POS)
    g2num = intersect((a2-2*mu)*nu, POS)
    denom = intersect(a2*v - 2*mn, POS)
    
    f2 = intersect( sqrt( intersect(f2num / denom, GEQ_ONE)), GEQ_ONE)
    g2 = intersect( sqrt( intersect(g2num / denom, UNIT_INT)), UNIT_INT)
    
    if g_pos
        return f2, g2
    end
    return f2, -g2
end

function fg4_assume2N4(mu, nu, a2, f2, g2, g_pos = true)
    
    f4 = intersect((1-a2/mu)*f2, UNIT_INT)
    g4 = ((1-a2/nu)*g2)
    
    if g_pos
        return f4, intersect(g4, GEQ_ONE)
    end
    return f4, intersect(g4, -GEQ_ONE)
end

# below are methods to directly determine feasibility
# as a convention:
    # avec is filled with the ai's, including 0's for the 
    # missing vertices
    # for fvec and gvec, a missing vertex is indicated by None


# checks that fvec and gvec can satisfy the eigen-equations

function fg_row_feasible(mu, nu, fvec, gvec, avec)

    J_INDICES_TO_COMPUTE = [
        1:7,
        [1,2,3,5,6,7],
        [1,2,5,6,7],
        [1,5,6,7],
        [1,2,3,4,5,6],
        [1,2,3,4,5],
        [1,2,3,4]
    ]

    for index = 1:7
        if fvec[index] != nothing
            fsum = interval(0)
            gsum = interval(0)
            for j = J_INDICES_TO_COMPUTE[index]
                if fvec[j] != nothing
                    fsum += avec[j]*fvec[j]
                    gsum += avec[j]*gvec[j]
                end
            end
            if isdisjoint(fsum, mu*fvec[index])
                return false
            end
            if isdisjoint(gsum, nu*gvec[index])
                return false
            end
        end
    end

    return true
end


# checks that the edge density is not less than 
# the sum of the squares of the two known eigenvalues

function density_feasible(mu, nu, avec)
    
    d = 1-(avec[3]+avec[4])^2-(avec[6]+avec[7])^2 - 2*(avec[2]*avec[4]+avec[5]*avec[7])
    
    if isdisjoint(d-mu^2-nu^2, POS)
        return false
    end

    return true
end

# checks that fvec and gvec can have norm 1

function norm_feasible(fvec, gvec, avec)
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i = 1:7
        if fvec[i] != nothing
            fnorm += avec[i]*fvec[i]^2
            gnorm += avec[i]*gvec[i]^2
        end
    end
    
    if isdisjoint(fnorm, interval(1))
        return false
    end
    
    if isdisjoint(gnorm, interval(1))
        return false
    end
    
    return true
end


# checks the ellipse equations can be satisfied

function ellipse_feasible(mu, nu, fvec, gvec, u)
    
    for i = 1:length(fvec)
        if fvec[i] != nothing
            if isdisjoint(mu*fvec[i]^2 - nu*gvec[i]^2, u)
                return false
            end
        end
    end
    
    return true
end
    


