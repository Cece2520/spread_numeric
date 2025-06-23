
using IntervalArithmetic, Dates

SPR_MAX = sqrt(interval(4)/3)
SPR_UPP = (interval(1) + sqrt(interval(2)))/interval(2)

NULL_INT = emptyinterval()
UNIT_INT = interval(0,1)
GEQ_ONE = interval(1, Inf)
POS = interval(0,Inf)

# These are helper methods meant to be used when 
# handling cases that share formulas 

# Every method begins with the variable to be returned 
# and ends with a string indicating assumptions made
# most often, the string is a list of indices w/ weight
# assumed to be positive
# when weights are assumed 0, the string corresponds to 
# an interval in increasing order; in the interval: 
    # any number in {1, ..., 7} indicates a positive weight
    # N indicates the weight is assumed 0
    # X (wildcard) indicates no assumption on the weight

# Always, we assume that:
    # u = mu - nu, 
    # v = mu + nu, and 
    # mn = mu*nu

# All other variables are self-evident
# Formulas are written to both minimize FLOPs and accumulated 
# error when possible


function a2_assume234(a3, mn, v)
    
    a2num = 2*a3*(mn)^2
    a2denom = 2*( mn + a3*v )^2 + a3^3*v
    
    return intersect(a2num / a2denom, UNIT_INT)
end
    
    
function a4_assume1234(a3, mn, v)
    
    a4num = a3*(2*(mn + a3*v)^2 + a3^3*v)^2
    a4denom = 4*mn^2*(mn + a3*v)^2
    a4denom -= 2*a3^3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom -= a3^5*v*( 2*mn + a3*v )
    
    return intersect(a4num / a4denom, UNIT_INT)
end


function a47_assume2457(a1,a2,a3,a5,a6,mu,nu,fvec,gvec)

    asum = intersect(interval(1)-a1-a2-a3-a5-a6, UNIT_INT)
    
    a4num_f = asum*fvec[7] - mu*(fvec[2]-fvec[5])
    a4num_g = asum*gvec[7] - nu*(gvec[2]-gvec[5])

    a7num_f = asum*fvec[4] - mu*(fvec[5]-fvec[2])
    a7num_g = asum*gvec[4] - nu*(gvec[5]-gvec[2])

    denom_f = fvec[4] + fvec[7]
    denom_g = gvec[4] + gvec[7]

    a4 = intersect(a4num_f / denom_f, a4num_g / denom_g)
    a7 = intersect(a7num_f / denom_f, a7num_g / denom_g)

    return intersect(a4, UNIT_INT), intersect(a7, UNIT_INT)
end
    

function a4_assumeN2347(a3, mn, v)
    a4num = 4*((3*a3*v + mn)*(2*a3*v + mn) - a3*mn*v)*mn^2*a3 
    a4num += 4*v*a3^4*((mn + a3*v)^2 + v^2*(4*mn + a3*v))
    a4num += v^2*a3^7

    a4denom = 4*mn^2*(mn + a3*v)^2
    a4denom -= 2*a3^3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom -= a3^5*v*( 2*mn + a3*v )
    
    return intersect(a4num / a4denom, UNIT_INT)
end


function a4_assume12N4(a2, mn, v)
    
    a4num = 2*a2*mn^2
    a4denom = 2*(a2*v - mn)^2 - a2^3*v
    
    return intersect(a4num / a4denom, UNIT_INT)
end


function fg3_assume23(mu, nu, a3, mn, c, v, s, g_pos = true)
    
    f3num = intersect(s*(a3+2*nu)*mu, -POS)
    g3num = intersect(s*(a3+2*mu)*nu, -POS)
    denom = intersect((mu-nu)*(a3*v + 2*mn), -POS)
    
    f3 = intersect( sqrt( intersect(f3num / (c * denom), UNIT_INT)), UNIT_INT)
    g3 = intersect( sqrt( intersect(g3num / ((interval(1)-c) * denom), GEQ_ONE)), GEQ_ONE)
    
    if g_pos
        return f3, g3
    end
    return f3, -g3
end


function fg2_assume23(mu, nu, a3, f3, g3, g_pos = true)
    
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


function fg2_assume2N4(mu, nu, a2, mn, c, v, s, g_pos = true)
    
    f2num = intersect(s*(a2-2*nu)*mu, POS)
    g2num = intersect(s*(a2-2*mu)*nu, POS)
    denom = intersect((mu-nu)*(a2*v - 2*mn), POS)
    
    f2 = intersect( sqrt( intersect(f2num / (c * denom), GEQ_ONE)), GEQ_ONE)
    g2 = intersect( sqrt( intersect(g2num / ((interval(1) - c) * denom), UNIT_INT)), UNIT_INT)
    
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


# Below are methods to directly determine feasibility
# As a convention:
    # avec is filled with the ai's, including 0's for the 
    # missing vertices
    # for fvec and gvec, a missing vertex is indicated by nothing


# Checks that mu and nu satisfy the necessary inequalities

function mu_nu_feasible(mu,nu,u,s,c)

    # if isdisjoint(u - SPR_MAX, POS)
    #     return false
    # end

    if isdisjoint(SPR_UPP - u, POS)
        return false
    end

    # if isdisjoint( sqrt(mu*(interval(1)-mu)) + nu, POS)
    #     return false
    # end

    if isdisjoint( sqrt(interval(1) - interval(2) * c + interval(2) * c^2) - s, POS)
        return false
    end

    return true
end


# Checks that fvec and gvec satisfy the necessary inequalities

function fg_ineq_feasible(fvec,gvec)

    if (fvec[1] != nothing) && (fvec[2] != nothing) && isdisjoint(fvec[1]-fvec[2], POS)
        return false
    end

    if (gvec[1] != nothing) && (gvec[2] != nothing) && isdisjoint(gvec[2]-gvec[1], POS)
        return false
    end

    if (fvec[3] != nothing) && (fvec[4] != nothing) && isdisjoint(fvec[3]-fvec[4], POS)
        return false
    end  

    if (gvec[3] != nothing) && (gvec[4] != nothing) && isdisjoint(gvec[4]-gvec[3], POS)
        return false
    end

    if (fvec[6] != nothing) && (fvec[7] != nothing) && isdisjoint(fvec[6]-fvec[7], POS)
        return false
    end

    if (gvec[6] != nothing) && (gvec[7] != nothing) && isdisjoint(gvec[6]-gvec[7], POS)
        return false
    end

    return true
end


# Checks that fvec and gvec can satisfy the eigen-equations

function fg_row_feasible(mu, nu, fvec, gvec, avec)

    J_INDICES_TO_COMPUTE = (
        (1,2,3,4,5,6,7),
        (1,2,3,5,6,7),
        (1,2,5,6,7),
        (1,5,6,7),
        (1,2,3,4,5,6),
        (1,2,3,4,5),
        (1,2,3,4)
    )
    
    for i = 1:7
        if fvec[i] != nothing
            fsum = interval(0)
            gsum = interval(0)
            for j = J_INDICES_TO_COMPUTE[i]
                if fvec[j] != nothing
                    fsum += avec[j]*fvec[j]
                    gsum += avec[j]*gvec[j]
                end
            end
            if isdisjoint(fsum, mu*fvec[i])
                return false
            end
            if isdisjoint(gsum, nu*gvec[i])
                return false
            end
        end
    end
    
    return true
end



# Checks that fvec and gvec can satisfy the eigen-equations
# Uses the equation a_1 = 1 - sum_{i \ne 1} a_i

function fg_row_feasible1(mu, nu, fvec, gvec, avec)

    J_INDICES_TO_COMPUTE = (
        (1,2,3,4,5,6,7),
        (1,2,3,5,6,7),
        (1,2,5,6,7),
        (1,5,6,7),
        (1,2,3,4,5,6),
        (1,2,3,4,5),
        (1,2,3,4)
    )

    J_INDICES_TO_COMPUTE_2 = (
        (2,3,4,5,6,7),
        (2,3,5,6,7),
        (2,5,6,7),
        (5,6,7),
        (2,3,4,5,6),
        (2,3,4,5),
        (2,3,4)
    )

    sum1_coeffs = (
        1,
        (1-avec[4]), ### not interval(1)?
        (1-avec[3]-avec[4]),
        (1-avec[2]-avec[3]-avec[4]),
        (1-avec[7]),
        (1-avec[6]-avec[7]),
        (1-avec[5]-avec[6]-avec[7])
    )
    
    for i = 1:7
        if fvec[i] != nothing
            fsum = interval(0)
            gsum = interval(0)
            for j = J_INDICES_TO_COMPUTE[i]
                if fvec[j] != nothing
                    fsum += avec[j]*fvec[j]
                    gsum += avec[j]*gvec[j]
                end
            end
            fsum1 = sum1_coeffs[i] * fvec[1]
            gsum1 = sum1_coeffs[i] * gvec[1]
            for j = J_INDICES_TO_COMPUTE_2[i]
                if fvec[j] != nothing
                    fsum1 += avec[j]*(fvec[j]-fvec[1])
                    gsum1 += avec[j]*(gvec[j]-gvec[1])
                end
            end
            if isdisjoint( intersect(fsum, fsum1), mu*fvec[i])
                return false
            end
            if isdisjoint( intersect(gsum, gsum1), nu*gvec[i])
                return false
            end
        end
    end
    
    return true
end


# Checks that fvec and gvec can satisfy the eigen-equations
# Uses the formulas a_4 = 1 - sum_{i \ne 4} a_i
# and a_7 = 1 - sum_{i \ne 7} a_i

function fg_row_feasible47(mu, nu, fvec, gvec, avec)

    J_INDICES_TO_COMPUTE = (
        (1,2,3,4,5,6,7),
        (1,2,3,5,6,7),
        (1,2,5,6,7),
        (1,5,6,7),
        (1,2,3,4,5,6),
        (1,2,3,4,5),
        (1,2,3,4)
    )

    J_INDICES_TO_COMPUTE_4 = (
        (1,2,3,5,6,7),
        (),
        (),
        (),
        (1,2,3,5,6),
        (1,2,3,5),
        (1,2,3)
    )

    J_INDICES_TO_COMPUTE_7 = (
        (1,2,3,4,5,6),
        (1,2,3,5,6),
        (1,2,5,6),
        (1,5,6),
        (),
        (),
        ()
    )

    sum4_coeffs = (
        1,
        0,
        0,
        0,
        1-avec[7],
        1-avec[6]-avec[7],
        1-avec[5]-avec[6]-avec[7]
    )

    sum7_coeffs = (
        1,
        1-avec[4],
        1-avec[3]-avec[4],
        1-avec[2]-avec[3]-avec[4],
        0,
        0,
        0
    )

    SHOULD_SUM4 = (1,0,0,0,1,1,1)
    SHOULD_SUM7 = (1,1,1,1,0,0,0)
    
    for i = 1:7
        if fvec[i] != nothing

            fsum = interval(0)
            gsum = interval(0)
            for j = J_INDICES_TO_COMPUTE[i]
                if fvec[j] != nothing
                    fsum += avec[j]*fvec[j]
                    gsum += avec[j]*gvec[j]
                end
            end
            if isdisjoint(fsum, mu*fvec[i]) || isdisjoint(gsum, nu*gvec[i])
                return false
            end

            if SHOULD_SUM4[i] == 1
                fsum4 = sum4_coeffs[i] * fvec[4]
                gsum4 = sum4_coeffs[i] * gvec[4]
                for j = J_INDICES_TO_COMPUTE_4[i]
                    if fvec[j] != nothing
                        fsum4 += avec[j]*(fvec[j]-fvec[4])
                        gsum4 += avec[j]*(gvec[j]-gvec[4])
                    end
                end
                if isdisjoint(fsum4, mu*fvec[i]) || isdisjoint(gsum4, nu*gvec[i])
                    return false
                end
            end
            
            if SHOULD_SUM7[i] == 1
                fsum7 = sum7_coeffs[i] * fvec[7]
                gsum7 = sum7_coeffs[i] * gvec[7]
                for j = J_INDICES_TO_COMPUTE_7[i]
                    if fvec[j] != nothing
                        fsum7 += avec[j]*(fvec[j]-fvec[7])
                        gsum7 += avec[j]*(gvec[j]-gvec[7])
                    end
                end
                if isdisjoint(fsum7, mu*fvec[i]) || isdisjoint(gsum7, nu*gvec[i])
                    return false
                end
            end
        end
    end
    
    return true
end


# Checks that the edge density is not less than 
# the sum of the squares of the two known eigenvalues

function density_feasible(mu, nu, avec)
    
    d = 1-(avec[3]+avec[4])^2-(avec[6]+avec[7])^2 - 2*(avec[2]*avec[4]+avec[5]*avec[7])
    
    if isdisjoint(d-mu^2-nu^2, POS)
        return false
    end

    return true
end


# Checks that fvec and gvec can have norm 1

function norm_feasible(fvec, gvec, avec)
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i = 1:7
        if fvec[i] != nothing
            fnorm += avec[i]*fvec[i]^2
            gnorm += avec[i]*gvec[i]^2
        end
    end
    
    if isdisjoint(fnorm, interval(1)) || isdisjoint(gnorm, interval(1))
        return false
    end
    
    return true
end


# Checks that fvec and gvec can have norm 1
# Uses the formula a_1 = 1 - sum_{i \ne 1} a_i

function norm_feasible1(fvec, gvec, avec)
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i = 1:7
        if fvec[i] != nothing
            fnorm += avec[i]*fvec[i]^2
            gnorm += avec[i]*gvec[i]^2
        end
    end

    fnorm1 = fvec[1]^2
    gnorm1 = gvec[1]^2
    
    for i = 2:7
        if fvec[i] != nothing
            fnorm1 += avec[i]*(fvec[i]^2 - fvec[1]^2)
            gnorm1 += avec[i]*(gvec[i]^2 - gvec[1]^2)
        end
    end
    
    if isdisjoint( intersect(fnorm, fnorm1), interval(1))
        return false
    end
    
    if isdisjoint( intersect(gnorm, gnorm1), interval(1))
        return false
    end
    
    return true
end


# Checks that fvec and gvec can satisfy the eigen-equations
# Uses the formulas a_4 = 1 - sum_{i \ne 4} a_i
# and a_7 = 1 - sum_{i \ne 7} a_i

function norm_feasible47(fvec, gvec, avec)
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i = 1:7
        if fvec[i] != nothing
            fnorm += avec[i]*fvec[i]^2
            gnorm += avec[i]*gvec[i]^2
        end
    end

    fnorm4 = fvec[4]^2
    gnorm4 = gvec[4]^2
    
    for i = (1,2,3,5,6,7)
        if fvec[i] != nothing
            fnorm4 += avec[i]*(fvec[i]^2 - fvec[4]^2)
            gnorm4 += avec[i]*(gvec[i]^2 - gvec[4]^2)
        end
    end
    
    fnorm7 = fvec[7]^2
    gnorm7 = gvec[7]^2
    
    for i = 1:6
        if fvec[i] != nothing
            fnorm7 += avec[i]*(fvec[i]^2 - fvec[7]^2)
            gnorm7 += avec[i]*(gvec[i]^2 - gvec[7]^2)
        end
    end

    if isdisjoint( intersect(fnorm, intersect(fnorm4, fnorm7)), interval(1))
        return false
    end
    
    if isdisjoint( intersect(gnorm, intersect(gnorm4, gnorm7)), interval(1))
        return false
    end
    
    return true
end


# Checks the ellipse equations can be satisfied

function ellipse_feasible(mu, nu, fvec, gvec, c, s)
    
    for i = 1:length(fvec)
        if fvec[i] != nothing && isdisjoint(c*mu*fvec[i]^2 - (interval(1)-c)*nu*gvec[i]^2, s)
            return false
        end
    end
    
    return true
end

function print_specs(case)
    println("============================== Case $case ==============================")

    ### TODO HALP ###
    # println("=== System Information ===")
    # uname = platform.uname()
    # println("System: $(Sys.KERNEL)")
    # println(f"Node Name: {uname.node}")
    # println(f"Release: {uname.release}")
    # println(f"Version: {uname.version}")
    # println("Machine: $(Sys.MACHINE)")
    # println(f"Processor: {uname.processor}\n\n")

    now = Dates.now()
    println("=== Initialization Date and Time ===")
    println("Current Date and Time: $now \n")
end