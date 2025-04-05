# Case 1|7

using IntervalArithmetic, DataStructures, Dates
include("utils.jl")

# import datetime, platform, queue


# This program numerically attempts to rule out solutions to
# case 1|7 via interval arithmetic and a divide and conquer
# algorithm

function is_feasible(mu, nu, a1, a7)
    
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
    
    asum = intersect(a1+a7, interval(1))
    
    if isempty(asum)
        return false
    end
    
    # Edge density equals the sum of the squares of eigenvalues

    avec = [a1, 0, 0, 0, 0, 0, a7]
    if !density_feasible(mu, nu, avec)
        return false
    end
    
    # Apply the relevant formulas for a_i, f_i, g_i
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, false)
    f1, g1 = fg2_assume23(mu, nu, a7, f7, g7)
    
    if isempty(f1)
        return false
    end
    if isempty(g1)
        return false
    end
    
    # Check the inequalities and eigenvector equations
    
    fvec = [f1, nothing, nothing, nothing, nothing, nothing, f7]
    gvec = [g1, nothing, nothing, nothing, nothing, nothing, g7]
    
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


# Divide-and-Conquer Algorithm! We construct a grid over (mu, nu, a1, a7) 
# in the box [.65, 1] x [-.5, -.15] x [0, 1] x [0, 1]
# Our grid has initial stepsizes .05, .05, .1, .1, respectively
# We queue each box as a separate case, stored with depth term
# If a case is not shown to be infeasible, we divide it in half along 
# a dimension, queueing each half of the box
# The halved dimension is chosen according to the 
# congruence mod 4 of the depth

function test17()

    print_specs("1|7")
    println("=== Divide and Conquer! ===")

    case_queue = Queue{NTuple{9, Int64}}()

    Mdenom = 20
    Ndenom = 20
    A1denom = 10
    A7denom = 10

    for M = 7:19
        for N = -10:-2
            for A1 in 0:9
                for A7 = 0:(9-A1)
                    enqueue!(case_queue, (M,Mdenom, N,Ndenom, A1,A1denom, A7,A7denom, 0))
                end
            end
        end
    end

    curr_depth = -1
    curr_size = 0
    next_size = length(case_queue)
    MAX_DEPTH = 50

    ctr = 0

    println("Attempting Case 1|7")

    while !isempty(case_queue) && curr_depth < MAX_DEPTH
        (M,Mdenom, N,Ndenom, A1,A1denom, A7,A7denom, depth) = dequeue!(case_queue)
        if depth != curr_depth
            curr_depth = depth
            curr_size = next_size
            ctr += curr_size
            next_size = 0
            println("\t Current Depth is $(lpad(curr_depth,3))... There are $(lpad(curr_size,7)) Boxes Remaining...")
        end
        
        mu = interval(M, M+1) / interval(Mdenom)
        nu = interval(N, N+1) / interval(Ndenom)
        a1 = interval(A1, A1+1) / interval(A1denom)
        a7 = interval(A7, A7+1) / interval(A7denom)
        
        
        if is_feasible(mu, nu, a1, a7)
            next_size += 2
            
            if depth % 4 == 0
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, 2*A1, 2*A1denom, A7, A7denom, depth+1) )
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, 2*A1+1, 2*A1denom, A7, A7denom, depth+1) )
            end
            
            if depth % 4 == 1
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, A1, A1denom, 2*A7, 2*A7denom, depth+1) )
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, A1, A1denom, 2*A7+1, 2*A7denom, depth+1) )
            end
            
            if depth % 4 == 2
                enqueue!(case_queue, (2*M,2*Mdenom, N,Ndenom, A1, A1denom, A7, A7denom, depth+1) )
                enqueue!(case_queue, (2*M+1,2*Mdenom, N,Ndenom, A1, A1denom, A7, A7denom, depth+1) )
            end
            
            if depth % 4 == 3
                enqueue!(case_queue, (M,Mdenom, 2*N,2*Ndenom, A1, A1denom, A7, A7denom, depth+1) )
                enqueue!(case_queue, (M,Mdenom, 2*N+1,2*Ndenom, A1, A1denom, A7, A7denom, depth+1) )
            end
        end
    end

    curr_depth += 1
    next_size = length(case_queue)

    if isempty(case_queue)
        print("\t Current Depth is $(lpad(curr_depth,3)) ... There are $(lpad(next_size,7)) Boxes Remaining...\n")
        println("Case 1|7 is Not Yet Infeasible! A Total of $ctr Boxes Were Considered...")
    else
        print("\t Current Depth is $(lpad(curr_depth,3)) ... There are $(lpad(next_size,7)) Boxes Remaining...\n")
        println("Case 1|7 is Infeasible! A Total of $ctr Boxes Were Considered...")
    end

    println("=== Termination Date and Time ===")
    now = Dates.now()
    println("Current Date and Time: $now")
end

test17()