using IntervalArithmetic, DataStructures
include("utils.jl")

# numerically attempts to rule out solutions to 
# the constraints, falling into the given intervals

function is_feasible_17_old(mu, nu, a1, a7)
    
    # first, we ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3)
    # some helper variables are also used
    
    u = mu-nu
    if isdisjoint(u - SPR_MAX, POS)
        return false
    end
    
    v = mu+nu
    mn = mu*nu
    
    
    # ignore cases where weight sum isn't 1
    
    asum = intersect(a1+a7, interval(1))
    
    if isempty(asum)
        return false
    end
    
    # formulas for fi's, gi's
    
    avec = [a1, 0, 0, 0, 0, 0, a7]
    
    
    # the sum of squares of eigenvalues equals 
    # the graph edge density
    
    if !density_feasible(mu, nu, avec)
        return false
    end
    
    
    # formulas for f1, g1, ..., f7, g7 must hold
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, false)
    f1, g1 = fg2_assume23(mu, nu, a7, f7, g7)
    
    if isempty(f1)
        return false
    end
    if isempty(g1)
        return false
    end
    
    
    # check eigenvector equations
    
    fvec = [f1, nothing, nothing, nothing, nothing, nothing, f7]
    gvec = [g1, nothing, nothing, nothing, nothing, nothing, g7]
    
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


# divide-and-conquer!  begin with a grid over (mu, nu, a1, a7)
# in the box [.65, 1] x [-.5, -.15] x [0, 1] x [0, 1]
# subdivide by the stepsizes .05, .05, .1, .1, respective
# queue up each box as a separate case, stored with depth term
# if a case cannot be ruled infeasible, split it in half along 
# one dimension, queueing each half of the box
# the halved dimension is chosen according to the 
# congruence mod 4 of the depth

function test_case17()
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
    next_size = 7*7*55
    MAX_DEPTH = 50

    ctr = 0

    print("trying case 1|7 ...")

    while !isempty(case_queue) && curr_depth < MAX_DEPTH
        (M,Mdenom, N,Ndenom, A1,A1denom, A7,A7denom, depth) = dequeue!(case_queue)
        if depth != curr_depth
            curr_depth = depth
            curr_size = next_size
            ctr += curr_size
            next_size = 0
            print("\ton depth = $curr_depth, size = $curr_size, so far $ctr ...\n")
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

    if isempty(case_queue)
        print("infeasible\n")
    else
        print("feasible\n")
    end
end

test_case17()