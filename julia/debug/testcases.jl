using IntervalArithmetic, DataStructures, Dates

include("case123457-test.jl")

function test_test_case(case, feas_func)
    print_specs(case)
    println("=== Divide and Conquer! ===")

    case_queue = Queue{NTuple{9, Int64}}()
    why_infeas = Dict()

    Mdenom = 20
    Ndenom = 20
    Alow_denom = 10
    Ahigh_denom = 10

    for M = 7:19
        for N = -10:-4
            for Alow in 0:9
                for Ahigh = 0:(9-Alow)
                    enqueue!(case_queue, (M,Mdenom, N,Ndenom, Alow,Alow_denom, Ahigh,Ahigh_denom, 0))
                end
            end
        end
    end

    curr_depth = -1
    curr_size = 0
    next_size = length(case_queue)
    MAX_DEPTH = 5

    ctr = 0

    println("Attempting Case $case")

    while !isempty(case_queue) && curr_depth < MAX_DEPTH
        (M,Mdenom, N,Ndenom, Alow,Alow_denom, Ahigh,Ahigh_denom, depth) = dequeue!(case_queue)
        if depth != curr_depth
            curr_depth = depth
            curr_size = next_size
            ctr += curr_size
            next_size = 0
            println("\t Current Depth is $(lpad(curr_depth,3))... There are $(lpad(curr_size,7)) Boxes Remaining...")
            println(why_infeas)
        end
        
        mu = interval(M, M+1) / interval(Mdenom)
        nu = interval(N, N+1) / interval(Ndenom)
        alow = interval(Alow, Alow+1) / interval(Alow_denom)
        ahigh = interval(Ahigh, Ahigh+1) / interval(Ahigh_denom)
        
        (is_feas, why) = feas_func(mu, nu, alow, ahigh)
        if haskey(why_infeas, why)
            why_infeas[why] += 1
        else
            why_infeas[why] = 1
        end


        if is_feas
            next_size += 2
            
            if depth % 4 == 0
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, 2*Alow, 2*Alow_denom, Ahigh, Ahigh_denom, depth+1) )
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, 2*Alow+1, 2*Alow_denom, Ahigh, Ahigh_denom, depth+1) )
            end
            
            if depth % 4 == 1
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, Alow, Alow_denom, 2*Ahigh, 2*Ahigh_denom, depth+1) )
                enqueue!(case_queue, (M,Mdenom, N,Ndenom, Alow, Alow_denom, 2*Ahigh+1, 2*Ahigh_denom, depth+1) )
            end
            
            if depth % 4 == 2
                enqueue!(case_queue, (2*M,2*Mdenom, N,Ndenom, Alow, Alow_denom, Ahigh, Ahigh_denom, depth+1) )
                enqueue!(case_queue, (2*M+1,2*Mdenom, N,Ndenom, Alow, Alow_denom, Ahigh, Ahigh_denom, depth+1) )
            end
            
            if depth % 4 == 3
                enqueue!(case_queue, (M,Mdenom, 2*N,2*Ndenom, Alow, Alow_denom, Ahigh, Ahigh_denom, depth+1) )
                enqueue!(case_queue, (M,Mdenom, 2*N+1,2*Ndenom, Alow, Alow_denom, Ahigh, Ahigh_denom, depth+1) )
            end
        end
    end

    curr_depth += 1
    next_size = length(case_queue)

    if !isempty(case_queue)
        print("\t Current Depth is $(lpad(curr_depth,3)) ... There are $(lpad(next_size,7)) Boxes Remaining...\n")
        println("Case $case is Not Yet Infeasible! A Total of $ctr Boxes Were Considered...")
    else
        print("\t Current Depth is $(lpad(curr_depth,3)) ... There are $(lpad(next_size,7)) Boxes Remaining...\n")
        println("Case $case is Infeasible! A Total of $ctr Boxes Were Considered...")
    end

    println("=== Termination Date and Time ===")
    now = Dates.now()
    println("Current Date and Time: $now")
end

test_test_case("TESTING", is_feasible_123457_test)