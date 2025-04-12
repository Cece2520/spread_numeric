using IntervalArithmetic, DataStructures, Dates
include("utils.jl")

include("case17.jl")

include("case147.jl")
include("case157.jl")
include("case457.jl")

include("case1247.jl")
include("case1457.jl")
include("case1567.jl")
include("case2457.jl")
include("case4567.jl")

include("case12347.jl")
include("case12457.jl")
include("case14567.jl")
include("case24567.jl")

include("case123457.jl")
include("case124567.jl")
include("case234567.jl")

include("case1234567.jl")

# Divide-and-Conquer Algorithm! We construct a grid over (mu, nu, alow, ahigh)
# in the box [.65, 1] x [-.5, -.15] x [0, 1] x [0, 1]
# Our grid has initial stepsizes .05, .05, .1, .1, respectively
# We queue each box as a separate case, stored with depth term
# If a case is not shown to be infeasible, we divide it in half along 
# a dimension, queueing each half of the box
# The halved dimension is chosen according to the 
# congruence mod 4 of the depth

function test_case(case, feas_func)

    print_specs(case)
    println("=== Divide and Conquer! ===")

    case_queue = Queue{NTuple{9, Int64}}()

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
    MAX_DEPTH = 50

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
        end
        
        mu = interval(M, M+1) / interval(Mdenom)
        nu = interval(N, N+1) / interval(Ndenom)
        alow = interval(Alow, Alow+1) / interval(Alow_denom)
        ahigh = interval(Ahigh, Ahigh+1) / interval(Ahigh_denom)
        
        
        if feas_func(mu, nu, alow, ahigh)
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

    is_feasible = !isempty(case_queue)
    return is_feasible
end

CASES = (
    ("1|7", is_feasible_17),
    ("1|4|7", is_feasible_147),
    ("1|57", is_feasible_157),
    ("4|57", is_feasible_457),
    ("1|24|7", is_feasible_1247),
    ("1|4|57", is_feasible_1457),
    ("1|567", is_feasible_1567),
    ("24|57", is_feasible_2457),
    ("4|567", is_feasible_4567),
    ("1|234|7", is_feasible_12347),
    ("1|24|57", is_feasible_12457),
    ("1|4|567", is_feasible_14567),
    ("24|567", is_feasible_24567),
    ("1|234|57", is_feasible_123457),
    ("1|24|567", is_feasible_124567),
    ("234|567", is_feasible_234567),
    ("1|2|3|4|5|6|7", is_feasible_1234567)
)

is_case_feasible = zeros(length(CASES))
for ((case_name, case_func), i) = zip(CASES, 1:length(CASES))
    is_case_feasible[i] = test_case(case_name, case_func)
end

println()
println(" =========================== RESULTS =========================== ")
for i = 1:length(CASES)
    println("Case $(lpad(CASES[i][1],13)) is considered $(lpad(is_case_feasible[i] == 1 ? "FEASIBLE" : "INFEASIBLE", 10))")
end
println(" =============================================================== ")


# Replacement strings for changing py to jl

# if (?<arg>.*) == NULL_INT:
# if isempty($1)

# if not (?<arg>.*):
# if !$1

# = \((?<arg1>.*)\) & (?<arg2>.*)
# = intersect($1, $2)

# g_pos = 

# None -> nothing, False -> false, True -> true