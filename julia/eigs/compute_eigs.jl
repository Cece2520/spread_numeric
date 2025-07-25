
using IntervalArithmetic, IntervalLinearAlgebra


# Returns [max(intervals), min(intervals)] i.e. [mu, nu]
function find_extreme_intervals(intervals)
    upper = reduce((x,y) -> (sup(x) > sup(y) ? x : y), intervals)
    lower = reduce((x,y) -> (inf(x) < inf(y) ? x : y), intervals)
    return [upper, lower]
end


function comp_eigs_17(a1)
    a2 = interval(1) - a1
    A = [a1 a2; a1 0]
    eigs, eigvecs, cert = IntervalLinearAlgebra.verify_eigen(A)

    if cert[1]
        return find_extreme_intervals(real(eigs))
    else
        throw("Case 17: Unable to verify eigenvalues")
    end
end

function comp_eigs_147(a1, a2)
    a3 = interval(1) - a1 - a2
    A = [a1 a2 a3; a1 0 a3; a1 a2 0]
    eigs, eigvecs, cert = IntervalLinearAlgebra.verify_eigen(A)

    if cert[1]
        return find_extreme_intervals(real(eigs))
    else
        throw("Case 147: Unable to verify eigenvalues")
    end
end

function comp_eigs_457(a1, a2)
    a3 = interval(1) - a1 - a2
    A = [0 a2 a3; a1 a2 0; a1 0 0]
    eigs, eigvecs, cert = IntervalLinearAlgebra.verify_eigen(A)

    if cert[1]
        return find_extreme_intervals(real(eigs))
    else
        throw("Case 457: Unable to verify eigenvalues")
    end
end

a1 = interval(1)/interval(10)
a2 = interval(2)/interval(10)
eigs = comp_eigs_457(a1, a2)

println(eigs)