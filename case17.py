
from interval import interval, inf, imath, fpu
from casework_helper import *

import queue


# numerically attempts to rule out solutions to 
# the constraints, falling into the given intervals

def is_feasible(mu, nu, a1, a7):
    
    # first, we ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3)
    # some helper variables are also used
    
    u = mu-nu
    if ((u - SPR_MAX) & POS) == NULL_INT:
        return False
    
    v = mu+nu
    mn = mu*nu
    
    
    # ignore cases where weight sum isn't 1
    
    asum = (a1+a7) & interval(1)
    
    if asum == NULL_INT:
        return False
    
    
    # formulas for fi's, gi's
    
    avec = [a1, 0, 0, 0, 0, 0, a7]
    
    
    # the sum of squares of eigenvalues equals 
    # the graph edge density
    
    if not density_feasible(mu, nu, avec):
        return False
    
    
    # formulas for f1, g1, ..., f7, g7 must hold
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, g_pos = False)
    f1, g1 = fg2_assume23(mu, nu, a7, f7, g7)
    
    if f1 == NULL_INT:
        return False
    if g1 == NULL_INT:
        return False
    
    
    # check eigenvector equations
    
    fvec = [f1, None, None, None, None, None, f7]
    gvec = [g1, None, None, None, None, None, g7]
    
    if not fg_row_feasible(mu, nu, fvec, gvec, avec):
        return False
    
    
    # might as well also check the norms and ellipse equations
    
    if not norm_feasible(fvec, gvec, avec):
        return False
    
    if not ellipse_feasible(mu, nu, fvec, gvec, u):
        return False
    
    return True


# divide-and-conquer!  begin with a grid over (mu, nu, a1, a7)
# in the box [.65, 1] x [-.5, -.15] x [0, 1] x [0, 1]
# subdivide by the stepsizes .05, .05, .1, .1, respective
# queue up each box as a separate case, stored with depth term
# if a case cannot be ruled infeasible, split it in half along 
# one dimension, queueing each half of the box
# the halved dimension is chosen according to the 
# congruence mod 4 of the depth

case_queue = queue.Queue()

Mdenom = 20
Ndenom = 20
A1denom = 10
A7denom = 10

for M in range(7, 20):
    for N in range(-10, -3):
        for A1 in range(0, 10):
            for A7 in range(0, 10-A1):
                case_queue.put( (M,Mdenom, N,Ndenom, A1,A1denom, A7,A7denom, 0) )

curr_depth = -1
curr_size = 0
next_size = 7*7*55

ctr = 0

print 'trying case 1|7 ...'

while not case_queue.empty() and curr_depth < MAX_DEPTH:
    (M,Mdenom, N,Ndenom, A1,A1denom, A7,A7denom, depth) = case_queue.get()
    if depth != curr_depth:
        curr_depth = depth
        curr_size = next_size
        ctr += curr_size
        next_size = 0
        print '\ton depth =', curr_depth, '...', 'size =', curr_size, '...', 'so far', ctr, '...'
        
    
    mu = interval[M, M+1] / interval(Mdenom)
    nu = interval[N, N+1] / interval(Ndenom)
    a1 = interval[A1, A1+1] / interval(A1denom)
    a7 = interval[A7, A7+1] / interval(A7denom)
    
    
    if is_feasible(mu, nu, a1, a7):
        next_size += 2
        
        if depth % 4 == 0:
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A1, 2*A1denom, A7, A7denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A1+1, 2*A1denom, A7, A7denom, depth+1) )
        
        if depth % 4 == 1:
            case_queue.put( (M,Mdenom, N,Ndenom, A1, A1denom, 2*A7, 2*A7denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, A1, A1denom, 2*A7+1, 2*A7denom, depth+1) )
        
        if depth % 4 == 2:
            case_queue.put( (2*M,2*Mdenom, N,Ndenom, A1, A1denom, A7, A7denom, depth+1) )
            case_queue.put( (2*M+1,2*Mdenom, N,Ndenom, A1, A1denom, A7, A7denom, depth+1) )
        
        if depth % 4 == 3:
            case_queue.put( (M,Mdenom, 2*N,2*Ndenom, A1, A1denom, A7, A7denom, depth+1) )
            case_queue.put( (M,Mdenom, 2*N+1,2*Ndenom, A1, A1denom, A7, A7denom, depth+1) )

if case_queue.empty():
    print 'infeasible\n'
else:
    print 'feasible\n'

