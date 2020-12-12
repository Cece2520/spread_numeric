# case 1|234|57

from interval import interval, inf, imath, fpu
from casework_helper import *

import queue


# numerically attempts to rule out solutions to 
# the constraints, falling into the given intervals

def is_feasible(mu, nu, a3, a5):
    
    # first, we ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3)
    # some helper variables are also used
    
    u = mu-nu
    if ((u - SPR_MAX) & POS) == NULL_INT:
        return False
    
    v = mu+nu
    mn = mu*nu
    
    
    # ignore cases where weight sum exceeds 1
    
    asum = (a3+a5) & UNIT_INT
    
    if asum == NULL_INT:
        return False
    
    # again, weight sum cannot exceed 1; 
    # apply formulas for other ai's, fi's, gi's
    
    a2 = a2_assume234(a3, mn, v)
    asum = (asum + a2) & UNIT_INT
    if asum == NULL_INT:
        return False
    
    a4 = a4_assume1234(a3, mn, v)
    asum = (asum + a4) & UNIT_INT
    if asum == NULL_INT:
        return False
    
    f5, g5 = fg2_assume2N4(mu, nu, a5, mn, v, g_pos = False)
    f7, g7 = fg4_assume2N4(mu, nu, a5, f5, g5, g_pos = False)
    
    if f7 == NULL_INT:
        return False
    if g7 == NULL_INT:
        return False
    
    f3, g3 = fg3_assume23(mu, nu, a3, mn, v)
    f2, g2 = fg2_assume23(mu, nu, a3, f3, g3)
    f4, g4 = fg4_assume234(mu, nu, a2, f2, f3, g2, g3)
    f1, g1 = fg1_assume124(mu, nu, a4, f2, f4, g2, g4)
    
    if f1 == NULL_INT:
        return False
    if g1 == NULL_INT:
        return False
    
    a7 = a4_assume12N4(a5, mn, v)
    asum = (asum + a7) & UNIT_INT
    a1 = (1 - asum) & UNIT_INT
    if a1 == NULL_INT:
        return False
    
    avec = [a1, a2, a3, a4, a5, 0, a7]
    # double-check the eigenvector eq'ns
    
    fvec = [f1, f2, f3, f4, f5, None, f7]
    gvec = [g1, g2, f3, g4, g5, None, g7]
    
    if not fg_row_feasible(mu, nu, fvec, gvec, avec):
        return False
    
    
    # might as well also check the norms and ellipse equations
    
    if not norm_feasible(fvec, gvec, avec):
        return False
    
    if not ellipse_feasible(mu, nu, fvec, gvec, u):
        return False
    
    return True


# divide-and-conquer!  begin with a grid over (mu, nu, a3, a5)
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
A3denom = 10
A5denom = 10

for M in range(7, 20):
    for N in range(-10, -3):
        for A3 in range(0, 10):
            for A5 in range(0, 10-A3):
                case_queue.put( (M,Mdenom, N,Ndenom, A3,A3denom, A5,A5denom, 0) )

curr_depth = -1
curr_size = 0
next_size = 7*7*55

ctr = 0

while not case_queue.empty():
    (M,Mdenom, N,Ndenom, A3,A3denom, A5,A5denom, depth) = case_queue.get()
    if depth != curr_depth:
        curr_depth = depth
        curr_size = next_size
        ctr += curr_size
        next_size = 0
        print 'on depth =', curr_depth, '...', 'size =', curr_size, '...', 'so far', ctr, '...'
        
    
    mu = interval[(M+0.0)/Mdenom, (M+1.0)/Mdenom]
    nu = interval[(N+0.0)/Ndenom, (N+1.0)/Ndenom]
    a3 = interval[(A3+0.0)/A3denom, (A3+1.0)/A3denom]
    a5 = interval[(A5+0.0)/A5denom, (A5+1.0)/A5denom]
    
    if is_feasible(mu, nu, a3, a5):
        next_size += 2
        
        if depth % 4 == 0:
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A3, 2*A3denom, A5, A5denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A3+1, 2*A3denom, A5, A5denom, depth+1) )
        
        if depth % 4 == 1:
            case_queue.put( (M,Mdenom, N,Ndenom, A3, A5denom, 2*A5, 2*A5denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, A3, A5denom, 2*A5+1, 2*A5denom, depth+1) )
        
        if depth % 4 == 2:
            case_queue.put( (2*M,2*Mdenom, N,Ndenom, A3, A3denom, A5, A5denom, depth+1) )
            case_queue.put( (2*M+1,2*Mdenom, N,Ndenom, A3, A3denom, A5, A5denom, depth+1) )
        
        if depth % 4 == 3:
            case_queue.put( (M,Mdenom, 2*N,2*Ndenom, A3, A3denom, A5, A5denom, depth+1) )
            case_queue.put( (M,Mdenom, 2*N+1,2*Ndenom, A3, A3denom, A5, A5denom, depth+1) )

print 'done with case 1|234|57'

