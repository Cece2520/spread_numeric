# case 1|57

from interval import interval, inf, imath, fpu
from casework_helper import *

import queue


# numerically attempts to rule out solutions to 
# the constraints, falling into the given intervals

def is_feasible(mu, nu, a1, a5):
    
    # first, we ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3)
    # some helper variables are also used
    
    u = mu-nu
    if ((u - SPR_MAX) & POS) == NULL_INT:
        return False
    
    v = mu+nu
    mn = mu*nu
    
    
    # ignore cases where weight sum exceeds 1
    # again, weight sum cannot exceed 1; 
    # apply formulas for other ai's, fi's, gi's
    
    f5, g5 = fg2_assume2N4(mu, nu, a5, mn, v, g_pos = False)
    f7, g7 = fg4_assume2N4(mu, nu, a5, f5, g5, g_pos = False)
    
    if f7 == NULL_INT:
        return False
    if g7 == NULL_INT:
        return False
    
    a7 = a4_assume12N4(a5, mn, v)
    asum = (a5+a7) & UNIT_INT
    a1 = a1 & (1-asum) & UNIT_INT
    
    if a1 == NULL_INT:
        return False
    
    f1, g1 = fg1_assume124(mu, nu, a7, f5, f7, g5, g7)
    
    # density equals sum of square of eigenvalues
    
    avec = [a1, 0, 0, 0, a5, 0, a7]
    if not density_feasible(mu, nu, avec):
        return False
    
    
    # double-check the eigenvector eq'ns
    
    fvec = [f1, None, None, None, f5, None, f7]
    gvec = [g1, None, None, None, g5, None, g7]
    
    if not fg_row_feasible(mu, nu, fvec, gvec, avec):
        return False
    
    
    # might as well also check the norms and ellipse equations
    
    if not norm_feasible(fvec, gvec, avec):
        return False
    
    if not ellipse_feasible(mu, nu, fvec, gvec, u):
        return False
    
    return True


# divide-and-conquer!  begin with a grid over (mu, nu, a1, a5)
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
A5denom = 10

for M in range(7, 20):
    for N in range(-10, -3):
        for A1 in range(0, 10):
            for A5 in range(0, 10-A1):
                case_queue.put( (M,Mdenom, N,Ndenom, A1,A1denom, A5,A5denom, 0) )

curr_depth = -1
curr_size = 0
next_size = case_queue.qsize()

ctr = 0

print('trying case 1|57 ...')

while not case_queue.empty() and curr_depth < MAX_DEPTH:
    (M,Mdenom, N,Ndenom, A1,A1denom, A5,A5denom, depth) = case_queue.get()
    if depth != curr_depth:
        curr_depth = depth
        curr_size = next_size
        ctr += curr_size
        next_size = 0
        print('\ton depth =', curr_depth, '...', 'size =', curr_size, '...', 'so far', ctr, '...')
        
    
    mu = interval[M, M+1] / interval(Mdenom)
    nu = interval[N, N+1] / interval(Ndenom)
    a1 = interval[A1, A1+1] / interval(A1denom)
    a5 = interval[A5, A5+1] / interval(A5denom)
    
    
    if is_feasible(mu, nu, a1, a5):
        next_size += 2
        
        if depth % 4 == 0:
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A1, 2*A1denom, A5, A5denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A1+1, 2*A1denom, A5, A5denom, depth+1) )
        
        if depth % 4 == 1:
            case_queue.put( (M,Mdenom, N,Ndenom, A1, A5denom, 2*A5, 2*A5denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, A1, A5denom, 2*A5+1, 2*A5denom, depth+1) )
        
        if depth % 4 == 2:
            case_queue.put( (2*M,2*Mdenom, N,Ndenom, A1, A1denom, A5, A5denom, depth+1) )
            case_queue.put( (2*M+1,2*Mdenom, N,Ndenom, A1, A1denom, A5, A5denom, depth+1) )
        
        if depth % 4 == 3:
            case_queue.put( (M,Mdenom, 2*N,2*Ndenom, A1, A1denom, A5, A5denom, depth+1) )
            case_queue.put( (M,Mdenom, 2*N+1,2*Ndenom, A1, A1denom, A5, A5denom, depth+1) )

if case_queue.empty():
    print('infeasible\n')
else:
    print('feasible\n')

