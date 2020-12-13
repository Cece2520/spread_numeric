# case 1|4|7

from interval import interval, inf, imath, fpu
from casework_helper import *

import queue


# numerically attempts to rule out solutions to 
# the constraints, falling into the given intervals

def is_feasible(mu, nu, a4, a7):
    
    # first, we ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3)
    # some helper variables are also used
    
    u = mu-nu
    if ((u - SPR_MAX) & POS) == NULL_INT:
        return False
    
    v = mu+nu
    mn = mu*nu
    
    
    # ignore cases where weight sum exceeds 1
    
    asum = (a4+a7) & UNIT_INT
    
    if asum == NULL_INT:
        return False
    
    # again, weight sum cannot exceed 1; 
    # apply formulas for other ai's, fi's, gi's
    
    f4, g4 = fg3_assume23(mu, nu, a4, mn, v)
    f1, g1 = fg2_assume23(mu, nu, a4, f4, g4)
    
    if f4 == NULL_INT:
        return False
    if g4 == NULL_INT:
        return False
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, g_pos = False)
    f1bot, g1bot = fg2_assume23(mu, nu, a7, f7, g7)
    
    f1 = f1 & f1bot
    g1 = g1 & g1bot
    
    if f1 == NULL_INT:
        return False
    if g1 == NULL_INT:
        return False
    
    a1 = (1-(a4+a7)) & UNIT_INT
    
    # graph density equals sum of squares of eigenvalues
    
    avec = [a1, 0, 0, a4, 0, 0, a7]
    if not density_feasible(mu, nu, avec):
        return False
    
    
    # double-check the eigenvector eq'ns
    
    fvec = [f1, None, None, f4, None, None, f7]
    gvec = [g1, None, None, g4, None, None, g7]
    
    if not fg_row_feasible(mu, nu, fvec, gvec, avec):
        return False
    
    # might as well also check the norms and ellipse equations
    
    if not norm_feasible(fvec, gvec, avec):
        return False
    
    if not ellipse_feasible(mu, nu, fvec, gvec, u):
        return False
    
    return True


# divide-and-conquer!  begin with a grid over (mu, nu, a4, a7)
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
A4denom = 10
A7denom = 10

for M in range(7, 20):
    for N in range(-10, -3):
        for A4 in range(0, 10):
            for A7 in range(0, 10-A4):
                case_queue.put( (M,Mdenom, N,Ndenom, A4,A4denom, A7,A7denom, 0) )

curr_depth = -1
curr_size = 0
next_size = case_queue.qsize()

ctr = 0

print 'trying case 1|4|7 ...'

while not case_queue.empty():
    (M,Mdenom, N,Ndenom, A4,A4denom, A7,A7denom, depth) = case_queue.get()
    if depth != curr_depth:
        curr_depth = depth
        curr_size = next_size
        ctr += curr_size
        next_size = 0
        print '\ton depth =', curr_depth, '...', 'size =', curr_size, '...', 'so far', ctr, '...'
        
    
    mu = interval[(M+0.0)/Mdenom, (M+1.0)/Mdenom]
    nu = interval[(N+0.0)/Ndenom, (N+1.0)/Ndenom]
    a4 = interval[(A4+0.0)/A4denom, (A4+1.0)/A4denom]
    a7 = interval[(A7+0.0)/A7denom, (A7+1.0)/A7denom]
    
    if is_feasible(mu, nu, a4, a7):
        next_size += 2
        
        if depth % 4 == 0:
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A4, 2*A4denom, A7, A7denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A4+1, 2*A4denom, A7, A7denom, depth+1) )
        
        if depth % 4 == 1:
            case_queue.put( (M,Mdenom, N,Ndenom, A4, A7denom, 2*A7, 2*A7denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, A4, A7denom, 2*A7+1, 2*A7denom, depth+1) )
        
        if depth % 4 == 2:
            case_queue.put( (2*M,2*Mdenom, N,Ndenom, A4, A4denom, A7, A7denom, depth+1) )
            case_queue.put( (2*M+1,2*Mdenom, N,Ndenom, A4, A4denom, A7, A7denom, depth+1) )
        
        if depth % 4 == 3:
            case_queue.put( (M,Mdenom, 2*N,2*Ndenom, A4, A4denom, A7, A7denom, depth+1) )
            case_queue.put( (M,Mdenom, 2*N+1,2*Ndenom, A4, A4denom, A7, A7denom, depth+1) )

print 'infeasible\n'

