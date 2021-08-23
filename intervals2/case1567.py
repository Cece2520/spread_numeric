# Case 1|567

from __future__ import print_function
from interval import *
from casework_helper import *

import datetime, platform, queue


# This program numerically attempts to rule out solutions to
# case 1|567 via interval arithmetic and a divide and conquer
# algorithm

def is_feasible(mu, nu, a1, a6):
    
    # We ignore cases that cannot exceed
    # the conjectured optimum of 2/sqrt(3),
    # cannot be less than the upper bound
    # (1+sqrt(2))/2, or cannot satisfy the 
    # inequality mu + sqrt(mu-mu^2) >= spread
    # Some helper variables are defined
    
    u = mu-nu
    v = mu+nu
    mn = mu*nu

    if not mu_nu_feasible(mu, nu, u):
        return False
    
    
    # We ignore cases where the weight sum exceeds 1
    # Apply the relevant formulas for a_i, f_i, g_i
    
    asum = (a1+a6) & UNIT_INT
    
    if asum == NULL_INT:
        return False
    
    a5 = a2_assume234(a6, mn, v)
    asum = (asum + a5) & UNIT_INT
    if asum == NULL_INT:
        return False
    
    a7 = a4_assume1234(a6, mn, v)
    asum = (asum + a7) & UNIT_INT
    if asum == NULL_INT:
        return False
    
    f6, g6 = fg3_assume23(mu, nu, a6, mn, v, g_pos = False)
    f5, g5 = fg2_assume23(mu, nu, a6, f6, g6, g_pos = False)
    f7, g7 = fg4_assume234(mu, nu, a5, f5, f6, g5, g6, g_pos = False)
    f1, g1 = fg1_assume124(mu, nu, a7, f5, f7, g5, g7)
    
    if f1 == NULL_INT:
        return False
    if g1 == NULL_INT:
        return False
    
    a1 = (1-asum) & UNIT_INT
    if a1 == NULL_INT:
        return False
    
    avec = [a1, 0, 0, 0, a5, a6, a7]
    
    
    # Edge density equals the sum of the squares of eigenvalues
    
    if not density_feasible(mu, nu, avec):
        return False
    
    
    # Check the eigenvector equations
    
    fvec = [f1, None, None, None, f5, f6, f7]
    gvec = [g1, None, None, None, g5, g6, g7]

    if not fg_ineq_feasible(fvec,gvec):
        return False
    
    if not fg_row_feasible1(mu, nu, fvec, gvec, avec):
        return False
    
    
    # Check the norm and ellipse equations
    
    if not norm_feasible1(fvec, gvec, avec):
        return False
    
    if not ellipse_feasible(mu, nu, fvec, gvec, u):
        return False
    
    return True


# Divide-and-Conquer Algorithm! We construct a grid over (mu, nu, a1, a6) 
# in the box [.65, 1] x [-.5, -.15] x [0, 1] x [0, 1]
# Our grid has initial stepsizes .05, .05, .1, .1, respectively
# We queue each box as a separate case, stored with depth term
# If a case is not shown to be infeasible, we divide it in half along 
# a dimension, queueing each half of the box
# The halved dimension is chosen according to the 
# congruence mod 4 of the depth


print("="*30, "Case 1|567", "="*30,"\n\n")

print("="*3, "System Information", "="*3)
uname = platform.uname()
print(f"System: {uname.system}")
print(f"Node Name: {uname.node}")
print(f"Release: {uname.release}")
print(f"Version: {uname.version}")
print(f"Machine: {uname.machine}")
print(f"Processor: {uname.processor}\n\n")

now = datetime.datetime.now()
print("="*3, "Initialization Date and Time", "="*3)
print ("Current Date and Time: ")
print (now.strftime("%Y-%m-%d %H:%M:%S"),"\n\n")

print("="*3, "Divide and Conquer!", "="*3)

case_queue = queue.Queue()

Mdenom = 20
Ndenom = 20
A1denom = 10
A6denom = 10

for M in range(7, 20):
    for N in range(-10, -3):
        for A1 in range(0, 10):
            for A6 in range(0, 10-A1):
                case_queue.put( (M,Mdenom, N,Ndenom, A1,A1denom, A6,A6denom, 0) )

curr_depth = -1
curr_size = 0
next_size = case_queue.qsize()

ctr = 0

print("Attempting Case 1|567")

while not case_queue.empty() and curr_depth < 50:
    (M,Mdenom, N,Ndenom, A1,A1denom, A6,A6denom, depth) = case_queue.get()
    if depth != curr_depth:
        curr_depth = depth
        curr_size = next_size
        ctr += curr_size
        next_size = 0
        print(f"\t Current Depth is {curr_depth:3} ... There are {curr_size:7} Boxes Remaining...")
    
    mu = interval[M, M+1] / interval(Mdenom)
    nu = interval[N, N+1] / interval(Ndenom)
    a1 = interval[A1, A1+1] / interval(A1denom)
    a6 = interval[A6, A6+1] / interval(A6denom)
    
    
    if is_feasible(mu, nu, a1, a6):
        next_size += 2
        
        if depth % 4 == 0:
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A1, 2*A1denom, A6, A6denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A1+1, 2*A1denom, A6, A6denom, depth+1) )
        
        if depth % 4 == 1:
            case_queue.put( (M,Mdenom, N,Ndenom, A1, A1denom, 2*A6, 2*A6denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, A1, A1denom, 2*A6+1, 2*A6denom, depth+1) )
        
        if depth % 4 == 2:
            case_queue.put( (2*M,2*Mdenom, N,Ndenom, A1, A1denom, A6, A6denom, depth+1) )
            case_queue.put( (2*M+1,2*Mdenom, N,Ndenom, A1, A1denom, A6, A6denom, depth+1) )
        
        if depth % 4 == 3:
            case_queue.put( (M,Mdenom, 2*N,2*Ndenom, A1, A1denom, A6, A6denom, depth+1) )
            case_queue.put( (M,Mdenom, 2*N+1,2*Ndenom, A1, A1denom, A6, A6denom, depth+1) )

curr_depth += 1
next_size = case_queue.qsize()

if not case_queue.empty():
    print(f"\t Current Depth is {curr_depth:3} ... There are {next_size:7} Boxes Remaining...\n")
    print(f"Case 1|567 is Not Yet Infeasible! A Total of {ctr} Boxes Were Considered...\n\n")
else:
    print(f"\t Current Depth is {curr_depth:3} ... There are {next_size:7} Boxes Remaining...\n")
    print(f"Case 1|567 is Infeasible! A Total of {ctr} Boxes Were Considered...\n\n")

print("="*3, "Termination Date and Time", "="*3)
now = datetime.datetime.now()
print ("Current Date and Time: ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))
