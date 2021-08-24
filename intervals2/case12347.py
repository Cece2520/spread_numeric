# Case 1|234|7

from __future__ import print_function
from interval import *
from casework_helper import *

import datetime, platform, queue


# This program numerically attempts to rule out solutions to
# case 1|234|7 via interval arithmetic and a divide and conquer
# algorithm

def is_feasible(mu, nu, a3, a7):
    
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
    
    asum = (a3+a7) & UNIT_INT
    
    if asum == NULL_INT:
        return False

    a2 = a2_assume234(a3, mn, v)
    asum = (asum + a2) & UNIT_INT
    if asum == NULL_INT:
        return False
    
    a4 = a4_assume1234(a3, mn, v)
    asum = (asum + a4) & UNIT_INT
    if asum == NULL_INT:
        return False
    
    f3, g3 = fg3_assume23(mu, nu, a3, mn, v)
    f2, g2 = fg2_assume23(mu, nu, a3, f3, g3)
    f4, g4 = fg4_assume234(mu, nu, a2, f2, f3, g2, g3)
    f1, g1 = fg1_assume124(mu, nu, a4, f2, f4, g2, g4)
    
    if f1 == NULL_INT:
        return False
    if g1 == NULL_INT:
        return False
    
    a1 = (1-asum) & UNIT_INT
    avec = [a1, a2, a3, a4, 0, 0, a7]
    
    
    f7, g7 = fg3_assume23(mu, nu, a7, mn, v, g_pos = False)
    f1bot, g1bot = fg2_assume23(mu, nu, a3, f3, g3)
    f1 = f1 & f1bot
    g1 = g1 & g1bot
    

    # Edge density equals the sum of the squares of eigenvalues

    if not density_feasible(mu, nu, avec):
        return False
    
    
    # Check the inequalities and eigenvector equations
    
    fvec = [f1, f2, f3, f4, None, None, f7]
    gvec = [g1, g2, g3, g4, None, None, g7]

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


# Divide-and-Conquer Algorithm! We construct a grid over (mu, nu, a3, a7) 
# in the box [.65, 1] x [-.5, -.15] x [0, 1] x [0, 1]
# Our grid has initial stepsizes .05, .05, .1, .1, respectively
# We queue each box as a separate case, stored with depth term
# If a case is not shown to be infeasible, we divide it in half along 
# a dimension, queueing each half of the box
# The halved dimension is chosen according to the 
# congruence mod 4 of the depth


print("="*30, "Case 1|234|7", "="*30,"\n\n")

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
A3denom = 10
A7denom = 10

for M in range(7, 20):
    for N in range(-10, -3):
        for A3 in range(0, 10):
            for A7 in range(0, 10-A3):
                case_queue.put( (M,Mdenom, N,Ndenom, A3,A3denom, A7,A7denom, 0) )

curr_depth = -1
curr_size = 0
next_size = case_queue.qsize()

ctr = 0

print("Attempting Case 1|234|7")

while not case_queue.empty() and curr_depth < 50:
    (M,Mdenom, N,Ndenom, A3,A3denom, A7,A7denom, depth) = case_queue.get()
    if depth != curr_depth:
        curr_depth = depth
        curr_size = next_size
        ctr += curr_size
        next_size = 0
        print(f"\t Current Depth is {curr_depth:3} ... There are {curr_size:7} Boxes Remaining...")
    
    mu = interval[M, M+1] / interval(Mdenom)
    nu = interval[N, N+1] / interval(Ndenom)
    a3 = interval[A3, A3+1] / interval(A3denom)
    a7 = interval[A7, A7+1] / interval(A7denom)
    
    
    if is_feasible(mu, nu, a3, a7):
        next_size += 2
        
        if depth % 4 == 0:
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A3, 2*A3denom, A7, A7denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, 2*A3+1, 2*A3denom, A7, A7denom, depth+1) )
        
        if depth % 4 == 1:
            case_queue.put( (M,Mdenom, N,Ndenom, A3, A3denom, 2*A7, 2*A7denom, depth+1) )
            case_queue.put( (M,Mdenom, N,Ndenom, A3, A3denom, 2*A7+1, 2*A7denom, depth+1) )
        
        if depth % 4 == 2:
            case_queue.put( (2*M,2*Mdenom, N,Ndenom, A3, A3denom, A7, A7denom, depth+1) )
            case_queue.put( (2*M+1,2*Mdenom, N,Ndenom, A3, A3denom, A7, A7denom, depth+1) )
        
        if depth % 4 == 3:
            case_queue.put( (M,Mdenom, 2*N,2*Ndenom, A3, A3denom, A7, A7denom, depth+1) )
            case_queue.put( (M,Mdenom, 2*N+1,2*Ndenom, A3, A3denom, A7, A7denom, depth+1) )


curr_depth += 1
next_size = case_queue.qsize()

if not case_queue.empty():
    print(f"\t Current Depth is {curr_depth:3} ... There are {next_size:7} Boxes Remaining...\n")
    print(f"Case 1|234|7 is Not Yet Infeasible! A Total of {ctr} Boxes Were Considered...\n\n")
else:
    print(f"\t Current Depth is {curr_depth:3} ... There are {next_size:7} Boxes Remaining...\n")
    print(f"Case 1|234|7 is Infeasible! A Total of {ctr} Boxes Were Considered...\n\n")

print("="*3, "Termination Date and Time", "="*3)
now = datetime.datetime.now()
print ("Current Date and Time: ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))
