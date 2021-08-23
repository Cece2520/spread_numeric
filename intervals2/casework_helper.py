
from interval import interval, inf, imath, fpu

SPR_MAX = imath.sqrt(interval[4]/3)
SPR_UPP = (interval(1)+imath.sqrt(interval(2)))/interval(2)

NULL_INT = interval()
UNIT_INT = interval[0,1]
GEQ_ONE = interval[1, inf]
POS = interval[0,inf]


# These are helper methods meant to be used when 
# handling cases that share formulas 

# Every method begins with the variable to be returned 
# and ends with a string indicating assumptions made
# most often, the string is a list of indices w/ weight
# assumed to be positive
# when weights are assumed 0, the string corresponds to 
# an interval in increasing order; in the interval: 
    # any number in {1, ..., 7} indicates a positive weight
    # N indicates the weight is assumed 0
    # X (wildcard) indicates no assumption on the weight

# Always, we assume that:
    # u = mu - nu, 
    # v = mu + nu, and 
    # mn = mu*nu

# All other variables are self-evident
# Formulas are written to both minimize FLOPs and accumulated 
# error when possible


def a2_assume234(a3, mn, v):
    
    a2num = 2*a3*(mn)**2
    a2denom = 2*( mn + a3*v )**2 + a3**3*v
    
    return (a2num / a2denom) & UNIT_INT
    
    
def a4_assume1234(a3, mn, v):
    
    a4num = a3*(2*(mn + a3*v)**2 + a3**3*v)**2
    a4denom = 4*mn**2*(mn + a3*v)**2
    a4denom -= 2*a3**3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom -= a3**5*v*( 2*mn + a3*v )
    
    return (a4num / a4denom) & UNIT_INT


def a47_assume2457(a1,a2,a3,a5,a6,mu,nu,fvec,gvec):

    asum = (interval[1]-a1-a2-a3-a5-a6) & UNIT_INT
    
    a4num_f = asum*fvec[6] - mu*(fvec[1]-fvec[4])
    a4num_g = asum*gvec[6] - nu*(gvec[1]-gvec[4])

    a7num_f = asum*fvec[3] - mu*(fvec[4]-fvec[1])
    a7num_g = asum*gvec[3] - nu*(gvec[4]-gvec[1])

    denom_f = fvec[3] + fvec[6]
    denom_g = gvec[3] + gvec[6]

    a4 = (a4num_f / denom_f) & (a4num_g / denom_g)
    a7 = (a7num_f / denom_f) & (a7num_g / denom_g)

    return a4 & UNIT_INT, a7 & UNIT_INT
    

def a4_assumeN2347(a3, mn, v):
    a4num = 4*((3*a3*v + mn)*(2*a3*v + mn) - a3*mn*v)*mn**2*a3 
    a4num += 4*v*a3**4*((mn + a3*v)**2 + v**2*(4*mn + a3*v))
    a4num += v**2*a3**7

    a4denom = 4*mn**2*(mn + a3*v)**2
    a4denom -= 2*a3**3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom -= a3**5*v*( 2*mn + a3*v )
    
    return (a4num / a4denom) & UNIT_INT


def a4_assume12N4(a2, mn, v):
    
    a4num = 2*a2*mn**2
    a4denom = 2*(a2*v - mn)**2 - a2**3*v
    
    return (a4num / a4denom) & UNIT_INT


def fg3_assume23(mu, nu, a3, mn, v, g_pos = True):
    
    f3num = ((a3+2*nu)*mu) & (-POS)
    g3num = ((a3+2*mu)*nu) & (-POS)
    denom = (a3*v + 2*mn) & (-POS)
    
    f3 = imath.sqrt((f3num / denom) & UNIT_INT) & UNIT_INT
    g3 = imath.sqrt((g3num / denom) & GEQ_ONE) & GEQ_ONE
    
    if g_pos:
        return f3, g3
    return f3, -g3


def fg2_assume23(mu, nu, a3, f3, g3, g_pos = True):
    
    f2 = ((1+a3/mu)*f3) & GEQ_ONE
    g2 = ((1+a3/nu)*g3)
    
    if g_pos:
        return f2, g2 & UNIT_INT
    return f2, g2 & (-UNIT_INT)

    
def fg4_assume234(mu, nu, a2, f2, f3, g2, g3, g_pos = True):
    
    f4 = (f3-a2*f2/mu) & UNIT_INT
    g4 = (g3-a2*g2/nu)
    
    if g_pos:
        return f4, g4 & GEQ_ONE
    return f4, g4 & (-GEQ_ONE)


def fg1_assume124(mu, nu, a4, f2, f4, g2, g4):
    
    f1 = (f2+a4*f4/mu) & GEQ_ONE
    g1 = (g2+a4*g4/nu) & UNIT_INT
    
    return f1, g1


def fg2_assume2N4(mu, nu, a2, mn, v, g_pos = True):
    
    f2num = ((a2-2*nu)*mu) & POS
    g2num = ((a2-2*mu)*nu) & POS
    denom = (a2*v - 2*mn) & POS
    
    f2 = imath.sqrt((f2num / denom) & GEQ_ONE) & GEQ_ONE
    g2 = imath.sqrt((g2num / denom) & UNIT_INT) & UNIT_INT
    
    if g_pos:
        return f2, g2
    return f2, -g2


def fg4_assume2N4(mu, nu, a2, f2, g2, g_pos = True):
    
    f4 = ((1-a2/mu)*f2) & UNIT_INT
    g4 = ((1-a2/nu)*g2)
    
    if g_pos:
        return f4, g4 & GEQ_ONE
    return f4, g4 & (-GEQ_ONE)


# Below are methods to directly determine feasibility
# As a convention:
    # avec is filled with the ai's, including 0's for the 
    # missing vertices
    # for fvec and gvec, a missing vertex is indicated by None


# Checks that mu and nu satisfy the necessary inequalities

def mu_nu_feasible(mu,nu,u):

    if ((u - SPR_MAX) & POS) == NULL_INT:
        return False

    if ((SPR_UPP - u) & POS) == NULL_INT:
        return False

    if ((imath.sqrt(mu*(interval(1)-mu)) + nu) & POS) == NULL_INT:
        return False

    return True


# Checks that fvec and gvec satisfy the necessary inequalities

def fg_ineq_feasible(fvec,gvec):

    if (not fvec[0] == None) and (not fvec[1] == None):
        if (fvec[0]-fvec[1]) & POS == NULL_INT:
            return False

    if (not gvec[0] == None) and (not gvec[1] == None):
        if (gvec[1]-gvec[0]) & POS == NULL_INT:
            return False

    if (not fvec[2] == None) and (not fvec[3] == None):
        if (fvec[2]-fvec[3]) & POS == NULL_INT:
            return False    

    if (not gvec[2] == None) and (not gvec[3] == None):
        if (gvec[3]-gvec[2]) & POS == NULL_INT:
            return False 

    if (not fvec[5] == None) and (not fvec[6] == None):
        if (fvec[5]-fvec[6]) & POS == NULL_INT:
            return False 

    if (not gvec[5] == None) and (not gvec[6] == None):
        if (gvec[5]-gvec[6]) & POS == NULL_INT:
            return False     

    return True


# Checks that fvec and gvec can satisfy the eigen-equations

def fg_row_feasible(mu, nu, fvec, gvec, avec):
    
    if not fvec[0] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in range(7):
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[0]) == NULL_INT:
            return False
        if gsum & (nu*gvec[0]) == NULL_INT:
            return False

    
    if not fvec[1] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[1]) == NULL_INT:
            return False
        if gsum & (nu*gvec[1]) == NULL_INT:
            return False

    
    if not fvec[2] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[2]) == NULL_INT:
            return False
        if gsum & (nu*gvec[2]) == NULL_INT:
            return False
    

    if not fvec[3] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[3]) == NULL_INT:
            return False
        if gsum & (nu*gvec[3]) == NULL_INT:
            return False
    

    if not fvec[4] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4,5]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[4]) == NULL_INT:
            return False
        if gsum & (nu*gvec[4]) == NULL_INT:
            return False
    

    if not fvec[5] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[5]) == NULL_INT:
            return False
        if gsum & (nu*gvec[5]) == NULL_INT:
            return False
    

    if not fvec[6] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        if fsum & (mu*fvec[6]) == NULL_INT:
            return False
        if gsum & (nu*gvec[6]) == NULL_INT:
            return False
    
    return True



# Checks that fvec and gvec can satisfy the eigen-equations
# Uses the equation a_1 = 1 - sum_{i \ne 1} a_i

def fg_row_feasible1(mu, nu, fvec, gvec, avec):
    
    if not fvec[0] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = fvec[0]
        gsum1 = gvec[0]
        for j in [1,2,3,4,5,6]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[0]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[0]) == NULL_INT:
            return False

    
    if not fvec[1] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = (1-avec[3])*fvec[0]
        gsum1 = (1-avec[3])*gvec[0]
        for j in [1,2,4,5,6]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[1]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[1]) == NULL_INT:
            return False

    
    if not fvec[2] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = (1-avec[2]-avec[3])*fvec[0]
        gsum1 = (1-avec[2]-avec[3])*gvec[0]
        for j in [1,4,5,6]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[2]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[2]) == NULL_INT:
            return False
    

    if not fvec[3] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = (1-avec[1]-avec[2]-avec[3])*fvec[0]
        gsum1 = (1-avec[1]-avec[2]-avec[3])*gvec[0]
        for j in [4,5,6]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[3]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[3]) == NULL_INT:
            return False
    

    if not fvec[4] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4,5]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = (1-avec[6])*fvec[0]
        gsum1 = (1-avec[6])*gvec[0]
        for j in [1,2,3,4,5]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[4]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[4]) == NULL_INT:
            return False
    

    if not fvec[5] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = (1-avec[5]-avec[6])*fvec[0]
        gsum1 = (1-avec[5]-avec[6])*gvec[0]
        for j in [1,2,3,4]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[5]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[5]) == NULL_INT:
            return False
    

    if not fvec[6] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum1 = (1-avec[4]-avec[5]-avec[6])*fvec[0]
        gsum1 = (1-avec[4]-avec[5]-avec[6])*gvec[0]
        for j in [1,2,3]:
            if not fvec[j] == None:
                fsum1 += avec[j]*(fvec[j]-fvec[0])
                gsum1 += avec[j]*(gvec[j]-gvec[0])
        if (fsum & fsum1) & (mu*fvec[6]) == NULL_INT:
            return False
        if (gsum & gsum1) & (nu*gvec[6]) == NULL_INT:
            return False

    
    return True


# Checks that fvec and gvec can satisfy the eigen-equations
# Uses the formulas a_4 = 1 - sum_{i \ne 4} a_i
# and a_7 = 1 - sum_{i \ne 7} a_i

def fg_row_feasible47(mu, nu, fvec, gvec, avec):
    
    if not fvec[0] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum4 = fvec[3]
        gsum4 = gvec[3]
        for j in [0,1,2,4,5,6]:
            if not fvec[j] == None:
                fsum4 += avec[j]*(fvec[j]-fvec[3])
                gsum4 += avec[j]*(gvec[j]-gvec[3])
        fsum7 = fvec[6]
        gsum7 = gvec[6]
        for j in [0,1,2,3,4,5]:
            if not fvec[j] == None:
                fsum7 += avec[j]*(fvec[j]-fvec[6])
                gsum7 += avec[j]*(gvec[j]-gvec[6])
        if (fsum & (fsum4 & fsum7)) & (mu*fvec[0]) == NULL_INT:
            return False
        if (gsum & (gsum4 & gsum7)) & (nu*gvec[0]) == NULL_INT:
            return False

    
    if not fvec[1] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum7 = (1-avec[3])*fvec[6]
        gsum7 = (1-avec[3])*gvec[6]
        for j in [0,1,2,4,5]:
            if not fvec[j] == None:
                fsum7 += avec[j]*(fvec[j]-fvec[6])
                gsum7 += avec[j]*(gvec[j]-gvec[6])
        if (fsum & fsum7) & (mu*fvec[1]) == NULL_INT:
            return False
        if (gsum & gsum7) & (nu*gvec[1]) == NULL_INT:
            return False

    
    if not fvec[2] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum7 = (1-avec[2]-avec[3])*fvec[6]
        gsum7 = (1-avec[2]-avec[3])*gvec[6]
        for j in [0,1,4,5]:
            if not fvec[j] == None:
                fsum7 += avec[j]*(fvec[j]-fvec[6])
                gsum7 += avec[j]*(gvec[j]-gvec[6])
        if (fsum & fsum7) & (mu*fvec[2]) == NULL_INT:
            return False
        if (gsum & gsum7) & (nu*gvec[2]) == NULL_INT:
            return False
    

    if not fvec[3] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,4,5,6]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum7 = (1-avec[1]-avec[2]-avec[3])*fvec[6]
        gsum7 = (1-avec[1]-avec[2]-avec[3])*gvec[6]
        for j in [0,4,5]:
            if not fvec[j] == None:
                fsum7 += avec[j]*(fvec[j]-fvec[6])
                gsum7 += avec[j]*(gvec[j]-gvec[6])
        if (fsum & fsum7) & (mu*fvec[3]) == NULL_INT:
            return False
        if (gsum & gsum7) & (nu*gvec[3]) == NULL_INT:
            return False
    

    if not fvec[4] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4,5]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum4 = (1-avec[6])*fvec[3]
        gsum4 = (1-avec[6])*gvec[3]
        for j in [0,1,2,4,5]:
            if not fvec[j] == None:
                fsum4 += avec[j]*(fvec[j]-fvec[3])
                gsum4 += avec[j]*(gvec[j]-gvec[3])
        if (fsum & fsum4) & (mu*fvec[4]) == NULL_INT:
            return False
        if (gsum & gsum4) & (nu*gvec[4]) == NULL_INT:
            return False
    

    if not fvec[5] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3,4]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum4 = (1-avec[5]-avec[6])*fvec[3]
        gsum4 = (1-avec[5]-avec[6])*gvec[3]
        for j in [0,1,2,4]:
            if not fvec[j] == None:
                fsum4 += avec[j]*(fvec[j]-fvec[3])
                gsum4 += avec[j]*(gvec[j]-gvec[3])
        if (fsum & fsum4) & (mu*fvec[5]) == NULL_INT:
            return False
        if (gsum & gsum4) & (nu*gvec[5]) == NULL_INT:
            return False
    

    if not fvec[6] == None:
        fsum = interval(0)
        gsum = interval(0)
        for j in [0,1,2,3]:
            if not fvec[j] == None:
                fsum += avec[j]*fvec[j]
                gsum += avec[j]*gvec[j]
        fsum4 = (1-avec[4]-avec[5]-avec[6])*fvec[3]
        gsum4 = (1-avec[4]-avec[5]-avec[6])*gvec[3]
        for j in [0,1,2]:
            if not fvec[j] == None:
                fsum4 += avec[j]*(fvec[j]-fvec[3])
                gsum4 += avec[j]*(gvec[j]-gvec[3])
        if (fsum & fsum4) & (mu*fvec[6]) == NULL_INT:
            return False
        if (gsum & gsum4) & (nu*gvec[6]) == NULL_INT:
            return False
    
    return True


# Checks that the edge density is not less than 
# the sum of the squares of the two known eigenvalues

def density_feasible(mu, nu, avec):
    
    d = 1-(avec[2]+avec[3])**2-(avec[5]+avec[6])**2 - 2*(avec[1]*avec[3]+avec[4]*avec[6])
    
    if (d-mu**2-nu**2) & POS == NULL_INT:
        return False

    return True


# Checks that fvec and gvec can have norm 1

def norm_feasible(fvec, gvec, avec):
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i in range(7):
        if not fvec[i] == None:
            fnorm += avec[i]*fvec[i]**2
            gnorm += avec[i]*gvec[i]**2
    
    if fnorm & interval(1) == NULL_INT:
        return False
    
    if gnorm & interval(1) == NULL_INT:
        return False
    
    return True


# Checks that fvec and gvec can have norm 1
# Uses the formula a_1 = 1 - sum_{i \ne 1} a_i

def norm_feasible1(fvec, gvec, avec):
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i in [0,1,2,3,4,5,6]:
        if not fvec[i] == None:
            fnorm += avec[i]*fvec[i]**2
            gnorm += avec[i]*gvec[i]**2

    fnorm1 = fvec[0]**2
    gnorm1 = gvec[0]**2
    
    for i in [1,2,3,4,5,6]:
        if not fvec[i] == None:
            fnorm1 += avec[i]*(fvec[i]**2 - fvec[0]**2)
            gnorm1 += avec[i]*(gvec[i]**2 - gvec[0]**2)
    
    if (fnorm & fnorm1) & interval(1) == NULL_INT:
        return False
    
    if (gnorm & gnorm1) & interval(1) == NULL_INT:
        return False
    
    return True


# Checks that fvec and gvec can satisfy the eigen-equations
# Uses the formulas a_4 = 1 - sum_{i \ne 4} a_i
# and a_7 = 1 - sum_{i \ne 7} a_i

def norm_feasible47(fvec, gvec, avec):
    
    fnorm = interval(0)
    gnorm = interval(0)
    
    for i in [0,1,2,3,4,5,6]:
        if not fvec[i] == None:
            fnorm += avec[i]*fvec[i]**2
            gnorm += avec[i]*gvec[i]**2

    fnorm4 = fvec[3]**2
    gnorm4 = gvec[3]**2
    
    for i in [0,1,2,4,5,6]:
        if not fvec[i] == None:
            fnorm4 += avec[i]*(fvec[i]**2 - fvec[3]**2)
            gnorm4 += avec[i]*(gvec[i]**2 - gvec[3]**2)
    
    fnorm7 = fvec[6]**2
    gnorm7 = gvec[6]**2
    
    for i in [0,1,2,3,4,5]:
        if not fvec[i] == None:
            fnorm7 += avec[i]*(fvec[i]**2 - fvec[6]**2)
            gnorm7 += avec[i]*(gvec[i]**2 - gvec[6]**2)

    if (fnorm & (fnorm4 & fnorm7)) & interval(1) == NULL_INT:
        return False
    
    if (gnorm & (gnorm4 & gnorm7)) & interval(1) == NULL_INT:
        return False
    
    return True


# Checks the ellipse equations can be satisfied

def ellipse_feasible(mu, nu, fvec, gvec, u):
    
    for i in range(len(fvec)):
        if not fvec[i] == None:
            if (mu*fvec[i]**2-nu*gvec[i]**2) & u == NULL_INT:
                return False
    
    return True
    


