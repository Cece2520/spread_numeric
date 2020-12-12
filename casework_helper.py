
from interval import interval, inf, imath, fpu

SPR_MAX = imath.sqrt(interval[4]/3)

NULL_INT = interval()
UNIT_INT = interval[0,1]
GEQ_ONE = interval[1, inf]
POS = interval[0,inf]


# these are helper methods meant to be used when 
# handling cases that share formulas 

# every method begins with the variable to be returned 
# and ends with a string indicating assumptions made
# most often, the string is a list of indices w/ weight
# assumed to be positive
# when weights are assumed 0, the string corresponds to 
# an interval in increasing order; in the interval: 
    # any number in {1, ..., 7} indicates a positive weight
    # N indicates the weight is assumed 0
    # X (wildcard) indicates no assumption on the weight


# always, we assume that:
    # u = mu - nu, 
    # v = mu + nu, and 
    # mn = mu*nu

# all other variables are self-evident
# formulas are written to minimize FLOPs and accumulated 
# error when possible

def a2_assume234(a3, mn, v):
    
    a2num = 2*a3*(mn)**2
    a2denom = 2*( mn + a3*v )**2 + a3**3*v
    
    return (a2num / a2denom) & UNIT_INT
    
    
def a4_assume1234(a3, mn, v):
    
    a4num = -a3*(2*(mn + a3*v)**2 + a3**3*v)**2
    a4denom = -4*mn**2*(mn + a3*v)**2
    a4denom += 2*a3**3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom += a3**5*v*( 2*mn + a3*v )
    
    return (a4num / a4denom) & UNIT_INT
    

# can we avoid using this one? so complicated... 


def a4_assumeN2347(a3, mn, v):
    a4num = -4*((3*a3*v + mn)*(2*a3*v + mn) - a3*mn*v)*mn**2*a3 
    a4num -= 4*v*a3**4*((mn + a3*v)**2 + v**2*(4*mn + a3*v))
    a4num -= v**2*a3**7

    a4denom = -4*mn**2*(mn + a3*v)**2
    a4denom += 2*a3**3 * (a3*v + mn)*(a3*v + 3*mn)*v
    a4denom += a3**5*v*( 2*mn + a3*v )
    
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
        return f2, g2 & GEQ_ONE
    return f2, g2 & (-GEQ_ONE)


# below are methods to directly determine feasibility
# as a convention:
    # avec is filled with the ai's, including 0's for the 
    # missing vertices
    # for fvec and gvec, a missing vertex is indicated by None


# checks that fvec and gvec can satisfy the eigen-equations

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


# checks that the edge density is not less than 
# the sum of the squares of the two known eigenvalues

def density_feasible(mu, nu, avec):
    
    d = 1-(avec[2]+avec[3])**2-(avec[5]+avec[6])**2 - 2*(avec[1]*avec[3]+avec[4]*avec[6])
    
    if (d-mu**2-nu**2) & POS == NULL_INT:
        return False
    return True


# checks that fvec and gvec can have norm 1

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


# checks the ellipse equations can be satisfied

def ellipse_feasible(mu, nu, fvec, gvec, u):
    
    for i in range(len(fvec)):
        if not fvec[i] == None:
            if (mu*fvec[i]**2-nu*gvec[i]**2) & u == NULL_INT:
                return False
    
    return True
    


