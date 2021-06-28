from ring_traffic_model import *
import cvxpy as cp
import mosek as msk


def lqr_sdp(N,s_star,gamma_s,gamma_v,gamma_u,AV_number):
    
    A, B, Q, R = ring_traffic_model(N,s_star,gamma_s,gamma_v,gamma_u,AV_number)
    
    n = 2*N
    m = AV_number
    
    epsilon = 1e-5
    
    H = np.identity(n)
    H[1:n:2][1:n:2] = 0
    
    
    X = cp.Variable((n,n), symmetric=True)
    Z = cp.Variable((m,n))
    Y = cp.Variable((m,m), symmetric=True)
    W = cp.Variable((m+n,m+n), symmetric=True)
    
    objective = cp.Minimize( cp.trace( Q@X ) + cp.trace( R@Y ) )  
    constraints = [-((A@X - B@Z) + (A@X - B@Z).T + H@(H).T) >= 0, 
                   X >= epsilon*np.identity(n),
                   W[0:m,0:m] == Y,
                   W[m:,0:m] == Z,
                   W[m:,m:] == X,
                   W >= 0]
    
    prob = cp.Problem(objective, constraints)
    
    prob.solve(solver=cp.MOSEK)
    Xd = X.value
    Zd = Z.value
    Pd = prob.value

    K = Zd*np.linalg.inv(Xd)   

    return K
    