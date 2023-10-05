# Import libraries
import numpy as np
import inputs as inp
import constants as ct
import parameters as para
import matplotlib.pyplot as plt

# Define alternative equations
#def pend(y, t, b, c):
#    theta, omega = y
#    dydt = [omega, -b*omega - c*np.sin(theta)]
#    return dydt

def capsid_q_2(c0,t,k1,a2):
    c1, c2 = c0
    dc1dt = -2*k1*c1*c1 + a2*c2
    dc2dt = +k1*c1*c1 - a2*c2 
    dcdt = [dc1dt, dc2dt]

    return dcdt

def test(c0,t,v):
    c1, c2 = c0
    dc1dt = -2*v[0]*c1*c1 + v[1]*c2
    dc2dt = +v[0]*c1*c1 - v[1]*c2
    dcdt = [dc1dt, dc2dt]

    return dcdt


def capsid_q(q,c,k,a):
    # q: number of subunits
    # k: vector with the association rate constant of a monomer to each partial capsid
    # a: vector with dissociation rate constant of each partical capsid into a monmer and a partial capsid
    q=inp.q
 
    # Declare vector containing rates
    dcdt =[0]*q

    # Compact system of differential equations
    tempbnd=c[0]*np.sum(np.multiply(k[1:q-1],c[1:q-1]))
    tempunbnd=np.sum(np.multiply(a[2:],c[2:]))

    dcdt[0] = (-2*k[0]*c[0]*c[0]) + (2*a[1]*c[1]) -(tempbnd)  + (tempunbnd)
    for n in range(1,q-1):
        dcdt[n] = (c[0]*k[n-1]*c[n-1]) - (c[0]*k[n]*c[n]) - (a[n]*c[n]) + (a[n+1]*c[n+1])
    dcdt[q-1] = (c[0]*k[q-2]*c[q-2]) - (a[q-1]*c[q-1])
    # x=np.linspace(1,q+1,q)
    # plt.plot(x,dcdt)
    return dcdt
