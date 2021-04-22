# nucleation


In this program we will be using scypi's odeint to solve for the concentrations of capsid formation in a solution as a functon of time.

# dCn/dt=J(n-1)-J(n)

We are building a capsid of size q=72 subunits
with a critical concentration Ccmc=0.2*10**(-6)Moles
and a line tension alpha=2kT

The solution is water at room temperature 293k

The Viscosity of water is 0.9544*10**(-3)Pa

The hydrodynamic radius of our subunit is 0.619*10**(-6)m, estimated from a 300 amino acid protein.
