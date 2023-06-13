# nucleation



In this program we will be using scypi's odeint to solve for the concentrations of capsid formation in a solution as a functon of time.
We use a modular programming schema for our simulation and label the main file as Virca.py, short for Viral capsid. An environment.py module is provided that contains all the dependencies used in our program for ease of use. There are 5 other modules included for a total of 7. Known values such as the boltzman constant, pi and a standard temperature of 293.15 K are contanined in the constants.py module. Dimensions.py contains reference units for S.I. fundamental dimensions to scale our simulation and contains Length, Energy, time, and concentration scales.inputs.py is the module that houses required user inputsto specify the sysm: number of subunits in the viral capsid q, mean binding energy of subunits in the completed capsid Dg, Temperatuer of the system T, the effective radius of a subunit of the capsid Rs, the initial concentration of subunit proteins c0, the the estimated hydrodynamic radius of the free subunit Rh, the radius of the competed capsid Rc, as well as the viscoity of th esolution nu. within this module we convert these values into Si units and from them we derive the subunit diffusion value D1 and the effective diffusion of the subunit-partial capsid pair D_rel.
in the parameters.py module we determine the values of the binding and unbinding rates. We first define an interaction length dL as 2 times the hydrodynamic radius Rh. WE then define the total interaction area A between subunits for each partial capsid. then we define a binding kernel k as the flux of subunits through A. We use a placeholder value for the unbinding rate. mastereqs.py holds the equations for the change in concentrations over time for our partial capsids and holds the main functions called in Virca.py. 
Virca.py imports all the other modules, defines the overal time length and number of data points for the integrating schema used in odeint. The equations in mastereqs are then solved for our inputs and are output into a plot as concentration v time.

We are building a capsid of size q=72 subunits
with a critical concentration Ccmc=0.2*10**(-6)Moles
and a line tension alpha=2kT

The solution is water at room temperature 293k

The Viscosity of water is 0.9544*10**(-3)Pa

The hydrodynamic radius of our subunit is 0.619*10**(-6)m, estimated from a 300 amino acid protein.
