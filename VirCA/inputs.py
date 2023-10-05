# List of inputs associated with the physical system to be modeled
import constants as ct
from scipy import constants as cts
import numpy as np

# Fundamental inputs in "human readable (hr)" (intuitive units) for capsid assembly
q_hr = 30		### Number of subunits forming the capsid.							Units: Dimensionless
Dg_hr = -4	### Mean binding energy of subunit full capsid. 					Units: kB*T0 where T0 is normal temperature (293.15 K).
T_hr = 20 		### Temperature of the system. 										Units: Celsius.
Rs_hr = 3.188 		### Subunit effective radius. 										Units: nm.
c0_hr = 10	### Initial concentration of subunits. 									Units: microMol/L
cs_hr=56 		### Solvent molarity 												Units:M (Mol/L)
Rh_hr=2.61 		### hydrodynamic radius for folded globular protein 				Units:(Ms2,mw=13.7kDa)
Rcap_hr=13.5 		### capsid radius (Ms2) 												Units: nM
### Note: Rational for working in nM: 1 molecule in a bacteria (1 um^3 in volume) would be 1.66 nM (so ~1 nM is an insightful unit)

### Warning!!! Water viscosity given at T=25C. Module needed to include temperature dependence.
nu_hr = 0.890 ### Water viscosity. Units: cP (centiPoise).

# Converstion of inputs to SI units
# These are the values used for consistency in operations.
import constants as ct
from scipy.constants import convert_temperature
import matplotlib.pyplot as plt

q = q_hr
Dg = Dg_hr*ct.kB*ct.T0 ## See Dg_hr. Units: Joules (J)
T = convert_temperature(T_hr, 'Celsius', 'Kelvin') ## See T_hr. Units: Kelvin (K)
Rs = Rs_hr*1E-9 ## See Rs_hr. Units: meters (m)
nu = nu_hr*1E-3 ## See nu_hr. Units: Pascal x seconds (Pa*s)
c0 = c0_hr*1E-6 *cts.N_A ## See c0. Units: 1/meter^3 (1/m^3)
cs=cs_hr#*cts.N_A ##See cs_hr. Units: 1/meter^3 (1/m^3)
Rh = Rh_hr*1E-9 ## See Rh_hr. Units: meters (m)
Rcap= Rcap_hr*1E-9 ## See Rc_hr. Units: meters (m)
# cq=c0#test for starting with populated cq

D1 = ct.kB*T/(6*ct.pi*nu*Rs) ### Subunit diffusion. Based on Stokes-Einstein equation. Units: m^2/sec

D=[D1/(i+1) for i in range(q)] ### Approximates diffusion for partial capsids decreasing by 1/n

D_rel=[D[n]+D[0]for n in range(q)] ### Approximates relative diffusion as D1+Dn

# DelG = [n*Dg_hr - n*np.log(c0/cs) + ct.pi/2*np.abs(Dg_hr)/np.sqrt(q)*np.sqrt(n*(q-n)) for n in range(q)] ### Total free energy of partial capsids

# def pentamers():
# 	q = q_hr/5
# 	Dg = Dg_hr*ct.kB*ct.T0*5 ## See Dg_hr. Units: Joules (J)
# 	T = convert_temperature(T_hr, 'Celsius', 'Kelvin') ## See T_hr. Units: Kelvin (K)
# 	Rs = 5.413*1E-9 ## See Rs_hr. Units: meters (m)
# 	nu = nu_hr*1E-3 ## See nu_hr. Units: Pascal x seconds (Pa*s)
# 	c0 = c0_hr*1E-6#*cts.N_A ## See c0. Units: 1/meter^3 (1/m^3)
# 	cs=cs_hr#*cts.N_A ##See cs_hr. Units: 1/meter^3 (1/m^3)
# 	Rh = 3.52*1E-9 ## See Rh_hr. Units: meters (m)
# 	Rcap= Rcap_hr*1E-9 ## See Rc_hr. Units: meters (m)
# 	return [q, Dg, T, Rs, nu, c0, cs, Rh, Rcap]
# def trimers():
# 	q = q_hr/5
# 	Dg = Dg_hr*ct.kB*ct.T0*5 ## See Dg_hr. Units: Joules (J)
# 	T = convert_temperature(T_hr, 'Celsius', 'Kelvin') ## See T_hr. Units: Kelvin (K)
# 	Rs = 4.101*1E-9 ## See Rs_hr. Units: meters (m)
# 	nu = nu_hr*1E-3 ## See nu_hr. Units: Pascal x seconds (Pa*s)
# 	c0 = c0_hr*1E-6#*cts.N_A ## See c0. Units: 1/meter^3 (1/m^3)
# 	cs=cs_hr#*cts.N_A ##See cs_hr. Units: 1/meter^3 (1/m^3)
# 	Rh = 2.98*1E-9 ## See Rh_hr. Units: meters (m)
# 	Rcap= Rcap_hr*1E-9 ## See Rc_hr. Units: meters (m)
# 	return[q,Dg,T,Rs,nu,c0,cs,Rh,Rcap]
# def dimers():
# 	q = q_hr/2
# 	Dg = Dg_hr*ct.kB*ct.T0*2 ## See Dg_hr. Units: Joules (J)
# 	T = convert_temperature(T_hr, 'Celsius', 'Kelvin') ## See T_hr. Units: Kelvin (K)
# 	Rs = 3.988*1E-9 ## See Rs_hr. Units: meters (m)
# 	nu = nu_hr*1E-3 ## See nu_hr. Units: Pascal x seconds (Pa*s)
# 	c0 = c0_hr*1E-6#*cts.N_A ## See c0. Units: 1/meter^3 (1/m^3)
# 	cs=cs_hr#*cts.N_A ##See cs_hr. Units: 1/meter^3 (1/m^3)
# 	Rh = 2.61*1E-9 ## See Rh_hr. Units: meters (m)
# 	Rcap= Rcap_hr*1E-9 ## See Rc_hr. Units: meters (m)
# 	return[q,Dg,T,Rs,nu,c0,cs,Rh,Rcap]
# def monomers():
# 	q = q_hr
# 	Dg = Dg_hr*ct.kB*ct.T0 ## See Dg_hr. Units: Joules (J)
# 	T = convert_temperature(T_hr, 'Celsius', 'Kelvin') ## See T_hr. Units: Kelvin (K)
# 	Rs = 3.165*1E-9 ## See Rs_hr. Units: meters (m)
# 	nu = nu_hr*1E-3 ## See nu_hr. Units: Pascal x seconds (Pa*s)
# 	c0 = c0_hr*1E-6#*cts.N_A ## See c0. Units: 1/meter^3 (1/m^3)
# 	cs = cs_hr#*cts.N_A ##See cs_hr. Units: 1/meter^3 (1/m^3)
# 	Rh = 2.08*1E-9 ## See Rh_hr. Units: meters (m)
# 	Rcap= Rcap_hr*1E-9 ## See Rc_hr. Units: meters (m)
# 	return[q,Dg,T,Rs,nu,c0,cs,Rh,Rcap]