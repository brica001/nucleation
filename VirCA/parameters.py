# Parameter values for the model
## The parameters are determined from physical constants and physical inputs from the system.
import inputs as inp
import constants as ct
import matplotlib.pyplot as plt
import numpy as np


Drel= inp.D_rel
I=5.2 ##Geometric factor, midsphere radius of icosah2edron

dL=inp.Rh*2###interaction distance defined as 2 times hydrodynamic radius of subunits 
###along interaction axis (units: m)
rim=[4*ct.pi * inp.Rcap / inp.q * np.sqrt(n*(inp.q-n)) for n in range(1,inp.q+1)] # rim length list of pcs

est=[(rim[n]/(inp.Rs)/I) for n in range(inp.q)] #estimation of number of subunits on rim

# Theta_eff=[2*(2*ct.pi - n/inp.q*np.pi) for n in range(inp.q)] ###Length inthe phi direction, minus interior (units*m)
# L=[4*ct.pi * inp.Rcap/q * (np.sqrt((n+1)*(q-n+1))) for n in range(q)] ###Rim length of partial capsids (units: m)
# Cap_area=np.multiply(Theta_eff,L) ###Total available area on the rim for capture through diffusion (units: m^2)

## Binding rates (1/s)
### Model based on effective Smoluchowski rate
k=[(est[n]/3* Drel[n] * dL * (1-np.exp(inp.Dg_hr/2))*ct.NA) for n in range(inp.q)] #how many subunits pass through interaction area (units:1/s)

#Decay rates (1/s)
a=[(est[n] * Drel[n] / (dL**2) * np.exp(inp.Dg_hr/2)) for n in range(1,inp.q-1)]
a.insert(0,0)
a.append(round(Drel[inp.q-1] / (dL**2) *np.exp(inp.Dg_hr)))


# x=np.linspace(1,inp.q,inp.q)
# plt.plot(x,k,label='Binding kernel')
# plt.plot(x,a,label='Unbinding kernel')
# plt.title('Steady state kernels')
# plt.ylabel('1/sec')
# plt.xlim(0,90)
# plt.xlabel('n')
# plt.legend()
# plt.show()
