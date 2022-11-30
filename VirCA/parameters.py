# Parameter values for the model
## The parameters are determined from physical constants and physical inputs from the system.
import inputs as inp
import constants as ct
import matplotlib.pyplot as plt
import numpy as np
q = inp.q
D= inp.D
dL=inp.Rh*2###interaction distance defined as 2 times hydrodynamic radius of subunits 
##along interaction axis (units: m)
L=4*ct.pi*inp.Rc/q###interaction length perpindicular to interaction axis (units:m)
A=[L*dL*(np.sqrt(n*(q-n))) for n in range(q)] ###total interaction area between subunits for each partial capsid
## Binding rates
### Model based on effective Smoluchowski rate

k=[A[n] for n in range(q)] #Relative to partial capsid, how many subunits pass through interaction area A (units of counts/second)

#

plt.plot(k)
plt.title('binding kernel')
plt.ylabel('1/sec')
plt.xlabel('n')
plt.grid()
plt.show()
## Decay rates

a=np.exp(inp.Dg_hr)*(4*ct.pi*inp.Rc**2)

a[0]=0

# plt.plot(a)
# plt.title('unbinding')
# plt.ylabel('1/sec')
# plt.xlabel('n')
# plt.grid()
# plt.show()

