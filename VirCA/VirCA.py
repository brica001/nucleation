# Main VirCA's python script
# Date (start): 2022-03-31
# Date (end): xxxx-xx-xx

# Standard libraries
import numpy as np
import matplotlib.pyplot as plt
import time

# Physical constants
import constants as ct

# Physical inputs
import inputs as inp

# Reference dimensions
import dimensions as dim
L = dim.L ### Length
E = dim.E ### Energy
t = dim.t ### Time (use lower case to avoid confusion with Temperature)
C = dim.C ### Concentration
# Capsid geometry


# Model parameters
q = inp.q
# (Defined based on q, Deltag, T, subunit characteristics, capsid geometry...) 
import parameters as para
# Set initial concentrations
initCn = [0] * q
initCn[0] = inp.c0

# define time step vector
t_start=0
t_end=20
t_steps=10001
t = np.linspace(t_start, t_end, t_steps)

# Obtain dynamics
from mastereqs import capsid_q as syseq

def Con(c,t):
## Number of subunits 

	## Binding rates
	k=para.k

	## Decay rates
	a =para.a

	dcdt=syseq(q,c,k,a)
	for n in range(0,q):
		if dcdt[n] <1:
			if dcdt[n]<-1:
				dcdt[n] = inp.c0/(n+1)
				print(dcdt[n])
			else: 
				dcdt[n]=0
				print(dcdt[n])
	return dcdt
##RK4################
h=t_end/len(t)
st3=time.time()
i=1
y3=np.zeros((len(t),q))
Cq3=[0]
y3[0]=initCn
while i < len(t):
    k1=Con(y3[i-1,:],t)
    k1e=np.multiply(h,k1)
    k2=Con(y3[i-1]+k1e,t+h/2)
    k2e=np.multiply(h/2,k2)
    k3=Con(y3[i-1]+k2e,t+h/2)
    k3e=np.multiply(h,k3)
    k4=Con(y3[i-1]+k3e,t+h/2)

    tempk2 =np.multiply(2,k2)
    tempk3 =np.multiply(2,k3)
    tempest1= np.add(k1,tempk2)
    tempest2=np.add(tempk3,k4)
    tempest3=np.add(tempest1,tempest2)
    est=np.multiply(tempest3,h/6)

    y3[i,:]=np.add(y3[i-1],est)
    Cq3.append(y3[i][-1])
    i=i+1
et3=time.time()
et=et3-st3
print(f'rk4 done in {et}s!')

####AB-4#####
st4=time.time()
y4=np.zeros((len(t),q))
Cq4=[0]
y4[0]=initCn
i=0
while i<4:
    y4[i,:]=y3[i,:]
    Cq4.append(y4[i][-1])
    i=i+1
del y3
while i< len(t)-1:
    k1=Con(y4[i-1,:],t+3*h)
    k1=np.multiply(k1,55)
    k2=Con(y4[i-2,:],t+2*h),
    k2=np.multiply(k2,59)
    k3=Con(y4[i-3,:],t+h)
    k3=np.multiply(k3,37)
    k4=Con(y4[i-4,:],t)
    k4=np.multiply(k4,9)
    temp1=np.subtract(k1,k2)
    temp2=np.subtract(k3,k4)
    temp3=np.add(temp1,temp2)
    temp4=np.multiply(temp3,h/24)
    y4[i,:]=np.add(y4[i-1],temp4)
    Cq4.append(y4[i][-1])
    i=i+1
et4=time.time()
et=et4-st4
sol=y4
del y4
print(f'ab4 done in {et}s!')
print(np.shape(Cq4))
#from mastereqs import test
#sol = odeint(pend, y0, t, args=(b, c))
#sol = odeint(capsid_q_2, y0, t, args=param)
# sol=BDF(func,t_start,
# print(np.shape(sol))

# sol = odeint(func, c_t0, t)

# print(np.shape(sol))
# print(sol,np.ones(len(sol)))
# f= open('Cqm1dis.txt', 'w')
# for element in sol[:,1]:
#     f.write(str(element)+"\n")
# f.close()


### Total subunits remains constant check
# nvec = np.arange(1,q+1)
# ctot = np.zeros((len(Cq4)))
# for m in range(0,len(Cq4)):
#     ctot[m] = sum(np.multiply(nvec,Cq4[m]))

# plt.plot(t,ctot, label='ctot')


#####Plot solution
plt.plot(t,sol[:,0],'-.', label='c1')
plt.plot(t,sol[:,1],label='c2')
# # plt.plot(t,sol[:,2],label='c3')
# # plt.plot(t,sol[:,9],label='c10')
# # plt.plot(t,sol[:,29],label='c30')
# # plt.plot(t,sol[:,59],label='c60')
# n=1
# while n < q:
#     plt.plot(t, sol[:, n], label=f'c{n+1}(t)')
#     n=n+10

# plt.plot(t,sol[:,0],'-.', label='c1')
# # plt.plot(t,sol[:,1],label='c2')
# plt.plot(t, sol[:, -2],'g', label='c89')
plt.plot(t, sol[:, -1],'-.', label='cq')

# plt.legend(loc='center')
plt.xlabel('t')
plt.ylabel('Concentration [mols]')
plt.title('initial concentration v time')
plt.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
plt.legend(loc='lower right')
# plt.yscale('log')
plt.legend()
plt.grid()
plt.show()



quit()

