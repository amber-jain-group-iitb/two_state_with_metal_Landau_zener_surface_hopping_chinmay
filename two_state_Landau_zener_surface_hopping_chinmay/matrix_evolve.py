import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

inputs=np.genfromtxt("AFSSH.inp",delimiter="!!",dtype=str,autostrip=True)


     # temperature   #Gamma        #gh           #dG           #U
text=[inputs[18][0],inputs[20][0],inputs[25][0],inputs[26][0],inputs[29][0]]

ptext=[]
for tex in text:
    ptext.append(float(tex.replace(".d","E")))

KT=ptext[0]
Gamma=ptext[1]
Er=ptext[2]
dG=ptext[3]
U=ptext[4]
print(KT,Gamma,Er,dG,U)

def fermi(x):
  return 1/(1+np.exp(x/KT))

def forward(x,exo):
  return Gamma*(1-fermi(x))*np.exp(-(exo-x-Er)**2/(4*Er*KT))/np.sqrt(4*np.pi*Er*KT)

def backward(x,exo):
  return Gamma*fermi(x)*np.exp(-(exo+x-Er)**2/(4*Er*KT))/np.sqrt(4*np.pi*Er*KT)

def integrate_function_forward(a):
    result, error = integrate.quad(lambda x: forward(x, a), -0.1, 0.1)
    return result, error

def integrate_function_backward(a):
    result, error = integrate.quad(lambda x: backward(x, a), -0.1, 0.1)
    return result, error


kf=integrate_function_forward(-dG)[0]
kb=integrate_function_backward(dG)[0]

exo=[-dG-Er,dG+Er,Er-dG-U,dG+U-Er]
exoi=[0,dG+Er,dG+Er,2*dG+U]
k1=integrate_function_forward(exo[0])[0]
k2=integrate_function_backward(exo[1])[0]
k3=integrate_function_forward(exo[2])[0]
k4=integrate_function_backward(exo[3])[0]





total_time=600000
A=np.array([[-kf,kb],[kf,-kb]])
Ai=np.array([[-2*k1,k2,k2,0],[k1,-(k2+k3),0,k4],[k1,0,-(k3+k2),k4],[0,k3,k3,-2*k4]])

x=np.array([1,0])
pop=np.array([0,1,0,0])
dt=100


nsteps=int(total_time/dt)

U0_pop=[]
U1_pop=[]
time=[]
for i in range(0,nsteps):
  x=x+np.dot(A,x)*dt
  U0_pop.append(x[0])
  U1_pop.append(x[1])
  time.append(i*dt)

#plt.plot(time,U0_pop)
#plt.plot(time,U1_pop)
#plt.show()

pop1=[]
pop2=[]
pop3=[]
pop4=[]
time=[]
for i in range(0,nsteps):
  pop=pop+np.dot(Ai,pop)*dt
  pop1.append(pop[0])
  pop2.append(pop[1])
  pop3.append(pop[2])
  pop4.append(pop[3])
  time.append(i*dt)

#plt.plot(time,pop1,label='pop1')
#plt.plot(time,pop2,label='pop2')
#plt.plot(time,pop3,label='pop3')
#plt.plot(time,pop4,label='pop4')
#plt.legend()
#plt.show()
np.savetxt("pops4.txt", np.c_[time,pop1,pop2,pop3,pop4],delimiter="  ")

Qi=np.exp(-exoi[0]/KT)+np.exp(-exoi[1]/KT)+np.exp(-exoi[2]/KT)+np.exp(-exoi[3]/KT)
pi=[np.exp(-i/KT)/Qi for i in exoi]
print(pi)
