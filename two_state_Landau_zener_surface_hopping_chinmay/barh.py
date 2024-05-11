import numpy as np
import matplotlib.pyplot as plt

gh=10.0
dG=-0.006
U=0.01
mass=2000
omega=0.0002
KT=0.001

Er=0.5*mass*omega**2*gh**2

exo_ba=-(Er+dG)
exo_bc=dG+U-Er
h_ba=dG**2/(4*Er)
h_bc=(dG+U)**2/(4*Er)

print(h_ba,h_bc)
print(exo_ba,exo_bc)

def curve_A(x):
    return 0.5*mass*omega**2*x**2

def curve_B(x):
    return 0.5*mass*omega**2*(x**2-2*gh*x+2*gh**2)+dG

def curve_C(x):
    return 0.5*mass*omega**2*(x-2*gh)**2+2*dG+U


pos=np.linspace(-40,40,100)

#plt.plot(pos,[curve_A(i) for i in pos])
#plt.plot(pos,[curve_B(i) for i in pos])
#plt.plot(pos,[curve_C(i) for i in pos])
#plt.show()
exo_a=0.0
exo_b=Er+dG
exo_c=2*dG+U

Qt=np.exp(-exo_a/KT)+np.exp(-exo_b/KT)+np.exp(-exo_c/KT)

Pa=np.exp(-exo_a/KT)/Qt

Pb=np.exp(-exo_b/KT)/Qt

Pc=np.exp(-exo_c/KT)/Qt

print(Pa,Pb,Pc,Qt)
