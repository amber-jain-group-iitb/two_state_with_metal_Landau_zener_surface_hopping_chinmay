import numpy as np

avr=np.zeros(2)

nfold=30
for i in range(1,nfold+1):
    avr=avr+np.loadtxt("./"+str(i)+"/fort.1400")

print("how often is pLZ>0.2:",(avr[0]/avr[1])*100,"%")


