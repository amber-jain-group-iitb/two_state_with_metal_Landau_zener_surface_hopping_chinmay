import numpy as np
import matplotlib.pyplot as plt

file="pop.out"
data=np.loadtxt(file)

plt.plot(data[1,:],data[0,:])

plt.show()
